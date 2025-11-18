"""
Pipeline Runner for Phredator

Executes the complete QC workflow: parse → analyze → fix → execute → verify
Shows before/after comparisons and improvement metrics.
"""

import os
import json
import subprocess
import tempfile
import shutil
from typing import Dict, Any, List, Optional
from dataclasses import dataclass, field, asdict
from pathlib import Path

from phredator.parser.fastqc_parser import FastQCParser
from phredator.analyzer.qc_analyzer import Analyzer
from phredator.fixer.qc_fixer import Fixer


@dataclass
class PipelineStep:
    """Represents a single pipeline step"""
    name: str
    status: str  # "pending", "running", "success", "failed", "skipped"
    command: Optional[str] = None
    output: Optional[str] = None
    error: Optional[str] = None
    duration_seconds: Optional[float] = None


@dataclass
class PipelineComparison:
    """Before/after comparison metrics"""
    metric: str
    before_status: str
    after_status: str
    before_value: Any
    after_value: Any
    improved: bool
    description: str


@dataclass
class PipelineResult:
    """Complete pipeline execution result"""
    input_file: str
    output_file: str
    steps: List[Dict[str, Any]] = field(default_factory=list)
    comparisons: List[Dict[str, Any]] = field(default_factory=list)
    overall_improvement: bool = False
    fixes_executed: int = 0
    metrics_improved: int = 0
    metrics_degraded: int = 0
    
    def to_json(self) -> str:
        """Convert to JSON string"""
        return json.dumps(asdict(self), indent=2)


class PipelineRunner:
    """
    Runs the complete Phredator pipeline with verification.
    
    Workflow:
    1. Parse initial FastQC
    2. Analyze and generate fixes
    3. Execute fixes (if --execute)
    4. Run FastQC on fixed file
    5. Parse and analyze fixed results
    6. Compare before/after metrics
    """
    
    def __init__(
        self,
        input_fastq: str,
        output_dir: str = "pipeline_output",
        organism: Optional[str] = None,
        experiment_type: Optional[str] = None,
        check_tools: bool = True,
        dry_run: bool = False,
        verbose: bool = False
    ):
        """
        Initialize pipeline runner.
        
        Args:
            input_fastq: Path to input FASTQ file (must have FastQC results)
            output_dir: Directory for pipeline outputs
            organism: Organism profile (human, mouse, etc.)
            experiment_type: Experiment type (wgs, rnaseq, etc.)
            check_tools: Check tool availability before running
            dry_run: If True, show what would be done without executing
            verbose: Verbose output
        """
        self.input_fastq = input_fastq
        self.output_dir = Path(output_dir)
        self.organism = organism
        self.experiment_type = experiment_type
        self.check_tools = check_tools
        self.dry_run = dry_run
        self.verbose = verbose
        
        # Derived paths
        self.sample_name = Path(input_fastq).stem.replace('.fastq', '').replace('.fq', '')
        self.fastqc_dir = self._find_fastqc_dir()
        
        # Output paths
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.before_parsed = self.output_dir / "before_parsed.json"
        self.before_analysis = self.output_dir / "before_analysis.json"
        self.fixes_json = self.output_dir / "fixes.json"
        self.fixed_fastq = self.output_dir / f"{self.sample_name}_fixed.fastq.gz"
        self.after_fastqc_dir = self.output_dir / f"{self.sample_name}_fixed_fastqc"
        self.after_parsed = self.output_dir / "after_parsed.json"
        self.after_analysis = self.output_dir / "after_analysis.json"
        self.comparison_json = self.output_dir / "comparison.json"
        
        # Pipeline state
        self.steps: List[PipelineStep] = []
        self.result: Optional[PipelineResult] = None
    
    def _find_fastqc_dir(self) -> Optional[Path]:
        """Find FastQC directory for input file"""
        base = Path(self.input_fastq).parent
        patterns = [
            f"{self.sample_name}_fastqc",
            f"{self.sample_name}_fastqc.zip"
        ]
        
        for pattern in patterns:
            candidate = base / pattern
            if candidate.exists():
                return candidate
        
        return None
    
    def _log(self, message: str):
        """Log message if verbose"""
        if self.verbose:
            print(f"[PIPELINE] {message}")
    
    def _add_step(self, step: PipelineStep):
        """Add step to pipeline"""
        self.steps.append(step)
        if self.verbose:
            status_icon = {"success": "✓", "failed": "✗", "skipped": "⊘", "running": "→"}.get(step.status, "•")
            print(f"{status_icon} {step.name}: {step.status.upper()}")
    
    def run(self) -> PipelineResult:
        """Execute the complete pipeline"""
        self._log(f"Starting pipeline for: {self.input_fastq}")
        
        # Step 1: Parse initial FastQC
        self._parse_initial_fastqc()
        
        # Step 2: Analyze initial data
        self._analyze_initial()
        
        # Step 3: Generate fixes
        self._generate_fixes()
        
        # Step 4: Execute fixes (if not dry-run)
        if not self.dry_run:
            self._execute_fixes()
            
            # Step 5: Run FastQC on fixed file
            self._run_fastqc_on_fixed()
            
            # Step 6: Parse fixed FastQC
            self._parse_fixed_fastqc()
            
            # Step 7: Analyze fixed data
            self._analyze_fixed()
            
            # Step 8: Compare results
            self._compare_results()
        else:
            self._add_step(PipelineStep(
                name="Execute Fixes",
                status="skipped",
                output="Dry-run mode - fixes not executed"
            ))
        
        # Create result
        self.result = PipelineResult(
            input_file=self.input_fastq,
            output_file=str(self.fixed_fastq) if not self.dry_run else "N/A (dry-run)",
            steps=[asdict(step) for step in self.steps]
        )
        
        # Save comparison
        comparison_path = self.output_dir / "pipeline_result.json"
        with open(comparison_path, 'w') as f:
            f.write(self.result.to_json())
        
        return self.result
    
    def _parse_initial_fastqc(self):
        """Parse initial FastQC results"""
        step = PipelineStep(name="Parse Initial FastQC", status="running")
        
        try:
            if not self.fastqc_dir:
                raise FileNotFoundError(f"FastQC directory not found for {self.input_fastq}")
            
            parser = FastQCParser(str(self.fastqc_dir))
            parsed = parser.parse()
            
            with open(self.before_parsed, 'w') as f:
                f.write(parsed.to_json())
            
            step.status = "success"
            step.output = str(self.before_parsed)
        except Exception as e:
            step.status = "failed"
            step.error = str(e)
        
        self._add_step(step)
        if step.status == "failed":
            raise RuntimeError(f"Failed to parse initial FastQC: {step.error}")
    
    def _analyze_initial(self):
        """Analyze initial parsed data"""
        step = PipelineStep(name="Analyze Initial QC", status="running")
        
        try:
            analyzer = Analyzer(
                str(self.before_parsed),
                organism=self.organism,
                experiment_type=self.experiment_type
            )
            analysis = analyzer.run()
            
            with open(self.before_analysis, 'w') as f:
                f.write(analysis.to_json())
            
            step.status = "success"
            step.output = str(self.before_analysis)
        except Exception as e:
            step.status = "failed"
            step.error = str(e)
        
        self._add_step(step)
        if step.status == "failed":
            raise RuntimeError(f"Failed to analyze initial data: {step.error}")
    
    def _generate_fixes(self):
        """Generate fix suggestions"""
        step = PipelineStep(name="Generate Fixes", status="running")
        
        try:
            fixer = Fixer(
                str(self.before_analysis),
                input_reads=self.input_fastq,
                check_tools=self.check_tools
            )
            fixes = fixer.run()
            
            with open(self.fixes_json, 'w') as f:
                f.write(fixes.to_json())
            
            step.status = "success"
            step.output = f"Generated {len(fixes.fixes_applied)} fixes"
            
            if self.result:
                self.result.fixes_executed = len(fixes.fixes_applied)
        except Exception as e:
            step.status = "failed"
            step.error = str(e)
        
        self._add_step(step)
        if step.status == "failed":
            raise RuntimeError(f"Failed to generate fixes: {step.error}")
    
    def _execute_fixes(self):
        """Execute the first (highest priority) fix"""
        step = PipelineStep(name="Execute Fixes", status="running")
        
        try:
            # Load fixes
            with open(self.fixes_json, 'r') as f:
                fixes_data = json.load(f)
            
            fixes = fixes_data.get('fixes_applied', [])
            
            if not fixes:
                step.status = "skipped"
                step.output = "No fixes to execute"
                self._add_step(step)
                return
            
            # Execute first fix (highest priority)
            fix = fixes[0]
            command = fix['command']
            
            # Replace placeholder with actual paths
            command = command.replace(fixes_data['input_file'], self.input_fastq)
            command = command.replace(f"{fixes_data['sample_name']}_trimmed.fastq.gz", str(self.fixed_fastq))
            command = command.replace(f"{fixes_data['sample_name']}_R1_trimmed.fastq.gz", str(self.fixed_fastq))
            
            step.command = command
            self._log(f"Executing: {command}")
            
            # Execute
            result = subprocess.run(
                command,
                shell=True,
                capture_output=True,
                text=True,
                timeout=600  # 10 minute timeout
            )
            
            if result.returncode == 0:
                step.status = "success"
                step.output = f"Executed: {fix['description']}"
            else:
                step.status = "failed"
                step.error = result.stderr
                
        except Exception as e:
            step.status = "failed"
            step.error = str(e)
        
        self._add_step(step)
        if step.status == "failed":
            raise RuntimeError(f"Failed to execute fixes: {step.error}")
    
    def _run_fastqc_on_fixed(self):
        """Run FastQC on the fixed file"""
        step = PipelineStep(name="Run FastQC on Fixed File", status="running")
        
        try:
            command = f"fastqc {self.fixed_fastq} -o {self.output_dir}"
            step.command = command
            
            self._log(f"Running: {command}")
            
            result = subprocess.run(
                command,
                shell=True,
                capture_output=True,
                text=True,
                timeout=600
            )
            
            if result.returncode == 0:
                step.status = "success"
                step.output = f"FastQC results: {self.after_fastqc_dir}"
            else:
                step.status = "failed"
                step.error = result.stderr
                
        except Exception as e:
            step.status = "failed"
            step.error = str(e)
        
        self._add_step(step)
        if step.status == "failed":
            raise RuntimeError(f"Failed to run FastQC: {step.error}")
    
    def _parse_fixed_fastqc(self):
        """Parse fixed FastQC results"""
        step = PipelineStep(name="Parse Fixed FastQC", status="running")
        
        try:
            parser = FastQCParser(str(self.after_fastqc_dir))
            parsed = parser.parse()
            
            with open(self.after_parsed, 'w') as f:
                f.write(parsed.to_json())
            
            step.status = "success"
            step.output = str(self.after_parsed)
        except Exception as e:
            step.status = "failed"
            step.error = str(e)
        
        self._add_step(step)
        if step.status == "failed":
            raise RuntimeError(f"Failed to parse fixed FastQC: {step.error}")
    
    def _analyze_fixed(self):
        """Analyze fixed data"""
        step = PipelineStep(name="Analyze Fixed QC", status="running")
        
        try:
            analyzer = Analyzer(
                str(self.after_parsed),
                organism=self.organism,
                experiment_type=self.experiment_type
            )
            analysis = analyzer.run()
            
            with open(self.after_analysis, 'w') as f:
                f.write(analysis.to_json())
            
            step.status = "success"
            step.output = str(self.after_analysis)
        except Exception as e:
            step.status = "failed"
            step.error = str(e)
        
        self._add_step(step)
        if step.status == "failed":
            raise RuntimeError(f"Failed to analyze fixed data: {step.error}")
    
    def _compare_results(self):
        """Compare before and after results"""
        step = PipelineStep(name="Compare Results", status="running")
        
        try:
            # Load analyses
            with open(self.before_analysis, 'r') as f:
                before = json.load(f)
            with open(self.after_analysis, 'r') as f:
                after = json.load(f)
            
            comparisons = []
            improved = 0
            degraded = 0
            
            # Compare each metric
            metrics = [
                'per_base_quality',
                'per_sequence_quality',
                'gc_content',
                'adapter_content',
                'duplication_levels',
                'overrepresented_sequences'
            ]
            
            for metric in metrics:
                before_metric = before['metrics'].get(metric, {})
                after_metric = after['metrics'].get(metric, {})
                
                if not before_metric or not after_metric:
                    continue
                
                before_status = before_metric.get('status', 'unknown')
                after_status = after_metric.get('status', 'unknown')
                
                # Check if improved
                status_order = {'pass': 3, 'warn': 2, 'fail': 1, 'unknown': 0}
                is_improved = status_order[after_status] > status_order[before_status]
                is_degraded = status_order[after_status] < status_order[before_status]
                
                if is_improved:
                    improved += 1
                elif is_degraded:
                    degraded += 1
                
                comparison = PipelineComparison(
                    metric=metric,
                    before_status=before_status,
                    after_status=after_status,
                    before_value=before_metric.get('summary', 'N/A'),
                    after_value=after_metric.get('summary', 'N/A'),
                    improved=is_improved,
                    description=f"{metric.replace('_', ' ').title()}: {before_status} → {after_status}"
                )
                comparisons.append(comparison)
            
            # Update result
            if self.result:
                self.result.comparisons = [asdict(c) for c in comparisons]
                self.result.metrics_improved = improved
                self.result.metrics_degraded = degraded
                self.result.overall_improvement = improved > degraded
            
            # Save comparison
            with open(self.comparison_json, 'w') as f:
                json.dump({
                    'before': before,
                    'after': after,
                    'comparisons': [asdict(c) for c in comparisons],
                    'summary': {
                        'improved': improved,
                        'degraded': degraded,
                        'unchanged': len(comparisons) - improved - degraded,
                        'overall_improvement': improved > degraded
                    }
                }, f, indent=2)
            
            step.status = "success"
            step.output = f"Improved: {improved}, Degraded: {degraded}"
            
        except Exception as e:
            step.status = "failed"
            step.error = str(e)
        
        self._add_step(step)
        if step.status == "failed":
            raise RuntimeError(f"Failed to compare results: {step.error}")
