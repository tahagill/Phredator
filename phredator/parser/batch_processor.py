

import os
import json
from typing import List, Dict, Any, Optional
from pathlib import Path
from dataclasses import dataclass, asdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime

from phredator.parser.fastqc_parser import FastQCParser
from phredator.analyzer.qc_analyzer import Analyzer
from phredator.fixer.qc_fixer import Fixer


@dataclass
class BatchSampleResult:
    """Result for a single sample in batch processing"""
    sample_name: str
    fastq_path: str
    status: str  # "success", "failed", "skipped"
    error_message: Optional[str] = None
    
    # Analysis results
    overall_status: Optional[str] = None
    issues_found: int = 0
    fixes_suggested: int = 0
    
    # Paths
    parsed_json: Optional[str] = None
    analysis_json: Optional[str] = None
    fixes_json: Optional[str] = None
    
    def to_dict(self) -> Dict[str, Any]:

        return asdict(self)


@dataclass
class BatchStatistics:
    """Aggregate statistics across all samples"""
    pass_count: int
    warn_count: int
    fail_count: int
    
    gc_mean: Optional[float] = None
    gc_std: Optional[float] = None
    
    quality_mean: Optional[float] = None
    quality_std: Optional[float] = None
    
    duplication_mean: Optional[float] = None
    duplication_std: Optional[float] = None
    
    outliers: List[Dict[str, Any]] = None
    
    def __post_init__(self):
        if self.outliers is None:
            self.outliers = []
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class BatchReport:
    """Overall batch processing report"""
    total_samples: int
    successful: int
    failed: int
    skipped: int
    start_time: str
    end_time: str
    duration_seconds: float
    organism: Optional[str] = None
    experiment_type: Optional[str] = None
    
    statistics: Optional[BatchStatistics] = None
    sample_results: List[BatchSampleResult] = None
    
    def __post_init__(self):
        if self.sample_results is None:
            self.sample_results = []
    
    def to_dict(self) -> Dict[str, Any]:

        result = {
            "total_samples": self.total_samples,
            "successful": self.successful,
            "failed": self.failed,
            "skipped": self.skipped,
            "start_time": self.start_time,
            "end_time": self.end_time,
            "duration_seconds": self.duration_seconds,
            "organism": self.organism,
            "experiment_type": self.experiment_type,
            "sample_results": [s.to_dict() for s in self.sample_results]
        }
        
        if self.statistics:
            result["statistics"] = self.statistics.to_dict()
        
        return result


class BatchProcessor:
    """Process multiple samples in batch with optional parallelization"""
    
    def __init__(self, 
                 sample_list: List[str],
                 output_dir: str,
                 organism: Optional[str] = None,
                 experiment_type: Optional[str] = None,
                 check_tools: bool = True,
                 parallel: int = 1,
                 dry_run: bool = False,
                 verbose: bool = False):
        pass  # docstring removed
        self.sample_list = sample_list
        self.output_dir = Path(output_dir)
        self.organism = organism
        self.experiment_type = experiment_type
        self.check_tools = check_tools
        self.parallel = parallel
        self.dry_run = dry_run
        self.verbose = verbose
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.results: List[BatchSampleResult] = []
    
    def process_sample(self, fastq_path: str, sample_idx: int) -> BatchSampleResult:
        pass  # docstring removed
        sample_name = Path(fastq_path).stem
        
        if self.verbose:
            print(f"[{sample_idx}/{len(self.sample_list)}] Processing: {sample_name}")
        sample_dir = self.output_dir / sample_name
        sample_dir.mkdir(exist_ok=True)
        
        result = BatchSampleResult(
            sample_name=sample_name,
            fastq_path=fastq_path,
            status="processing"
        )
        
        try:
            # Step 1: Check if FastQC data exists
            fastqc_dir = self._find_fastqc_dir(fastq_path)
            if not fastqc_dir:
                result.status = "skipped"
                result.error_message = "No FastQC data found - run FastQC first"
                if self.verbose:
                    print(f"  ⚠️  {result.error_message}")
                return result
            
            # Step 2: Parse FastQC
            parsed_path = sample_dir / "parsed.json"
            parser = FastQCParser(fastqc_dir)
            parsed_report = parser.parse()
            with open(parsed_path, 'w') as f:
                f.write(parsed_report.to_json())
            result.parsed_json = str(parsed_path)
            
            if self.verbose:
                print(f"  ✓ Parsed FastQC data")
            
            # Step 3: Analyze with profiles
            analysis_path = sample_dir / "analysis.json"
            analyzer = Analyzer(
                str(parsed_path),
                organism=self.organism,
                experiment_type=self.experiment_type
            )
            analysis_result = analyzer.run()
            
            with open(analysis_path, 'w') as f:
                f.write(analysis_result.to_json())
            result.analysis_json = str(analysis_path)
            result.overall_status = analysis_result.overall_status
            
            # Count issues
            result.issues_found = sum(
                1 for metric in analysis_result.metrics.values()
                if metric.get('status') in ['WARN', 'FAIL']
            )
            
            if self.verbose:
                print(f"  ✓ Analysis complete: {result.overall_status} ({result.issues_found} issues)")
            
            # Step 4: Generate fixes
            fixes_path = sample_dir / "fixes.json"
            fixer = Fixer(
                str(analysis_path),
                input_reads=fastq_path,
                check_tools=self.check_tools
            )
            fixes_result = fixer.run()
            
            with open(fixes_path, 'w') as f:
                f.write(fixes_result.to_json())
            result.fixes_json = str(fixes_path)
            result.fixes_suggested = len(fixes_result.fixes_applied)
            
            if self.verbose:
                print(f"  ✓ Generated {result.fixes_suggested} fix suggestions")
            
            result.status = "success"
            
        except Exception as e:
            result.status = "failed"
            result.error_message = str(e)
            if self.verbose:
                print(f"  ❌ Failed: {e}")
        
        return result
    
    def _find_fastqc_dir(self, fastq_path: str) -> Optional[str]:
        """Find FastQC output - accepts FastQC dirs/zips OR FASTQ files"""
        fastq_path = Path(fastq_path)
        
        # Case 1: Already a FastQC directory (e.g., sample_fastqc/)
        if fastq_path.is_dir() and fastq_path.name.endswith('_fastqc'):
            if (fastq_path / "fastqc_data.txt").exists():
                return str(fastq_path)
        
        # Case 2: FastQC zip file (e.g., sample_fastqc.zip)
        if fastq_path.is_file() and fastq_path.name.endswith('_fastqc.zip'):
            return str(fastq_path)
        
        # Case 3: FASTQ file - look for corresponding FastQC directory
        base_name = fastq_path.stem
        
        # Remove common FASTQ extensions
        for ext in ['.fastq', '.fq', '.gz']:
            if base_name.endswith(ext):
                base_name = base_name[:-len(ext)]
        
        # Look for FastQC directory
        possible_dirs = [
            fastq_path.parent / f"{base_name}_fastqc",
            fastq_path.parent / "fastqc_output" / f"{base_name}_fastqc",
            Path(f"{base_name}_fastqc")
        ]
        
        for dir_path in possible_dirs:
            if dir_path.exists() and dir_path.is_dir():
                return str(dir_path)
        
        return None
    
    def _calculate_statistics(self) -> BatchStatistics:
        """Calculate aggregate statistics from all results"""
        successful_results = [r for r in self.results if r.status == "success"]
        
        if not successful_results:
            return BatchStatistics(pass_count=0, warn_count=0, fail_count=0)
        
        # Count PASS/WARN/FAIL (case-insensitive)
        pass_count = sum(1 for r in successful_results if r.overall_status and r.overall_status.upper() == "PASS")
        warn_count = sum(1 for r in successful_results if r.overall_status and r.overall_status.upper() == "WARN")
        fail_count = sum(1 for r in successful_results if r.overall_status and r.overall_status.upper() == "FAIL")
        
        # Collect metrics from analysis JSONs
        gc_values = []
        quality_values = []
        duplication_values = []
        
        for result in successful_results:
            if not result.analysis_json:
                continue
            
            try:
                with open(result.analysis_json, 'r') as f:
                    analysis_data = json.load(f)
                
                # Extract GC content (supports both old and new format)
                gc_val = None
                if 'metrics' in analysis_data:
                    # New format: gc_content with summary
                    if 'gc_content' in analysis_data['metrics']:
                        gc_metric = analysis_data['metrics']['gc_content']
                        summary = gc_metric.get('summary', '')
                        # Parse "Normal GC content: 49.7% (expected ~52.0%)"
                        import re
                        match = re.search(r'(\d+\.?\d*)%', summary)
                        if match:
                            gc_val = float(match.group(1))
                    # Old format: GC Content with details
                    elif 'GC Content' in analysis_data['metrics']:
                        gc_metric = analysis_data['metrics']['GC Content']
                        if 'details' in gc_metric and 'actual_gc' in gc_metric['details']:
                            gc_val = gc_metric['details']['actual_gc']
                
                if gc_val is not None:
                    gc_values.append(gc_val)
                
                # Extract quality score
                quality_val = None
                if 'metrics' in analysis_data:
                    # New format: per_base_quality with summary
                    if 'per_base_quality' in analysis_data['metrics']:
                        quality_metric = analysis_data['metrics']['per_base_quality']
                        summary = quality_metric.get('summary', '')
                        # Parse "Excellent quality: mean Q=39.4, median Q=40.3"
                        import re
                        match = re.search(r'mean Q=([\d.]+)', summary)
                        if match:
                            quality_val = float(match.group(1))
                    # Old format: Per Base Sequence Quality with details
                    elif 'Per Base Sequence Quality' in analysis_data['metrics']:
                        quality_metric = analysis_data['metrics']['Per Base Sequence Quality']
                        if 'details' in quality_metric and 'mean_quality' in quality_metric['details']:
                            quality_val = quality_metric['details']['mean_quality']
                
                if quality_val is not None:
                    quality_values.append(quality_val)
                
                # Extract duplication rate
                dup_val = None
                if 'metrics' in analysis_data:
                    # New format: duplication_levels with summary
                    if 'duplication_levels' in analysis_data['metrics']:
                        dup_metric = analysis_data['metrics']['duplication_levels']
                        summary = dup_metric.get('summary', '')
                        # Parse "High duplication: 86.1% (acceptable for RNA-seq/ChIP-seq)"
                        import re
                        match = re.search(r'(\d+\.?\d*)%', summary)
                        if match:
                            dup_val = float(match.group(1))
                    # Old format: Sequence Duplication Levels with details
                    elif 'Sequence Duplication Levels' in analysis_data['metrics']:
                        dup_metric = analysis_data['metrics']['Sequence Duplication Levels']
                        if 'details' in dup_metric and 'percent_duplicates' in dup_metric['details']:
                            dup_val = dup_metric['details']['percent_duplicates']
                
                if dup_val is not None:
                    duplication_values.append(dup_val)
                
            except (json.JSONDecodeError, FileNotFoundError, KeyError):
                continue
        
        # Calculate statistics
        import statistics
        
        gc_mean = statistics.mean(gc_values) if gc_values else None
        gc_std = statistics.stdev(gc_values) if len(gc_values) > 1 else None
        
        quality_mean = statistics.mean(quality_values) if quality_values else None
        quality_std = statistics.stdev(quality_values) if len(quality_values) > 1 else None
        
        duplication_mean = statistics.mean(duplication_values) if duplication_values else None
        duplication_std = statistics.stdev(duplication_values) if len(duplication_values) > 1 else None
        
        # Detect outliers (values > 2 standard deviations from mean)
        outliers = []
        
        if gc_mean and gc_std and gc_std > 0:
            for result in successful_results:
                if not result.analysis_json:
                    continue
                try:
                    with open(result.analysis_json, 'r') as f:
                        analysis_data = json.load(f)
                    
                    if 'metrics' in analysis_data and 'GC Content' in analysis_data['metrics']:
                        gc_metric = analysis_data['metrics']['GC Content']
                        if 'details' in gc_metric and 'actual_gc' in gc_metric['details']:
                            gc_val = gc_metric['details']['actual_gc']
                            deviation = abs(gc_val - gc_mean) / gc_std
                            if deviation > 2.0:
                                outliers.append({
                                    'sample': result.sample_name,
                                    'metric': 'GC content',
                                    'value': f"{gc_val:.1f}%",
                                    'reason': 'possible contamination' if gc_val < gc_mean else 'unusual GC distribution'
                                })
                except:
                    pass
        
        if quality_mean and quality_std and quality_std > 0:
            for result in successful_results:
                if not result.analysis_json:
                    continue
                try:
                    with open(result.analysis_json, 'r') as f:
                        analysis_data = json.load(f)
                    
                    if 'metrics' in analysis_data and 'Per Base Sequence Quality' in analysis_data['metrics']:
                        quality_metric = analysis_data['metrics']['Per Base Sequence Quality']
                        if 'details' in quality_metric and 'mean_quality' in quality_metric['details']:
                            q_val = quality_metric['details']['mean_quality']
                            deviation = abs(q_val - quality_mean) / quality_std
                            if deviation > 2.0 and q_val < 30:
                                outliers.append({
                                    'sample': result.sample_name,
                                    'metric': 'Quality',
                                    'value': f"Q{q_val:.1f}",
                                    'reason': 'resequence recommended' if q_val < 28 else 'below average quality'
                                })
                except:
                    pass
        
        return BatchStatistics(
            pass_count=pass_count,
            warn_count=warn_count,
            fail_count=fail_count,
            gc_mean=gc_mean,
            gc_std=gc_std,
            quality_mean=quality_mean,
            quality_std=quality_std,
            duplication_mean=duplication_mean,
            duplication_std=duplication_std,
            outliers=outliers
        )

    
    def process_all(self) -> BatchReport:

        start_time = datetime.now()
        
        if self.verbose:
            print(f"\n{'='*70}")
            print(f"🚀 Batch Processing: {len(self.sample_list)} samples")
            if self.organism:
                print(f"   Organism: {self.organism}")
            if self.experiment_type:
                print(f"   Experiment: {self.experiment_type}")
            print(f"   Parallel: {self.parallel} process(es)")
            if self.dry_run:
                print(f"   Mode: DRY-RUN (no fixes executed)")
            print(f"{'='*70}\n")
        
        if self.parallel > 1:
            # Parallel processing
            with ProcessPoolExecutor(max_workers=self.parallel) as executor:
                futures = {
                    executor.submit(self.process_sample, sample, idx+1): sample
                    for idx, sample in enumerate(self.sample_list)
                }
                
                for future in as_completed(futures):
                    result = future.result()
                    self.results.append(result)
        else:
            # Sequential processing
            for idx, sample in enumerate(self.sample_list):
                result = self.process_sample(sample, idx+1)
                self.results.append(result)
        
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        # Generate report
        successful = sum(1 for r in self.results if r.status == "success")
        failed = sum(1 for r in self.results if r.status == "failed")
        skipped = sum(1 for r in self.results if r.status == "skipped")
        
        # Calculate aggregate statistics
        statistics = self._calculate_statistics()
        
        report = BatchReport(
            total_samples=len(self.sample_list),
            successful=successful,
            failed=failed,
            skipped=skipped,
            start_time=start_time.isoformat(),
            end_time=end_time.isoformat(),
            duration_seconds=duration,
            organism=self.organism,
            experiment_type=self.experiment_type,
            statistics=statistics,
            sample_results=self.results
        )
        
        # Save batch report
        report_path = self.output_dir / "batch_report.json"
        with open(report_path, 'w') as f:
            json.dump(report.to_dict(), f, indent=2)
        
        if self.verbose:
            print(f"\n{'='*70}")
            print(f"✅ Batch Processing Complete!")
            print(f"   Duration: {duration:.1f}s ({duration/len(self.sample_list):.1f}s per sample)")
            print(f"   Successful: {successful}/{len(self.sample_list)}")
            if failed > 0:
                print(f"   Failed: {failed}")
            if skipped > 0:
                print(f"   Skipped: {skipped}")
            print(f"   Report: {report_path}")
            print(f"{'='*70}\n")
        
        return report
