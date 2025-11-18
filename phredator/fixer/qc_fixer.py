

import json
import os
import re
from typing import Dict, List, Any, Optional
from dataclasses import dataclass, field, asdict

from phredator.utils.tool_checker import ToolChecker


@dataclass
class FixSuggestion:
    """Represents a single fix suggestion"""
    category: str  # e.g., "quality_trimming", "adapter_removal", "duplicate_removal"
    priority: str  # "high", "medium", "low"
    description: str
    command: str  # Actual command to run
    reason: str
    tool_required: Optional[str] = None  # Tool key required to run this fix


@dataclass
class QCFixResult:
    """Results from QC fixer"""
    sample_name: str
    input_file: str  # Placeholder for actual input file
    fixes_applied: List[Dict[str, Any]] = field(default_factory=list)
    suggested_pipeline: List[str] = field(default_factory=list)
    tool_availability: Optional[Dict[str, Any]] = None  # Tool availability info
    read_length: Optional[int] = None
    is_paired_end: Optional[bool] = None
    
    def to_json(self) -> str:

        return json.dumps(asdict(self), indent=4)
    
    def to_dict(self) -> Dict[str, Any]:

        return asdict(self)


class Fixer:
    
    def __init__(self, input_path: str, input_reads: str = None, check_tools: bool = False):
        self.input_path = input_path
        self.input_reads = input_reads or "INPUT_READS.fastq.gz"
        self.analysis_data = None
        self.sample_name = "Unknown"
        self.fixes = []
        self.check_tools = check_tools
        self.tool_checker = ToolChecker() if check_tools else None
        
        self.read_length = None
        self.is_paired_end = False
        self.quality_threshold = 20
        
    def load_analysis(self) -> Dict[str, Any]:
        if not os.path.exists(self.input_path):
            raise FileNotFoundError(f"Analysis file not found: {self.input_path}")
        
        with open(self.input_path, 'r') as f:
            self.analysis_data = json.load(f)
        
        if 'sample_name' in self.analysis_data:
            self.sample_name = self.analysis_data['sample_name']
        
        self._detect_parameters(self.analysis_data)
        
        return self.analysis_data
    
    def _detect_parameters(self, data=None):
        if data is None:
            data = self.analysis_data
            
        if not data:
            return
            
        basic_stats = data.get('basic_statistics', data.get('parsed_data', {}).get('Basic Statistics', []))
        
        if isinstance(basic_stats, list):
            stats_dict = {item[0]: item[1] for item in basic_stats if len(item) >= 2}
        else:
            stats_dict = basic_stats
        
        if 'Sequence length' in stats_dict:
            seq_len = stats_dict['Sequence length']
            if '-' in str(seq_len):
                parts = str(seq_len).split('-')
                try:
                    self.read_length = int(parts[1])  # Use max value from range
                except:
                    pass
            else:
                try:
                    self.read_length = int(seq_len)
                except:
                    pass
        
        filename = stats_dict.get('Filename', self.input_reads)
        if filename:
            self.is_paired_end = bool(re.search(r'[_\.]R?[12][_\.]|[_\.]forward|[_\.]reverse', filename, re.I))
        
        profile_info = data.get('profile_info', '')
        if 'COVID' in profile_info or 'SARS-CoV-2' in profile_info:
            self.quality_threshold = 30
        elif 'metagenomics' in profile_info.lower():
            self.quality_threshold = 15
        else:
            self.quality_threshold = 20
    
    def _calculate_minlen(self) -> int:
        if self.read_length is None:
            return 36
        
        read_len = self.read_length
        if read_len < 75:
            return int(read_len * 0.5)
        elif read_len <= 150:
            return 36
        else:
            return int(read_len * 0.4)
    
    def _sort_and_select_fixes(self, all_fixes):
        priority_order = {"high": 0, "medium": 1, "low": 2}
        category_order = {
            "quality_trimming": 0,
            "adapter_removal": 1,
            "contamination_screening": 2,
            "contamination_removal": 3,
            "duplicate_removal": 4
        }
        
        sorted_fixes = sorted(
            all_fixes,
            key=lambda x: (priority_order.get(x.priority, 3), category_order.get(x.category, 5))
        )
        
        seen_categories = set()
        pipeline = []
        
        for fix in sorted_fixes:
            if fix.category not in seen_categories:
                pipeline.append(f"# {fix.description}")
                pipeline.append(fix.command)
                pipeline.append("")  # Blank line for readability
                seen_categories.add(fix.category)
        
        # Add final QC check
        pipeline.append("# Re-run FastQC to verify improvements")
        pipeline.append(f"fastqc {self.sample_name}_trimmed.fastq.gz -o fastqc_output/")
        
        return pipeline
    
    def generate_quality_trim_fixes(self, metrics: Dict) -> List[FixSuggestion]:
        fixes = []
        quality_metric = metrics.get("per_base_quality", {})
        
        if quality_metric.get("status") in ["warn", "fail"]:
            minlen = self._calculate_minlen()
            quality_threshold = self.quality_threshold if self.quality_threshold else 20
            input_file = self.input_reads if self.input_reads else "INPUT_READS.fastq.gz"
            
            if self.check_tools and self.tool_checker:
                available_tools = self.tool_checker.suggest_alternatives("quality_trim")
                if not available_tools or not isinstance(available_tools, list) or len(available_tools) == 0:
                    available_tools = ["fastp", "trimmomatic", "cutadapt"]
            else:
                available_tools = ["fastp", "trimmomatic", "cutadapt"]
            
            if "fastp" in available_tools:
                if self.is_paired_end:
                    r1_file = input_file
                    r2_file = input_file.replace("_R1", "_R2").replace("_1.", "_2.")
                    cmd = f"fastp -i {r1_file} -I {r2_file} -o {self.sample_name}_R1_trimmed.fastq.gz -O {self.sample_name}_R2_trimmed.fastq.gz -q {quality_threshold} -l {minlen}"
                else:
                    cmd = f"fastp -i {input_file} -o {self.sample_name}_trimmed.fastq.gz -q {quality_threshold} -l {minlen}"
                
                fixes.append(FixSuggestion(
                    category="quality_trimming",
                    priority="high",
                    description="Quality trim using fastp",
                    command=cmd,
                    reason=f"Low quality bases detected (threshold Q{quality_threshold})",
                    tool_required="fastp"
                ))
            
            if "trimmomatic" in available_tools:
                if self.is_paired_end:
                    r1_file = input_file
                    r2_file = input_file.replace("_R1", "_R2").replace("_1.", "_2.")
                    cmd = f"trimmomatic PE -phred33 {r1_file} {r2_file} {self.sample_name}_R1_paired.fastq.gz {self.sample_name}_R1_unpaired.fastq.gz {self.sample_name}_R2_paired.fastq.gz {self.sample_name}_R2_unpaired.fastq.gz LEADING:{quality_threshold} TRAILING:{quality_threshold} SLIDINGWINDOW:4:{quality_threshold} MINLEN:{minlen}"
                else:
                    cmd = f"trimmomatic SE -phred33 {input_file} {self.sample_name}_trimmed.fastq.gz LEADING:{quality_threshold} TRAILING:{quality_threshold} SLIDINGWINDOW:4:{quality_threshold} MINLEN:{minlen}"
                
                fixes.append(FixSuggestion(
                    category="quality_trimming",
                    priority="high",
                    description="Quality trim using Trimmomatic",
                    command=cmd,
                    reason=f"Low quality bases detected (threshold Q{quality_threshold})",
                    tool_required="trimmomatic"
                ))
        
        return fixes
    
    def generate_adapter_trim_fixes(self, metrics: Dict) -> List[FixSuggestion]:
        fixes = []
        adapter_metric = metrics.get("adapter_content", {})
        
        if adapter_metric.get("status") in ["warn", "fail"]:
            minlen = self._calculate_minlen()
            quality_threshold = self.quality_threshold if self.quality_threshold else 20
            input_file = self.input_reads if self.input_reads else "INPUT_READS.fastq.gz"
            
            if self.check_tools and self.tool_checker:
                available_tools = self.tool_checker.suggest_alternatives("adapter_removal")
                if not available_tools or not isinstance(available_tools, list) or len(available_tools) == 0:
                    available_tools = ["cutadapt", "fastp", "trimmomatic"]
            else:
                available_tools = ["cutadapt", "fastp", "trimmomatic"]
            
            if "cutadapt" in available_tools:
                if self.is_paired_end:
                    r1_file = input_file
                    r2_file = input_file.replace("_R1", "_R2").replace("_1.", "_2.")
                    cmd = f"cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -q {quality_threshold} --minimum-length {minlen} -o {self.sample_name}_R1_trimmed.fastq.gz -p {self.sample_name}_R2_trimmed.fastq.gz {r1_file} {r2_file}"
                else:
                    cmd = f"cutadapt -a AGATCGGAAGAG -q {quality_threshold} --minimum-length {minlen} -o {self.sample_name}_trimmed.fastq.gz {input_file}"
                
                fixes.append(FixSuggestion(
                    category="adapter_removal",
                    priority="medium",
                    description="Remove Illumina adapters using Cutadapt",
                    command=cmd,
                    reason="Adapter contamination detected in reads",
                    tool_required="cutadapt"
                ))
            
            if "fastp" in available_tools:
                if self.is_paired_end:
                    r1_file = input_file
                    r2_file = input_file.replace("_R1", "_R2").replace("_1.", "_2.")
                    cmd = f"fastp -i {r1_file} -I {r2_file} -o {self.sample_name}_R1_trimmed.fastq.gz -O {self.sample_name}_R2_trimmed.fastq.gz --detect_adapter_for_pe -q {quality_threshold} -l {minlen}"
                else:
                    cmd = f"fastp -i {input_file} -o {self.sample_name}_trimmed.fastq.gz --detect_adapter_for_pe -q {quality_threshold} -l {minlen}"
                
                fixes.append(FixSuggestion(
                    category="adapter_removal",
                    priority="medium",
                    description="Auto-detect and remove adapters using fastp",
                    command=cmd,
                    reason="Adapter contamination detected in reads",
                    tool_required="fastp"
                ))
            
            if "trimmomatic" in available_tools:
                if self.is_paired_end:
                    r1_file = input_file
                    r2_file = input_file.replace("_R1", "_R2").replace("_1.", "_2.")
                    cmd = f"trimmomatic PE -phred33 {r1_file} {r2_file} {self.sample_name}_R1_paired.fastq.gz {self.sample_name}_R1_unpaired.fastq.gz {self.sample_name}_R2_paired.fastq.gz {self.sample_name}_R2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 MINLEN:{minlen}"
                else:
                    cmd = f"trimmomatic SE -phred33 {input_file} {self.sample_name}_trimmed.fastq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 MINLEN:{minlen}"
                
                fixes.append(FixSuggestion(
                    category="adapter_removal",
                    priority="medium",
                    description="Remove adapters using Trimmomatic with adapter file",
                    command=cmd,
                    reason="Adapter contamination detected in reads",
                    tool_required="trimmomatic"
                ))
        
        return fixes
    
    def generate_deduplication_fixes(self, metrics: Dict) -> List[FixSuggestion]:
        fixes = []
        dup_metric = metrics.get("sequence_duplication", {})
        
        if dup_metric.get("status") in ["warn", "fail"]:
            profile_info = ""
            if self.analysis_data:
                profile_info = self.analysis_data.get("profile_info", "")
            is_rnaseq = "rna" in profile_info.lower() or "rnaseq" in str(self.input_reads).lower()
            
            if is_rnaseq:
                fixes.append(FixSuggestion(
                    category="duplicate_removal",
                    priority="low",
                    description="Note: High duplication is normal for RNA-seq",
                    command="# RNA-seq samples naturally have high duplication from highly expressed genes",
                    reason="RNA-seq duplication: Do not remove duplicates - they represent biological signal from gene expression",
                    tool_required=None
                ))
            else:
                if self.check_tools and self.tool_checker:
                    available_tools = self.tool_checker.suggest_alternatives("deduplication")
                    if not available_tools or not isinstance(available_tools, list) or len(available_tools) == 0:
                        available_tools = ["picard", "samtools"]
                else:
                    available_tools = ["picard", "samtools"]
                
                if "picard" in available_tools:
                    fixes.append(FixSuggestion(
                        category="duplicate_removal",
                        priority="medium",
                        description="Remove PCR duplicates using Picard",
                        command=f"picard MarkDuplicates I=aligned.bam O={self.sample_name}_dedup.bam M=metrics.txt REMOVE_DUPLICATES=true",
                        reason=f"High duplication level detected: {dup_metric.get('duplication_level', 'unknown')}%",
                        tool_required="picard"
                    ))
        
        return fixes
    
    def generate_contamination_fixes(self, metrics: Dict) -> List[FixSuggestion]:
        fixes = []
        return fixes
    
    def generate_pipeline_suggestion(self, all_fixes: List[FixSuggestion]) -> List[str]:
        return self._sort_and_select_fixes(all_fixes)
    
    def generate_fixes(self) -> QCFixResult:

        if self.analysis_data is None:
            self.load_analysis()
        
        metrics = self.analysis_data.get("metrics", {})
        
        all_fixes = []
        
        # Generate fixes for each category
        all_fixes.extend(self.generate_quality_trim_fixes(metrics))
        all_fixes.extend(self.generate_adapter_trim_fixes(metrics))
        all_fixes.extend(self.generate_deduplication_fixes(metrics))
        all_fixes.extend(self.generate_contamination_fixes(metrics))
        
        # Generate suggested pipeline
        pipeline = self.generate_pipeline_suggestion(all_fixes)
        
        # Convert fixes to dict format
        fixes_dict = [asdict(fix) for fix in all_fixes]
        
        # Prepare tool availability info
        tool_availability = None
        if self.check_tools and self.tool_checker:
            installed = self.tool_checker.get_installed_tools()
            all_tools = list(self.tool_checker.tools.keys())
            tool_availability = {
                'installed': installed,
                'missing': [t for t in all_tools if t not in installed],
                'total_installed': len(installed),
                'total_missing': len(all_tools) - len(installed)
            }
        
        result = QCFixResult(
            sample_name=self.sample_name,
            input_file=self.input_reads,
            fixes_applied=fixes_dict,
            suggested_pipeline=pipeline,
            tool_availability=tool_availability,
            read_length=self.read_length,
            is_paired_end=self.is_paired_end
        )
        
        return result
    
    def run(self) -> QCFixResult:

        self.load_analysis()
        return self.generate_fixes()
