"""
Batch parser for multi-sample FastQC analysis.
"""
import os
import json
from typing import List, Dict
from pathlib import Path
from phredator.parser.fastqc_parser import FastQCParser, FastQCReport


class BatchParser:
    """Parse multiple FastQC outputs and generate aggregate statistics."""
    
    def __init__(self, input_paths: List[str]):
        """
        Initialize batch parser.
        
        Args:
            input_paths: List of paths to FastQC folders or zip files
        """
        self.input_paths = input_paths
        self.reports: List[FastQCReport] = []
    
    def parse_all(self) -> List[FastQCReport]:
        """Parse all FastQC outputs."""
        for path in self.input_paths:
            try:
                parser = FastQCParser(path)
                report = parser.parse()
                self.reports.append(report)
            except Exception as e:
                print(f"[WARNING] Failed to parse {path}: {e}")
                continue
        
        return self.reports
    
    def get_aggregate_statistics(self) -> Dict:
        """Calculate aggregate statistics across all samples."""
        if not self.reports:
            return {}
        
        # Aggregate statistics
        total_samples = len(self.reports)
        total_sequences = sum(r.total_sequences for r in self.reports)
        
        # Quality statistics
        mean_qualities = []
        for report in self.reports:
            if report.per_base_quality:
                sample_mean = sum(
                    pos_data['mean'] 
                    for pos_data in report.per_base_quality.values()
                ) / len(report.per_base_quality)
                mean_qualities.append(sample_mean)
        
        avg_quality_across_samples = (
            sum(mean_qualities) / len(mean_qualities) if mean_qualities else 0
        )
        
        # GC content statistics
        gc_contents = [r.gc_content_mean for r in self.reports if r.gc_content_mean > 0]
        avg_gc = sum(gc_contents) / len(gc_contents) if gc_contents else 0
        min_gc = min(gc_contents) if gc_contents else 0
        max_gc = max(gc_contents) if gc_contents else 0
        
        # Duplication statistics
        duplication_rates = [
            100 - r.total_deduplicated_percentage 
            for r in self.reports
        ]
        avg_duplication = sum(duplication_rates) / len(duplication_rates)
        
        # Adapter contamination
        samples_with_adapters = sum(
            1 for r in self.reports if r.adapter_content
        )
        
        return {
            "total_samples": total_samples,
            "total_sequences": total_sequences,
            "average_quality": avg_quality_across_samples,
            "gc_content": {
                "mean": avg_gc,
                "min": min_gc,
                "max": max_gc,
                "values": gc_contents
            },
            "duplication": {
                "mean": avg_duplication,
                "rates": duplication_rates
            },
            "adapter_contamination": {
                "samples_affected": samples_with_adapters,
                "percentage": (samples_with_adapters / total_samples * 100) if total_samples > 0 else 0
            },
            "sample_names": [r.sample_name for r in self.reports]
        }
    
    def save_batch_report(self, output_path: str):
        """Save batch analysis report."""
        batch_data = {
            "aggregate_statistics": self.get_aggregate_statistics(),
            "individual_samples": [
                {
                    "sample_name": r.sample_name,
                    "total_sequences": r.total_sequences,
                    "gc_content_mean": r.gc_content_mean,
                    "duplication_pct": 100 - r.total_deduplicated_percentage,
                    "has_adapters": len(r.adapter_content) > 0
                }
                for r in self.reports
            ]
        }
        
        with open(output_path, 'w') as f:
            json.dump(batch_data, f, indent=4)
