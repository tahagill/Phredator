"""
MultiQC parser for aggregate QC analysis.
"""
import json
from typing import Dict, List
from pathlib import Path


class MultiQCParser:
    """Parse MultiQC JSON output for aggregate analysis."""
    
    def __init__(self, multiqc_json_path: str):
        """
        Initialize MultiQC parser.
        
        Args:
            multiqc_json_path: Path to multiqc_data.json file
        """
        self.multiqc_json_path = multiqc_json_path
        self.data = None
    
    def parse(self) -> Dict:
        """Parse MultiQC JSON and extract FastQC data."""
        with open(self.multiqc_json_path, 'r') as f:
            self.data = json.load(f)
        
        # Extract FastQC-specific data
        fastqc_data = {}
        
        # MultiQC stores data in report_general_stats_data
        if 'report_general_stats_data' in self.data:
            for sample_data in self.data['report_general_stats_data']:
                for sample_name, metrics in sample_data.items():
                    if sample_name not in fastqc_data:
                        fastqc_data[sample_name] = {}
                    
                    # Extract relevant metrics
                    if 'percent_gc' in metrics:
                        fastqc_data[sample_name]['gc_content'] = metrics['percent_gc']
                    if 'percent_duplicates' in metrics:
                        fastqc_data[sample_name]['duplication'] = metrics['percent_duplicates']
                    if 'avg_sequence_length' in metrics:
                        fastqc_data[sample_name]['sequence_length'] = metrics['avg_sequence_length']
                    if 'total_sequences' in metrics:
                        fastqc_data[sample_name]['total_sequences'] = metrics['total_sequences']
        
        return {
            'samples': fastqc_data,
            'total_samples': len(fastqc_data),
            'multiqc_version': self.data.get('config_version', 'unknown')
        }
    
    def get_summary_statistics(self) -> Dict:
        """Calculate summary statistics from MultiQC data."""
        if not self.data:
            self.parse()
        
        parsed = self.parse()
        samples = parsed['samples']
        
        if not samples:
            return {}
        
        # Calculate aggregate statistics
        gc_values = [s.get('gc_content', 0) for s in samples.values() if 'gc_content' in s]
        dup_values = [s.get('duplication', 0) for s in samples.values() if 'duplication' in s]
        
        return {
            'total_samples': len(samples),
            'gc_content': {
                'mean': sum(gc_values) / len(gc_values) if gc_values else 0,
                'min': min(gc_values) if gc_values else 0,
                'max': max(gc_values) if gc_values else 0
            },
            'duplication': {
                'mean': sum(dup_values) / len(dup_values) if dup_values else 0,
                'samples_high_dup': sum(1 for d in dup_values if d > 50)
            },
            'sample_names': list(samples.keys())
        }
