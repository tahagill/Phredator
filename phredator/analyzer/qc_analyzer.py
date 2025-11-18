

import json
import os
from typing import Dict, List, Any, Optional
from dataclasses import dataclass, field, asdict

from phredator.rules.qc_rules import QCRulesEngine, QCStatus
from phredator.utils.profile_loader import ProfileLoader


@dataclass
class QCAnalysisResult:
    """Results from QC analysis"""
    sample_name: str
    overall_status: str
    overall_summary: str
    metrics: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    all_recommendations: List[str] = field(default_factory=list)
    organism: Optional[str] = None
    experiment_type: Optional[str] = None
    profile_info: Optional[str] = None
    
    def to_json(self) -> str:

        return json.dumps(asdict(self), indent=4)
    
    def to_dict(self) -> Dict[str, Any]:

        return asdict(self)


class Analyzer:
    
    def __init__(self, input_path: str, expected_gc: float = 50.0, organism: Optional[str] = None, experiment_type: Optional[str] = None):
        self.input_path = input_path
        self.expected_gc = expected_gc
        self.organism = organism
        self.experiment_type = experiment_type
        self.parsed_data = None
        self.sample_name = "Unknown"
        self.profile_loader = ProfileLoader()
        self.thresholds = self.profile_loader.get_combined_thresholds(organism, experiment_type)
        
        # Initialize rules engine with custom thresholds
        self.rules_engine = QCRulesEngine(thresholds=self.thresholds)
        if organism and 'gc_content' in self.thresholds:
            self.expected_gc = self.thresholds['gc_content'].get('mean', expected_gc)
        
    def load_data(self) -> Dict[str, Any]:

        if not os.path.exists(self.input_path):
            raise FileNotFoundError(f"Input file not found: {self.input_path}")
        
        with open(self.input_path, 'r') as f:
            data = json.load(f)
        
        # Validate required fields
        required_fields = ["sample_name"]
        missing_fields = [f for f in required_fields if f not in data]
        if missing_fields:
            raise ValueError(f"Missing required fields in input data: {missing_fields}")
        
        self.sample_name = data.get("sample_name", "Unknown")
        self.parsed_data = data
        
        return data
    
    def analyze(self) -> QCAnalysisResult:

        if self.parsed_data is None:
            self.load_data()
        
        if self.parsed_data is None:
            raise RuntimeError("Failed to load data")
        
        # Dictionary to store individual metric results
        individual_results = {}
        all_recommendations = []
        
        # Analyze per-base quality
        per_base_quality = self.parsed_data.get("per_base_quality", {})
        if per_base_quality:
            status, summary, recommendations = self.rules_engine.evaluate_per_base_quality(per_base_quality)
            individual_results["per_base_quality"] = (status, summary, recommendations)
            all_recommendations.extend(recommendations)
        
        # Analyze GC content (use mean from distribution)
        gc_content_mean = self.parsed_data.get("gc_content_mean", 0)
        if gc_content_mean > 0:
            status, summary, recommendations = self.rules_engine.evaluate_gc_content(gc_content_mean, self.expected_gc)
            individual_results["gc_content"] = (status, summary, recommendations)
            all_recommendations.extend(recommendations)
        
        # Analyze duplication levels
        duplication_levels = self.parsed_data.get("duplication_levels", {})
        total_dedup_pct = self.parsed_data.get("total_deduplicated_percentage", 100.0)
        if duplication_levels:
            # Calculate actual duplication percentage from levels
            # total_deduplicated_percentage tells us what % of unique library remains
            # So duplication % = 100 - total_deduplicated_percentage
            duplication_percentage = 100.0 - total_dedup_pct
            status, summary, recommendations = self.rules_engine.evaluate_duplication(duplication_levels)
            individual_results["duplication_levels"] = (status, summary, recommendations)
            all_recommendations.extend(recommendations)
        
        # Analyze adapter content
        adapter_content = self.parsed_data.get("adapter_content", {})
        status, summary, recommendations = self.rules_engine.evaluate_adapter_content(adapter_content)
        individual_results["adapter_content"] = (status, summary, recommendations)
        all_recommendations.extend(recommendations)
        
        # Analyze overrepresented sequences
        overrepresented_sequences = self.parsed_data.get("overrepresented_sequences", [])
        status, summary, recommendations = self.rules_engine.evaluate_overrepresented_sequences(overrepresented_sequences)
        individual_results["overrepresented_sequences"] = (status, summary, recommendations)
        all_recommendations.extend(recommendations)
        
        # Generate overall assessment
        overall_status, overall_summary = self.rules_engine.generate_overall_assessment(individual_results)
        
        # Format metrics for output
        metrics = {}
        for metric_name, (status, summary, recommendations) in individual_results.items():
            metrics[metric_name] = {
                "status": status.value,
                "summary": summary,
                "recommendations": recommendations
            }
        profile_info = []
        if self.organism:
            profile_info.append(f"Organism: {self.thresholds.get('organism_name', self.organism)}")
        if self.experiment_type:
            profile_info.append(f"Experiment: {self.thresholds.get('experiment_name', self.experiment_type)}")
        
        result = QCAnalysisResult(
            sample_name=self.sample_name,
            overall_status=overall_status.value,
            overall_summary=overall_summary,
            metrics=metrics,
            all_recommendations=list(dict.fromkeys(all_recommendations)),  # Remove duplicates, preserve order
            organism=self.organism,
            experiment_type=self.experiment_type,
            profile_info=" | ".join(profile_info) if profile_info else None
        )
        
        return result
    
    def run(self) -> QCAnalysisResult:

        self.load_data()
        return self.analyze()
