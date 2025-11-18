

from typing import Dict, List, Tuple
from dataclasses import dataclass
from enum import Enum


class QCStatus(Enum):
    """QC status levels"""
    PASS = "pass"
    WARN = "warn"
    FAIL = "fail"


@dataclass
class QCRule:
    """Represents a single QC rule with thresholds"""
    name: str
    description: str
    pass_threshold: float
    warn_threshold: float
    higher_is_better: bool = True  # If False, lower values are better


# Industry-standard QC thresholds based on FastQC and best practices
QC_RULES = {
    # Per-base quality score thresholds (Phred scores)
    "per_base_quality_mean": QCRule(
        name="Per Base Quality (Mean)",
        description="Average quality score across all bases should be high",
        pass_threshold=28.0,  # Q28 = 99.84% accuracy
        warn_threshold=20.0,  # Q20 = 99% accuracy
        higher_is_better=True
    ),
    
    "per_base_quality_median": QCRule(
        name="Per Base Quality (Median)",
        description="Median quality score should be consistently high",
        pass_threshold=30.0,  # Q30 = 99.9% accuracy
        warn_threshold=25.0,
        higher_is_better=True
    ),
    
    # GC content thresholds (organism-dependent, but these are general)
    "gc_content": QCRule(
        name="GC Content",
        description="GC content should be within expected range for organism",
        pass_threshold=65.0,  # Upper bound
        warn_threshold=70.0,
        higher_is_better=False  # Deviation from expected is bad
    ),
    
    "gc_content_lower": QCRule(
        name="GC Content (Lower Bound)",
        description="GC content should not be too low",
        pass_threshold=35.0,  # Lower bound
        warn_threshold=30.0,
        higher_is_better=True
    ),
    
    # Sequence duplication levels (percentage)
    "duplication_level": QCRule(
        name="Sequence Duplication Level",
        description="High duplication may indicate PCR over-amplification or low complexity",
        pass_threshold=20.0,  # < 20% duplication is good
        warn_threshold=50.0,  # > 50% is concerning
        higher_is_better=False
    ),
    
    # Adapter content (percentage)
    "adapter_content": QCRule(
        name="Adapter Content",
        description="Adapter sequences should be minimal or absent",
        pass_threshold=5.0,   # < 5% is acceptable
        warn_threshold=10.0,  # > 10% needs trimming
        higher_is_better=False
    ),
    
    # Overrepresented sequences
    "overrepresented_sequences": QCRule(
        name="Overrepresented Sequences",
        description="Number of overrepresented sequences (contamination indicator)",
        pass_threshold=5.0,   # < 5 sequences is fine
        warn_threshold=10.0,  # > 10 may indicate contamination
        higher_is_better=False
    )
}


class QCRulesEngine:
    
    def __init__(self, custom_rules: Dict[str, QCRule] = None, thresholds: Dict[str, any] = None):
        self.rules = QC_RULES.copy()
        if custom_rules:
            self.rules.update(custom_rules)
        
        self.thresholds = thresholds or {}
        
        if thresholds:
            self._update_rules_from_thresholds(thresholds)
    
    def _update_rules_from_thresholds(self, thresholds: Dict):

        # Update duplication rule
        if 'duplication' in thresholds:
            dup = thresholds['duplication']
            self.rules["duplication_level"] = QCRule(
                name="Sequence Duplication Level",
                description="High duplication may indicate PCR over-amplification or low complexity",
                pass_threshold=dup.get('acceptable', 20.0),
                warn_threshold=dup.get('critical', 50.0),
                higher_is_better=False
            )
        
        # Update adapter rule
        if 'adapters' in thresholds:
            adp = thresholds['adapters']
            self.rules["adapter_content"] = QCRule(
                name="Adapter Content",
                description="Adapter sequences should be minimal or absent",
                pass_threshold=adp.get('acceptable', 5.0),
                warn_threshold=adp.get('critical', 15.0),
                higher_is_better=False
            )
        
        # Update quality rules
        if 'quality' in thresholds:
            qual = thresholds['quality']
            mean_min = qual.get('mean_quality_min', 28)
            self.rules["per_base_quality_mean"] = QCRule(
                name="Per Base Quality (Mean)",
                description="Average quality score across all bases should be high",
                pass_threshold=mean_min,
                warn_threshold=max(20.0, mean_min - 8),
                higher_is_better=True
            )
    
    def evaluate_per_base_quality(self, per_base_quality: Dict[str, Dict[str, float]]) -> Tuple[QCStatus, str, List[str]]:
        pass  # docstring removed
        if not per_base_quality:
            return QCStatus.FAIL, "No per-base quality data available", ["Check if FastQC data is complete"]
        
        mean_rule = self.rules["per_base_quality_mean"]
        median_rule = self.rules["per_base_quality_median"]
        
        # Calculate overall statistics
        all_means = [v.get("mean", 0) for v in per_base_quality.values()]
        all_medians = [v.get("median", 0) for v in per_base_quality.values()]
        
        if not all_means or not all_medians:
            return QCStatus.FAIL, "Invalid quality data", ["Verify FastQC output format"]
        
        avg_mean = sum(all_means) / len(all_means)
        avg_median = sum(all_medians) / len(all_medians)
        quality_drop = False
        if len(all_means) > 10:
            last_10_pct = all_means[-max(1, len(all_means)//10):]
            if sum(last_10_pct) / len(last_10_pct) < mean_rule.warn_threshold:
                quality_drop = True
        
        recommendations = []
        
        # Determine status
        if avg_mean >= mean_rule.pass_threshold and avg_median >= median_rule.pass_threshold:
            status = QCStatus.PASS
            summary = f"Excellent quality: mean Q={avg_mean:.1f}, median Q={avg_median:.1f}"
        elif avg_mean >= mean_rule.warn_threshold and avg_median >= median_rule.warn_threshold:
            status = QCStatus.WARN
            summary = f"Acceptable quality: mean Q={avg_mean:.1f}, median Q={avg_median:.1f}"
            recommendations.append("Consider quality filtering or trimming low-quality bases")
        else:
            status = QCStatus.FAIL
            summary = f"Poor quality: mean Q={avg_mean:.1f}, median Q={avg_median:.1f}"
            recommendations.append("Quality trimming strongly recommended")
            recommendations.append("Consider discarding this sample or re-sequencing")
        
        if quality_drop:
            recommendations.append("Quality drops at read ends - trim last 5-10 bases")
        
        return status, summary, recommendations
    
    def evaluate_gc_content(self, gc_content: float, expected_gc: float = 50.0) -> Tuple[QCStatus, str, List[str]]:
        pass  # docstring removed
        if gc_content <= 0 or gc_content > 100:
            return QCStatus.FAIL, f"Invalid GC content: {gc_content}%", ["Check FastQC data integrity"]
        if 'gc_content' in self.thresholds:
            gc_profile = self.thresholds['gc_content']
            gc_range = gc_profile.get('range', [35, 65])
            expected_gc = gc_profile.get('mean', expected_gc)
            tolerance = gc_profile.get('tolerance', 5.0)
            
            lower_bound = gc_range[0]
            upper_bound = gc_range[1]
        else:
            # Fallback to rules
            upper_rule = self.rules["gc_content"]
            lower_rule = self.rules["gc_content_lower"]
            lower_bound = lower_rule.pass_threshold
            upper_bound = upper_rule.pass_threshold
            tolerance = 10.0
        
        deviation = abs(gc_content - expected_gc)
        recommendations = []
        if lower_bound <= gc_content <= upper_bound and deviation < tolerance:
            status = QCStatus.PASS
            summary = f"Normal GC content: {gc_content:.1f}% (expected ~{expected_gc:.1f}%)"
        elif deviation < tolerance * 2:
            status = QCStatus.WARN
            summary = f"GC content {gc_content:.1f}% deviates {deviation:.1f}% from expected {expected_gc:.1f}%"
            recommendations.append("GC content outside typical range")
            if deviation > tolerance * 1.5:
                recommendations.append("May indicate contamination or adapter content")
        else:
            status = QCStatus.FAIL
            summary = f"GC content {gc_content:.1f}% deviates {deviation:.1f}% from expected {expected_gc:.1f}%"
            recommendations.append("Severe GC bias detected")
            recommendations.append("Check for contamination, adapter dimers, or wrong organism")
        
        return status, summary, recommendations
    
    def evaluate_duplication(self, duplication_levels: Dict[str, float]) -> Tuple[QCStatus, str, List[str]]:
        pass  # docstring removed
        if not duplication_levels:
            return QCStatus.WARN, "No duplication data available", []
        
        rule = self.rules["duplication_level"]
        
        # Calculate percentage of library that is duplicated (duplication level > 1)
        # duplication_levels contains: "1" -> % of unique seqs, "2" -> % with 2 copies, etc.
        total_duplication = sum(
            pct for level, pct in duplication_levels.items()
            if level != "1"  # Exclude level 1 (unique sequences)
        )
        
        recommendations = []
        check_duplicates = True
        allow_high_dup = False
        if 'special' in self.thresholds:
            check_duplicates = self.thresholds['special'].get('check_duplicates', True)
            allow_high_dup = self.thresholds['special'].get('allow_high_duplication', False)
        
        # For RNA-seq or experiments that allow high duplication
        if allow_high_dup or not check_duplicates:
            if total_duplication <= rule.warn_threshold:
                status = QCStatus.PASS
                summary = f"Duplication: {total_duplication:.1f}% (normal for this experiment type)"
            else:
                status = QCStatus.WARN
                summary = f"High duplication: {total_duplication:.1f}% (acceptable for RNA-seq/ChIP-seq)"
                recommendations.append("High duplication is normal for RNA-seq (abundant transcripts)")
                recommendations.append("DO NOT remove duplicates for RNA-seq - they are real biological signal")
        else:
            # Standard DNA-seq evaluation
            if total_duplication <= rule.pass_threshold:
                status = QCStatus.PASS
                summary = f"Low duplication: {total_duplication:.1f}%"
            elif total_duplication <= rule.warn_threshold:
                status = QCStatus.WARN
                summary = f"Moderate duplication: {total_duplication:.1f}%"
                recommendations.append("Duplication may indicate PCR over-amplification")
                recommendations.append("Consider using duplicate removal tools (e.g., Picard MarkDuplicates)")
            else:
                status = QCStatus.FAIL
                summary = f"High duplication: {total_duplication:.1f}%"
                recommendations.append("Severe duplication detected")
                recommendations.append("Strong recommendation: remove duplicates before downstream analysis")
                recommendations.append("May indicate low library complexity or PCR issues")
        
        return status, summary, recommendations
    
    def evaluate_adapter_content(self, adapter_content: Dict[str, float]) -> Tuple[QCStatus, str, List[str]]:
        pass  # docstring removed
        if not adapter_content:
            return QCStatus.PASS, "No adapter content detected", []
        
        rule = self.rules["adapter_content"]
        
        # Find maximum adapter content
        max_adapter = max(adapter_content.values()) if adapter_content else 0
        max_adapter_pos = max(adapter_content, key=adapter_content.get) if adapter_content else "Unknown"
        
        recommendations = []
        
        if max_adapter <= rule.pass_threshold:
            status = QCStatus.PASS
            summary = f"Minimal adapter content: {max_adapter:.2f}%"
        elif max_adapter <= rule.warn_threshold:
            status = QCStatus.WARN
            summary = f"Adapter content detected: {max_adapter:.2f}% at position {max_adapter_pos}"
            recommendations.append(f"Trim adapters using tools like Cutadapt or Trimmomatic")
            recommendations.append(f"Adapters start appearing at position {max_adapter_pos}")
        else:
            status = QCStatus.FAIL
            summary = f"High adapter content: {max_adapter:.2f}% at position {max_adapter_pos}"
            recommendations.append("Adapter trimming is essential before analysis")
            recommendations.append(f"Use Cutadapt or Trimmomatic to remove adapters")
            recommendations.append("High adapter content may reduce mappability")
        
        return status, summary, recommendations
    
    def evaluate_overrepresented_sequences(self, overrepresented_sequences: List[str]) -> Tuple[QCStatus, str, List[str]]:
        pass  # docstring removed
        rule = self.rules["overrepresented_sequences"]
        
        num_overrep = len(overrepresented_sequences)
        recommendations = []
        
        if num_overrep == 0:
            status = QCStatus.PASS
            summary = "No overrepresented sequences detected"
        elif num_overrep <= rule.pass_threshold:
            status = QCStatus.PASS
            summary = f"Few overrepresented sequences: {num_overrep}"
            recommendations.append("Monitor for potential contamination")
        elif num_overrep <= rule.warn_threshold:
            status = QCStatus.WARN
            summary = f"Multiple overrepresented sequences: {num_overrep}"
            recommendations.append("Check for contamination or adapter dimers")
            recommendations.append("Run BLAST on overrepresented sequences to identify source")
        else:
            status = QCStatus.FAIL
            summary = f"Many overrepresented sequences: {num_overrep}"
            recommendations.append("Severe contamination likely")
            recommendations.append("BLAST overrepresented sequences to identify contaminants")
            recommendations.append("Consider re-library prep or sample cleanup")
        
        return status, summary, recommendations
    
    def generate_overall_assessment(self, individual_results: Dict[str, Tuple[QCStatus, str, List[str]]]) -> Tuple[QCStatus, str]:
        pass  # docstring removed
        status_counts = {QCStatus.PASS: 0, QCStatus.WARN: 0, QCStatus.FAIL: 0}
        
        for status, _, _ in individual_results.values():
            status_counts[status] += 1
        
        total_metrics = len(individual_results)
        
        # Overall status logic - more lenient approach
        if status_counts[QCStatus.FAIL] > total_metrics / 2:
            # More than half the metrics failed
            overall_status = QCStatus.FAIL
            overall_summary = f"Sample FAILED QC: {status_counts[QCStatus.FAIL]}/{total_metrics} metrics failed"
        elif status_counts[QCStatus.FAIL] > 0 or status_counts[QCStatus.WARN] > 0:
            # Some failures or warnings
            overall_status = QCStatus.WARN
            fails = status_counts[QCStatus.FAIL]
            warns = status_counts[QCStatus.WARN]
            if fails > 0:
                overall_summary = f"Sample quality concerns: {fails} failed, {warns} warned"
            else:
                overall_summary = f"Sample passed with warnings: {warns}/{total_metrics} metrics need attention"
        else:
            # All passed
            overall_status = QCStatus.PASS
            overall_summary = f"Sample PASSED QC: all {total_metrics} metrics passed"
        
        return overall_status, overall_summary
