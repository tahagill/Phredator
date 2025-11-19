"""
Adaptive Quality Threshold Calibration (AQTC)

Uses statistical analysis to calibrate QC thresholds based on actual data
distribution rather than fixed universal cutoffs.
"""

import statistics
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import math


@dataclass
class ThresholdCalibration:
    """Results from adaptive threshold calibration"""
    mean_quality: float
    std_dev: float
    median: float
    mad: float  # Median Absolute Deviation
    recommended_threshold: float
    confidence_level: str  # "high", "medium", "low"
    outlier_positions: List[int]
    trend: str  # "stable", "degrading", "improving"
    trend_rate: float  # Quality change per position


class AdaptiveThresholdCalibrator:
    """
    Calibrates QC thresholds based on data distribution.
    
    Algorithm:
    1. Analyze quality score distribution (mean, std, MAD)
    2. Detect outliers using MAD (robust to extreme values)
    3. Detect trend using simple linear regression
    4. Adjust thresholds based on overall data quality
    5. Estimate confidence using data consistency
    """
    
    def __init__(self):
        # Standard baseline thresholds
        self.baseline_thresholds = {
            "excellent": 35,
            "good": 30,
            "acceptable": 25,
            "marginal": 20
        }
    
    def calculate_mad(self, values: List[float]) -> float:
        """
        Calculate Median Absolute Deviation (MAD).
        
        MAD is more robust than standard deviation for outlier detection.
        Formula: MAD = median(|x_i - median(x)|)
        """
        if not values:
            return 0.0
        
        median = statistics.median(values)
        deviations = [abs(x - median) for x in values]
        mad = statistics.median(deviations)
        return mad
    
    def detect_outliers_mad(self, values: List[float], threshold: float = 3.0) -> List[int]:
        """
        Detect outlier positions using MAD method.
        
        A point is an outlier if |x_i - median| > threshold * MAD
        Using threshold=3 is standard in statistics (similar to 3-sigma rule)
        """
        if len(values) < 3:
            return []
        
        median = statistics.median(values)
        mad = self.calculate_mad(values)
        
        if mad == 0:  # All values identical
            return []
        
        outliers = []
        for i, val in enumerate(values):
            if abs(val - median) > threshold * mad:
                outliers.append(i)
        
        return outliers
    
    def calculate_trend(self, values: List[float]) -> Tuple[str, float]:
        """
        Detect quality trend using simple linear regression.
        
        Returns:
            trend: "stable", "degrading", or "improving"
            rate: slope (quality change per position)
        """
        if len(values) < 5:
            return "stable", 0.0
        
        n = len(values)
        x = list(range(n))
        y = values
        
        # Simple linear regression: y = mx + b
        x_mean = sum(x) / n
        y_mean = sum(y) / n
        
        numerator = sum((x[i] - x_mean) * (y[i] - y_mean) for i in range(n))
        denominator = sum((x[i] - x_mean) ** 2 for i in range(n))
        
        if denominator == 0:
            return "stable", 0.0
        
        slope = numerator / denominator
        
        # Classify trend based on slope magnitude
        if abs(slope) < 0.05:  # Less than 0.05 quality units per position
            trend = "stable"
        elif slope < 0:
            trend = "degrading"
        else:
            trend = "improving"
        
        return trend, slope
    
    def calibrate_from_per_base_quality(self, per_base_quality: Dict) -> ThresholdCalibration:
        """
        Main calibration algorithm.
        
        Steps:
        1. Extract mean quality at each position
        2. Calculate distribution statistics (mean, std, median, MAD)
        3. Detect outlier positions
        4. Detect quality trend
        5. Calibrate threshold based on overall quality level
        6. Estimate confidence based on data consistency
        """
        if not per_base_quality:
            return self._default_calibration()
        
        # Step 1: Extract quality values
        quality_values = []
        for pos_data in per_base_quality.values():
            if isinstance(pos_data, dict) and "mean" in pos_data:
                quality_values.append(pos_data["mean"])
        
        if not quality_values:
            return self._default_calibration()
        
        # Step 2: Distribution statistics
        mean_q = statistics.mean(quality_values)
        std_q = statistics.stdev(quality_values) if len(quality_values) > 1 else 0.0
        median_q = statistics.median(quality_values)
        mad = self.calculate_mad(quality_values)
        
        # Step 3: Outlier detection
        outliers = self.detect_outliers_mad(quality_values)
        
        # Step 4: Trend detection
        trend, trend_rate = self.calculate_trend(quality_values)
        
        # Step 5: Adaptive threshold calibration
        # If overall quality is high, use stricter thresholds
        # If overall quality is low, adjust expectations
        if mean_q >= 35:
            # High quality data - use strict thresholds
            recommended = self.baseline_thresholds["excellent"]
            adjustment = 0
        elif mean_q >= 30:
            # Good quality - standard thresholds
            recommended = self.baseline_thresholds["good"]
            adjustment = 0
        elif mean_q >= 25:
            # Acceptable quality - slightly relaxed
            recommended = self.baseline_thresholds["acceptable"]
            adjustment = 0
        else:
            # Low quality - adaptive thresholds
            # Don't fail everything just because baseline is low
            recommended = max(20, mean_q - std_q)
            adjustment = self.baseline_thresholds["acceptable"] - recommended
        
        # Step 6: Confidence estimation
        # High confidence: low std, few outliers, stable trend
        # Low confidence: high std, many outliers, degrading trend
        confidence = self._estimate_confidence(std_q, len(outliers), len(quality_values), trend)
        
        return ThresholdCalibration(
            mean_quality=mean_q,
            std_dev=std_q,
            median=median_q,
            mad=mad,
            recommended_threshold=recommended,
            confidence_level=confidence,
            outlier_positions=outliers,
            trend=trend,
            trend_rate=trend_rate
        )
    
    def _estimate_confidence(self, std: float, num_outliers: int, 
                           total_positions: int, trend: str) -> str:
        """
        Estimate confidence in the calibration.
        
        High confidence: consistent data with few outliers
        Low confidence: variable data with many outliers
        """
        # Calculate coefficient of variation (CV = std/mean ratio)
        # and outlier proportion
        outlier_ratio = num_outliers / total_positions if total_positions > 0 else 0
        
        # High confidence: low variability, few outliers, stable trend
        if std < 3.0 and outlier_ratio < 0.1 and trend == "stable":
            return "high"
        # Low confidence: high variability, many outliers, or degrading trend
        elif std > 8.0 or outlier_ratio > 0.3 or trend == "degrading":
            return "low"
        else:
            return "medium"
    
    def _default_calibration(self) -> ThresholdCalibration:
        """Return default calibration when no data available"""
        return ThresholdCalibration(
            mean_quality=0.0,
            std_dev=0.0,
            median=0.0,
            mad=0.0,
            recommended_threshold=self.baseline_thresholds["acceptable"],
            confidence_level="low",
            outlier_positions=[],
            trend="unknown",
            trend_rate=0.0
        )
    
    def interpret_calibration(self, calibration: ThresholdCalibration) -> str:
        """
        Generate human-readable interpretation of calibration results.
        """
        interpretation = []
        
        # Overall quality assessment
        if calibration.mean_quality >= 35:
            interpretation.append("High quality data (mean Q={:.1f})".format(calibration.mean_quality))
        elif calibration.mean_quality >= 30:
            interpretation.append("Good quality data (mean Q={:.1f})".format(calibration.mean_quality))
        elif calibration.mean_quality >= 25:
            interpretation.append("Acceptable quality data (mean Q={:.1f})".format(calibration.mean_quality))
        else:
            interpretation.append("Low quality data (mean Q={:.1f})".format(calibration.mean_quality))
        
        # Consistency assessment
        if calibration.std_dev < 3.0:
            interpretation.append("Very consistent quality across positions")
        elif calibration.std_dev < 6.0:
            interpretation.append("Moderate variation in quality")
        else:
            interpretation.append("High variation in quality (std={:.1f})".format(calibration.std_dev))
        
        # Trend information
        if calibration.trend == "degrading":
            interpretation.append("Quality degrades along read (rate={:.3f}/position)".format(abs(calibration.trend_rate)))
        elif calibration.trend == "improving":
            interpretation.append("Quality improves along read (unusual)")
        else:
            interpretation.append("Stable quality across read length")
        
        # Outlier information
        if len(calibration.outlier_positions) > 0:
            interpretation.append("{} outlier positions detected".format(len(calibration.outlier_positions)))
        
        # Confidence
        interpretation.append("Confidence: {}".format(calibration.confidence_level))
        
        return ". ".join(interpretation) + "."
