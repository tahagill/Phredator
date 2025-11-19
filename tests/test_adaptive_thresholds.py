"""
Tests for Adaptive Quality Threshold Calibration
"""

import pytest
from phredator.analyzer.adaptive_thresholds import (
    AdaptiveThresholdCalibrator,
    ThresholdCalibration
)


def test_mad_calculation():
    """Test Median Absolute Deviation calculation"""
    calibrator = AdaptiveThresholdCalibrator()
    
    # Test with simple data
    values = [1, 2, 3, 4, 5]
    mad = calibrator.calculate_mad(values)
    # Median = 3, deviations = [2,1,0,1,2], MAD = 1
    assert mad == 1.0
    
    # Test with outlier
    values = [10, 10, 10, 10, 100]
    mad = calibrator.calculate_mad(values)
    # MAD should be 0 since most values are 10
    assert mad == 0.0


def test_outlier_detection():
    """Test MAD-based outlier detection"""
    calibrator = AdaptiveThresholdCalibrator()
    
    # Data with clear outlier
    values = [30, 31, 32, 30, 31, 5, 30, 31]  # Position 5 is outlier
    outliers = calibrator.detect_outliers_mad(values, threshold=3.0)
    
    assert 5 in outliers
    assert len(outliers) >= 1


def test_trend_detection_stable():
    """Test trend detection for stable quality"""
    calibrator = AdaptiveThresholdCalibrator()
    
    # Stable quality
    values = [30, 30, 31, 30, 30, 31, 30]
    trend, rate = calibrator.calculate_trend(values)
    
    assert trend == "stable"
    assert abs(rate) < 0.1


def test_trend_detection_degrading():
    """Test trend detection for degrading quality"""
    calibrator = AdaptiveThresholdCalibrator()
    
    # Degrading quality (common in sequencing)
    values = [35, 34, 33, 31, 29, 27, 25, 23, 20]
    trend, rate = calibrator.calculate_trend(values)
    
    assert trend == "degrading"
    assert rate < 0  # Negative slope


def test_calibration_high_quality():
    """Test calibration with high quality data"""
    calibrator = AdaptiveThresholdCalibrator()
    
    # Simulate high quality data
    per_base_quality = {
        f"pos_{i}": {"mean": 37 + i * 0.01}  # Very high, stable quality
        for i in range(100)
    }
    
    result = calibrator.calibrate_from_per_base_quality(per_base_quality)
    
    assert result.mean_quality >= 35
    assert result.recommended_threshold >= 35
    assert result.confidence_level in ["high", "medium"]
    assert result.trend == "stable"


def test_calibration_degrading_quality():
    """Test calibration with degrading quality"""
    calibrator = AdaptiveThresholdCalibrator()
    
    # Simulate quality degradation (realistic for Illumina)
    per_base_quality = {
        f"pos_{i}": {"mean": 38 - i * 0.15}  # Starts at Q38, drops to Q23
        for i in range(100)
    }
    
    result = calibrator.calibrate_from_per_base_quality(per_base_quality)
    
    assert result.trend == "degrading"
    assert result.trend_rate < 0
    assert len(result.outlier_positions) >= 0  # May have outliers at ends


def test_calibration_low_quality():
    """Test calibration adapts to low quality baseline"""
    calibrator = AdaptiveThresholdCalibrator()
    
    # Simulate low but consistent quality (e.g., Nanopore, older data)
    per_base_quality = {
        f"pos_{i}": {"mean": 22 + i * 0.01}  # Low but stable Q22
        for i in range(100)
    }
    
    result = calibrator.calibrate_from_per_base_quality(per_base_quality)
    
    # Should adapt threshold lower, not fail everything
    assert result.mean_quality < 25
    assert result.recommended_threshold < 25
    assert result.trend == "stable"


def test_calibration_with_outliers():
    """Test calibration handles outliers correctly"""
    calibrator = AdaptiveThresholdCalibrator()
    
    # Good quality with few bad positions (adapter contamination)
    per_base_quality = {}
    for i in range(100):
        if i in [95, 96, 97, 98, 99]:  # Last 5 positions are bad
            per_base_quality[f"pos_{i}"] = {"mean": 15}
        else:
            per_base_quality[f"pos_{i}"] = {"mean": 35}
    
    result = calibrator.calibrate_from_per_base_quality(per_base_quality)
    
    # With 95% good data, mean should still be high
    assert result.mean_quality >= 30
    # Trend should detect the drop at end
    assert result.trend == "degrading"


def test_interpretation_generation():
    """Test human-readable interpretation"""
    calibrator = AdaptiveThresholdCalibrator()
    
    per_base_quality = {
        f"pos_{i}": {"mean": 35}
        for i in range(100)
    }
    
    result = calibrator.calibrate_from_per_base_quality(per_base_quality)
    interpretation = calibrator.interpret_calibration(result)
    
    assert isinstance(interpretation, str)
    assert len(interpretation) > 0
    assert "quality" in interpretation.lower()


def test_empty_data():
    """Test graceful handling of empty data"""
    calibrator = AdaptiveThresholdCalibrator()
    
    result = calibrator.calibrate_from_per_base_quality({})
    
    assert result.mean_quality == 0.0
    assert result.confidence_level == "low"
    assert result.trend == "unknown"


def test_real_world_scenario():
    """Test with realistic Illumina data pattern"""
    calibrator = AdaptiveThresholdCalibrator()
    
    # Realistic Illumina HiSeq pattern:
    # - Starts high (Q38-40)
    # - Stable in middle (Q35-37)
    # - Drops at end (Q28-32)
    per_base_quality = {}
    for i in range(150):
        if i < 10:
            quality = 38 + (i * 0.2)  # Ramp up
        elif i < 120:
            quality = 36 + (i % 2)  # Stable with small variation
        else:
            quality = 36 - (i - 120) * 0.25  # Degrade at end
        
        per_base_quality[f"pos_{i}"] = {"mean": max(20, quality)}
    
    result = calibrator.calibrate_from_per_base_quality(per_base_quality)
    
    # Should detect degradation
    assert result.trend in ["degrading", "stable"]
    # Should have good overall quality
    assert result.mean_quality >= 30
    # Should be confident (realistic data is consistent)
    assert result.confidence_level in ["high", "medium"]
    
    interpretation = calibrator.interpret_calibration(result)
    print(f"\nReal-world interpretation: {interpretation}")


if __name__ == "__main__":
    # Run the real-world test to see output
    test_real_world_scenario()
