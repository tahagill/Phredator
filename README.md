# PHREDATOR

**Phred-based Quality Assessment Tool for Sequencing Data**

A professional, rule-based QC toolkit for NGS (Next-Generation Sequencing) data analysis. Phredator provides intelligent quality control analysis with context-aware recommendations for FastQC reports.

[![PyPI version](https://badge.fury.io/py/phredator.svg)](https://badge.fury.io/py/phredator)
[![PyPI downloads](https://img.shields.io/pypi/dm/phredator.svg)](https://pypi.org/project/phredator/)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![GitHub stars](https://img.shields.io/github/stars/tahagill/Phredator.svg?style=social&label=Star)](https://github.com/tahagill/Phredator)
[![GitHub issues](https://img.shields.io/github/issues/tahagill/Phredator.svg)](https://github.com/tahagill/Phredator/issues)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

---

## 🌟 Features

- **🔍 Intelligent QC Analysis**: Context-aware quality assessment using organism and experiment-specific profiles
- **⚡ Batch Processing**: Process hundreds of samples with parallel execution
- **🧬 25+ Organism Profiles**: Pre-configured thresholds for human, mouse, E. coli, yeast, and more
- **🧪 Experiment-Specific Rules**: Different expectations for WGS, WES, RNA-seq, ChIP-seq, and metagenomics
- **📊 Aggregate Statistics**: Mean ± SD for GC%, quality, duplication with outlier detection
- **🎯 Smart Recommendations**: Context-aware advice (e.g., "DON'T remove duplicates for RNA-seq")
- **🔤 Fuzzy Matching**: Natural input like "Human", "RNA seq", "chip-seq" automatically matched
- **🎨 Professional Output**: Subread-style terminal display with color-coded status indicators
- **📦 MultiQC Support**: Parse and analyze MultiQC aggregate reports
- **🧪 Automated Fixes**: Generate trimming/filtering commands with smart parameter detection

---

## 📋 Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
  - [Single File Analysis](#single-file-analysis)
  - [Batch Processing](#batch-processing)
  - [MultiQC Integration](#multiqc-integration)
- [Command Reference](#command-reference)
- [Parameters](#parameters)
- [Output Format](#output-format)
- [Organism Profiles](#organism-profiles)
- [Experiment Types](#experiment-types)
- [Examples](#examples)
- [Development](#development)
- [Citation](#citation)
- [License](#license)

---

## 🚀 Installation

### From PyPI (Recommended)

```bash
pip install phredator
```

### From Source

```bash
git clone https://github.com/tahagill/Phredator.git
cd phredator
pip install -e .
```

### Requirements

- Python 3.8+
- NumPy
- PyYAML

---

## ⚡ Quick Start

```bash
# Single file analysis
phredator parse sample_fastqc.zip --output parsed.json
phredator analyze parsed.json --organism human --experiment-type rnaseq --output analysis.json

# Batch processing (multiple samples)
phredator batch *_fastqc.zip --organism human --experiment-type rnaseq --parallel 4

# List available organisms
phredator list-organisms

# Show usage examples
phredator examples
```

---

## 📖 Usage

### Single File Analysis

**Step 1: Parse FastQC output**
```bash
phredator parse sample1_fastqc.zip --output sample1.json
```

**Step 2: Analyze with organism/experiment profiles**
```bash
phredator analyze sample1.json \
    --organism human \
    --experiment-type rnaseq \
    --output analysis1.json
```

**Step 3: Generate fix suggestions (optional)**
```bash
phredator fix analysis1.json \
    --input-reads sample1.fastq.gz \
    --output fixes1.json
```

**Example Output:**
```
  ██████╗ ██╗  ██╗██████╗ ███████╗██████╗  ██████╗ ████████╗ ██████╗ ██████╗ 
  ██╔══██╗██║  ██║██╔══██╗██╔════╝██╔══██╗██╔═══██╗╚══██╔══╝██╔═══██╗██╔══██╗
  ██████╔╝███████║██████╔╝█████╗  ██║  ██║███████║   ██║   ██║   ██║██████╔╝
  ██╔═══╝ ██╔══██║██╔══██╗██╔══╝  ██║  ██║██╔══██║   ██║   ██║   ██║██╔══██╗
  ██║     ██║  ██║██║  ██║███████╗██████╔╝██║  ██║   ██║   ╚██████╔╝██║  ██║
  ╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝╚═════╝ ╚═╝  ╚═╝   ╚═╝    ╚═════╝ ╚═╝  ╚═╝

//================================\\
||       QC Analysis Report       ||
\\================================//

         Sample : SRR29487748
         Status : WARN
        Profile : Organism: Human (Homo sapiens) | Experiment: RNA Sequencing

//================================\\
||    Quality Assessments         ||
\\================================//

   [PASS] Per Base Quality               : Excellent quality: mean Q=39.4, median Q=40.3
   [PASS] Gc Content                     : Normal GC content: 49.7% (expected ~52.0%)
   [WARN] Duplication Levels             : High duplication: 86.1% (acceptable for RNA-seq)
   [PASS] Adapter Content                : Minimal adapter content: 0.55%
   [PASS] Overrepresented Sequences      : No overrepresented sequences detected

//================================\\
||       Recommendations          ||
\\================================//

   1. High duplication is normal for RNA-seq (abundant transcripts)
   2. DO NOT remove duplicates for RNA-seq - they are real biological signal
```

### Batch Processing

Process multiple samples in one command with aggregate statistics:

**Method 1: Direct file list**
```bash
phredator batch sample1.zip sample2.zip sample3.zip \
    --organism human \
    --experiment-type rnaseq \
    --output-dir batch_results/ \
    --parallel 4
```

**Method 2: Using wildcards**
```bash
phredator batch results/qc/*_fastqc.zip \
    --organism human \
    --experiment-type chipseq \
    --output-dir chipseq_qc/ \
    --parallel 8
```

**Method 3: From a list file**
```bash
# Create sample list
ls /path/to/*_fastqc.zip > samples.txt

# Process all
phredator batch samples.txt \
    --organism mouse \
    --experiment-type wgs \
    --parallel 4
```

**Batch Output:**
```
//================================\\
||   Batch Processing Summary     ||
\\================================//

   Processed 8/8 samples:
   [PASS] 5 (62%)
   [WARN] 3 (38%)
   [FAIL] 0 (0%)

//================================\\
||      Sample Details            ||
\\================================//

   [PASS] Sample1_fastqc
      [PASS] Per Base Quality   : Excellent quality: mean Q=39.4
      [PASS] Gc Content          : Normal GC content: 49.7%
      [PASS] Duplication Levels  : Duplication: 81.8% (normal for RNA-seq)
      ...

//================================\\
||      Average Metrics           ||
\\================================//

   GC content   : 49.6% ± 2.6%
   Quality      : Q39.3 ± 0.1
   Duplication  : 81.3% ± 11.1% (normal for RNA-seq)

//================================\\
||      Recommendations           ||
\\================================//

   1. DO NOT remove duplicates for RNA-seq
   2. High duplication is normal for RNA-seq (abundant transcripts)
   3. Monitor for potential contamination
```

### MultiQC Integration

Parse MultiQC aggregate reports:

```bash
phredator parse multiqc_data/multiqc_data.json --output multiqc_parsed.json
```

---

## 📚 Command Reference

### `phredator parse`
Parse FastQC output or MultiQC JSON.

```bash
phredator parse <input> [OPTIONS]
```

**Arguments:**
- `input`: Path to FastQC zip/folder or MultiQC JSON file

**Options:**
- `--output FILE`: Output JSON file (default: `parsed.json`)
- `--verbose`: Enable verbose logging

### `phredator analyze`
Analyze parsed QC data with organism/experiment profiles.

```bash
phredator analyze <input_json> [OPTIONS]
```

**Arguments:**
- `input_json`: Path to parsed JSON from `parse` command

**Options:**
- `--organism NAME`: Organism profile (see `list-organisms`)
- `--experiment-type TYPE`: Experiment type (`wgs`, `wes`, `rnaseq`, `chipseq`, `metagenomics`)
- `--expected-gc FLOAT`: Override expected GC% (default: from profile)
- `--output FILE`: Output JSON file (default: `analysis.json`)
- `--verbose`: Enable verbose logging

### `phredator batch`
Process multiple samples with parallel execution.

```bash
phredator batch <samples...> [OPTIONS]
```

**Arguments:**
- `samples`: FastQC files/zips OR path to sample list file

**Options:**
- `--organism NAME`: Organism profile for all samples
- `--experiment-type TYPE`: Experiment type for all samples
- `--output-dir DIR`: Output directory (default: `batch_output`)
- `--parallel N`: Number of parallel processes (default: 1)
- `--verbose`: Enable verbose logging

### `phredator list-organisms`
List available organism profiles.

```bash
phredator list-organisms [OPTIONS]
```

**Options:**
- `--detailed`: Show full profile details (GC%, assembly, etc.)

### `phredator examples`
Show comprehensive usage examples.

```bash
phredator examples
```

---

## 🔧 Parameters

### Organism Profiles

25+ pre-configured organism profiles with expected GC content and quality thresholds:

| Organism | Code | GC% | Notes |
|----------|------|-----|-------|
| Human | `human` | 41% | Homo sapiens |
| Mouse | `mouse` | 42% | Mus musculus |
| Rat | `rat` | 42% | Rattus norvegicus |
| Zebrafish | `zebrafish` | 36% | Danio rerio |
| Drosophila | `drosophila` | 42% | D. melanogaster |
| C. elegans | `celegans` | 36% | Worm model |
| Yeast | `yeast` | 38% | S. cerevisiae |
| E. coli | `ecoli` | 51% | Bacteria |
| Arabidopsis | `arabidopsis` | 36% | Plant model |
| COVID-19 | `covid19` | 38% | SARS-CoV-2 |

**View all organisms:**
```bash
phredator list-organisms --detailed
```

### Experiment Types

Different QC expectations for each experiment type:

| Type | Code | Key Characteristics |
|------|------|---------------------|
| Whole Genome Sequencing | `wgs` | Low duplication, uniform coverage |
| Whole Exome Sequencing | `wes` | Moderate duplication, targeted regions |
| RNA Sequencing | `rnaseq` | High duplication OK (abundant transcripts) |
| ChIP Sequencing | `chipseq` | Variable duplication, enriched regions |
| Metagenomics | `metagenomics` | Variable GC, mixed organisms |

**Fuzzy matching supported:**
- "RNA seq" → `rnaseq`
- "chip-seq" → `chipseq`
- "WGS" → `wgs`
- "Human" → `human`

---

## 📊 Output Format

### JSON Output

All commands produce structured JSON output:

```json
{
  "sample_name": "SRR29487748",
  "overall_status": "WARN",
  "profile_info": "Organism: Human | Experiment: RNA-seq",
  "metrics": {
    "per_base_quality": {
      "status": "PASS",
      "summary": "Excellent quality: mean Q=39.4, median Q=40.3",
      "details": {...}
    },
    "gc_content": {
      "status": "PASS",
      "summary": "Normal GC content: 49.7% (expected ~52.0%)",
      "details": {...}
    }
  },
  "all_recommendations": [
    "High duplication is normal for RNA-seq",
    "DO NOT remove duplicates for RNA-seq"
  ]
}
```

### Log Files

Automatic `.log` files created alongside JSON outputs with timestamp and full analysis details.

---

## 💡 Examples

### Example 1: RNA-seq Quality Control

```bash
# Parse 8 RNA-seq samples
phredator batch ~/rnaseq/results/*_fastqc.zip \
    --organism human \
    --experiment-type rnaseq \
    --output-dir rnaseq_qc \
    --parallel 4

# Output shows:
# - 5 PASS, 3 WARN, 0 FAIL
# - Average GC: 49.6% ± 2.6%
# - Average Quality: Q39.3 ± 0.1
# - Duplication: 81.3% ± 11.1% (normal for RNA-seq)
# - Recommendations: DON'T remove duplicates
```

### Example 2: ChIP-seq with Adapter Detection

```bash
phredator batch ~/chipseq/qc/*_fastqc.zip \
    --organism human \
    --experiment-type chipseq \
    --output-dir chipseq_qc \
    --parallel 8

# Output includes:
# - Adapter content warnings
# - Recommendations: "Trim adapters using Cutadapt"
# - Fix commands generated automatically
```

### Example 3: COVID-19 Variant Surveillance

```bash
phredator batch sarscov2_samples/*_fastqc.zip \
    --organism covid19 \
    --experiment-type metagenomics \
    --parallel 16

# High quality thresholds (Q30+) enforced
# GC% ~38% expected for SARS-CoV-2
```

### Example 4: MultiQC Integration

```bash
# Parse MultiQC aggregate report
phredator parse multiqc_data/multiqc_data.json --output multiqc_summary.json

# Extract metrics for all 100+ samples at once
```

---

## 🧪 Development

### Running Tests

```bash
# Run comprehensive integration test
pytest tests/test_integration.py -v

# Run with sample data
pytest tests/ -v
```

### Project Structure

```
phredator/
├── phredator/
│   ├── __init__.py
│   ├── __main__.py
│   ├── cli/
│   │   ├── cli.py              # Command-line interface
│   │   └── commands.py
│   ├── parser/
│   │   ├── fastqc_parser.py    # Parse FastQC output
│   │   ├── multiqc_parser.py   # Parse MultiQC JSON
│   │   └── batch_processor.py  # Batch processing engine
│   ├── analyzer/
│   │   └── qc_analyzer.py      # QC analysis engine
│   ├── rules/
│   │   └── qc_rules.py         # Rule-based evaluation
│   ├── fixer/
│   │   └── qc_fixer.py         # Fix suggestion generator
│   ├── config/
│   │   ├── organisms/          # 25+ organism profiles
│   │   └── experiment_types/   # Experiment-specific rules
│   └── utils/
│       └── helpers.py
├── tests/
│   ├── test_integration.py     # Comprehensive integration test
│   └── test_data/              # Sample FastQC data
├── README.md
├── setup.py
├── requirements.txt
└── LICENSE
```

---

## 📖 Citation

If you use Phredator in your research, please cite:

```bibtex
@software{phredator2025,
  title={Phredator: Intelligent Quality Control for NGS Data},
  author={Ahmad, Taha},
  year={2025},
  publisher={Zenodo},
  doi={10.5281/zenodo.XXXXXXX},
  url={https://github.com/tahagill/Phredator}
}
```

**Note**: After creating a release on Zenodo, replace `XXXXXXX` with your actual DOI number in both the badge above and this citation.

---

## 📝 License

MIT License - see [LICENSE](LICENSE) file for details.

---

## 🤝 Contributing

Contributions welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

---

## 🐛 Issues

Found a bug? Have a feature request? Please open an issue on [GitHub](https://github.com/tahagill/Phredator/issues).

---

## 📧 Contact

- **Author**: Taha Ahmad
- **Email**: tahagill99@gmail.com
- **GitHub**: [@tahagill](https://github.com/tahagill)

---

## 🙏 Acknowledgments

- FastQC team for the excellent QC tool
- MultiQC for aggregate reporting
- Subread for terminal output inspiration
