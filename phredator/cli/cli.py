# phredator/cli/cli.py
import argparse
import os
import sys
import json
from pathlib import Path

from phredator.parser.fastqc_parser import FastQCParser
from phredator.parser.batch_parser import BatchParser
from phredator.parser.batch_processor import BatchProcessor
from phredator.analyzer.qc_analyzer import Analyzer
from phredator.fixer.qc_fixer import Fixer
from phredator.reporter.report_generator import Reporter
from phredator.pipeline.pipeline_runner import PipelineRunner
from phredator.utils.profile_loader import ProfileLoader

def main():
    parser = argparse.ArgumentParser(
        prog="phredator",
        description="Phredator: Rule-based QC toolkit for sequencing data"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # ---------------- Examples ----------------
    examples_parser = subparsers.add_parser("examples", help="Show usage examples for single file and batch processing")

    # ---------------- List Organisms ----------------
    list_orgs_parser = subparsers.add_parser("list-organisms", help="List all available organism profiles")
    list_orgs_parser.add_argument("--detailed", action="store_true", help="Show detailed GC content and thresholds")

    # ---------------- Parse ----------------
    parse_parser = subparsers.add_parser("parse", help="Parse FastQC/MultiQC data")
    parse_parser.add_argument("input_path", type=str, help="Path to FastQC zip/folder or MultiQC JSON (multiqc_data.json)")
    parse_parser.add_argument("--output", type=str, default="parsed.json", help="Output file")
    parse_parser.add_argument("--verbose", action="store_true", help="Verbose logging")

    # ---------------- Analyze ----------------
    analyze_parser = subparsers.add_parser("analyze", help="Analyze QC metrics using rule-based assessment")
    analyze_parser.add_argument("input_path", type=str, help="Path to parsed data (JSON)")
    analyze_parser.add_argument("--output", type=str, default="analysis.json", help="Output file")
    analyze_parser.add_argument("--expected-gc", type=float, default=None, help="Expected GC content for organism (overrides profile)")
    analyze_parser.add_argument("--organism", type=str,
                                help="Organism profile (use 'phredator list-organisms' to see all 25 available)")
    analyze_parser.add_argument("--experiment-type", type=str,
                                help="Experiment type (wgs, wes, rnaseq, chipseq, metagenomics) - fuzzy matching supported!")
    analyze_parser.add_argument("--verbose", action="store_true", help="Verbose logging")

    # ---------------- Fix ----------------
    fix_parser = subparsers.add_parser("fix", help="Generate automated fix suggestions")
    fix_parser.add_argument("input_path", type=str, help="Path to analysis data (JSON)")
    fix_parser.add_argument("--output", type=str, default="fixes.json", help="Output file")
    fix_parser.add_argument("--input-reads", type=str, help="Path to raw reads file for concrete commands")
    fix_parser.add_argument("--check-tools", action="store_true", help="Check tool availability and filter suggestions")
    fix_parser.add_argument("--show-tool-status", action="store_true", help="Display tool availability status")
    fix_parser.add_argument("--verbose", action="store_true", help="Verbose logging")

    # ---------------- Report ----------------
    report_parser = subparsers.add_parser("report", help="Generate comprehensive QC report")
    report_parser.add_argument("input_path", type=str, help="Path to data to report on (parsed/analysis/fixes JSON)")
    report_parser.add_argument("--output", type=str, default="report.json", help="Output file")
    report_parser.add_argument("--format", choices=["json", "csv", "summary"], default="json", help="Output format")
    report_parser.add_argument("--verbose", action="store_true", help="Verbose logging")

    # ---------------- Batch ----------------
    batch_parser = subparsers.add_parser("batch", help="Process multiple samples in batch")
    batch_parser.add_argument("samples", type=str, nargs='+', help="FastQC files/zips OR path to sample list file (one per line)")
    batch_parser.add_argument("--organism", type=str, help="Organism profile for all samples")
    batch_parser.add_argument("--experiment-type", type=str,
                             help="Experiment type (wgs, wes, rnaseq, chipseq, metagenomics) - fuzzy matching supported!")
    batch_parser.add_argument("--output-dir", type=str, default="batch_output", help="Output directory")
    batch_parser.add_argument("--parallel", type=int, default=1, help="Number of parallel processes (default: 1)")
    batch_parser.add_argument("--check-tools", action="store_true", default=True, help="Check tool availability")
    batch_parser.add_argument("--dry-run", action="store_true", help="Show what would be done without executing")
    batch_parser.add_argument("--verbose", action="store_true", help="Verbose logging")

    # ---------------- Pipeline ----------------
    pipeline_parser = subparsers.add_parser("pipeline", help="Run complete QC workflow with verification")
    pipeline_parser.add_argument("input_fastq", type=str, help="Path to input FASTQ file (must have FastQC results)")
    pipeline_parser.add_argument("--output-dir", type=str, default="pipeline_output", help="Output directory")
    pipeline_parser.add_argument("--organism", type=str,
                                 help="Organism profile (use 'phredator list-organisms' to see all 25 available)")
    pipeline_parser.add_argument("--experiment-type", type=str,
                                 help="Experiment type (wgs, wes, rnaseq, chipseq, metagenomics) - fuzzy matching supported!")
    pipeline_parser.add_argument("--check-tools", action="store_true", default=True, help="Check tool availability (default: True)")
    pipeline_parser.add_argument("--dry-run", action="store_true", help="Show what would be done without executing")
    pipeline_parser.add_argument("--verbose", action="store_true", help="Verbose logging")

    args = parser.parse_args()

    # Friendly fuzzy matching for organism and experiment-type
    fuzzy_loader = ProfileLoader()
    
    if hasattr(args, 'organism') and args.organism:
        original_input = args.organism
        available = fuzzy_loader.list_organisms()
        matched = fuzzy_loader._find_closest_match(args.organism, available)
        
        if matched and matched != original_input:
            print(f"\nINFO: Using organism profile '{matched}' (matched from '{original_input}')")
            args.organism = matched
        elif matched:
            args.organism = matched
    
    if hasattr(args, 'experiment_type') and args.experiment_type:
        original_input = args.experiment_type
        available = fuzzy_loader.list_experiment_types()
        matched = fuzzy_loader._find_closest_match(args.experiment_type, available)
        
        if matched and matched != original_input:
            print(f"INFO: Using experiment type '{matched}' (matched from '{original_input}')")
            args.experiment_type = matched
        elif matched:
            args.experiment_type = matched

    try:
        if args.command == "examples":
            print("""
//================================================================\\\\
||              PHREDATOR USAGE EXAMPLES                          ||
\\\\================================================================//

SINGLE FILE PROCESSING (parse → analyze → fix)
-----------------------------------------------
Process one sample at a time through the pipeline.
Good for: Exploratory analysis, debugging, testing

Example workflow:
  # Step 1: Parse FastQC output
  $ phredator parse sample1_fastqc.zip --output sample1.json

  # Step 2: Analyze with organism/experiment profiles
  $ phredator analyze sample1.json \\
      --organism human \\
      --experiment-type chipseq \\
      --output analysis1.json

  # Step 3: Generate fix suggestions
  $ phredator fix analysis1.json \\
      --input-reads sample1.fastq.gz \\
      --output fixes1.json

  Each command shows detailed output with [PASS]/[WARN]/[FAIL] status
  and saves a .log file alongside the .json file.


BATCH PROCESSING (all at once)
-------------------------------
Process multiple samples in one command.
Good for: Production runs, many samples, reproducibility

Method 1: From a list file
  # Create a sample list file
  $ cat > samples.txt << EOF
/path/to/sample1_fastqc.zip
/path/to/sample2_fastqc.zip
/path/to/sample3_fastqc.zip
EOF

  # Process all samples
  $ phredator batch samples.txt \\
      --organism human \\
      --experiment-type rnaseq \\
      --output-dir batch_results/ \\
      --parallel 4

Method 2: Multiple files directly
  # Pass files directly on command line
  $ phredator batch sample1_fastqc.zip sample2_fastqc.zip sample3_fastqc.zip \\
      --organism human \\
      --experiment-type rnaseq \\
      --output-dir batch_results/ \\
      --parallel 4

Method 3: Using wildcards (shell expansion)
  # Process all FastQC files in directory
  $ phredator batch results/qc/*_fastqc.zip \\
      --organism human \\
      --experiment-type rnaseq \\
      --output-dir batch_results/ \\
      --parallel 4

  Output:
    Processed 96 samples:
    [PASS] 78 (81%)
    [WARN] 12 (13%)
    [FAIL] 6 (6%)

    Average metrics:
    - GC content: 41.2% ± 2.3%
    - Quality: Q38.5 ± 1.2
    - Duplication: 65.3% ± 18.9% (normal for RNA-seq)

    Outliers:
    - Sample_23: GC content = 35.0% (possible contamination)
    - Sample_67: Quality = Q22.0 (resequence recommended)


FUZZY MATCHING
--------------
Phredator accepts natural input for organism and experiment types:

  "Human" → human
  "chip-seq" → chipseq
  "RNA seq" → rnaseq

See all available profiles:
  $ phredator list-organisms
  $ phredator list-organisms --detailed


REAL-WORLD EXAMPLE
------------------
You have 16 ChIP-seq samples in ~/MEF2C_ChIPseq/results/qc/:

Method 1: Quick batch with wildcards
  $ phredator batch ~/MEF2C_ChIPseq/results/qc/*_fastqc.zip \\
      --organism human \\
      --experiment-type chipseq \\
      --output-dir chipseq_qc/ \\
      --parallel 4

Method 2: Using a list file
  $ ls ~/MEF2C_ChIPseq/results/qc/*_fastqc.zip > chipseq_samples.txt
  $ phredator batch chipseq_samples.txt \\
      --organism human \\
      --experiment-type chipseq \\
      --output-dir chipseq_qc/ \\
      --parallel 4

Method 3: Individual sample inspection
  $ phredator parse ~/MEF2C_ChIPseq/results/qc/SRR35220282_1_fastqc.zip \\
      --output sample.json
  $ phredator analyze sample.json \\
      --organism human \\
      --experiment-type chipseq \\
      --output analysis.json

//================================================================\\\\
||  Need more help? Use --help with any command for details      ||
\\\\================================================================//
""")
        
        elif args.command == "list-organisms":
            profile_loader = ProfileLoader()
            organisms = profile_loader.list_organisms()
            
            print("\n🧬 Available Organism Profiles ({}):".format(len(organisms)))
            print("="*70)
            
            if args.detailed:
                for org in organisms:
                    profile = profile_loader.load_organism_profile(org)
                    if profile:
                        gc_mean = profile.gc_content.get('mean', 'N/A')
                        gc_range = profile.gc_content.get('range', ['N/A', 'N/A'])
                        print(f"\n  {org:15s} - {profile.name}")
                        print(f"    GC%:      {gc_mean}% (range: {gc_range[0]}-{gc_range[1]}%)")
                        print(f"    Assembly: {profile.assembly}")
                    else:
                        print(f"\n  {org:15s} - Could not load profile")
            else:
                # Simple 3-column list
                for i in range(0, len(organisms), 3):
                    row = organisms[i:i+3]
                    print("  " + "".join(f"{org:25s}" for org in row))
            
            print("\n" + "="*70)
            print(f"💡 Usage: phredator analyze ... --organism <name>")
            print(f"💡 Details: phredator list-organisms --detailed\n")
        
        elif args.command == "parse":
            if args.verbose:
                print(f"[INFO] Parsing data from: {args.input_path}")
            
            # Auto-detect MultiQC vs FastQC
            if args.input_path.endswith('.json') or 'multiqc_data.json' in args.input_path:
                # MultiQC JSON file
                from phredator.parser.multiqc_parser import MultiQCParser
                parser_obj = MultiQCParser(args.input_path)
                multiqc_data = parser_obj.parse()
                with open(args.output, "w") as f:
                    json.dump(multiqc_data, f, indent=2)
                if args.verbose:
                    print(f"[INFO] MultiQC parsing complete. {multiqc_data['total_samples']} samples found")
                    print(f"[INFO] Output saved to {args.output}")
            else:
                # FastQC file/zip
                parser_obj = FastQCParser(args.input_path)
                report = parser_obj.parse()
                with open(args.output, "w") as f:
                    f.write(report.to_json())
                if args.verbose:
                    print(f"[INFO] FastQC parsing complete. Output saved to {args.output}")

        elif args.command == "analyze":
            if args.verbose:
                print(f"[INFO] Analyzing parsed data: {args.input_path}")
                if args.organism:
                    print(f"[INFO] Using organism profile: {args.organism}")
                if args.experiment_type:
                    print(f"[INFO] Using experiment type: {args.experiment_type}")
            
            # Determine expected_gc
            expected_gc = args.expected_gc if args.expected_gc is not None else 50.0
            
            analyzer = Analyzer(
                args.input_path, 
                expected_gc=expected_gc,
                organism=args.organism,
                experiment_type=args.experiment_type
            )
            analysis = analyzer.run()
            
            # Save JSON output
            with open(args.output, "w") as f:
                f.write(analysis.to_json())
            
            # ANSI color codes
            RED = '\033[31m'  # Dark red instead of bright red
            BOLD = '\033[1m'
            RESET = '\033[0m'
            
            # Generate and display terminal summary (Subread-style)
            summary_lines = []
            summary_lines.append("")
            summary_lines.append(f"{RED}{BOLD}")
            summary_lines.append("  ██████╗ ██╗  ██╗██████╗ ███████╗██████╗  ██████╗ ████████╗ ██████╗ ██████╗ ")
            summary_lines.append("  ██╔══██╗██║  ██║██╔══██╗██╔════╝██╔══██╗██╔═══██╗╚══██╔══╝██╔═══██╗██╔══██╗")
            summary_lines.append("  ██████╔╝███████║██████╔╝█████╗  ██║  ██║███████║   ██║   ██║   ██║██████╔╝")
            summary_lines.append("  ██╔═══╝ ██╔══██║██╔══██╗██╔══╝  ██║  ██║██╔══██║   ██║   ██║   ██║██╔══██╗")
            summary_lines.append("  ██║     ██║  ██║██║  ██║███████╗██████╔╝██║  ██║   ██║   ╚██████╔╝██║  ██║")
            summary_lines.append("  ╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝╚═════╝ ╚═╝  ╚═╝   ╚═╝    ╚═════╝ ╚═╝  ╚═╝")
            summary_lines.append(f"{RESET}")
            summary_lines.append("")
            summary_lines.append("//================================\\\\")
            summary_lines.append("||       QC Analysis Report       ||")
            summary_lines.append("\\\\================================//")
            summary_lines.append("")
            summary_lines.append(f"         Sample : {analysis.sample_name}")
            summary_lines.append(f"         Status : {analysis.overall_status.upper()}")
            if analysis.profile_info:
                summary_lines.append(f"        Profile : {analysis.profile_info}")
            summary_lines.append("")
            summary_lines.append("//================================\\\\")
            summary_lines.append("||    Quality Assessments         ||")
            summary_lines.append("\\\\================================//")
            summary_lines.append("")
            
            for module, result in analysis.metrics.items():
                status = result['status'].upper()
                status_symbol = "[PASS]" if status == 'PASS' else "[WARN]" if status == 'WARN' else "[FAIL]"
                summary_lines.append(f"   {status_symbol} {module:30s} : {result.get('summary', '')}")
            
            if hasattr(analysis, 'all_recommendations') and analysis.all_recommendations:
                summary_lines.append("")
                summary_lines.append("//================================\\\\")
                summary_lines.append("||       Recommendations          ||")
                summary_lines.append("\\\\================================//")
                summary_lines.append("")
                for i, rec in enumerate(analysis.all_recommendations[:8], 1):
                    summary_lines.append(f"   {i}. {rec}")
                if len(analysis.all_recommendations) > 8:
                    summary_lines.append(f"   ... ({len(analysis.all_recommendations) - 8} more recommendations)")
            
            summary_lines.append("")
            summary_lines.append("//================================\\\\")
            summary_lines.append(f"   Output saved to : {args.output}")
            summary_lines.append(f"   Log saved to    : {args.output.replace('.json', '.log')}")
            summary_lines.append("\\\\================================//")
            summary_lines.append("")
            
            # Display in terminal
            summary_text = "\n".join(summary_lines)
            print(summary_text)
            
            # Save to log file
            log_file = args.output.replace('.json', '.log')
            with open(log_file, "w") as f:
                f.write(summary_text)

        elif args.command == "fix":
            if args.verbose:
                print(f"[INFO] Generating fix suggestions from: {args.input_path}")
            
            # Show tool status if requested
            if args.show_tool_status:
                from phredator.utils.tool_checker import ToolChecker
                checker = ToolChecker()
                print("\n" + "="*60)
                print("Tool Availability Status".center(60))
                print("="*60)
                checker.print_tool_status(verbose=True)
                print("="*60 + "\n")
            
            fixer = Fixer(args.input_path, input_reads=args.input_reads, check_tools=args.check_tools)
            fixes = fixer.run()
            
            # Save JSON output
            with open(args.output, "w") as f:
                f.write(fixes.to_json())
            
            # Generate and display terminal summary (Subread-style)
            summary_lines = []
            summary_lines.append("")
            summary_lines.append("//================================\\\\")
            summary_lines.append("||     Fix Suggestions Report     ||")
            summary_lines.append("\\\\================================//")
            summary_lines.append("")
            summary_lines.append(f"   Total fixes suggested : {len(fixes.fixes_applied)}")
            summary_lines.append("")
            
            if fixes.fixes_applied:
                summary_lines.append("//================================\\\\")
                summary_lines.append("||    Recommended Actions         ||")
                summary_lines.append("\\\\================================//")
                summary_lines.append("")
                for i, fix in enumerate(fixes.fixes_applied, 1):
                    priority_tag = f"[{fix['priority'].upper():6s}]"
                    summary_lines.append(f"   {i}. {priority_tag} {fix['description']}")
                    if fix.get('reason'):
                        summary_lines.append(f"      Reason  : {fix['reason']}")
                    if fix.get('tool_required'):
                        summary_lines.append(f"      Tool    : {fix['tool_required']}")
                    if fix.get('command'):
                        summary_lines.append(f"      Command : {fix['command'][:80]}")
                        if len(fix['command']) > 80:
                            summary_lines.append(f"                {fix['command'][80:]}")
                    summary_lines.append("")
            else:
                summary_lines.append("   No fixes needed - sample quality is acceptable.")
                summary_lines.append("")
            
            if hasattr(fixes, 'suggested_pipeline') and fixes.suggested_pipeline:
                summary_lines.append("//================================\\\\")
                summary_lines.append("||      Suggested Pipeline        ||")
                summary_lines.append("\\\\================================//")
                summary_lines.append("")
                for line in fixes.suggested_pipeline[:15]:
                    summary_lines.append(f"   {line}")
                if len(fixes.suggested_pipeline) > 15:
                    summary_lines.append(f"   ... ({len(fixes.suggested_pipeline) - 15} more lines)")
                summary_lines.append("")
            
            summary_lines.append("//================================\\\\")
            summary_lines.append(f"   Output saved to : {args.output}")
            summary_lines.append(f"   Log saved to    : {args.output.replace('.json', '.log')}")
            summary_lines.append("\\\\================================//")
            summary_lines.append("")
            
            # Display in terminal
            summary_text = "\n".join(summary_lines)
            print(summary_text)
            
            # Save to log file
            log_file = args.output.replace('.json', '.log')
            with open(log_file, "w") as f:
                f.write(summary_text)

                if args.check_tools and fixes.tool_availability:
                    ta = fixes.tool_availability
                    print(f"[INFO] Tools available: {ta['total_installed']}/{ta['total_installed'] + ta['total_missing']}")
                    if ta['missing']:
                        print(f"[WARN] Missing tools: {', '.join(ta['missing'])}")
                
                if fixes.read_length:
                    print(f"[INFO] Detected read length: {fixes.read_length}bp")
                if fixes.is_paired_end:
                    print(f"[INFO] Detected paired-end reads")

        elif args.command == "report":
            if args.verbose:
                print(f"[INFO] Generating report from: {args.input_path}")
            reporter = Reporter(args.input_path)
            reporter.generate(args.output, fmt=args.format)
            if args.verbose:
                print(f"[INFO] Report generated: {args.output}")
                print(f"[INFO] Format: {args.format.upper()}")

        elif args.command == "batch":
            # Handle multiple input formats:
            # 1. Single file path (list file)
            # 2. Multiple FastQC files/zips directly
            # 3. Glob pattern results
            
            sample_list = []
            
            # If single argument and it's a text file, read it as list
            if len(args.samples) == 1 and Path(args.samples[0]).is_file():
                sample_path = Path(args.samples[0])
                
                # Check if it's a text file (list) or a FastQC file
                if sample_path.suffix in ['.txt', '.list']:
                    # Read from list file (one sample per line)
                    with open(sample_path, 'r') as f:
                        sample_list = [line.strip() for line in f if line.strip() and not line.startswith('#')]
                else:
                    # It's a single FastQC file
                    sample_list = [str(sample_path)]
            else:
                # Multiple files provided directly
                sample_list = args.samples
            
            if args.verbose:
                print(f"[INFO] Loaded {len(sample_list)} samples")
            
            # Initialize batch processor
            batch_processor = BatchProcessor(
                sample_list=sample_list,
                output_dir=args.output_dir,
                organism=args.organism,
                experiment_type=args.experiment_type,
                check_tools=args.check_tools,
                parallel=args.parallel,
                dry_run=args.dry_run,
                verbose=args.verbose
            )
            
            # Process all samples
            report = batch_processor.process_all()
            
            # ANSI color codes
            RED = '\033[31m'  # Dark red instead of bright red
            BOLD = '\033[1m'
            RESET = '\033[0m'
            
            # Generate professional terminal summary with red Phredator title
            summary_lines = []
            summary_lines.append("")
            summary_lines.append(f"{RED}{BOLD}")
            summary_lines.append("  ██████╗ ██╗  ██╗██████╗ ███████╗██████╗  ██████╗ ████████╗ ██████╗ ██████╗ ")
            summary_lines.append("  ██╔══██╗██║  ██║██╔══██╗██╔════╝██╔══██╗██╔═══██╗╚══██╔══╝██╔═══██╗██╔══██╗")
            summary_lines.append("  ██████╔╝███████║██████╔╝█████╗  ██║  ██║███████║   ██║   ██║   ██║██████╔╝")
            summary_lines.append("  ██╔═══╝ ██╔══██║██╔══██╗██╔══╝  ██║  ██║██╔══██║   ██║   ██║   ██║██╔══██╗")
            summary_lines.append("  ██║     ██║  ██║██║  ██║███████╗██████╔╝██║  ██║   ██║   ╚██████╔╝██║  ██║")
            summary_lines.append("  ╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝╚═════╝ ╚═╝  ╚═╝   ╚═╝    ╚═════╝ ╚═╝  ╚═╝")
            summary_lines.append(f"{RESET}")
            summary_lines.append("")
            summary_lines.append("//================================\\\\")
            summary_lines.append("||   Batch Processing Summary     ||")
            summary_lines.append("\\\\================================//")
            summary_lines.append("")
            summary_lines.append(f"   Processed {report.successful}/{report.total_samples} samples:")
            
            if report.statistics:
                stats = report.statistics
                total_processed = stats.pass_count + stats.warn_count + stats.fail_count
                
                if total_processed > 0:
                    pass_pct = (stats.pass_count / total_processed) * 100
                    warn_pct = (stats.warn_count / total_processed) * 100
                    fail_pct = (stats.fail_count / total_processed) * 100
                    
                    summary_lines.append(f"   [PASS] {stats.pass_count} ({pass_pct:.0f}%)")
                    summary_lines.append(f"   [WARN] {stats.warn_count} ({warn_pct:.0f}%)")
                    summary_lines.append(f"   [FAIL] {stats.fail_count} ({fail_pct:.0f}%)")
            
            if report.failed > 0:
                summary_lines.append(f"   Failed  : {report.failed} samples")
            if report.skipped > 0:
                summary_lines.append(f"   Skipped : {report.skipped} samples")
            
            # Show per-sample details
            summary_lines.append("")
            summary_lines.append("//================================\\\\")
            summary_lines.append("||      Sample Details            ||")
            summary_lines.append("\\\\================================//")
            summary_lines.append("")
            
            for sample_result in report.sample_results[:20]:  # Show first 20 samples
                if sample_result.status == "success" and sample_result.analysis_json:
                    try:
                        with open(sample_result.analysis_json, 'r') as f:
                            analysis_data = json.load(f)
                        
                        status_symbol = "[PASS]" if sample_result.overall_status.upper() == "PASS" else "[WARN]" if sample_result.overall_status.upper() == "WARN" else "[FAIL]"
                        summary_lines.append(f"   {status_symbol} {sample_result.sample_name}")
                        summary_lines.append("")
                        
                        # Show ALL metrics with status symbols
                        if 'metrics' in analysis_data:
                            metrics = analysis_data['metrics']
                            
                            # Define metric display order
                            metric_order = [
                                'per_base_quality',
                                'gc_content', 
                                'duplication_levels',
                                'adapter_content',
                                'overrepresented_sequences',
                                'per_base_n_content',
                                'sequence_length_distribution',
                                'per_sequence_quality',
                                'per_tile_quality',
                                'per_sequence_gc_content',
                                'kmer_content'
                            ]
                            
                            for metric_key in metric_order:
                                if metric_key in metrics:
                                    metric = metrics[metric_key]
                                    status = metric.get('status', '').upper()
                                    summary = metric.get('summary', '')
                                    
                                    # Status symbol
                                    if status == 'PASS':
                                        symbol = "[PASS]"
                                    elif status == 'WARN':
                                        symbol = "[WARN]"
                                    elif status == 'FAIL':
                                        symbol = "[FAIL]"
                                    else:
                                        symbol = "[INFO]"
                                    
                                    # Format metric name nicely
                                    metric_name = metric_key.replace('_', ' ').title()
                                    
                                    summary_lines.append(f"      {symbol} {metric_name:30s} : {summary}")
                        
                        summary_lines.append("")
                    except:
                        pass
            
            if len(report.sample_results) > 20:
                summary_lines.append(f"   ... and {len(report.sample_results) - 20} more samples")
                summary_lines.append("")
            
            summary_lines.append("")
            summary_lines.append("//================================\\\\")
            summary_lines.append("||      Average Metrics           ||")
            summary_lines.append("\\\\================================//")
            summary_lines.append("")
            
            if report.statistics:
                stats = report.statistics
                
                if stats.gc_mean is not None:
                    gc_display = f"{stats.gc_mean:.1f}%"
                    if stats.gc_std is not None:
                        gc_display += f" ± {stats.gc_std:.1f}%"
                    summary_lines.append(f"   GC content   : {gc_display}")
                
                if stats.quality_mean is not None:
                    q_display = f"Q{stats.quality_mean:.1f}"
                    if stats.quality_std is not None:
                        q_display += f" ± {stats.quality_std:.1f}"
                    summary_lines.append(f"   Quality      : {q_display}")
                
                if stats.duplication_mean is not None:
                    dup_display = f"{stats.duplication_mean:.1f}%"
                    if stats.duplication_std is not None:
                        dup_display += f" ± {stats.duplication_std:.1f}%"
                    
                    # Add context for duplication based on experiment type
                    context = ""
                    if args.experiment_type and 'rna' in args.experiment_type.lower():
                        context = " (normal for RNA-seq)"
                    elif args.experiment_type and 'chip' in args.experiment_type.lower():
                        context = " (expected for ChIP-seq)"
                    
                    summary_lines.append(f"   Duplication  : {dup_display}{context}")
                
                # Show outliers if any
                if stats.outliers:
                    summary_lines.append("")
                    summary_lines.append("//================================\\\\")
                    summary_lines.append("||         Outliers               ||")
                    summary_lines.append("\\\\================================//")
                    summary_lines.append("")
                    
                    for outlier in stats.outliers[:10]:
                        summary_lines.append(f"   {outlier['sample']:20s} : {outlier['metric']} = {outlier['value']} ({outlier['reason']})")
                    
                    if len(stats.outliers) > 10:
                        summary_lines.append(f"   ... and {len(stats.outliers) - 10} more outliers")
            
            # Collect all recommendations from all samples
            all_recommendations = set()
            for sample_result in report.sample_results:
                if sample_result.analysis_json and sample_result.status == "success":
                    try:
                        with open(sample_result.analysis_json, 'r') as f:
                            analysis_data = json.load(f)
                            if 'all_recommendations' in analysis_data:
                                for rec in analysis_data['all_recommendations']:
                                    all_recommendations.add(rec)
                    except:
                        pass
            
            # Show common recommendations
            if all_recommendations:
                summary_lines.append("")
                summary_lines.append("//================================\\\\")
                summary_lines.append("||      Recommendations           ||")
                summary_lines.append("\\\\================================//")
                summary_lines.append("")
                
                for i, rec in enumerate(sorted(all_recommendations)[:10], 1):
                    summary_lines.append(f"   {i}. {rec}")
                
                if len(all_recommendations) > 10:
                    summary_lines.append(f"   ... and {len(all_recommendations) - 10} more recommendations")
            
            summary_lines.append("")
            summary_lines.append("//================================\\\\")
            summary_lines.append(f"   Output : {args.output_dir}")
            summary_lines.append(f"   Report : {args.output_dir}/batch_report.json")
            summary_lines.append("\\\\================================//")
            summary_lines.append("")
            
            # Display in terminal
            summary_text = "\n".join(summary_lines)
            print(summary_text)
            
            # Save to log file
            log_file = Path(args.output_dir) / "batch_summary.log"
            with open(log_file, "w") as f:
                f.write(summary_text)


        elif args.command == "pipeline":
            if args.verbose:
                print(f"[INFO] Starting pipeline for: {args.input_fastq}")
                if args.dry_run:
                    print(f"[INFO] DRY-RUN MODE - No fixes will be executed")
            
            runner = PipelineRunner(
                input_fastq=args.input_fastq,
                output_dir=args.output_dir,
                organism=args.organism,
                experiment_type=args.experiment_type,
                check_tools=args.check_tools,
                dry_run=args.dry_run,
                verbose=args.verbose
            )
            
            result = runner.run()
            
            # Print summary
            print("\n" + "="*70)
            print("PIPELINE EXECUTION SUMMARY".center(70))
            print("="*70)
            
            if not args.dry_run:
                print(f"\n📊 Metrics Comparison:")
                print(f"   ✓ Improved:  {result.metrics_improved}")
                print(f"   ✗ Degraded:  {result.metrics_degraded}")
                print(f"   = Unchanged: {len(result.comparisons) - result.metrics_improved - result.metrics_degraded}")
                
                if result.overall_improvement:
                    print(f"\n✅ OVERALL: Quality IMPROVED")
                else:
                    print(f"\n⚠️  OVERALL: No significant improvement")
                
                # Show specific improvements
                if result.comparisons:
                    print(f"\n📈 Detailed Changes:")
                    for comp in result.comparisons:
                        if comp['improved']:
                            print(f"   ✓ {comp['metric'].replace('_', ' ').title()}: {comp['before_status'].upper()} → {comp['after_status'].upper()}")
                
                print(f"\n📁 Output:")
                print(f"   Fixed file: {result.output_file}")
                print(f"   Results:    {args.output_dir}/pipeline_result.json")
            else:
                print(f"\n📋 Dry-run completed - {result.fixes_executed} fixes identified")
                print(f"   Run without --dry-run to execute fixes")
            
            print("="*70 + "\n")
            
            if args.verbose:
                print(f"[INFO] Complete results saved to: {args.output_dir}/")

        elif args.command == "batch_old":  # Keep old handler temporarily
            if args.verbose:
                print(f"[INFO] Batch parsing {len(args.input_paths)} samples")
            
            batch = BatchParser(args.input_paths)
            batch.parse_all()
            batch.save_batch_report(args.output)
            
            if args.verbose:
                print(f"[INFO] Batch report saved to {args.output}")
                stats = batch.get_aggregate_statistics()
                print(f"[INFO] Total samples: {stats['total_samples']}")
                print(f"[INFO] Average GC: {stats['gc_content']['mean']:.2f}%")

    except FileNotFoundError as e:
        print(f"[ERROR] File not found: {e}", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"[ERROR] Invalid JSON file: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"[ERROR] Invalid input: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] Unexpected error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    sys.exit(main())