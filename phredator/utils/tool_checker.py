

import subprocess
import shutil
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass


@dataclass
class ToolInfo:
    """Information about a bioinformatics tool"""
    name: str
    command: str
    version_flag: str
    is_installed: bool = False
    version: Optional[str] = None
    install_conda: Optional[str] = None
    install_pip: Optional[str] = None
    description: str = ""


class ToolChecker:
    
    def __init__(self):
        self._checked = {}
        self.tools = {
            'fastqc': ToolInfo(
                name="FastQC",
                command="fastqc",
                version_flag="--version",
                install_conda="conda install -c bioconda fastqc",
                description="Quality control tool for high throughput sequence data"
            ),
            'cutadapt': ToolInfo(
                name="Cutadapt",
                command="cutadapt",
                version_flag="--version",
                install_conda="conda install -c bioconda cutadapt",
                install_pip="pip install cutadapt",
                description="Tool to remove adapters from high-throughput sequencing reads"
            ),
            'trimmomatic': ToolInfo(
                name="Trimmomatic",
                command="trimmomatic",
                version_flag="-version",
                install_conda="conda install -c bioconda trimmomatic",
                description="Flexible read trimming tool for Illumina NGS data"
            ),
            'fastp': ToolInfo(
                name="fastp",
                command="fastp",
                version_flag="--version",
                install_conda="conda install -c bioconda fastp",
                description="Ultra-fast all-in-one FASTQ preprocessor"
            ),
            'picard': ToolInfo(
                name="Picard",
                command="picard",
                version_flag="MarkDuplicates --version",
                install_conda="conda install -c bioconda picard",
                description="Java-based command-line utilities for manipulating HTS data"
            ),
            'samtools': ToolInfo(
                name="Samtools",
                command="samtools",
                version_flag="--version",
                install_conda="conda install -c bioconda samtools",
                description="Tools for manipulating alignments in SAM/BAM format"
            ),
            'bbduk': ToolInfo(
                name="BBDuk",
                command="bbduk.sh",
                version_flag="-h",
                install_conda="conda install -c bioconda bbmap",
                description="Decontamination and quality filtering tool"
            )
        }
    
    def check_tool(self, tool_key: str) -> bool:
        if tool_key in self._checked:
            return self._checked[tool_key]
        
        if tool_key not in self.tools:
            self._checked[tool_key] = False
            return False
        
        tool = self.tools[tool_key]
        command_path = shutil.which(tool.command)
        
        if command_path is None:
            tool.is_installed = False
            self._checked[tool_key] = False
            return False
        try:
            # Some tools print version to stderr, some to stdout
            result = subprocess.run(
                f"{tool.command} {tool.version_flag}",
                shell=True,
                capture_output=True,
                text=True,
                timeout=5,
                errors='replace'  # Handle encoding errors gracefully
            )
            version_output = (result.stdout or "") + (result.stderr or "")
            
            # Extract version (first line usually)
            if version_output:
                version_line = version_output.strip().split('\n')[0]
                tool.version = version_line[:100]  # Limit length
            else:
                tool.version = "installed (version unknown)"
            
            tool.is_installed = True
            self._checked[tool_key] = True
            return True
            
        except (subprocess.TimeoutExpired, subprocess.SubprocessError, FileNotFoundError):
            # Command exists but version check failed - assume installed
            tool.is_installed = True
            tool.version = "installed (version check failed)"
            self._checked[tool_key] = True
            return True
    
    def check_all_tools(self) -> Dict[str, bool]:
        results = {}
        for tool_key in self.tools.keys():
            results[tool_key] = self.check_tool(tool_key)
        return results
    
    def get_tool_info(self, tool_key: str) -> Optional[ToolInfo]:
        return self.tools.get(tool_key)
    
    def get_tool_version(self, tool_key: str) -> Optional[str]:
        if tool_key not in self.tools:
            return None
        self.check_tool(tool_key)
        return self.tools[tool_key].version
    
    def get_installed_tools(self) -> List[str]:
        installed = []
        for tool_key in self.tools.keys():
            if self.check_tool(tool_key):
                installed.append(tool_key)
        return installed
    
    def print_tool_status(self, verbose: bool = False) -> None:
        self.check_all_tools()
        
        print("\n" + "="*70)
        print("Tool Availability Check")
        print("="*70)
        
        installed = []
        missing = []
        
        for key, tool in self.tools.items():
            if tool.is_installed:
                installed.append((key, tool))
            else:
                missing.append((key, tool))
        
        # Print installed tools
        if installed:
            print("\n✓ Installed Tools:")
            for key, tool in installed:
                if verbose:
                    print(f"  • {tool.name:15} {tool.version}")
                else:
                    print(f"  ✓ {tool.name}")
        
        # Print missing tools
        if missing:
            print("\n✗ Missing Tools:")
            for key, tool in missing:
                print(f"  ✗ {tool.name:15} - {tool.description}")
        
        # Print installation suggestions
        if missing:
            print("\n" + "-"*70)
            print("Installation Suggestions:")
            print("-"*70)
            for key, tool in missing:
                print(f"\n{tool.name}:")
                if tool.install_conda:
                    print(f"  Conda:  {tool.install_conda}")
                if tool.install_pip:
                    print(f"  Pip:    {tool.install_pip}")
        
        print("\n" + "="*70)
        print(f"Summary: {len(installed)} installed, {len(missing)} missing")
        print("="*70 + "\n")
    
    def filter_tool_categories(self, required_tools: List[str]) -> Tuple[List[str], List[str]]:
        pass  # docstring removed
        available = []
        missing = []
        
        for tool_key in required_tools:
            if self.check_tool(tool_key):
                available.append(tool_key)
            else:
                missing.append(tool_key)
        
        return available, missing
    
    def suggest_alternatives(self, tool_category: str) -> List[str]:
        pass  # docstring removed
        # Map categories to tool alternatives (in order of preference)
        alternatives = {
            'adapter_removal': ['cutadapt', 'fastp', 'trimmomatic', 'bbduk'],
            'quality_trim': ['fastp', 'trimmomatic', 'cutadapt'],
            'deduplication': ['picard', 'samtools'],
            'contamination': ['bbduk', 'fastqc']
        }
        
        if tool_category not in alternatives:
            return []
        
        # Return installed tools in preference order
        installed = []
        for tool_key in alternatives[tool_category]:
            if self.check_tool(tool_key):
                installed.append(tool_key)
        
        return installed
