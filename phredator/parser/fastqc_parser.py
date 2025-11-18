import os
import zipfile
import json
from typing import Dict, List
from dataclasses import dataclass, field, asdict

@dataclass
class FastQCReport:
    sample_name: str
    # Basic Statistics
    filename: str = ""
    file_type: str = ""
    encoding: str = ""
    total_sequences: int = 0
    total_bases: str = ""
    sequences_flagged_poor_quality: int = 0
    sequence_length: str = ""
    percent_gc: int = 0
    # Per base sequence quality
    per_base_quality: Dict[str, Dict[str, float]] = field(default_factory=dict)
    # Per base sequence content (ATCG% at each position)
    per_base_sequence_content: Dict[str, Dict[str, float]] = field(default_factory=dict)
    # Per sequence quality scores (distribution of mean quality per read)
    per_sequence_quality_scores: Dict[int, float] = field(default_factory=dict)
    # Per base N content (% of N calls at each position)
    per_base_n_content: Dict[str, float] = field(default_factory=dict)
    # Sequence Length Distribution
    sequence_length_distribution: Dict[str, float] = field(default_factory=dict)
    # GC content as distribution (GC% -> count) plus calculated mean
    gc_content_mean: float = 0.0
    gc_content_distribution: Dict[int, float] = field(default_factory=dict)
    # Duplication: total % deduplicated library and level distribution
    total_deduplicated_percentage: float = 100.0
    duplication_levels: Dict[str, float] = field(default_factory=dict)
    adapter_content: Dict[str, float] = field(default_factory=dict)
    overrepresented_sequences: List[str] = field(default_factory=list)
    
    def to_json(self) -> str:

        return json.dumps(asdict(self), indent=4)

class FastQCParser:
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.sample_name = os.path.basename(filepath).split('_fastqc')[0]
        self.data = FastQCReport(sample_name=self.sample_name)

    def parse(self) -> FastQCReport:

        if not os.path.exists(self.filepath):
            raise FileNotFoundError(f"File not found: {self.filepath}")

        # Handle zip file
        if self.filepath.endswith('.zip'):
            with zipfile.ZipFile(self.filepath, 'r') as z:
                # Find fastqc_data.txt (might be in root or subdirectory)
                data_file = None
                for name in z.namelist():
                    if name.endswith('fastqc_data.txt'):
                        data_file = name
                        break
                
                if not data_file:
                    raise ValueError("fastqc_data.txt missing in zip")
                
                with z.open(data_file) as f:
                    self._parse_fastqc_data(f.read().decode('utf-8'))
        else:
            # Handle folder with fastqc_data.txt
            data_path = os.path.join(self.filepath, 'fastqc_data.txt')
            if not os.path.exists(data_path):
                raise FileNotFoundError(f"{data_path} missing")
            with open(data_path, 'r') as f:
                self._parse_fastqc_data(f.read())

        return self.data

    def _parse_fastqc_data(self, content: str):
        """Internal method to parse fastqc_data.txt content."""
        lines = content.splitlines()
        section = None

        for line in lines:
            line = line.strip()
            if line.startswith('>>Basic Statistics'):
                section = 'basic_statistics'
                continue
            elif line.startswith('>>Per base sequence quality'):
                section = 'per_base_quality'
                continue
            elif line.startswith('>>Per base sequence content'):
                section = 'per_base_sequence_content'
                continue
            elif line.startswith('>>Per sequence quality scores'):
                section = 'per_sequence_quality_scores'
                continue
            elif line.startswith('>>Per base N content'):
                section = 'per_base_n_content'
                continue
            elif line.startswith('>>Sequence Length Distribution'):
                section = 'sequence_length_distribution'
                continue
            elif line.startswith('>>Per sequence GC content'):
                section = 'gc_content'
                continue
            elif line.startswith('>>Sequence Duplication Levels'):
                section = 'duplication_levels'
                continue
            elif line.startswith('>>Adapter Content'):
                section = 'adapter_content'
                continue
            elif line.startswith('>>Overrepresented sequences'):
                section = 'overrepresented_sequences'
                continue
            elif line.startswith('>>END_MODULE'):
                section = None
                continue

            if section == 'basic_statistics' and line and not line.startswith('#'):
                # Format: "Measure\tValue"
                parts = line.split('\t')
                if len(parts) >= 2:
                    measure = parts[0].strip()
                    value = parts[1].strip()
                    
                    if measure == 'Filename':
                        self.data.filename = value
                    elif measure == 'File type':
                        self.data.file_type = value
                    elif measure == 'Encoding':
                        self.data.encoding = value
                    elif measure == 'Total Sequences':
                        try:
                            self.data.total_sequences = int(value)
                        except ValueError:
                            pass
                    elif measure == 'Total Bases':
                        self.data.total_bases = value
                    elif measure == 'Sequences flagged as poor quality':
                        try:
                            self.data.sequences_flagged_poor_quality = int(value)
                        except ValueError:
                            pass
                    elif measure == 'Sequence length':
                        self.data.sequence_length = value
                    elif measure == '%GC':
                        try:
                            self.data.percent_gc = int(value)
                        except ValueError:
                            pass

            elif section == 'per_base_quality' and line and not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 3:
                    base_range = parts[0]
                    mean_quality = float(parts[1])
                    median_quality = float(parts[2])
                    self.data.per_base_quality[base_range] = {
                        'mean': mean_quality,
                        'median': median_quality
                    }

            elif section == 'per_base_sequence_content' and line and not line.startswith('#'):
                # Format: "#Base\tG\tA\tT\tC"
                # Data: "1\t22.36\t28.37\t27.67\t21.6"
                parts = line.split()
                if len(parts) >= 5:
                    try:
                        base_position = parts[0]  # e.g., "1", "10-14"
                        g_percent = float(parts[1])
                        a_percent = float(parts[2])
                        t_percent = float(parts[3])
                        c_percent = float(parts[4])
                        self.data.per_base_sequence_content[base_position] = {
                            'G': g_percent,
                            'A': a_percent,
                            'T': t_percent,
                            'C': c_percent
                        }
                    except ValueError:
                        continue

            elif section == 'per_sequence_quality_scores' and line and not line.startswith('#'):
                # Format: "#Quality\tCount"
                # Data: "33\t4.0", "34\t8883.0"
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        quality = int(parts[0])
                        count = float(parts[1])
                        self.data.per_sequence_quality_scores[quality] = count
                    except ValueError:
                        continue

            elif section == 'per_base_n_content' and line and not line.startswith('#'):
                # Format: "#Base\tN-Count"
                # Data: "1\t0.0", "10-14\t0.0"
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        base_position = parts[0]
                        n_count = float(parts[1])
                        self.data.per_base_n_content[base_position] = n_count
                    except ValueError:
                        continue

            elif section == 'sequence_length_distribution' and line and not line.startswith('#'):
                # Format: "#Length\tCount"
                # Data: "150\t10000.0" or "40-49\t123.0"
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        length = parts[0]  # Can be "150" or "40-49"
                        count = float(parts[1])
                        self.data.sequence_length_distribution[length] = count
                    except ValueError:
                        continue

            elif section == 'gc_content' and line and not line.startswith('#'):
                # Format: "#GC Content\tCount"
                # Parse full distribution: GC% (0-100) -> count of sequences
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        gc_percent = int(parts[0])  # GC percentage (0-100)
                        count = float(parts[1])      # Number of sequences
                        if count > 0:
                            self.data.gc_content_distribution[gc_percent] = count
                    except ValueError:
                        continue  # skip malformed lines

            elif section == 'duplication_levels' and line:
                # Special case: "#Total Deduplicated Percentage\t0.03"
                if line.startswith('#Total Deduplicated Percentage'):
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        try:
                            self.data.total_deduplicated_percentage = float(parts[1])
                        except ValueError:
                            continue
                # Skip other header lines
                elif line.startswith('#'):
                    continue
                # Data lines: "1\t100.0", "2\t0.0", ">10\t5.0", etc.
                else:
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            duplication_level = parts[0]  # e.g., "1", "2", ">10", ">10k+"
                            percentage = float(parts[1])   # Percentage of total sequences at this level
                            self.data.duplication_levels[duplication_level] = percentage
                        except ValueError:
                            continue

            elif section == 'adapter_content' and line and not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 2:
                    adapter = parts[0]
                    try:
                        fraction = float(parts[1])
                        self.data.adapter_content[adapter] = fraction
                    except ValueError:
                        continue

            elif section == 'overrepresented_sequences' and line and not line.startswith('#'):
                seq = line.split()[0]  # sequence is the first column
                self.data.overrepresented_sequences.append(seq)

        # Calculate mean GC content from distribution
        self._calculate_gc_mean()

    def _calculate_gc_mean(self):

        if not self.data.gc_content_distribution:
            return
        
        total_sequences = sum(self.data.gc_content_distribution.values())
        if total_sequences == 0:
            return
        
        weighted_sum = sum(gc_pct * count 
                          for gc_pct, count in self.data.gc_content_distribution.items())
        self.data.gc_content_mean = weighted_sum / total_sequences
