import os
import yaml
from typing import Dict, Any, Optional, List, Tuple
from dataclasses import dataclass
from difflib import get_close_matches


@dataclass
class OrganismProfile:
    name: str
    gc_content: Dict[str, Any]
    quality: Dict[str, Any]
    duplication: Dict[str, Any]
    adapters: Dict[str, Any]
    n_content: Dict[str, Any]
    overrepresented: Dict[str, Any]
    read_length: Dict[str, Any]
    notes: str = ""
    assembly: str = ""
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'OrganismProfile':
        return cls(
            name=data.get('name', 'Unknown'),
            gc_content=data.get('gc_content', {}),
            quality=data.get('quality', {}),
            duplication=data.get('duplication', {}),
            adapters=data.get('adapters', {}),
            n_content=data.get('n_content', {}),
            overrepresented=data.get('overrepresented', {}),
            read_length=data.get('read_length', {}),
            notes=data.get('notes', ''),
            assembly=data.get('assembly', '')
        )


@dataclass
class ExperimentProfile:
    name: str
    abbrev: str
    description: str
    quality: Dict[str, Any]
    duplication: Dict[str, Any]
    adapters: Dict[str, Any]
    special: Dict[str, Any]
    gc_content: Optional[Dict[str, Any]] = None
    notes: str = ""
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'ExperimentProfile':
        return cls(
            name=data.get('name', 'Unknown'),
            abbrev=data.get('abbrev', 'unknown'),
            description=data.get('description', ''),
            quality=data.get('quality', {}),
            duplication=data.get('duplication', {}),
            adapters=data.get('adapters', {}),
            special=data.get('special', {}),
            gc_content=data.get('gc_content'),
            notes=data.get('notes', '')
        )


class ProfileLoader:
    
    def __init__(self):
        current_dir = os.path.dirname(os.path.abspath(__file__))
        self.config_dir = os.path.join(os.path.dirname(current_dir), 'config')
        self.organism_dir = os.path.join(self.config_dir, 'organisms')
        self.experiment_dir = os.path.join(self.config_dir, 'experiment_types')
    
    def _normalize_input(self, value: str) -> str:
        # rna-seq, RNA_SEQ, RnaSeq all become rnaseq
        return value.strip().lower().replace('-', '').replace('_', '').replace(' ', '')
    
    def _find_closest_match(self, input_value: str, available_values: List[str]) -> Optional[str]:
        normalized_input = self._normalize_input(input_value)
        normalized_available = {self._normalize_input(v): v for v in available_values}
        
        if normalized_input in normalized_available:
            return normalized_available[normalized_input]
        
        matches = get_close_matches(normalized_input, normalized_available.keys(), n=1, cutoff=0.6)
        if matches:
            return normalized_available[matches[0]]
        
        return None
    
    def _get_suggestions(self, input_value: str, available_values: List[str], top_n: int = 3) -> List[str]:
        normalized_input = self._normalize_input(input_value)
        normalized_available = {self._normalize_input(v): v for v in available_values}
        
        matches = get_close_matches(normalized_input, normalized_available.keys(), n=top_n, cutoff=0.4)
        return [normalized_available[m] for m in matches]
    
    def load_organism_profile(self, organism: str) -> Optional[OrganismProfile]:
        available_organisms = self.list_organisms()
        
        matched_organism = self._find_closest_match(organism, available_organisms)
        
        if matched_organism:
            organism_normalized = matched_organism.lower()
        else:
            organism_normalized = organism.strip().lower()
        
        profile_path = os.path.join(self.organism_dir, f"{organism_normalized}.yaml")
        
        if not os.path.exists(profile_path):
            suggestions = self._get_suggestions(organism, available_organisms)
            print(f"\nError: Organism '{organism}' not found.")
            if suggestions:
                print(f"Did you mean: {', '.join(suggestions)}?")
            print(f"\nAvailable organisms: {', '.join(available_organisms[:10])}")
            if len(available_organisms) > 10:
                print(f"... and {len(available_organisms) - 10} more")
            print("Use 'phredator list-organisms' to see all available organisms.\n")
            return None
            
        try:
            with open(profile_path, 'r') as f:
                data = yaml.safe_load(f)
            return OrganismProfile.from_dict(data)
        except Exception as e:
            print(f"Warning: Could not load organism profile '{organism}': {e}")
            return None
    
    def load_experiment_profile(self, experiment_type: str) -> Optional[ExperimentProfile]:
        available_experiments = self.list_experiment_types()
        
        matched_experiment = self._find_closest_match(experiment_type, available_experiments)
        
        if matched_experiment:
            experiment_normalized = matched_experiment.lower()
        else:
            experiment_normalized = experiment_type.strip().lower()
        
        profile_path = os.path.join(self.experiment_dir, f"{experiment_normalized}.yaml")
        
        if not os.path.exists(profile_path):
            suggestions = self._get_suggestions(experiment_type, available_experiments)
            print(f"\nError: Experiment type '{experiment_type}' not found.")
            if suggestions:
                print(f"Did you mean: {', '.join(suggestions)}?")
            print(f"\nAvailable experiment types: {', '.join(available_experiments)}")
            print("Use 'phredator list-organisms --detailed' to see experiment types.\n")
            return None
            
        try:
            with open(profile_path, 'r') as f:
                data = yaml.safe_load(f)
            return ExperimentProfile.from_dict(data)
        except Exception as e:
            print(f"Warning: Could not load experiment profile '{experiment_type}': {e}")
            return None
    
    def list_organisms(self) -> list:
        if not os.path.exists(self.organism_dir):
            return []
        
        profiles = []
        for filename in os.listdir(self.organism_dir):
            if filename.endswith('.yaml'):
                profiles.append(filename.replace('.yaml', ''))
        return sorted(profiles)
    
    def list_experiment_types(self) -> list:
        if not os.path.exists(self.experiment_dir):
            return []
        
        profiles = []
        for filename in os.listdir(self.experiment_dir):
            if filename.endswith('.yaml'):
                profiles.append(filename.replace('.yaml', ''))
        return sorted(profiles)
    
    def get_combined_thresholds(self, 
                                organism: Optional[str] = None,
                                experiment_type: Optional[str] = None) -> Dict[str, Any]:
        thresholds = self._get_default_thresholds()
        
        if organism:
            org_profile = self.load_organism_profile(organism)
            if org_profile:
                thresholds['gc_content'] = org_profile.gc_content
                thresholds['quality'].update(org_profile.quality)
                thresholds['duplication'].update(org_profile.duplication)
                thresholds['adapters'].update(org_profile.adapters)
                thresholds['n_content'] = org_profile.n_content
                thresholds['organism_name'] = org_profile.name
        
        if experiment_type:
            exp_profile = self.load_experiment_profile(experiment_type)
            if exp_profile:
                thresholds['quality'].update(exp_profile.quality)
                thresholds['duplication'].update(exp_profile.duplication)
                thresholds['adapters'].update(exp_profile.adapters)
                thresholds['special'] = exp_profile.special
                thresholds['experiment_name'] = exp_profile.name
                
                if exp_profile.gc_content:
                    thresholds['gc_content'] = exp_profile.gc_content
        
        return thresholds
    
    def _get_default_thresholds(self) -> Dict[str, Any]:
        return {
            'gc_content': {
                'mean': 50.0,
                'range': [35, 65],
                'tolerance': 5.0
            },
            'quality': {
                'q20_threshold': 0.80,
                'q28_threshold': 0.80,
                'q30_threshold': 0.75,
                'mean_quality_min': 28
            },
            'duplication': {
                'acceptable': 20,
                'warning': 35,
                'critical': 50,
                'check_duplicates': True
            },
            'adapters': {
                'acceptable': 5,
                'warning': 10,
                'critical': 15,
                'required': True
            },
            'n_content': {
                'max_per_base': 5,
                'max_total': 1
            },
            'special': {},
            'organism_name': 'Generic',
            'experiment_name': 'Generic'
        }
