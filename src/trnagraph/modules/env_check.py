import os
import sys
import shutil
import subprocess
import re
import importlib.metadata
from typing import Dict, Optional, Tuple

def get_requirements_path() -> str:
    """
    Locate requirements.yaml relative to this file.
    Assumes this file is in src/trnagraph/modules/
    and requirements.yaml is in the project root.
    """
    current_dir = os.path.dirname(os.path.abspath(__file__))

    project_root = os.path.dirname(os.path.dirname(os.path.dirname(current_dir)))
    return os.path.join(project_root, 'requirements.yaml')

def parse_requirements(req_path: str) -> Dict[str, Tuple[str, str]]:
    """
    Parse requirements.yaml to extract package names and versions.
    Returns a dict {package: (operator, version)}.
    """
    requirements = {}
    if not os.path.exists(req_path):
        # Fallback or warning if file not found
        return requirements

    with open(req_path, 'r') as f:
        in_dependencies = False
        for line in f:
            line = line.strip()
            if line.startswith('dependencies:'):
                in_dependencies = True
                continue
            
            if in_dependencies and line.startswith('-'):
                # Remove '- ' and quotes
                dep = line.lstrip('- ').strip('"\'')
                
                # Match package name followed by operator and version
                # Operators: =, ==, >=, <=, >, <
                match = re.match(r'^([a-zA-Z0-9\-_]+)\s*([<>=!]+)\s*(.+)$', dep)
                if match:
                    pkg = match.group(1)
                    op = match.group(2)
                    ver = match.group(3)
                    # Normalize single = to == for comparison logic if needed, 
                    # but conda uses = for pinning. We'll handle it.
                    requirements[pkg] = (op, ver)
                else:
                    # Maybe just package name?
                    pass
    return requirements

def compare_versions(installed: str, required: str, operator: str) -> bool:
    """
    Compare two version strings using the given operator.
    """
    def parse_ver(v):
        # Remove build info if present (e.g. 1.2.3=py39_0 -> 1.2.3)
        v = v.split('=')[0]
        # Split by dots and convert to int where possible
        parts = []
        for part in re.split(r'[.-]', v):
            if part.isdigit():
                parts.append(int(part))
            else:
                parts.append(part)
        return parts

    v1 = parse_ver(installed)
    v2 = parse_ver(required)
    
    if operator == '=' or operator == '==':
        # For equality, we might want to be loose if installed has more detail?
        # But usually exact match or prefix match.
        # Let's do prefix match for '='
        if operator == '=':
             return str(installed).startswith(str(required))
        return v1 == v2
    elif operator == '>=':
        return v1 >= v2
    elif operator == '>':
        return v1 > v2
    elif operator == '<=':
        return v1 <= v2
    elif operator == '<':
        return v1 < v2
    
    return False

def check_python_package(package: str, requirement: Tuple[str, str]) -> Tuple[bool, str]:
    """
    Check if a python package is installed and matches the version.
    """
    op, expected_version = requirement
    
    # Mapping for package names that differ from import names or metadata names
    # Most match, but some might need adjustment
    pkg_map = {
        "scikit-learn": "scikit-learn",
        "umap-learn": "umap-learn",
        "sra-tools": None, # Binary
        "bowtie2": None, # Binary
        "samtools": None, # Binary
        "bedtools": None, # Binary
        "fastp": None, # Binary
        "gffread": None, # Binary
        "git": None, # Binary
        "infernal": None, # Binary
        "ucsc-bedgraphtobigwig": None, # Binary
        "umi_tools": None, # Binary
        "python": "python", # Special case
    }
    
    if package == "python":
        import platform
        current_ver = platform.python_version()
        if compare_versions(current_ver, expected_version, op):
            return True, f"Found {current_ver}"
        return False, f"Expected {op}{expected_version}, found {current_ver}"

    if package in pkg_map:
        if pkg_map[package] is None:
            return True, "Binary check skipped here" # Handled separately
        pkg_name = pkg_map[package]
    else:
        pkg_name = package

    try:
        installed_version = importlib.metadata.version(pkg_name)
        if compare_versions(installed_version, expected_version, op):
            return True, f"Found {installed_version}"
        else:
            return False, f"Expected {op}{expected_version}, found {installed_version}"
    except importlib.metadata.PackageNotFoundError:
        return False, "Not installed"

def get_binary_version(cmd: str, flag: str, regex: str) -> Optional[str]:
    """
    Try to get version of a binary.
    """
    if not shutil.which(cmd):
        return None
    
    try:
        # Capture both stdout and stderr as some tools print version to stderr (e.g. fastp)
        args = [cmd]
        if flag:
            args.append(flag)
            
        result = subprocess.run(args, capture_output=True, text=True, timeout=5)
        output = result.stdout + result.stderr
        match = re.search(regex, output)
        if match:
            return match.group(1)
    except Exception:
        pass
    return None

def check_binary_package(package: str, requirement: Tuple[str, str]) -> Tuple[bool, str]:
    """
    Check if a binary package is installed and matches the version.
    """
    op, expected_version = requirement
    
    # Special handling for ucsc-bedgraphtobigwig which has mismatched conda/binary versions
    if package == "ucsc-bedgraphtobigwig":
        expected_version = "2.10"
        op = ">="

    bin_map = {
        "bedtools": {"cmd": "bedtools", "flag": "--version", "regex": r"bedtools v([\d\.]+)"},
        "bowtie2": {"cmd": "bowtie2", "flag": "--version", "regex": r"version ([\d\.]+)"},
        "fastp": {"cmd": "fastp", "flag": "--version", "regex": r"fastp ([\d\.]+)"},
        "gffread": {"cmd": "gffread", "flag": "--version", "regex": r"([\d\.]+)"},
        "git": {"cmd": "git", "flag": "--version", "regex": r"git version ([\d\.]+)"},
        "samtools": {"cmd": "samtools", "flag": "--version", "regex": r"samtools ([\d\.]+)"},
        "sra-tools": {"cmd": "fastq-dump", "flag": "--version", "regex": r": ([\d\.]+)"},
        "ucsc-bedgraphtobigwig": {"cmd": "bedGraphToBigWig", "flag": "", "regex": r"v ([\d\.]+)"},
        "infernal": {"cmd": "cmscan", "flag": "-h", "regex": r"INFERNAL ([\d\.]+)"},
        "umi_tools": {"cmd": "umi_tools", "flag": "--version", "regex": r": ([\d\.]+)"},
    }

    if package not in bin_map:
        return True, "Not a known binary"

    info = bin_map[package]
    cmd = info["cmd"]
    
    if not shutil.which(cmd):
        return False, f"Command '{cmd}' not found"

    found_version = get_binary_version(cmd, info["flag"], info["regex"])
    
    if found_version:
        if compare_versions(found_version, expected_version, op):
            return True, f"Found {found_version}"
        return False, f"Expected {op}{expected_version}, found {found_version}"
    
    return False, f"Command '{cmd}' found but version could not be determined"

def validate_environment():
    """
    Validates that the environment matches requirements.yaml.
    """
    req_path = get_requirements_path()
    if not os.path.exists(req_path):
        # If we can't find requirements.yaml, we can't validate against it.
        # This might happen in installed package context if file not included.
        # We'll skip validation or print a warning.
        # print(f"Warning: Could not find requirements.yaml at {req_path}")
        return

    requirements = parse_requirements(req_path)
    
    errors = []
    
    # Define which are binaries
    binaries = [
        "bedtools", "bowtie2", "fastp", "gffread", "git", 
        "samtools", "sra-tools", "ucsc-bedgraphtobigwig", "infernal", "umi_tools"
    ]

    print("Validating environment...")
    
    for pkg, req in requirements.items():
        if pkg in binaries:
            ok, msg = check_binary_package(pkg, req)
        else:
            ok, msg = check_python_package(pkg, req)
            if msg == "Binary check skipped here":
                continue

        if not ok:
            errors.append(f"{pkg}: {msg}")

    if errors:
        print("\n\033[91mError: Environment validation failed.\033[0m")
        print("The following dependencies are missing or have incorrect versions:")
        for err in errors:
            print(f"  - {err}")
        
        in_conda = os.environ.get("CONDA_PREFIX") is not None
        # Try to determine the active conda env name
        conda_env = os.environ.get("CONDA_DEFAULT_ENV")
        if not conda_env and os.environ.get("CONDA_PREFIX"):
            conda_env = os.path.basename(os.environ.get("CONDA_PREFIX").rstrip(os.sep))

        required_env = "trnagraph"

        if not in_conda:
            print("\n\033[93mWarning: No Conda environment detected.\033[0m")
            print(f"Please activate/install the {required_env} environment:")
            print(f"  conda activate {required_env}")
            sys.exit(1)
        elif conda_env != required_env:
            print("\n\033[91mError: Incorrect Conda environment.\033[0m")
            print(f"Current conda env: {conda_env or os.environ.get('CONDA_PREFIX')}")
            print("Please switch to the required environment:")
            print(f"  conda activate {required_env}")
            sys.exit(1)
        else:
            print(f"\nUsing Conda environment '{required_env}'.")
            
        sys.exit(1)
    else:
        print("Environment validation passed.")
        pass
