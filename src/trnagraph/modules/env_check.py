import os
import sys
import json
import time
import shutil
import subprocess
import re
import importlib.metadata
from typing import Dict, Optional, Tuple

from .. import __version__

def _is_trnagraph_source_root(path: str) -> bool:
    """
    True only if `path` is *this project's* source checkout: it must contain both the package
    source (`src/trnagraph/__init__.py`) and a `pyproject.toml` naming tRNAgraph.

    Checking for a `.git` directory would not do. Callers pass the result to git as a working
    directory, and git searches upward from there, so any path nested inside an unrelated
    repository looks like a valid checkout to it -- which is exactly the failure this guards.
    """
    if not os.path.isfile(os.path.join(path, 'src', 'trnagraph', '__init__.py')):
        return False
    pyproject = os.path.join(path, 'pyproject.toml')
    if not os.path.isfile(pyproject):
        return False
    try:
        import tomllib
        with open(pyproject, 'rb') as f:
            name = tomllib.load(f).get('project', {}).get('name', '')
    except Exception:
        # An unreadable or malformed pyproject.toml is not evidence of anything; treat it the
        # same as an absent one rather than guessing.
        return False
    return name.lower() == 'trnagraph'

def get_project_root() -> Optional[str]:
    """
    Absolute path of the tRNAgraph source checkout this package was installed from, or None when
    there isn't one.

    The path is derived from this file's location (assumed to be `src/trnagraph/modules/`) and
    then positively identified via _is_trnagraph_source_root(). The identification step is the
    point: under a non-editable `pip install .` the same walk lands on `<venv>/lib/pythonX.Y`,
    and returning that unchecked let `trnagraph update` run git commands whose working directory
    resolved to whatever repository happened to enclose the install location.

    None means "installed as a distribution, no source tree available" -- requirements.yaml,
    the git history and the branch/tag are all absent, and every caller must say what it does
    without them rather than proceeding against a path that merely exists.
    """
    current_dir = os.path.dirname(os.path.abspath(__file__))
    candidate = os.path.dirname(os.path.dirname(os.path.dirname(current_dir)))
    return candidate if _is_trnagraph_source_root(candidate) else None

def get_requirements_path() -> Optional[str]:
    """
    Path to requirements.yaml in the source checkout, or None when there is no checkout.
    requirements.yaml is not part of the built distribution, so it exists only alongside source.
    """
    project_root = get_project_root()
    return os.path.join(project_root, 'requirements.yaml') if project_root else None

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

def get_binary_version(cmd: str, flag: str, regex: str, timeout: float = 20.0) -> Tuple[Optional[str], Optional[str]]:
    """
    Try to get the version of a binary. Returns (version, failure_reason) -- exactly one is
    None. Distinguishing WHY a lookup failed (timed out vs no regex match vs some other
    subprocess error) instead of collapsing every failure into one generic message matters in
    practice: a Python-based CLI tool (e.g. umi_tools) has to boot a full interpreter and import
    numpy/pandas/pysam just to print --version, which can push past a short timeout when the
    machine is under heavy concurrent load (e.g. a `trnagraph` pipeline already running on it) --
    that failure mode looks identical to "the output format doesn't match" unless it's reported
    separately. `timeout` defaults to a generous 20s for exactly that reason -- a fast C binary
    (bedtools, bowtie2, samtools) costs nothing extra from the higher ceiling, but a Python-based
    one needs the headroom.
    """
    if not shutil.which(cmd):
        return None, "not found"

    args = [cmd]
    if flag:
        args.append(flag)

    try:
        # Capture both stdout and stderr as some tools print version to stderr (e.g. fastp)
        result = subprocess.run(args, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None, f"timed out after {timeout}s"
    except Exception as e:
        return None, str(e)

    output = result.stdout + result.stderr
    match = re.search(regex, output)
    if match:
        return match.group(1), None
    return None, "version string not found in command output"

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

    found_version, failure_reason = get_binary_version(cmd, info["flag"], info["regex"])

    if found_version:
        if compare_versions(found_version, expected_version, op):
            return True, f"Found {found_version}"
        return False, f"Expected {op}{expected_version}, found {found_version}"

    return False, f"Command '{cmd}' found but version could not be determined ({failure_reason})"

def validate_environment():
    """
    Validates that the environment matches requirements.yaml.
    """
    req_path = get_requirements_path()
    if req_path is None or not os.path.exists(req_path):
        # requirements.yaml ships with the source, not with the built distribution, so a
        # non-editable install has nothing to validate against. Skip silently rather than
        # warning: the environment may well be correct, we simply cannot confirm it here.
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


def get_tracking_remote_name(project_root: str) -> Optional[str]:
    """
    Best-effort lookup of the git remote the current branch tracks (e.g. 'github' for a branch
    tracking 'github/dev'). Returns None on any failure -- detached HEAD, no upstream configured,
    git unavailable, or `project_root` not a git checkout at all -- rather than raising; callers
    decide how to handle a missing remote (e.g. `update`'s own multi-remote fallback).
    """
    try:
        result = subprocess.run(
            ['git', '-C', project_root, 'rev-parse', '--abbrev-ref', '--symbolic-full-name', '@{u}'],
            capture_output=True, text=True, timeout=5,
        )
    except Exception:
        return None
    if result.returncode != 0:
        return None
    upstream = result.stdout.strip()
    if '/' not in upstream:
        return None
    return upstream.split('/')[0]


def get_remote_url(project_root: str, remote_name: str) -> Optional[str]:
    """Best-effort lookup of a configured remote's URL. Returns None on any failure."""
    try:
        result = subprocess.run(
            ['git', '-C', project_root, 'remote', 'get-url', remote_name],
            capture_output=True, text=True, timeout=5,
        )
    except Exception:
        return None
    if result.returncode != 0:
        return None
    return result.stdout.strip() or None


def get_current_branch(project_root: str) -> Optional[str]:
    """
    Best-effort lookup of the currently checked-out branch name. Returns None on any failure
    (detached HEAD, git unavailable, or `project_root` not a git checkout at all) rather than
    raising -- callers decide how to handle "no current branch" (e.g. get_version_channel()'s
    exact-tag/nightly fallback, or `trnagraph update`'s own hard error for its branch default).
    """
    try:
        result = subprocess.run(
            ['git', '-C', project_root, 'rev-parse', '--abbrev-ref', 'HEAD'],
            capture_output=True, text=True, timeout=5,
        )
    except Exception:
        return None
    if result.returncode != 0:
        return None
    branch = result.stdout.strip()
    return None if branch in ('', 'HEAD') else branch


def _get_exact_tag(project_root: str) -> Optional[str]:
    """The tag HEAD is exactly at, if any (e.g. after `trnagraph update --tag v1.9.0`, which
    leaves a detached HEAD). None if HEAD isn't exactly at a tag, or on any failure."""
    try:
        result = subprocess.run(
            ['git', '-C', project_root, 'describe', '--tags', '--exact-match'],
            capture_output=True, text=True, timeout=5,
        )
    except Exception:
        return None
    if result.returncode != 0:
        return None
    return result.stdout.strip() or None


def _get_short_hash(project_root: str) -> Optional[str]:
    try:
        result = subprocess.run(
            ['git', '-C', project_root, 'rev-parse', '--short', 'HEAD'],
            capture_output=True, text=True, timeout=5,
        )
    except Exception:
        return None
    if result.returncode != 0:
        return None
    return result.stdout.strip() or None


def get_version_channel(project_root: Optional[str]) -> str:
    """
    Short label describing where this checkout's code is actually coming from, since
    __version__ alone doesn't distinguish 'the official main-branch release' from 'a stale/ahead
    dev or feature checkout with the same version string bumped in' -- see this project's own
    `main` being deliberately behind `dev` during stabilization (docs/roadmap.md's update-tool
    item), which is exactly the kind of thing this label makes visible at a glance:
      - 'main' branch -> 'stable'
      - 'dev' branch -> 'prerelease'
      - any other branch -> 'nightly @ <short-hash>' (a fork/feature branch, not one of the two
        known release channels, so the hash is the only precise identifier)
      - detached HEAD (e.g. from `trnagraph update --tag`) exactly at a release tag -> that tag
      - detached HEAD not at a tag -> 'nightly @ <short-hash>', same reasoning as other branches
      - no source checkout at all (a non-editable install) -> 'installed'
    """
    if project_root is None:
        # A built distribution carries no git history, so there is no branch, tag or hash to
        # report. 'stable' would be a lie -- a wheel built from a dirty dev branch would claim
        # it -- and 'nightly' implies a hash we do not have.
        return 'installed'
    branch = get_current_branch(project_root)
    if branch == 'main':
        return 'stable'
    if branch == 'dev':
        # 'prerelease' rather than 'beta': it pairs with 'stable' as a release status, whereas
        # 'beta' implies a maturity grade this project does not otherwise use. It describes the
        # branch's position relative to main -- not semantic-release's `prerelease` setting,
        # which is false for dev (see pyproject.toml), so a release cut here is a plain x.y.z.
        return 'prerelease'
    if branch is None:
        tag = _get_exact_tag(project_root)
        if tag:
            return tag
    short_hash = _get_short_hash(project_root)
    return f'nightly @ {short_hash}' if short_hash else 'nightly'


_UPDATE_CHECK_INTERVAL_SECONDS = 24 * 60 * 60

def _update_check_cache_path() -> str:
    cache_home = os.environ.get('XDG_CACHE_HOME') or '~/.cache'
    return os.path.join(os.path.expanduser(cache_home), 'trnagraph', 'update_check.json')

def _latest_stable_tag(ls_remote_output: str) -> Optional[str]:
    """Pick the highest non-prerelease semver tag out of a `git ls-remote --tags` listing."""
    from packaging.version import Version, InvalidVersion

    best_version = None
    best_tag = None
    for line in ls_remote_output.splitlines():
        if 'refs/tags/' not in line:
            continue
        tag = line.rsplit('refs/tags/', 1)[1].strip()
        if tag.endswith('^{}'):
            tag = tag[:-3]
        raw = tag[1:] if tag[:1] in ('v', 'V') else tag
        try:
            version = Version(raw)
        except InvalidVersion:
            continue
        if version.is_prerelease or version.is_devrelease:
            continue
        if best_version is None or version > best_version:
            best_version, best_tag = version, tag
    return best_tag

def _print_update_notice_if_newer(latest_tag: Optional[str]) -> None:
    if not latest_tag:
        return
    from packaging.version import Version, InvalidVersion

    raw_latest = latest_tag[1:] if latest_tag[:1] in ('v', 'V') else latest_tag
    try:
        if Version(raw_latest) > Version(__version__):
            print(f"A newer version ({latest_tag}) is available -- run `trnagraph update` to upgrade.")
    except InvalidVersion:
        pass

def warn_if_output_capture_suspected(isatty_fn=None, environ=None) -> Optional[str]:
    """
    Returns a warning message if stdout/stderr look non-interactive despite plausibly being
    inside a human-attended screen/tmux session -- a strong signal of an output-capturing
    wrapper (most commonly `conda run` without --no-capture-output, which defaults to buffering
    all output until the wrapped process exits, making a long-running command look completely
    silent for its entire duration regardless of the actual terminal). Returns None when
    genuinely interactive, or when non-interactive with no $STY/$TMUX -- that's an ordinary
    background/cron/CI invocation, nothing suspicious about it. Doesn't print directly, so
    callers can log/print/test the message independently.
    """
    if isatty_fn is None:
        isatty_fn = lambda: sys.stdout.isatty() or sys.stderr.isatty()
    if isatty_fn():
        return None
    if environ is None:
        environ = os.environ
    if environ.get('STY') or environ.get('TMUX'):
        return (
            "NOTE: output appears non-interactive even though this looks like a screen/tmux "
            "session. If this was launched via `conda run`, live progress display needs "
            "`conda run --no-capture-output` -- by default conda run buffers all output until "
            "the command finishes, which can look like nothing is happening for a long time."
        )
    return None


def check_for_updates() -> None:
    """
    Best-effort, non-blocking notice if a newer tRNAgraph release exists on the configured git
    remote (the remote the current branch tracks). Compares tags via
    `packaging.version.Version`, the same logic `update`'s own version guard uses -- never string
    comparison. Must never raise or block a command: every failure mode (no network, git
    unavailable, no remote configured, an unreachable/timing-out remote, not a git checkout at
    all) is swallowed silently. Cached to a per-user timestamp file
    (~/.cache/trnagraph/update_check.json, respecting $XDG_CACHE_HOME) so the network is only
    touched once every 24h, not on every invocation.
    """
    try:
        cache_path = _update_check_cache_path()
        now = time.time()

        if os.path.exists(cache_path):
            with open(cache_path, 'r') as f:
                cached = json.load(f)
            if now - float(cached.get('checked_at', 0)) < _UPDATE_CHECK_INTERVAL_SECONDS:
                _print_update_notice_if_newer(cached.get('latest_tag'))
                return

        project_root = get_project_root()
        if project_root is None:
            # No checkout means no configured remote to ask, so there is nothing to check
            # against. Return before touching the network at all.
            return
        remote_name = get_tracking_remote_name(project_root)
        if not remote_name:
            return
        remote_url = get_remote_url(project_root, remote_name)
        if not remote_url:
            return

        result = subprocess.run(
            ['git', 'ls-remote', '--tags', remote_url],
            capture_output=True, text=True, timeout=5,
        )
        if result.returncode != 0:
            return

        latest_tag = _latest_stable_tag(result.stdout)

        os.makedirs(os.path.dirname(cache_path), exist_ok=True)
        with open(cache_path, 'w') as f:
            json.dump({'checked_at': now, 'latest_tag': latest_tag}, f)

        _print_update_notice_if_newer(latest_tag)
    except Exception:
        pass
