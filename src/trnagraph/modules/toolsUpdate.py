#!/usr/bin/env python3

import os
import re
import sys
import shutil
import logging
import subprocess
from typing import Optional, Tuple

from packaging.version import Version, InvalidVersion

from . import env_check
from .. import __version__

# The tRNAgraph version `trnagraph update` was itself introduced in -- a release older than this
# has no `update` command to fall back on, so checking out a tag from before it would leave the
# user with a package that (post git-checkout) can no longer run the very command they just ran.
MIN_TAG_VERSION = "1.9.0"


class UpdateManager:
    '''
    Implements `trnagraph update`: updates this git checkout (the `main` branch by default, or an
    explicit --branch/--tag), then re-syncs the environment (`conda env update -f
    requirements.yaml --prune`) and the package's own install (`pip install -e .`). Mirrors the
    manual update procedure a developer would otherwise run by hand -- see docs/roadmap.md's
    "update CLI command" item.
    '''
    def __init__(self, args):
        self.args = args
        self.logger = logging.getLogger(__name__)
        self.project_root = env_check.get_project_root()
        if self.project_root is None:
            # Raised here, not in run(), because project_root is the `cwd` of every git command
            # this class issues -- the object has no meaning without it. Failing at construction
            # guarantees no git command runs at all, which matters: git searches upward from its
            # working directory, so under a non-editable install (where the guessed root sits
            # inside the Python installation) `git status`/`fetch`/`checkout`/`pull` would
            # resolve against whatever repository happens to enclose it -- commonly the user's
            # own project, with a `python -m venv .venv` layout.
            raise ValueError(
                "`trnagraph update` requires a source checkout: it updates the working tree with "
                "git and re-syncs the conda environment from requirements.yaml, neither of which "
                "exists in an installed distribution. This looks like a non-editable install "
                f"(the package is at {os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}). "
                "Clone the repository and `pip install -e .` from it to use this command, or "
                "upgrade with pip directly."
            )

    def run(self) -> None:
        self._check_git_available()
        self._check_clean_working_tree()

        tag_ref: Optional[str] = None
        if self.args.tag:
            tag_ref, tag_version = self._resolve_tag(self.args.tag)
            min_version = Version(MIN_TAG_VERSION)
            if tag_version < min_version:
                raise ValueError(
                    f"Refusing to update to '{self.args.tag}': `trnagraph update` was introduced in "
                    f"v{MIN_TAG_VERSION}, so it cannot fetch a release older than that. Check out the "
                    f"tag manually instead: git checkout {tag_ref}"
                )

        remote = self._resolve_remote()
        self.logger.info(f"Fetching latest refs from '{remote}'...")
        self._fetch(remote)

        if tag_ref:
            self.logger.info(f"Checking out tag '{tag_ref}'...")
            self._run(['git', 'checkout', tag_ref])
            self.logger.info(
                f"\nNote: '{tag_ref}' is a tag, not a branch, so git has left this checkout in "
                "'detached HEAD' state -- standard git behavior for checking out a tag, not an "
                "error. Commits made here won't belong to any branch and are easy to lose track "
                "of. To go back to tracking a branch, run e.g.:\n"
                "  git checkout main\n"
                "  git checkout dev"
            )
        else:
            branch = self._resolve_target_branch()
            self._check_branch_not_older_than_installed(remote, branch)
            self._checkout_branch(remote, branch)
            self.logger.info(f"Pulling latest '{branch}' from '{remote}'...")
            self._run(['git', 'pull', remote, branch])

        env_tool = 'mamba' if shutil.which('mamba') else 'conda'
        self.logger.info(f"Syncing conda environment ({env_tool} env update -f requirements.yaml --prune)...")
        self._conda_env_update(env_tool)

        self.logger.info("Reinstalling tRNAgraph package (pip install -e .)...")
        self._run([sys.executable, '-m', 'pip', 'install', '-e', '.'])

        self.logger.info("\nUpdate complete.")

    def _resolve_target_branch(self) -> str:
        '''
        Default to the currently checked-out branch when --branch isn't given, so `update`
        continues tracking whatever branch was last explicitly selected (via --branch) instead of
        always snapping back to a hardcoded 'main' -- git's own branch state IS the memory, no
        separate persisted preference needed. Detached HEAD (e.g. after a previous `trnagraph
        update --tag ...`) has no "current branch" to continue from, so that's a hard error
        rather than a silent guess.
        '''
        if self.args.branch:
            return self.args.branch
        result = self._run(['git', 'rev-parse', '--abbrev-ref', 'HEAD'], log_output=False)
        current = result.stdout.strip()
        if current == 'HEAD':
            raise ValueError(
                "Cannot determine a default branch to update: this checkout is in detached HEAD "
                "state (e.g. from a previous `trnagraph update --tag ...`). Pass --branch explicitly."
            )
        return current

    # -- git/environment helpers --

    def _check_git_available(self) -> None:
        if shutil.which('git') is None:
            raise ValueError("Error: 'git' is not installed or not in PATH.")

    def _check_clean_working_tree(self) -> None:
        result = self._run(['git', 'status', '--porcelain'], log_output=False)
        # '??' entries are untracked files, not uncommitted changes to anything git already
        # tracks -- those are left alone (git wouldn't touch them on checkout/pull either).
        dirty = [line for line in result.stdout.splitlines() if line and not line.startswith('??')]
        if dirty:
            raise ValueError(
                "Refusing to update: this checkout has uncommitted local changes to tracked files:\n"
                + '\n'.join(dirty)
                + "\nCommit, stash, or discard them first (e.g. `git stash`), then re-run `trnagraph update`."
            )

    def _resolve_remote(self) -> str:
        result = self._run(['git', 'remote'], log_output=False)
        remotes = [r for r in result.stdout.split() if r]
        if not remotes:
            raise ValueError('Error: no git remote is configured for this repository.')
        tracking = env_check.get_tracking_remote_name(self.project_root)
        if tracking and tracking in remotes:
            return tracking
        return remotes[0]

    def _fetch(self, remote: str) -> None:
        '''
        `git fetch <remote> --tags --prune` fails outright if a *local* tag ref has diverged from
        the same-named tag on the remote (git treats tags as immutable pointers and refuses to
        silently move one -- "! [rejected] vX.Y.Z -> vX.Y.Z (would clobber existing tag)"). This
        can happen if a release tag was ever recreated upstream, or an old checkout cached a tag
        before that happened. Since `trnagraph update`'s whole purpose is to make this checkout
        match the remote, that's exactly the case where forcing the local tag to match the
        remote's is correct -- so on that specific rejection, retry once with `--force` (which
        only affects tags/remote-tracking refs here, not local branches).
        '''
        cmd = ['git', 'fetch', remote, '--tags', '--prune']
        result = self._run(cmd, check=False)
        if result.returncode != 0:
            if 'clobber existing tag' not in result.stderr:
                raise ValueError(f"Command failed ({' '.join(cmd)}): {result.stderr.strip() or result.stdout.strip()}")
            self.logger.info(
                "A local tag has diverged from the one on the remote (likely recreated "
                f"upstream since this checkout last fetched) -- re-fetching from '{remote}' with "
                "tags forced to match..."
            )
            self._run(cmd + ['--force'])

    def _check_branch_not_older_than_installed(self, remote: str, branch: str) -> None:
        '''
        Refuse to update to a branch whose OWN `__version__` is older than the version currently
        installed -- this is what actually catches `main` being stale (main is deliberately not
        fast-forwarded to dev yet, see docs/roadmap.md's "update CLI command" item), rather than
        a hardcoded floor like `--tag`'s MIN_TAG_VERSION. Peeks at `<remote>/<branch>`'s own
        `src/trnagraph/__init__.py` via `git show` -- no checkout happens, so nothing is mutated
        if this refuses.

        Two different "can't tell" cases are handled deliberately differently:
        - The ref itself (`<remote>/<branch>`) doesn't exist at all -- fails open (does nothing)
          and lets the real checkout step below surface git's own clear error instead of a
          confusing one here.
        - The ref exists but has no `src/trnagraph/__init__.py` at all -- fails CLOSED (refuses):
          on this project's actual history, that means the target commit predates the versioned
          package structure entirely (added well after v1.0.0), so it's certainly older, not
          merely "unknown". This is exactly the real-world case that motivated this check: `main`
          currently points at pre-package-restructure history with no `__init__.py` at all, so
          treating "file missing" as "unknown, proceed anyway" would silently reproduce the exact
          bug this guard exists to catch.
        A malformed/unparsable `__version__` string, on the other hand, is genuinely ambiguous
        (the file exists, so this isn't pre-versioning history) and fails open.
        '''
        ref_check = self._run(['git', 'rev-parse', '--verify', '--quiet', f'{remote}/{branch}'], log_output=False, check=False)
        if ref_check.returncode != 0:
            return

        result = self._run(
            ['git', 'show', f'{remote}/{branch}:src/trnagraph/__init__.py'], log_output=False, check=False,
        )
        installed_version = Version(__version__)
        if result.returncode != 0:
            raise ValueError(
                f"Refusing to update to '{branch}': couldn't find src/trnagraph/__init__.py on "
                f"'{remote}/{branch}' at all, which means it predates this project's versioned "
                f"package structure -- it's almost certainly older than the version currently "
                f"installed ({installed_version}). If you really want it anyway, check it out "
                f"manually instead: git checkout {branch}"
            )
        match = re.search(r'__version__\s*=\s*[\'"]([^\'"]+)[\'"]', result.stdout)
        if not match:
            return
        try:
            branch_version = Version(match.group(1))
        except InvalidVersion:
            return
        if branch_version < installed_version:
            raise ValueError(
                f"Refusing to update to '{branch}': its own version ({branch_version}) is older than "
                f"the version currently installed ({installed_version}) -- '{branch}' on '{remote}' is "
                f"behind. If you really want it anyway, check it out manually instead: git checkout {branch}"
            )

    def _checkout_branch(self, remote: str, branch: str) -> None:
        local_branches = self._run(['git', 'branch', '--list', branch], log_output=False).stdout
        if local_branches.strip():
            self.logger.info(f"Checking out '{branch}'...")
            self._run(['git', 'checkout', branch])
        else:
            self.logger.info(f"Creating local branch '{branch}' tracking '{remote}/{branch}'...")
            self._run(['git', 'checkout', '-b', branch, '--track', f'{remote}/{branch}'])

    def _resolve_tag(self, tag_arg: str) -> Tuple[str, Version]:
        raw = tag_arg[1:] if tag_arg[:1] in ('v', 'V') else tag_arg
        try:
            version = Version(raw)
        except InvalidVersion as e:
            raise ValueError(f"'{tag_arg}' is not a valid version tag (expected e.g. v1.9.0 or 1.9.0).") from e
        ref = tag_arg if tag_arg[:1] in ('v', 'V') else f'v{tag_arg}'
        return ref, version

    def _conda_env_update(self, env_tool: str) -> None:
        # Prefer mamba's much faster dependency resolver when it's available (e.g. any
        # miniforge/mambaforge install) -- same env update semantics, drop-in CLI-compatible
        # with conda. Falls back to conda when mamba isn't installed rather than requiring it.
        if shutil.which(env_tool) is None:
            raise ValueError(f"Error: '{env_tool}' is not installed or not in PATH.")
        # Cannot be None in practice -- the constructor refuses to build without a source
        # checkout -- but get_requirements_path() is Optional, so handle it as one path.
        req_path = env_check.get_requirements_path()
        if req_path is None or not os.path.exists(req_path):
            raise ValueError(f'Error: requirements.yaml not found at {req_path}.')
        # Streamed rather than self._run()'s buffered subprocess.run(capture_output=True): a
        # conda/mamba solve can legitimately take a long time, and with output fully buffered
        # until the process exits, that's indistinguishable from an actual hang. Live output
        # (conda/mamba's own progress lines) makes it visible that it's still working.
        self._run_streaming([env_tool, 'env', 'update', '-f', req_path, '--prune'])

    def _run_streaming(self, cmd) -> None:
        '''
        Like _run(), but streams stdout/stderr line-by-line as the subprocess produces them
        (into the logger, so both the console -- unless --quiet -- and the log file see it in
        near-real-time) instead of buffering everything until the process exits. Reserved for
        steps that can legitimately run long (currently just the conda/mamba env sync); the git
        commands elsewhere in this class are normally fast enough that buffering is fine.
        '''
        proc = subprocess.Popen(cmd, cwd=self.project_root, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        output_lines = []
        for line in proc.stdout:
            line = line.rstrip()
            output_lines.append(line)
            self.logger.info(line)
        returncode = proc.wait()
        if returncode != 0:
            raise ValueError(f"Command failed ({' '.join(cmd)}): {output_lines[-1] if output_lines else '(no output)'}")

    def _run(self, cmd, log_output: bool = True, check: bool = True) -> subprocess.CompletedProcess:
        # check=False -- failures are surfaced via our own `check`/return-code handling below,
        # not subprocess's CalledProcessError, so we can attach a friendlier message.
        result = subprocess.run(cmd, cwd=self.project_root, capture_output=True, text=True, check=False)
        if log_output:
            if result.stdout.strip():
                self.logger.info(result.stdout.rstrip())
            if result.stderr.strip():
                self.logger.info(result.stderr.rstrip())
        if check and result.returncode != 0:
            raise ValueError(f"Command failed ({' '.join(cmd)}): {result.stderr.strip() or result.stdout.strip()}")
        return result
