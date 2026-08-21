#!/usr/bin/env python3

import os
import sys
import shutil
import logging
import subprocess
from typing import Optional, Tuple

from packaging.version import Version, InvalidVersion

from . import env_check

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
        self._run(['git', 'fetch', remote, '--tags', '--prune'])

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
            branch = self.args.branch or "main"
            self._checkout_branch(remote, branch)
            self.logger.info(f"Pulling latest '{branch}' from '{remote}'...")
            self._run(['git', 'pull', remote, branch])

        env_tool = 'mamba' if shutil.which('mamba') else 'conda'
        self.logger.info(f"Syncing conda environment ({env_tool} env update -f requirements.yaml --prune)...")
        self._conda_env_update(env_tool)

        self.logger.info("Reinstalling tRNAgraph package (pip install -e .)...")
        self._run([sys.executable, '-m', 'pip', 'install', '-e', '.'])

        self.logger.info("\nUpdate complete.")

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
        req_path = env_check.get_requirements_path()
        if not os.path.exists(req_path):
            raise ValueError(f'Error: requirements.yaml not found at {req_path}.')
        self._run([env_tool, 'env', 'update', '-f', req_path, '--prune'])

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
