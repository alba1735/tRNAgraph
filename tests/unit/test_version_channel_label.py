"""Regression test for roadmap.md's update-tool item: show which branch/channel (main=stable,
dev=prerelease, anything else=nightly) tRNAgraph is running from, in --version output and every
command's log, so the confusion this project's error_log_update_server.txt shows (main being far
behind dev, invisibly) is easier to spot at a glance."""
import subprocess

from trnagraph.modules.env_check import get_current_branch, get_version_channel


def _run(cwd, *args):
    subprocess.run(['git', *args], cwd=str(cwd), check=True, capture_output=True)


def _init_repo(repo):
    repo.mkdir()
    _run(repo, 'init', '-q', '-b', 'main')
    _run(repo, 'config', 'user.email', 'test@example.com')
    _run(repo, 'config', 'user.name', 'Test')
    (repo / 'f.txt').write_text('1')
    _run(repo, 'add', 'f.txt')
    _run(repo, 'commit', '-q', '-m', 'init')


def test_get_current_branch_reports_main(tmp_path):
    repo = tmp_path / 'repo'
    _init_repo(repo)

    assert get_current_branch(str(repo)) == 'main'


def test_get_current_branch_reports_dev(tmp_path):
    repo = tmp_path / 'repo'
    _init_repo(repo)
    _run(repo, 'checkout', '-q', '-b', 'dev')

    assert get_current_branch(str(repo)) == 'dev'


def test_get_current_branch_returns_none_on_detached_head(tmp_path):
    repo = tmp_path / 'repo'
    _init_repo(repo)
    sha = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=str(repo), capture_output=True, text=True, check=True).stdout.strip()
    _run(repo, 'checkout', '-q', sha)

    assert get_current_branch(str(repo)) is None


def test_version_channel_main_is_stable(tmp_path):
    repo = tmp_path / 'repo'
    _init_repo(repo)

    assert get_version_channel(str(repo)) == 'stable'


def test_version_channel_dev_is_prerelease(tmp_path):
    repo = tmp_path / 'repo'
    _init_repo(repo)
    _run(repo, 'checkout', '-q', '-b', 'dev')

    assert get_version_channel(str(repo)) == 'prerelease'


def test_version_channel_other_branch_is_nightly_with_hash(tmp_path):
    repo = tmp_path / 'repo'
    _init_repo(repo)
    _run(repo, 'checkout', '-q', '-b', 'my-feature')
    short_hash = subprocess.run(
        ['git', 'rev-parse', '--short', 'HEAD'], cwd=str(repo), capture_output=True, text=True, check=True
    ).stdout.strip()

    channel = get_version_channel(str(repo))

    assert channel == f'nightly @ {short_hash}'


def test_version_channel_detached_head_at_exact_tag_shows_the_tag(tmp_path):
    repo = tmp_path / 'repo'
    _init_repo(repo)
    _run(repo, 'tag', 'v1.9.0')
    sha = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=str(repo), capture_output=True, text=True, check=True).stdout.strip()
    _run(repo, 'checkout', '-q', sha)

    assert get_version_channel(str(repo)) == 'v1.9.0'


def test_version_channel_detached_head_not_at_a_tag_is_nightly_with_hash(tmp_path):
    repo = tmp_path / 'repo'
    _init_repo(repo)
    sha = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=str(repo), capture_output=True, text=True, check=True).stdout.strip()
    short_hash = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=str(repo), capture_output=True, text=True, check=True).stdout.strip()
    _run(repo, 'checkout', '-q', sha)

    assert get_version_channel(str(repo)) == f'nightly @ {short_hash}'
