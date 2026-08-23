"""Regression test for roadmap.md's tempfile/TMPDIR item. toolsMap.py's map_sample() and
toolsTrackHub.py's convertbam() both built their `samtools sort -T` prefix under
tempfile.gettempdir() (system /tmp unless $TMPDIR is set). tRAX historically hit a server
/tmp-fills-up failure and fixed it by writing temp files into the invocation directory instead --
tRNAgraph inherited the same class of risk. Fix: a shared helper that defaults the sort-temp
directory to the output file's own directory instead of the system temp dir."""
import os

from trnagraph.modules.toolsTG import sort_temp_dir


def test_sort_temp_dir_uses_the_outputs_own_directory():
    assert sort_temp_dir("/some/dir/output.bam") == "/some/dir"


def test_sort_temp_dir_resolves_relative_paths_to_an_absolute_directory():
    result = sort_temp_dir("output.bam")
    assert os.path.isabs(result)
    assert result == os.getcwd()


def test_sort_temp_dir_resolves_relative_subdirectories():
    result = sort_temp_dir(os.path.join("nested", "output.bam"))
    assert result == os.path.join(os.getcwd(), "nested")
