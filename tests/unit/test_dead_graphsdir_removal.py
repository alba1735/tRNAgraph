"""Regression tests for roadmap.md Phase 2's dead `graphsdir` creation removal. Both
`preprocess map` (MapSamples.main()) and `analyze build` (AnalysisPipeline.run()) unconditionally
created a `graphsdir` that neither command ever writes into (`build` already had a comment
marking this as dead legacy scaffolding) -- confirmed via toolsMap.py/adataBuild.py inspection
that nothing in either method's code path touches `self.expinfo.graphsdir`/`graph_path(...)`
paths for writing. These tests confirm the directory is no longer created by either command,
while `resultsdir` (which both commands genuinely write result files into) still is."""
from types import SimpleNamespace
from unittest.mock import patch

from trnagraph.modules.adataBuild import AnalysisPipeline
from trnagraph.modules.toolsMap import MapSamples


def _make_map_args(tmp_path, samplenames):
    fastqs = []
    for name in samplenames:
        fq = tmp_path / f"{name}.fastq"
        fq.write_text("@read\nACGT\n+\nIIII\n")
        fastqs.append(fq)
    manifest = tmp_path / "manifest.txt"
    lines = ["fastq\tsample\tgroup"]
    lines += [f"{fq}\t{name}\t{name}" for fq, name in zip(fastqs, samplenames)]
    manifest.write_text("\n".join(lines) + "\n")

    return SimpleNamespace(
        mode='map', output=str(tmp_path / "exp"), database=str(tmp_path / "db"),
        input=str(manifest), lazy=False, minnontrnasize=20,
        local=False, threads=1, skipcheck=True, bamdir=str(tmp_path / "bam"), quiet=True,
    )


def test_mapsamples_main_does_not_create_graphsdir(tmp_path):
    args = _make_map_args(tmp_path, ["sample0"])
    mapper = MapSamples(args)

    with patch.object(MapSamples, "mapsamples", lambda self: None):
        mapper.main()

    import os
    assert not os.path.exists(mapper.expinfo.graphsdir)
    assert os.path.exists(mapper.expinfo.resultsdir)


def _make_build_args(tmp_path):
    return SimpleNamespace(
        database=str(tmp_path / "db"), output=str(tmp_path / "exp" / "exp.h5ad"),
        input=str(tmp_path / "metadata.txt"), gtf=None, bed=[],
        bamdir=str(tmp_path / "bam"), threads=1,
        minnontrnasize=20, maxmismatches=None, mincoverage=None, filtermultimapped=False,
        pairs=None, hubonly=False, hub=False, filterother=False, quiet=True,
    )


def test_analysispipeline_run_does_not_create_graphsdir(tmp_path):
    args = _make_build_args(tmp_path)
    pipeline = AnalysisPipeline(args)

    patches = [
        patch.object(AnalysisPipeline, name, lambda self, *a, **k: None)
        for name in ["makefeaturebed", "countfeatures", "run_deseq2", "counttypes",
                     "run_unique_deseq2", "gettrnacoverage", "write_runinfo"]
    ]
    with patches[0], patches[1], patches[2], patches[3], patches[4], patches[5], patches[6], \
         patch.object(pipeline.trnainfo, "getorgtype", lambda: "euk"):
        pipeline.run()

    import os
    assert not os.path.exists(pipeline.expinfo.graphsdir)
    assert os.path.exists(pipeline.expinfo.resultsdir)
