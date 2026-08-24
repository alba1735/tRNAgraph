"""Organism-mode handling must be explicit for euk, arch, bact and mito alike.

The Sprinzl position tables are defined twice -- literally in
`toolsGetCoverage.py` and compositionally in `toolsTDatabase._init_positions()` --
and nothing tied the two together, so a change to one (the planned e19->e27
variable-arm extension is the obvious candidate) would silently desynchronise the
database from the coverage step.

Separately, an unrecognised organism mode used to fall through to eukaryotic
behaviour silently in three places: `gettnanums()`'s else branch,
`pos_maps.get(mode, pos_maps["euk"])`, and `getorgtype()` when the database has no
recorded mode. Since `--orgmode` was an unvalidated free string, `-s bacteria` or
`-s Bact` would build a eukaryotic database without complaint.
"""
import pytest

from trnagraph.modules import toolsGetCoverage as gc

MODES = ("euk", "arch", "bact", "mito")


def _db_tables():
    from trnagraph.modules.toolsTDatabase import tRNADatabaseBuilder as TRNADatabase
    obj = TRNADatabase.__new__(TRNADatabase)
    obj._init_positions()
    return obj.pos_maps


def test_the_two_position_tables_stay_in_sync():
    db = _db_tables()
    literal = {"euk": gc.eukpositions, "arch": gc.archpositions,
               "bact": gc.bactpositions, "mito": gc.mitopositions}
    for mode in MODES:
        assert db[mode] == literal[mode], f"{mode} tables have diverged"


def test_every_supported_mode_has_a_table_in_both_places():
    db = _db_tables()
    assert set(db) == set(MODES)
    for mode in MODES:
        assert gc.positions_for(mode), f"no coverage table for {mode}"


def test_unknown_organism_mode_is_rejected_not_silently_treated_as_euk():
    with pytest.raises(ValueError, match="bacteria"):
        gc.positions_for("bacteria")


def test_getorgtype_survives_a_blank_line_and_reports_unknown_modes(tmp_path):
    from trnagraph.modules.toolsMap import trnadatabase
    db = trnadatabase(str(tmp_path / "db"))

    # a trailing blank line must not raise IndexError on fields[0]
    open(db.dbinfo, "w").write("orgmode\tbact\n\n")
    assert db.getorgtype() == "bact"

    # an unrecognised recorded mode is an error, not a silent euk fallback
    open(db.dbinfo, "w").write("orgmode\tbacteria\n")
    with pytest.raises(ValueError, match="bacteria"):
        db.getorgtype()


def test_makedb_rejects_an_unknown_organism_mode(tmp_path):
    """`-s bacteria` used to build a eukaryotic database without complaint."""
    from types import SimpleNamespace
    from trnagraph.modules.toolsTDatabase import tRNADatabaseBuilder

    args = SimpleNamespace(
        output=str(tmp_path / "db"), genome=None, trnaout=None, trnafa=None,
        namemap=None, addtrna=None, addseqs=None, forcecca=False,
        orgmode="bacteria", threads=1,
    )
    with pytest.raises(ValueError, match="bacteria"):
        tRNADatabaseBuilder(args)


def test_makedb_accepts_every_supported_mode(tmp_path):
    from types import SimpleNamespace
    from trnagraph.modules.toolsTDatabase import tRNADatabaseBuilder

    for mode in MODES:
        args = SimpleNamespace(
            output=str(tmp_path / f"db_{mode}"), genome=None, trnaout=None,
            trnafa=None, namemap=None, addtrna=None, addseqs=None,
            forcecca=False, orgmode=mode, threads=1,
        )
        assert tRNADatabaseBuilder(args).orgmode == mode
