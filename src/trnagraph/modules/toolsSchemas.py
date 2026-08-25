'''
Pydantic data models shared across tRNAgraph's variant-selection and size-split-merge
machinery (see toolsTG.py's parse_variant()/build_variant_view() and adataBuild.py's
compute_variant_contribution()/merge_variant_into_adata()). Kept in a dedicated module,
following the tools*.py convention, so these schemas have one canonical, validated definition
instead of being duplicated as plain dataclasses in the modules that use them.
'''

import re
from typing import Any, Dict, List, Optional, Literal, Tuple
import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict, RootModel, field_validator, model_validator

_TAG_PATTERN = re.compile(r'^(full|[uo]\d+)$')


class VariantTag(BaseModel):
    '''
    Parsed representation of a `--variant` CLI string, e.g. "raw:u60" -> norm='raw', tag='u60'.
    '''
    model_config = ConfigDict(frozen=True)

    raw: str
    norm: Literal['norm', 'raw', 'allfeatures', 'vst']
    tag: str

    @field_validator('tag')
    @classmethod
    def _validate_tag_format(cls, v: str) -> str:
        if not _TAG_PATTERN.match(v):
            raise ValueError(f"Split tag '{v}' must be 'full' or match '<u|o><cutoff>' (e.g. 'u60', 'o50').")
        return v


class VariantContribution(BaseModel):
    '''
    Coverage matrices, numeric per-obs columns, and count/sizefactor summaries for one
    size-split variant, computed by AnnDataBuilder.compute_variant_contribution() and merged
    into an existing target AnnData object by adataBuild.merge_variant_into_adata().
    '''
    model_config = ConfigDict(arbitrary_types_allowed=True)

    x_raw: pd.DataFrame
    x_norm: pd.DataFrame
    # None for a read-length split variant: all-feature normalization is computed only for the
    # complete variant, since a split has its non-tRNA features excluded entirely.
    x_norm_allfeatures: Optional[pd.DataFrame] = None
    x_vst: Optional[np.ndarray] = None
    obsm_counts: pd.DataFrame
    sizefactors_trna: Dict[str, Any]
    sizefactors_allfeatures: Optional[Dict[str, Any]] = None
    type_counts: pd.DataFrame
    type_real_counts: pd.DataFrame
    amino_counts: pd.DataFrame
    anticodon_counts: pd.DataFrame
    nontrna_counts: pd.DataFrame


class GraphFilterConfig(BaseModel):
    '''
    Validates the `--config` filter JSON passed to `trnagraph graph`, at the point it's read
    in adataGraph.py -- previously a bare `json.load()` with no structural checks, which let a
    `var_r`-without-`var` config crash later as an uncaught `KeyError` on `config['var']`
    instead of a clear, immediate validation error.
    '''
    model_config = ConfigDict(extra='forbid')

    name: str
    obs: Optional[Dict[str, List[Any]]] = None
    obs_r: Optional[Dict[str, List[Any]]] = None
    var: Optional[Dict[str, List[Any]]] = None
    var_r: Optional[Dict[str, List[Any]]] = None

    @field_validator('name')
    @classmethod
    def _validate_name(cls, v: str) -> str:
        if not v.strip():
            raise ValueError('Config "name" must not be empty.')
        # "name" is concatenated directly into the output path (self.args.output += '/' + name)
        # -- reject path separators/traversal rather than letting a config write outside output.
        if '/' in v or '\\' in v or '..' in v:
            raise ValueError(f'Config "name" must not contain path separators or "..": got "{v}".')
        return v


class ColormapFile(RootModel[Dict[str, Dict[str, str]]]):
    '''
    Validates the `--colormap` JSON passed to `trnagraph graph`/`preprocess trim`: a top-level
    dict keyed by grouping-column name (e.g. 'group', or 'trimtype' for the trim-stats plot),
    each mapping category value -> color string. Previously a bare `json.load()` -- a
    non-string color value (or a flat, non-nested dict) would only surface much later as an
    opaque matplotlib error inside `mplcolors.to_rgb()`.
    '''
    root: Dict[str, Dict[str, str]]


class MetadataFile(BaseModel):
    '''
    Validates the metadata/samples file passed to `analyze build -i`, after the fastq column
    has been dropped and any 'sample'/'group' header row separated out, at the point it's
    read in adataBuild.py -- previously just a column-count check, so a duplicate sample name
    silently produced duplicate `adata.obs` rows instead of failing here with a clear error.
    '''
    path: str
    header: List[str]
    rows: List[List[Optional[str]]]

    @field_validator('header')
    @classmethod
    def _check_min_columns(cls, v: List[str]) -> List[str]:
        if len(v) < 2:
            raise ValueError(f'Metadata file must have at least 2 columns after the fastq column (sample, group): got {len(v)}.')
        return v

    @model_validator(mode='after')
    def _check_rows(self) -> 'MetadataFile':
        if not self.rows:
            raise ValueError(f'Metadata file {self.path} has no data rows.')
        bad_lengths = {len(r) for r in self.rows if len(r) != len(self.header)}
        if bad_lengths:
            raise ValueError(
                f'Metadata file {self.path} has row(s) with {sorted(bad_lengths)} column(s), '
                f'expected {len(self.header)} (matching the header).'
            )
        samples = [r[0] for r in self.rows]
        if any(not s for s in samples):
            raise ValueError(f'Metadata file {self.path} has one or more rows with an empty sample name.')
        dupes = sorted({s for s in samples if samples.count(s) > 1})
        if dupes:
            raise ValueError(f'Metadata file {self.path} has duplicate sample name(s): {dupes}.')
        return self


class PairsFile(BaseModel):
    '''
    Validates the `--pairs` file (whitespace-separated Sample1/Sample2 rows used for build-time
    DESeq2 pairwise comparisons) at read time, in adataBuild.py.run_pairwise_de -- previously a
    read failure was caught and swallowed into an empty result with only a printed message, so
    a mistyped `--pairs` path silently produced a build with zero pairwise-DE results.
    '''
    path: str
    pairs: List[Tuple[str, str]]

    @model_validator(mode='after')
    def _check_non_empty(self) -> 'PairsFile':
        if not self.pairs:
            raise ValueError(f'Pairs file {self.path} contains no valid Sample1/Sample2 pairs.')
        return self
