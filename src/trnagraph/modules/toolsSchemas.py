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
    # Read-level mismatch histogram, raw counts. None when the variant's build produced no
    # mismatch count file. A split variant carries tRNA rows only, since non-tRNA features
    # are excluded from read-length splits entirely.
    mismatch_counts: Optional[pd.DataFrame] = None


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


# Graph types that accept their own style block. Kept here rather than imported from
# adataGraph to avoid a circular import; the test suite pins the two lists together.
STYLE_GRAPH_TYPES = ('cluster', 'compare', 'correlation', 'count', 'coverage', 'heatmap',
                     'logo', 'mismatch', 'pca', 'radar', 'volcano', 'trimming')

OUTPUT_FORMATS = ('pdf', 'svg', 'png')

# Keys every graph type can honour.
UNIVERSAL_STYLE_KEYS = frozenset({'format', 'dpi', 'font_size', 'figsize'})

# Extra keys each graph type consumes, beyond the universal set. A key set in a graph type's
# own block but absent here is rejected at load rather than silently ignored -- a style file
# that appears to do nothing is the failure mode this whole schema exists to avoid. The
# `defaults` block is exempt: it broadcasts, so a type that cannot use a key just skips it.
GRAPH_STYLE_SUPPORT = {
    'volcano':     frozenset({'marker_size', 'alpha', 'rasterize_over'}),
    'pca':         frozenset({'marker_size', 'alpha', 'rasterize_over'}),
    'cluster':     frozenset({'marker_size', 'alpha', 'rasterize_over'}),
    'mismatch':    frozenset({'marker_size', 'alpha', 'rasterize_over'}),
    # coverage/radar/logo deliberately do NOT accept alpha: their alpha values are structural
    # (shaded arm bands, fill translucency scaled by how many series overlay each other) rather
    # than a point-opacity knob, so a global override would break tuned visuals rather than
    # restyle them.
    'coverage':    frozenset(),
    'radar':       frozenset(),
    'logo':        frozenset(),
    'compare':     frozenset(),
    'correlation': frozenset(),
    'count':       frozenset(),
    'heatmap':     frozenset(),
    'trimming':    frozenset(),
}


class StyleBlock(BaseModel):
    '''
    Presentation settings, valid either as the file's `defaults` or as one graph type's
    override. `extra='forbid'` so a typo like "marker_sizes" fails at load instead of being
    silently ignored while the user wonders why the figure did not change -- the same
    reasoning that put extra='forbid' on GraphFilterConfig.
    '''
    model_config = ConfigDict(extra='forbid')

    # Applied to individual plots only. Combined/multi-page pages compute their own page
    # geometry from how many panels they are laying out, so a fixed figsize there would
    # either clip panels or leave dead space; those ignore it with a warning.
    figsize: Optional[Tuple[float, float]] = None
    marker_size: Optional[float] = None
    font_size: Optional[float] = None
    dpi: Optional[int] = None
    alpha: Optional[float] = None
    # Above this many points, a scatter layer is rasterized while text and axes stay vector.
    # A vector PDF carrying tens of thousands of points is slow to open and is often rejected
    # by journal submission systems.
    rasterize_over: Optional[int] = None
    format: Optional[Literal['pdf', 'svg', 'png']] = None

    @field_validator('figsize')
    @classmethod
    def _positive_figsize(cls, v):
        if v is not None and (v[0] <= 0 or v[1] <= 0):
            raise ValueError(f'figsize must be two positive numbers, got {v}.')
        return v

    @field_validator('marker_size', 'font_size', 'dpi')
    @classmethod
    def _positive(cls, v, info):
        if v is not None and v <= 0:
            raise ValueError(f'{info.field_name} must be positive, got {v}.')
        return v

    @field_validator('alpha')
    @classmethod
    def _unit_interval(cls, v):
        if v is not None and not 0 <= v <= 1:
            raise ValueError(f'alpha must be between 0 and 1, got {v}.')
        return v


class StyleFile(BaseModel):
    '''
    Validates the `--style` JSON passed to `trnagraph graph` and `preprocess trim`.

    Supersedes the former `--colormap` file, so that one file carries both the palette and
    the presentation settings for a figure set rather than the user juggling two. A file in
    the old colormap shape (a bare mapping of column -> value -> color) is still accepted and
    read as `colors`, which is what keeps existing colormap.json files working.
    '''
    model_config = ConfigDict(extra='forbid')

    colors: Optional[Dict[str, Dict[str, str]]] = None
    defaults: Optional[StyleBlock] = None
    cluster: Optional[StyleBlock] = None
    compare: Optional[StyleBlock] = None
    correlation: Optional[StyleBlock] = None
    count: Optional[StyleBlock] = None
    coverage: Optional[StyleBlock] = None
    heatmap: Optional[StyleBlock] = None
    logo: Optional[StyleBlock] = None
    mismatch: Optional[StyleBlock] = None
    pca: Optional[StyleBlock] = None
    radar: Optional[StyleBlock] = None
    volcano: Optional[StyleBlock] = None
    trimming: Optional[StyleBlock] = None

    @model_validator(mode='before')
    @classmethod
    def _accept_bare_colormap(cls, data):
        '''
        A legacy colormap file has no `colors`/`defaults`/graph-type keys at all -- every top
        level key is a grouping column whose value is a color mapping. Detect that shape and
        lift it into `colors` rather than failing on extra='forbid', so an existing
        colormap.json keeps working under --style unchanged.
        '''
        if not isinstance(data, dict) or not data:
            return data
        known = {'colors', 'defaults', *STYLE_GRAPH_TYPES}
        if any(key in known for key in data):
            return data
        if all(isinstance(v, dict) and all(isinstance(c, str) for c in v.values())
               for v in data.values()):
            return {'colors': data}
        return data

    @model_validator(mode='after')
    def _reject_unsupported_per_type_keys(self) -> 'StyleFile':
        for graph_type in STYLE_GRAPH_TYPES:
            block = getattr(self, graph_type, None)
            if block is None:
                continue
            allowed = UNIVERSAL_STYLE_KEYS | GRAPH_STYLE_SUPPORT.get(graph_type, frozenset())
            unsupported = sorted(set(block.model_dump(exclude_none=True)) - allowed)
            if unsupported:
                raise ValueError(
                    f"'{graph_type}' does not use {unsupported}. Supported for this graph "
                    f"type: {sorted(allowed)}. Put it in 'defaults' if you meant it to apply "
                    f"wherever it is relevant."
                )
        return self

    def resolve(self, graph_type: str) -> Dict[str, Any]:
        '''
        Settings for one graph type: `defaults` with that type's own block laid over it.
        Only keys the file actually set appear, so a caller can distinguish "not configured"
        from "configured to the same value as the built-in default".
        '''
        merged: Dict[str, Any] = {}
        if self.defaults is not None:
            merged.update(self.defaults.model_dump(exclude_none=True))
        override = getattr(self, graph_type, None) if graph_type in STYLE_GRAPH_TYPES else None
        if override is not None:
            merged.update(override.model_dump(exclude_none=True))
        return merged

    def colors_for(self, column: Optional[str]) -> Optional[Dict[str, str]]:
        '''The palette for one grouping column, or None when the file does not set it.'''
        if not self.colors or column is None:
            return None
        return self.colors.get(column)
