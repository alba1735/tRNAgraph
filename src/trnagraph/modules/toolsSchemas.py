'''
Pydantic data models shared across tRNAgraph's variant-selection and size-split-merge
machinery (see toolsTG.py's parse_variant()/build_variant_view() and adataBuild.py's
compute_variant_contribution()/merge_variant_into_adata()). Kept in a dedicated module,
following the tools*.py convention, so these schemas have one canonical, validated definition
instead of being duplicated as plain dataclasses in the modules that use them.
'''

import re
from typing import Any, Dict, List, Optional, Literal, Tuple, Union
import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict, field_validator, model_validator

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


class VennSet(BaseModel):
    '''
    One circle of a Venn diagram: which population of features it contains.

    A circle is identified by where its counts come from -- a read-length variant (`tag`) and a
    read type -- and optionally by one level of a grouping column, which is what lets a complex
    Venn cross timepoint with variant. `label` is what the reader sees on the figure.
    '''
    model_config = ConfigDict(extra='forbid')

    label: str
    readtype: str
    tag: str = 'full'
    level: Optional[str] = None


class VennPlan(BaseModel):
    '''One Venn diagram to draw: its output name, its title, and its circles.'''
    model_config = ConfigDict(extra='forbid')

    name: str
    title: str
    sets: List[VennSet]


class MultivariateConfig(BaseModel):
    '''
    The `multivariate` block: which analysis the set-membership plots describe.

    Top level rather than a `flags.graph` key, for the reason the filters and `order` are: it
    says what the experiment IS. It is also what GATES `-g venn` and `-g agreement`, which are
    excluded from `-g all` -- the sets and thresholds are choices about a specific design, and
    figures made from unchosen ones invite wrong conclusions.

    Thresholds default to the project-wide pair (see plotsThresholds). Membership uses the
    stricter of the two, since a set-overlap claim is stronger than a single volcano's.
    '''
    model_config = ConfigDict(extra='forbid')

    #: The obs column whose levels the analysis is taken over.
    grouping: str = 'group'
    #: The level every contrast is measured against. Defaults to the first category of an
    #: ordered `order` declaration, and to the first level present otherwise.
    reference: Optional[str] = None
    #: Minimum mean normalized count, per group, for a feature to count as PRESENT.
    presence_cutoff: float = 20.0
    #: Significance thresholds for the DE-hit membership mode.
    log2fc: float = 1.5
    padj: float = 0.001


class RunConfig(BaseModel):
    '''
    Validates the `--config` JSON, which describes one saved analysis: which subset of the
    data, and which options each command runs with.

    The top level holds what is global to the file -- its `name`, and the obs/var filters that
    say which subset of the object every AnnData-consuming command sees. Per-command options
    live under `flags`, one block per command. Filters are not CLI flags and are not scoped to
    a single command, which is why they stay at the top level rather than moving inside a block.
    '''
    model_config = ConfigDict(extra='forbid')

    name: str
    obs: Optional[Dict[str, List[Any]]] = None
    obs_r: Optional[Dict[str, List[Any]]] = None
    var: Optional[Dict[str, List[Any]]] = None
    var_r: Optional[Dict[str, List[Any]]] = None
    # Explicit category order per obs column, e.g. {"timepoint": ["Day 0", "Day 35", "Day 70"]}.
    # Top level rather than a `flags.build` key for the same reason the filters are: it says what
    # the experiment IS, not how one command ran. It is also the single source of truth for two
    # different consumers -- plot legend/axis order, and the DESeq2 reference level (the FIRST
    # listed level) -- which a per-command flag could not be, since it would have to be
    # re-supplied identically on every plotting and DE call to stay consistent.
    order: Optional[Dict[str, List[str]]] = None
    # Declares the multivariate analysis, and gates `-g venn`/`-g agreement`. Top level for the
    # same reason `order` is: it describes the experiment, not one command's invocation.
    multivariate: Optional['MultivariateConfig'] = None
    # Per-command options this config pins, one block per command. They live here rather than
    # in --style because these are selection/analysis choices, not presentation: a style file
    # is meant to be shared across differently-parameterized runs, which it could not be if it
    # also fixed the cutoffs.
    flags: Optional['CommandFlags'] = None

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


# Options NO command's flags block may set. What is excluded is where a run writes and how it
# is driven -- an output directory, a thread count, verbosity -- plus the config and style
# files themselves, which a file able to name them could redirect. `format` is excluded
# because --style already owns it, which leaves exactly zero keys settable from two files and
# so no precedence rule to document between them.
#
# INPUT paths are deliberately settable. For `build`, `map` and `makedb` the reference, GTF,
# bam directory and metadata ARE the analysis; withholding them would leave those blocks
# nearly empty and a config file could not describe a run end to end, which is the point of
# having them. The line is between what an analysis IS and where its output goes.
COMMAND_FLAG_EXCLUSIONS = frozenset({
    'anndata', 'anndata_path', 'output', 'config', 'style', 'format',
    'threads', 'quiet', 'verbose',
})

#: Kept as the old name for `graph`'s own exclusions, which are the same set.
GRAPH_FLAG_EXCLUSIONS = COMMAND_FLAG_EXCLUSIONS

# Flags whose value is a list. Declared once so the empty-list guard and the tests that pin
# replace-not-append semantics agree on the set.
GRAPH_LIST_FLAGS = ('graphtypes', 'diffrts', 'pcareadtypes', 'radarmethod', 'logomanualgrp')


class GraphFlags(BaseModel):
    '''
    `graph` options settable from a --config file, so one file can carry a whole saved
    analysis -- which subset of the data, under which grouping columns and cutoffs -- instead
    of a shell line that has to be retyped correctly every time.

    Every field is Optional and defaults to None, which is what makes "not configured"
    distinguishable from "configured to the same value as the CLI default". A CLI flag the
    user actually typed always wins; see adataGraph's merge and cli.py's ParameterSource
    capture for how that is detected.

    Written out by hand rather than derived from the typer command, because importing cli.py
    here would invert the package's dependency direction (cli -> modules) and defeat
    lazy_imports by dragging the whole CLI in whenever a schema is touched. A unit test
    asserts this model's fields equal the graph command's eligible parameters, so the two
    cannot drift apart.
    '''
    model_config = ConfigDict(extra='forbid')

    # Whole-command selection.
    graphtypes: Optional[List[str]] = None
    variant: Optional[str] = None
    allreads: Optional[bool] = None
    regen_uns: Optional[bool] = None

    # cluster
    clustergrp: Optional[str] = None
    clusterlabels: Optional[str] = None
    clusteroverview: Optional[bool] = None
    clusternumeric: Optional[bool] = None
    clustermask: Optional[bool] = None

    # compare
    comparegrp1: Optional[str] = None
    comparegrp2: Optional[str] = None

    # correlation
    corrmethod: Optional[str] = None
    corrgroup: Optional[str] = None
    corrmask: Optional[Literal['none', 'upper', 'lower']] = None

    # coverage
    covgrp: Optional[str] = None
    covobs: Optional[str] = None
    covtype: Optional[str] = None
    covgap: Optional[bool] = None
    covmethod: Optional[str] = None
    combinedpdfonly: Optional[bool] = None

    # heatmap (diffrts is shared with volcano)
    heatgrp: Optional[str] = None
    diffrts: Optional[List[str]] = None
    heatcutoff: Optional[int] = None
    heatbound: Optional[int] = None
    heatsubplots: Optional[bool] = None
    heatorient: Optional[Literal['vertical', 'horizontal']] = None

    # pca
    pcamarkers: Optional[str] = None
    pcacolors: Optional[str] = None
    pcareadtypes: Optional[List[str]] = None

    # radar
    radargrp: Optional[str] = None
    radarmethod: Optional[List[str]] = None
    radarscaled: Optional[bool] = None

    # logo
    logogrp: Optional[str] = None
    logomanualgrp: Optional[List[str]] = None
    logomanualname: Optional[str] = None
    logopseudocount: Optional[int] = None
    logosize: Optional[str] = None
    ccatail: Optional[bool] = None
    pseudogenes: Optional[bool] = None
    logornamode: Optional[bool] = None

    # mismatch
    mismatchpseudocount: Optional[int] = None

    # volcano
    volgrp: Optional[str] = None
    volcutoff: Optional[int] = None
    shrink: Optional[Literal['apeGLM', 'none']] = None
    volxlim: Optional[float] = None
    vollabels: Optional[int] = None

    @field_validator(*GRAPH_LIST_FLAGS)
    @classmethod
    def _reject_empty_list(cls, v, info):
        '''
        A list REPLACES the default rather than extending it -- replacing is the only
        semantics that lets a config narrow a list, which is the point. That makes an empty
        list a request to draw nothing, so it is rejected: `"graphtypes": []` would otherwise
        produce a run that renders no figures and reports success.
        '''
        if v is not None and len(v) == 0:
            raise ValueError(
                f"'{info.field_name}' is empty. A list here replaces the default rather than "
                f"adding to it, so an empty list would select nothing at all. Remove the key "
                f"to keep the default."
            )
        return v


class MakedbFlags(BaseModel):
    '''`preprocess makedb` options settable from a --config file.'''
    model_config = ConfigDict(extra='forbid')

    genome: Optional[str] = None
    trnaout: Optional[str] = None
    trnafa: Optional[str] = None
    namemap: Optional[str] = None
    addtrna: Optional[str] = None
    addseqs: Optional[str] = None
    orgmode: Optional[str] = None
    forcecca: Optional[bool] = None


class TrimFlags(BaseModel):
    '''`preprocess trim` options settable from a --config file.'''
    model_config = ConfigDict(extra='forbid')

    input: Optional[str] = None
    adapter1: Optional[str] = None
    adapter2: Optional[str] = None
    length: Optional[int] = None
    umilength: Optional[int] = None
    umi3: Optional[bool] = None


class MapFlags(BaseModel):
    '''`preprocess map` options settable from a --config file.'''
    model_config = ConfigDict(extra='forbid')

    database: Optional[str] = None
    input: Optional[str] = None
    force_remap: Optional[bool] = None
    minnontrnasize: Optional[int] = None
    local: Optional[bool] = None
    skipcheck: Optional[bool] = None
    bamdir: Optional[str] = None
    dedup: Optional[bool] = None
    keep_prededup: Optional[bool] = None
    dedup_method: Optional[str] = None


class BuildFlags(BaseModel):
    '''`analyze build` options settable from a --config file.'''
    model_config = ConfigDict(extra='forbid')

    input: Optional[str] = None
    database: Optional[str] = None
    gtf: Optional[str] = None
    pairs: Optional[str] = None
    bed: Optional[List[str]] = None
    maxmismatches: Optional[str] = None
    minfeaturereads: Optional[str] = None
    minnontrnasize: Optional[int] = None
    hub: Optional[bool] = None
    hubonly: Optional[bool] = None
    filterother: Optional[bool] = None
    bamdir: Optional[str] = None
    dispfittype: Optional[str] = None
    readlengthsplit: Optional[int] = None
    overwritebams: Optional[bool] = None
    savesplitbams: Optional[bool] = None
    vst: Optional[str] = None


class AddsplitFlags(BaseModel):
    '''`analyze addsplit` options settable from a --config file.'''
    model_config = ConfigDict(extra='forbid')

    readlengthsplit: Optional[int] = None
    metadata: Optional[str] = None
    bamdir: Optional[str] = None
    database: Optional[str] = None
    gtf: Optional[str] = None
    dispfittype: Optional[str] = None
    vst: Optional[str] = None
    minfeaturereads: Optional[str] = None
    overwritebams: Optional[bool] = None
    savesplitbams: Optional[bool] = None
    overwrite: Optional[bool] = None
    force: Optional[bool] = None


class ClusterFlags(BaseModel):
    '''`analyze cluster` options settable from a --config file.'''
    model_config = ConfigDict(extra='forbid')

    randomstate: Optional[int] = None
    readcutoff: Optional[int] = None
    coveragetype: Optional[List[str]] = None
    ncomponentsmp: Optional[int] = None
    ncomponentgrp: Optional[int] = None
    neighborclusmp: Optional[int] = None
    neighborclusgrp: Optional[int] = None
    neighborstdsmp: Optional[int] = None
    neighborstdgrp: Optional[int] = None
    hdbscanminsampsmp: Optional[int] = None
    hdbscanminsampgrp: Optional[int] = None
    hdbscanminclusmp: Optional[int] = None
    hdbscanminclugrp: Optional[int] = None
    mindist: Optional[float] = None
    variancethreshold: Optional[float] = None
    umapstatsmetrics: Optional[str] = None
    hdbstatsmetrics: Optional[str] = None
    clusterobsexperimental: Optional[List[str]] = None
    variant: Optional[str] = None
    overwrite: Optional[bool] = None


class Log2fcFlags(BaseModel):
    '''`tools log2fc` options settable from a --config file.'''
    model_config = ConfigDict(extra='forbid')

    group: Optional[str] = None
    readtypes: Optional[List[str]] = None
    cutoff: Optional[List[int]] = None
    variant: Optional[str] = None


class CommandFlags(BaseModel):
    '''
    One block of pinned options per command, so a single --config file can carry a whole run
    rather than each stage's settings living in a shell script.

    Blocks exist for the commands whose options are worth saving between runs. `tools csv`,
    `merge`, `info`, `test` and `template` have none: their options are paths, output
    destinations and one-shot selectors, with nothing an analysis would want to fix.
    '''
    model_config = ConfigDict(extra='forbid')

    graph: Optional['GraphFlags'] = None
    build: Optional[BuildFlags] = None
    map: Optional[MapFlags] = None
    trim: Optional[TrimFlags] = None
    makedb: Optional[MakedbFlags] = None
    cluster: Optional[ClusterFlags] = None
    addsplit: Optional[AddsplitFlags] = None
    log2fc: Optional[Log2fcFlags] = None


#: {command name: its flags model}. Drives the drift tests, the template generator and the
#: per-command merge, so a new block is added in exactly one place.
COMMAND_FLAG_MODELS = {
    'graph': None,  # bound below, once GraphFlags is defined
    'build': BuildFlags,
    'map': MapFlags,
    'trim': TrimFlags,
    'makedb': MakedbFlags,
    'cluster': ClusterFlags,
    'addsplit': AddsplitFlags,
    'log2fc': Log2fcFlags,
}


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
# `line_width` goes to the types that draw a data trace or a bar edge. The scatter types are
# served by marker_size instead, and heatmap/correlation/logo have no stroke a user would set.
GRAPH_STYLE_SUPPORT = {
    'volcano':     frozenset({'marker_size', 'alpha', 'rasterize_over'}),
    'pca':         frozenset({'marker_size', 'alpha', 'rasterize_over'}),
    'cluster':     frozenset({'marker_size', 'alpha', 'rasterize_over'}),
    'mismatch':    frozenset({'marker_size', 'alpha', 'rasterize_over', 'line_width'}),
    # coverage/radar/logo deliberately do NOT accept alpha: their alpha values are structural
    # (shaded arm bands, fill translucency scaled by how many series overlay each other) rather
    # than a point-opacity knob, so a global override would break tuned visuals rather than
    # restyle them.
    'coverage':    frozenset({'line_width'}),
    'radar':       frozenset({'line_width'}),
    'logo':        frozenset(),
    'compare':     frozenset({'line_width'}),
    'correlation': frozenset(),
    'count':       frozenset({'line_width'}),
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
    # Stroke width for data traces and bar edges, in points, replacing the module's own tuned
    # default. Shrinking figsize scales geometry but not strokes, so a small figure otherwise
    # renders with lines far too heavy for the panel they sit in.
    line_width: Optional[float] = None
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

    @field_validator('marker_size', 'line_width', 'font_size', 'dpi')
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


# A palette value in a style file is either a registered colormap/palette NAME or an explicit
# list of color tokens. One shape for both the ordered ramps and the categorical fallback, so
# a user learns it once.
PaletteValue = Union[str, List[str]]

# Ordered/continuous scales a style file can set, named by what they ENCODE rather than by
# appearance, matching plotsPalette's role naming. Deliberately not keyed by graph type: the
# heatmap draws two of these at once (`lfc` and `significance`), the seqlogo draws two more
# (`score` and `sequence`), and `ordered` is shared by coverage and cluster -- so a per-graph
# type `cmap` key could not express what the plots actually do.
GRADIENT_ROLES = ('correlation', 'significance', 'score', 'sequence', 'ordered', 'lfc')


def _check_colors(values, field):
    '''Every entry must be something matplotlib can interpret as a color.'''
    from matplotlib.colors import is_color_like

    bad = [c for c in values if not (isinstance(c, str) and is_color_like(c))]
    if bad:
        raise ValueError(
            f"'{field}' contains {bad}, which matplotlib does not recognize as colors. Use "
            f"hex ('#1f77b4'), a named color ('tab:blue', 'red'), or any other matplotlib "
            f"color token."
        )


def _validate_colormap_name(name, field):
    '''
    Reject an unknown colormap at FILE LOAD rather than at draw time.

    seaborn is imported first because mako/rocket/crest/flare/vlag/icefire are registered by
    seaborn, not shipped with matplotlib -- checking matplotlib alone would reject the very
    ramps tRNAgraph itself defaults to.
    '''
    import seaborn as sns  # noqa: F401  -- registers seaborn's colormaps as a side effect
    import matplotlib.pyplot as plt

    try:
        plt.get_cmap(name)
    except (ValueError, KeyError):
        # Colormap names are case-sensitive ('Blues', not 'blues'), which is the easiest way to
        # get this wrong, so point at the near match rather than leaving the user to guess.
        import difflib

        registered = list(plt.colormaps)
        # A case-only difference is the common mistake, and difflib compares case-sensitively,
        # so check that exactly before falling back to fuzzy matching.
        folded = {c.lower(): c for c in registered}
        match = folded.get(name.lower())
        if match is None:
            close = difflib.get_close_matches(name.lower(), list(folded), n=1, cutoff=0.7)
            match = folded[close[0]] if close else None
        suggestion = f" Did you mean '{match}'?" if match else ''
        raise ValueError(
            f"'{field}' names colormap '{name}', which is not registered.{suggestion} Names are "
            f"case-sensitive. Use a matplotlib or seaborn colormap name (e.g. 'mako_r', 'vlag', "
            f"'Blues'), or give a list of two or more colors to build your own ramp."
        ) from None


def _validate_gradient(v, field):
    if v is None:
        return v
    if isinstance(v, str):
        _validate_colormap_name(v, field)
        return v
    if len(v) < 2:
        raise ValueError(
            f"'{field}' needs at least two colors to interpolate a ramp, got {len(v)}."
        )
    _check_colors(v, field)
    return v


class GradientBlock(BaseModel):
    '''
    The ordered/continuous scales, each settable to a colormap name or a list of colors.

    extra='forbid' so a misspelled role (`significants`) fails at load rather than leaving the
    user staring at an unchanged figure -- the same reasoning behind StyleBlock's forbid.
    '''
    model_config = ConfigDict(extra='forbid')

    # correlation R^2 heatmap.
    correlation: Optional[PaletteValue] = None
    # -log10(p) panel; drawn beside `lfc`, so the two must stay mutually distinguishable.
    significance: Optional[PaletteValue] = None
    # seqlogo per-position score heatmap and its colorbar.
    score: Optional[PaletteValue] = None
    # seqlogo sequence heatmap background.
    sequence: Optional[PaletteValue] = None
    # coverage specificity partition and cluster numeric coloring: light = least, dark = most.
    ordered: Optional[PaletteValue] = None
    # diverging scale for log2 fold change, centered on zero -- wants an odd number of stops
    # with a neutral middle when given as a list.
    lfc: Optional[PaletteValue] = None

    @field_validator(*GRADIENT_ROLES)
    @classmethod
    def _valid_gradient(cls, v, info):
        return _validate_gradient(v, info.field_name)


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
    # Ordered/continuous scales, keyed by the quantity they encode. A sibling of `colors`
    # rather than a per-graph-type key, because the roles are shared across graph types.
    gradients: Optional[GradientBlock] = None
    # The fallback palette for unordered categories a `colors` entry does not name. Its own
    # key rather than a `gradients` role: it is a seaborn palette namespace, not a colormap
    # one, so 'husl' is valid here and meaningless there.
    categorical: Optional[PaletteValue] = None
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

    @field_validator('categorical')
    @classmethod
    def _valid_categorical(cls, v):
        '''
        A seaborn palette NAME or an explicit list of colors.

        Validated separately from the gradients because the namespaces differ: 'husl' is a
        legitimate categorical palette and not a colormap, while 'mako_r' is the reverse.
        '''
        if v is None or not isinstance(v, str):
            if v is not None:
                if not v:
                    raise ValueError("'categorical' must name at least one color.")
                _check_colors(v, 'categorical')
            return v
        import seaborn as sns

        try:
            sns.color_palette(v, 2)
        except Exception:
            raise ValueError(
                f"'categorical' names palette '{v}', which seaborn does not recognize. Use a "
                f"seaborn palette name (e.g. 'colorblind', 'husl', 'tab10'), or give an "
                f"explicit list of colors."
            ) from None
        return v

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
        known = {'colors', 'gradients', 'categorical', 'defaults', *STYLE_GRAPH_TYPES}
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

def _blocks_owning(key: str):
    """Which command blocks declare `key`, in the order they appear in the file."""
    return [command for command, model in COMMAND_FLAG_MODELS.items()
            if model is not None and key in model.model_fields]


def explain_rejected_keys(exc, file_kind: str):
    '''
    Turn pydantic's `extra_forbidden` reports into a sentence naming where the key belongs.

    Three mistakes are worth naming, and none of them is a typo. A key can be in the wrong
    FILE (`--config` and `--style` deliberately share no keys, which is what lets a style file
    be reused across differently-parameterized runs). It can be in the wrong BLOCK, an easy
    slip once `flags` holds eight of them. Or it can be at the wrong DEPTH: `flags` used to be
    graph's options directly, so every config file written before the per-command blocks has
    its keys one level too shallow. Pydantic reports all three as "Extra inputs are not
    permitted", which is true and useless; only a typo falls through to a spelling guess.

    `file_kind` is 'style' or 'config'. Returns a list of lines to append to the error.
    '''
    import difflib

    style_keys = set(StyleBlock.model_fields)
    lines = []
    for error in exc.errors():
        if error.get('type') != 'extra_forbidden':
            continue
        loc = [str(part) for part in error['loc']]
        key = loc[-1]
        # Where the key was written: the block it sat in, or None for a bare `flags` key.
        placed_in = loc[1] if len(loc) >= 3 and loc[0] == 'flags' else None

        if file_kind == 'config' and key in style_keys:
            lines.append(f"  '{key}' is a presentation setting, not a command option: put it "
                         f"in your --style file instead.")
            continue

        owners = _blocks_owning(key)
        if owners and (file_kind == 'style' or placed_in not in owners):
            where = ' or '.join(f'`flags.{command}`' for command in owners)
            if file_kind == 'style':
                lines.append(f"  '{key}' is a {owners[0]} option, not a presentation setting: "
                             f"put it under {where} in your --config file instead.")
            elif placed_in is None:
                # The commonest case by far: a file written when `flags` WAS graph's options.
                lines.append(f"  '{key}' now belongs one level deeper, under {where}: `flags` "
                             f"takes one block per command rather than graph's options directly.")
            else:
                lines.append(f"  '{key}' is not a `{placed_in}` option; it belongs under {where}.")
            continue

        # Nothing structural -- fall back to the nearest valid spelling at this location.
        if file_kind == 'style':
            valid = style_keys
        elif placed_in in COMMAND_FLAG_MODELS and COMMAND_FLAG_MODELS[placed_in] is not None:
            valid = set(COMMAND_FLAG_MODELS[placed_in].model_fields)
        elif loc[:1] == ['flags']:
            valid = set(CommandFlags.model_fields)
        else:
            valid = set(RunConfig.model_fields)
        close = difflib.get_close_matches(key, sorted(valid), n=3, cutoff=0.7)
        if close:
            lines.append(f"  '{key}' is not a valid key; did you mean: {', '.join(close)}?")
    return lines


#: The previous name for RunConfig, from when the file carried only `graph`'s options.
GraphFilterConfig = RunConfig

COMMAND_FLAG_MODELS['graph'] = GraphFlags

# RunConfig is declared before CommandFlags, so its forward reference is resolved here.
RunConfig.model_rebuild()
