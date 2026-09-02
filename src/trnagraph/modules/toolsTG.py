#!/usr/bin/env python3

import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats
import contextlib
import importlib.resources
import itertools
import os
import sys
import logging
import gzip
import re
import tempfile
import subprocess
import pysam
from shutil import which
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from typing import Generator, Iterable, Tuple, Union, TextIO, Dict, Optional, List, Any, Callable
from collections import defaultdict, deque
from pathlib import Path
from dataclasses import dataclass, field

from rich.live import Live
from rich.progress import Progress, BarColumn, TextColumn, TaskProgressColumn
from rich.spinner import Spinner
from rich.console import Group, Console
from rich.text import Text

logger = logging.getLogger(__name__)


def assets_dir() -> str:
    '''
    Absolute path to the package's bundled `assets/` directory -- the Sprinzl covariance models
    under `cm/`, the vibrChol1 demo manifest/metadata/pair files, and the default colormap.

    Resolved through importlib.resources rather than by counting directories up from `__file__`,
    because `assets/` is a *package* subdirectory shipped by [tool.setuptools.package-data], not
    a source-tree sibling. The directory walk this replaces assumed a src-layout checkout and so
    landed on `<venv>/lib/pythonX.Y` under a non-editable `pip install .`, where nothing exists.

    Returns a plain `str`: callers interpolate it into shell commands (toolsTestSuite's
    `cp {assets_dir}/*.txt config/.`), and importlib.resources yields a real filesystem path for
    any non-zip install -- which this package is by construction, since it shells out to
    bowtie2/samtools/infernal against real files on disk.
    '''
    return os.fspath(importlib.resources.files('trnagraph') / 'assets')


def builder(directory: Union[str, Path]) -> str:
    '''
    Function to create output directory if it does not already exist
    '''
    dir_path = Path(directory)
    # Fix for potential relative path issues if needed, though Path handles most
    if not dir_path.is_absolute() and not str(dir_path).startswith('.'):
         dir_path = Path('./') / dir_path

    if not dir_path.exists():
        output = f'Creating output directory: {dir_path}'
        dir_path.mkdir(parents=True, exist_ok=True)
    else:
        output = f'Output directory already exists: {dir_path}'
    return output


class UsageError(ValueError):
    '''
    The command as typed cannot work, for a reason the user can fix by retyping it.

    A dedicated hierarchy rather than ValueError because cli.py renders these as a one-line
    usage error and exits, the way toolsTestSuite.WorkspaceNotOwnedError is already handled --
    a mistake in a flag is a usage problem, and a traceback buries the one sentence that tells
    the user what to type instead. usage_error_guard() catches this base, so a new kind of
    usage mistake is rendered correctly by subclassing rather than by editing the CLI.

    Subclasses ValueError so that callers (and tests) written before this hierarchy existed,
    which catch ValueError from the config/style loaders, keep working unchanged.
    '''


class UnknownLabelError(UsageError):
    '''A parameter names something the AnnData object does not contain.'''


class InvalidParameterError(UsageError):
    '''
    Parameters that each name something real, in a combination that cannot work.

    Distinct from UnknownLabelError because the failure is different in kind: nothing is
    missing, so there is no near match to suggest and no vocabulary to list.
    '''


def _near_matches(value: str, candidates, limit: int = 3):
    '''
    The candidates a mistyped label most plausibly meant, best first.

    Case-only differences are checked before fuzzy matching, exactly as the --style colormap
    validator does: difflib compares case-sensitively, so 'Group' against 'group' would
    otherwise score as a near miss rather than the certain match it is.
    '''
    import difflib

    candidates = list(candidates)
    folded = {str(c).lower(): c for c in candidates}
    exact_fold = folded.get(str(value).lower())
    if exact_fold is not None:
        return [exact_fold]
    # 0.7, not difflib's 0.6 default, and the same threshold the --style colormap validator
    # uses. Measured on real labels from the vibrChol1 object, every genuine typo scores >=
    # 0.75 ('grp'/'group' 0.750, 'goup'/'group' 0.889) and every misleading pairing <= 0.667
    # ('treatement'/'fragment' 0.667) -- and a confident suggestion of an unrelated column is
    # worse than offering none.
    close = difflib.get_close_matches(str(value).lower(), list(folded), n=limit, cutoff=0.7)
    return [folded[c] for c in close]


#: Values listed when a label has no near match, before eliding. Bounded so a column with
#: hundreds of values cannot turn one error line into a wall of text.
_AVAILABLE_CAP = 12


def _label_candidates(adata: ad.AnnData, domain: str):
    '''
    The valid values for one kind of label, sorted.

    A label is not always an obs column: --covtype names an adata.var['coverage'] value, so
    checking it against obs columns would reject every valid coverage type and suggest
    nonsense for the invalid ones.
    '''
    if domain == 'obs':
        return sorted(str(c) for c in adata.obs.columns)
    if domain == 'column':
        # Any column of the object, either frame -- `tools info --column` accepts both, since
        # a user asking "what is in this column" has no reason to know which frame holds it.
        return sorted({str(c) for c in adata.obs.columns} | {str(c) for c in adata.var.columns})
    if domain == 'coverage':
        if 'coverage' not in adata.var.columns:
            return []
        return sorted({str(v) for v in adata.var['coverage']})
    if domain == 'readtype':
        # Bare readtypes, recovered from the count columns this object actually carries --
        # the same source resolve_readtype() consults, so a readtype that validates here is
        # one it can resolve. The basis component is stripped because --diffrts/--pcareadtypes
        # deliberately do not carry it; --allreads sets it once for the whole command.
        readtypes = set()
        for column in adata.obs.columns:
            match = re.fullmatch(r'nreads_(.+?)(_unique)?_norm', str(column))
            if match:
                readtypes.add(match.group(1))
        return sorted(readtypes)
    if domain == 'readtype_with_basis':
        # `tools log2fc -r` builds 'nreads_{rt}_norm' directly, so its readtypes DO carry the
        # basis ('total_unique') -- the opposite of `graph --diffrts`, where --allreads sets
        # the basis once for the command and a basis-carrying readtype is rejected outright.
        # Two vocabularies, deliberately: validating either against the other's would reject
        # every correct value.
        readtypes = set()
        for column in adata.obs.columns:
            match = re.fullmatch(r'nreads_(.+)_norm', str(column))
            if match:
                readtypes.add(match.group(1))
        return sorted(readtypes)
    raise ValueError(f"Unknown label domain '{domain}'.")


def validate_labels(adata: ad.AnnData, requests, extra_problems=()) -> None:
    '''
    Check every CLI label against the object up front, reporting all failures at once.

    `requests` is an iterable of (param_name, value, domain) triples. Batching matters: a
    command carries a dozen label-valued options, and aborting on the first bad one turns
    fixing a command line into as many round trips as there are typos.

    `extra_problems` carries already-worded problems that are not unknown labels -- parameters
    that name real things in a combination that cannot work -- so those are reported in the
    same pass rather than after the user has fixed their typos and re-run.
    '''
    problems = []
    for param_name, value, domain in requests:
        candidates = _label_candidates(adata, domain)
        if str(value) in candidates:
            continue
        problems.append((param_name, value, _near_matches(value, candidates), candidates))
    extra_problems = list(extra_problems)
    if not problems and not extra_problems:
        return
    total = len(problems) + len(extra_problems)
    if extra_problems:
        # Mixed batch: "unknown labels" would misdescribe half of it.
        lines = [f'{total} problem{"s" if total > 1 else ""} with this command:']
    else:
        lines = [f'{total} unknown label{"s" if total > 1 else ""} for this AnnData object:']
    for param_name, value, matches, candidates in problems:
        if matches:
            hint = f'  did you mean: {", ".join(matches)}?'
        else:
            # No near match leaves the message saying only "that is wrong", with nothing to
            # act on -- so fall back to the vocabulary itself, capped so a 436-value column
            # cannot turn one error line into a wall of text.
            available = candidates
            shown = ', '.join(available[:_AVAILABLE_CAP])
            if len(available) > _AVAILABLE_CAP:
                shown += f' ... +{len(available) - _AVAILABLE_CAP} more'
            hint = f'  available: {shown}' if available else '  (this object carries none)'
        lines.append(f"  --{param_name} '{value}'{hint}")
    lines.extend(f'  {problem}' for problem in extra_problems)
    lines.append('Run `trnagraph tools info -i <object.h5ad>` to list every obs/var column and its values.')
    # The type names what the batch actually contains, so a caller catching the narrower
    # UnknownLabelError is not handed something else under that name.
    error = InvalidParameterError if extra_problems else UnknownLabelError
    raise error('\n'.join(lines))


def resolve_grp_column(adata: ad.AnnData, grp: str, param_name: str) -> str:
    '''
    Validate that a user-specified grouping column exists in adata.obs, raising if it does not.

    This used to fall back to 'sample' with a warning, on the reasoning that a typo should
    degrade a command rather than abort it. That was the wrong trade: the fallback produced a
    complete, plausible set of figures grouped by a column the user never asked for, and a
    warning on stderr does not survive a redirected log or a long run. A grouping column is
    not something to guess at.

    anndataGrapher validates every label up front (see GRAPH_LABEL_PARAMS), so in the CLI
    path this is a second line of defence. It still matters for the plots*.py modules used
    standalone and for the Python API, which reach them without the grapher.
    '''
    if grp in adata.obs.columns:
        return grp
    validate_labels(adata, [(param_name, grp, 'obs')])
    return grp


def is_trna_feature(feature_name: str) -> bool:
    '''
    Whether a feature name refers to a tRNA/tRX feature rather than a non-tRNA GTF feature.

    The single definition of that distinction. It decides both which features control DESeq2
    size-factor estimation (the `use_trna_control` path in adataBuild.run_deseq2_on_file) and
    which rows survive into a read-length split variant's outputs -- keeping those two on one
    rule is what stops a split being filtered by one definition and normalized by another.

    Matches on the name rather than on genetypes.txt's type column because the name is available
    everywhere the decision has to be made, including on a counts matrix that carries no type
    annotation at all. The two agree exactly: verified as 0 disagreements across all 8343
    features of a real hg38 build, and pinned by a unit test.
    '''
    return 'tRNA' in feature_name or 'tRX' in feature_name


def resolve_nontrna_counts(adata: ad.AnnData, is_full_variant: bool, feature_label: str):
    '''
    Shared gate for plotsVolcano.py/plotsPca.py/plotsCorrelation.py's non-tRNA/combined plots.
    Returns (nontrna_df, skip_message) -- exactly one is non-None. A length-based split cutoff
    partitions tRNAs by design, but non-tRNA reads aren't classified by that criterion at all, so
    a split-variant's non-tRNA data doesn't represent a meaningful subset of anything -- non-tRNA
    plots must only ever be generated for the complete (non-split) variant, regardless of whether
    adata.uns['nontRNA_counts'] happens to be present on this view. `feature_label` names what's
    being skipped for the log message (e.g. 'PCA plots', 'correlation matrices').
    '''
    if not is_full_variant:
        return None, (
            f"Skipping non-tRNA {feature_label}: only ever generated for the complete (non-split) "
            f"variant, since non-tRNA reads can't be meaningfully partitioned by a tRNA length split."
        )
    nontrna_df = adata.uns.get('nontRNA_counts')
    if nontrna_df is None or nontrna_df.empty:
        return None, (
            f"No non-tRNA feature counts found in AnnData object (uns['nontRNA_counts'] missing or "
            f"empty). Skipping non-tRNA {feature_label}. Re-run `trnagraph analyze build` with --gtf "
            f"to enable these {feature_label}."
        )
    return nontrna_df, None


def sort_temp_dir(output_path: str) -> str:
    '''
    Directory to use for a `samtools sort -T` prefix, given the sort's final output file path.
    Defaults to the output file's own directory rather than the system temp dir
    (tempfile.gettempdir(), i.e. /tmp unless $TMPDIR is set) -- tRAX historically hit a server
    /tmp-fills-up failure on large samples and fixed it by writing temp files next to the
    invocation/output instead, and tRNAgraph inherited the same class of risk by relying on the
    system temp dir here.
    '''
    return os.path.dirname(os.path.abspath(output_path))


def variant_dir_names(args, tag: Optional[str] = None) -> Tuple[str, str]:
    '''
    Resolve (results_dir_name, graphs_dir_name) for one `analyze build`/`analyze addsplit`
    variant: variants nest as subfolders of a shared `results`/`graphs` root instead of parallel
    sibling directories (the old `results_u60`/`graphs_u60`). The single source of truth for
    this naming, shared by
    `adataBuild.py` (the AnnData build pipeline) and `toolsSplit.py` (BAM splitting, which
    writes its own per-variant `mapinfo.txt`/`trnamapinfo.txt` into the same directories) --
    keeping both in step matters, since one writes what the other later reads.

    A specific split-variant `tag` (e.g. 'u60'/'o60') always nests as `results/<tag>`/
    `graphs/<tag>`. The default/full variant (`tag=None`) only nests under `results/complete`/
    `graphs/complete` when `args.readlengthsplit` is actually set for this build -- i.e. only
    when there's more than one variant on disk that needs disambiguating. A plain build with no
    split keeps writing straight to `results`/`graphs`, unchanged, so the common (non-split)
    case has no behavior change at all.
    '''
    if tag is not None:
        return os.path.join('results', tag), os.path.join('graphs', tag)
    if getattr(args, 'readlengthsplit', None):
        return os.path.join('results', 'complete'), os.path.join('graphs', 'complete')
    return 'results', 'graphs'


class _TailCaptureHandler(logging.Handler):
    '''
    Captures the last N formatted log messages into a bounded deque, for a live-updating
    scrolling "tail" panel underneath a rich progress bar/spinner.
    '''
    def __init__(self, maxlen: int = 10):
        super().__init__()
        self.lines: deque = deque(maxlen=maxlen)

    def emit(self, record: logging.LogRecord) -> None:
        self.lines.append(self.format(record))


def progress_iterator(
    iterable: Iterable[Any], total: int, desc: str, logger: logging.Logger,
    quiet: bool = False, isatty_fn: Optional[Callable[[], bool]] = None,
) -> Generator[Any, None, None]:
    '''
    Shared progress-reporting helper for long-running per-item loops (trimming, mapping,
    counting, graphing).

    - `quiet=True`: never starts the interactive rich display (it writes straight to the real
      terminal, bypassing the logging system entirely) -- but percentage-milestone logging still
      runs below, since file persistence via the logging system is unconditional. Console
      visibility is controlled separately, by whether cli.py's configure_logging() attached a
      console handler under `--quiet`.
    - A real interactive terminal and not quiet: a live rich spinner (switching to a determinate
      bar once progress starts) plus a scrolling tail of the most recent log messages, driven by
      temporarily swapping out the shared 'trnagraph' logger's console-tagged handler for a local
      _TailCaptureHandler on `logger` itself. Any FileHandler on 'trnagraph' is left untouched, so
      the persisted .log/ file keeps receiving every message regardless of the live display.
    - Otherwise (non-tty, or quiet): percentage-milestone INFO log lines at ~10% increments
      (falling back to every item for totals under 10), in the exact "N/total (P%) complete"
      format that toolsTestSuite.py's own live-box display parses to drive its own progress bar.
    '''
    if isatty_fn is None:
        isatty_fn = sys.stderr.isatty
    is_tty = isatty_fn()

    if is_tty and not quiet:
        tail_handler = _TailCaptureHandler()
        tail_handler.setFormatter(logging.Formatter('%(message)s'))
        logger.addHandler(tail_handler)

        trnagraph_logger = logging.getLogger('trnagraph')
        removed_console_handlers = [h for h in trnagraph_logger.handlers if getattr(h, '_is_console_handler', False)]
        for h in removed_console_handlers:
            trnagraph_logger.removeHandler(h)

        spinner = Spinner('dots', text=desc, style='green')
        progress = Progress(
            TextColumn('[green]{task.description}'),
            BarColumn(complete_style='green', finished_style='bright_green'),
            TaskProgressColumn(),
        )
        task_id = progress.add_task(desc, total=total)
        started_bar = False

        def render():
            body = progress if started_bar else spinner
            tail_text = '\n'.join(tail_handler.lines)
            return Group(body, Text(tail_text)) if tail_text else body

        try:
            # console=Console(stderr=True) is required, not cosmetic: cli.py's handle_output()
            # wraps every command in contextlib.redirect_stdout(tee) so plain print()s reach both
            # the console and the persisted .log/ file. rich.live.Live's internal Console
            # defaults to sys.stdout unless told otherwise -- left implicit, every spinner/
            # progress-bar frame's raw ANSI escape codes get captured into `tee` right along with
            # everything else, flooding the log file. redirect_stdout only ever touches stdout,
            # so pointing Live at stderr instead makes it write straight to the real terminal,
            # bypassing that capture entirely -- matching isatty_fn's own use of stderr above.
            with Live(get_renderable=render, refresh_per_second=4, transient=True, console=Console(stderr=True)):
                for i, item in enumerate(iterable):
                    yield item
                    started_bar = True
                    progress.update(task_id, completed=i + 1)
        finally:
            logger.removeHandler(tail_handler)
            for h in removed_console_handlers:
                trnagraph_logger.addHandler(h)
        return

    milestone_step = max(1, total // 10)
    for i, item in enumerate(iterable):
        yield item
        completed = i + 1
        if completed % milestone_step == 0 or completed == total:
            pct = int(completed / total * 100)
            logger.info(f"{desc}: {completed}/{total} ({pct}%) complete")


class PhaseTracker:
    '''
    Shared outer-progress tracker for a fixed, named sequence of build/graph phases. Each phase
    is announced via one INFO log line the moment it *completes*, weighted by an optional
    per-phase `weight` (default 1, i.e. equal weighting) -- the outer percentage is cumulative
    completed weight over total weight, so callers that want proportional progress (e.g. a future
    graphing command weighting each graph-type phase by its own plot count) get it for free by
    just passing that count as the weight, with no extra bookkeeping.

    Deliberately NOT a rich Live/spinner owner like progress_iterator: phases that have their own
    inner per-item loop (e.g. toolsCountReads.py's per-sample counting, already wired via
    progress_iterator) keep using that exactly as before for fine-grained feedback -- this class
    only reports coarse, phase-level progress on top, so a real terminal running a phase-tracked
    command directly sees both: a plain log line at each phase boundary from this class, plus
    progress_iterator's own rich bar/milestones whenever the current phase happens to wrap one.

    The message format ("<desc> phase N/Total (P%) complete: <label>") is deliberately NOT the
    same shape as progress_iterator's bare "N/Total (P%) complete" -- note the literal word
    "phase" immediately before the fraction. toolsTestSuite.py's _LiveBoxHandler uses that to tell
    the two apart, so a phase-tracked command's box only ever reflects genuine phase-level
    progress and never gets overridden mid-phase by an inner per-item milestone. That was the
    actual bug this was built to fix: Stage 3's per-sample counting milestones reached "10/10
    (100%) complete" almost immediately during a build, pinning the box at 100% for everything
    that ran afterward (DESeq2 fitting, coverage generation, VST, writing the h5ad).

    Some phases wrap their own inner per-item loop whose item count happens to equal that phase's
    own weight (e.g. a graphing command's "coverage" phase, weighted by its expected plot count,
    generating one plot at a time) -- call `advance()` once per completed inner item from within
    that phase's `with phase():` block to tick the outer percentage in lockstep, rather than
    waiting for the phase to complete atomically. `advance()` throttles its own milestone logging
    to ~10%-of-that-phase's-weight steps (mirroring progress_iterator's own convention), so a
    phase with a large weight doesn't flood the log with one line per item.
    '''
    def __init__(
        self, phases: List[str], logger: logging.Logger, desc: str = "Build",
        weights: Optional[List[float]] = None,
    ):
        if weights is None:
            weights = [1] * len(phases)
        if len(weights) != len(phases):
            raise ValueError("weights must be the same length as phases")
        self.phases = list(phases)
        self.weights = list(weights)
        self.total_weight = sum(self.weights) or 1
        self.logger = logger
        self.desc = desc
        self._done_weight: float = 0
        self._index = 0
        self._active: Optional[Dict[str, Any]] = None

    @contextlib.contextmanager
    def phase(self, variant: Optional[str] = None) -> Generator[None, None, None]:
        '''
        Wrap one phase's work. Advances to the next phase in `phases` automatically on each call
        -- pass `variant` (e.g. "Under60") to fold a split-build variant's name into the log line
        without treating it as a separate nesting level (the phase sequence just restarts/repeats
        per variant, labeled, rather than adding a third rendered level).

        `_index` auto-increments assuming phases complete in the SAME order they were declared in
        `phases=[...]` -- a caller that fans work out to a `multiprocessing.Pool` must use ordered
        `imap()`, not `imap_unordered()`, if the pooled work is itself wrapped in phase tracking.
        '''
        if self._index >= len(self.phases):
            raise IndexError(f"{self.desc}: no more phases registered (declared {len(self.phases)})")
        index = self._index
        self._index += 1
        name = self.phases[index]
        weight = self.weights[index]
        full_name = f"[{variant}] {name}" if variant else name
        self._active = {'index': index, 'weight': weight, 'full_name': full_name, 'progress': 0}
        try:
            yield
        finally:
            progress = self._active['progress']
            self._active = None
            remaining = weight - progress
            if remaining:
                self._done_weight += remaining
        pct = int(round(100 * self._done_weight / self.total_weight))
        self.logger.info(f"{self.desc} phase {index + 1}/{len(self.phases)} ({pct}%) complete: {full_name}")

    def advance(self, amount: float = 1) -> None:
        '''
        Partial progress within the CURRENT phase -- must be called from inside a `with
        phase():` block. See the class docstring for the lockstep use case.
        '''
        if self._active is None:
            raise RuntimeError(f"{self.desc}: advance() called outside of an active phase")
        self._active['progress'] += amount
        self._done_weight += amount
        weight = self._active['weight']
        progress = self._active['progress']
        milestone_step = max(1, weight // 10)
        if progress % milestone_step == 0 or progress >= weight:
            pct = int(round(100 * self._done_weight / self.total_weight))
            self.logger.info(f"{self.desc} phase {self._active['index'] + 1}/{len(self.phases)} ({pct}%): {self._active['full_name']}")


from .toolsSchemas import VariantTag

# ---------------------------------------------------------------------------
# Read basis: unique (transcript-specific) vs. all reads
# ---------------------------------------------------------------------------
# "Unique" means TRANSCRIPT-SPECIFIC throughout both tRNAgraph and tRAX: a read whose
# bowtie2 YR tag shows it aligned to exactly one mature-tRNA transcript. It is NOT
# genome-level uniqueness. toolsCountReads.getbamcounts() gates adduniquecount() on
# isuniquetrnamapping() (the YR tag), and toolsGetCoverage fills 'uniquecoverage' from the
# same predicate, so obs['nreads_*_unique_*'] and var 'uniquecoverage' are the same concept.
# The genome MAPQ >= 2 filter is a separate, always-on prefilter sitting beneath every read
# basis and every coverage category -- see tests/unit/test_filtermultimapped_default.py.
#
# Every `trnagraph graph` plot selects its counts through resolve_readtype(), so a single
# --allreads switches the entire command at once. Readtypes therefore arrive WITHOUT a
# '_unique' component: accepting one would let a caller reintroduce, per graph type, exactly
# the silent cross-plot inconsistency --allreads exists to remove.

def load_style_file(path, logger=None):
    '''
    Read and validate a `--style` JSON file, returning a StyleFile.

    Shared by `analyze graph` and `preprocess trim` so both accept the same file format --
    trim only reads the colors block, but there is no reason for it to need a different file.
    A file in the old `--colormap` shape is accepted unchanged; see StyleFile.
    '''
    import json
    from pydantic import ValidationError
    from .toolsSchemas import StyleFile

    if logger is not None:
        logger.info(f'Loading style file: {path}')
    with open(path, 'r') as handle:
        raw = json.load(handle)
    try:
        return StyleFile.model_validate(raw)
    except ValidationError as exc:
        from .toolsSchemas import explain_rejected_keys

        # The pydantic report is kept in full -- it names the exact location, which matters
        # when several keys are wrong at once -- and the guidance is appended to it.
        detail = '\n'.join([f'Invalid style file {path}:', str(exc)]
                           + explain_rejected_keys(exc, 'style'))
        raise InvalidParameterError(detail) from exc


# Built-in presentation defaults, used when a --style file does not set a key. dpi/fonttype
# were previously set as module-level rcParams in each plots*.py file independently.
PLOT_STYLE_DEFAULTS = {
    'format': 'pdf',
    'dpi': 300,
    'alpha': None,          # None = each plot keeps its own tuned opacity
    'rasterize_over': None, # None = never rasterize
}


def resolve_plot_style(style, graph_type, **overrides):
    '''
    Presentation settings for one graph type: built-in defaults, then the style file's
    `defaults`, then that graph type's block, then any CLI flag passed in `overrides`.

    A CLI flag always wins over the file -- passing None for an override means "not
    specified", so a flag left at its default never masks the file.
    '''
    from . import plotsPalette

    resolved = dict(PLOT_STYLE_DEFAULTS)
    if style is not None:
        resolved.update(style.resolve(graph_type))
    resolved.update({k: v for k, v in overrides.items() if v is not None})
    # The palette rides in `settings` rather than being threaded as a new parameter: every
    # plot module already receives this dict, and tests/unit/test_plot_settings_scope.py
    # already statically guards the failure mode a NEW threaded parameter reintroduces (a
    # helper reading it without being given it, invisible until it runs inside a pool).
    # Every graph type gets every gradient role -- see plotsPalette.resolve_gradients.
    resolved['gradients'] = plotsPalette.resolve_gradients(
        getattr(style, 'gradients', None) if style is not None else None)
    resolved['categorical'] = getattr(style, 'categorical', None) if style is not None else None
    return resolved


def figsize_for(settings, default):
    '''
    The --style figure size for an INDIVIDUAL plot, or the module's own tuned default.

    style_context() puts `figsize` into rcParams, but every plot module passes an explicit
    figsize= to plt.figure/plt.subplots, and an explicit argument beats the rcParam -- so the
    setting silently did nothing anywhere. Individual plots go through here instead.

    Combined and multi-page pages deliberately do NOT: they compute their geometry from how
    many panels they are laying out, so a fixed size would clip panels or leave dead space.
    '''
    configured = (settings or {}).get('figsize')
    return tuple(configured) if configured else default


def linewidth_for(settings, default):
    '''
    The --style line width for an INDIVIDUAL plot, or the module's own tuned default.

    Shrinking `figsize` scales a figure's geometry but not its stroke widths, so a small
    coverage plot renders with traces disproportionately heavy for the panel they sit in.
    This is an ABSOLUTE point value replacing the module default, not a multiplier scaling
    it: a user asking for hairlines wants hairlines, not "whatever this module chose, times
    0.4".

    Data traces and bar edges only. The dashed modification guides, grid lines and arm-band
    shading keep their own widths -- that furniture is structural rather than a presentation
    knob, the same reasoning that withholds `alpha` from coverage.

    Individual plots only, exactly as with figsize_for(): combined and multi-page pages
    compute their own geometry, are already unaffected by `figsize`, and stay static.
    '''
    configured = (settings or {}).get('line_width')
    return configured if configured else default


@contextlib.contextmanager
def style_context(settings):
    '''
    Apply the rcParam-expressible style settings for the duration of one graph type's run.

    font_size and figsize have to be in force while a figure is *created*, not when it is
    saved, so they go through rcParams here rather than through save_current(). figsize is
    applied as the default figure size, which individual plots pick up; combined/multi-page
    pages call plt.subplots/figure with an explicit size computed from their panel count and
    so are unaffected, which is the intended split.
    '''
    import matplotlib.pyplot as plt

    settings = settings or {}
    overrides = {}
    if settings.get('font_size'):
        overrides['font.size'] = settings['font_size']
    if settings.get('figsize'):
        overrides['figure.figsize'] = list(settings['figsize'])
    if settings.get('dpi'):
        overrides['savefig.dpi'] = settings['dpi']
    if not overrides:
        yield
        return
    with plt.rc_context(overrides):
        yield


def save_current(path, settings=None, **savefig_kwargs):
    '''
    Save the current matplotlib figure to `path`, swapping its extension for the configured
    output format and applying the configured dpi. Returns the path actually written.

    Individual plots go through here so `--format svg` reaches all of them. Multi-page
    combined outputs deliberately do NOT: they are built with PdfPages, and neither SVG nor
    PNG has a multi-page concept, so those stay PDF whatever the format setting says.
    '''
    import matplotlib.pyplot as plt

    settings = settings or PLOT_STYLE_DEFAULTS
    fmt = settings.get('format') or 'pdf'
    root, _ = os.path.splitext(path)
    out = f'{root}.{fmt}'
    kwargs = {'bbox_inches': 'tight'}
    kwargs.update(savefig_kwargs)
    if settings.get('dpi'):
        kwargs['dpi'] = settings['dpi']
    plt.savefig(out, **kwargs)
    return out


def save_figure(fig, path_without_extension, settings=None, **savefig_kwargs):
    '''
    Write `fig` in the configured format, returning the path written.

    Centralises what 35 call sites used to hardcode as '.pdf'. Format and dpi come from the
    resolved style settings, so `--format svg` reaches every plot without each module
    knowing about the flag.
    '''
    settings = settings or PLOT_STYLE_DEFAULTS
    fmt = settings.get('format') or 'pdf'
    path = f'{path_without_extension}.{fmt}'
    kwargs = {'bbox_inches': 'tight'}
    kwargs.update(savefig_kwargs)
    if settings.get('dpi'):
        kwargs['dpi'] = settings['dpi']
    fig.savefig(path, **kwargs)
    return path


READ_BASIS_UNIQUE = 'unique'
READ_BASIS_ALL = 'allreads'

#: Readtypes carrying BOTH bases in obs (nreads_<rt>_unique_norm and nreads_<rt>_norm).
#: Verified against a built object; these five are the ones adataBuild populates twice.
DUAL_BASIS_READTYPES = ('wholecounts', 'fiveprime', 'threeprime', 'other', 'total')

#: Readtypes adataBuild only ever populates in the all-reads basis -- pre-tRNA and antisense
#: categories, which tRAX never counted transcript-specifically. Requesting one under the
#: unique basis falls back to its all-reads column rather than failing, mirroring how the
#: non-tRNA panels stay all-reads by structural necessity.
ALL_READS_ONLY_READTYPES = ('wholeprecounts', 'partialprecounts', 'trailercounts', 'antisense')

#: --covtype default per basis. tRAX's four coverage categories partition 'coverage', so the
#: all-reads basis is that total and the unique basis is its transcript-specific part.
COVTYPE_DEFAULTS = {READ_BASIS_UNIQUE: 'uniquecoverage', READ_BASIS_ALL: 'coverage'}

# tRAX bins each read into exactly one of four categories by how specifically it could be
# assigned, using bowtie2's YM/YA/YR tags (counts of aminos/anticodons/transcripts the read
# hit). getcoverage.py's getsamplecoverage() does this as an if/elif chain, so the four are
# mutually exclusive and sum to 'coverage'. They are a PARTITION, not a filter -- the
# genome MAPQ >= 2 prefilter sits beneath all four equally.
#
# The display labels are tRAX's own, from newcoverageplots.R's column renaming, and read as
# "the finest level this read resolves to". They are kept verbatim so a tRAX user recognises
# them; the CLI aliases are the short forms of the same idea.

#: --covtype alias -> adata.var['coverage'] value. Raw var values stay accepted so nothing
#: that already scripted --covtype uniquecoverage breaks.
COVERAGE_CATEGORY_ALIASES = {
    'unique': 'uniquecoverage',
    'transcript': 'uniquecoverage',
    'isodecoder': 'multitrnacoverage',
    'isotype': 'multianticodoncoverage',
    'notamino': 'multiaminocoverage',
    'total': 'coverage',
}

#: adata.var['coverage'] value -> the alias used for its output directory.
COVERAGE_CATEGORY_DIRS = {
    'uniquecoverage': 'unique',
    'multitrnacoverage': 'isodecoder',
    'multianticodoncoverage': 'isotype',
    'multiaminocoverage': 'notamino',
    'coverage': 'total',
}

#: Display labels, verbatim from newcoverageplots.R.
COVERAGE_CATEGORY_LABELS = {
    'uniquecoverage': 'Transcript Specific',
    'multitrnacoverage': 'Isodecoder Specific',
    'multianticodoncoverage': 'Isotype Specific',
    'multiaminocoverage': 'Not Amino Specific',
}

#: Stacking order, least specific first, matching the factor levels newcoverageplots.R sets
#: so the stacked overview reads the same way round as tRAX's.
COVERAGE_PARTITION = (
    'multiaminocoverage', 'multianticodoncoverage', 'multitrnacoverage', 'uniquecoverage',
)


def coverage_category_dir(covtype: str) -> str:
    '''
    Output subdirectory for a resolved --covtype. Partition members get their short alias;
    anything else (readstarts, mismatchedbases, ...) uses its own name, so every coverage
    run lands in exactly one directory named for what it plotted.
    '''
    return COVERAGE_CATEGORY_DIRS.get(covtype, covtype)


def read_basis(allreads: bool) -> str:
    '''Map the --allreads flag onto a basis token. The default (flag absent) is unique.'''
    return READ_BASIS_ALL if allreads else READ_BASIS_UNIQUE


def resolve_readtype(readtype: str, basis: str, adata: Optional[ad.AnnData] = None) -> str:
    '''
    Turn a bare readtype (e.g. 'total') into the obs column for `basis`, e.g.
    'nreads_total_unique_norm' or 'nreads_total_norm'.

    Rejects a readtype that already names a basis ('total_unique'), because --diffrts and
    --pcareadtypes deliberately no longer carry that dimension -- honouring it would let one
    graph type sit on a different denominator than the rest without saying so.

    Falls back to the all-reads column, with a warning, for a readtype that has no unique
    counterpart at all (see ALL_READS_ONLY_READTYPES). When `adata` is supplied the fallback
    is decided by what obs actually contains rather than by the static list, so an object
    built by a future adataBuild that does populate them needs no change here.
    '''
    logger = logging.getLogger(__name__)
    if readtype.startswith('nreads_'):
        raise ValueError(
            f"Readtype '{readtype}' is already a full obs column name. Pass a bare readtype "
            f"(one of {list(DUAL_BASIS_READTYPES)}) and let --allreads choose the basis."
        )
    if '_unique' in readtype or readtype.endswith('unique'):
        raise ValueError(
            f"Readtype '{readtype}' names a read basis. Basis is set once for the whole "
            f"`graph` command by --allreads, so pass '{readtype.replace('_unique', '')}' "
            f"instead of '{readtype}'."
        )
    if basis not in (READ_BASIS_UNIQUE, READ_BASIS_ALL):
        raise ValueError(f"Unknown read basis '{basis}'. Expected one of "
                         f"{[READ_BASIS_UNIQUE, READ_BASIS_ALL]}.")

    all_reads_column = f'nreads_{readtype}_norm'
    if basis == READ_BASIS_ALL:
        return all_reads_column

    unique_column = f'nreads_{readtype}_unique_norm'
    if adata is not None:
        available = unique_column in adata.obs.columns
    else:
        available = readtype not in ALL_READS_ONLY_READTYPES
    if available:
        return unique_column
    logger.warning(
        f"Readtype '{readtype}' has no transcript-specific (unique) counts; plotting it from "
        f"all reads ({all_reads_column}) instead. Pre-tRNA and antisense categories are only "
        f"ever counted across all reads."
    )
    return all_reads_column


def resolve_covtype(covtype: Optional[str], basis: str) -> str:
    '''
    Resolve --covtype. An explicit value is always honoured -- the four tRAX coverage
    categories are a partition a user may legitimately want to inspect in either mode -- so
    the basis only supplies the default when none was given.
    '''
    if covtype:
        return COVERAGE_CATEGORY_ALIASES.get(covtype, covtype)
    return COVTYPE_DEFAULTS[basis]


_VARIANT_LAYER_MAP = {'norm': 'norm', 'raw': 'raw', 'allfeatures': 'norm_allfeatures', 'vst': 'vst'}

# Maps a `uns['size_splits'][tag]` key to the default/full uns key it stands in for
# once overlaid onto a resolved variant view (see build_variant_view()). 'sample_cluster_umap'
# is NOT here -- it's obs-aligned data stored in obsm (as f'sample_cluster_umap_{tag}' for a
# split variant), overlaid separately below, not through this uns-only map.
_VARIANT_UNS_KEY_MAP = {
    'sizefactors_trna': 'deseq2_sizefactors_trna',
    'sizefactors_allfeatures': 'deseq2_sizefactors_allfeatures',
    'type_counts': 'type_counts',
    'type_real_counts': 'type_real_counts',
    'amino_counts': 'amino_counts',
    'anticodon_counts': 'anticodon_counts',
    'nontRNA_counts': 'nontRNA_counts',
    'mismatch_counts': 'mismatch_counts',
    'log2FC': 'log2FC',
    'cluster_runinfo': 'cluster_runinfo',
    'group_cluster_umap': 'group_cluster_umap',
}


def parse_variant(adata: ad.AnnData, variant: str = 'norm:full') -> 'VariantTag':
    '''
    Parse a `--variant` string of the form "<norm>:<tag>" (tag defaults to 'full' if omitted).
    <norm> must be one of 'norm', 'raw', 'allfeatures', 'vst'. If <tag> is not 'full', it must
    already exist in adata.uns['size_splits'] (i.e. a split variant added via `analyze build
    --readlengthsplit`/`analyze addsplit`).
    '''
    norm, sep, tag = variant.partition(':')
    tag = tag if sep else 'full'
    if norm not in _VARIANT_LAYER_MAP:
        raise ValueError(f"Unknown normalization '{norm}' in --variant '{variant}'. Expected one of: {list(_VARIANT_LAYER_MAP)}.")
    if tag != 'full':
        available = list(adata.uns.get('size_splits', {}).keys())
        if tag not in available:
            raise ValueError(f"Split tag '{tag}' not found in this AnnData object (in --variant '{variant}'). Available split tags: {available or 'none (no splits added yet)'}.")
    return VariantTag(raw=variant, norm=norm, tag=tag)


def build_variant_view(adata: ad.AnnData, spec: 'VariantTag') -> ad.AnnData:
    '''
    Return a working COPY of `adata` with .X swapped to the layer `spec` selects, and, for a
    split variant (spec.tag != 'full'), the split's adata.obsm[f'size_split_{tag}'] columns
    overlaid onto the copy's .obs, and applicable adata.uns['size_splits'][tag] entries
    overlaid onto the copy's .uns -- all under the same (unsuffixed) names used by the
    default/full variant, so downstream plotting/clustering/DE code needs no changes at all.

    IMPORTANT: this returns a disposable working copy. Never write it back to `adata`'s
    original h5ad path -- doing so would overwrite the real full/default variant's data with
    the split variant's overlaid values. Callers that compute new results against this view
    (e.g. a fresh log2FC cache entry, new cluster output) must write those results back onto
    the ORIGINAL (unresolved) adata, into its namespaced location (uns['size_splits'][tag],
    obsm[f'size_split_{tag}']), not onto this view.
    '''
    view = adata.copy()
    if spec.tag == 'full':
        if spec.norm != 'norm':
            view.X = view.layers[_VARIANT_LAYER_MAP[spec.norm]]
        return view

    layer_name = f'{_VARIANT_LAYER_MAP[spec.norm]}_{spec.tag}' if spec.norm != 'norm' else f'norm_{spec.tag}'
    if layer_name not in view.layers:
        # All-feature normalization is deliberately complete-variant only: a split variant has
        # its non-tRNA features excluded, so there is no all-feature size-factor set to
        # normalize against. Say so, rather than sending the reader looking for a build flag
        # that would produce the layer -- there isn't one.
        if spec.norm == 'allfeatures' and spec.tag != 'full':
            raise ValueError(
                f"--variant '{spec.raw}' is not available: all-feature normalization is only "
                f"computed for the complete (non-split) variant. A read-length split excludes "
                f"non-tRNA features entirely, so a split has no all-feature size factors to "
                f"normalize against. Use 'norm:{spec.tag}', 'raw:{spec.tag}' or 'vst:{spec.tag}' "
                f"for this split, or 'allfeatures:full' for the complete variant."
            )
        raise ValueError(f"--variant '{spec.raw}' resolves to layer '{layer_name}', which is not present in this AnnData object (was this normalization computed for this split?).")
    view.X = view.layers[layer_name]

    split_obs = adata.obsm.get(f'size_split_{spec.tag}')
    if split_obs is not None:
        for col in split_obs.columns:
            view.obs[col] = split_obs[col].reindex(view.obs.index).values

    split_sample_umap = adata.obsm.get(f'sample_cluster_umap_{spec.tag}')
    if split_sample_umap is not None:
        view.obsm['sample_cluster_umap'] = split_sample_umap.reindex(view.obs.index)

    split_uns = adata.uns.get('size_splits', {}).get(spec.tag, {})
    for split_key, default_key in _VARIANT_UNS_KEY_MAP.items():
        if split_key in split_uns:
            view.uns[default_key] = split_uns[split_key]

    return view


def scatter_subset_graph_to_full(subset_graph, subset_index, full_index):
    '''
    Embed a sparse pairwise graph computed over a SUBSET of observations (e.g. UMAP's
    fuzzy-simplicial-set `.graph_`, sized [n_subset, n_subset]) into a full [n_obs, n_obs]
    sparse matrix aligned to `full_index`, so it can be stored in adata.obsp. Rows/columns
    for observations not present in `subset_index` (e.g. excluded by --readcutoff or the
    'Und' amino-acid filter) are left as all-zero connectivity.
    '''
    import scipy.sparse as sp
    subset_graph = sp.csr_matrix(subset_graph)
    n_full = len(full_index)
    full_pos = {name: i for i, name in enumerate(full_index)}
    rows = [full_pos[name] for name in subset_index if name in full_pos]
    if len(rows) != len(subset_index):
        missing = [name for name in subset_index if name not in full_pos]
        raise ValueError(f"scatter_subset_graph_to_full: {len(missing)} subset observations not found in full_index (e.g. {missing[:3]}).")
    row_map = np.array(rows)
    # Embed subset_graph[i, j] at full[row_map[i], row_map[j]] via a sparse selector matrix.
    selector = sp.csr_matrix((np.ones(len(row_map)), (np.arange(len(row_map)), row_map)), shape=(len(row_map), n_full))
    full_graph = selector.T @ subset_graph @ selector
    return full_graph.tocsr()


#: Shrinkage estimators for the log2 fold changes. PyDESeq2 implements apeGLM only;
#: tRAX's own DESeq2 betaPrior has no equivalent here, so it is not offered.
SHRINK_METHODS = ('apeGLM', 'none')


def _as_pair_list(pairs) -> List[Tuple[str, str]]:
    '''
    Coerce a `pairs` value to a list of plain string tuples, whichever path produced it.

    A freshly computed value is a list of tuples straight out of itertools.combinations, but the
    same value read back from a written .h5ad is a numpy array of shape (n, 2): anndata stores a
    list of same-length sequences as an array and has no way to restore the original container.
    Consumers therefore saw a different type on a warm object than on a cold one, and one of them
    (plotsVolcano's `if pairs:`) raised on the array form -- so a combined volcano overview was
    only ever written on the first run for a given config_name/cutoff.

    Normalizing here, at the single boundary both paths pass through, is what keeps that
    asymmetry from reaching any consumer at all.
    '''
    return [tuple(str(level) for level in pair) for pair in pairs]


class adataLog2FC:
    def __init__(self, adata: ad.AnnData, compare: str, readtype: str, readcount_cutoff: int = 80, config_name: str = 'default', overwrite: bool = False, n_cpus: Optional[int] = 1, shrink: str = 'apeGLM'):
        self.adata = adata
        self.compare = compare
        self.readtype = readtype
        self.readcount_cutoff = str(readcount_cutoff)
        self.config_name = config_name
        self.overwrite = overwrite
        # Default of 1 is a safety requirement, not a performance tweak: this class is called
        # from inside adataGraph.py's own multiprocessing.Pool worker processes (plotsVolcano.py,
        # plotsHeatmap.py), and PyDESeq2 defaults to using ALL available CPUs via joblib's loky
        # backend -- spawning a second real process pool from inside an already-forked worker,
        # which deadlocks. joblib's n_jobs=1 is the one setting that skips creating a nested pool
        # entirely (it runs sequentially in-process) rather than just shrinking it, so only a
        # caller that KNOWS it isn't running inside a pool (e.g. adataBuild.py's or
        # adataGraph.py's own pre-pool precompute step) should ever pass something else.
        self.n_cpus = n_cpus
        # How the log2 fold changes are shrunk: 'apeGLM' (Zhu, Ibrahim & Love 2019,
        # doi:10.1093/bioinformatics/bty895) or 'none'. A named method rather than a boolean
        # because the choice of estimator is a real one -- tRAX shrinks via DESeq2's
        # betaPrior=TRUE (analyzecounts.R:104), which PyDESeq2 does not implement, so a third
        # value can be added here without changing the flag's shape. Shrinking costs one
        # DESeq2 fit per distinct baseline level instead of one overall -- see log2fc_df().
        if shrink not in SHRINK_METHODS:
            raise InvalidParameterError(
                f"Unknown --shrink method {shrink!r}. Expected one of: "
                + ', '.join(SHRINK_METHODS))
        self.shrink = shrink
        self.log2fc_dict: Dict[str, Any] = {}
        # Set by main(): whether THIS call actually computed a fresh fit (cache miss) rather than
        # reusing a cached df (cache hit). Callers that precompute several combos and want to
        # know whether anything new was added (e.g. to decide whether to persist to disk) should
        # check this instead of diffing the uns['log2FC'] dict before/after -- a naive dict
        # comparison silently never detects new entries under an already-existing config_name/
        # compare key when using a shallow copy (shared nested dict references short-circuit
        # equality via identity), and raises ValueError ("truth value of a DataFrame is
        # ambiguous") when using a deep copy instead, since dict equality falls through to
        # comparing the nested DataFrame values directly.
        self.computed_fresh = False

    def main(self) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        # Search nested adata.uns for log2FC for compare, readtype, the readcutoff
        self.log2fc_dict = self.adata.uns.get('log2FC', {})

        # Ensure nested structure exists
        self.log2fc_dict.setdefault(self.config_name, {}) \
                        .setdefault(self.compare, {}) \
                        .setdefault(self.readtype, {}) \
                        .setdefault(self.readcount_cutoff, {})

        entry = self.log2fc_dict[self.config_name][self.compare][self.readtype][self.readcount_cutoff]
        df = entry.get('df', pd.DataFrame())
        # Which shrinkage produced a cached fold change is part of what that value IS, so it
        # belongs in the cache key. An entry written before this key existed reads as None,
        # which correctly misses against any named method and recomputes -- otherwise an object
        # built under the old boolean key would keep serving values from the wrong estimator.
        if entry.get('shrink') != self.shrink:
            df = pd.DataFrame()

        # `pairs` is read by consumers out of the RETURNED dict, not from a local -- so both
        # branches normalize the stored value rather than a variable of their own.
        compare_entry = self.log2fc_dict[self.config_name][self.compare]
        if df.empty or self.overwrite:
            df, pairs = self.log2fc_df()
            entry['df'] = df
            entry['shrink'] = self.shrink
            compare_entry['pairs'] = _as_pair_list(pairs)
            self.adata.uns['log2FC'] = self.log2fc_dict
            self.computed_fresh = True
        else:
            # Normalized in place: a cached entry that came back from a written .h5ad holds an
            # ndarray, and leaving it there would hand the next reader the other type again.
            compare_entry['pairs'] = _as_pair_list(compare_entry['pairs'])

        return df, self.log2fc_dict

    def log2fc_df(self, index_col: str = 'trna') -> Tuple[pd.DataFrame, List[Tuple[Any, Any]]]:
        '''
        Compute pairwise log2FC/significance between every pair of self.compare's levels,
        using PyDESeq2's own dispersion/GLM model rather than a manual two-sample t-test on
        precomputed per-group mean/std -- this gets a real negative-binomial fit (appropriate
        for count data, unlike a t-test's normality assumption) and BH-adjusted p-values
        (multiple-testing correction the old t-test never applied) essentially for free, since
        adataBuild.py already runs the equivalent PyDESeq2 flow for the build-time --pairs
        comparisons. Every caller of adataLog2FC passes a NORMALIZED ('..._norm') readtype
        (used for the readcount_cutoff filter, matching prior behavior exactly); this method
        derives the matching RAW ('..._raw') column to actually feed PyDESeq2, since it does
        its own internal normalization and dispersion modeling -- feeding it already-normalized
        values would double-normalize and bias the fit.
        '''
        if not self.readtype.endswith('_norm') or self.readtype.replace('_norm', '_raw') not in self.adata.obs.columns:
            raise ValueError(
                f"log2fc_df expects a '..._norm' readtype with a matching '..._raw' obs column "
                f"(got '{self.readtype}')."
            )
        raw_readtype = self.readtype.replace('_norm', '_raw')
        obs = self.adata.obs

        # Same readcount-cutoff filter as before: mean of the per-compare-group NORMALIZED
        # readtype average must clear the cutoff. `index_col` is the FEATURE AXIS the fit is
        # taken over -- 'trna' for the heatmap and volcano, 'amino'/'iso' for `-g compare`,
        # which plots per amino acid and per isoacceptor and so could not use this method at
        # all while both pivots below were hardwired to 'trna'.
        mdf = obs.pivot_table(index=index_col, columns=self.compare, values=self.readtype, aggfunc='mean', observed=True)
        keep_features = mdf.index[mdf.mean(axis=1) >= int(self.readcount_cutoff)]

        # Create permutations of pairings of groups for heatmap
        pairs = list(itertools.combinations(mdf.columns, 2))
        # Always create the log2_<pair>/pval_<pair> columns (pairs come from self.compare's
        # column labels, independent of which/how many trna rows pass the cutoff below) even
        # if no rows end up surviving -- callers (e.g. plotsVolcano.py) index these columns
        # unconditionally, so an empty-but-column-less df breaks them with a KeyError.
        pair_columns = sorted(c for p in pairs for c in (f'log2_{p[0]}-{p[1]}', f'pval_{p[0]}-{p[1]}'))
        df_pairs = pd.DataFrame(index=keep_features, columns=pair_columns, dtype=float)
        if not pairs or len(keep_features) == 0:
            return df_pairs, pairs

        # Build the RAW per-(trna, sample) counts matrix PyDESeq2 needs, and each sample's
        # condition (self.compare is a per-sample covariate: samples sharing a value are
        # replicates of that group).
        # Counts aggregate additively, so a coarser axis SUMS its tRNAs into one feature -- an
        # amino acid's count in a sample is the total over its isodecoders. The 'trna' axis
        # keeps 'first' rather than being folded into the same branch: one row per
        # (trna, sample) makes the two identical, and leaving that path untouched means adding
        # this parameter cannot move a single existing number.
        wide_raw = obs.pivot_table(index=index_col, columns='sample', values=raw_readtype,
                                   aggfunc='first' if index_col == 'trna' else 'sum').loc[keep_features]
        sample_df = obs.drop_duplicates('sample').set_index('sample')
        # self.compare == 'sample' (the Parameter Fallback default) is a degenerate but valid
        # case: 'sample' becomes the index via set_index above, so it's no longer a column to
        # select -- each sample is trivially its own condition, i.e. its own index value.
        if self.compare == 'sample':
            sample_condition = pd.Series(sample_df.index, index=sample_df.index)
        else:
            sample_condition = sample_df[self.compare]

        counts_df = wide_raw.T.fillna(0).clip(lower=0).round().astype(int)  # samples x feature, as PyDESeq2 expects
        meta_df = pd.DataFrame({'condition': sample_condition.reindex(counts_df.index)}).dropna()
        counts_df = counts_df.loc[meta_df.index]

        # size_factors_fit_type='poscounts' avoids PyDESeq2's "iterative" size-factor fallback
        # (a scipy.optimize Powell search over one parameter per sample) that the default
        # 'ratio' method silently falls into on zero-heavy count data like this -- see the VST
        # hang fix in adataBuild.py._compute_vst_ for the full story on why that path is
        # pathologically slow at anything beyond ~50-100 samples.
        # n_cpus=1 is required, not just a performance tweak: this runs inside adataGraph.py's own
        # multiprocessing.Pool worker processes, and PyDESeq2 defaults to using ALL available CPUs
        # via joblib's loky backend -- spawning a second real process pool from inside an already-
        # forked worker, which deadlocks (confirmed live: orphaned trnagraph graph processes found
        # hung 24+ hours later with orphaned loky resource-tracker children attached). joblib's
        # n_jobs=1 is the one setting that skips creating a nested pool entirely rather than just
        # shrinking it, so this must stay pinned at 1 even though DESeq2 could otherwise parallelize.
        def fit(reference=None):
            meta = meta_df
            if reference is not None:
                # apeGLM shrinks a design-matrix COEFFICIENT, not an arbitrary contrast, and
                # the coefficients a fit exposes are all relative to its reference level. So a
                # pair can only be shrunk from a fit whose reference is that pair's own
                # baseline. Ordering the condition categories puts the wanted level first;
                # PyDESeq2's ref_level argument is deprecated.
                levels = [reference] + [lv for lv in sorted(set(meta_df['condition'])) if lv != reference]
                meta = meta_df.assign(condition=pd.Categorical(meta_df['condition'], categories=levels))
            model = DeseqDataSet(counts=counts_df, metadata=meta, design_factors='condition',
                                 size_factors_fit_type='poscounts', quiet=True, n_cpus=self.n_cpus)
            model.deseq2()
            return model

        # Group pairs by the level acting as their baseline, so each distinct reference costs
        # one fit rather than each pair costing one. Without shrinkage a single fit serves
        # every contrast, which is why this is opt-outable.
        if self.shrink != 'none':
            by_reference = defaultdict(list)
            for pair in pairs:
                by_reference[pair[0]].append(pair)
            fits = {reference: fit(reference) for reference in by_reference}
        else:
            by_reference = {None: list(pairs)}
            fits = {None: fit()}

        for reference, reference_pairs in by_reference.items():
            dds = fits[reference]
            for pair in reference_pairs:
                col_name = f'{pair[0]}-{pair[1]}'
                # contrast=[factor, test_level, ref_level] -> log2FoldChange = log2(test/ref);
                # test=pair[1], ref=pair[0] matches the previous log2(mean(pair[1])/mean(pair[0])) convention.
                stat_res = DeseqStats(dds, contrast=['condition', pair[1], pair[0]], quiet=True)
                stat_res.summary()
                if self.shrink != 'none':
                    # Leaves p-values untouched by construction -- shrinkage moves the effect
                    # size, not the significance call.
                    stat_res.lfc_shrink(coeff=f'condition[T.{pair[1]}]')
                res = stat_res.results_df
                df_pairs[f'log2_{col_name}'] = res['log2FoldChange'].reindex(keep_features)
                # BH-adjusted p-value (padj), not the raw per-comparison p-value -- the old
                # t-test never applied any multiple-testing correction.
                df_pairs[f'pval_{col_name}'] = res['padj'].reindex(keep_features)

        # sort the columns alphabetically so log2FC are followed by pvals
        df_pairs = df_pairs.reindex(sorted(df_pairs.columns), axis=1)
        return df_pairs, pairs



def require_options(args, command, required):
    """
    Fail with one sentence when an option a run cannot proceed without is missing.

    These options used to be typer-required, which made them impossible to supply from a
    --config file: click rejected the invocation before the file was ever read. They are now
    optional at the parser and checked here instead, AFTER the config has been merged, so
    either source satisfies them and the error names both.
    """
    missing = [name for name in required if getattr(args, name, None) is None]
    if missing:
        flags = ', '.join(f"--{name.replace('_', '-')}" for name in missing)
        keys = ', '.join(f"'{name}'" for name in missing)
        raise InvalidParameterError(
            f"{command} needs {flags}. Give {'them' if len(missing) > 1 else 'it'} on the "
            f"command line, or set {keys} in your --config file's `flags.{command}` block."
        )


def read_run_config(path):
    """
    Parse and validate a --config file, or return None when none was given.

    Deliberately silent: some commands need a value out of the config (a manifest path, say)
    before their log file exists, so parsing is separated from the logging that
    `apply_config_flags` does once a logger is available.
    """
    import json

    from pydantic import ValidationError

    from .toolsSchemas import RunConfig, explain_rejected_keys

    if not path:
        return None
    try:
        with open(path, 'r') as handle:
            return RunConfig.model_validate(json.load(handle))
    except FileNotFoundError:
        raise InvalidParameterError(f'Config file not found: {path}')
    except json.JSONDecodeError as e:
        raise InvalidParameterError(f'Config file {path} is not valid JSON: {e}')
    except ValidationError as e:
        detail = '\n'.join([f'Invalid config file {path}:', str(e)]
                           + explain_rejected_keys(e, 'config'))
        raise InvalidParameterError(detail)


def load_run_config(path, command, args, logger):
    """
    Read a --config file and apply that command's `flags` block to `args`.

    Shared by every command that takes --config so all eight apply a block identically:
    a value typed on the command line always beats the file, detected through
    `args.cli_specified` (cli.py builds it from click's ParameterSource -- comparing against
    the default cannot work, since nearly every option has a real default rather than a None
    sentinel, so "set to 80" and "left alone" look identical).

    `cli_specified` is absent when a namespace is built directly (the Python API, and tests),
    in which case nothing counts as typed and the file applies in full.
    """
    config = read_run_config(path)
    if config is None:
        return None
    logger.info(f'Loading config file: {path}')
    apply_config_flags(config, command, args, logger)
    return config


def apply_config_flags(config, command, args, logger):
    '''Lay one command's `flags` block over `args`, skipping anything typed on the CLI.'''
    if config is None:
        return
    # Named on every command that reads the file, so a run is traceable to what drove it.
    logger.info(f'Using config "{config.name}"')
    block = getattr(config.flags, command, None) if config.flags else None
    if block is None:
        return
    typed = getattr(args, 'cli_specified', None) or frozenset()
    for key, value in block.model_dump(exclude_none=True).items():
        if key in typed:
            logger.info(f'Config sets {key}, but --{key} was given on the command line; '
                        f'keeping the command-line value.')
            continue
        setattr(args, key, value)
        logger.info(f'Config sets {command}.{key} = {value!r}')


def log2fc_from_wide_df(df: pd.DataFrame, sample_group_map: Dict[str, str], readcount_cutoff: int = 80) -> Tuple[pd.DataFrame, List[Tuple[Any, Any]]]:
    '''
    Compute pairwise log2FC/p-value statistics from a wide (feature x sample) dataframe using an
    explicit sample->group mapping. Mirrors adataLog2FC.log2fc_df()'s statistics (mean/std/count
    per group, ttest_ind_from_stats) but for data that isn't stored per-tRNA in adata.obs, e.g.
    adata.uns['nontRNA_counts'] or a combined tRNA+non-tRNA dataframe, as used by the non-tRNA and
    combined volcano plots in plotsVolcano.py.
    '''
    tdf = df.T
    tdf.index = tdf.index.astype(str)
    groups = tdf.index.map(sample_group_map)

    mdf = tdf.groupby(groups).mean().T
    sdf = tdf.groupby(groups).std().T
    cdf = tdf.groupby(groups).count().T

    mean_drop_list = mdf.mean(axis=1) >= int(readcount_cutoff)
    sdf = sdf[mean_drop_list].dropna()
    mdf = mdf[mean_drop_list].dropna()
    cdf = cdf[mean_drop_list].dropna()

    mdf = mdf.replace(0, 1e-20)

    pairs = list(itertools.combinations(mdf.columns, 2))
    df_pairs = pd.DataFrame()
    for pair in pairs:
        col_name = f'{pair[0]}-{pair[1]}'
        df_pairs[f'log2_{col_name}'] = np.log2(mdf[pair[1]]) - np.log2(mdf[pair[0]])
        _, pval = stats.ttest_ind_from_stats(
            mdf[pair[0]].values, sdf[pair[0]].values, cdf[pair[0]].values,
            mdf[pair[1]].values, sdf[pair[1]].values, cdf[pair[1]].values
        )
        df_pairs[f'pval_{col_name}'] = pval

    df_pairs = df_pairs.reindex(sorted(df_pairs.columns), axis=1)
    return df_pairs, pairs

def read_multi_fasta(fa_file: str) -> Generator[Tuple[str, str], None, None]:
    """
    Iterates over a FASTA file and yields (header, sequence) tuples.
    
    Handles:
    - Standard text files
    - Gzipped files (.gz)
    - Standard Input ('stdin')
    
    Args:
        fa_file (str): Path to the fasta file or 'stdin'.

    Yields:
        Tuple[str, str]: A tuple containing (sequence_id, sequence_string).
    """
    
    # Internal helper to handle file opening context
    def get_file_handle(filename: str) -> TextIO:
        if filename == "stdin":
            return sys.stdin
        elif filename.endswith(".gz"):
            return gzip.open(filename, "rt", encoding='utf-8')
        else:
            return open(filename, "r", encoding='utf-8')

    # Main parsing logic
    # Using a context manager ensures file handles (except stdin) are closed properly
    file_handle = get_file_handle(fa_file)
    
    try:
        current_header = None
        seq_buffer = []

        for line in file_handle:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # If we have a previous record, yield it
                if current_header is not None:
                    # Join buffer list for performance (faster than string += string)
                    yield current_header, "".join(seq_buffer)
                
                # Parse new header
                # Logic matches original regex: r"\>([^\s\,]+)"
                # Splits at first space, then splits at first comma, removing '>'
                current_header = line[1:].split()[0].split(',')[0]
                seq_buffer = []
            else:
                seq_buffer.append(line)

        # Yield the final sequence after loop finishes
        if current_header is not None and seq_buffer:
            yield current_header, "".join(seq_buffer)

    finally:
        # Close file unless it is stdin
        if fa_file != "stdin" and file_handle is not sys.stdin:
            file_handle.close()

def read_rna_stk(stk_file: Iterable[str]) -> Generator['RnaAlignment', None, None]:
    """
    Parses a Stockholm formatted RNA alignment file.

    Args:
        stk_file (Iterable[str]): An iterable yielding lines (e.g., open file handle).

    Yields:
        RnaAlignment: An object containing parsed sequences, structure, and consensus.
    """
    
    sequences: Dict[str, List[str]] = defaultdict(list)
    structure_buffer: List[str] = []
    consensus_buffer: List[str] = []

    for line in stk_file:
        line = line.strip()
        
        # Skip empty lines
        if not line:
            continue

        if line.startswith("//"):
            # End of a record: Process buffers and yield
            
            # Join lists into final strings
            final_consensus = "".join(consensus_buffer) if consensus_buffer else None
            final_structure = "".join(structure_buffer)
            
            # Convert sequence buffers to strings
            final_sequences = {k: "".join(v) for k, v in sequences.items()}

            yield RnaAlignment(final_sequences, final_structure, consensus=final_consensus)

            # Reset buffers for the next record
            sequences = defaultdict(list)
            structure_buffer = []
            consensus_buffer = []

        elif line.startswith("#=GC SS_cons"):
            # Parse Secondary Structure Consensus
            # Format: #=GC SS_cons <structure>
            parts = line.split()
            if len(parts) >= 3:
                structure_buffer.append(parts[2])

        elif line.startswith("#=GC RF"):
            # Parse Reference / Consensus Sequence
            # Format: #=GC RF <sequence>
            parts = line.split()
            if len(parts) >= 3:
                consensus_buffer.append(parts[2])

        elif not line.startswith("#"):
            # Parse Sequence Data
            # Format: <seq_name> <sequence>
            parts = line.split()
            if len(parts) >= 2:
                seq_id = parts[0]
                seq_data = parts[1]
                sequences[seq_id].append(seq_data)

class RnaAlignment:
    """
    Represents an RNA sequence alignment including structure and consensus information.
    """
    def __init__(self, aligned_sequences: Dict[str, str], structure: str, 
                 consensus: Optional[str] = None, energies: Optional[Any] = None):
        """
        Initialize the RNA Alignment object.

        Args:
            aligned_sequences (Dict[str, str]): Dictionary mapping sequence IDs to aligned strings.
            structure (str): Secondary structure string (e.g., dot-bracket notation).
            consensus (Optional[str]): The consensus sequence string.
            energies (Optional[Any]): Energy scoring data (optional).
        """
        self.aligned_sequences = aligned_sequences
        self.structure = structure
        self.energies = energies
        self.consensus = consensus
        
        # Calculate max length efficiently without creating a full intermediate list
        self.alignment_length = max((len(seq) for seq in aligned_sequences.values()), default=0)

    def add_upstream(self, prefix_sequences: Dict[str, str], prefix_structure: Optional[str] = None) -> 'RnaAlignment':
        """
        Prepends sequences to the aligned sequences.

        Args:
            prefix_sequences (Dict[str, str]): Dictionary of sequences to prepend. 
                                               Keys must match existing sequences.
            prefix_structure (Optional[str]): Structure string corresponding to the prefix. 
                                              If None, ':' padding is used.

        Returns:
            RnaAlignment: A new instance with updated sequences and structure.
        """
        # Combine prefix + existing sequence
        new_sequences = {
            seq_id: prefix_sequences[seq_id] + curr_seq
            for seq_id, curr_seq in self.aligned_sequences.items()
        }

        if prefix_structure is None:
            # Determine padding length based on the longest prefix provided
            max_len = max((len(seq) for seq in prefix_sequences.values()), default=0)
            new_structure = (max_len * ":") + self.structure
        else:
            new_structure = prefix_structure + self.structure

        return RnaAlignment(new_sequences, new_structure)

    def add_downstream(self, suffix_sequences: Dict[str, str], suffix_structure: Optional[str] = None) -> 'RnaAlignment':
        """
        Appends sequences to the aligned sequences.

        Args:
            suffix_sequences (Dict[str, str]): Dictionary of sequences to append.
                                               Keys must match existing sequences.
            suffix_structure (Optional[str]): Structure string corresponding to the suffix.
                                              If None, ':' padding is used.

        Returns:
            RnaAlignment: A new instance with updated sequences and structure.
        """
        # Combine existing sequence + suffix
        new_sequences = {
            seq_id: curr_seq + suffix_sequences[seq_id]
            for seq_id, curr_seq in self.aligned_sequences.items()
        }

        if suffix_structure is None:
            # Determine padding length based on the longest suffix provided
            max_len = max((len(seq) for seq in suffix_sequences.values()), default=0)
            new_structure = self.structure + (max_len * ":")
        else:
            new_structure = self.structure + suffix_structure

        return RnaAlignment(new_sequences, new_structure)

    def add_margin(self, length: int) -> 'RnaAlignment':
        """
        Adds 'N' padding to sequences and ':' padding to structure/consensus on both sides.

        Args:
            length (int): The number of padding characters to add to each side.

        Returns:
            RnaAlignment: A new instance with padded sequences, structure, and consensus.
        """
        seq_padding = "N" * length
        struct_padding = ":" * length

        new_sequences = {
            seq_id: seq_padding + curr_seq + seq_padding
            for seq_id, curr_seq in self.aligned_sequences.items()
        }

        new_structure = struct_padding + self.structure + struct_padding

        new_consensus = None
        if self.consensus is not None:
            new_consensus = seq_padding + self.consensus + seq_padding

        return RnaAlignment(new_sequences, new_structure, consensus=new_consensus)

def invertstrand(strand: str) -> str:
    if strand == "+":
        return "-"
    elif strand == "-":
        return "+"
    return strand

def revcom(sequence: str) -> str:
    # Use translation table for better performance
    complement = str.maketrans({
        'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', 'N': 'N', 
        'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'u': 'a', 'n': 'n',
        '-': '-'
    })
    return sequence.translate(complement)[::-1]

@dataclass(slots=True)
class GenomeRange:
    dbname: str
    chrom: str
    start: int
    end: int
    strand: Optional[str] = None
    name: Optional[str] = None
    orderstrand: bool = False
    data: Optional[Any] = None
    fastafile: Optional[str] = None

    def __post_init__(self):
        self.start = int(self.start)
        self.end = int(self.end)
        if self.orderstrand and self.start > self.end:
            self.start, self.end = self.end, self.start
            self.strand = invertstrand(self.strand)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, GenomeRange):
            return NotImplemented
        return (self.strand == other.strand and 
                self.chrom == other.chrom and 
                self.start == other.start and 
                self.end == other.end)

    def __hash__(self) -> int:
        return self.start + self.end + hash(self.chrom) + hash(self.strand)

    def coverage(self, other: 'GenomeRange') -> int:
        if self.strand == other.strand and self.chrom == other.chrom:
            start = max(self.start, other.start)
            end = min(self.end, other.end)
            return max(0, end - start)
        return 0
            
    def bedstring(self, name: Optional[str] = None, score: int = 1000) -> str:
        strand = self.strand if self.strand is not None else "+"
        
        if name is None:
            name = self.name if self.name is not None else "FEAT"
            
        return "\t".join([self.chrom, str(self.start), str(self.end), name, str(score), strand])

    def length(self) -> int:
        return self.end - self.start

    def addmargin(self, dist: int = 50, name: Optional[str] = None) -> 'GenomeRange':
        newname = name if name is not None else self.name
        return GenomeRange(self.dbname, self.chrom, self.start - dist, self.end + dist, self.strand, name=newname, fastafile=self.fastafile)

    def getupstream(self, dist: int = 50, name: Optional[str] = None) -> 'GenomeRange':
        newname = name if name is not None else self.name
        return GenomeRange(self.dbname, self.chrom, self.start - dist, self.start, self.strand, name=newname, fastafile=self.fastafile)

    def getdownstream(self, dist: int = 50, name: Optional[str] = None) -> 'GenomeRange':
        newname = name if name is not None else self.name
        if self.strand != '-':
            return GenomeRange(self.dbname, self.chrom, self.end, self.end + dist, self.strand, name=newname, fastafile=self.fastafile)
        else:
            return GenomeRange(self.dbname, self.chrom, self.start - dist, self.start, self.strand, name=newname, fastafile=self.fastafile)

    def getfirst(self, dist: int = 50, name: Optional[str] = None) -> 'GenomeRange':
        newname = name if name is not None else self.name
        if self.strand == "-":
            return GenomeRange(self.dbname, self.chrom, self.end - dist, self.end, self.strand, name=newname, fastafile=self.fastafile)
        else:
            return GenomeRange(self.dbname, self.chrom, self.start, self.start + dist, self.strand, name=newname, fastafile=self.fastafile)

    def getlast(self, dist: int = 50, name: Optional[str] = None) -> 'GenomeRange':
        newname = name if name is not None else self.name
        if self.strand == "-":
            return GenomeRange(self.dbname, self.chrom, self.start, self.start + dist, self.strand, name=newname, fastafile=self.fastafile)
        else:
            return GenomeRange(self.dbname, self.chrom, self.end - dist, self.end, self.strand, name=newname, fastafile=self.fastafile)

    def getbase(self, basenum: int, name: Optional[str] = None) -> 'GenomeRange':
        newname = name if name is not None else self.name
        if self.strand == "-":
            return GenomeRange(self.dbname, self.chrom, self.end - basenum - 1, self.end - basenum, self.strand, name=newname, fastafile=self.fastafile)
        else:
            return GenomeRange(self.dbname, self.chrom, self.start + basenum, self.start + basenum + 1, self.strand, name=newname, fastafile=self.fastafile)

    def antisense(self) -> 'GenomeRange':
        newstrand = "-" if self.strand == "+" else "+"
        return GenomeRange(self.dbname, self.chrom, self.start, self.end, newstrand, name=self.name, fastafile=self.fastafile)

    def bamseq(self) -> str:
        if self.data is None or "seq" not in self.data:
             return ""
        if self.strand == "+":
            return self.data["seq"]
        else:
            return revcom(self.data["seq"])

    def getgc(self) -> int:
        seq = self.bamseq()
        return sum(1 if currbase in {"G", "C"} else 0 for currbase in seq)

@dataclass(slots=True, unsafe_hash=True)
class tRNAlocus:
    loc: GenomeRange
    seq: str
    score: float
    amino: str
    anticodon: str
    intronseq: str
    intron: Optional[Tuple[int, int]]
    rawseq: Optional[str] = None
    name: Optional[str] = field(init=False)

    def __post_init__(self):
        self.name = self.loc.name

@dataclass(slots=True)
class tRNAtranscript:
    seq: str
    score: Union[float, set]
    amino: str
    anticodon: str
    loci: Union[Tuple[tRNAlocus, ...], set, list]
    intronseqs: Union[str, set]
    name: Optional[str] = None
    rawseq: Optional[str] = None
    artificialtrna: bool = False

    def __post_init__(self):
        if isinstance(self.loci, (set, list)):
            self.loci = tuple(self.loci)

    def getmatureseq(self, addcca: bool = True) -> str:
        prefix = ""
        if self.amino == "His":
            prefix = "G"
        end = ""
        if addcca and not self.artificialtrna:
            end = "CCA"
        return prefix + self.seq + end 

def getuniquetRNAs(trnalist: List[tRNAlocus]) -> Generator[tRNAtranscript, None, None]:
    sequencedict: Dict[str, List[tRNAlocus]] = defaultdict(list)
    
    for curr in trnalist:
        sequencedict[curr.seq].append(curr)
        
    for currtrans, loci_list in sequencedict.items():
        scores = {curr.score for curr in loci_list}
        anticodon = {curr.anticodon for curr in loci_list}
        amino = {curr.amino for curr in loci_list}
        introns = {curr.intronseq for curr in loci_list}
        
        #remove psuedo if there's a real one somewheres
        if len(anticodon) > 1:
            anticodon.discard('Xxx')
            amino.discard('X')
        if introns == {""}:
            introns = set()
            
        if len(scores) > 1:
            # Multiple scores found
            pass
        if len(anticodon) > 1:
            logger.error("tRNA file contains identical tRNAs with seperate anticodons, cannot continue")
            sys.exit(1)
            
        yield tRNAtranscript(currtrans, scores, list(amino)[0], list(anticodon)[0], loci_list, introns)

def striplocus(trnaname: str) -> str:
    return re.sub(r"\-\d+$", "", trnaname)

def readtRNAscan(scanfile: Union[str, TextIO], genomefile: str, mode: Optional[str] = None) -> Generator[tRNAlocus, None, None]:
    trnalist: List[GenomeRange] = []
    orgname = "genome"
    
    if hasattr(scanfile, 'read'):
        trnascan = scanfile
    else:
        trnascan = open(scanfile, "r")
        
    trnascore: Dict[str, float] = {}
    trnaanticodon: Dict[str, str] = {}
    trnaamino: Dict[str, str] = {}
    tRNAintron: Dict[str, Tuple[int, int]] = {}
    trnas: Dict[str, GenomeRange] = {}
    
    try:
        for currline in trnascan:
            if currline.startswith("Sequence") or currline.startswith("Name") or currline.startswith("------"):
                continue
            
            fields = currline.split()
            
            if mode == "gtRNAdb":
                del fields[6:8]
                
            curramino = fields[4]
            currac = fields[5]
            
            if currac == "???":
                currac = 'Xxx'
                
            start = int(fields[2])
            end = int(fields[3])
            
            if start > end:
                end -= 1
            else:
                start -= 1
                
            currchrom = fields[0]
            trnanum = fields[1]
            
            # Create GenomeRange
            currtRNA = GenomeRange(
                dbname=orgname, 
                chrom=currchrom, 
                start=start, 
                end=end, 
                name=f"{currchrom}.tRNA{trnanum}-{curramino}{currac}", 
                strand="+", 
                orderstrand=True,
                fastafile=genomefile
            )
            
            currtrans = currtRNA

            trnaamino[currtrans.name] = curramino
            trnaanticodon[currtrans.name] = currac
            trnascore[currtrans.name] = float(fields[8])
            trnas[currtrans.name] = currtrans
        
            trnalist.append(currtRNA)
            
            if int(fields[6]) != 0:
                if currtRNA.strand == "-":
                    intronstart = int(fields[2]) - int(fields[6]) 
                    intronend = int(fields[2]) - int(fields[7]) + 1
                else:
                    intronstart = int(fields[6]) - int(fields[2]) - 1
                    intronend = int(fields[7]) - int(fields[2])
                tRNAintron[currtRNA.name] = (intronstart, intronend)
    finally:
        if not hasattr(scanfile, 'read'):
            trnascan.close()
            
    #should add check to make sure all trnas are grabbed
    trnaseqs = getseqdict(trnalist, faifiles={orgname: genomefile + ".fai"})
    intronseq: Dict[str, str] = defaultdict(str)
    
    for curr in trnaseqs.keys():
        currintron = None
        if curr in tRNAintron:
            start_intron, end_intron = tRNAintron[curr]
            intronseq[curr] = trnaseqs[curr][start_intron:end_intron]
            trnaseqs[curr] = trnaseqs[curr][:start_intron] + trnaseqs[curr][end_intron:]
            currintron = tRNAintron[curr]
            
        yield tRNAlocus(trnas[curr], trnaseqs[curr], trnascore[curr], trnaamino[curr], trnaanticodon[curr], intronseq[curr], currintron)

def readtRNAdb(scanfile: Union[str, TextIO], genomefile: str, trnamap: Dict[str, str]) -> Generator[tRNAtranscript, None, None]:
    trnalist: List[GenomeRange] = []
    orgname = "genome"
    
    if hasattr(scanfile, 'read'):
        trnascan = scanfile
    else:
        trnascan = open(scanfile, "r")
        
    trnascore: Dict[str, float] = {}
    trnaanticodon: Dict[str, str] = {}
    trnaamino: Dict[str, str] = {}
    tRNAintron: Dict[str, Tuple[int, int]] = {}
    trnas: Dict[str, GenomeRange] = {}
    
    try:
        for linenum, currline in enumerate(trnascan):
            if currline.startswith("Sequence") or currline.startswith("Name") or currline.startswith("------"):
                continue
            if len(currline) < 5:
                logger.warning(f"cannot read line: {linenum} of {scanfile}")
                continue
                
            fields = currline.split()
            curramino = fields[4]
            currac = fields[5]
            
            start = int(fields[2])
            end = int(fields[3])
            
            if start > end:
                end -= 1
            else:
                start -= 1
                
            currchrom = fields[0]
            trnanum = fields[1]
            
            trnascanname = f"{currchrom}.trna{trnanum}-{curramino}{currac}"
            shorttrnascanname = f"{currchrom}.trna{trnanum}"
            
            if trnascanname in trnamap:
                currtRNA = GenomeRange(orgname, currchrom, start, end, name=trnamap[trnascanname], strand="+", orderstrand=True)
            elif shorttrnascanname in trnamap:
                currtRNA = GenomeRange(orgname, currchrom, start, end, name=trnamap[shorttrnascanname], strand="+", orderstrand=True)
            else:
                logger.warning(f"Skipping {trnascanname}, has no transcript name")
                continue
                
            currtrans = currtRNA

            trnaamino[currtrans.name] = curramino
            trnaanticodon[currtrans.name] = currac
            trnascore[currtrans.name] = float(fields[8])
            trnas[currtrans.name] = currtrans
        
            currtRNA.fastafile = genomefile
            trnalist.append(currtRNA)
            
            if int(fields[6]) != 0:
                if currtRNA.strand == "-":
                    intronstart = int(fields[2]) - int(fields[6])
                    intronend = int(fields[2]) - int(fields[7]) + 1
                else:
                    intronstart = int(fields[6]) - int(fields[2]) 
                    intronend = int(fields[7]) - int(fields[2]) + 1
                tRNAintron[currtRNA.name] = (intronstart, intronend)
    finally:
        if not hasattr(scanfile, 'read'):
            trnascan.close()
            
    trnaseqs = getseqdict(trnalist, faifiles={orgname: genomefile + ".fai"})
    intronseq: Dict[str, str] = defaultdict(str)
    trnaloci: List[tRNAlocus] = []
    
    for curr in trnaseqs.keys():
        currintron = None
        if curr in tRNAintron:
            start_intron, end_intron = tRNAintron[curr]
            intronseq[curr] = trnaseqs[curr][start_intron:end_intron]
            trnaseqs[curr] = trnaseqs[curr][:start_intron] + trnaseqs[curr][end_intron:]
            currintron = tRNAintron[curr]

        trnaloci.append(tRNAlocus(trnas[curr], trnaseqs[curr], trnascore[curr], trnaamino[curr], trnaanticodon[curr], intronseq[curr], currintron))
    
    trnaloci.sort(key=lambda x: striplocus(x.name))
    
    for transname, currloci in itertools.groupby(trnaloci, lambda x: striplocus(x.name)):
        currlocus = list(currloci)
        yield tRNAtranscript(
            currlocus[0].seq, 
            {curr.score for curr in currlocus}, 
            currlocus[0].amino, 
            currlocus[0].anticodon, 
            set(currlocus), 
            currlocus[0].intronseq, 
            name=transname
        )

def tempmultifasta(allseqs) -> tempfile.NamedTemporaryFile:
    fafile = tempfile.NamedTemporaryFile(suffix=".fa", mode="w+")

    for seqname, seq in allseqs:
        print(f">{seqname}\n{seq}", file=fafile)

    fafile.flush()
    return fafile

def getseqdict(genelist, faifiles = None) -> Dict[str, str]:
    allorgs = {currgene.dbname for currgene in genelist}
    dbdict = {currorg: {} for currorg in allorgs}
    fastafiles = {}
    
    for currgene in genelist:
        dbdict[currgene.dbname][currgene.name] = currgene
        if currgene.fastafile is not None:
            fastafiles[currgene.dbname] = currgene.fastafile
        else:
            pass 
            
    seqdict = {}
    
    for currorg in allorgs:
        if faifiles is not None and currorg in faifiles:
            currseqs = getseqs(fastafiles[currorg], dbdict[currorg], faindex=faifiles[currorg])
        else:
            currseqs = getseqs(fastafiles[currorg], dbdict[currorg])
        seqdict.update(currseqs)
    return seqdict

def getnamedict(genelist) -> Dict[str, GenomeRange]:
    return {currgene.name: currgene for currgene in genelist}

class FastqSeqException(Exception):
    pass

def getseqs(fafile: str, rangedict: Dict[str, GenomeRange], faindex: Optional[str] = None) -> Dict[str, str]:
    if faindex is not None:
        try:
            faifile_obj = fastaindex(fafile, faindex)
        except IOError:
            logger.error(f"Cannot read fasta file {fafile}")
            logger.error(f"Ensure that file {fafile} exits and generate fastaindex {faindex} with samtools faidx")
            sys.exit(1)
        return faifile_obj.getseqs(rangedict)
    
    with open(fafile, "r") as genomefile:
        reheader = re.compile(r"\>([^]+)")
        allseqs: Dict[str, str] = defaultdict(str)
        currloc = 0
        currseq = ""
        
        for line in genomefile:
            line = line.rstrip("\n")
            currheader = reheader.match(line)
            if currheader:
                currseq = currheader.group(1)
                currloc = 0
            else:
                line_len = len(line)
                for currname, location in rangedict.items():
                    if currseq == location.chrom:
                        chromstart = location.start
                        chromend = location.end
                        
                        overlap = False
                        if currloc <= chromstart < currloc + line_len:
                            overlap = True
                        elif currloc < chromend <= currloc + line_len:
                            overlap = True
                        elif chromstart < currloc and chromend > currloc + line_len:
                            overlap = True
                            
                        if overlap:
                            start_in_line = max(0, chromstart - currloc)
                            end_in_line = min(line_len, chromend - currloc)
                            
                            if start_in_line < end_in_line:
                                allseqs[currname] += line[start_in_line:end_in_line]
                    
                currloc += line_len

    finalseqs: Dict[str, str] = {}
    
    comp_trans_full = str.maketrans({
        'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', 'N': 'N', 
        'R': 'Y', 'Y': 'R', 'S': 'W', 'W': 'S', 'K': 'M', 'M': 'K'
    })

    for currname, seq_str in allseqs.items():
        if rangedict[currname].strand == "-":
            finalseqs[currname] = seq_str.upper().translate(comp_trans_full)[::-1]
        else:
            finalseqs[currname] = seq_str.upper()
            
    for currseq in rangedict.keys():
        if currseq not in finalseqs:
            logger.warning(f"No sequence extracted for {rangedict[currseq].dbname}.{rangedict[currseq].chrom}:{rangedict[currseq].start}-{rangedict[currseq].end}")
            
    return finalseqs        

class fastaindex:
    def __init__(self, fafile: str, faifile: str):
        self.fafile = fafile
        self.chromsize: Dict[str, int] = {}
        self.chromoffset: Dict[str, int] = {}
        self.seqlinesize: Dict[str, int] = {}
        self.seqlinebytes: Dict[str, int] = {}
        
        with open(faifile) as fai:
            for line in fai:
                fields = line.split("\t")
                self.chromsize[fields[0]] = int(fields[1])
                self.chromoffset[fields[0]] = int(fields[2])
                self.seqlinesize[fields[0]] = int(fields[3])
                self.seqlinebytes[fields[0]] = int(fields[4])

    def getchrombed(self, dbname: str = 'genome') -> Generator[GenomeRange, None, None]:
        for curr in self.chromsize.keys():
            yield GenomeRange(dbname, curr, 0, self.chromsize[curr], name=curr, strand="+")

    def getseek(self, currchrom: str, loc: int) -> int:
        if currchrom not in self.chromsize:
            raise FastqSeqException(f"sequence {currchrom} not found in index for {self.fafile}")
        
        return self.chromoffset[currchrom] + loc + int(loc / self.seqlinesize[currchrom]) * (self.seqlinebytes[currchrom] - self.seqlinesize[currchrom])

    def getfullseqs(self, names: List[str]) -> Generator[Tuple[str, str], None, None]:
        with open(self.fafile, "r") as genomefile:
            for currchrom in names:
                genomefile.seek(self.getseek(currchrom, 0))
                seq = genomefile.read(self.getseek(currchrom, self.chromsize[currchrom]) - self.getseek(currchrom, 0))
                yield currchrom, seq 

    def getseqs(self, rangedict: Dict[str, GenomeRange]) -> Dict[str, str]:
        allseqs: Dict[str, Optional[str]] = {}
        
        with open(self.fafile, "r") as genomefile:
            for currname, currregion in rangedict.items():
                try:
                    currchrom = currregion.chrom
                    start_seek = self.getseek(currchrom, currregion.start)
                    end_seek = self.getseek(currchrom, currregion.end)
                    
                    genomefile.seek(start_seek)
                    seq = genomefile.read(end_seek - start_seek)
                    seq = seq.replace("\n", "")
                    allseqs[currname] = seq
                except FastqSeqException:
                    allseqs[currname] = None
                    pass
        
        finalseqs: Dict[str, str] = {}
        
        comp_trans_full = str.maketrans({
            'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', 'N': 'N', 
            'R': 'Y', 'Y': 'R', 'S': 'W', 'W': 'S', 'K': 'M', 'M': 'K'
        })

        for currname, seq_val in allseqs.items():
            if seq_val is None:
                continue
            
            if rangedict[currname].strand == "-":
                finalseqs[currname] = seq_val.upper().translate(comp_trans_full)[::-1]
            else:
                finalseqs[currname] = seq_val.upper()
                
        return finalseqs    

@dataclass(slots=True)
class GenomeRead(GenomeRange):
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, GenomeRead):
            return NotImplemented
        return self.name == other.name
    
    def __hash__(self) -> int:
        return hash(self.name)

@dataclass(slots=True)
class BamRead(GenomeRange):
    bamline: Any = None

    def getmismatches(self) -> Optional[int]:
        if self.bamline.has_tag("XM"):
            return int(self.bamline.get_tag("XM"))
        else:
            return None
    def isuniqueaminomapping(self) -> Optional[bool]:
        if self.bamline.has_tag("YM"):
            return not self.bamline.get_tag("YM") > 1
        else:
            return None
    def isuniqueacmapping(self) -> Optional[bool]:
        if self.bamline.has_tag("YA"):
            return not self.bamline.get_tag("YA") > 1
        else:
            return None
    def isuniquetrnamapping(self) -> Optional[bool]:
        if self.bamline.has_tag("YR"):
            return not self.bamline.get_tag("YR") > 1
        else:
            return None

    def hasindel(self) -> bool:
        return len(self.bamline.cigar) > 1
    def getlength(self) -> int:
        return len(self.bamline.seq)
    def getseq(self) -> str:
        if self.strand == '-':
            return revcom(self.bamline.seq)
        else:
            return self.bamline.seq
            
    def getgc(self) -> int:
        seq = self.getseq()
        return sum(1 if currbase in set(["G","C"]) else 0 for currbase in seq) 
    def issinglemapped(self) -> bool:
        return self.bamline.mapq >= 2
    def getcigar(self) -> List[Tuple[int, int]]:
        return self.bamline.cigar

def isprimarymapping(mapping: Any) -> bool:
    return not (mapping.flag & 0x0100 > 0)        
def issinglemapping(mapping: Any) -> bool:
    return mapping.mapq > 2
    
def mismatchnum(mapping: Any) -> int:
    return int(mapping.get_tag("XM"))

def getbamrangeshortseq(bamfile: Any, chromrange: Optional[GenomeRange] = None, primaryonly: bool = False, singleonly: bool = False, maxmismatches: Optional[int] = None, allowindels: bool = True, skiptags: bool = False) -> Generator[GenomeRead, None, None]:
    try:
        if chromrange is not None:
            bamiter = bamfile.fetch(chromrange.chrom, chromrange.start, chromrange.end)
        else:
            bamiter = bamfile.fetch()
   
        for currline in bamiter: 
            if primaryonly and not isprimarymapping(currline):
                continue

            if singleonly and not issinglemapping(currline):
                continue
            if not allowindels and len(currline.cigar) > 1:
                continue
            if maxmismatches is not None and mismatchnum(currline) > maxmismatches:
                continue
            rname = bamfile.getrname(currline.rname)
            strand = ifelse(currline.is_reverse, '-','+')
            yield GenomeRead( "genome",rname,currline.pos,currline.aend,strand, name = currline.qname, data = {"seq":currline.seq} )

    except ValueError:
        if chromrange is not None:
            pass

def getbamrangeshort(bamfile: Any, chromrange: Optional[GenomeRange] = None, primaryonly: bool = False, singleonly: bool = False, maxmismatches: Optional[int] = None, allowindels: bool = True, skiptags: bool = False) -> Generator[GenomeRead, None, None]:
    try:
        if chromrange is not None:
            bamiter = bamfile.fetch(chromrange.chrom, chromrange.start, chromrange.end)
        else:
            bamiter = bamfile.fetch()
   
        for currline in bamiter:
            if primaryonly and not isprimarymapping(currline):
                continue

            if singleonly and not issinglemapping(currline):
                continue
            if not allowindels and len(currline.cigar) > 1:
                continue
            if maxmismatches is not None and mismatchnum(currline) > maxmismatches:
                continue
            rname = bamfile.getrname(currline.rname)
            strand = ifelse(currline.is_reverse, '-','+')
            
            yield GenomeRead( "genome",rname,currline.pos,currline.aend,strand, name = currline.qname)

    except ValueError:
        if chromrange is not None:
            pass

def getbamrange(bamfile: Any, chromrange: Optional[GenomeRange] = None, primaryonly: bool = False, singleonly: bool = False, maxmismatches: Optional[int] = None, allowindels: bool = True, skiptags: bool = False) -> Generator[GenomeRead, None, None]:
    try:
        if chromrange is not None:
            bamiter = bamfile.fetch(chromrange.chrom, chromrange.start, chromrange.end)
        else:
            bamiter = bamfile.fetch()
   
        for currline in bamiter:
            if primaryonly and not isprimarymapping(currline):
                continue

            if singleonly and not issinglemapping(currline):
                continue
            if maxmismatches is not None and mismatchnum(currline) > maxmismatches:
                continue
            rname = bamfile.getrname(currline.rname)
            strand = ifelse(currline.is_reverse, '-','+')
            
            seq = currline.seq
            
            uniqueac = True
            uniqueamino = True
            uniquetrna = True
            # mismatches = None
            alignscore = None
            secondbestscore = None
            uniquemapping = False

            if not skiptags:
                if currline.has_tag("YA") and currline.get_tag("YA") > 1:
                    uniqueac = False
                if currline.has_tag("YM") and currline.get_tag("YM") > 1:
                    uniqueamino = False
                if currline.has_tag("YR") and currline.get_tag("YR") > 1:
                    uniquetrna = False
                if currline.has_tag("XM"):
                    pass # mismatches = currline.get_tag("XM")
                if currline.has_tag("XS"):
                    secondbestscore = float(currline.get_tag("XS"))
                if currline.has_tag("AS"):
                    alignscore = float(currline.get_tag("AS"))

            if secondbestscore is None or alignscore > secondbestscore:
                uniquemapping = True
            
            if not allowindels and len(currline.cigar) > 1:
                continue
                
            yield GenomeRead( "genome",rname,currline.pos,currline.aend,strand, name = currline.qname , data = {"score":currline.mapq, "CIGAR":currline.cigar,"CIGARstring":currline.cigarstring, "seq":seq, "flags": currline.flag, "qual":currline.qual,"bamline":currline,'uniqueac':uniqueac,"uniqueamino":uniqueamino,"uniquetrna":uniquetrna,"uniquemapping":uniquemapping})
    except ValueError:
        if chromrange is not None:
            pass

def getbam(bamfile: Any, chromrange: Optional[GenomeRange] = None, primaryonly: bool = False, singleonly: bool = False, allowindels: bool = True) -> Generator[BamRead, None, None]:
    try:
        if chromrange is not None:
            bamiter = bamfile.fetch(chromrange.chrom, chromrange.start, chromrange.end)
        else:
            bamiter = bamfile.fetch()
   
        for currline in bamiter:
            if primaryonly and not isprimarymapping(currline):
                continue

            if singleonly and not issinglemapping(currline):
                continue
           
            rname = bamfile.getrname(currline.rname)
            strand = ifelse(currline.is_reverse, '-','+')

            if currline.aend is None:
                continue
            if not allowindels and len(currline.cigar) > 1:
                continue
            yield BamRead("genome", rname, currline.pos, currline.aend, strand=strand, name=currline.qname, bamline=currline)
    except ValueError:
        if chromrange is not None:
            pass

def isuniquetrnamapping(read: GenomeRead) -> bool:
    return read.data["uniquetrna"]
def isuniqueaminomapping(read: GenomeRead) -> bool:
    return read.data["uniqueamino"]
def isuniqueacmapping(read: GenomeRead) -> bool:
    return read.data["uniqueac"]
def issinglemapped(read: GenomeRead) -> bool:
    return (read.data["score"] >= 2)  

def getfragtype(currfeat: GenomeRange, currread: GenomeRange, maxoffset: int = 10) -> str:
    if currread.start < currfeat.start + maxoffset and currread.end > currfeat.end - maxoffset:
        return "Whole"
    elif currread.start < currfeat.start + maxoffset:
        if currfeat.strand == "+":
            return "Fiveprime"
        else:
            return "Threeprime"
    elif currread.end > currfeat.end - maxoffset:
        if currfeat.strand == "+":
            return "Threeprime"
        else:
            return "Fiveprime"
    return "Other"

def getendtype(currfeat: GenomeRange, currread: GenomeRange, maxoffset: int = 10) -> Optional[str]:
    endtype = None
    if currread.end == currfeat.end:
        endtype = "CCA"
    elif currread.end == currfeat.end -1:
        endtype = "CC"
    elif currread.end == currfeat.end -2:
        endtype = "C"
    elif currread.end == currfeat.end -3:
        endtype = ""
    return endtype

class RangeBin:     
    def __init__(self, rangelist: Iterable[GenomeRange], binfactor: int = 10000):
        self.binfactor = binfactor
        self.bins: List[set] = []
        self.length = 0
        for curr in rangelist:
            self.additem(curr)
        
    def __len__(self) -> int:
        return self.length
    def __iter__(self) -> Generator[GenomeRange, None, None]:
        for currbin in self.bins:
            for currgene in currbin:
                yield currgene
    def additem(self, item: GenomeRange):
        binstart = int(item.start / self.binfactor)
        binend = int(item.end / self.binfactor) + 1
        while (binend + 2 >= len(self.bins)):
            self.bins.append(set())
        self.bins[binstart].add(item)
        for i in range(binstart, binend):
            self.bins[i].add(item)
        
        self.length += 1
        
    def getrange(self, item: GenomeRange) -> Generator[GenomeRange, None, None]:
        for i in range(int(item.start / self.binfactor)-1,int(item.end / self.binfactor)+1):
            if i < len(self.bins):
                for currrange in self.bins[i]:
                    if currrange.start >= item.start and currrange.end <= item.end:
                            yield currrange
                            
    def getbin(self, item: GenomeRange) -> Generator[GenomeRange, None, None]:
        outputregions = set()
        for i in range(int(item.start / self.binfactor)-1,int(item.end / self.binfactor)+1):
            if i < len(self.bins) and i >= 0:
                for currrange in self.bins[i]:
                    regionid = (currrange.name, currrange.start, currrange.end)
                    if regionid not in outputregions:
                        outputregions.add(regionid)
                        yield currrange
                    
    def getbinpos(self, item: int) -> Generator[GenomeRange, None, None]:
        for i in range(int(item / self.binfactor)-1,int(item / self.binfactor)+1):
            if i < len(self.bins) and i >= 0:
                for currrange in self.bins[i]:
                    yield currrange
                    
    def getfeatbin(self, name: str) -> Generator[int, None, None]:
        for i, currbin in enumerate(self.bins):
            for currgene in currbin:
                if currgene.name == name:
                    yield i
    def getbinnums(self, item: GenomeRange) -> Generator[int, None, None]:
        for i in range(int(item.start / self.binfactor)-1,int(item.end / self.binfactor)+1):
            if i < len(self.bins) and i >= 0:
                yield i

def uniqueorder(inp: Iterable[Any]) -> Generator[Any, None, None]:
    alldata = set()
    for curr in inp:
        if curr not in alldata:
            yield curr
            alldata.add(curr)

class transcriptfile:
    def __init__(self, trnafilename: str):
        self.locustranscript: Dict[str, str] = {}
        self.transcripts: List[str] = []
        self.loci: List[str] = []
        self.amino: Dict[str, str] = {}
        self.anticodon: Dict[str, str] = {}
        self.transcriptdict: Dict[str, set] = defaultdict(set)
        aminoorder: List[str] = []
        anticodonorder: List[str] = []
        
        with open(trnafilename, 'r') as trnafile:
            for i, line in enumerate(trnafile):
                fields = line.split()
                if len(fields) < 2:
                    continue
                self.transcripts.append(fields[0])
                self.amino[fields[0]] = fields[2]
                self.anticodon[fields[0]] = fields[3]
                aminoorder.append(fields[2])
                anticodonorder.append(fields[3])
                for currlocus in fields[1].split(','):
                    self.locustranscript[currlocus] = fields[0]
                    self.loci.append(currlocus)
                    self.transcriptdict[fields[0]].add(currlocus)

        self.aminoorder = tuple(uniqueorder(aminoorder))             
        self.anticodonorder = tuple(uniqueorder(anticodonorder))

    def gettranscripts(self) -> set:
        return set(self.transcripts)
    def getlocustranscript(self, locus: str) -> str:
        return self.locustranscript[locus]
    def getloci(self) -> List[str]:
        return self.loci
    def getamino(self, trna: str) -> str:
        return self.amino[trna]
    def getanticodon(self, trna: str) -> str:
        return self.anticodon[trna]
        
    def allaminos(self) -> Tuple[str, ...]:
        return self.aminoorder
    def allanticodons(self) -> Tuple[str, ...]:
        return self.anticodonorder
    def getaminotranscripts(self, trnaamino: str) -> set:
        return set(curr for curr in self.transcripts if trnaamino == self.amino[curr])
    def getanticodontranscripts(self, trnaanticodon: str) -> set:
        return set(curr for curr in self.transcripts if trnaanticodon == self.anticodon[curr])

def getpairfile(pairfilename: str) -> Generator[Tuple[str, str], None, None]:
    with open(pairfilename, 'r') as pairfile:
        for currline in pairfile:
            fields = currline.split()
            if len(fields) > 1:
                yield fields[0], fields[1]

class extraseqfile:
    def __init__(self, extraseqfilename: str):
        self.logger = logging.getLogger(__name__)
        self.seqlist: List[str] = []
        self.seqfasta: Dict[str, str] = {}
        self.seqbed: Dict[str, str] = {}
        
        if extraseqfilename is None or not os.path.exists(extraseqfilename):
            return

        try:
            with open(extraseqfilename, 'r') as extrafile:
                directory = os.path.dirname(extraseqfilename)
                for i, line in enumerate(extrafile):
                    fields = line.split()
                    if len(fields) < 2:
                        continue
                    self.seqfasta[fields[0]] = directory+"/"+fields[1]
                    self.seqbed[fields[0]] = directory+"/"+fields[2]
                    self.seqlist.append(fields[0])
        except IOError as e:
            self.seqlist = []
            self.seqfasta = {}
            self.seqbed = {}
            self.logger.error(f"extraseqfile I/O error({e.errno}): {e.strerror}")

    def getseqnames(self) -> Dict[str, set]:
        seqnamedict = defaultdict(set)
        for currseq in self.seqlist:
            seqnamedict[currseq] = set(curr.name for curr in readbed(self.seqbed[currseq]))
        return seqnamedict
    def getseqbeds(self) -> Dict[str, str]:
        return self.seqbed

def getsizefactors(sizefactorfilename: str) -> Dict[str, float]:
    try:
        with open(sizefactorfilename, 'r') as sizefactorfile:
            lines = sizefactorfile.readlines()
    except IOError:
        logger.error(f"Cannot read size factor file {sizefactorfilename}")
        logger.error("check Rlog.txt")
        sys.exit(1)
        
    sizefactors = {}
    bamheaders = []
    sizes = []
    
    for i, line in enumerate(lines):
        if i == 0:
            bamheaders = [curr.strip("\"\n") for curr in line.split()]
        elif i == 1:
            sizes = [float(curr.strip("\"\n")) for curr in line.split()]
            
    for i in range(len(bamheaders)):
        sizefactors[bamheaders[i]] = sizes[i]
    return sizefactors

def ifelse(arg: bool, trueres: Any, falseres: Any) -> Any:
    if arg:
        return trueres
    else:
        return falseres

class samplefile:
    def __init__(self, samplefilename, bamdir = "./"):
        self.logger = logging.getLogger(__name__)
        try:
            samplefile = open(samplefilename)
            samplelist = list()
            samplefiles = dict()
            replicatename = dict()
            
            
            replicatelist = list()
            for i, line in enumerate(samplefile):
                fields = line.split()
                if len(fields) < 3:
                    continue
                # Skip header if present
                if i == 0 and fields[0].lower() == 'fastq':
                    continue
                    
                # New format: fastq (0), sample (1), group (2)
                sample_name = fields[1]
                fastq_path = fields[0]
                group_name = fields[2]
                
                samplefiles[sample_name] = fastq_path
                replicatename[sample_name] = group_name
                
                samplelist.append(sample_name)
                if group_name not in set(replicatelist):
                    replicatelist.append(group_name)
            
            #bamlist = list(curr + "_sort.bam" for curr in samplefiles.iterkeys())
            #samplenames = list(curr  for curr in samplefiles.keys())
            if bamdir is None:
                bamdir = "./"
            self.bamdir = bamdir
            self.samplelist = samplelist
            self.samplefiles = samplefiles
            self.replicatename = replicatename
            self.replicatelist = replicatelist
            #self.bamlist = list(curr+ "_sort.bam" for curr in samplelist)
        except IOError as e:
            self.logger.error(f"Cannot read sample file {samplefilename}: {e}")
            raise
    def getsamples(self):
        return self.samplelist
    def getbamlist(self):
        return list(curr+ ".bam" for curr in self.samplelist)
    def getbam(self, sample):
        return self.bamdir + "/" + sample + ".bam" 
    def getmergebam(self, sample):
        return self.bamdir + "/" + sample + "-merge.bam" 
    def getfastq(self, sample):
        return self.samplefiles[sample]
    def getreplicatename(self, sample):
        return self.replicatename[sample]
    def allreplicates(self):
        return self.replicatelist
    def getrepsamples(self, replicate):
        return list(currsample for currsample in self.samplelist if self.replicatename[currsample] == replicate)
    def getfastqsample(self, fastqname):
        for currsample in self.samplefiles.keys():
            if self.samplefiles[currsample] == fastqname:
                return currsample

def readfeatures(filename, orgdb="genome", seqfile= None, removepseudo = False):
    if filename.endswith(".bed") or filename.endswith(".bed.gz"):
        return readbed(filename, orgdb, seqfile)
    elif filename.endswith(".gtf") or filename.endswith(".gtf.gz") or filename.endswith(".gff") or filename.endswith(".gff.gz"):
        #print >>sys.stderr, removepseudo
        return (curr for curr in readgtf(filename, orgdb, seqfile, filterpsuedo = removepseudo, filtertypes =set(['retained_intron','antisense','lincRNA']) ))
    else:
        logger.error(filename+" not valid feature file")
        sys.exit()


def readgtf(filename, orgdb="genome", seqfile= None, filterpsuedo = False, replacename = False, filtertypes = set(['retained_intron','antisense','lincRNA'])):
    bedfile = None
    skippedlines = 0
    #print >>sys.stderr, "****"
    if filename == "stdin":
        bedfile = sys.stdin
    elif filename.endswith(".gz"):
        bedfile = gzip.open(filename, 'rt', encoding='utf-8')
    else:
        bedfile = open(filename, "r")
    
    for currline in bedfile:
        #print currline
        if currline.startswith('track') or currline.startswith('#'):
            continue
        fields = currline.rstrip().split("\t")
        if len(fields) > 2:
            biotype = "Unknown"
            featname = None
            genename = None
            #print >>sys.stderr, len(fields)
            genesource = fields[1]  
            #retained introns are often other things as well, so I skip em
            if fields[2] != "transcript" or genesource in filtertypes:
                continue

              
            for currattr in fields[8].rstrip(";").split(";"):
                #print >>sys.stderr,  currattr
                currname = currattr.strip().split()[0]
                currvalue = currattr.strip().split()[1]
                if currname == "transcript_name":
                    featname = currvalue.strip('"')
                elif currname == "name" or currname == "gene_id" and featname is None:
                    featname = currvalue.strip('"')

                elif currname == "gene_biotype":
                    biotype = currvalue.strip('"')
                elif currname == "gene_name":
                    genename = currvalue.strip('"')
                #print >>sys.stderr, "***||"
            if genename is None:
                #print >>sys.stderr, "No name for gtf entry "+featname
                genename = featname
                pass
            if filterpsuedo and biotype == "pseudogene":
                #print >>sys.stderr, "*******"
                continue
                
            if genesource == 'ensembl':
                #print >>sys.stderr, biotype
                genesource = biotype
            if not (fields[6] == "+" or fields[6] == "-"):
                logger.warning("strand error in "+filename)
                skippedlines += 1
            elif not (fields[3].isdigit() and fields[4].isdigit()):
                logger.warning("non-number coordinates in "+filename)
                logger.warning(currline)

                skippedlines += 1
            else:
                if replacename:
                    featname = genename
                yield GenomeRange( orgdb, fields[0],fields[3],fields[4],fields[6], name = featname, fastafile = seqfile, data = {"biotype":biotype, "source":genesource, "genename":genename,"feature":fields[2]})
            
def readbed(filename, orgdb="genome", seqfile= None, includeintrons = False):
    bedfile = None
    if filename == "stdin":
        bedfile = sys.stdin
    elif filename.endswith(".gz"):
        bedfile = gzip.open(filename, 'rt', encoding='utf-8')
    else:
        bedfile = open(filename, "r")
    skippedlines = 0
    
    for currline in bedfile:
        #print currline
        data = dict()
        if currline.startswith('track') or currline.startswith('#'):
            continue
        fields = currline.rstrip().split()
        if len(fields) > 2:
            if len(fields) < 5:
                strand = "+"
            else:
                strand = fields[5]
            if not (strand == "+" or strand == "-"):
                logger.warning("strand error in "+filename)
                skippedlines += 1
            elif not (fields[1].isdigit() and fields[2].isdigit()):
                logger.warning("non-number coordinates in "+filename)
                logger.warning(currline)
                skippedlines += 1
            else:
                if includeintrons and len(fields) > 7:
                    data["blockcount"] = int(fields[9])
                    data["blocksizes"] = tuple(int(curr) for curr in fields[10].rstrip(",").split(","))
                    data["blockstarts"] = tuple(int(curr) for curr in fields[11].rstrip(",").split(",")) 
                yield GenomeRange( orgdb, fields[0],fields[1],fields[2],strand, name = fields[3], fastafile = seqfile, data = data)
    
    if skippedlines > 0:
        logger.warning("skipped "+str(skippedlines)+" in "+filename)



# Column carrying the per-feature counts replicate correlation is computed from. Normalized
# rather than raw, so a depth difference between samples is not read as a biological one.
REPLICATE_CORR_VALUE = 'nreads_total_unique_norm'


def replicate_correlation(obs, group_col='group', value_col=REPLICATE_CORR_VALUE,
                          feature_col='trna', sample_col='sample'):
    '''
    Pearson r-squared between every pair of samples, summarised by whether the pair shares a group.

    The headline number is `separation` -- mean within-group r2 minus mean between-group r2. It
    says whether the experiment's own grouping is visible in its counts at all, and it is the
    objective test of whether a processing step helped: removing a technical artifact should widen
    it, because replicates become more alike while genuinely different conditions do not.

    Counts are log1p-transformed first. tRNA abundance spans orders of magnitude, and on raw
    counts Pearson r would be dominated by whichever handful of transcripts happen to be most
    abundant, making every pair of samples look near-identical regardless of grouping.

    Returns per-sample and per-pair tables alongside the summary, since the summary alone cannot
    identify *which* sample is the problem.
    '''
    import numpy as np
    import pandas as pd

    matrix = obs.pivot_table(index=feature_col, columns=sample_col, values=value_col,
                             aggfunc='sum', observed=True).dropna()
    matrix = np.log1p(matrix)
    r2 = matrix.corr(method='pearson') ** 2

    groups = obs.drop_duplicates(sample_col).set_index(sample_col)[group_col].to_dict()
    samples = list(r2.columns)

    pairs = []
    for i, s1 in enumerate(samples):
        for s2 in samples[i + 1:]:
            pairs.append({'sample_1': s1, 'sample_2': s2,
                          'group_1': groups.get(s1), 'group_2': groups.get(s2),
                          'same_group': groups.get(s1) == groups.get(s2),
                          'r2': float(r2.loc[s1, s2])})
    pairs_df = pd.DataFrame(pairs)

    within = pairs_df.loc[pairs_df['same_group'], 'r2'] if len(pairs_df) else pd.Series(dtype=float)
    between = pairs_df.loc[~pairs_df['same_group'], 'r2'] if len(pairs_df) else pd.Series(dtype=float)

    per_sample = []
    for s in samples:
        own = [p['r2'] for p in pairs
               if (p['sample_1'] == s or p['sample_2'] == s) and p['same_group']]
        other = [p['r2'] for p in pairs
                 if (p['sample_1'] == s or p['sample_2'] == s) and not p['same_group']]
        per_sample.append({
            'sample': s,
            'group': groups.get(s),
            'mean_r2_within_group': float(np.mean(own)) if own else float('nan'),
            'mean_r2_other_groups': float(np.mean(other)) if other else float('nan'),
            'n_replicates': len(own),
        })

    return {
        'summary': {
            'within_mean': float(within.mean()) if len(within) else float('nan'),
            'between_mean': float(between.mean()) if len(between) else float('nan'),
            'separation': (float(within.mean()) - float(between.mean()))
                          if len(within) and len(between) else float('nan'),
            'n_within_pairs': int(len(within)),
            'n_between_pairs': int(len(between)),
        },
        'per_sample': pd.DataFrame(per_sample),
        'pairs': pairs_df,
    }
