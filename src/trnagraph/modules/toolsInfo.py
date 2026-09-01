#!/usr/bin/env python3
'''
Report what an AnnData object actually contains, so a user can see what to type.

Every `graph`, `analyze cluster` and `tools log2fc` option that names a grouping column, a
coverage type, a readtype or a variant is drawn from a vocabulary that lives inside the .h5ad
and nowhere else. Before this command the only way to read that vocabulary was to open the
object in Python, which is why a mistyped label mattered as much as it did -- there was no
cheap way to check one. toolsTG.validate_labels() points here by name when it rejects a
label, so what this reports has to be exactly the vocabulary it validates against.

The structured summary is the product; the human rendering and the --json payload are two
views of it. Keeping it that way means the text a user reads and the JSON a script parses
cannot drift apart.
'''

import json
import logging

import anndata as ad

from . import toolsTG

logger = logging.getLogger(__name__)

#: Unique values printed inline per column before truncating. A cap is required rather than
#: nice to have: 'trna' alone is 436 values on a human build, which would bury every column
#: worth reading under one that is never typed by hand. --column lifts it for one column.
VALUE_CAP = 20

#: Characters shown per value before eliding. The count cap says how MANY values are printed
#: and says nothing about how WIDE each is: obs['refseq'] holds 75-90nt sequences, so twenty
#: of them ran to about 2000 characters on one line. --column lifts this too.
VALUE_WIDTH = 32


class AnnDataInspector:
    '''Summarise one AnnData object's columns, keys and vocabularies.'''

    def __init__(self, args):
        self.args = args
        self.logger = logging.getLogger(__name__)

    def _describe_frame(self, frame):
        '''
        One entry per column of an obs/var frame: what it is called, and what is in it.

        `--column` narrows the frame to that one column and lifts the cap for it. Narrowing
        as well as uncapping is the point: lifting the cap across every column would
        reproduce exactly the unreadable output the cap exists to prevent.
        '''
        wanted = getattr(self.args, 'column', None)
        described = []
        for name in frame.columns:
            if wanted is not None and str(name) != str(wanted):
                continue
            series = frame[name]
            uniques = sorted({str(v) for v in series.dropna().unique()})
            entry = {
                'name': str(name),
                'dtype': str(series.dtype),
                # The count is reported from the FULL set even when the values are capped:
                # it is what tells the reader the printed list is not the whole story.
                'n_unique': len(uniques),
            }
            if wanted is None and self._is_continuous(series, uniques):
                # A measurement is not a label. obs carries ~30 float count columns on a real
                # object, and enumerating them buried 'group' and 'amino' -- the columns
                # anyone actually types -- under hundreds of values like 0.8767311244266289.
                numeric = series.dropna()
                entry['values'] = []
                entry['truncated'] = False
                entry['range'] = [float(numeric.min()), float(numeric.max())] if len(numeric) else []
            else:
                capped = wanted is None and len(uniques) > VALUE_CAP
                shown = uniques[:VALUE_CAP] if capped else uniques
                if wanted is None:
                    shown = [v if len(v) <= VALUE_WIDTH else v[:VALUE_WIDTH - 3] + '...'
                             for v in shown]
                entry['values'] = shown
                entry['truncated'] = capped
            described.append(entry)
        return described

    def _describe_uns(self, adata):
        '''
        uns keys with enough shape to tell a populated slot from an empty one.

        Values are deliberately not rendered: uns holds whole DataFrames (the log2FC cache,
        the non-tRNA counts) and printing them would drown the report. What a reader needs
        here is whether a slot exists and has anything in it.
        '''
        described = []
        for key in sorted(adata.uns):
            value = adata.uns[key]
            entry = {'name': str(key), 'type': type(value).__name__}
            shape = getattr(value, 'shape', None)
            if shape is not None:
                entry['shape'] = tuple(shape)
            elif hasattr(value, '__len__'):
                entry['shape'] = (len(value),)
            described.append(entry)
        return described

    def _variants(self, adata):
        '''
        The --variant strings this object can actually resolve, as ready-to-paste values.

        Built from the layers and split tags that are present rather than from the four
        normalizations parse_variant() accepts: 'vst:full' resolves only if a vst layer was
        built and 'norm:u60' only if that split was added, so listing all four unconditionally
        would advertise variants that abort the moment they are used. 'norm:full' is always
        available -- it is .X itself, which every object has.
        '''
        norms = ['norm'] + [name for name, layer in toolsTG._VARIANT_LAYER_MAP.items()
                            if name != 'norm' and layer in adata.layers]
        tags = ['full'] + sorted(adata.uns.get('size_splits', {}))
        variants = []
        for tag in tags:
            for norm in norms:
                if tag != 'full':
                    layer = f'norm_{tag}' if norm == 'norm' else f'{toolsTG._VARIANT_LAYER_MAP[norm]}_{tag}'
                    if layer not in adata.layers:
                        continue
                variants.append(f'{norm}:{tag}')
        return variants

    def _is_continuous(self, series, uniques):
        '''
        Whether a column is a measurement rather than a set of labels.

        Numeric AND high-cardinality, both parts required: a low-cardinality numeric column is
        a label (an ordered timepoint or dose is exactly what --volgrp gets pointed at), and
        its values have to stay visible. --column overrides this entirely -- a user who names
        a column has asked for its contents whatever they are.
        '''
        import pandas as pd

        return (pd.api.types.is_numeric_dtype(series)
                and not isinstance(series.dtype, pd.CategoricalDtype)
                and len(uniques) > VALUE_CAP)

    def summary(self, adata):
        '''The structured report: the single source both the text and --json views render.'''
        column = getattr(self.args, 'column', None)
        if column is not None:
            # Reporting nothing for a mistyped --column would be the same silently-wrong
            # answer the strict validator exists to remove, in the very command its error
            # message tells the user to run.
            toolsTG.validate_labels(adata, [('column', column, 'column')])
        return {
            'shape': tuple(adata.shape),
            'obs': self._describe_frame(adata.obs),
            'var': self._describe_frame(adata.var),
            'uns': self._describe_uns(adata),
            'layers': sorted(adata.layers),
            'obsm': sorted(adata.obsm),
            'obsp': sorted(adata.obsp),
            'variants': self._variants(adata),
        }

    def _render_columns(self, entries, heading):
        '''One frame's columns as aligned text, widest-first alignment on the name.'''
        if not entries:
            return []
        width = max(len(e['name']) for e in entries)
        lines = [f'{heading} ({len(entries)} column{"s" if len(entries) != 1 else ""})']
        for entry in entries:
            if 'range' in entry:
                low, high = (entry['range'] + ['', ''])[:2]
                values = f'range {low} to {high}' if entry['range'] else 'no values'
            else:
                values = ', '.join(entry['values'])
            if entry['truncated']:
                # Say how many were withheld and how to see them: a bare ellipsis would leave
                # the reader with no way to get the rest.
                hidden = entry['n_unique'] - len(entry['values'])
                values += f' ... +{hidden} more (--column {entry["name"]} for all)'
            lines.append(f'  {entry["name"]:<{width}}  [{entry["dtype"]}, {entry["n_unique"]}]  {values}')
        return lines

    def render(self, summary):
        '''The human view of the summary.'''
        lines = [f'AnnData object: {summary["shape"][0]} observations x {summary["shape"][1]} variables', '']
        lines += self._render_columns(summary['obs'], 'obs')
        if summary['obs'] and summary['var']:
            lines.append('')
        lines += self._render_columns(summary['var'], 'var')
        # A --column report is a focused dump of one column; the object-wide slots below would
        # bury the thing that was actually asked for.
        if getattr(self.args, 'column', None) is not None:
            return '\n'.join(lines)
        for heading, key in (('layers', 'layers'), ('obsm', 'obsm'), ('obsp', 'obsp')):
            if summary[key]:
                lines += ['', f'{heading}: {", ".join(summary[key])}']
        if summary['uns']:
            lines += ['', 'uns']
            for entry in summary['uns']:
                shape = f' {tuple(entry["shape"])}' if 'shape' in entry else ''
                lines.append(f'  {entry["name"]} [{entry["type"]}{shape}]')
        lines += ['', f'--variant: {", ".join(summary["variants"])}']
        return '\n'.join(lines)

    def run(self):
        '''Read the object named by --input and return the requested view of it.'''
        adata = ad.read_h5ad(self.args.anndata)
        summary = self.summary(adata)
        if getattr(self.args, 'json', False):
            return json.dumps(summary, indent=2, default=str)
        return self.render(summary)
