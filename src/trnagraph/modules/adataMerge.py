import logging

import anndata as ad
import numpy as np
import pandas as pd
from . import toolsTG

class anndataMerger():
    '''
    Class for merging multiple AnnData objects into a single object. This is useful for combining multiple tRAX runs into a single object for analysis.
    '''
    def __init__(self, args):
        self.logger = logging.getLogger(__name__)
        self.adata1 = ad.read_h5ad(args.anndata1)
        self.adata2 = ad.read_h5ad(args.anndata2)
        self.args = args

    def merge(self):
        # Conflicting-run-info validation: refuse (unless --force) if the two objects' own
        # recorded build provenance (database/gtf) disagrees -- mirrors `analyze addsplit`'s
        # --force-gated check (adataBuild.py.add_split), previously only implemented there and
        # missing entirely from `tools merge`.
        flags1 = self.adata1.uns.get('trnagraphruninfo', {}).get('flags', {})
        flags2 = self.adata2.uns.get('trnagraphruninfo', {}).get('flags', {})
        conflicts = []
        for flag_key, label in [('database', 'database'), ('gtf', 'gtf')]:
            value1 = flags1.get(flag_key)
            value2 = flags2.get(flag_key)
            if value1 is not None and value2 is not None and value1 != value2:
                conflicts.append(f"object 1 was built with {label}='{value1}', but object 2 was built with {label}='{value2}'")
        if conflicts:
            message = "Detected conflicting build provenance between the two AnnData objects:\n  " + "\n  ".join(conflicts)
            if not getattr(self.args, 'force', False):
                raise ValueError(message + "\nPass --force to proceed anyway.")
            self.logger.warning(f"WARNING: {message}\nProceeding anyway due to --force.")

        if len(np.intersect1d(self.adata1.obs.index, self.adata2.obs.index))>0:
            raise Exception('WARNING: The two AnnData objects have overlapping indices. This will cause issues with merging as duplicate groups/samples may occur. \
                  Regenerate your input annData objects with unique sample names across all merged objects.\n')
        # print diff in columns betwen adata1 and adata2
        if len(np.setdiff1d(self.adata1.obs.columns, self.adata2.obs.columns))>0 or len(np.setdiff1d(self.adata2.obs.columns, self.adata1.obs.columns))>0:
            self.logger.info('The following columns are not present in each AnnData object and will be dropped:\n'+ \
                  ' '.join(set(np.setdiff1d(self.adata1.obs.columns, self.adata2.obs.columns)) | set(np.setdiff1d(self.adata2.obs.columns, self.adata1.obs.columns))))
        # ad.concat's uns_merge='same' (below) silently drops any uns key -- including
        # uns['size_splits'] -- that isn't byte-identical across both objects, and drops
        # layers/obsm/obsp entries not present in both. Size-split variants added to only one
        # side (e.g. u60/o60 added via `analyze addsplit` to just one of these two objects)
        # will be silently lost from the merged result -- warn so that isn't a silent surprise.
        splits1 = set(self.adata1.uns.get('size_splits', {}).keys())
        splits2 = set(self.adata2.uns.get('size_splits', {}).keys())
        if splits1 != splits2:
            self.logger.info('WARNING: The two AnnData objects have different size-split variants present '
                  f'(object 1: {sorted(splits1) or "none"}, object 2: {sorted(splits2) or "none"}). '
                  'ad.concat only keeps uns/layers/obsm/obsp entries that are identical/present in both objects, '
                  'so mismatched split-variant data will be dropped from the merged result.\n')
        # Merge uns data
        amino_counts = pd.concat([self.adata1.uns['amino_counts'], self.adata2.uns['amino_counts']], axis=1).fillna(0)
        self.adata1.uns['amino_counts'] = amino_counts
        self.adata2.uns['amino_counts'] = amino_counts
        anticodon_counts = pd.concat([self.adata1.uns['anticodon_counts'], self.adata2.uns['anticodon_counts']], axis=1).fillna(0)
        self.adata1.uns['anticodon_counts'] = anticodon_counts
        self.adata2.uns['anticodon_counts'] = anticodon_counts
        # Correct for non-tRNA counts having different genes
        if self.args.dropno:
            nontRNA_counts = pd.concat([self.adata1.uns['nontRNA_counts'], self.adata2.uns['nontRNA_counts']], axis=1).dropna()
        else:
            self.logger.info(str(len(set(np.setdiff1d(self.adata1.uns['nontRNA_counts'].index, self.adata2.uns['nontRNA_counts'].index))))+' genes are not present in AnnData object 1 and '+ \
                  str(len(set(np.setdiff1d(self.adata2.uns['nontRNA_counts'].index, self.adata1.uns['nontRNA_counts'].index))))+' genes are not present in AnnData object 2. '+\
                  'They will be filled with 0 counts where no overlap occurs. If you wish to remove these genes instead, please use the --dropno option.\n')
            nontRNA_counts = pd.concat([self.adata1.uns['nontRNA_counts'], self.adata2.uns['nontRNA_counts']], axis=1).fillna(0)
        self.adata1.uns['nontRNA_counts'] = nontRNA_counts
        self.adata2.uns['nontRNA_counts'] = nontRNA_counts
        # Correct for missing RNA categories in uns type_counts
        if self.args.droprna:
            type_counts = pd.concat([self.adata1.uns['type_counts'], self.adata2.uns['type_counts']], axis=1).dropna()
        else:
            if len(np.setdiff1d(self.adata1.uns['type_counts'].index, self.adata2.uns['type_counts'].index))>0 or len(np.setdiff1d(self.adata2.uns['type_counts'].index, self.adata1.uns['type_counts'].index))>0:
                self.logger.info('The following RNA categories are not present in each AnnData object and will be filled with 0 counts:\n'+ \
                      ' '.join(set(np.setdiff1d(self.adata1.uns['type_counts'].index, self.adata2.uns['type_counts'].index)) | set(np.setdiff1d(self.adata2.uns['type_counts'].index, self.adata1.uns['type_counts'].index)))+'\n'+ \
                      'If you wish to remove these RNA categories instead, please use the --droprna option.\n')
            type_counts = pd.concat([self.adata1.uns['type_counts'], self.adata2.uns['type_counts']], axis=1).fillna(0)      
        self.adata1.uns['type_counts'] = type_counts
        self.adata2.uns['type_counts'] = type_counts
        # Merge AnnData objects
        self.adata = ad.concat([self.adata1, self.adata2], merge='unique', uns_merge='same')
        # Recompute the 'group' log2FC value/pval to uns since it is the default for the volcano/heatmap and saves time later and needs to exist for some plots
        for i in [20,40,80,100,200]: # These are common read cutoffs for tRNAseq
            toolsTG.adataLog2FC(self.adata, 'group', 'nreads_total_unique_norm', readcount_cutoff=i, config_name='default', overwrite=True).main()
            toolsTG.adataLog2FC(self.adata, 'group', 'nreads_total_norm', readcount_cutoff=i, config_name='default', overwrite=True).main()
        # Save merged AnnData object
        self.adata.write(f'{self.args.output}')
        self.logger.info(f'Writing h5ad database object to {self.args.output}')
