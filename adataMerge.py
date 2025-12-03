import anndata as ad
import numpy as np
import pandas as pd
from . import toolsTG

class anndataMerger():
    '''
    Class for merging multiple AnnData objects into a single object. This is useful for combining multiple tRAX runs into a single object for analysis.
    '''
    def __init__(self, args):
        self.adata1 = ad.read_h5ad(args.anndata1)
        self.adata2 = ad.read_h5ad(args.anndata2)
        self.args = args

    def merge(self):
        if len(np.intersect1d(self.adata1.obs.index, self.adata2.obs.index))>0:
            raise Exception('WARNING: The two AnnData objects have overlapping indices. This will cause issues with merging as duplicate groups/samples may occur. \
                  Regenerate your input annData objects with unique sample names across all merged objects.\n')
        # print diff in columns betwen adata1 and adata2
        if len(np.setdiff1d(self.adata1.obs.columns, self.adata2.obs.columns))>0 or len(np.setdiff1d(self.adata2.obs.columns, self.adata1.obs.columns))>0:
            print('The following columns are not present in each AnnData object and will be dropped:\n'+ \
                  ' '.join(set(np.setdiff1d(self.adata1.obs.columns, self.adata2.obs.columns)) | set(np.setdiff1d(self.adata2.obs.columns, self.adata1.obs.columns))))
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
            print(str(len(set(np.setdiff1d(self.adata1.uns['nontRNA_counts'].index, self.adata2.uns['nontRNA_counts'].index))))+' genes are not present in AnnData object 1 and '+ \
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
                print('The following RNA categories are not present in each AnnData object and will be filled with 0 counts:\n'+ \
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
        print(f'Writing h5ad database object to {self.args.output}')
