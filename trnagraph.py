#!/usr/bin/env python3

import numpy as np
import pandas as pd
import anndata as ad
import argparse
import os
import sys
import json
import contextlib
import multiprocessing
# Custom functions
import toolsTG
import plotsBar
import plotsCount
import plotsCluster
import plotsCompare
import plotsCorrelation
import plotsCoverage
import plotsHeatmap
import plotsSeqlogo
import plotsPca
import plotsRadar
import plotsVolcano
# Cluster functions
import umap
from sklearn.preprocessing import RobustScaler
# from sklearn.preprocessing import Normalizer
import hdbscan

class trax2anndata():
    '''
    Create h5ad AnnData object from a trax run
    '''
    def __init__(self, traxdir, metadata, output):
        '''
        Initialize trax2anndata object by loading various files from tRAX run
        '''
        # Add unique feature column to coverage file for alignment and sorting
        traxcoverage = traxdir + '/' + traxdir.split('/')[-1] + '-coverage.txt'
        self.coverage = pd.read_csv(traxcoverage, sep='\t', header=0)
        self.coverage['uniquefeat'] = self.coverage['Feature'] + '_' + self.coverage['Sample']
        self.positions = pd.unique(self.coverage['position'])
        # Add size factors to coverage file
        sizefactors = traxdir + '/' + traxdir.split('/')[-1] + '-SizeFactors.txt'
        self.size_factors = pd.read_csv(sizefactors, sep=" ", header=0).to_dict('index')[0]
        self.size_factors_list = None
        # For adding unique counts to coverage file
        trnauniquecounts = traxdir + '/unique/' + traxdir.split('/')[-1] + '-unique-trnas.txt' #'-trnauniquecounts.txt'
        self.unique_counts = pd.read_csv(trnauniquecounts, sep='\t', header=0).to_dict('index')
        # For adding normalized read counts to coverage file split by read type
        normalizedreadcounts = traxdir + '/' + traxdir.split('/')[-1] + '-normalizedreadcounts.txt'
        normalized_read_counts = pd.read_csv(normalizedreadcounts, sep='\t', header=0)
        # Clean all non tRNAs from normalized read counts by removing all rows that do not have a tRNA in the feature column
        non_trna_read_counts = normalized_read_counts[~(normalized_read_counts.index.str.contains('tRNA'))]
        self.non_trna_read_counts = non_trna_read_counts[~(non_trna_read_counts.index.str.contains('tRX'))]
        normalized_read_counts = normalized_read_counts[(normalized_read_counts.index.str.contains('tRNA')) | (normalized_read_counts.index.str.contains('tRX'))]
        self.read_types = pd.unique(normalized_read_counts.index.str.split('_').str[1])
        self.normalized_read_counts = normalized_read_counts.to_dict('index')
        # For adding anticoodon counts to coverage file
        anticodoncounts = traxdir + '/' + traxdir.split('/')[-1] + '-anticodoncounts.txt'
        self.anticodon_counts = pd.read_csv(anticodoncounts, sep='\t')
        # For adding type counts to coverage file
        typecounts = traxdir + '/' + traxdir.split('/')[-1] + '-typecounts.txt'
        self.type_counts = pd.read_csv(typecounts, sep='\t')
        typerealcounts = traxdir + '/' + traxdir.split('/')[-1] + '-typerealcounts.txt'
        self.type_real_counts = pd.read_csv(typerealcounts, sep='\t')
        # For adding amino acid counts to coverage file
        aminoacidcounts = traxdir + '/' + traxdir.split('/')[-1] + '-aminocounts.txt'
        self.amino_counts = pd.read_csv(aminoacidcounts, sep='\t')
        # For adding metadata to adata object
        metadata_type = '\t' if metadata.endswith('.tsv') else ',' if metadata.endswith('.csv') else ' '
        try:
            self.metadata = pd.read_csv(metadata, sep=metadata_type, header=None, index_col=None)
        except:
            raise ValueError(f'Could not read metadata file, check to make sure it is formated correctly: {metadata}')
        # turn first row into a list then remove it from the dataframe
        self.observations = self.metadata.iloc[0].dropna().values.tolist()
        # self.metadata.columns = self.observations
        self.metadata = self.metadata.iloc[1:]
        # Create list of reference sequences from actualbase column of coverage file - skips gap positions
        self.seqs = self._seq_build_()
        self.seqs_full = self._seq_build_(gap=True)
        # Make sure that observations are not going to be generated automatically
        auto_obs = ['trna', 'iso', 'amino', 'deseq2_sizefactor', 'refseq', 'dataset', 'pseudogene']
        if any(x in self.observations for x in auto_obs):
            raise ValueError(f'The following observation categories will automatically be generated please remove these if you included them: {auto_obs}')
        # Make sure that observations are unique
        if len(self.observations) != len(set(self.observations)):
            raise ValueError(f'Observation categories must be unique, please remove duplicates from the observation catgories you wish to generate: {self.observations}')
        # Make sure that sample and group are the first two observations
        if self.observations[0] != 'sample' or self.observations[1] != 'group':
            raise ValueError(f'The first two observation categories must be "sample" and "group" please reorder your observation categories to match the following: ["sample", "group", ...]: {self.observations}')
        # Add manual observations to obs list if they are not provided or if the length of the observations list does not match the number of parameters in the trax coverage file
        if len(self.observations) != len(self.metadata.columns):
            diff_obs_count = len(self.metadata.columns)-len(self.observations)
            print(f'Number of observations does not match number of parameters in trax coverage file by {diff_obs_count}. To create a more specific database object, please provide the correct number of observations.')
            if diff_obs_count > 0:
                print(f'Adding {diff_obs_count} observations to the end of the list')
                self.observations += ['obs_' + str(x) for x in range(diff_obs_count)]
            if diff_obs_count < 0:
                print(f'Removing {abs(diff_obs_count)} observations from the end of the list')
                self.observations = self.observations[:diff_obs_count]
        # Add align observation categories to metadata
        self.metadata.columns = self.observations
        self.metadata.set_index('sample', inplace=True)
        self.metadata = self.metadata.to_dict()
        # For adding the tRAX runinfo
        with open(traxdir + '/' + traxdir.split('/')[-1] + '-runinfo.txt', 'r') as f:
            runinfo = f.readlines()
        # Split the strings and convert to dictionary
        runinfo = [x.rstrip().split('\t') for x in runinfo]
        self.trax_run_info = {x[0]: x[1] for x in runinfo if len(x) == 2}
        if self.trax_run_info['git version'] == 'Cannot find git version':
            print('WARNING: Could not find git version in trax runinfo file. This is likely because the version of trax was not from a git repository or downloaded directly.\n'+ \
                  'Please make sure that the version of trax is at least v1.1.0-beta or later. Merging datasets may not work properly without matching git versions.\n')
        # Search self.trax_run_info['command'] for '--nofrag' and raise a warning if present
        if '--nofrag' in self.trax_run_info['command']:
            raise ValueError('The --nofrag option was used in the trax run, please rerun trax without this flag. This option is not recommended for generating a database object for use with trnagraph.\n' + \
                             'The --nofrag option combines fragment types into whole-reads and will not generate the necessary observations for a full database object.\n')
        # Add trnagraph version to trax run info bashed on github hash - Will be changed to git describe once the package is deployed
        trnagraphdir = os.path.dirname(os.path.abspath(__file__))
        self.trnagraph_run_info = {'expname': traxdir.split('/')[-1], 
                                   'time': os.popen('date').read().rstrip(),
                                   'trax_directory': traxdir,
                                   'git version': os.popen('git --git-dir='+trnagraphdir+'/.git describe').read().rstrip(),
                                   'git version hash': os.popen('git --git-dir='+trnagraphdir+'/.git rev-parse HEAD').read().rstrip()}
        # Output file name
        self.output = output
        # Names of coverage types to add to adata object from coverage file
        self.cov_types = ['coverage', 'readstarts', 'readends', 'uniquecoverage', 'multitrnacoverage',
                          'multianticodoncoverage', 'multiaminocoverage','tRNAreadstotal', 'mismatchedbases',
                          'deletedbases', 'adenines', 'thymines', 'cytosines', 'guanines', 'deletions']

    def create(self):
        '''
        Create h5ad database object
        '''
        # Build obs and x dataframes
        x_dfs = []
        for cov_type in self.cov_types:
            x_df = self._x_build_(cov_type)
            # Build size factors list if it does not exist
            if not self.size_factors_list:
                self.size_factors_list = [self.size_factors.get(i) for i in ['_'.join(x.split('_')[1:]) for x in x_df.index.values]]
            # 'adenines', 'thymines', 'cytosines', 'guanines', 'deletions' are already raw counts so they need to be normalized by size factor first
            if cov_type in ['adenines', 'thymines', 'cytosines', 'guanines', 'deletions']:
                x_df = x_df.div(self.size_factors_list, axis=0)
            x_dfs.append(x_df)
        # Build obs dataframe
        obs_df = self._obs_build_(x_df)
        x_df = pd.concat(x_dfs, axis=1, sort=False)
        x_df = x_df.astype('float64')  # Not sure if this is the dtype I want to use defaults to float64
        # Rename columns of x_df to include position and coverage type
        clist = [[p + '_' + cov for p in self.positions] for cov in self.cov_types]
        clist = [a for l in clist for a in l]
        x_df.columns = clist
        # obs_df,x_df = self._group_sort_(obs_df,x_df) # Not sure if I need this function
        # Check that the index of the obs and x dataframes are the same
        if not obs_df.index.equals(x_df.index):
            raise ValueError('The index of the obs and x dataframes are not the same. This means somthing went wrong in the sorting process.')
        # Build adata object
        adata = self._adata_build_(obs_df, x_df)
        # Add size factors to adata object as raw layer
        adata.layers['raw'] = adata.X * adata.obs['deseq2_sizefactor'].values[:,None]
        # Quality check adata by dropping NaN values and printing summary
        if adata.obs.isna().any(axis=0).any():
            print('WARNING: NaN values found in obs dataframe this is commonly caused by missing samples in your metadata file or havinge a different number of observations per sample.\n' + \
                  'Another consideration is to make sure your .tsv/.csv is using tabs or commas appropriately.\n' + \
                  'It can also be caused by metadata containing NaN or None values. Please check your metadata file to make sure the following are correct\n' + \
                  f'Observations:\n{str(adata.obs.columns[adata.obs.isna().any(axis=0)].tolist())}\n' + \
                  f'Samples with NaN:\n{str(list(set(adata.obs["sample"][adata.obs.isna().any(axis=1)].tolist())))}\n')
        # Add output name to adata object index
        adata.obs.index = [os.path.basename(self.output).split('.')[0] + '_' + str(x) for x in adata.obs.index]
        # Save adata object
        adata.write(self.output)
        print(f'Writing h5ad database object to {self.output}')

    def _seq_build_(self, gap=False):
        # Build reference sequence dataframe
        seq_df = self._x_build_('actualbase')
        if not gap:
            # Drop gap positions from reference sequence dataframe
            seq_df = seq_df.loc[:,~seq_df.columns.str.contains('a')]
            seq_df = seq_df.loc[:,~seq_df.columns.str.contains('b')]
            seq_df = seq_df.loc[:,~seq_df.columns.str.contains('e')]
            seq_df = seq_df.loc[:,~seq_df.columns.str.contains('-1')]
        seqs = [''.join(x) for x in seq_df.values.tolist()]

        return seqs

    def _x_build_(self, cov_type):
        '''
        Build x dataframe from trax coverage file
        '''
        x_df = self.coverage.pivot(index='uniquefeat', values=[cov_type], columns='position')
        cols = x_df.columns.get_level_values(1).values
        x_df = x_df.T.reset_index(drop=True).T
        x_df.columns = cols
        x_df = x_df.reindex(columns=self.positions)
        
        return x_df
    
    def _obs_build_(self, x_df):
        '''
        Build obs dataframe from trax coverage file derived x dataframe
        '''
        # Create obs dataframe from trax coverage file derived x dataframe
        obs_df = pd.DataFrame([[x.split('_')[0], '_'.join(x.split('_')[1:])] for x in x_df.index.values], columns=['trna','sample'], index=x_df.index)
        # Add metadata to obs dataframe
        for i in self.metadata:
            obs_df[i] = [self.metadata.get(i).get(x) for x in obs_df['sample'].values]
        # Add tRNA type, amino acid, iso, refseq, sizefactor, and nreads to obs dataframe
        trna_obs = obs_df['trna'].str.split('-',n=4,expand=True)
        obs_df['pseudogene'] = trna_obs[0]
        obs_df['amino'] = trna_obs[1]
        obs_df['iso'] = trna_obs[2]
        obs_df['refseq'] = self.seqs
        obs_df['refseq_full'] = self.seqs_full
        obs_df['sizefactor'] = self.size_factors_list # [self.size_factors.get(i) for i in ['_'.join(x.split('_')[1:]) for x in x_df.index.values]]
        # obs_df['nreads_unique_raw'] = [self.unique_counts.get(i[0]).get('_'.join(i[1:])) if self.unique_counts.get(i[0]) else 0 for i in [x.split('_') for x in x_df.index.values]] # Some samples may not have any reads for a given tRNA in the unique_counts dictionary might want to double check this
        # Create unique counts for each tRNA split by type
        for rt in ['fiveprime','threeprime','wholecounts','other']: # Changed whole to wholecounts in recent trax versions
            obs_df[f'nreads_{rt}_unique_raw'] = [self.unique_counts.get(i[0]+f'_{rt}').get('_'.join(i[1:])) if self.unique_counts.get(i[0]+f'_{rt}') else 0 for i in [x.split('_') for x in x_df.index.values]]
            obs_df[f'nreads_{rt}_unique_norm'] = obs_df[f'nreads_{rt}_unique_raw']/obs_df['sizefactor']
        # Add total unique read counts to obs dataframe
        obs_df['nreads_total_unique_raw'] = obs_df[[f'nreads_{rt}_unique_raw' for rt in ['fiveprime','threeprime','wholecounts','other']]].sum(axis=1)
        obs_df['nreads_total_unique_norm'] = obs_df['nreads_total_unique_raw']/obs_df['sizefactor']
        # Add read counts (from bowtie2) for each tRNA split by type
        for rt in self.read_types:
            obs_df['nreads_' + rt + '_norm'] = [self.normalized_read_counts.get(j[0] + '_' + rt).get('_'.join(j[1:])) if self.normalized_read_counts.get(j[0] + '_' + rt) else 0 for j in [x.split('_') for x in x_df.index.values]]
            obs_df['nreads_' + rt + '_raw'] = obs_df['nreads_' + rt + '_norm']*obs_df['sizefactor']
        # Add total read counts (from bowtie2) to obs dataframe - excluding partial precounts, trailer counts, and antisense counts so that it matches unique counts
        obs_df['nreads_total_raw'] = obs_df[[f'nreads_{rt}_raw' for rt in ['fiveprime','threeprime','wholecounts','other']]].sum(axis=1)
        obs_df['nreads_total_norm'] = obs_df['nreads_total_raw']/obs_df['sizefactor']
        # Add unique feature column to obs dataframe - this is the index of the x dataframe - not sure if this is necessary
        obs_df['uniquefeat'] = obs_df.index.values
        return obs_df

    def _adata_build_(self, obs_df, x_df):
        '''
        Build AnnData object from obs and x dataframes
        '''
        positions = [i.split('_')[0] for i in x_df.columns.values]
        coverage = [i.split('_')[-1] for i in x_df.columns.values]
        # Generate gap list
        gap_list = pd.Series(positions)
        gap_list = gap_list[gap_list.str.contains('a') | gap_list.str.contains('b') | gap_list.str.contains('e') | gap_list.str.contains('-1')]
        gap_list = pd.Series(positions).isin(gap_list).tolist()
        # Build AnnData object
        adata = ad.AnnData(x_df)
        adata.var['gap'] = gap_list
        adata.var['positions'] = positions
        adata.var['coverage'] = coverage
        # Create sprinzl position and fragment type information
        loc_dict, loc_half_dict = {}, {}
        # Define the location of acceptor stem
        loc_dict.update({i:'fiveprime_acceptorstem' for i in [str(i) for i in range(-1,8)]})
        loc_dict.update({i:'threeprime_acceptorstem' for i in [str(i) for i in range(66,77)]})
        loc_half_dict.update({i:'fiveprime' for i in [str(i) for i in range(-1,8)]})
        loc_half_dict.update({i:'threeprime' for i in [str(i) for i in range(66,77)]})
        # Define the location of acceptor stem to d stem
        loc_a_to_d_internal = ['8','9']
        loc_dict.update({i:'a_to_d_internal' for i in loc_a_to_d_internal})
        loc_half_dict.update({i:'fiveprime' for i in loc_a_to_d_internal})
        # Define the location of d stem and loop
        loc_dstem = [str(i) for i in range(10,14)] + [str(i) for i in range(22,26)]
        loc_dict.update({i:'dstem' for i in loc_dstem})
        loc_half_dict.update({i:'fiveprime' for i in loc_dstem})
        loc_dloop = [str(i) for i in range(14,22)] + ['17a','20a','20b']
        loc_dict.update({i:'dloop' for i in loc_dloop})
        loc_half_dict.update({i:'fiveprime' for i in loc_dloop})
        # Define the location of d stem to anticodon stem
        loc_d_to_anticodon_internal = ['26']
        loc_dict.update({i:'d_to_anticodon_internal' for i in loc_d_to_anticodon_internal})
        loc_half_dict.update({i:'fiveprime' for i in loc_d_to_anticodon_internal})
        # Define the location of anticodon stem and loop
        loc_dict.update({i:'fiveprime_anticodonstem' for i in range(27,32)})
        loc_dict.update({i:'threeprime_anticodonstem' for i in range(39,44)})
        loc_half_dict.update({i:'center' for i in range(27,32)})
        loc_half_dict.update({i:'center' for i in range(39,44)})
        loc_anticodonloop = [str(i) for i in range(32,39)]
        loc_dict.update({i:'anticodonloop' for i in loc_anticodonloop})
        loc_half_dict.update({i:'center' for i in loc_anticodonloop})
        # Define the location of anticodon stem to t stem
        loc_anticodon_to_t_internal = [str(i) for i in range(44,49)]
        loc_dict.update({i:'anticodon_to_t_internal' for i in loc_anticodon_to_t_internal})
        loc_half_dict.update({i:'threeprime' for i in loc_anticodon_to_t_internal})
        # Define the location of extension loop
        loc_e = ['e' + str(i) for i in range(1,20)]
        loc_dict.update({i:'extensionloop' for i in loc_e})
        loc_half_dict.update({i:'threeprime' for i in loc_e})
        # Define the location of t stem and loop
        loc_tstem = [str(i) for i in range(49,54)] + [str(i) for i in range(61,66)]
        loc_dict.update({i:'tstem' for i in loc_tstem})
        loc_half_dict.update({i:'threeprime' for i in loc_tstem})
        loc_tloop = [str(i) for i in range(54,61)]
        loc_dict.update({i:'tloop' for i in loc_tloop})
        loc_half_dict.update({i:'threeprime' for i in loc_tloop})
        # Add the location data to adata object
        adata.var['location'] = adata.var['positions'].map(loc_dict)
        adata.var['half'] = adata.var['positions'].map(loc_half_dict)
        # Add metadata dataframe
        adata.obs['dataset'] = os.path.basename(self.output) # Name of the dataset output file - Usefull for combining multiple datasets if the merge function is used
        adata.obs['trna'] = obs_df['trna'].values
        adata.obs['iso'] = obs_df['iso'].values
        adata.obs['amino'] = obs_df['amino'].values
        # Add sample and group metadata
        adata.obs['sample'] = obs_df['sample'].values
        adata.obs['group'] = obs_df['group'].values
        adata.obs['pseudogene'] = obs_df['pseudogene'].values
        # Add custom dataframe obs
        for ob in self.observations:
            adata.obs[ob] = obs_df[ob].values
        # Add the numer of reads per tRNA as observations-annotation to adata from trna unique counts file
        for rt in ['wholecounts','fiveprime','threeprime','other']:
            adata.obs['nreads_' + rt + '_unique_raw'] = obs_df['nreads_' + rt + '_unique_raw'].values
            adata.obs['nreads_' + rt + '_unique_norm'] = obs_df['nreads_' + rt + '_unique_norm'].values
        adata.obs['nreads_total_unique_raw'] = obs_df['nreads_total_unique_raw'].values
        adata.obs['nreads_total_unique_norm'] = obs_df['nreads_total_unique_norm'].values
        # Add the numer of reads per tRNA as observations-annotation to adata from trna normalized read counts file
        for rt in self.read_types:
            adata.obs['nreads_' + rt + '_raw'] = obs_df['nreads_' + rt + '_raw'].values
            adata.obs['nreads_' + rt + '_norm'] = obs_df['nreads_' + rt + '_norm'].values
        adata.obs['nreads_total_raw'] = obs_df['nreads_total_raw'].values
        adata.obs['nreads_total_norm'] = obs_df['nreads_total_norm'].values
        # Add fragment type classification to the adata object
        fragtype = []
        for i in range(adata.shape[0]):
            fiveprimereadends = adata.obs['nreads_fiveprime_unique_norm'].iloc[i]
            fiveprimemean = adata.X[i][adata.var['half']=='fiveprime'].mean()
            fiveprimestd = np.std(adata.X[i][adata.var['half']=='fiveprime'])
            threeprimereadends = adata.obs['nreads_threeprime_unique_norm'].iloc[i]
            threeprimemean = adata.X[i][adata.var['half']=='threeprime'].mean()
            threeprimestd = np.std(adata.X[i][adata.var['half']=='threeprime'])
            totalreads = adata.obs['nreads_total_unique_norm'].iloc[i]
            totalstd = np.std(adata.X[i])
            centermean = adata.X[i][adata.var['half']=='center'].mean()
            centerstd = np.std(adata.X[i][adata.var['half']=='center'])
            # If the readend counts are more than 2 standard deviations away from the opposite end then it is a fragment
            if totalreads <= 20:
                fragtype.append('low_coverage')
            elif fiveprimereadends > np.abs(threeprimereadends + 2*totalstd):
                if np.abs(centermean - adata.X[i][adata.var['location']=='fiveprime_acceptorstem'].mean()) < fiveprimestd:
                    fragtype.append('fiveprime_half')
                else:
                    fragtype.append('fiveprime_fragment')
            elif threeprimereadends > np.abs(fiveprimereadends + 2*totalstd):
                if np.abs(centermean - adata.X[i][adata.var['location']=='threeprime_acceptorstem'].mean()) < threeprimestd:
                    fragtype.append('threeprime_half')
                else:
                    fragtype.append('threeprime_fragment')
            else:
                if np.abs(fiveprimemean - threeprimemean) > totalstd:
                    fragtype.append('other_fragment')
                elif np.abs(adata.X[i][adata.var['location']=='fiveprime_acceptorstem'].mean() - \
                            adata.X[i][adata.var['location']=='threeprime_acceptorstem'].mean()) > centermean + centerstd:
                    fragtype.append('multiple_fragment')
                else:
                    fragtype.append('whole')
        adata.obs['fragment'] = fragtype
        # Add size factor
        adata.obs['deseq2_sizefactor'] = obs_df['sizefactor'].values
        # Add aligned reference sequence based on values in coverage file
        adata.obs['refseq'] = obs_df['refseq'].values
        adata.obs['refseq_full'] = obs_df['refseq_full'].values
        # Add anticodon counts as uns
        adata.uns['anticodon_counts'] = self.anticodon_counts
        # Add amino acid counts as uns
        adata.uns['amino_counts'] = self.amino_counts
        # Add type counts as uns
        adata.uns['type_counts'] = self.type_counts
        adata.uns['type_real_counts'] = self.type_real_counts
        # Add non tRNA counts as uns
        adata.uns['nontRNA_counts'] = self.non_trna_read_counts
        # Add runinfo as uns
        adata.uns['traxruninfo'] = self.trax_run_info
        adata.uns['trnagraphruninfo'] = self.trnagraph_run_info
        # Add 'group' log2FC value/pval to uns since it is the default for the volcano/heatmap and saves time later
        for i in [20,40,80,100,200]: # These are common read cutoffs for tRNAseq
            toolsTG.adataLog2FC(adata, 'group', 'nreads_total_unique_norm', readcount_cutoff=i, config_name='default', overwrite=True).main()
            toolsTG.adataLog2FC(adata, 'group', 'nreads_total_norm', readcount_cutoff=i, config_name='default', overwrite=True).main()
        
        return adata

class anndataGrapher:
    '''
    Class to generate graphs from an AnnData object by calling the appropriate graphing functions
    '''
    def __init__(self, args):
        self.args = args
        self.adata = ad.read_h5ad(self.args.anndata)
        self.config_name = 'default'
        # Load cmap dict for each graph type
        self.cmap_dict = {'bar':self.args.bargrp, 'cluster':self.args.clustergrp, 'compare':self.args.comparegrp1, \
                          'coverage':self.args.covgrp, 'pca':self.args.pcacolors, 'radar':self.args.radargrp}
        # Load all graph types if specified
        if self.args.graphtypes == 'all' or 'all' in self.args.graphtypes:
            self.args.graphtypes = ['bar', 'cluster', 'correlation', 'count', 'coverage', 'heatmap', 'logo', 'pca', 'radar', 'volcano']
            self.args.clusteroverview = True
        # Load max threads available unless specified
        if self.args.threads == 0:
            try:
                # This is a linux only function but is less likely to cause problems than multiprocessing.cpu_count()
                self.args.threads = len(os.sched_getaffinity(0))
            except:
                self.args.threads = multiprocessing.cpu_count()
        # Load config file if specified
        if self.args.config:
            print('Loading config file: ' + self.args.config)
            with open(self.args.config, 'r') as f:
                self.args.config = json.load(f)
            if 'name' in self.args.config:
                self.args.output += '/' + self.args.config['name']
                self.config_name = self.args.config['name']
                print(toolsTG.builder(self.args.output))
            else:
                raise ValueError('Config file must contain a "name" field')
            if 'obs' in self.args.config or 'obs_r' in self.args.config:
                # Dictionary of uns columns and values to filter by as groups and samples since the coulmns are different from the main adata obs
                obs_dict = {i:True for i in self.adata.uns['amino_counts'].columns.values}
                obs_dict.update({i:True for i in self.adata.uns['type_counts'].columns.values})
                filter_dict = self.args.config['obs']
                if 'obs_r' in self.args.config:
                    # Add the inverse of the obs_r filter to the filter_dict
                    for k,v in self.args.config['obs_r'].items():
                        filter_dict[k] = [i for i in self.adata.obs[k].unique() if i not in v]
                for k,v in filter_dict.items():
                    print('Filtering AnnData object by observation: ' + k + ' , ' + str(v))
                    # Filter all uns columns by the observation and update the obs_dict
                    sub_obs_dict = dict(zip(self.adata.obs['sample'], self.adata.obs[k]))
                    sub_obs_dict.update(dict(zip(self.adata.obs['group'], self.adata.obs[k])))
                    obs_dict = {i:False if sub_obs_dict.get(i,False) not in v else True for i,j in obs_dict.items()}
                    # Filter the adata object by the observation
                    self.adata = self.adata[self.adata.obs[k].isin(v), :]
                # Filter uns columns by the obs_dict
                uns_dict = self.adata.uns.copy() # This was the only way I could get the uns columns to update without an implicit copy warning
                for uns_key in ['amino_counts', 'anticodon_counts', 'nontRNA_counts', 'type_counts']:
                    uns_value = self.adata.uns[uns_key].loc[:, [i for i in self.adata.uns[uns_key].columns.values if obs_dict[i]]].copy()
                    uns_dict[uns_key] = uns_value
                self.adata.uns = uns_dict
            if 'var' in self.args.config or 'var_r' in self.args.config:
                filter_dict = self.args.config['var']
                if 'var_r' in self.args.config:
                    # Add the inverse of the var_r filter to the filter_dict
                    for k,v in self.args.config['var_r'].items():
                        filter_dict[k] = [i for i in self.adata.var[k].unique() if i not in v]
                for k,v in filter_dict.items():
                    print('Filtering AnnData object by variable: ' + k + ' , ' + str(v))
                    self.adata = self.adata[:, self.adata.var[k].isin(v)]
            print('Config file loaded.\n')

        else:
            self.args.config = {}
        # Load the colormap if specified
        if self.args.colormap:
            print('Loading colormap file: ' + self.args.colormap)
            with open(self.args.colormap, 'r') as f:
                self.args.colormap = json.load(f)
            print('Colormap loaded.\n')
        # Check for heatmap or volcano in graph types and if present check for readcount_cutoff in log2FC - Will precompute this and save it back to uns
        # This is done now to prevent saving issues later if multiprocessing is used
        log2FC_dict = self.adata.uns['log2FC'].copy()
        if 'heatmap' in self.args.graphtypes or 'volcano' in self.args.graphtypes:
            for grp in list(set([self.args.heatgrp, self.args.volgrp])):
                for readtype in [f'nreads_{i}_norm' for i in self.args.diffrts]: #list(set(self.args.heatrts+self.args.volrts))]:
                    for cutoff in list(set([self.args.heatcutoff, self.args.volcutoff])):
                        toolsTG.adataLog2FC(self.adata, grp, readtype, readcount_cutoff=cutoff, config_name=self.config_name, overwrite=self.args.regen_uns).main()
        if log2FC_dict != self.adata.uns['log2FC'] or self.args.regen_uns:
            print('The log2FC uns dictionary has been updated.\n')
            self.adata.write(self.args.anndata)

    def main(self):
        # Generate graphs
        if self.args.verbose:
            print('Generating graphs with the following parameters:\n')
            for i in self.args.__dict__: print(f'{i}: {self.args.__dict__[i]}')
            print('')
        # Remove bar and coverage from self.args.graphtypes and add them to non_pooled_graphs if they are present
        # This is because the bar and coverage plots already implement multiprocessing
        non_pooled_graphs = []
        if 'bar' in self.args.graphtypes:
            self.args.graphtypes.remove('bar')
            non_pooled_graphs.append('bar')
        if 'coverage' in self.args.graphtypes:
            self.args.graphtypes.remove('coverage')
            non_pooled_graphs.append('coverage')
        # Pool if applicable
        if self.args.threads > 1 and len(self.args.graphtypes) > 1:
            print(f'Multithreading enabled with {self.args.threads} threads to generate plots\n')
            # Create a multiprocessing pool
            pool = multiprocessing.Pool(self.args.threads)
            # Generate graphs
            pool_output = pool.starmap(self.plot, [(gt, True) for gt in self.args.graphtypes])
            # Close the pool and wait for the work to finish
            pool.close()
            pool.join()
            for gt in non_pooled_graphs:
                self.plot(gt)
            for po in pool_output:
                if po:
                    print(po + '\n')
        else:
            # Combine non_pooled_graphs with self.args.graphtypes
            self.args.graphtypes += non_pooled_graphs
            for gt in self.args.graphtypes:
                self.plot(gt)

    def plot(self, gt, threaded=None):
        if threaded:
            threaded = f'Generating {gt} plots...\n'
        else:
            print(f'Generating {gt} plots...')
        adata_c = self.adata.copy()
        # Define the colormap to use for the graph type
        colormap, colormap_tc, colormap_bg = None, None, None # For counts plot
        cmapgrp = gt
        if gt == 'coverage':
            if not self.args.covgrp:
                cmapgrp = 'group'
        cmappar = self.cmap_dict.get(cmapgrp, None)
        if self.args.colormap:
            if cmappar in self.args.colormap:
                colormap = self.args.colormap[cmappar]
            if gt == 'count':
                if 'sample' in self.args.colormap:
                    colormap_bg = self.args.colormap['sample']
                if 'group' in self.args.colormap:
                    colormap_tc = self.args.colormap['group']
        # Create the output directory
        output = self.args.output + '/' + gt + '/'
        if threaded:
            threaded += toolsTG.builder(output) + '\n'
        else:
            print(toolsTG.builder(output))
        # Plot specific parameters
        if gt == 'bar':
            if self.args.barsubgrp:
                plotsBar.visualizer(adata_c, self.args.threads, self.args.barcol, self.args.bargrp, self.args.barsubgrp, self.args.barsort, self.args.barlabel, colormap, output).generate_subplots()
            else:
                plotsBar.visualizer(adata_c, self.args.threads, self.args.barcol, self.args.bargrp, self.args.barsubgrp, self.args.barsort, self.args.barlabel, colormap, output).generate_plots()
        if gt == 'cluster':
            if not 'cluster_runinfo' in self.adata.uns:
                if threaded:
                    threaded += 'No cluster run information found in AnnData object. Please run the cluster command first.\n'
                else:
                    print('No cluster run information found in AnnData object. Please run the cluster command first.\n')
            else:
                threaded = plotsCluster.visualizer(adata_c, self.args.clustergrp, self.args.clusteroverview, self.args.clusternumeric, self.args.clusterlabels, self.args.clustermask, colormap, output, threaded=threaded).generate_plots()
        if gt == 'compare':
            threaded = plotsCompare.visualizer(adata_c, self.args.comparegrp1, self.args.comparegrp2, colormap, output, threaded=threaded)
        if gt == 'correlation':
            threaded = plotsCorrelation.visualizer(adata_c, self.args.corrmethod, self.args.corrgroup, output, threaded=threaded)
        if gt == 'count':
            threaded = plotsCount.visualizer(adata_c, colormap_tc, colormap_bg, output, threaded=threaded) # Need to add better threading to this
        if gt == 'coverage':
            pcV = plotsCoverage.visualizer(adata_c, self.args.threads, self.args.covgrp, self.args.covobs, self.args.covtype, self.args.covgap, self.args.covmethod, colormap, output)
            # Generate folders/subfolders if coveragecombine is specified
            print(toolsTG.builder(f'{output}{self.args.covobs}/'))
            print(toolsTG.builder(f'{output}{self.args.covobs}/low_coverage/'))
            # Generate coverage plots with combine or split pdfs
            if not self.args.combinedpdfonly:
                print('Generating individual coverage plots pdfs...')
                pcV.generate_split()
            print('Generating combined coverage plots pdf...')
            pcV.generate_combine()
        if gt == 'heatmap':
            threaded = plotsHeatmap.visualizer(adata_c, self.args.heatgrp, self.args.diffrts, self.args.heatcutoff, self.args.heatbound, self.args.heatsubplots, output, threaded=threaded, config_name=self.config_name, overwrite=self.args.regen_uns)
        if gt == 'logo':
            plotsSeqlogo.visualizer(adata_c, self.args.logogrp, self.args.logomanualgrp, self.args.logomanualname, self.args.logopseudocount, self.args.logosize, self.args.ccatail, self.args.pseudogenes, self.args.logornamode, output).generate_plots()
        if gt == 'pca':
            threaded = plotsPca.visualizer(adata_c, self.args.pcamarkers, self.args.pcacolors, self.args.pcareadtypes, colormap, output, threaded=threaded)
        if gt == 'radar':
            if 'all' in self.args.radarmethod:
                self.args.radarmethod = ['mean', 'median', 'max', 'sum']
            for radarmethod in self.args.radarmethod:
                pRd = plotsRadar.visualizer(adata_c, self.args.radargrp, radarmethod, self.args.radarscaled, colormap, output, threaded=threaded)
                threaded = pRd.isotype_plots()
        if gt == 'volcano':
            threaded = plotsVolcano.visualizer(adata_c, self.args.volgrp, self.args.diffrts, self.args.volcutoff, output, threaded=threaded, config_name=self.config_name, overwrite=self.args.regen_uns)
        # Return threaded output  
        if threaded:
            threaded += f'{gt.capitalize()} plots generated!\n'
            return threaded
        else:
            print(f'{gt.capitalize()} plots generated!\n')

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
        # Save merged AnnData object
        self.adata.write(f'{self.args.output}')
        print(f'Writing h5ad database object to {self.args.output}')

class anndataCluster():
    '''
    Class for performing UMAP clustering on an AnnData object
    '''
    def __init__(self, args):
        self.adata = ad.read_h5ad(args.anndata)
        self.overwrite = args.overwrite
        self.output = args.output
        self.randomstate = args.randomstate
        self.readcutoff = args.readcutoff
        self.coveragetype = args.coveragetype
        self.group_n_components = args.ncomponentgrp
        self.group_neighbors_cluster = args.neighborclusgrp
        self.group_neighbors_plot = args.neighborstdgrp
        self.group_hdbscan_min_samples = args.hdbscanminsampgrp
        self.group_hdbscan_min_cluster_size = args.hdbscanminclugrp
        self.sample_n_components = args.ncomponentsmp
        self.sample_neighbors_cluster = args.neighborclusmp
        self.sample_neighbors_plot = args.neighborstdsmp
        self.sample_hdbscan_min_samples = args.hdbscanminsampsmp
        self.sample_hdbscan_min_cluster_size = args.hdbscanminclusmp
        self.cluster_obs = args.clusterobsexperimental

    def main(self):
        # Check if the output file already exists
        if os.path.isfile(self.output) and self.overwrite == False:
            try:
                if 'cluster_runinfo' in ad.read_h5ad(self.output).uns:
                    print('Cluster information already present in AnnData object. No new clustering will be performed. If you wish to overwrite the existing clustering information, please use the --overwrite option.')
                    exit()
            except:
                print('Output file already exists but not in AnnData format. Please remove the file or use the --overwrite option.')
                exit()
        # Subset the AnnData object to only include samples with a minimum number of reads
        print('Performing UMAP clustering...')
        # Preprocess AnnData object
        print('Preprocessing AnnData object...')
        adata_sub_sample = self.adataPreprocess(self.adata.copy())
        adata_sub_group = self.adataPreprocess(self.adata.copy(), grpby='group')
        # # Cluster the data
        print('Clustering AnnData object...')
        sample_df = self.adataCluster(adata_sub_sample, self.sample_neighbors_plot, self.sample_neighbors_cluster, self.sample_hdbscan_min_samples, self.sample_hdbscan_min_cluster_size, self.sample_n_components)
        group_df = self.adataCluster(adata_sub_group, self.group_neighbors_plot, self.group_neighbors_cluster, self.group_hdbscan_min_samples, self.group_hdbscan_min_cluster_size, self.group_n_components)
        # Add the cluster information to the original AnnData object
        print('Adding cluster information to original AnnData object...')
        self.adata = self.adataCombine(self.adata, sample_df, 'sample')
        self.adata = self.adataCombine(self.adata, group_df, 'group')
        # Save all the variables as a dictionary in adata.uns
        if not self.randomstate:
            self.randomstate = -1
        cluster_runinfo = {'sample_neighbors_cluster':self.sample_neighbors_cluster,'sample_neighbors_plot':self.sample_neighbors_plot,'sample_hdbscan_min_samples':self.sample_hdbscan_min_samples,\
                           'sample_hdbscan_min_cluster_size':self.sample_hdbscan_min_cluster_size,'sample_n_components':self.sample_n_components,\
                           'group_neighbors_cluster':self.group_neighbors_cluster,'group_neighbors_plot':self.group_neighbors_plot,'group_hdbscan_min_samples':self.group_hdbscan_min_samples,\
                           'group_hdbscan_min_cluster_size':self.group_hdbscan_min_cluster_size,'group_n_components':self.group_n_components,\
                           'readcutoff':self.readcutoff,'randomstate':self.randomstate}
        # Convert the dictionary to a dataframe then transform it to a single column dataframe
        self.adata.uns['cluster_runinfo'] = cluster_runinfo
        # Save the AnnData object
        print(f'Writing h5ad database object to: {self.output}')
        self.adata.write(f'{self.output}')

    def adataPreprocess(self, adata, grpby=None):
        '''
        Preprocess AnnData object for clustering
        '''
        # Clean up the data by removing gaps from var and subsetting the var to uniquecoverage, readstarts, readends, mismatchedbases, and deletions
        adata = adata[:, ~adata.var['gap']]
        adata = adata[:, adata.var['coverage'].isin(self.coveragetype)]
        # Take the columns from adata.obs[self.cluster_obs] and add them to the adata.X and adata.var
        if self.cluster_obs:
            # Create a dataframe from the cluster_obs column
            df_X = pd.DataFrame(adata.obs[self.cluster_obs])
            # replace NaN with 0
            df_X = df_X.fillna(0)
            df_X = np.concatenate((adata.X, df_X.values), axis=1)
            # Add the cluster_obs column to the adata.var by copying a row from the adata.var then filling with NaN
            df_var = pd.DataFrame(adata.var.iloc[0:len(self.cluster_obs),:])
            # Replace various var data with arbitrary values
            df_var['positions'] = [1000]*len(self.cluster_obs)
            df_var['coverage'] = ['custom']*len(self.cluster_obs)
            df_var['location'] = ['custom']*len(self.cluster_obs)
            df_var['half'] = ['custom']*len(self.cluster_obs)
            df_var.index = [self.cluster_obs]
            df_var = pd.concat([adata.var, df_var])
            # Convert the dataframes to an AnnData object
            adata = ad.AnnData(X=df_X, obs=adata.obs, var=df_var)
        # Filter out Und samples from that amino acid category
        adata = adata[~(adata.obs['amino'] == 'Und'), :]
        if grpby:
            # Combine the mean of the samples in each group for the entire dataset
            df = pd.DataFrame(adata.X, index=adata.obs.index, columns=adata.var.index)
            df['trna'] = adata.obs['trna']
            df[grpby] = adata.obs[grpby]
            df = df.groupby(['trna', grpby], observed=True).mean()
            df_obs = adata.obs.groupby(['trna', grpby], observed=True).first()
            # Remove multiindex
            df.index = df.index.map('_'.join)
            df_obs.index = df_obs.index.map('_'.join)
            df_obs.index.all() == df.index.all()
            # Convert the dataframe to an AnnData object
            adata = ad.AnnData(X=df, obs=df_obs, var=adata.var)
         # Filter out samples with low coverage
        adata = adata[adata.obs['nreads_total_unique_raw'] >= self.readcutoff, :]
        # Normalize the each read at each position by the total coverage - This would collapse the variation between positons so should only be used in certain cases depending on the var slice used
        # adata.X = Normalizer().fit_transform(adata.X)
        # Scale the data - Seems to perform well with the robust scaler compared to the standard scaler
        adata.X = RobustScaler().fit_transform(adata.X)

        return adata
    
    def adataCluster(self, adata, neighbors_plot, neighbors_cluster, min_samples, min_cluster_size, n_components):
        # Apply a standardscaler to the data and reduce dimensions
        standard_embedding = umap.UMAP(random_state=self.randomstate, n_neighbors=neighbors_plot).fit_transform(adata.X)
        cluster_embedding = umap.UMAP(random_state=self.randomstate, n_neighbors=neighbors_cluster, min_dist=0.0, n_components=n_components).fit_transform(adata.X)
        # Perform clustering with HDBSCAN
        hdbscan_results = hdbscan.HDBSCAN(min_samples=min_samples, min_cluster_size=min_cluster_size).fit_predict(cluster_embedding)
        # Create a dataframe of the cluster information
        df = pd.DataFrame(standard_embedding, index=adata.obs.index, columns=['standard_umap1','standard_umap2'])
        df_c = pd.DataFrame(cluster_embedding, index=adata.obs.index, columns=['cluster_umap'+str(i) for i in range(1,n_components+1)])
        df = pd.concat([df, df_c], axis=1)
        df['cluster_hdbscan'] = hdbscan_results
        
        return df
    
    def adataCombine(self, adata, df, group):
        # Add the uns information to the original AnnData object for reference by using the 3: column of the dataframe
        adata.uns['_'.join([group,'cluster_umap'])] = df
        # Create dictionaries to map the cluster information to the original AnnData object obs for convience
        for i in [('cluster_hdbscan','cluster'), ('standard_umap1','umap1'), ('standard_umap2','umap2')]:
            # Create dictionaries to map the cluster information to the original AnnData object
            temp_dict = dict(zip(df.index, df[i[0]]))
            # Add the cluster information to the original AnnData object
            if group == 'sample':
                adata.obs['_'.join([group,i[1]])] = adata.obs.index.map(temp_dict)
            else:
                adata.obs['_'.join([group,i[1]])] = adata.obs['trna'].astype('str') + '_' + adata.obs[group].astype('str')
                adata.obs['_'.join([group,i[1]])] = adata.obs['_'.join([group,i[1]])].map(temp_dict)
        
        return adata

def main(args):
    '''
    Main function for running argparse and calling the appropriate class
    '''
        # Read database object or create one from trax run if none provided
    if args.mode == 'build':
        # Clean the path to the trax directory
        args.traxdir = os.path.abspath(args.traxdir)
        # Raise exception if trax directory is empty or doesn't exist
        if not os.path.isdir(args.traxdir):
            raise Exception('Error: trax directory does not exist.')
        # Raise exception if metadata file is empty or doesn't exist
        if args.metadata:
            if not os.path.isfile(args.metadata):
                raise Exception('Error: metadata file does not exist.')
        print('Building AnnData object...')
        # Create path to output directory if it doesn't exist
        print(toolsTG.builder(args.output))
        # Create AnnData object
        trax2anndata(args.traxdir, args.metadata, args.output).create()
        print('Done!\n')
    elif args.mode == 'merge':
        # Raise exception if h5ad file is empty or doesn't exist
        if not os.path.isfile(args.anndata1):
            raise Exception('Error: first h5ad file does not exist.')
        if not os.path.isfile(args.anndata2):
            raise Exception('Error: second h5ad file does not exist.')
        # Create path to output directory if it doesn't exist
        print(toolsTG.builder(args.output))
        # Merge AnnData objects
        print('Merging database objects...\n')
        anndataMerger(args).merge()
        print('Done!\n')
    elif args.mode == 'cluster':
        # Raise exception if h5ad file is empty or doesn't exist
        if not os.path.isfile(args.anndata):
            raise Exception('Error: h5ad file does not exist.')
        # Create path to output directory if it doesn't exist
        print(toolsTG.builder(args.output))
        print('Clustering data from database object...\n')
        anndataCluster(args).main()
        print('Done!\n')
    elif args.mode == 'graph':
        # Create output directory if it doesn't exist
        args.output = os.path.abspath(args.output)
        print(toolsTG.builder(args.output))
        # Raise exception if h5ad file is empty or doesn't exist
        if not os.path.isfile(args.anndata):
            raise Exception('Error: h5ad file does not exist.')
        print('Graphing data from database object...\n')
        anndataGrapher(args).main()
        print('Done!\n')
    elif args.mode == 'log2fc':
        # Raise exception if h5ad file is empty or doesn't exist
        if not os.path.isfile(args.anndata):
            raise Exception('Error: h5ad file does not exist.')
        # Load the AnnData object
        adata = ad.read_h5ad(args.anndata)
        # Load config file for name if specified
        config_name = 'default'
        if args.config:
            with open(args.config, 'r') as f:
                args.config = json.load(f)
            if 'name' in args.config:
                # self.args.output += '/' + self.args.config['name']
                config_name = args.config['name']
                # print(toolsTG.builder(self.args.output))
            else:
                raise ValueError('Config file must contain a "name" field')
        print('Calculating log2FC for database object...\n')
        adata_copy = adata.copy()
        log2FC_dict = adata.uns['log2FC']
        for readtype in [f'nreads_{i}_norm' for i in args.readtypes]:
            for cutoff in args.cutoff:
                toolsTG.adataLog2FC(adata_copy, args.group, readtype, readcount_cutoff=cutoff, config_name=config_name, overwrite=True).main()
        # if log2FC_dict.items() != adata_copy.uns['log2FC'].items(): # Can fix this to be more efficient later
        print('The log2FC uns dictionary has been updated.\nWriting h5ad database object to: ' + args.anndata)
        adata_copy.write(args.anndata)
        # else:
        # print('The log2FC uns dictionary has not been updated.\n')
        print('Done!\n')
    elif args.mode == 'csv':
        args.output = os.path.abspath(args.output)
        # Add the name of the h5ad file to the output directory minus the extension account for periods in the name removing just .h5ad
        args.output += '/' + '.'.join(os.path.basename(args.anndata).split('.')[:-1]) + '/'
        print(toolsTG.builder(args.output))
        adata = ad.read_h5ad(args.anndata)
        print('Writing csv files to: ' + args.output)
        adata.write_csvs(args.output, skip_data=False)
        print('Done!\n')
    else:
        print('Invalid operating mode. Exiting...')
        parser.print_help()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='tRNAgraph',
        description='tRNAgraph is a tool for analyzing tRNA-seq data generated from tRAX. It can be used to create an AnnData database object from \
            a trax output folder, or to analyze an existing database object and generate expanded visulizations. The database object can also be used to \
            perform further multivariate analysis such as clustering and classification of readcoverages.',
        allow_abbrev=False
    )

    subparsers = parser.add_subparsers(
        title='Operating modes',
        description='Choose between building a database object, mergeing two database objects together, running dimensionality reduction and clustering of coverage, \
            graphing data from an existing database object or using extra utilities to manipulate the database object.',
        dest='mode',
        required=True
    )

    # Tools parser
    parser_tools = subparsers.add_parser("tools", help="Extra utilities for working with tRNAgraph objects")
    tools_subparsers = parser_tools.add_subparsers(
        title='Operating modes',
        description='Extra utilities for working with tRNAgraph objects',
        dest='mode',
        required=True
    )

    # Build parser
    parser_build = subparsers.add_parser("build", help="Build a h5ad AnnData object from a tRAX run")
    parser_build.add_argument('-i', '--traxdir', help='Specify location of trax directory (required)', required=True)
    parser_build.add_argument('-m', '--metadata', help='Specify a metadata file to create annotations, you can also use the sample file used to generate tRAX DB (required)', required=True)
    parser_build.add_argument('-o', '--output', help='Specify output h5ad file (default: h5ad/trnagraph.h5ad) (optional)', default='h5ad/trnagraph.h5ad')
    parser_build.add_argument('--log', help='Log output to file (optional)', default=None)
    parser_build.add_argument('-q', '--quiet', help='Suppress output to stdout (optional)', action='store_true')

    # Merge parser
    parser_merge = subparsers.add_parser("merge", help="Merge data from two existing h5ad AnnData objects")
    parser_merge.add_argument('-i1', '--anndata1', help='Specify location of first h5ad object (required)', required=True)
    parser_merge.add_argument('-i2', '--anndata2', help='Specify location of second h5ad object (required)', required=True)
    parser_merge.add_argument('--dropno', help='Drop non tRNAs genes that are not present in both AnnData objects (optional)', action='store_true')
    parser_merge.add_argument('--droprna', help='Drop RNA categories that are not present in both AnnData objects (optional)', action='store_true')
    parser_merge.add_argument('-o', '--output', help='Specify output h5ad file (default: h5ad/trnagraph.merge.h5ad) (optional)', default='h5ad/trnagraph.merge.h5ad')
    parser_merge.add_argument('--log', help='Log output to file (optional)', default=None)
    parser_merge.add_argument('-q', '--quiet', help='Suppress output to stdout (optional)', action='store_true')

    # Cluster parser
    parser_cluster = subparsers.add_parser("cluster", help="Cluster data from an existing h5ad AnnData object")
    parser_cluster.add_argument('-i', '--anndata', help='Specify location of h5ad object (required)', required=True)
    parser_cluster.add_argument('-r', '--randomstate', help='Specify random state for UMAP if you want to have a static seed (default: None) (optional)', default=None, type=int)
    parser_cluster.add_argument('-t', '--readcutoff', help='Specify readcount cutoff to use for clustering (default: 20) (optional)', default=20, type=int)
    parser_cluster.add_argument('-v', '--coveragetype', help='Specify coverage types for umap clustering treated as features (default: uniquecoverage, readstarts, readends, mismatchedbases, deletions) (optional)', \
                                choices=['coverage', 'readstarts', 'readends', 'uniquecoverage', 'multitrnacoverage', 'multianticodoncoverage', 'multiaminocoverage','tRNAreadstotal', 'mismatchedbases', \
                                         'deletedbases', 'adenines', 'thymines', 'cytosines', 'guanines', 'deletions'], nargs='+', default=['uniquecoverage', 'readstarts', 'readends', 'mismatchedbases', 'deletions'])
    parser_cluster.add_argument('-c1', '--ncomponentsmp', help='Specify number of components to use for UMAP clustering of samples (default: 2) (optional)', default=2, type=int)
    parser_cluster.add_argument('-c2', '--ncomponentgrp', help='Specify number of components to use for UMAP clustering of groups (default: 2) (optional)', default=2, type=int)
    parser_cluster.add_argument('-l1', '--neighborclusmp', help='Specify number of neighbors to use for UMAP clustering of samples (default: 150) (optional)', default=150, type=int)
    parser_cluster.add_argument('-l2', '--neighborclusgrp', help='Specify number of neighbors to use for UMAP clustering of groups (default: 40) (optional)', default=40, type=int)
    parser_cluster.add_argument('-n1', '--neighborstdsmp', help='Specify number of neighbors to use for UMAP projection plotting of samples (default: 75) (optional)', default=75, type=int)
    parser_cluster.add_argument('-n2', '--neighborstdgrp', help='Specify number of neighbors to use for UMAP projection plotting of groups (default: 20) (optional)', default=20, type=int)
    parser_cluster.add_argument('-d1', '--hdbscanminsampsmp', help='Specify minsamples size to use for HDBSCAN clustering of samples (default: 6) (optional)', default=6, type=int)
    parser_cluster.add_argument('-d2', '--hdbscanminsampgrp', help='Specify minsamples size to use for HDBSCAN clustering of groups (default: 3) (optional)', default=3, type=int)
    parser_cluster.add_argument('-b1', '--hdbscanminclusmp', help='Specify min cluster size to use for HDBSCAN clustering of samples (default: 30) (optional)', default=30, type=int)
    parser_cluster.add_argument('-b2', '--hdbscanminclugrp', help='Specify min cluster size to use for HDBSCAN clustering of groups (default: 10) (optional)', default=10, type=int)
    parser_cluster.add_argument('--clusterobsexperimental', help='This is an experimental feature to add columns from adata.obs to the adata.var and adata.X to be used for clustering (optional)', nargs='+', default=[])
    parser_cluster.add_argument('-w', '--overwrite', help='Overwrite existing cluster information in AnnData object (optional)', action='store_true')
    parser_cluster.add_argument('-o', '--output', help='Specify output directory (default: h5ad/trnagraph.cluster.h5ad) (optional)', default='h5ad/trnagraph.cluster.h5ad')
    parser_cluster.add_argument('--log', help='Log output to file (optional)', default=None)
    parser_cluster.add_argument('-q', '--quiet', help='Suppress output to stdout (optional)', action='store_true')

    # Graph parser
    parser_graph = subparsers.add_parser("graph", help="Graph data from an existing h5ad AnnData object")
    parser_graph.add_argument('-i', '--anndata', help='Specify location of h5ad object (required)', required=True)
    parser_graph.add_argument('-o', '--output', help='Specify output directory (optional)', default='figures')
    parser_graph.add_argument('-g', '--graphtypes', choices=['all','bar','cluster','compare','correlation','count','coverage','heatmap','logo','pca','radar','volcano'], \
                              help='Specify graphs to create, if not specified it will default to "all" (optional)', nargs='+', default='all')
    # Add argument to filter parameters from AnnData object
    parser_graph.add_argument('--config', help='Specify a json file containing observations/variables to filter out and other config options (optional)', default=None)
    parser_graph.add_argument('--colormap', help='Specify a json file containing colormaps for the graphs (optional)', default=None)
    # Options to imporve speed or log output
    parser_graph.add_argument('--regen_uns', help='Force regenerate uns log2fc data if it would be generated again (optional)', action='store_true', default=False)
    parser_graph.add_argument('-n', '--threads', help='Specify number of threads to use (default: cpu_max) (optional)', default=0, type=int)
    parser_graph.add_argument('--log', help='Log output to file (optional)', default=None)
    parser_graph.add_argument('-q', '--quiet', help='Suppress output to stdout (optional)', action='store_true')
    parser_graph.add_argument('-v', '--verbose', help='Print verbose output to stdout (optional)', action='store_true')
    # Bar options
    parser_graph.add_argument('--barcol', help='Specify AnnData column to of what the individal stacks of bars will be (default: group) (optional)', default='group', required=False)
    parser_graph.add_argument('--bargrp', help='Specify AnnData column to of what will stack in bar columns (default: amino) (optional)', default='amino', required=False)
    parser_graph.add_argument('--barsubgrp', help='Specify AnnData column for secondary spliting of bars into subplots (default: None) (optional)', default=None, required=False)
    parser_graph.add_argument('--barsort', help='Specify AnnData column to sort the bars by (default: None) (optional)', default=None, required=False)
    parser_graph.add_argument('--barlabel', help='Specify wether to label the bars using a different AnnData column (default: None) (optional)', default=None, required=False)
    # Cluster options
    parser_graph.add_argument('--clustergrp', help='Specify AnnData column to group by (default: amino) (optional)', default='amino', required=False)
    parser_graph.add_argument('--clusterlabels', help='Specify a AnnData column of names to use for the clusters instead of the default and will place them on the plot (optional)', default=None, required=False)
    parser_graph.add_argument('--clusteroverview', help='Specify wether to generate an overview of the clusters (default: False) (optional)', default=False, action='store_true', required=False)
    parser_graph.add_argument('--clusternumeric', help='Specify wether to the cluster category is numeric (default: False) (optional)', default=False, action='store_true', required=False)
    parser_graph.add_argument('--clustermask', help='Specify wether to mask the cluster plots to annotated HDBSCAN clusters (default: False) (optional)', default=False, action='store_true', required=False)
    # Compare options
    parser_graph.add_argument('--comparegrp1', help='Specify AnnData column as main comparative group (default: group) (optional)', default='group', required=False)
    parser_graph.add_argument('--comparegrp2', help='Specify AnnData column to group by (default: group) (optional)', default='group', required=False)
    # Correlation options
    parser_graph.add_argument('--corrmethod', choices=['pearson', 'spearman', 'kendall'], help='Specify correlation method (default: pearson) (optional)', default='pearson', required=False)
    parser_graph.add_argument('--corrgroup', help='Specify a grouping variable to generate correlation matrices for (default: sample) (optional)', default='sample', required=False)
    # Coverage options
    parser_graph.add_argument('--covgrp', help='Specify a grouping variable to generate coverage plots for (default: group) (optional)', default='group', required=False)
    parser_graph.add_argument('--covobs', help='Specify the basis for each individual coverage plot (default: trna) (optional)', default='trna', required=False)
    parser_graph.add_argument('--covtype', help='Specify a coverage type for coverage plots corresponding to trax coverage file outputs (default: uniquecoverage) (optional)', \
                              choices=['coverage', 'readstarts', 'readends', 'uniquecoverage', 'multitrnacoverage', 'multianticodoncoverage', 'multiaminocoverage','tRNAreadstotal', 'mismatchedbases', \
                                       'deletedbases', 'adenines', 'thymines', 'cytosines', 'guanines', 'deletions'], default='uniquecoverage')
    parser_graph.add_argument('--covgap', help='Specify wether to include gaps in coverage plots (default: False) (optional)', default=False)
    parser_graph.add_argument('--covmethod', help='Specify method to use for coverage plots when combining multiple groups (default: mean) (optional)', choices=['mean','median','max','min','sum'], default='mean', required=False)
    parser_graph.add_argument('--combinedpdfonly', help='Do not generate single tRNA coverage plot PDFs for every tRNA, only keep the combined output (optional)', action='store_true', required=False)
    # Heatmap options
    parser_graph.add_argument('--heatgrp', help='Specify group to use for heatmap', default='group', required=False)
    parser_graph.add_argument('--diffrts', choices=['wholecounts_unique', 'fiveprime_unique', 'threeprime_unique', 'other_unique', 'total_unique', \
                                                    'wholecounts', 'fiveprime', 'threeprime', 'other', 'total', 'all'], \
                             help='Specify readtypes to use for heatmap/volcano (default: wholecounts_unique, fiveprime_unique, threeprime_unique, other_unique, total_unique) (optional)', \
                             nargs='+', default=['wholecounts_unique', 'fiveprime_unique', 'threeprime_unique', 'other_unique', 'total_unique'], required=False)
    parser_graph.add_argument('--heatcutoff', help='Specify readcount cutoff to use for heatmap', default=80, required=False, type=int)
    parser_graph.add_argument('--heatbound', help='Specify range to use for bounding the heatmap to top and bottom counts', default=25, required=False)
    parser_graph.add_argument('--heatsubplots', help='Specify wether to generate subplots for each comparasion in addition to the sum (default: False)', action='store_true', default=False, required=False)
    # PCA options
    parser_graph.add_argument('--pcamarkers', help='Specify AnnData column to use for PCA markers (default: sample) (optional)', default='sample')
    parser_graph.add_argument('--pcacolors', help='Specify AnnData column to color PCA markers by (default: group) (optional)', default='group')
    parser_graph.add_argument('--pcareadtypes', choices=['wholecounts_unique', 'fiveprime_unique', 'threeprime_unique', 'other_unique', 'total_unique', \
                                                         'wholecounts', 'fiveprime', 'threeprime', 'other', 'total', 'all'], \
                             help='Specify read types to use for PCA markers (default: total_unique, total) (optional)', nargs='+', default=['total_unique', 'total'])
    # Radar options
    parser_graph.add_argument('--radargrp', help='Specify AnnData column to group by (default: group) (optional)', default='group', required=False)
    parser_graph.add_argument('--radarmethod', help='Specify method to use for radar plots (default: mean) (optional)', choices=['mean','median','max','sum','all'], default=['mean'], nargs='+', required=False)
    parser_graph.add_argument('--radarscaled', help='Specify wether to scale the radar plots to 100%% (optional)', action='store_true', default=False, required=False)
    # Seqlogo options
    parser_graph.add_argument('--logogrp', help='Specify AnnData column to group sequences by (default: amino) (optional)', default='amino', required=False)
    parser_graph.add_argument('--logomanualgrp', help='Specify a manual group of tRNAs to use for seqlogo plots instead of using the AnnData column (optional)', nargs='+', default=None)
    parser_graph.add_argument('--logomanualname', help='Specify a name for the manual group of tRNAs output file, will be ignored and timestamped if not specified (optional)', default=None)
    parser_graph.add_argument('--logopseudocount', help='Specify the number of pseudocounts to add to each position when calculating as ratio of the bases in the pool of RNAs (default: 20) (optional)', default=20, required=False, type=int)
    parser_graph.add_argument('--logosize', help='Specify the sequence size to use for the logo plots from presets (default: noloop)', choices=['sprinzl', 'noloop', 'full'], default='noloop', required=False)
    parser_graph.add_argument('--ccatail', help='Specify wether to keep the CCA tail from the sequences (optional)', action='store_false', default=True, required=False)
    parser_graph.add_argument('--pseudogenes', help='Specify wether to keep the pseudo-tRNAs (tRX) (optional)', action='store_false', default=True, required=False)
    parser_graph.add_argument('--logornamode', help='Specify wether to print the output as RNA rather than DNA (optional)', action='store_true', default=False, required=False)
    # Volcano options
    parser_graph.add_argument('--volgrp', help='Specify group to use for volcano plot', default='group', required=False)
    parser_graph.add_argument('--volcutoff', help='Specify readcount cutoff to use for volcano plot', default=80, required=False)
    # Log2fc parser
    parser_tools_log2fc = tools_subparsers.add_parser("log2fc", help="Compute log2fc data from an existing h5ad AnnData object")
    parser_tools_log2fc.add_argument('-i', '--anndata', help='Specify location of h5ad object (required)', required=True)
    parser_tools_log2fc.add_argument('-g', '--group', help='Specify group to use for log2fc from obs (default: group) (optional)', default='group', required=False)
    parser_tools_log2fc.add_argument('-r', '--readtypes', choices=['wholecounts_unique', 'fiveprime_unique', 'threeprime_unique', 'other_unique', 'total_unique', \
                                                                   'wholecounts', 'fiveprime', 'threeprime', 'other', 'total', 'all'], \
                                    help='Specify readtypes to generate log2fc for (default: wholecounts_unique, fiveprime_unique, threeprime_unique, other_unique, total_unique) (optional)', \
                                    nargs='+', default=['wholecounts_unique', 'fiveprime_unique', 'threeprime_unique', 'other_unique', 'total_unique'], required=False)
    parser_tools_log2fc.add_argument('-x', '--cutoff', help='Specify readcounts cutoff to use for log2fc (default: 80) (optional)', default=80, required=False, type=int, nargs='+')
    parser_tools_log2fc.add_argument('-c', '--config', help='Specify a json file containing observations/variables to filter out and other config options (optional)', default=None)
    parser_tools_log2fc.add_argument('--log', help='Log output to file (optional)', default=None)
    parser_tools_log2fc.add_argument('-q', '--quiet', help='Suppress output to stdout (optional)', action='store_true')
    # CSV parser
    parser_tools_csv = tools_subparsers.add_parser("csv", help="Output .h5ad to CSV")
    parser_tools_csv.add_argument('-i', '--anndata', help='Specify location of h5ad object (required)', required=True)
    parser_tools_csv.add_argument('-o', '--output', help='Specify output directory (optional)', default='csv')
    parser_tools_csv.add_argument('--log', help='Log output to file (optional)', default=None)
    parser_tools_csv.add_argument('-q', '--quiet', help='Suppress output to stdout (optional)', action='store_true')

    args = parser.parse_args()

    # Set log file if specified
    sys.stdout = open(args.log, 'w') if args.log else sys.stdout
    # Run main function
    if args.quiet:
        with open(os.devnull, 'w') as f, contextlib.redirect_stdout(f):
            main(args)
    else:
        main(args)