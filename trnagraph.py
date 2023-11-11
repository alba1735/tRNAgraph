#!/usr/bin/env python3

import numpy as np
import pandas as pd
import anndata as ad
import argparse
import os
import sys
import json
# Custom functions
import directory_tools
import coverage_tools
import heatmap_tools
import pca_tools
import correlation_tools
import volcano_tools
import radar_tools
import bar_tools
import compare_tools
# Cluster functions
import umap
from sklearn.preprocessing import RobustScaler
from sklearn.preprocessing import Normalizer
import hdbscan
# Plotting functions
import matplotlib.pyplot as plt
import seaborn as sns

class trax2anndata():
    '''
    Create h5ad AnnData object from a trax run
    '''
    def __init__(self, traxdir, metadata, observations, output):
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
        # For adding unique counts to coverage file
        trnauniquecounts = traxdir + '/unique/' + traxdir.split('/')[-1] + '-trnauniquecounts.txt'
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
        # For adding amino acid counts to coverage file
        aminoacidcounts = traxdir + '/' + traxdir.split('/')[-1] + '-aminocounts.txt'
        self.amino_counts = pd.read_csv(aminoacidcounts, sep='\t')
        # For adding metadata to adata object
        metadata_type = '\t' if metadata.endswith('.tsv') else ',' if metadata.endswith('.csv') else ' '
        self.metadata = pd.read_csv(metadata, sep=metadata_type, header=None, index_col=None)
        # For adding categories to adata object

        self.observations = observations
        if self.observations:
            # Make sure that observations are not going to be generated automatically
            auto_obs = ['trna', 'iso', 'amino', 'sample', 'group', 'deseq2_sizefactor', 'refseq', 'dataset']
            if any(x in self.observations for x in auto_obs):
                raise ValueError('The following observation categories will automatically be generated please remove these if you included them: {}'.format(auto_obs))
            # Make sure that observations are unique
            if len(self.observations) != len(set(self.observations)):
                raise ValueError('Observation categories must be unique, please remove duplicates from the observation catgories you wish to generate: {}'.format(self.observations))
            self.observations.insert(0,'sample')
            self.observations.insert(1,'group')
            # Add manual observations to obs list if they are not provided or if the length of the observations list does not match the number of parameters in the trax coverage file
            if len(self.observations) != len(self.metadata.columns):
                diff_obs_count = len(self.metadata.columns)-len(self.observations)
                print('Number of observations does not match number of parameters in trax coverage file by {}. To create a more specific database object, please provide the correct number of observations.'.format(diff_obs_count))
                if diff_obs_count > 0:
                    print('Adding {} observations to the end of the list'.format(diff_obs_count))
                    self.observations += ['obs_' + str(x) for x in range(diff_obs_count)]
                if diff_obs_count < 0:
                    print('Removing {} observations from the end of the list'.format(abs(diff_obs_count)))
                    self.observations = self.observations[:diff_obs_count]
        else:
            diff_obs_count = abs(len(self.metadata.columns)-2)
            self.observations = ['sample', 'group'] + ['obs_' + str(x) for x in range(diff_obs_count)]
        # Add align observation categories to metadata
        self.metadata.columns = self.observations
        self.metadata.set_index('sample', inplace=True)
        self.metadata = self.metadata.to_dict()
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
        # Build reference sequence dataframe
        self._ref_seq_build_()
        # Build obs and x dataframes
        x_dfs = []
        for cov_type in self.cov_types:
            x_df = self._x_build_(cov_type)
            x_dfs.append(x_df)
        obs_df = self._obs_build_(x_df)
        x_df = pd.concat(x_dfs,axis=1,sort=False)
        # Rename columns of x_df to include position and coverage type
        clist = [[p + '_' + cov for p in self.positions] for cov in self.cov_types]
        clist = [a for l in clist for a in l]
        x_df.columns = clist
        # obs_df,x_df = self._group_sort_(obs_df,x_df) # Not sure if I need this function
        # Check that the index of the obs and x dataframes are the same
        if not obs_df.index.equals(x_df.index):
            raise ValueError('The index of the obs and x dataframes are not the same. This means somthing went wrong in the sorting process.')
        # Build adata object
        adata = self._adata_build_(obs_df,x_df)
        # Add size factors to adata object as raw layer
        adata.layers['raw'] = adata.X*adata.obs['deseq2_sizefactor'].values[:,None]
        # Quality check adata by dropping NaN values and printing summary
        if adata.obs.isna().any(axis=0).any():
            print('WARNING: NaN values found in obs dataframe this is commonly caused by missing samples in your metadata file or have a different number of observations per sample. ' \
                  'Another consideration is to make sure your .tsv/.csv is using tabs or commas appropriately. ' \
                  'It can also be caused by metadata containing NaN or None values. Please check your metadata file to make sure the following are correct:\n{}'.format(adata.obs.columns[adata.obs.isna().any(axis=0)].tolist()))
        # Add output name to adata object index
        adata.obs.index = [self.output.split('.')[0] + '_' + str(x) for x in adata.obs.index] 
        # Save adata object
        adata.write(f'{self.output}')
        print(f'Writing h5ad database object to {self.output}')
        
    def _ref_seq_build_(self):
        '''
        Build reference sequence from trax coverage file
        '''
        x_df = self._x_build_('actualbase')
        seqs = x_df.values.tolist()
        self.seqs = [''.join(x) for x in seqs]

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
        obs_df['trna_type'] = trna_obs[0]
        obs_df['amino'] = trna_obs[1]
        obs_df['iso'] = trna_obs[2]
        obs_df['refseq'] = self.seqs
        obs_df['sizefactor'] = [self.size_factors.get(i) for i in ['_'.join(x.split('_')[1:]) for x in x_df.index.values]]
        # obs_df['nreads_unique_raw'] = [self.unique_counts.get(i[0]).get('_'.join(i[1:])) if self.unique_counts.get(i[0]) else 0 for i in [x.split('_') for x in x_df.index.values]] # Some samples may not have any reads for a given tRNA in the unique_counts dictionary might want to double check this
        # Create unique counts for each tRNA split by type
        for rt in ['fiveprime','threeprime','whole','other']:
            obs_df[f'nreads_{rt}_unique_raw'] = [self.unique_counts.get(i[0]+f'_{rt}').get('_'.join(i[1:])) if self.unique_counts.get(i[0]+f'_{rt}') else 0 for i in [x.split('_') for x in x_df.index.values]]
            obs_df[f'nreads_{rt}_unique_norm'] = obs_df[f'nreads_{rt}_unique_raw']/obs_df['sizefactor']
        # Add total unique read counts to obs dataframe
        obs_df['nreads_total_unique_raw'] = obs_df[[f'nreads_{rt}_unique_raw' for rt in ['fiveprime','threeprime','whole','other']]].sum(axis=1)
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
        gap_list = [True if len(i.split('.')) > 1 else False for i in positions]
        gap_list = [True if positions[i][-1] == 'a' else gap_list[i] for i in range(len(gap_list))]
        gap_list = [True if positions[i][-1] == 'b' else gap_list[i] for i in range(len(gap_list))]
        gap_list = [True if positions[i][0] == 'e' else gap_list[i] for i in range(len(gap_list))]
        # Build AnnData object
        adata = ad.AnnData(x_df, dtype=x_df.dtypes[0]) # Not sure if this is the dtype I want to use defaults to float64
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
        # Add custom dataframe obs
        for ob in self.observations:
            adata.obs[ob] = obs_df[ob].values
        # Add the numer of reads per tRNA as observations-annotation to adata from trna unique counts file
        for rt in ['whole','fiveprime','threeprime','other']:
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
        # Add anticodon counts as uns
        adata.uns['anticodon_counts'] = self.anticodon_counts
        # Add amino acid counts as uns
        adata.uns['amino_counts'] = self.amino_counts
        # Add type counts as uns
        adata.uns['type_counts'] = self.type_counts
        # Add non tRNA counts as uns
        adata.uns['nontRNA_counts'] = self.non_trna_read_counts
        
        return adata

class anndataGrapher:
    '''
    Class to generate graphs from an AnnData object by calling the appropriate graphing functions
    '''
    def __init__(self, args):
        self.adata = ad.read_h5ad(args.anndata)
        self.args = args

    def graph(self):
        if self.args.graphtypes == 'all' or 'all' in self.args.graphtypes:
            self.args.graphtypes = ['coverage', 'heatmap', 'pca', 'correlation', 'volcano', 'radar', 'bar']

        if self.args.config:
            with open(self.args.config, 'r') as f:
                d_config = json.load(f)
            if 'name' in d_config:
                self.args.output += '/' + d_config['name']
                directory_tools.builder(self.args.output)
            if 'obs' in d_config:
                # Dictionary of uns columns and values to filter by as groups and samples since the coulmns are different from the main adata obs
                uns_dict = {i:True for i in self.adata.uns['amino_counts'].columns.values}
                uns_dict.update({i:True for i in self.adata.uns['type_counts'].columns.values})
                for k,v in d_config['obs'].items():
                    print('Filtering AnnData object by observation: ' + k + ' , ' + str(v))
                    # Filter all uns columns by the observation and update the uns_dict
                    sub_uns_dict = dict(zip(self.adata.obs['sample'], self.adata.obs[k]))
                    sub_uns_dict.update(dict(zip(self.adata.obs['group'], self.adata.obs[k])))
                    uns_dict = {i:False if sub_uns_dict.get(i,False) not in v else True for i,j in uns_dict.items()}
                    # Filter the adata object by the observation
                    self.adata = self.adata[self.adata.obs[k].isin(v), :]
                # Filter uns columns by the uns_dict
                for uns in self.adata.uns:
                    self.adata.uns[uns] = self.adata.uns[uns].loc[:, [i for i in self.adata.uns[uns].columns.values if uns_dict[i]]]
            if 'var' in d_config:
                for k,v in d_config['var'].items():
                    print('Filtering AnnData object by variable: ' + k + ' , ' + str(v))
                    self.adata = self.adata[:, self.adata.var[k].isin(v)]
        else:
            d_config = {}

        if 'coverage' in self.args.graphtypes:
            output = self.args.output + '/coverage'
            directory_tools.builder(output)

            colormap = None
            if 'colormap' in d_config:
                if self.args.coveragegrp in d_config['colormap']:
                    colormap = d_config['colormap'][self.args.coveragegrp]

            if self.args.combineonly:
                print('Generating combined coverage plots...')
                coverage_tools.visualizer(self.adata.copy(), self.args.threads, self.args.coveragegrp, self.args.coverageobs, self.args.coveragetype, self.args.coveragegap, colormap, output).generate_combine()
            else:
                print('Generating individual coverage plots...')
                directory_tools.builder(output+'/single')
                directory_tools.builder(output+'/single/low_coverage')
                coverage_tools.visualizer(self.adata.copy(), self.args.threads, self.args.coveragegrp, self.args.coverageobs, self.args.coveragetype, self.args.coveragegap, colormap, output).generate_split()
                print('Generating combined coverage plots...')
                coverage_tools.visualizer(self.adata.copy(), self.args.threads, self.args.coveragegrp, self.args.coverageobs, self.args.coveragetype, self.args.coveragegap, colormap, output).generate_combine()
            print('Coverage plots generated.\n')

        if 'heatmap' in self.args.graphtypes:
            print('Generating heatmaps...')
            output = self.args.output + '/heatmap'
            directory_tools.builder(output)
            heatmap_tools.visualizer(self.adata.copy(), self.args.heatgrp, self.args.heatrts, self.args.heatcutoff, self.args.heatbound, self.args.heatsubplots, output)
            print('Heatmaps generated.\n')

        if 'pca' in self.args.graphtypes:
            colormap = None
            if 'colormap' in d_config:
                if self.args.pcacolors in d_config['colormap']:
                    colormap = d_config['colormap'][self.args.pcacolors]

            print('Generating pca plots...')
            output = self.args.output + '/pca'
            directory_tools.builder(output)
            pca_tools.visualizer(self.adata.copy(), self.args.pcamarkers, self.args.pcacolors, self.args.pcareadtypes, colormap, output)
            print('PCA plots generated.\n')

        if 'correlation' in self.args.graphtypes:
            print('Generating correlation plots...')
            output = self.args.output + '/correlation'
            directory_tools.builder(output)
            correlation_tools.visualizer(self.adata.copy(), output, self.args.corrmethod, self.args.corrgroup)
            print('Correlation plots generated.\n')

        if 'volcano' in self.args.graphtypes:
            print('Generating volcano plots...')
            output = self.args.output + '/volcano'
            directory_tools.builder(output)
            volcano_tools.visualizer(self.adata.copy(), self.args.volgrp, self.args.volrt, self.args.volcutoff, output)
            print('Volcano plots generated.\n')

        if 'radar' in self.args.graphtypes:
            colormap = None
            if 'colormap' in d_config:
                if self.args.radargrp in d_config['colormap']:
                    colormap = d_config['colormap'][self.args.radargrp]

            print('Generating radar plots...')
            output = self.args.output + '/radar'
            directory_tools.builder(output)
            radar_tools.visualizer(self.adata.copy(), self.args.radargrp, colormap, output)
            print('Radar plots generated.\n')

        if 'bar' in self.args.graphtypes:
            colormap_tc, colormap_bg = None, None
            if 'colormap' in d_config:
                if self.args.bargrp in d_config['colormap']:
                    colormap_bg = d_config['colormap'][self.args.bargrp]
                if 'group' in d_config['colormap']:
                    colormap_tc = d_config['colormap']['group']

            print('Generating bar plots...')
            output = self.args.output + '/bar'
            directory_tools.builder(output)
            bar_tools.visualizer(self.adata.copy(), self.args.bargrp, colormap_tc, colormap_bg, output)
            print('Bar plots generated.\n')

        if 'compare' in self.args.graphtypes:
            colormap = None
            if 'colormap' in d_config:
                if self.args.comparegrp1 in d_config['colormap']:
                    colormap = d_config['colormap'][self.args.comparegrp1]

            print('Generating comparison plots...')
            output = self.args.output + '/compare'
            directory_tools.builder(output)
            compare_tools.visualizer(self.adata.copy(), self.args.comparegrp1, self.args.comparegrp2, colormap, output)
            print('Comparison plots generated.\n')

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
        self.randomstate = args.randomstate
        self.neighbors_cluster = args.neighborsclu
        self.neighbors_plot = args.neighborsstd
        self.output = args.output + '/clustering'
        directory_tools.builder(self.output)

    def main(self):
        print('Performing UMAP clustering...')
        # Preprocess AnnData object
        adata_sub = self.adataPreprocess(self.adata.copy())
        # Apply a standardscaler to the data and reduce dimensions
        standard_embedding = umap.UMAP(random_state=self.randomstate, n_neighbors=self.neighbors_plot).fit_transform(adata_sub.X)
        cluster_embedding = umap.UMAP(random_state=self.randomstate, n_neighbors=self.neighbors_cluster, min_dist=0.0, n_components=12).fit_transform(adata_sub.X)
        # Perform clustering with HDBSCAN
        cluster_labels = hdbscan.HDBSCAN(min_samples=5, min_cluster_size=40).fit_predict(cluster_embedding)
        clustered = (cluster_labels >= 0)
        # Save the clustered data to csv
        adata_sub.obs['umap1'] = standard_embedding[:,0]
        adata_sub.obs['umap2'] = standard_embedding[:,1]
        adata_sub.obs['cluster'] = cluster_labels
        adata_sub.obs['clustered'] = clustered
        adata_sub.obs.to_csv(self.output + '/umap.csv')
        # Plot the UMAP projection
        self.clusterPlot(standard_embedding, cluster_labels, clustered, adata_sub)

    def adataPreprocess(self, adata):
        '''
        Preprocess AnnData object for clustering
        '''
        # Clean up the data by removing gaps from var and subsetting the var to uniquecoverage, readstarts, readends, mismatchedbases, and deletions
        adata = adata[:, ~adata.var['gap']]
        adata = adata[:, adata.var['coverage'].isin(['uniquecoverage', 'readstarts', 'readends', 'mismatchedbases', 'deletions'])]
        # Filter out Und samples from that amino acid category
        adata = adata[~(adata.obs['amino'] == 'Und'), :]
        # Filter out samples with low coverage
        adata = adata[adata.obs['nreads_total_unique_raw'] > 20, :]
        # Normalize the each read at each position by the total coverage
        # adata.X = Normalizer().fit_transform(adata.X)
        # Scale the data
        adata.X = RobustScaler().fit_transform(adata.X)

        return adata
    
    def clusterPlot(self, embedding, cluster_labels, clustered, adata):
        # Create a 2 x 2 subplot with the umap projection and the cluster labels as the last subplot
        fig, axs = plt.subplots(3, 3, figsize=(20,20))
        # Plot first through eigth subplots
        plot_list = [('Amino Acid','amino',0,0), ('Tissue Type','tissuetype',0,1), ('Treatment','treatment',0,2), \
                     ('Timepoint','timepoint',1,0), ('Fragment Type','fragment',1,1), ('Group','group',1,2), \
                     ('Sample','sample',2,0), ('Dataset','dataset',2,1)]
        for i in plot_list:
            pal = sns.color_palette("hls", len(pd.unique(adata.obs[i[1]][clustered])))
            sns.scatterplot(x=embedding[~clustered, 0], y=embedding[~clustered, 1], s=1, ax=axs[i[2],i[3]], hue=adata.obs[i[1]][~clustered], palette=pal)
            axs[i[2],i[3]].set_title(i[0])
            axs[i[2],i[3]].set_xlabel('UMAP 1')
            axs[i[2],i[3]].set_ylabel('UMAP 2')
            axs[i[2],i[3]].set_xticks([])
            axs[i[2],i[3]].set_yticks([])
            # Create legend from pal adding 
            axs[i[2],i[3]].legend(pd.unique(adata.obs[i[1]]), loc='upper right')
        # Plot ninth subplot colored by cluster
        pal = sns.color_palette("hls", len(pd.unique(cluster_labels[clustered])))
        sns.scatterplot(x=embedding[~clustered, 0], y=embedding[~clustered, 1], s=1, ax=axs[2,2], c=(0.5,0.5,0.5), alpha=0.5)
        sns.scatterplot(x=embedding[clustered, 0], y=embedding[clustered, 1], s=1, ax=axs[2,2], hue=cluster_labels[clustered], palette=pal)
        axs[2,2].set_title('HDBSCAN')
        axs[2,2].set_xlabel('UMAP 1')
        axs[2,2].set_ylabel('UMAP 2')
        axs[2,2].set_xticks([])
        axs[2,2].set_yticks([])
        axs[2,2].legend(pd.unique(cluster_labels[clustered]), loc='upper right')
        # Set aspect ratio to 'equal' so that umap projection is not distorted
        plt.gca().set_aspect('equal', 'datalim')
        # Add title
        plt.title('UMAP projection of the tRNA dataset')
        # Save figure
        plt.tight_layout()
        plt.savefig(self.output + '/umap.pdf')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='tRNAgraph',
        description='tRNAgraph is a tool for analyzing tRNA-seq data generated from tRAX. It can be used to create a AnnData database object from \
            a trax coverage file, or to analyze an existing database object and generate expanded visulizations. The database object can also be used to \
            perform further analysis beyond the scope of what tRAX can do.',
        allow_abbrev=False
    )
    # parser.add_argument('-o', '--output', help='Specify output directory', default='trnagraph')

    subparsers = parser.add_subparsers(
        title='Operating modes',
        description='Choose between building a database object, mergeing two database objects together, or graphing data from an existing database object',
        dest='mode',
        required=True
    )

    parser_build = subparsers.add_parser("build", help="Build a h5ad AnnData object from a tRAX run")
    parser_build.add_argument('-i', '--traxdir', help='Specify location of trax directory (required)', required=True)
    parser_build.add_argument('-m', '--metadata', help='Specify a metadata file to create annotations, you can also use the sample file used to generate tRAX DB (required)', required=True)
    parser_build.add_argument('-l', '--observationslist', help='Specify the observations of sample names in order (optional)', nargs='*', default=None)
    parser_build.add_argument('-f', '--observationsfile', help='Specify a file containing the observations of sample names in order as tab seperated file (optional)', default=None)
    parser_build.add_argument('-o', '--output', help='Specify output h5ad file (default: trnagraph.h5ad) (optional)', default='trnagraph.h5ad')
    parser_build.add_argument('--log', help='Log output to file (optional)', default=None)

    parser_merge = subparsers.add_parser("merge", help="Merge data from two existing h5ad AnnData objects")
    parser_merge.add_argument('-i1', '--anndata1', help='Specify location of first h5ad object (required)', required=True)
    parser_merge.add_argument('-i2', '--anndata2', help='Specify location of second h5ad object (required)', required=True)
    parser_merge.add_argument('--dropno', help='Drop non tRNAs genes that are not present in both AnnData objects (optional)', action='store_true')
    parser_merge.add_argument('--droprna', help='Drop RNA categories that are not present in both AnnData objects (optional)', action='store_true')
    parser_merge.add_argument('-o', '--output', help='Specify output h5ad file (default: trnagraph_merge.h5ad) (optional)', default='trnagraph_merge.h5ad')
    parser_merge.add_argument('--log', help='Log output to file (optional)', default=None)

    parser_cluster = subparsers.add_parser("cluster", help="Cluster data from an existing h5ad AnnData object")
    parser_cluster.add_argument('-i', '--anndata', help='Specify location of h5ad object (required)', required=True)
    parser_cluster.add_argument('-n1', '--neighborsstd', help='Specify number of neighbors to use for UMAP Plotting (default: 25) (optional)', default=25, type=int)
    parser_cluster.add_argument('-n2', '--neighborsclu', help='Specify number of neighbors to use for UMAP Clustering (default: 200) (optional)', default=100, type=int)
    parser_cluster.add_argument('-r', '--randomstate', help='Specify random state for UMAP (default: 19) (optional)', default=19, type=int)
    parser_cluster.add_argument('-o', '--output', help='Specify output directory (optional)', default='trnagraph')
    parser_cluster.add_argument('--log', help='Log output to file (optional)', default=None)

    parser_graph = subparsers.add_parser("graph", help="Graph data from an existing h5ad AnnData object")
    parser_graph.add_argument('-i', '--anndata', help='Specify location of h5ad object (required)', required=True)
    parser_graph.add_argument('-o', '--output', help='Specify output directory (optional)', default='trnagraph')
    parser_graph.add_argument('-g', '--graphtypes', choices=['all','coverage','heatmap','pca','correlation','volcano','radar','bar','compare'], \
        help='Specify graphs to create, if not specified it will default to "all" (optional)', nargs='+', default='all')
    # Add argument to filter parameters from AnnData object
    parser_graph.add_argument('--config', help='Specify a json file containing observations/variables to filter out and other config options (optional)', default=None)
    # Options to imporve speed or log output
    parser_graph.add_argument('-n', '--threads', help='Specify number of threads to use (default: 1) (optional)', default=1, type=int)
    parser_graph.add_argument('--log', help='Log output to file (optional)', default=None)
    # Coverage options
    parser_graph.add_argument('--coveragegrp', help='Specify a grouping variable to generate coverage plots for (default: group) (optional)', default='group')
    parser_graph.add_argument('--coverageobs', help='Specify a observation subsetting for coverage plots (optional)', nargs='+', default=None)
    parser_graph.add_argument('--coveragetype', help='Specify a coverage type for coverage plots corresponding to trax coverage file outputs (default: uniquecoverage) (optional)', \
        choices=['coverage', 'readstarts', 'readends', 'uniquecoverage', 'multitrnacoverage', 'multianticodoncoverage', 'multiaminocoverage','tRNAreadstotal', 'mismatchedbases', \
            'deletedbases', 'adenines', 'thymines', 'cytosines', 'guanines', 'deletions'], default='uniquecoverage')
    parser_graph.add_argument('--coveragegap', help='Specify wether to include gaps in coverage plots (default: False) (optional)', default=False)
    parser_graph.add_argument('--combineonly', help='Do not generate single tRNA coverage plot PDFs for every tRNA, only keep the combined output (optional)', action='store_true', required=False)
    # Heatmap options
    parser_graph.add_argument('--heatgrp', help='Specify group to use for heatmap', default='group', required=False)
    parser_graph.add_argument('--heatrts', choices=['whole_unique', 'fiveprime_unique', 'threeprime_unique', 'other_unique', 'total_unique', \
         'wholecounts', 'fiveprime', 'threeprime', 'other', 'total', 'antisense', 'wholeprecounts', 'partialprecounts', 'trailercounts', 'all'], \
            help='Specify readtypes to use for heatmap (default: whole_unique, fiveprime_unique, threeprime_unique, other_unique, total_unique) (optional)', \
                nargs='+', default=['whole_unique', 'fiveprime_unique', 'threeprime_unique', 'other_unique', 'total_unique'], required=False)
    parser_graph.add_argument('--heatcutoff', help='Specify readcount cutoff to use for heatmap', default=80, required=False, type=int)
    parser_graph.add_argument('--heatbound', help='Specify range to use for bounding the heatmap to top and bottom counts', default=25, required=False)
    parser_graph.add_argument('--heatsubplots', help='Specify wether to generate subplots for each comparasion in addition to the sum (default: False)', action='store_true', default=False, required=False)
    # PCA options
    parser_graph.add_argument('--pcamarkers', help='Specify AnnData column to use for PCA markers (default: sample) (optional)', default='sample')
    parser_graph.add_argument('--pcacolors', help='Specify AnnData column to color PCA markers by (default: group) (optional)', default='group')
    parser_graph.add_argument('--pcareadtypes', choices=['whole_unique', 'fiveprime_unique', 'threeprime_unique', 'other_unique', 'total_unique', \
         'wholecounts', 'fiveprime', 'threeprime', 'other', 'total', 'antisense', 'wholeprecounts', 'partialprecounts', 'trailercounts', 'all'], \
            help='Specify read types to use for PCA markers (default: total_unique, total) (optional)', nargs='+', default=['total_unique', 'total'])
    # Correlation options
    parser_graph.add_argument('--corrmethod', choices=['pearson', 'spearman', 'kendall'], help='Specify correlation method (default: pearson) (optional)', default='pearson', required=False)
    parser_graph.add_argument('--corrgroup', help='Specify a grouping variable to generate correlation matrices for (default: sample) (optional)', default='sample', required=False)
    # Volcano options
    parser_graph.add_argument('--volgrp', help='Specify group to use for volcano plot', default='group', required=False)
    parser_graph.add_argument('--volrt', help='Specify readtype to use for volcano plot', default='nreads_total_unique_norm', required=False)
    parser_graph.add_argument('--volcutoff', help='Specify readcount cutoff to use for volcano plot', default=80, required=False)
    # Radar options
    parser_graph.add_argument('--radargrp', help='Specify AnnData column to group by (default: group) (optional)', default='group', required=False)
    # Bar options
    parser_graph.add_argument('--bargrp', help='Specify AnnData column to group by (default: sample) (optional)', default='sample', required=False)
    # Compare options
    parser_graph.add_argument('--comparegrp1', help='Specify AnnData column as main comparative group (default: group) (optional)', default='group', required=False)
    parser_graph.add_argument('--comparegrp2', help='Specify AnnData column to group by (default: group) (optional)', default='group', required=False)

    args = parser.parse_args()

    # Set log file if specified
    sys.stdout = open(args.log, 'w') if args.log else sys.stdout

    if args.mode == 'build':
        args.traxdir = directory_tools.path_cleaner(args.traxdir)
    if args.mode == 'graph':
        args.anndata = directory_tools.path_cleaner(args.anndata)

    # Read database object or create one from trax run if none provided
    if args.mode == 'build':
        # Raise exception if trax directory is empty or doesn't exist
        if not os.path.isdir(args.traxdir):
            raise Exception('Error: trax directory does not exist.')
        # Raise exception if metadata file is empty or doesn't exist
        if args.metadata:
            if not os.path.isfile(args.metadata):
                raise Exception('Error: metadata file does not exist.')
        # Raise exception if observations file is empty or doesn't exist
        if args.observationsfile:
            if not os.path.isfile(args.observationsfile):
                raise Exception('Error: observations file does not exist.')
        print('Building AnnData object...')
        if args.observationslist and args.observationsfile:
            print('Error: Only one of --observationslist or --observationsfile can be used. Defaulting to --observationsfile...')
        if not args.observationslist and not args.observationsfile:
            print('No observations provided. Defaulting to sample names from trax coverage file...')
        if args.observationsfile:
            args.observationslist = []
            with open(args.observationsfile, 'r') as f:
                for line in f:
                    args.observationslist += line.split()
        
        trax2anndata(args.traxdir, args.metadata, args.observationslist, args.output).create()
        print('Done!\n')
    elif args.mode == 'merge':
        # Raise exception if h5ad file is empty or doesn't exist
        if not os.path.isfile(args.anndata1):
            raise Exception('Error: first h5ad file does not exist.')
        if not os.path.isfile(args.anndata2):
            raise Exception('Error: second h5ad file does not exist.')
        print('Merging database objects...\n')
        anndataMerger(args).merge()
        print('Done!\n')
    elif args.mode == 'cluster':
        # Raise exception if h5ad file is empty or doesn't exist
        if not os.path.isfile(args.anndata):
            raise Exception('Error: h5ad file does not exist.')
        print('Clustering data from database object...\n')
        anndataCluster(args).main()
        print('Done!\n')
    elif args.mode == 'graph':
        # Create output directory if it doesn't exist
        directory_tools.builder(args.output)
        # Raise exception if h5ad file is empty or doesn't exist
        if not os.path.isfile(args.anndata):
            raise Exception('Error: h5ad file does not exist.')
        print('Graphing data from database object...\n')
        anndataGrapher(args).graph()
        print('Done!\n')
    else:
        print('Invalid operating mode. Exiting...')
        parser.print_help()