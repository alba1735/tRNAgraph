#!/usr/bin/env python3

import numpy as np
import pandas as pd
import anndata as ad
import os
import argparse

class trax2anndata():
    '''
    Create h5ad AnnData object from a trax run
    '''
    def __init__(self, traxcoverage, observations, output, database):
        self.coverage = pd.read_csv(traxcoverage, sep='\t', header=0)
        self.size_factors = pd.read_csv(traxcoverage.split('-')[0] + '-SizeFactors.txt',sep=" ",header=0).to_dict('index')[0]
        self.observations = observations
        self.cov_types = ['coverage', 'readstarts', 'readends', 'uniquecoverage', 'multitrnacoverage',
                          'multianticodoncoverage', 'multiaminocoverage','tRNAreadstotal', 'mismatchedbases',
                          'deletedbases', 'adenines', 'thymines', 'cytosines', 'guanines', 'deletions']
        self.positions = pd.unique(self.coverage['position'])
        self.output = output
        self.database = database

    def create(self):
        '''
        Create h5ad database object
        '''
        self._ref_seq_build_()

        x_dfs = []
        for cov_type in self.cov_types:
            x_df = self._x_build_(cov_type)
            x_dfs.append(x_df)

        obs_df = self._obs_build_(x_df)
        x_df = pd.concat(x_dfs,axis=1,sort=False)

        clist = [[p + '_' + cov for p in self.positions] for cov in self.cov_types]
        clist = [a for l in clist for a in l]
        x_df.columns = clist

        # if raw==True:
        #     x_raw_df = x_df.T*obs_df['sizefactor']
        #     x_raw_df = x_raw_df.T
        #     x_df = x_raw_df
        
        obs_df,x_df = self._group_sort_(obs_df,x_df)
        adata = self._adata_build_(obs_df,x_df)
        self.save(adata)

    def save(self, adata):
        '''
        Save h5ad database object
        '''
        print('Writing h5ad database object to {}.h5ad'.format(self.output))
        adata.write('{}.h5ad'.format(self.output))

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
        x_df = self.coverage
        x_df['uniquefeat'] = x_df['Feature'] + '_' + x_df['Sample']
        x_df = x_df.pivot(index='uniquefeat', values=[cov_type], columns='position')
        cols = x_df.columns.get_level_values(1).values
        x_df = x_df.T.reset_index(drop=True).T
        x_df.columns = cols
        
        x_df = x_df.reindex(columns=self.positions)
        
        return x_df
    
    def _obs_build_(self, x_df):
        '''
        Build obs dataframe from trax coverage file derived x dataframe
        '''
        obs = ['trna']
        if self.observations:
            obs += self.observations
        if len(obs) != len([x.split('_') for x in x_df.index.values][0]):
            diff_obs_count = abs(len(obs)-len([x.split('_') for x in x_df.index.values][0]))
            print('Number of observations does not match number of parameters in trax coverage file by {}.\nTo create a more specific database object, please provide proper the correct number of observations.\nAdding blank observations to database object.\n'.format(diff_obs_count))
            obs += ['obs_' + str(x) for x in range(diff_obs_count)]

        self.observations = obs[1:]

        obs_df = pd.DataFrame([x.split('_') for x in x_df.index.values], columns=obs, index=x_df.index)
        trna_obs = obs_df['trna'].str.split('-',n=4,expand=True)
        obs_df['trna_type'] = trna_obs[0]
        obs_df['amino'] = trna_obs[1]
        obs_df['iso'] = trna_obs[2]
        obs_df['refseq'] = self.seqs
        obs_df['sizefactor'] = [self.size_factors.get(i) for i in ['_'.join(x.split('_')[1:]) for x in x_df.index.values]]
        obs_df['sample_group'] = obs_df['trna'] + ('_' + obs_df[[y for y in self.observations if y != 'sample']]).sum(axis=1).str.strip() 

        return obs_df

    def _group_sort_(self, obs_df, x_df):
        '''
        Sort obs and x dataframes by sample group
        '''
        obs_dict = {k:'first' for k in obs_df.columns.values if k != 'sample_group'}
        x_dict = {k:'mean' for k in x_df.columns.values}
        combo_dict = {**obs_dict, **x_dict}
        
        combo_df = pd.concat([obs_df,x_df], axis=1, sort=False).reset_index(drop=True)
        combo_df = combo_df.groupby('sample_group', as_index=False).agg(combo_dict)
        
        obs_df = combo_df.iloc[:,:len(obs_df.columns)]
        x_df = combo_df.iloc[:,len(obs_df.columns):]
        
        obs_df.index = obs_df.sample_group.values
        x_df.index = obs_df.sample_group.values
        
        return obs_df, x_df

    def _adata_build_(self, obs_df, x_df):
        '''
        Build AnnData object from obs and x dataframes
        '''
        print('Building AnnData object...')
        positions = [i.split('_')[0] for i in x_df.columns.values]
        coverage = [i.split('_')[-1] for i in x_df.columns.values]
        
        gap_list = [True if len(i.split('.')) > 1 else False for i in positions]
        gap_list = [True if positions[i][-1] == 'a' else gap_list[i] for i in range(len(gap_list))]
        gap_list = [True if positions[i][-1] == 'b' else gap_list[i] for i in range(len(gap_list))]
        gap_list = [True if positions[i][0] == 'e' else gap_list[i] for i in range(len(gap_list))]
        
        adata = ad.AnnData(x_df, dtype='float64')
        adata.var['gap'] = gap_list
        adata.var['positions'] = positions
        adata.var['coverage'] = coverage

        # Add metadata dataframe
        adata.obs['trna'] = obs_df['trna'].values
        adata.obs['iso'] = obs_df['iso'].values
        adata.obs['amino'] = obs_df['amino'].values
        
        # Add custom dataframe obs
        for ob in [y for y in self.observations if y != 'sample']:
            adata.obs[ob] = obs_df[ob].values
            
        # Add the total counts per tRNA as observations-annotation to adata
        adata.obs['n_counts'] = adata.X.sum(axis=1)
        
        # Add the total reads per tRNA as observations-annotation to adata
        ## MAKE SURE THESE AREN'T GETTING HIGHEST READS FROM GAPS?
        mask_out = np.isin(adata.var.coverage, ['uniquecoverage'])
        adata_u = adata[:,mask_out]
        adata.obs['n_reads'] = adata_u.X.max(axis=1)
        
        # Add reference sequence
        adata.obs['refseq'] = obs_df['refseq'].values

        #if seq_build == True:
        #    adata = self._seq_builder(adata)
        #adata = self._mods_builder(adata) 
        #adata = self._epi_builder(adata)
        
        return adata

class anndataGrapher:
    def __init__(self, db, graph_types, graph_args, output_dir):
        self.adata = ad.read_h5ad(db)
        self.output_dir = output_dir
        self.graph_types = graph_types
        self.graph_args = graph_args

    def graph(self, graph_type, graph_args):
        pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='tRNAgraph',
        description='tRNAgraph is a tool for analyzing tRNA-seq data generated from tRAX. It can be used to create a AnnData database object from \
            a trax coverage file, or to analyze an existing database object and generate expanded visulizations. The database object can also be used to \
            perform further analysis beyond the scope of what tRAX can do.'
    )
    parser.add_argument('-o', '--output', help='Specify output directory', default='tRNAgraph')

    subparsers = parser.add_subparsers(
        title='Operating modes',
        description='Choose between building a database object or graphing data from an existing database object',
        dest='mode',
        required=True
    )

    parser_build = subparsers.add_parser("build", help="Build a h5ad AnnData object from a tRAX run")
    parser_build.add_argument('-tc', '--traxcoverage', help='Specify location of trax coverage file', required=True)
    parser_build.add_argument('-obs', '--observations', help='Specify the observations of sample names in order', nargs='*', default=None)
    parser_build.add_argument('-obsf', '--observations_file', help='Specify a file containing the observations of sample names in order as tab seperated file', default=None)

    parser_graph = subparsers.add_parser("graph", help="Graph data from an existing h5ad AnnData object")
    parser_graph.add_argument('-db', '--database', help='Specify location of h5ad object', required=True)
    parser_graph.add_argument('-gt', '--graph_types', help='Specify graphs to create', nargs='+', required=True)
    parser_graph.add_argument('-ga', '--graph_args', help='Specify arguments for graphs', nargs='+', required=False)

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    if args.output:
        if not os.path.exists(args.output):
            print('Creating output directory: {}'.format(args.output))
            os.mkdir(args.output)
        else:
            print('Output directory already exists: {}'.format(args.output))

    # Read database object or create one from trax run if none provided
    if args.mode == 'build':
        print('Building database object...\n')
        if args.observations and args.observations_file:
            print('Error: Only one of --observations or --observations_file can be used. Defaulting to observations_file...\n')
            trax2anndata(args.traxcoverage, args.observations_file, args.output).create()
        elif args.observations_file:
            trax2anndata(args.traxcoverage, args.observations_file, args.output).create()
        else:
            trax2anndata(args.traxcoverage, args.observations, args.output).create()
    elif args.mode == 'graph':
        print('Graphing data from database object...\n')
        anndataGrapher(args.database, args.graph_types, args.graph_args, args.output).graph()
    else:
        print('Invalid operating mode. Exiting...\n')
        parser.print_help()