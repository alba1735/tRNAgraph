#!/usr/bin/env python3

import pandas as pd
import anndata as ad
import argparse
import os

import directory_tools

# graphing functions
# import graphs_coverage
import graphs_pca

class trax2anndata():
    '''
    Create h5ad AnnData object from a trax run
    '''
    def __init__(self, traxdir, observations, output):
        traxcoverage = traxdir + '/' + traxdir.split('/')[-1] + '-coverage.txt'
        sizefactors = traxdir + '/' + traxdir.split('/')[-1] + '-SizeFactors.txt'
        trnauniquecounts = traxdir + '/unique/' + traxdir.split('/')[-1] + '-trnauniquecounts.txt'

        self.coverage = pd.read_csv(traxcoverage, sep='\t', header=0)
        self.size_factors = pd.read_csv(sizefactors, sep=" ", header=0).to_dict('index')[0]
        self.unique_counts = pd.read_csv(trnauniquecounts, sep='\t', header=0).to_dict('index')
        self.observations = observations
        self.cov_types = ['coverage', 'readstarts', 'readends', 'uniquecoverage', 'multitrnacoverage',
                          'multianticodoncoverage', 'multiaminocoverage','tRNAreadstotal', 'mismatchedbases',
                          'deletedbases', 'adenines', 'thymines', 'cytosines', 'guanines', 'deletions']
        self.positions = pd.unique(self.coverage['position'])
        self.output = output

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
        
        obs_df,x_df = self._group_sort_(obs_df,x_df)
        adata = self._adata_build_(obs_df,x_df)

        # Add size factors to adata object as raw layer
        adata.layers['raw'] = adata.X*adata.obs['deseq2_sizefactor'].values[:,None]

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
            print('Number of observations does not match number of parameters in trax coverage file by {}.\nTo create a more specific database object, please provide the correct number of observations.\nAdding blank observations to database object...\n'.format(diff_obs_count))
            obs += ['obs_' + str(x) for x in range(diff_obs_count)]

        self.observations = obs[1:]

        obs_df = pd.DataFrame([x.split('_') for x in x_df.index.values], columns=obs, index=x_df.index)
        trna_obs = obs_df['trna'].str.split('-',n=4,expand=True)
        obs_df['trna_type'] = trna_obs[0]
        obs_df['amino'] = trna_obs[1]
        obs_df['iso'] = trna_obs[2]
        obs_df['refseq'] = self.seqs
        obs_df['sizefactor'] = [self.size_factors.get(i) for i in ['_'.join(x.split('_')[1:]) for x in x_df.index.values]]
        obs_df['nreads_raw'] = [self.unique_counts.get(i[0]).get('_'.join(i[1:])) if self.unique_counts.get(i[0]) else 0 for i in [x.split('_') for x in x_df.index.values]] # Some samples may not have any reads for a given tRNA in the unique_counts dictionary might want to double check this
        obs_df['nreads'] = obs_df['nreads_raw']/obs_df['sizefactor']
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
        positions = [i.split('_')[0] for i in x_df.columns.values]
        coverage = [i.split('_')[-1] for i in x_df.columns.values]
        
        # Generate gap list
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
            
        # Add the numer of reads per tRNA as observations-annotation to adata from trna unique counts file
        adata.obs['nreads_raw'] = obs_df['nreads_raw'].values
        adata.obs['nreads'] = obs_df['nreads'].values

        # combine all columns of obs_df after trna, iso, and amino to make a unique sample column
        adata.obs['sample'] = obs_df[[y for y in self.observations if y != 'sample']].agg('_'.join, axis=1).values

        # Add size factor
        adata.obs['deseq2_sizefactor'] = obs_df['sizefactor'].values

        # Add aligned reference sequence based on values in coverage file
        adata.obs['refseq'] = obs_df['refseq'].values
        
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
            self.args.graphtypes = ['coverage', 'heatmap', 'pca', 'volcano', 'radar']

        # if 'coverage' in self.graph_types:
        #     print('Generating coverage graphs...')
        #     output = self.output + '/coverage'
        #     directory_builder(output)
        #     graphs_coverage.grapher(self.adata, output).create()
        #     print('Coverage graphs generated.\n')

        if 'pca' in self.args.graphtypes:
            print('Generating pca plots...')
            output = self.args.output + '/pca'
            directory_tools.builder(output)
            graphs_pca.grapher(self.adata, output, self.args.pcamarkers, self.args.pcacolors).create()
            print('Coverage pca plots.\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='tRNAgraph',
        description='tRNAgraph is a tool for analyzing tRNA-seq data generated from tRAX. It can be used to create a AnnData database object from \
            a trax coverage file, or to analyze an existing database object and generate expanded visulizations. The database object can also be used to \
            perform further analysis beyond the scope of what tRAX can do.',
        allow_abbrev=False
    )
    parser.add_argument('-o', '--output', help='Specify output directory', default='trnagraph')

    subparsers = parser.add_subparsers(
        title='Operating modes',
        description='Choose between building a database object or graphing data from an existing database object',
        dest='mode',
        required=True
    )

    parser_build = subparsers.add_parser("build", help="Build a h5ad AnnData object from a tRAX run")
    parser_build.add_argument('-i', '--traxdir', help='Specify location of trax directory (required)', required=True)
    parser_build.add_argument('-l', '--observationslist', help='Specify the observations of sample names in order (optional)', nargs='*', default=None)
    parser_build.add_argument('-f', '--observationsfile', help='Specify a file containing the observations of sample names in order as tab seperated file (optional)', default=None)

    parser_graph = subparsers.add_parser("graph", help="Graph data from an existing h5ad AnnData object")
    parser_graph.add_argument('-i', '--anndata', help='Specify location of h5ad object (required)', required=True)
    parser_graph.add_argument('-g', '--graphtypes', choices=['all','coverage','heatmap','pca','volcano','radar'], help='Specify graphs to create, if not specified it will default to "all" (optional)', nargs='+', default='all')
    # parser_graph.add_argument('-r', '--raw', help='Specify to use raw data for graphing (optional)', action='store_true') # Implement this later
    parser_graph.add_argument('--pcamarkers', help='Specify AnnData column to use for PCA markers (default: sample) (optional)', default='sample')
    parser_graph.add_argument('--pcacolors', help='Specify AnnData column to color PCA markers by (default: sample) (optional)', default='sample')

    args = parser.parse_args()

    if args.mode == 'build':
        args.traxdir = directory_tools.path_cleaner(args.traxdir)
    if args.mode == 'graph':
        args.anndata = directory_tools.path_cleaner(args.anndata)

    # Create output directory if it doesn't exist
    directory_tools.builder(args.output)

    # Read database object or create one from trax run if none provided
    if args.mode == 'build':
        # Raise exception if trax directory is empty or doesn't exist
        if not os.path.isdir(args.traxdir):
            raise Exception('Error: trax directory does not exist.')
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
        trax2anndata(args.traxdir, args.observationslist, args.output).create()
        print('Done!\n')
    elif args.mode == 'graph':
        # Raise exception if h5ad file is empty or doesn't exist
        if not os.path.isfile(args.anndata):
            raise Exception('Error: h5ad file does not exist.')
        print('Graphing data from database object...')
        anndataGrapher(args).graph()
        print('Done!\n')
    else:
        print('Invalid operating mode. Exiting...')
        parser.print_help()