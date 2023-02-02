#!/usr/bin/env python3

import pandas as pd
import anndata as ad
import argparse
import os

import directory_tools

# graphing functions
import coverage_tools
import pca_tools
import correlation_tools

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
        self.read_types = pd.unique(normalized_read_counts.index.str.split('_').str[1])
        self.normalized_read_counts = normalized_read_counts.to_dict('index')
        # For adding metadata to adata object
        metadata_type = '\t' if metadata.endswith('.tsv') else ',' if metadata.endswith('.csv') else ' '
        self.metadata = pd.read_csv(metadata, sep=metadata_type, header=None, index_col=None)
        # For adding categories to adata object
        self.observations = observations
        if self.observations:
            # Make sure that observations are not going to be generated automatically
            auto_obs = ['trna', 'iso', 'amino', 'sample', 'group', 'deseq2_sizefactor', 'refseq']
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
        # Save adata object
        adata.write('{}.h5ad'.format(self.output))
        print('Writing h5ad database object to {}.h5ad'.format(self.output))
        
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
        obs_df['nreads_unique_raw'] = [self.unique_counts.get(i[0]).get('_'.join(i[1:])) if self.unique_counts.get(i[0]) else 0 for i in [x.split('_') for x in x_df.index.values]] # Some samples may not have any reads for a given tRNA in the unique_counts dictionary might want to double check this
        obs_df['nreads_unique_norm'] = obs_df['nreads_unique_raw']/obs_df['sizefactor']
        for i in self.read_types:
            obs_df['nreads_' + i + '_norm'] = [self.normalized_read_counts.get(j[0] + '_' + i).get('_'.join(j[1:])) if self.normalized_read_counts.get(j[0] + '_' + i) else 0 for j in [x.split('_') for x in x_df.index.values]]
        obs_df['uniquefeat'] = obs_df.index.values
        return obs_df

    # def _group_sort_(self, obs_df, x_df):
    #     '''
    #     Sort obs and x dataframes by sample group to ensure that all dfs are aligned
    #     '''
    #     obs_dict = {k:'first' for k in obs_df.columns.values if k != 'sample_group'}
    #     x_dict = {k:'mean' for k in x_df.columns.values}

    #     combo_dict = {**obs_dict, **x_dict}
        
    #     combo_df = pd.concat([obs_df,x_df], axis=1, sort=False).reset_index(drop=True)
    #     combo_df = combo_df.groupby('group', as_index=False).agg(combo_dict)
        
    #     obs_df = combo_df.iloc[:,:len(obs_df.columns)]
    #     x_df = combo_df.iloc[:,len(obs_df.columns):]
        
    #     obs_df.index = obs_df['group'].values
    #     x_df.index = obs_df['group'].values

    #     # Check that the index of the obs and x dataframes are the same
    #     if not obs_df.index.equals(x_df.index):
    #         raise ValueError('The index of the obs and x dataframes are not the same. This means somthing went wrong in the sorting process.')

    #     return obs_df, x_df

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
        # Add metadata dataframe
        adata.obs['trna'] = obs_df['trna'].values
        adata.obs['iso'] = obs_df['iso'].values
        adata.obs['amino'] = obs_df['amino'].values
        # Add custom dataframe obs
        for ob in self.observations:
            adata.obs[ob] = obs_df[ob].values 
        # Add the numer of reads per tRNA as observations-annotation to adata from trna unique counts file
        adata.obs['nreads_unique_raw'] = obs_df['nreads_unique_raw'].values
        adata.obs['nreads_unique_norm'] = obs_df['nreads_unique_norm'].values
        for i in self.read_types:
            adata.obs['nreads_' + i + '_norm'] = obs_df['nreads_' + i + '_norm'].values
        # Add sample and group metadata
        adata.obs['sample'] = obs_df['sample'].values
        adata.obs['group'] = obs_df['group'].values
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
            self.args.graphtypes = ['coverage', 'heatmap', 'pca', 'correlation', 'volcano', 'radar']

        if 'coverage' in self.args.graphtypes:
            print('Generating coverage plots...')
            output = self.args.output + '/coverage'
            directory_tools.builder(output)
            coverage_tools.visualizer(self.adata.copy(), self.args.coveragegrp, self.args.coverageobs, self.args.coveragetype, self.args.coveragegap, self.args.coveragefill, output).generate_all()

        if 'pca' in self.args.graphtypes:
            print('Generating pca plots...')
            output = self.args.output + '/pca'
            directory_tools.builder(output)
            pca_tools.visualizer(self.adata.copy(), output, self.args.pcamarkers, self.args.pcacolors)

        if 'correlation' in self.args.graphtypes:
            print('Generating correlation plots...')
            output = self.args.output + '/correlation'
            directory_tools.builder(output)
            correlation_tools.visualizer(self.adata.copy(), output, self.args.corrmethod, self.args.corrgroup)

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
    parser_build.add_argument('-m', '--metadata', help='Specify a metadata file to create annotations, you can also use the sample file used to generate tRAX DB (required)', required=True)
    parser_build.add_argument('-l', '--observationslist', help='Specify the observations of sample names in order (optional)', nargs='*', default=None)
    parser_build.add_argument('-f', '--observationsfile', help='Specify a file containing the observations of sample names in order as tab seperated file (optional)', default=None)

    parser_graph = subparsers.add_parser("graph", help="Graph data from an existing h5ad AnnData object")
    parser_graph.add_argument('-i', '--anndata', help='Specify location of h5ad object (required)', required=True)
    parser_graph.add_argument('-g', '--graphtypes', choices=['all','coverage','heatmap','pca', 'correlation', 'volcano','radar'], help='Specify graphs to create, if not specified it will default to "all" (optional)', nargs='+', default='all')
    # Coverage options
    parser_graph.add_argument('--coveragegrp', help='Specify a grouping variable to generate coverage plots for (default: sample) (optional)', default='group')
    parser_graph.add_argument('--coverageobs', help='Specify a observation subsetting for coverage plots (optional)', nargs='+', default=None)
    parser_graph.add_argument('--coveragetype', help='Specify a coverage type for coverage plots corresponding to trax coverage file outputs (default: uniquecoverage) (optional)', choices=['coverage', 'readstarts', \
         'readends', 'uniquecoverage', 'multitrnacoverage', 'multianticodoncoverage', 'multiaminocoverage','tRNAreadstotal', 'mismatchedbases', 'deletedbases', 'adenines', 'thymines', 'cytosines', 'guanines', 'deletions'], default='uniquecoverage')
    parser_graph.add_argument('--coveragegap', help='Specify wether to include gaps in coverage plots (default: False) (optional)', default=False)
    parser_graph.add_argument('--coveragefill', choices=['fill', 'ci', 'none'], help='Specify wether to fill area under coverage plots or use confidence intervals (default: ci) (optional)', default='ci')
    # PCA options
    parser_graph.add_argument('--pcamarkers', help='Specify AnnData column to use for PCA markers (default: sample) (optional)', default='sample')
    parser_graph.add_argument('--pcacolors', help='Specify AnnData column to color PCA markers by (default: sample) (optional)', default='sample')
    # Correlation options
    parser_graph.add_argument('--corrmethod', choices=['pearson', 'spearman', 'kendall'], help='Specify correlation method (default: pearson) (optional)', default='pearson', required=False)
    parser_graph.add_argument('--corrgroup', help='Specify a grouping variable to generate correlation matrices for (optional)', required=False)

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