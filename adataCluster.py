import anndata as ad
import os
import pandas as pd
import numpy as np
from . import toolsTG
from .lazy_imports import umap, hdbscan

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
        self.variance_threshold = args.variancethreshold
        self.stats_metrics_umap = args.umapstatsmetrics
        self.stats_metrics_hdbscan = args.hdbstatsmetrics
        self.mindist = args.mindist

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
        from sklearn.preprocessing import RobustScaler
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
        # Convert view to copy to avoid ImplicitModificationWarning when modifying X
        if adata.is_view:
            adata = adata.copy()
        # Normalize the each read at each position by the total coverage - This would collapse the variation between positons so should only be used in certain cases depending on the var slice used
        # adata.X = Normalizer().fit_transform(adata.X)
        # Scale the data - Seems to perform well with the robust scaler compared to the standard scaler
        adata.X = RobustScaler().fit_transform(adata.X)

        return adata
    
    def adataCluster(self, adata, neighbors_plot, neighbors_cluster, min_samples, min_cluster_size, n_components):
        from sklearn.feature_selection import VarianceThreshold
        # Remove low variance features
        sel = VarianceThreshold(threshold=(self.variance_threshold))
        # Apply a standardscaler to the data and reduce dimensions
        standard_embedding = umap.UMAP(random_state=self.randomstate, n_neighbors=neighbors_plot, min_dist=self.mindist, metric=self.stats_metrics_umap).fit_transform(sel.fit_transform(adata.X))
        cluster_embedding = umap.UMAP(random_state=self.randomstate, n_neighbors=neighbors_cluster, min_dist=0.0, n_components=n_components, metric=self.stats_metrics_umap).fit_transform(sel.fit_transform(adata.X))
        # Perform clustering with HDBSCAN
        hdbscan_results = hdbscan.HDBSCAN(min_samples=min_samples, min_cluster_size=min_cluster_size, metric=self.stats_metrics_hdbscan).fit_predict(cluster_embedding)
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
