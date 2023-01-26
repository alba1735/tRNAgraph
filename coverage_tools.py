#!/usr/bin/env python3

import numpy as np
import pandas as pd
import anndata as ad
import argparse

import directory_tools

import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import seaborn as sns

class visualizer():
    '''
    Generate coverage plots for each sample in an AnnData object.
    '''
    def __init__(self, adata, coverage_grp, coverage_obs, coverage_type, coverage_gap, coverage_pal, output):
        self.coverage_grp = coverage_grp
        self.coverage_obs = coverage_obs
        self.coverage_type = coverage_type
        self.coverage_gap = coverage_gap
        self.coverage_pal = coverage_pal
        self.adata = self.clean_adata(adata)
        self.output = output

    def clean_adata(self, adata):
        '''
        Clean AnnData object for plotting.
        '''
        # Remove end positions
        #adata_test = adata_test[:,np.invert(np.isin(adata_test.var.positions, ['-1','76']))] 
        ## Need to determine if we should also remove 76/75/74 etc.
        adata = adata[:,np.invert(np.isin(adata.var.positions, ['-1']))]
        # Subset AnnData observations if specified
        if self.coverage_obs:
            for k,v in self.coverage_obs.items():
                adata = adata[np.isin(adata.obs[k].values, v),:]
        # Subset by coverage type from AnnData variables
        adata = adata[:,np.isin(adata.var.coverage, [self.coverage_type])]
        # Subset gaps from AnnData variables
        adata = adata[:,np.isin(adata.var.gap, self.coverage_gap)]
        
        return adata

    def generate_all(self):
        '''
        Generate coverage plots for all tRNAs.
        '''
        for trna in sorted(self.adata.obs.trna.unique()):
            self.generate_single(trna)

    def generate_single(self, trna):
        '''
        Generate coverage plots for a single tRNA.
        '''
        # Subset AnnData object to a single tRNA
        adata_t = self.adata[self.adata.obs.trna == trna].copy()
        max_reads = adata_t.obs['nreads_unique_norm'].values.max()
        df = pd.DataFrame(adata_t.X.T, columns=adata_t.obs[self.coverage_grp].values)
    
        print(df)
        exit()




if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='coverage_tools.py',
        description='Generate tRNA coverage plots'
    )

    parser.add_argument('-i', '--anndata', help='Specify AnnData input', required=True)
    parser.add_argument('-o', '--output', help='Specify output directory', default='coverage', required=False)
    parser.add_argument('--coveragegrp', help='Specify a grouping variable to generate coverage plots for (default: sample) (optional)', default='sample')
    parser.add_argument('--coverageobs', help='Specify a observation subsetting for coverage plots (optional)', nargs='+', default=None)
    parser.add_argument('--coveragetype', help='Specify a coverage type for coverage plots corresponding to trax coverage file outputs (default: uniquecoverage) (optional)', \
        choices=['coverage', 'readstarts', 'readends', 'uniquecoverage', 'multitrnacoverage', 'multianticodoncoverage', 'multiaminocoverage','tRNAreadstotal', 'mismatchedbases', \
            'deletedbases', 'adenines', 'thymines', 'cytosines', 'guanines', 'deletions'], default='uniquecoverage')
    parser.add_argument('--coveragegap', help='Specify wether to include gaps in coverage plots (default: False) (optional)', default=False)
    parser.add_argument('--coveragepal', help='Specify a palette for coverage plots (optional)', nargs='+', default=None)

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    directory_tools.builder(args.output)

    adata = ad.read_h5ad(args.anndata)

    visualizer(adata, args.coveragegrp, args.coverageobs, args.coveragetype, args.coveragegap, args.coveragepal, args.output).generate_all()

class coverage_plots():
    def __init__(self, adata_path, save_path='coverage_plots', key_feat='index_names', raw=False, subset_o={}, subset_v={}, label_dict={}, pal_dict=None):
        sns.set(font_scale=1)
        self.adata_path = adata_path
        self.allvalue = False
        self.raw = raw
        self.overlaptemp = False
        self.subset_o = subset_o
        self.subset_v = subset_v
        self.key_feat = key_feat
        self.label_dict = label_dict
        self.save_path = save_path
        self.adata = self.adata_to_load(adata_path)
        self.max_reads = 0
        
        os.makedirs(self.save_path, exist_ok=True)
        
        if pal_dict:
            self.pal = pal_dict
        else:
            lv = sorted(self.adata.obs[self.key_feat].unique())
            self.pal = {lv[i]:sns.hls_palette(len(lv))[i] for i in range(len(lv))}
        
    def louvain(self, overlaptemp=False):
        self.overlaptemp = overlaptemp
        self.adata = self.adata_to_load(self.adata_path)
        
        with sns.axes_style(style="white", rc={'axes.grid':False}):
            t = self.adata.obs['lclass'].unique()
            #t = sorted([int(i) for i in t])
            t = sorted([i for i in t])
            t = ['{}'.format(i) for i in t]
            pal = sns.color_palette("hls", len(self.adata.obs['lclass'].unique()))
            pal_dict = dict(zip(t,pal))

            for i in t:
                adata_t = self.adata[np.isin(self.adata.obs.lclass.values, [i]),:]
                if self.overlaptemp == True:
                    title ='Louvain Overlap {} {}'.format(i,' '.join(self.subset_o.values()))
                    self.overlap(adata_t,title,pal_dict[i])
                else:
                    title = '{}/louvain_{}_{}_by_{}.pdf'.format(self.save_path,i,'_'.join(self.subset_o.values()),self.key_feat)
                    self.multi_cplot(adata_t,title)
                    
    def isotype(self, overlaptemp=False):
        self.overlaptemp = overlaptemp
        self.adata = self.adata_to_load(self.adata_path)
        
        with sns.axes_style(style="white", rc={'axes.grid':False}):
            t = self.adata.obs['iso'].unique()
            t = sorted([i for i in t])
            t = ['{}'.format(i) for i in t]
            pal = sns.color_palette("hls", len(self.adata.obs['iso'].unique()))
            pal_dict = dict(zip(t,pal))

            for i in t:
                adata_t = self.adata[np.isin(self.adata.obs.iso.values, [i]),:]
                if self.overlaptemp == True:
                    title ='Iso Overlap {} {}'.format(i,' '.join(self.subset_o.values()))
                    self.overlap(adata_t,title,pal_dict[i])
                else:
                    title = '{}/iso_{}_{}_by_{}.pdf'.format(self.save_path,i,'_'.join(self.subset_o.values()),self.key_feat)
                    self.multi_cplot(adata_t,title)
                    
    def multi(self):
        with sns.axes_style(style="white", rc={'axes.grid':False}):
            self.multi_cplot(self.adata)

    def allvalues(self):
        self.allvalue = True
        for trna in sorted(self.adata.obs.trna.unique()):
            self.single(trna)
                    
    def single(self, trna, alt_lines=None):
        with sns.axes_style(style="white", rc={'axes.grid':False}):
            title_str = ' '.join(self.subset_o.values())
            title ='{} {}'.format(trna,title_str)
            df = self.adata_to_df(self.adata,trna)
            self.cplot(df,title,alt_lines)

    def multi_cplot(self, adata, title=None):
        if not title:
            title = '{}/multi_{}_by_{}.pdf'.format(self.save_path,'_'.join(self.subset_o.values()),self.key_feat)
        with PdfPages(title) as pdf:
            cutoff = 20
            adata_sub = adata[adata.obs.n_reads >= cutoff, :]
            self.multi_cplot_plot(adata_sub,pdf)
            adata_sub = adata[adata.obs.n_reads < cutoff, :]
            self.multi_cplot_plot(adata_sub,pdf)
            
    def overlap(self, adata, title=None, pal=False):
        trna_ls = sorted(adata.obs.trna.unique().to_list())
 
        jdf = []
        for trna in trna_ls:
            jdf.append(self.adata_to_df(adata,trna))
            
        df = pd.concat(jdf, axis=1)
        
        self.cplot(df,title,overlap=pal)    
            
    def multi_cplot_plot(self, adata, pdf, n = 20):
            trna_ls = sorted(adata.obs.trna.unique().to_list())
            trna_ls = [trna_ls[i * n:(i + 1) * n] for i in range((len(trna_ls) + n - 1) // n )]  
        
            for trna_l in trna_ls:
                fig = plt.figure(figsize=(18,18))
                for i in range(len(trna_l)):
                    trna = trna_l[i]
                    title_str = ' '.join(self.subset_o.values())
                    title ='{} {}'.format(trna,title_str)
                    df = self.adata_to_df(adata,trna)
                    ax,fig = self.cplot(df,title,subplot_pos=i+1,subplot_fig=fig)
                pdf.savefig(ax.get_figure())
                plt.close()
            
    def adata_to_load(self, adata_path):
        adata_test = sc.read(adata_path)
        
        if self.raw == True:
            adata_raw = adata_test.raw[:,np.isin(adata_test.raw.var.gap.values, [False])]
            adata_raw = adata_raw[:,np.isin(adata_raw.var.coverage.values, [self.subset_v['coverage']])]
            
            adata_test = adata_test[:,np.isin(adata_test.var.gap.values, [False])]
            adata_test = adata_test[:,np.isin(adata_test.var.coverage.values, [self.subset_v['coverage']])]

            adata_test.X = adata_raw.X
        
        # Add default name parameters
        adata_test.obs['index_names'] = adata_test.obs.index.values
        
        # Remove end positions
        #adata_test = adata_test[:,np.invert(np.isin(adata_test.var.positions, ['-1','76']))]
        adata_test = adata_test[:,np.invert(np.isin(adata_test.var.positions, ['-1']))]

        if self.overlaptemp != True:
            # Select AlkB pos only
            #adata_test = adata_test[np.isin(adata_test.obs.alkb.values, ['pAlkB']),:]

            for k,v in self.subset_o.items():
                adata_test = adata_test[np.isin(adata_test.obs[k].values, v),:]

            for k,v in self.subset_v.items():
                adata_test = adata_test[:,np.isin(adata_test.var[k].values, v)]

            #adata_test = adata_test[:,np.isin(adata_test.var.gap.values, [False])]
        
        return adata_test
    
    def adata_to_df(self, adata, trna):
        # Select specific tRNA only
        mask_out = np.isin(adata.obs.trna.values, [trna])
        adata = adata[mask_out,:]
        
        if adata.obs.n_reads.values.max() > self.max_reads:
            self.max_reads = adata.obs.n_reads.values.max()
        
        # Convert to DF
        df = pd.DataFrame(adata.X.T, columns=adata.obs[self.key_feat].values)
        
        return df

    def cplot(self, df, title, alt_lines=None, subplot_pos=None, subplot_fig=None, overlap=False):
        if subplot_fig:
            fig = subplot_fig
            ax = fig.add_subplot(5,4,subplot_pos)
        else:
            if overlap or self.subset_v['coverage'] != 'uniquecoverage':
                fig, ax = plt.subplots(figsize=(5.5,6)) 
            else:
                fig, ax = plt.subplots(figsize=(6,5.5))
        
        #try:
        y = np.array([[0],[self.max_reads]])
        hcut = 0.1
        if overlap == False:
            if self.subset_v['coverage'] == 'uniquecoverage':
                sns.lineplot(data=df, linewidth=2, dashes=False, ax=ax)
            else:
                x = df.values
                
                min_max_scaler = preprocessing.MinMaxScaler()
                y_scaled = min_max_scaler.fit_transform(y)
                x_scaled = min_max_scaler.transform(x)
                
                # If A,G,T,C
                #x_scaled = min_max_scaler.fit_transform(x)
                
                s_df = pd.DataFrame(x_scaled,columns=df.columns.values)
                
                s_df = s_df/s_df.sum()
                
                alt_lines=[0]
                for i in df:
                    ax.bar(list(df.index.values+0.5), s_df[i].values, 1, alpha=0.5, color=self.pal[i], linewidth=0, edgecolor='black', clip_on=False)
                    alt_lines += list(s_df[s_df[i] >= hcut].index.values)
                
                if len(list(s_df[s_df[i] >= hcut].index.values)) == 0:
                    alt_lines = None
        else:
            ci_val = 0.95
            x = df.values
            min_max_scaler = preprocessing.MinMaxScaler()
            
            if self.subset_v['coverage'] == 'uniquecoverage':
                # Deviation from mean at position method
                df = df.mean(axis=1)/df.mean(axis=1).max()
                ci = st.norm.interval(alpha=ci_val, loc=df, scale=df.sem())
                plt.plot(df, color=overlap, linewidth=2.5, alpha=.5,)
                
                # Min-Max Normalization Method
                ############################
                #x_scaled = min_max_scaler.fit_transform(x)
                #df = pd.DataFrame(x_scaled)
                #ci = st.norm.interval(alpha=ci_val, loc=np.mean(df.T), scale=st.sem(df.T))
                #plt.plot(df.mean(axis=1), color=overlap, linewidth=2.5, alpha=.5,)
                ############################

                #For all lines
                #plt.plot(df, color='black', alpha=0.10)

                # For CI
                ax.fill_between(df.index, (ci[0]), (ci[1]), color=overlap, alpha=.15)
                plt.text(10, .975, 'ci={}'.format(ci_val), fontsize=10, horizontalalignment='right', verticalalignment='center')
            else:
                x = x.mean(axis=1)
                x = [[i] for i in x]
                
                y_scaled = min_max_scaler.fit_transform(y)
                x_scaled = min_max_scaler.transform(x)
                
                x_scaled = min_max_scaler.fit_transform(x)
                
                s_df = pd.DataFrame(x_scaled)
                
                #Think this is important to scale between 0 and 1
                s_df = s_df/s_df.sum()
                
                alt_lines=[]
                for i in s_df:
                    ax.bar(list(df.index.values+0.5), s_df[i].values, 1, alpha=0.5, color=self.pal[i], linewidth=0, edgecolor='black', clip_on=False)
                    alt_lines += list(s_df[s_df[i] >= hcut].index.values)
                
                if len(list(s_df[s_df[i] >= hcut].index.values)) == 0:
                    alt_lines = None
            
        self.cplot_params(df,fig,ax,title,overlap,alt_lines)

        if subplot_fig:
            return ax,subplot_fig
        else:
            title_save = title.replace(' ','_').lower() + '_by_{}'.format(self.key_feat)
            plt.savefig('{}/{}.pdf'.format(self.save_path,title_save), bbox_inches = 'tight')
            if self.allvalue == True:
                plt.close()
        #except:
        #    print('Error with:\n',title,df)
        
    def cplot_params(self, df, fig, ax, title, overlap, alt_lines=None):
        if self.key_feat == 'index_names':
            self.pal = {df.columns[i]:sns.hls_palette(len(df.columns))[i] for i in range(len(df.columns))}
        
        ax.set_xlabel("Positions on tRNA")
        if overlap == False:
            ax.set_ylabel("Normalized Readcounts")
        else:
            ax.set_ylabel("Normalized Fragment Coverage")

        bottom, top = ax.get_ylim()
        ax.set_ylim(0, top)
        ax.set_xlim(0, 73)
        #ax.set_xlim(0, 75)
        #ax.set_xlim(0, 76)
        
        if overlap or self.subset_v['coverage'] != 'uniquecoverage':
            ax.set_ylim(0, 1.005)
            ax.set_yticks([0,0.25,0.5,0.75,1])
            ax.set_yticklabels(['0%','25%','50%','75%','100%'])
            bottom, top = ax.get_ylim()

        ax.fill_between([14,21],[top,top], color='#cacaca', alpha=0.35)
        ax.fill_between([32,38],[top,top], color='#cacaca', alpha=0.35)
        ax.fill_between([54,60],[top,top], color='#cacaca', alpha=0.35)
        
        ax.fill_between([10,25],[top,top], color='#cacaca', alpha=0.35)
        ax.fill_between([27,43],[top,top], color='#cacaca', alpha=0.35)
        ax.fill_between([49,65],[top,top], color='#cacaca', alpha=0.35)
        
        if overlap == False:
            if self.subset_v['coverage'] == 'uniquecoverage':
                ax = sns.lineplot(data=df, linewidth=2, palette=self.pal, dashes=False)

                for i in range(0,len(df.columns))[::-1]:
                    line = ax.lines[i]
                    x = line.get_xydata()[:,0]
                    y = line.get_xydata()[:,1]
                    ax.fill_between(x,y, color=self.pal.get(df.columns[i]), alpha=0.25)
            else:
                pass
                #ax.set_xlim(0, 76)
                #for i in df:
                #    ax.bar(list(df.index.values), df[i].values, 1, alpha=0.25, linewidth=0, edgecolor='black', clip_on=False)

        if alt_lines:
            for l in alt_lines:
                plt.plot([l,l],[0,top],linewidth=1,ls='--',color='black')
                #ax.fill_between([l,l+1],[top,top], color='white', alpha=0.5)
                
                #plt.plot([l,l],[0,top],linewidth=1,ls='--',color='black')
                #plt.plot([l+1,l+1],[0,top],linewidth=1,ls='--',color='black')
                
                #plt.plot([0,76],[0.1,0.1],linewidth=1,ls='--',color='black')
                #ax.annotate('10% Threshold', (1, 0.115), fontsize=6, horizontalalignment='left', verticalalignment='center')
                
            #plt.plot([37,37],[0,top],linewidth=1,ls='--',color='black')  
            #plt.plot([38,38],[0,top],linewidth=1,ls='--',color='black')  
            #plt.plot([58,58],[0,top],linewidth=1,ls='--',color='black')
            #plt.plot([59,59],[0,top],linewidth=1,ls='--',color='black')
                
            #plt.plot([0,76],[0.1,0.1],linewidth=1,ls='--',color='black')
            #ax.annotate('10% Threshold', (1, 0.115), fontsize=6, horizontalalignment='left', verticalalignment='center')   
                
            xtl = sorted(alt_lines + [18.01,35.01,57.01])  
            xtll = [i if i != 18.01 else '\nD-Arm' for i in xtl]
            xtll = [i if i != 35.01 else '\nA-Arm' for i in xtll]
            xtll = [i if i != 57.01 else '\nT-Arm' for i in xtll]
            ax.set_xticks(xtl)
            ax.set_xticklabels(xtll)
        else:
            plt.plot([37,37],[0,top],linewidth=1,ls='--',color='black')  
            #plt.plot([38,38],[0,top],linewidth=1,ls='--',color='black')  
            plt.plot([58,58],[0,top],linewidth=1,ls='--',color='black')
            #plt.plot([59,59],[0,top],linewidth=1,ls='--',color='black')
            
            #plt.plot([0,76],[0.1,0.1],linewidth=1,ls='--',color='black')
            #ax.annotate('10% Threshold', (1, 0.115), fontsize=6, horizontalalignment='left', verticalalignment='center')
            
            ax.set_xticks([18.01,35.01,37,57.01,58])
            ax.set_xticklabels(['\nD-Arm','\nA-Arm','37','\nT-Arm','58'])
                
        ax.set_title(title, fontsize='x-large')

        if overlap == False and self.subset_v['coverage'] == 'uniquecoverage':
            leg_handles, leg_labels = ax.get_legend_handles_labels()
            labels = [self.label_dict.get(label,label) for label in leg_labels][:len(leg_handles)//2]
            if self.key_feat=='index_names':
                labels = ['_'.join(label.split('_')[1:]) for label in labels]
            leg = ax.legend(handles=leg_handles[::-1], labels=labels[::-1], loc='center left', bbox_to_anchor=(1.025, 0.5),
                            ncol=1, frameon=False, fontsize='medium', title='{}'.format(self.key_feat.capitalize()), title_fontsize='large')
            leg._legend_box.align = "left"