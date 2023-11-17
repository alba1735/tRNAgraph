
#!/usr/bin/env python3

import seaborn as sns
import numpy as np
import pandas as pd
import anndata as ad
import argparse

import toolsDirectory
import toolsAnalysis

import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

import logomaker

class visualizer():
    def __init__(self,adata,save_path='seq_logos',cutoff=0,pal_dict=None,cov_mask=False):
        os.makedirs(save_path, exist_ok=True)
        self.adata = adata
        self.save_path = save_path
        self.cutoff = cutoff
        self.pal_dict = pal_dict
        self.cov_mask = cov_mask
        self.multi_plot = False
        self.divergence_plot = False
        
    def multi(self,category):
        self.multi_plot = True
        cat_list = sorted(adata.obs[category].values.unique())
        for i in cat_list:
            plt = self.single(category,i)
            if category == 'trna':
                plt2 = self.raw_sequence(category,i)
    
    def divergence(self,category1,unit1,category2,unit2):
        self.divergence_plot = True
        raw_1 = self.list_gen(category1,unit1)
        raw_2 = self.list_gen(category2,unit2)
        
        df_1,cov_df_1 = self.norm_counts(raw_1,category1,unit1,pseudo=20)
        df_2,cov_df_2 = self.norm_counts(raw_2,category2,unit2,pseudo=20)
        
        df_1 = logomaker.transform_matrix(df_1,from_type='counts',to_type='information')
        df_2 = logomaker.transform_matrix(df_2,from_type='counts',to_type='information')
        
        df = df_1 - df_2
        
        title = 'divergence_{}_{}_{}_{}'.format(category1,unit1,category2,unit2)
        
        # Masking of the DF with tDRs using cov_df_1
        if self.cov_mask:
            for x in range(len(df)):
                if cov_df_1[x] == False: df.iloc[x] = [0,0,0,0]
            title = 'divergence_{}_{}_{}_{}_cov_{}'.format(category1,unit1,category2,unit2,self.cov_mask)
   
        #df = self.freq_matrix(df)
        #df = df.clip_lower(0)
        #df = logomaker.transform_matrix(df,from_type='counts',to_type='information')
        
        self.logo_plot(df,cov_df_1,'counts',title)
        
    
    def single(self,category,unit):
        raw_l = self.list_gen(category,unit)
        
        df,cov_df = self.norm_counts(raw_l,category,unit,pseudo=20)
        df = logomaker.transform_matrix(df,from_type='counts',to_type='information')
        
        title = '{}_{}_normalized_bits'.format(category,unit)
        
        if self.cov_mask:
            for x in range(len(df)):
                if cov_df[x] == False: df.iloc[x] = [0,0,0,0]
            title = '{}_{}_normalized_bits_cov_{}'.format(category,unit,self.cov_mask)
        
        self.logo_plot(df,cov_df,'information',title)
        
        df,cov_df = self.norm_counts(raw_l,category,unit,pseudo=0)
        df = self.freq_matrix(df)
        
        title = '{}_{}_frequency_cutoff_{}'.format(category,unit,self.cutoff)
        
        if self.cov_mask:
            for x in range(len(df)):
                if cov_df[x] == False: df.iloc[x] = [0,0,0,0]
            title = '{}_{}_frequency_cutoff_{}_cov_{}'.format(category,unit,self.cutoff,self.cov_mask)
        
        self.logo_plot(df,cov_df,'counts',title)
        
    def raw_sequence(self,category,unit):
        raw_l = self.list_gen(category,unit)
        self.seq_plot(raw_l,category,unit,'counts','{}_{}_sequence'.format(category,unit))
        
    def list_gen(self,category,unit):
        raw_l = [i for i in adata.obs[adata.obs[category] == unit].refseq.values]
        read_l = [int(i) for i in adata.obs[adata.obs[category] == unit].n_reads.values]
        #norm_l = []

        #for i in range(len(read_l)):
        #    norm_l += [raw_l[i]] * read_l[i]
            
        return raw_l
    
    def freq_matrix(self,df):
        fsize = df.max().max()
        for k,v in df.iterrows():
            tl = []
            for i in v:
                if i/fsize <= self.cutoff:
                    tl.append(0)
                else:
                    tl.append(i*2/fsize)
            df.at[k,['A','C','G','T']] = tl
        return df
    
    def norm_counts(self,logo_list,category,unit,pseudo=0):
        logo_list = [i for i in self.adata.obs[self.adata.obs[category] == unit].refseq.values]
        
        pos_dict = {'A':0,'C':1,'G':2,'T':3}
        mxa = np.full((len(logo_list[0]),4),pseudo)

        mask_out = np.isin(self.adata.obs[category], [unit])
        mx = self.adata[mask_out,:]
        
        mx = mx.raw[:,np.isin(mx.raw.var.gap.values, [False])]
        mx = mx[:,np.isin(mx.var.coverage.values, ['uniquecoverage'])]
        
        for i in range(len(logo_list)):    
            for j in range(len(mx.X[i])):
                if pos_dict.get(logo_list[i][j]) != None:
                    mp = pos_dict.get(logo_list[i][j])
                    mxa[j][mp] += mx.X[i][j]          
        
        cov_df = pd.DataFrame(mx.X).mean(axis=0)
        cov_df = [True if i >= self.cov_mask else False for i in cov_df/cov_df.max()]
        
        return pd.DataFrame(mxa, index=list(range(0,77)), columns=['A','C','G','T']),cov_df
        #return pd.DataFrame(mxa, index=mx.var.positions.values, columns=['A','C','G','T'])

    def logo_plot(self,df,cov_df,seq_type,title='seq_logo'):
        fig, ax = plt.subplots(figsize=(9,3))
        if self.pal_dict:
            logo_ax = logomaker.Logo(df,color_scheme=pal_dict,ax=ax)
        else:
            logo_ax = logomaker.Logo(df,ax=ax)

        logo_ax.style_spines(visible=False)
        ax.set_xlim([-1, 77])
        ax.set_xticks([])
        ax.set_title(title, fontsize='x-large')
        
        if self.divergence_plot == True:
            logo_ax.style_spines(spines=['left'], visible=True, bounds=[-2,2])
            y = -2.25
            yys = 0.935
            ax.set_ylim([-2.5, 2])
            ax.set_yticks([-2,-1,0,1,2])
            ax.set_yticklabels(['-2','-1','0','1','2'])
            ax.set_ylabel("       Bits",fontsize=12,verticalalignment='center')
        else:
            logo_ax.style_spines(spines=['left'], visible=True, bounds=[0, 2])
            y = -.16
            yys = 0.6
            if seq_type == 'information':
                ax.set_ylim([-0.25, 2])
                ax.set_yticks([0,1,2])
                ax.set_yticklabels(['0','1','2'])
                ax.set_ylabel("       Bits",fontsize=12,verticalalignment='center')
            elif seq_type == 'counts':
                ax.set_ylim([-0.25, 2])
                ax.set_yticks([0,2])
                ax.set_yticklabels(['0%','100%'])
                ax.set_ylabel("       Frequency",fontsize=12,verticalalignment='center')    
                if self.cutoff > 0:
                    ax.plot([-1,77],[self.cutoff*2,self.cutoff*2],linewidth=1,ls='--',color='black')
            
        # set parameters for drawing gene
        
        xs = np.arange(-3, len(df),10)
        ys = y*np.ones(len(xs))
        
        # draw trna scheme
        ax.axhline(y, color='k', linewidth=1)
        ax.plot(xs, ys, marker='4', linewidth=0, markersize=7, color='k')
        ax.plot([9.5, 25.5], [y, y], color='k', linewidth=12, solid_capstyle='butt')
        ax.plot([26.5, 43.5], [y, y], color='k', linewidth=12, solid_capstyle='butt')
        ax.plot([48.5, 65.5], [y, y], color='k', linewidth=12, solid_capstyle='butt')
        if self.cov_mask:
            fvar = False
            cov_map = []
            for x in range(len(cov_df)):
                if cov_df[x] != fvar:
                    fvar = cov_df[x]
                    cov_map.append(x)
                    
            if len(cov_map)%2 == 1: cov_map.append(77)
            cov_map = [(cov_map[i],cov_map[i+1]) for i in range(0,len(cov_map),2)]
            for x1,x2 in cov_map:
                ax.plot([x1-0.5, x2-0.5], [y+0.35, y+0.35], color='orange', alpha=0.5, linewidth=12, solid_capstyle='butt')

            ax.text(np.mean(cov_map[0]),(yys*y)+0.35,'tDR Coverage',fontsize=12,verticalalignment='top',horizontalalignment='center',color='black')

        # annotate trna scheme
        ax.text(17.5,yys*y,'D-Arm/Loop',fontsize=12,verticalalignment='top',horizontalalignment='center',color='white')
        ax.text(35,yys*y,'A-Arm/Loop',fontsize=12,verticalalignment='top',horizontalalignment='center',color='white')
        ax.text(57,yys*y,'T-Arm/Loop',fontsize=12,verticalalignment='top',horizontalalignment='center',color='white')

        logo_ax.highlight_position_range(pmin=10, pmax=25, color='#f0f0f0')
        logo_ax.highlight_position_range(pmin=27, pmax=43, color='#f0f0f0')
        logo_ax.highlight_position_range(pmin=49, pmax=65, color='#f0f0f0')
        
        logo_ax.highlight_position_range(pmin=14, pmax=21, color='#cacaca')
        logo_ax.highlight_position_range(pmin=32, pmax=38, color='#cacaca')
        logo_ax.highlight_position_range(pmin=54, pmax=60, color='#cacaca')

        logo_ax.fig.tight_layout()
        
        title_save = title.replace(' ','_').replace('(','').replace(')','').lower()
            
        plt.savefig('{}/{}.pdf'.format(self.save_path,title_save),bbox_inches='tight')
        
        if self.multi_plot == True:
            plt.close()
            
    def seq_plot(self,logo_list,category,unit,seq_type,title='sequence'):
        df = self.norm_counts(logo_list,category,unit,pseudo=0)
        df[df > 0] = 1

        fig, ax = plt.subplots(figsize=(9,1.5))
        if self.pal_dict:
            logo_ax = logomaker.Logo(df,color_scheme=pal_dict,show_spines=False,ax=ax)
        else:
            logo_ax = logomaker.Logo(df,show_spines=False,ax=ax)

        # style using Axes methods
        ax.set_xlim([0, 77])
        ax.set_xticks([1,37,58,73])
        ax.set_xticklabels(['1','37','58','73'])
        ax.set_title(title, fontsize='x-large')
        ax.set_yticks([])
        ax.set_yticklabels([])

        logo_ax.highlight_position_range(pmin=10, pmax=25, color='#F0F0F0')
        logo_ax.highlight_position_range(pmin=27, pmax=43, color='#F0F0F0')
        logo_ax.highlight_position_range(pmin=49, pmax=65, color='#F0F0F0')
        
        logo_ax.fig.tight_layout()

        title_save = title.replace(' ','_').replace('(','').replace(')','').lower()
        plt.savefig('{}/{}.pdf'.format(self.save_path,title_save),bbox_inches='tight')
        
        if self.multi_plot == True:
            plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='logo_tools.py',
        description='Generate seqlogos from adata.',
    )

    parser.add_argument('-i', '--anndata', help='Specify AnnData input', required=True)
    parser.add_argument('-o', '--output', help='Specify output directory', default='barplot', required=False)
    parser.add_argument('--bargrp', help='Specify AnnData column to group by (default: group) (optional)', default='group', required=False)
    parser.add_argument('--colormap_tc', help='Specify a colormap for coverage plots for type counts (optional)', default=None)
    parser.add_argument('--colormap_bg', help='Specify a colormap for coverage plots for bargrp (optional)', default=None)

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    toolsDirectory.builder(args.output)

    adata = ad.read_h5ad(args.anndata)

    visualizer(adata, args.bargrp, args.colormap_tc, args.colormap_bg, args.output)