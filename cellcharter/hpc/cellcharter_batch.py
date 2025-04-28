import numpy as np
import pandas as pd
import os
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import cellcharter as cc
import squidpy as sq
from multiprocessing import Pool 


# Load the data
adata = ad.read_h5ad("adata.h5ad")

#load the autok
model_params = {
    'random_state': 42,
    'trainer_params': {
        'accelerator':'cpu',
        'enable_progress_bar': False
        },
    }

autok = cc.tl.ClusterAutoK(n_clusters=(4,25), model_class=cc.tl.GaussianMixture, 
                           model_params=model_params, 
                           max_runs=5)

# define function to run cluster stability analysis

def run_cs(exp):
    adata_exp = adata[adata.obs["exp_name"] == exp]
    sc.pp.scale(adata_exp)
    sq.gr.spatial_neighbors(adata_exp, library_key='exp_name', coord_type='generic', delaunay=True)
    adata_exp.uns["spatial"] = {roi: {} for roi in adata_exp.obs['exp_name'].unique()}
    cc.gr.remove_long_links(adata_exp)
    cc.gr.aggregate_neighbors(adata_exp, n_layers=3)
    autok.fit(adata_exp, use_rep='X_cellcharter')
    ax = cc.pl.autok_stability(autok, return_ax=True, save=f'figures/autok_stability_{exp}.png')
    ax.figure.set_size_inches(12, 6)
    return ax


# Run the function in parallel

with Pool() as pool:
     ax_list = pool.map(run_cs, adata.obs["exp_name"].unique())
