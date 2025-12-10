import numpy as np
import anndata as ad
import pandas as pd
import scanpy as sc
import scipy as sp
import os
import yaml
from pathlib import Path
from collections import OrderedDict

def extract_interaction(adata, selected_cell_type, excluded_cell_types = None, interaction_label = 'scimap_delaunay_conditional_zscore', direction = "from", cells_threshold = 10, save_data = True, results_dir = None):
    """
    Extract interaction data for a specific cell type from the adata object.
    
    Parameters:
    adata (AnnData): The annotated data object containing spatial interaction data.
    selected_cell_type (str): The cell type to extract interactions for
    excluded_cell_types (list): List of cell types to exclude from the interaction data.
    interaction_label (str): The label of the interaction data to extract.
    direction (str): The direction of interaction to consider ('from' or 'to').
    cells_threshold (int): The minimum number of cells required for an interaction to be considered valid.
        
    Returns:
    pd.DataFrame: A DataFrame containing the interaction data for the selected cell type.
    """
    # Extract interaction data for all cell types
    interaction_data = adata.uns[interaction_label]
    interaction_data.index = interaction_data["phenotype"] + "-" + interaction_data["neighbour_phenotype"]
    interaction_data= interaction_data.drop(columns=["phenotype", "neighbour_phenotype"])
    interaction_data = interaction_data.filter(like='zscore', axis=1)
    #arrange according to patient number for consistency 
    interaction_data.columns = interaction_data.columns.str.replace('zscore_', '')  # remove zscore_ from column names
    ordered_columns = adata.obs.sort_values(by = "patient_ID").loc[:,["exp_name", "patient_ID", "patient_exp"]].drop_duplicates().reset_index(drop=True)
    interaction_data = interaction_data[ordered_columns["patient_exp"]]
    interaction_data.columns = ordered_columns["patient_exp"].values
    # # Normalize by the total number of cells per sample
    total_cells = adata.obs.groupby("patient_exp")["cell_category"].count()[ordered_columns["patient_exp"].values]
    interaction_data = interaction_data.div(np.sqrt(total_cells), axis=1).multiply(1000)  # Scale to per 1000 cells
    #Filter for cell type of interest
    if direction == "from":
        interaction_data = interaction_data[interaction_data.index.str.contains(selected_cell_type + "-", regex=False)]
    elif direction == "to":
        interaction_data = interaction_data[interaction_data.index.str.contains("-" + selected_cell_type, regex = False)]
    #make data in long format
    interaction_data_long = interaction_data.reset_index().melt(id_vars='index', var_name='patient_exp', value_name='zscore')
    interaction_data_long.rename(columns={'index': 'cell_pair'}, inplace=True)
    #rearrange the columns
    interaction_data_long = interaction_data_long[['patient_exp', 'cell_pair', 'zscore']]
    #Get the names of the cell type without the neighbouring cell type
    if direction == "from":
        interaction_data_long['cell_category'] = interaction_data_long['cell_pair'].str.replace(selected_cell_type + "-", '', regex=False)
    elif direction == "to":
        interaction_data_long['cell_category'] = interaction_data_long['cell_pair'].str.replace('-' + selected_cell_type, '', regex=False)
    #Make any interactions with less than threshold cells NA
    total_cells = adata.obs.groupby(["patient_exp", "cell_category"]).size().reset_index(name='total_cells')
    interaction_data_long = interaction_data_long.merge(total_cells, on=['patient_exp', 'cell_category'], how='left')
    interaction_data_long.loc[interaction_data_long["total_cells"] <= cells_threshold, "zscore"] = np.nan
    interaction_data_long = interaction_data_long.drop(columns=["total_cells"])
    #make the values that are -Inf or Inf na
    interaction_data_long['zscore'] = interaction_data_long['zscore'].replace([np.inf, -np.inf], np.nan)
    #Take out excluded cell types if included
    if excluded_cell_types is not None:
        interaction_data_long = interaction_data_long[~interaction_data_long['cell_category'].isin(excluded_cell_types)]
    #Save the data if required
    if save_data:
        interaction_data_long.to_csv(results_dir / f"{selected_cell_type}_interaction_{direction}.csv", index=False)
    return interaction_data_long 