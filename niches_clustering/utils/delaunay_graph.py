### CREATES THE DELAUNAY GRAPH FOR THE DATA FOR SPATIAL NICHES ###
import numpy as np
from scipy.spatial import Delaunay
import pandas as pd


def delaunay_adata_internal(adata, coords_key = "spatial", phenotype = "cell_category"):
    points = adata.obsm[coords_key]
    tri = Delaunay(points)
    neighbours_dict = {i: set() for i in range(len(points))}
    for simplex in tri.simplices:
                for i in range(len(simplex)):
                    for j in range(i + 1, len(simplex)):
                        neighbours_dict[simplex[i]].add(simplex[j])
                        neighbours_dict[simplex[j]].add(simplex[i])
    neighbours_list = [list(neighbours) for neighbours in neighbours_dict.values()]
    # Ensure each list has the same number of elements by padding with -1 (assuming indices are non-negative)
    max_neigh_len = max(len(neigh) for neigh in neighbours_list)
    neighbours_list_padded = [neigh + [-1] * (max_neigh_len - len(neigh)) for neigh in neighbours_list]
    # Convert to numpy array for consistency with KNN method
    ind = np.array(neighbours_list_padded)
    # Convert to DataFrame for the same output format as the original function
    neighbours = pd.DataFrame(ind.tolist(), index=adata.obs.index)
    # Replace -1 with None
    neighbours.replace(-1, None, inplace=True)
    phenomap = dict(zip(list(range(len(ind))), adata.obs[phenotype])) # Used for mapping
    for i in neighbours.columns:
            neighbours[i] = neighbours[i].dropna().map(phenomap, na_action='ignore')
    # Collapse all the neighbours into a single column
    n = pd.DataFrame(neighbours.stack(), columns = ["neighbour_phenotype"])
    n.index = n.index.get_level_values(0) # Drop the multi index
    n = pd.DataFrame(n)
    n['order'] = list(range(len(n)))
    # Merge with real phenotype
    n_m = n.merge(adata.obs[phenotype], how='inner', left_index=True, right_index=True)
    n_m['neighbourhood'] = n_m.index
    n = n_m.sort_values(by=['order'])
    # Normalize based on total cell count
    k = n.groupby(['neighbourhood','neighbour_phenotype']).size().unstack().fillna(0)
    # k = k.div(k.sum(axis=1), axis=0)    
    return k
def delaunay_adata(adata, coords_key = "spatial", phenotype = "cell_category", image_id = "patient_exp", label = "spatial_count"):
    adata_list = [adata[adata.obs[image_id] == i] for i in adata.obs[image_id].unique()]
    r_delaunay_data = lambda x: delaunay_adata_internal(x, coords_key="spatial", phenotype="cell_category")
    all_data = list(map(r_delaunay_data, adata_list))
    result = []
    for i in range(len(all_data)):
        result.append(all_data[i])
    result = pd.concat(result, join='outer') 
    # Reindex the cells
    result = result.reindex(adata.obs.index)
    result = result.fillna(0)

    # Add to adata
    adata.uns[label] = result
    
    # Return        
    return adata
