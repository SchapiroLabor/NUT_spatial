#########################################################
##TUMOR CONTOUR DETECTION AND DISTANCE TO CONTOUR PLOTS##
#########################################################

#import necessary packages
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, distance_transform_edt
from scipy.ndimage import label as ndi_label
from skimage.measure import find_contours
from skimage.morphology import (
    binary_closing,
    remove_small_objects,
    disk,
)
import matplotlib.pyplot as plt


#rasterize the coordinates function

def _rasterize_counts(coords, x_min, y_min, nx, ny, grid_size):
    """Rasterize point counts onto a (ny, nx) grid."""
    if coords.size == 0:
        return np.zeros((ny, nx), dtype=float)
    gx = np.clip(((coords[:, 0] - x_min) / grid_size).astype(int), 0, nx - 1)
    gy = np.clip(((coords[:, 1] - y_min) / grid_size).astype(int), 0, ny - 1)
    grid = np.zeros((ny, nx), dtype=float)
    np.add.at(grid, (gy, gx), 1.0)
    return grid
   

#tumor mask creation and calculating distances
def tumor_mask_and_distance_to_contour(
    adata,
    *,
    img_id="patient_exp",
    x="Cell Center X",
    y="Cell Center Y",
    cell_type="cell_category",
    tumor_cell_name="Tumor cells",
    # raster / smoothing
    grid_size=30,
    gaussian_sigma=2,
    # fraction model
    alpha_non_tumor=1,
    frac_thresh=0.35,
    # tissue mask
    tissue_thresh=0.05,
    # minimum evidence gate
    tumor_present_min_count=1,
    tumor_present_dilate=1,
    # morphology
    closing_radius=2,
    min_obj_px=50,      
    min_islet_px=50,    
    # plotting
    plot=True,
    point_size=3,
    fig_size=(7, 7),
    contour_lw=2,
    resolution = 0.17,
    show_titles=False,
    show_legend=False,
    save_svg=False,
    svg_dir=None,
    dpi=300,
):
    """
    Outputs added to adata.obs (per cell):
      - tumor_territory (bool): inside cleaned tumor mask
      - tumor_islet_id (int): 0 outside, 1..K inside
      - tumor_signed_dist_to_contour (float): signed distance to contour in ORIGINAL coordinate units
          negative inside, positive outside
      - tumor_dist_to_contour (float): absolute distance to contour
      - tumor_total_area_um2 (float): total area of all tumor islets combined (in µm²), same value for all cells in the same image
      - non_tumor_area_um2 (float): total tissue area outside the tumor mask (in µm²), same value for all cells in the same image
    """

    # init output cols
    for col, default in [
        ("tumor_territory", False),
        ("tumor_islet_id", 0),
        ("tumor_signed_dist_to_contour", np.nan),
        ("tumor_dist_to_contour", np.nan),
        ("tumor_total_area_um2", np.nan),
        ("non_tumor_area_um2", np.nan),     # <-- NEW
    ]:
        if col not in adata.obs.columns:
            adata.obs[col] = default

    imgs = adata.obs[img_id].unique()

    for img_val in imgs:
        df = adata.obs.loc[adata.obs[img_id] == img_val, [x, y, cell_type]].copy()
        if df.shape[0] == 0:
            continue

        coords = df[[x, y]].to_numpy(float)
        labels = df[cell_type].to_numpy()
        is_tumor = labels == tumor_cell_name

        tumor_coords = coords[is_tumor]
        non_tumor_coords = coords[~is_tumor]

        # grid bounds
        x_min, y_min = coords.min(axis=0)
        x_max, y_max = coords.max(axis=0)
        nx = int(np.ceil((x_max - x_min) / grid_size)) + 1
        ny = int(np.ceil((y_max - y_min) / grid_size)) + 1
        extent = [x_min, x_max, y_max, y_min]

        # raster counts
        tumor_grid = _rasterize_counts(tumor_coords, x_min, y_min, nx, ny, grid_size)
        non_tumor_grid = _rasterize_counts(non_tumor_coords, x_min, y_min, nx, ny, grid_size)

        # smooth
        tumor_sm = gaussian_filter(tumor_grid, sigma=gaussian_sigma)
        non_tumor_sm = gaussian_filter(non_tumor_grid, sigma=gaussian_sigma)
        

        # fraction
        eps = 1e-6
        tumor_frac = tumor_sm / (tumor_sm + alpha_non_tumor * non_tumor_sm + eps)

        tumor_norm = tumor_sm / tumor_sm.max() if tumor_sm.max() > 0 else tumor_sm
        stroma_norm = non_tumor_sm / non_tumor_sm.max() if non_tumor_sm.max() > 0 else non_tumor_sm

        dominance = tumor_norm - 1.0 * stroma_norm
        tumor_mask = (tumor_norm >= 0.12) & (dominance >= 0) & (
            tumor_norm >= 1.2 * np.maximum(stroma_norm, 1e-6)
        ) & (tumor_frac >= frac_thresh)


        # morphology
        if min_obj_px and min_obj_px > 0:
            tumor_mask = remove_small_objects(tumor_mask.astype(bool), min_size=int(min_obj_px))
        if closing_radius and closing_radius > 0:
            tumor_mask = binary_closing(tumor_mask, footprint=disk(int(closing_radius)))


        if plot:
            plt.figure(figsize=fig_size)
            plt.imshow(tumor_mask, origin="lower", cmap="gray", extent=extent)
            plt.axis("off")
            if save_svg:
                if svg_dir is None:
                    raise ValueError("svg_dir must be provided if save_svg=True")

                roi_tag = str(img_val).replace("/", "_")
                filename = f"{svg_dir}/{roi_tag}_tumor_mask.pdf"
                plt.savefig(filename, format="pdf", dpi=dpi, bbox_inches="tight")
                print(f"[EXPORT] Saved tumor mask plot for {roi_tag} to {filename}")
            plt.show()
        
        # connected components (islets)
        comp_map, n_comp = ndi_label(tumor_mask, structure=np.ones((3, 3), dtype=int))
        if n_comp == 0:
            if plot:
                print(f"[{img_val}] No tumor mask after filtering.")
            continue

        comp_sizes = np.bincount(comp_map.ravel())
        valid = np.where(comp_sizes >= int(min_islet_px))[0]
        valid = valid[valid != 0]

        tumor_mask_clean = np.isin(comp_map, valid)
        comp_map_clean = comp_map * tumor_mask_clean

        # Each grid pixel covers (grid_size * resolution)² µm²
        pixel_area_um2 = (grid_size * resolution) ** 2
        total_tumor_px = int(tumor_mask_clean.sum())
        total_tumor_area_um2 = total_tumor_px * pixel_area_um2

        # for non-tumor (stroma) area
        # Tissue mask: any grid pixel that has at least one cell of any type.
        # This avoids counting empty background pixels at the edges of the bounding box.
        all_grid = _rasterize_counts(coords, x_min, y_min, nx, ny, grid_size)
        tissue_mask = all_grid > 0
        non_tumor_area_um2 = int((tissue_mask & ~tumor_mask_clean).sum()) * pixel_area_um2  # <-- NEW

        # signed distance to CONTOUR via EDT
        # distance to nearest outside pixel for inside points, and to nearest inside pixel for outside points
        dist_out = distance_transform_edt(~tumor_mask_clean) * grid_size   # outside -> distance to tumor
        dist_in  = distance_transform_edt(tumor_mask_clean) * grid_size    # inside -> distance to outside boundary

        
        # signed: inside negative (distance to boundary), outside positive
        signed_dist = dist_out.copy()
        signed_dist[tumor_mask_clean] = -dist_in[tumor_mask_clean]
        
        # Distance field (outside tumor)
        D_out = dist_out.astype(np.float32)
        
        # Unique ROI tag
        roi_tag = str(img_val).replace("/", "_")
        
        # Save coordinate metadata
        meta = {
            "x_min": float(x_min),
            "y_min": float(y_min),
            "grid_size": float(grid_size),
            "nx": int(nx),
            "ny": int(ny)
        }
        

###########

#mapping grid to cells

        gx = np.clip(((coords[:, 0] - x_min) / grid_size).astype(int), 0, nx - 1)
        gy = np.clip(((coords[:, 1] - y_min) / grid_size).astype(int), 0, ny - 1)

        out = pd.DataFrame(index=df.index)
        out["tumor_territory"] = tumor_mask_clean[gy, gx]
        out["tumor_islet_id"] = comp_map_clean[gy, gx].astype(int)
        out["tumor_signed_dist_to_contour"] = signed_dist[gy, gx]
        out["tumor_dist_to_contour"] = np.abs(out["tumor_signed_dist_to_contour"].to_numpy(float))
        out["tumor_signed_dist_to_contour_um"] = out["tumor_signed_dist_to_contour"] * resolution 
        out["tumor_dist_to_contour_um"] = out["tumor_dist_to_contour"] * resolution
        out["tumor_total_area_um2"] = total_tumor_area_um2
        out["non_tumor_area_um2"] = non_tumor_area_um2          # <-- NEW

        adata.obs.loc[out.index, out.columns] = out.values

        # store in adata.uns for easy image-level lookup
        adata.uns.setdefault("tumor_total_area_um2", {})[img_val] = total_tumor_area_um2
        adata.uns.setdefault("non_tumor_area_um2", {})[img_val] = non_tumor_area_um2      # <-- NEW

        # --- plot: tumor vs non-tumor + contour overlay ---
        if plot:
            # contour lines in grid coords -> convert to original coords
            contours = find_contours(tumor_mask_clean.astype(float), level=0.5)
            if "tumor_contours" not in adata.uns:
                adata.uns["tumor_contours"] = {}
            adata.uns["tumor_contours"][img_val] = contours
            adata.uns.setdefault("tumor_contours_grid_size", {})[img_val] = grid_size
            adata.uns.setdefault("tumor_contours_origin", {})[img_val] = (x_min, y_min)

            plt.figure(figsize=fig_size)

            # non-tumor first (background)
            plt.scatter(
                non_tumor_coords[:, 0], non_tumor_coords[:, 1],
                s=point_size, c="lightgrey", linewidth=0, label="Non-tumor"
            )
            # tumor on top
            plt.scatter(
                tumor_coords[:, 0], tumor_coords[:, 1],
                s=point_size, c="#44AA99", linewidth=0, label="Tumor"
            )

            # overlay contours
            for cont in contours:
                # cont is (row=gy, col=gx) in grid space
                yy = y_min + cont[:, 0] * grid_size
                xx = x_min + cont[:, 1] * grid_size
                plt.plot(xx, yy, linewidth=contour_lw, color="#1f77b4")

            # plt.gca().invert_yaxis()
            plt.gca().set_aspect("equal")
            if show_titles:
                plt.title(f"{img_val}: Tumor vs Non-tumor + Tumor contour")
            plt.axis("off")
            if show_legend:
                plt.legend(frameon=False, loc="upper right")
            if save_svg:
                if svg_dir is None:
                    raise ValueError("svg_dir must be provided if save_svg=True")

                roi_tag = str(img_val).replace("/", "_")
                filename = f"{svg_dir}/{roi_tag}_tumor_contour.pdf"
                plt.savefig(filename, format="pdf", dpi=dpi, bbox_inches="tight")
                print(f"[EXPORT] Saved tumor contour plot for {roi_tag} to {filename}")
            plt.show()

            plt.figure(figsize=fig_size)

            # non-tumor first (background)
            plt.scatter(
                non_tumor_coords[:, 0], non_tumor_coords[:, 1],
                s=point_size, c="lightgrey", linewidth=0, label="Non-tumor"
            )
            # tumor on top
            plt.scatter(
                tumor_coords[:, 0], tumor_coords[:, 1],
                s=point_size, c="#44AA99", linewidth=0, label="Tumor"
            )

            # plt.gca().invert_yaxis()
            plt.gca().set_aspect("equal")
            if show_titles:
                plt.title(f"{img_val}: Tumor vs Non-tumor")
            plt.axis("off")
            if show_legend:
                plt.legend(frameon=False, loc="upper right")
            if save_svg:
                if svg_dir is None:
                    raise ValueError("svg_dir must be provided if save_svg=True")

                roi_tag = str(img_val).replace("/", "_")
                filename = f"{svg_dir}/{roi_tag}_tumor_non_tumor.pdf"
                plt.savefig(filename, format="pdf", dpi=dpi, bbox_inches="tight")
                print(f"[EXPORT] Saved tumor vs non-tumor plot for {roi_tag} to {filename}")
            plt.show()

    return adata

def add_tumor_layers(
    adata,
    dist_col="tumor_signed_dist_to_contour_um",
    step=25,
    out_col="tumor_layer",
):
    d = adata.obs[dist_col].to_numpy()
    layer = np.empty(len(d), dtype=object)

    # OUTSIDE tumor (positive)
    outside = d > 0
    out_bin = np.floor((d[outside] - 1e-12) / step).astype(int)
    conditions = [out_bin >= 4, out_bin == 3, out_bin == 2, out_bin == 1, out_bin == 0]
    choices = ["stroma distal", "stroma 75-100 µm", "stroma 50-75 µm", "stroma 25-50 µm", "stroma 0-25 µm"]
    out_labels = np.select(conditions, choices, default="")
    layer[outside] = out_labels

    # INSIDE tumor (negative)
    inside = d < 0
    depth_in = -d[inside]
    in_bin = np.floor((depth_in - 1e-12) / step).astype(int)
    conditions = [in_bin >= 1, in_bin == 0]
    choices = ["tumor 25+ µm", "tumor 0-25 µm"]
    in_labels = np.select(conditions, choices, default="")
    layer[inside] = in_labels

    #stroma_tumor_interface
    layer[d == 0] = "stroma 0-25 µm"


    ordered = ["stroma distal", "stroma 75-100 µm", "stroma 50-75 µm", "stroma 25-50 µm", "stroma 0-25 µm", "tumor 0-25 µm", "tumor 25+ µm"]
    adata.obs[out_col] = pd.Categorical(layer, categories=ordered, ordered=True)

    return adata


def compute_patient_layer_enrichment(
    adata,
    patient_col = "patient_ID",
    exp_col="patient_exp",
    layer_col="tumor_layer",
    cell_col="cell_category"
):
    df = adata.obs[[patient_col, exp_col, layer_col, cell_col]].copy()

    results = []

    for patient, df_p in df.groupby(exp_col):

        # global cell-type frequency for this patient
        global_freq = df_p[cell_col].value_counts(normalize=True)

        for layer, df_l in df_p.groupby(layer_col):
            layer_freq = df_l[cell_col].value_counts(normalize=True)

            for cell in global_freq.index:
                layer_freq_cell = layer_freq.get(cell, 0)
                global_freq_cell = global_freq[cell]
                fe = layer_freq.get(cell, 0) / global_freq[cell]

                results.append({
                    "patient_ID" : df_p[patient_col].iloc[0],
                    "patient_exp": patient,
                    "tumor_layer": layer,
                    "cell_type": cell,
                    "layer_proportion" : layer_freq_cell,
                    "global_proportion" : global_freq_cell,
                    "fold_enrichment": fe,
                    "log2_fold_enrichment": np.log2(fe) if fe > 0 else np.nan,
                    "n_cells_layer": len(df_l),
                    "n_cells_patient": len(df_p)
                })

    return pd.DataFrame(results)

