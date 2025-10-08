from skimage import measure, morphology
from skimage.filters import threshold_otsu
from shapely.geometry import Polygon, Point
from scipy.ndimage import gaussian_filter, distance_transform_edt
import matplotlib.pyplot as plt
import numpy as np

def extract_all_tumor_cores_otsu(
    adata, 
    image_size=256, 
    sigma=1, 
    min_size=50, 
    resolution_um_per_px=0.17, 
    plot=True, 
    save=True,
    results_dir=None
):
    
    # get patient sample name
    patient_exp = adata.obs['patient_exp'].unique()[0]

    # Select tumor cells coordinates
    df = adata.obs[adata.obs['cell_type'] == 'Tumor_cells']
    x_coords = np.round(df['Cell Center X'].values)
    y_coords = np.round(df['Cell Center Y'].values)
    
    # determine the dimension of the image
    x_min, x_max = x_coords.min(), x_coords.max()
    y_min, y_max = y_coords.min(), y_coords.max()

    # Scale tthe tumor cell coordinates
    scale_x = (x_coords - x_min) / (x_max - x_min) * (image_size - 1)
    scale_y = (y_coords - y_min) / (y_max - y_min) * (image_size - 1)

    # Rasterize the scaled coordinates into binary image
    image = np.zeros((image_size, image_size), dtype=np.uint8)
    for x_pix, y_pix in zip(scale_x.astype(int), scale_y.astype(int)):
        image[y_pix, x_pix] = 255

    # Gaussian blur and threshold
    blurred_image = gaussian_filter(image, sigma=sigma)
    blurred_image = blurred_image * (blurred_image > np.percentile(blurred_image, 10))
    otsu_thresh = threshold_otsu(blurred_image) * 0.5 #a little threshold corrector
    binary_mask = blurred_image > otsu_thresh
    
    # Clean small objects
    cleaned_mask = morphology.remove_small_objects(binary_mask, min_size=min_size)

    # Pad mask to prevent edge artifacts in contour finding
    padded_mask = np.pad(cleaned_mask, pad_width=1, mode='constant', constant_values=0)

    # Find contours on padded mask
    contours = measure.find_contours(padded_mask.astype(float), 0.5)
    if len(contours) == 0:
        print("No contour found.")
        return None, None, None, None

    polygons = []
    for contour in contours:
        # Adjust coordinates due to padding (-1 pixel)
        contour_x = (contour[:, 1] - 1) / (image_size - 1) * (x_max - x_min) + x_min
        contour_y = (contour[:, 0] - 1) / (image_size - 1) * (y_max - y_min) + y_min
        poly = Polygon(zip(contour_x, contour_y))
        polygons.append(poly)

    n_cores = len(polygons)

    # Compute total area of polygons
    tumor_area = sum(poly.area for poly in polygons)

    # --- Compute distance map from tumor cores ---
    inverse_mask = ~cleaned_mask  # True outside tumor cores
    distance_map = distance_transform_edt(inverse_mask)  # distances in pixel units

    # Calculate distance for every cell in adata.obs
    all_x = np.round(adata.obs['Cell Center X'].values)
    all_y = np.round(adata.obs['Cell Center Y'].values)

    # Scale all cell coords
    scale_all_x = (all_x - x_min) / (x_max - x_min) * (image_size - 1)
    scale_all_y = (all_y - y_min) / (y_max - y_min) * (image_size - 1)

    pix_x = np.clip(scale_all_x.astype(int), 0, image_size - 1)
    pix_y = np.clip(scale_all_y.astype(int), 0, image_size - 1)

    # Lookup distance for each cell
    cell_distances_pixels = distance_map[pix_y, pix_x]

    # Convert pixel distances back to original coordinate units
    pixel_size_x = (x_max - x_min) / (image_size - 1)
    pixel_size_y = (y_max - y_min) / (image_size - 1)
    avg_pixel_size = (pixel_size_x + pixel_size_y) / 2
    cell_distances = cell_distances_pixels * avg_pixel_size
    cell_distances_um = cell_distances * resolution_um_per_px  # distances in microns

    # --- NEW: Assign tumor core ids + contour flags ---
    cells = [Point(x, y) for x, y in zip(all_x, all_y)]
    tumor_core_ids = []
    on_contour_flags = []

    # distance threshold for "on contour"
    threshold = 2 * ((x_max - x_min) / (image_size - 1))

    for cell in cells:
        found = False
        for i, poly in enumerate(polygons):
            if poly.contains(cell):
                tumor_core_ids.append(i)
                on_contour_flags.append(cell.distance(poly.exterior) <= threshold)
                found = True
                break
        if not found:
            tumor_core_ids.append(-1)  # outside all tumor cores
            on_contour_flags.append(False)


    # --- Plotting ---
    if plot:
        plt.rcParams['svg.fonttype'] = 'none'
        plt.rcParams['font.family'] = 'Arial'  
        fig, axs = plt.subplots(1, 4, figsize=(20, 5))

        axs[0].imshow(image, cmap='gray')
        axs[0].set_title('Tumor Cells')

        axs[1].imshow(blurred_image, cmap='gray')
        axs[1].set_title(f'Gaussian Blurred + Otsu Threshold')

        axs[2].imshow(cleaned_mask, cmap='gray')
        axs[2].set_title('Filled Mask after threshold + clean')

        axs[3].imshow(image, cmap='gray')
        for poly in polygons:
            x_img = (np.array(poly.exterior.xy[0]) - x_min) / (x_max - x_min) * (image_size - 1)
            y_img = (np.array(poly.exterior.xy[1]) - y_min) / (y_max - y_min) * (image_size - 1)
            axs[3].plot(x_img, y_img, color='red')
        axs[3].set_title('All Tumor Core Contours')

        for ax in axs:
            ax.axis('off')
        
        plt.show()
        if save:
            if results_dir is None:
                raise ValueError("Please provide `results_dir` if `save=True`.")
            else:
                fig.savefig(results_dir / f"{patient_exp}_tumor_cores.svg", format='svg', bbox_inches='tight')
                print(f"Saved tumor cores figure for {patient_exp} at {results_dir / f'{patient_exp}_tumor_cores.svg'}")
        else:
            print("Not saving figure.")

    return polygons, cell_distances_um, n_cores, tumor_area, on_contour_flags
