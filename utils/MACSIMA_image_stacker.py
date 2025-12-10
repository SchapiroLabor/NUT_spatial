import os
import numpy as np
import tifffile

#Stack up images function for MACSIMA data
def stack_images(image_dir, channels_list, output_path):
    #find the images that include the channels of interest
    image_files = sorted([f for f in os.listdir(image_dir) if f.endswith('.tif')])
    final_images_files = []
    final_channel_names = []
    #if among the channels there is dapi we only take one dapi image
    if "DAPI" in channels_list:
        one_dapi_file = [f for f in image_files if "-DAPI.tif" in f][1]
        final_images_files.append(one_dapi_file)
        final_channel_names.append("DAPI")
    for ch in channels_list:
        if ch != "DAPI":
            selected_image_files = [f for f in image_files if ch in f]
            #if there is more than one image of that channel append the channel name twice
            if len(selected_image_files) > 1:
                final_channel_names.extend([ch]*len(selected_image_files))
                final_images_files.extend(selected_image_files)
            else:
                final_channel_names.append(ch)
                final_images_files.extend(selected_image_files)
    #stack the images
    img_stacks = []
    for image in final_images_files:
        img_path = os.path.join(image_dir, image)
        img = tifffile.imread(img_path)
        img_stacks.append(img)
    # Stack images 
    stacked_image = np.stack(img_stacks, axis=0)
    # Add metadata
    metadata = {'channels': final_channel_names, "axes": "CYX"}
    # Save the final stacked image
    tifffile.imwrite(output_path, stacked_image)
    return stacked_image, final_images_files, final_channel_names