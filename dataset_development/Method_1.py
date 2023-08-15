from osgeo import gdal
import numpy as np

AC_grid_path = 'data/0_source/PctAnnChance_clip.tif'
DV_grid_path = 'data/0_source/composite_DxV_clip.tif'
#DV_output_path = 'data/1_interim/Method_1_DV.tif'
#AC_output_path = 'data/1_interim/Method_1_AC.tif'
raster_band = 5 #This is 1% AEP
output_path = 'data/2_output/Method_1_Scores.tif'

# FUNCTION: N/A
#
# DESCRIPTION: Calculates the 'Method 1' score for a  multi-band DxV raster
#   and a single band Annual Chance raster. Prior to running this script, the 
#   rasters were clipped to identical boundaries and cell-size, ensuring that 
#   they align into a numpy array. A single band raster of calculated scores
#   is returned. Any value not in the defined ranges is assigned a value of 
#   zero. Equal weight was given to the two rasters. The option exists to 
#   export the interim rasters, but is commented out for now.
#
# INPUTS: 
#   AC_grid_path - the path to the single band Annual Chance raster
#   DV_grid_path - the path to the multi-band DxV raster
#   raster_band - the band to use from the DxV raster
#   output_path - the single-band score raster to be created
#
# CODE:
#
# --------------------------
# Process Annual Chance Grid
# --------------------------
 # Open the composite DV multiband raster dataset
AC_input_dataset = gdal.Open(AC_grid_path, gdal.GA_ReadOnly)
if AC_input_dataset is None:
    print("Failed to open the raster dataset.")

# Get the raster dimensions
AC_width = AC_input_dataset.RasterXSize
AC_height = AC_input_dataset.RasterYSize

# Read data
AC_band = AC_input_dataset.GetRasterBand(1)
AC_array = AC_band.ReadAsArray()

# Define condition
condition = (AC_array >= 0.2) & (AC_array <= 50)

# Define calculation for condition
calculation = 1.597626 * np.log(AC_array) + 11.172

# Apply conditional calculation
AC_result_array = np.where(condition, calculation, 0)

# For exporting interim raster:
''' # Create a new single-band raster dataset for the output
driver = gdal.GetDriverByName('GTiff')
AC_out_dataset = driver.Create(AC_output_path, AC_width, AC_height, 1, gdal.GDT_Float32)

# Write the scores to the single-band raster dataset
AC_band = AC_out_dataset.GetRasterBand(1)
AC_band.WriteArray(AC_result_array, 0, 0)

# Set the geotransform and projection of the output raster same as the input
AC_out_dataset.SetGeoTransform(AC_input_dataset.GetGeoTransform())
AC_out_dataset.SetProjection(AC_input_dataset.GetProjection())

# Close the datasets
AC_input_dataset = None
AC_out_dataset = None
'''


# ----------------
# Process DxV Grid
# ----------------
 # Open the composite DV multiband raster dataset
DV_input_dataset = gdal.Open(DV_grid_path, gdal.GA_ReadOnly)
if DV_input_dataset is None:
    print("Failed to open the raster dataset.")

# Get the number of bands and dimensions of the raster
DV_width = DV_input_dataset.RasterXSize
DV_height = DV_input_dataset.RasterYSize
    
DV_band = DV_input_dataset.GetRasterBand(raster_band)
DV_array = DV_band.ReadAsArray()

# Define conditions
condition1 = (DV_array > 0) & (DV_array <= 3.2)
condition2 = (DV_array > 3.2) & (DV_array <= 43.1)
condition3 = DV_array > 43.1

# Define calculations for each condition
calculation1 = DV_array * 0.7813
calculation2 = 2.89249 * np.log(DV_array) - 0.8825
calculation3 = 10

# Apply conditional calculations
DV_result_array = np.where(condition1, calculation1,
                    np.where(condition2, calculation2, 
                      np.where(condition3, calculation3, 0)))

# For exporting interim raster:
'''# Create a new single-band raster dataset for the output
driver = gdal.GetDriverByName('GTiff')
DV_out_dataset = driver.Create(DV_output_path, DV_width, DV_height, 1, gdal.GDT_Float32)

# Write the scores to the single-band raster dataset
score_band = DV_out_dataset.GetRasterBand(1)
score_band.WriteArray(DV_result_array, 0, 0)

# Set the geotransform and projection of the output raster same as the input
DV_out_dataset.SetGeoTransform(DV_input_dataset.GetGeoTransform())
DV_out_dataset.SetProjection(DV_input_dataset.GetProjection())

# Close the datasets
DV_input_dataset = None
DV_out_dataset = None
'''

# --------------
# Combine Scores
# --------------
# Check that the DxV and Annual Chance grids align
if AC_width == DV_width and AC_height == DV_height:

    # Create an empty array to store the individual classification results
    raster_values = np.zeros((2, AC_height, AC_width), dtype=np.float32)
    raster_values[0, :, :] = AC_result_array
    raster_values[1, :, :] = DV_result_array

    # Calculate the mean across all bands as the final score
    final_score = np.nanmean(raster_values, axis = 0)

    # Create a new single-band raster dataset for the output
    driver = gdal.GetDriverByName('GTiff')
    score_dataset = driver.Create(output_path, DV_width, DV_height, 1, gdal.GDT_Float32)

    # Write the scores to the single-band raster dataset
    score_band = score_dataset.GetRasterBand(1)
    score_band.WriteArray(final_score, 0, 0)

    # Set the geotransform and projection of the output raster same as the input
    score_dataset.SetGeoTransform(DV_input_dataset.GetGeoTransform())
    score_dataset.SetProjection(DV_input_dataset.GetProjection())

    # Close the datasets
    dataset = None
    score_dataset = None
        
    print("Scores calculated and saved to:", output_path)
else:
    print("Raster shapes to not match. Unable to calculate score.")