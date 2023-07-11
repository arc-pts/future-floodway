import arcpy
import os
import shutil
from os import PathLike, path

def zScorePerChunkSingleband(value_raster: str, chunk_polygon: str, output_folder: PathLike) -> PathLike:
    '''
    A generic helper function to calculate z-score values per raster cell, per chunk
    A z-score refers to the number of standard deviations each data value is from the mean, with a z-score of zero indicating the exact mean
    Inputs: 
        value_raster: Single-band raster
        chunk_polygon: Polygon shapefile containing raster chunk polygons

    Output:
        Single-band with z-score values per cell
    '''
        
    # Environment Settings
    desc = arcpy.Describe(value_raster)
    spatialReference = desc.spatialReference
    cellSize = desc.meanCellHeight
    extent = desc.extent

    if output_folder:
        workspace = output_folder
    else:
        workspace = desc.path

    arcpy.env.snapRaster = value_raster
    arcpy.env.outputCoordinateSystem = spatialReference
    # arcpy.env.extent = extent
    arcpy.env.cellSize = cellSize
    arcpy.env.workspace = workspace
    arcpy.env.overwriteOutput = True

    # Populate IDs - to be used for join
    try:
        arcpy.AddField_management(chunk_polygon, "Join_ID", "SHORT")
    except:
        arcpy.AddMessage("Id field already exists")

    with arcpy.da.UpdateCursor(chunk_polygon, ['FID','Join_ID']) as cursor:
        for row in cursor:
            row[1] = row[0] # Populate Id field with FID
            cursor.updateRow(row)


    # Zonal stats as table to get mean and stdev for each polygon chunk. Assign values back to polygons
    arcpy.AddMessage("Calcuating Stats")
    try:
        arcpy.DeleteField_management(chunk_polygon, ["MEAN", "STD"])
    except:
        arcpy.AddMessage("Previous stats fields not present")

    statsTable = workspace + "\statsTable.dbf"
    arcpy.sa.ZonalStatisticsAsTable(chunk_polygon, "Join_ID", value_raster, statsTable)
    arcpy.JoinField_management(chunk_polygon, "Join_ID", statsTable, "Join_ID", ["MEAN", "STD"])

    # Create folder to hold individual z score rasters
    temp_folder = workspace + r"\TempFolder"
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)

    # To deal with overlapping polygons
    # For each polygon, developed the z-score raster to be mosaiced together using MAX operator at the end
    chunk_lyr = arcpy.MakeFeatureLayer_management(chunk_polygon,"chunk_lyr")
    arcpy.AddMessage("Calculating Z-Scores per chunk")
    with arcpy.da.SearchCursor(chunk_polygon, ['OID@', "MEAN", "STD"]) as cursor:
        for row in cursor:
            # Select feature and create MEAN and STD rasters
            arcpy.AddMessage("Computing feature {0}".format(row[0]))
            query = '"FID" = ' + str(row[0])
            arcpy.SelectLayerByAttribute_management(chunk_lyr, "NEW_SELECTION", query)
            mean_raster = arcpy.conversion.FeatureToRaster(chunk_lyr, "MEAN", workspace + r"\MEAN.tif")
            std_raster = arcpy.conversion.FeatureToRaster(chunk_lyr, "STD", workspace + r"\STD.tif")

            # Raster Calculator - calculate z-score
            # I cannot get this raster calculator function to work! Error below
            # The parameter is incorrect. Parameter 'Rasters' is missing or invalid. Bind failed in function 'Raster Calculator Function'
            # z_score_raster = arcpy.sa.RasterCalculator(rasters=[value_raster, mean_raster, std_raster], input_names=["VALUE", "MEAN", "STD"], expression="(VALUE-MEAN)/STD", extent_type="IntersectionOf")
            value_raster_Extract = arcpy.sa.ExtractByMask(value_raster, chunk_lyr)
            outMinus = arcpy.sa.Minus(value_raster_Extract, mean_raster)
            z_score_raster_temp = arcpy.sa.Divide(outMinus, std_raster)
            z_score_raster_temp.save(temp_folder + r"\z_score_raster" + str(row[0]) + r".tif")

    arcpy.AddMessage("Mosaicing rasters")
    rasters_to_mosaic = [temp_folder + r"\\" +  f for f in os.listdir(temp_folder) if f.endswith(".tif")]
    z_score_raster_mosaic = arcpy.MosaicToNewRaster_management(rasters_to_mosaic, workspace, "z_score_raster_mosaic.tif", spatialReference,
                                   "32_BIT_FLOAT", cellSize, number_of_bands=1, mosaic_method="MAXIMUM")
    z_score_raster = arcpy.sa.SetNull(z_score_raster_mosaic, z_score_raster_mosaic, "VALUE > 1000")
    z_score_raster.save(workspace + r"\z_score_raster.tif")

    # Delete Intermediate Data
    arcpy.Delete_management(statsTable)
    arcpy.Delete_management(mean_raster)
    arcpy.Delete_management(std_raster)
    arcpy.Delete_management(z_score_raster_mosaic)
    try:
        shutil.rmtree(temp_folder)
    except:
        arcpy.AddMessage("File in use")

    output = arcpy.Describe(z_score_raster).catalogPath
    arcpy.AddMessage("Output saved to {0}".format(output))
    return output


def addToMap(layer):
    arcpy.AddMessage("Adding to Map")
    m = arcpy.mp.ArcGISProject("CURRENT").listMaps('*')[0]
    m.addDataFromPath(arcpy.Describe(layer).catalogPath)


if __name__ == "__main__":

    value_raster = arcpy.GetParameterAsText(0)
    chunk_polygon = arcpy.GetParameterAsText(1)
    output_folder = arcpy.GetParameterAsText(2)

    z_score_raster = zScorePerChunkSingleband(value_raster, chunk_polygon, output_folder)
    addToMap(z_score_raster)