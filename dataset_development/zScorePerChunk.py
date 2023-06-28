# A generic helper function to calculate z-score values per raster cell, per chunk
# A z-score refers to the number of standard deviations each data value is from the mean, with a z-score of zero indicating the exact mean
# Inputs: 
#   Single-band raster
#   Polygon shapefile containing raster chunk polygons

# Output:
#   Single-band with z-score values per cell

import arcpy


# Inputs
value_raster = arcpy.GetParameterAsText(0)
chunk_polygon = arcpy.GetParameterAsText(1)


# Environment Settings
desc = arcpy.Describe(value_raster)
spatialReference = desc.spatialReference
cellSize = desc.meanCellHeight
extent = desc.extent
workspace = desc.path

arcpy.env.snapRaster = value_raster
arcpy.env.outputCoordinateSystem = spatialReference
arcpy.env.extent = extent
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
zonalStatChunks = arcpy.sa.ZonalStatisticsAsTable(chunk_polygon, "Join_ID", value_raster, statsTable)
arcpy.JoinField_management(chunk_polygon, "Join_ID", statsTable, "Join_ID", ["MEAN", "STD"])


# Create MEAN and STD rasters to use in Raster Calculator
arcpy.AddMessage("Stats to raster")
mean_raster = arcpy.conversion.FeatureToRaster(chunk_polygon, "MEAN", workspace + r"\MEAN.tif")
std_raster = arcpy.conversion.FeatureToRaster(chunk_polygon, "STD", workspace + r"\STD.tif")


# Raster Calculator - calculate z-score
arcpy.AddMessage("Calculating Z-Scores per chunk")
# I cannot get this raster calculator function to work! Error below
# The parameter is incorrect. Parameter 'Rasters' is missing or invalid. Bind failed in function 'Raster Calculator Function'
# z_score_raster = arcpy.sa.RasterCalculator(rasters=[value_raster, mean_raster, std_raster], input_names=["VALUE", "MEAN", "STD"], expression="(VALUE-MEAN)/STD", extent_type="IntersectionOf")
outMinus = arcpy.sa.Minus(value_raster, mean_raster)
z_score_raster = arcpy.sa.Divide(outMinus, std_raster)
z_score_raster.save(workspace + r"\z_score_raster.tif")


# Add to map
arcpy.AddMessage("Adding to Map")
m = arcpy.mp.ArcGISProject("CURRENT").listMaps('*')[0]
m.addDataFromPath(arcpy.Describe(z_score_raster).catalogPath)