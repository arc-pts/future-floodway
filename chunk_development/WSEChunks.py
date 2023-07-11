"""
Script to create polygon chunks based on WSE intervals
Inputs: 
  Single-band WSE raster (Horizontal render)
  scoped_stream_network - streams used to assign chunking extents
  Output folder - Optional. Will save in WSE raster directory if not provided
  out_name - name of output chunk dataset
  Interval (1ft, 2ft, etc). Good up to 2 decimal places

  Potential future input? - Scoped streamlines: Used to clip out polygons

Output:
  Polygon shapefile containing interval chunks. Holes filled in (except at polygon boundaries)
"""
import arcpy
import math
from os import PathLike, path


def WSEChunks(
  WSE_raster: PathLike, 
  scoped_stream_network: PathLike,
  output_folder: PathLike, 
  out_name: str, 
  interval: int = 2
) -> PathLike:
  # Environment Settings
  desc = arcpy.Describe(WSE_raster)
  spatialReference = desc.spatialReference
  cellSize = desc.meanCellHeight
  extent = desc.extent

  if output_folder:
    workspace = output_folder
  else:
    workspace = desc.path

  arcpy.AddMessage("Output will be saved in {0}".format(workspace))

  arcpy.env.snapRaster = WSE_raster
  arcpy.env.outputCoordinateSystem = spatialReference
  arcpy.env.extent = extent
  arcpy.env.cellSize = cellSize
  arcpy.env.workspace = workspace
  arcpy.env.overwriteOutput = True

  minWSE = math.floor(float(arcpy.GetRasterProperties_management(WSE_raster, "MINIMUM")[0]))
  maxWSE = math.ceil(float(arcpy.GetRasterProperties_management(WSE_raster, "MAXIMUM")[0]))

  # Reclassify. Set range to lower bound
  arcpy.AddMessage("Reclassifying raster")
  reclassify_range = []
  workingWSE = minWSE

  while workingWSE < maxWSE:
      reclassify_range.append([workingWSE, workingWSE+interval, int(workingWSE*100)])  # Remap*100 since it needs to be an integer and we want at least a few decimal places
      workingWSE += interval

  # arcpy.AddMessage(reclassify_range)
  WSE_Reclassify = arcpy.sa.Reclassify(WSE_raster, "Value", arcpy.sa.RemapRange(reclassify_range))
  # WSE_Reclassify.save(workspace + r"\WSE_Reclassify.tif")

  # Raster to Polygon and fill in holes
  arcpy.AddMessage("Raster to Polygon")
  WSE_chunks = arcpy.conversion.RasterToPolygon(WSE_Reclassify, path.join(workspace, "WSE_chunks.shp"))
  WSE_chunks_eliminate = arcpy.EliminatePolygonPart_management(WSE_chunks, path.join(workspace, out_name), "PERCENT", part_area_percent=90)

  # Add field with WSE range
  arcpy.AddField_management(WSE_chunks_eliminate, "WSE_Low", "FLOAT", field_alias='WSE Lower Bound')
  arcpy.AddField_management(WSE_chunks_eliminate, "WSE_High", "FLOAT", field_alias='WSE Upper Bound')

  with arcpy.da.UpdateCursor(WSE_chunks_eliminate, ['gridcode','WSE_Low', 'WSE_High']) as cursor:
      for row in cursor:
          row[1] = row[0]/100 
          row[2] = row[0]/100 + interval
          cursor.updateRow(row)
          
  # Clip polygons to scoped stream - future feature?

  # # Delete Intermediate Data
  arcpy.Delete_management(WSE_chunks)

  arcpy.AddMessage("Output saved to {0}".format(WSE_chunks_eliminate))
  return arcpy.Describe(WSE_chunks_eliminate).catalogPath


def addToMap(layer):
  arcpy.AddMessage("Adding to Map")
  m = arcpy.mp.ArcGISProject("CURRENT").listMaps('*')[0]
  m.addDataFromPath(arcpy.Describe(layer).catalogPath)


if __name__ == "__main__":

    WSE_raster = arcpy.GetParameterAsText(0)
    interval = arcpy.GetParameter(1)
    output_folder = arcpy.GetParameterAsText(2)

    wse_chunk_polygon = WSEChunks(WSE_raster, interval, output_folder)
    addToMap(wse_chunk_polygon)