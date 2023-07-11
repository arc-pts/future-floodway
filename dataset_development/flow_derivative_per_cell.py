from os import PathLike, path
import arcpy
from arcpy.sa import *
arcpy.env.overwriteOutput = True


def flow_derivative(
    base_d_raster: PathLike, 
    base_dv_raster: PathLike, 
    plus_d_raster: PathLike, 
    plus_dv_raster: PathLike,
    out_directory: PathLike, 
    out_name: str = "flow_derivative.tif",
    ) -> PathLike:

    out_raster = path.join(out_directory, out_name)
    base_dv_raster = Raster(base_dv_raster); plus_dv_raster = Raster(plus_dv_raster); base_d_raster = Raster(base_d_raster); plus_d_raster = Raster(plus_d_raster)
    flow_diff = Con(base_dv_raster, (base_dv_raster-plus_dv_raster)/base_dv_raster)
    flow_diff.save(path.join(out_directory, "flow_diff.tif"))
    depth_diff = Con(base_d_raster, (base_d_raster-plus_d_raster)/base_d_raster)
    depth_diff.save(path.join(out_directory, "depth_diff.tif"))
    fd_rast = Con(base_dv_raster, (flow_diff/depth_diff))
    fd_rast.save(out_raster)

    return out_raster