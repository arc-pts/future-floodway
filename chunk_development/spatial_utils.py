from os import path, PathLike
import numpy as np
from scipy.spatial import KDTree
from cell import cell
from face import face
import arcpy
arcpy.env.overwriteOutput = True
arcpy.env.addOutputsToMap = 0
arcpy.env.workspace = arcpy.env.scratchGDB
scratch_folder = arcpy.env.scratchFolder


class BadProjectionError(Exception):
    pass


def check_projection(scoped_stream_network: PathLike) -> None:
    sr = arcpy.Describe(scoped_stream_network).spatialReference
    if not sr.PCSCode or sr.PCSCode == 0:
        raise BadProjectionError("input stream network is not georeferenced in a state plane coordinate system. please reproject to state plane and try again.")


def streams_to_points_array(
    stream_network_layer: PathLike,
    mesh_perimeter_shape: PathLike,
    point_interval_ft: int = 100
) -> np.ndarray:
    start_pnts = arcpy.GeneratePointsAlongLines_management(
        stream_network_layer, 
        "start_pnts", 
        Distance=f"{point_interval_ft} Feet", 
        Include_End_Points="END_POINTS"
    )
    start_pnts_clipped = arcpy.SelectLayerByLocation_management(start_pnts, "INTERSECT", mesh_perimeter_shape)
    return np.array(tuple(row[0] for row in arcpy.da.SearchCursor(start_pnts_clipped, "Shape@XY")))


def perimeter_to_shapefile(
    mesh_perimeter_coordinates: np.ndarray, # np.ndarray[np.ndarray[float]]
    projection_file: PathLike, 
    out_directory: PathLike = arcpy.env.workspace, 
    out_name: str = "mesh_perimeter"
) -> PathLike:
    mesh_perimeter_shp = arcpy.CreateFeatureclass_management(
        out_directory, 
        out_name, 
        "POLYGON", 
        spatial_reference=arcpy.Describe(projection_file).spatialReference
    )
    with arcpy.da.InsertCursor(mesh_perimeter_shp, ["SHAPE@"]) as cur:
        cur.insertRow([arcpy.Polygon(arcpy.Array([arcpy.Point(*point) for point in mesh_perimeter_coordinates]))])
    return arcpy.Describe(mesh_perimeter_shp).catalogPath


def chunks_to_shape(
    chunk_dict: dict[int:set[cell]], 
    projection_file: PathLike, 
    out_directory: PathLike, 
    out_name: str = "chunks"
) -> PathLike:
    return PathLike


def copy_shape(
    source: PathLike, 
    out_directory: PathLike, 
    out_name: str = "out_shape"
) -> PathLike:
    return PathLike


def cells_to_shape(
    cell_faces: dict[int:list[face]], 
    projection_file: PathLike, 
    out_directory: PathLike = arcpy.env.workspace, 
    out_name: str = "cells"
) -> PathLike:
    """
    cell_faces: {cell_id : [face objects]}
    """
    cell_faces_shp = arcpy.CreateFeatureclass_management(
        arcpy.env.workspace, 
        "cell_faces_shp", 
        "POLYLINE", 
        spatial_reference=arcpy.Describe(projection_file).spatialReference
    )
    with arcpy.da.InsertCursor(cell_faces_shp, ["SHAPE@"]) as cur:
        for cell_id, faces in cell_faces.items():
            for face_obj in faces:
                cur.insertRow([arcpy.Polyline(arcpy.Array([arcpy.Point(*pnt) for pnt in face_obj.coordinates]))])
    cell_shape = arcpy.FeatureToPolygon_management(cell_faces_shp, path.join(out_directory, out_name))
    return arcpy.Describe(cell_shape).catalogPath


def merge_shapes(
    shapes: list[PathLike], 
    projection_file: PathLike, 
    out_directory: PathLike, 
    out_name: str = "chunks"
) -> PathLike:
    return PathLike


def point_spatial_index(input_points: np.ndarray = None, xy_array: np.ndarray = None):
    """
    input_points - this is an array of point xy coord arrays array(array(x,y), array(x,y), ...) that need to be associated with their nearest neighbor
    xy_array - an array of point xy coord arrays array(array(x,y), array(x,y), ...) of potential nearest neighbors to the input_points
    """
    mytree = KDTree(xy_array)
    nearest_neighbor_dists, nearest_neighbor_indexes = mytree.query(input_points)
    return nearest_neighbor_indexes


def scratch_me() -> None:
    scratch_db = arcpy.env.scratchGDB
    dir, filename = path.split(scratch_db)
    try:
        if arcpy.Exists(scratch_db):
            try: arcpy.Delete_management(scratch_db)
            except Exception as e: arcpy.AddMessage(f"could not rewrite scratch database: {e}")
        arcpy.CreateFileGDB_management(dir, filename)
    except: pass


if __name__ == "__main__":
    pass