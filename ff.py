from os import PathLike
import numpy as np
import datetime
from chunk_development.ras_mesh_data import fetch_mesh_data
from chunk_development.spatial_utils import merge_shapes, copy_shape
from chunk_development.WSEChunks import WSEChunks
from chunk_development.spatial_utils import *
from dataset_development.flow_derivative_per_cell import flow_derivative
from dataset_development.zscore_per_cell import zscore_per_cell
from dataset_development.cov_per_cell import cov_per_cell


def main() -> None:

    ras_plan_hdf_file = r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\testing_inputs\test_hdf_file\lbr1.p04.hdf"
    # chunk_layer = r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\testing_outputs\chunks\neighbor_chunks.shp"
    scoped_stream_network = r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\testing_inputs\streams_layer\streams_scope.shp"
    DV2_raster = r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\testing_inputs\RAS\100yr\D _ V^2 (Max).Terrain.hydroDEM.tif"
    percentiles = [30,40,50,60,70,80,90]
    out_directory = r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\testing_outputs\datasets\percentile_DV2_per_chunk"

    me = fetch_mesh_data(ras_plan_hdf_file)[0]
    cell_pnts_shp = me.cell_points_to_shapefile(scoped_stream_network, out_directory)
    me.trace_floodplains(scoped_stream_network)
    fp_mask = me.floodplain_cells_to_mask_shape(scoped_stream_network, out_directory)
    me.trace_neighbor_chunks(scoped_stream_network)
    chunk_polys = me.neighbor_chunks_to_polygons(scoped_stream_network)

    # chunks = arcpy.CreateFeatureclass_management(out_directory, "chunks", "polygon", spatial_reference=arcpy.Describe(scoped_stream_network).spatialReference)
    # with arcpy.da.InsertCursor(chunks, ["shape@"]) as ic:
    #     for chunk_shape in chunk_polys.values():
    #         ic.insertRow(chunk_shape)

    # arcpy.CopyFeatures_management(chunk_polys.values(), path.join(out_directory, "chunks"))

    start = datetime.datetime.now(); print(f"started getting raster percentiles for chunks at {start}...")
    rast = arcpy.Raster(DV2_raster)
    chunk_percentile_vals = dict() # {chunk_id : [percentile values]}
    cnt = 0
    total = len(chunk_polys.keys())
    for chunk_id, chunk_poly in chunk_polys.items():
        # cnt+=1; print(f"    getting raster percentiles for chunk {cnt} of {total}")
        chunk_rast = ExtractByMask(rast, chunk_poly)
        chunk_percentile_vals[chunk_id] = np.percentile(chunk_rast.read(), percentiles)
    print(list(chunk_percentile_vals.values())[0])
    print(f"completed raster percentiles for chunks in {datetime.datetime.now() - start}.")

    start = datetime.datetime.now(); print(f"started getting cell percentiles at {start}...")
    cell_percentile_vals = dict() # {cell_id : [percentile values]}
    cnt = 0
    total = len(me.cell_residence_neighbor_chunks.keys())
    for cell_id, chunk_ids in me.cell_residence_neighbor_chunks.items():
        # cnt+=1; print(f"    getting cell percentiles for cell {cnt} of {total}")
        cell_percentile_vals[cell_id] = np.stack([chunk_percentile_vals[chunk_id] for chunk_id in chunk_ids], axis=1).max(axis=1)
    print(list(cell_percentile_vals.values())[0])
    print(f"completed cell percentiles in {datetime.datetime.now() - start}.")

    start = datetime.datetime.now(); print(f"started creating tin points at {start}...")
    p_fields = list(f"per_{p}" for p in percentiles)
    [arcpy.AddField_management(cell_pnts_shp, pf, "FLOAT") for pf in p_fields]
    query = f"cell_id IN ({str(list(me.cell_residence_neighbor_chunks.keys())).strip('[]')})"
    cell_pnts_shp_cleaned = arcpy.SelectLayerByAttribute_management(cell_pnts_shp, where_clause=query)
    with arcpy.da.UpdateCursor(cell_pnts_shp_cleaned, ["cell_id"]+p_fields) as uc:
        for row in uc:
            row[1:] = cell_percentile_vals[row[0]]
            uc.updateRow(row)
    print(f"completed tin points in {datetime.datetime.now() - start}.")
    arcpy.CopyFeatures_management(cell_pnts_shp_cleaned, path.join(out_directory, "cell_pnts_shp_cleaned"))


    # """
    # here need to get the cell points and the floodplain cell polygon mask and then grab values from percentile raster and create tin raster -> stats_rast
    # """
    # stats_rast = None

    # floodway_from_rasters(
    #     DV2_raster,
    #     stats_rast,
    #     out_directory,
    #     out_name="floodway.tif"
    # )

    return


class chunks(object):

    def __init__(self) -> None:
        pass

    @staticmethod
    def from_backwater_cells(
        ras_plan_hdf_file: PathLike, 
        scoped_stream_network: PathLike,
        out_directory: PathLike, 
        out_name: str = "backwater_chunks"
    ) -> PathLike:
        mesh_areas = fetch_mesh_data(ras_plan_hdf_file)
        chunk_shapes = list()
        for me in mesh_areas:
            me.trace_backwater_chunks(scoped_stream_network)
            chunk_shapes.append(
                me.backwater_chunks_to_shapefile(
                    scoped_stream_network,
                    r"C:\temp",
                    out_name+f"_{me.name}"
                )
            )
        if len(chunk_shapes) > 1:
            chunk_shape = merge_shapes(
                chunk_shapes,
                scoped_stream_network,
                out_directory,
                out_name
            )
        else:
            chunk_shape = copy_shape(
                chunk_shapes[0],
                out_directory,
                out_name
            )
        return chunk_shape

    @staticmethod
    def from_neighbor_cells(
        ras_plan_hdf_file: PathLike, 
        scoped_stream_network: PathLike,
        out_directory: PathLike, 
        out_name: str = "neighbor_chunks"
    ) -> PathLike:
        mesh_areas = fetch_mesh_data(ras_plan_hdf_file)
        chunk_shapes = list()
        for me in mesh_areas:
            me.trace_neighbor_chunks(scoped_stream_network)
            chunk_shapes.append(
                me.neighbor_chunks_to_shapefile(
                    scoped_stream_network,
                    r"C:\temp",
                    out_name+f"_{me.name}"
                )
            )
        if len(chunk_shapes) > 1:
            chunk_shape = merge_shapes(
                chunk_shapes,
                scoped_stream_network,
                out_directory,
                out_name
            )
        else:
            chunk_shape = copy_shape(
                chunk_shapes[0],
                out_directory,
                out_name
            )
        return chunk_shape

    @staticmethod
    def from_wse_raster(
        WSE_raster: PathLike, 
        scoped_stream_network: PathLike,
        output_folder: PathLike, 
        out_name: str, 
        interval: int = 2
    ) -> PathLike:
        chunk_shape = WSEChunks(
            WSE_raster, 
            scoped_stream_network,
            output_folder, 
            out_name, 
            interval
        )
        return chunk_shape


class datasets(object):

    def __init__(self) -> None:
        pass

    @staticmethod
    def get_flow_derivative(
        base_d_raster: PathLike,
        base_v_raster: PathLike,
        base_dv_raster: PathLike, 
        plus_d_raster: PathLike, 
        plus_v_raster: PathLike, 
        plus_dv_raster: PathLike,
        out_directory: PathLike
    ) -> PathLike:
        out_raster = flow_derivative(
            base_d_raster,
            base_v_raster,
            base_dv_raster, 
            plus_d_raster,
            plus_v_raster,
            plus_dv_raster,
            out_directory
        )
        return out_raster

    @staticmethod
    def get_z_score_per_cell(
        multiband_raster: PathLike,
        target_raster_band: int,
        output_raster: PathLike,
        consider_nodata: bool = False
    ) -> PathLike:
        out_raster = zscore_per_cell(
            multiband_raster,
            target_raster_band,
            output_raster,
            consider_nodata
        )
        return out_raster

    @staticmethod
    def get_cov_per_cell(
        multiband_raster: PathLike,
        output_raster: PathLike,
        consider_nodata: bool = False
    ) -> PathLike:
        out_raster = cov_per_cell(
            multiband_raster,
            output_raster,
            consider_nodata
        )
        return out_raster



if __name__ == "__main__":
    main()