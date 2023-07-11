from os import path, PathLike
from ras_mesh_data import fetch_mesh_data
from spatial_utils import merge_shapes, copy_shape


def main() -> None:
    ras_plan_hdf_file = r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\model_testing\lbr1.p04.hdf"
    scoped_stream_network = r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\model_testing\streams_scope.shp"
    out_directory = r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\model_testing"

    mesh_areas = fetch_mesh_data(ras_plan_hdf_file)
    for me in mesh_areas:
        me.get_trace_sources(scoped_stream_network)
        # me.trace_floodplains()
        me.trace_neighbor_chunks()
        me.neighbor_chunks_to_shapefile(scoped_stream_network, out_directory)
        # me.cells_to_shapefile(me.floodplain_cell_ids, scoped_stream_network, out_directory, "floodplain_cells")
        # me.trace_neighbor_chunks(scoped_stream_network)
        # me.cells_to_shapefile(me.neighbor_chunks[22016], scoped_stream_network, out_directory, "neighbor_chunk_22016")
        # me.cells_to_shapefile(me.neighbor_chunks[22728], scoped_stream_network, out_directory, "neighbor_chunk_22728")
        # me.trace_backwater_chunks(scoped_stream_network)
        # me.cells_to_shapefile(me.backwater_chunks[22016], scoped_stream_network, out_directory, "backwater_chunk_22016")
        # me.cells_to_shapefile(me.backwater_chunks[22728], scoped_stream_network, out_directory, "backwater_chunk_22728")

    return


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


if __name__ == "__main__":
    main()