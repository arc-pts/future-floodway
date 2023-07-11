from flow_derivative_per_cell import flow_derivative
from zscore_per_cell import zscore_per_cell
from os import PathLike


def main() -> None:

    get_flow_derivative(
        r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\testing_inputs\RAS\100yr\Depth (Max).Terrain.hydroDEM.tif",
        r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\testing_inputs\RAS\100yr\D _ V (Max).Terrain.hydroDEM.tif",
        r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\testing_inputs\RAS\100Plus\Depth (Max).Terrain.hydroDEM.tif",
        r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\testing_inputs\RAS\100Plus\D _ V (Max).Terrain.hydroDEM.tif",
        r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\testing_outputs\datasets"
    )

    # get_z_score_per_cell(
    #     r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\testing_inputs\composite_DxV_raster\composite_DxV.tif",
    #     4,
    #     r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\testing_outputs\datasets\z_scored_exclude_nodata.tif",
    #     False
    # )

    return


def get_flow_derivative(
    base_d_raster: PathLike, 
    base_dv_raster: PathLike, 
    plus_d_raster: PathLike, 
    plus_dv_raster: PathLike,
    out_directory: PathLike, 
    out_name: str = "flow_derivative.tif",
) -> PathLike:
    out_raster = flow_derivative(
        base_d_raster, 
        base_dv_raster, 
        plus_d_raster, 
        plus_dv_raster,
        out_directory, 
        out_name,
    )
    return out_raster


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


if __name__ == "__main__":
    main()