from flow_derivative_per_cell import flow_derivative


def main() -> None:

    flow_derivative(
        r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\model_testing\RAS\100yr\Depth (Max).Terrain.hydroDEM.tif",
        r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\model_testing\RAS\100yr\D _ V (Max).Terrain.hydroDEM.tif",
        r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\model_testing\RAS\100Plus\Depth (Max).Terrain.hydroDEM.tif",
        r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\model_testing\RAS\100Plus\D _ V (Max).Terrain.hydroDEM.tif",
        r"C:\Users\jacob.bates\OneDrive - WSP O365\2d_floodway_testing\future_floodway\testing_outputs\datasets"
    )

    return


if __name__ == "__main__":
    main()