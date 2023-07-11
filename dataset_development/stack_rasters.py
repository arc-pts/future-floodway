from os import PathLike, path
import arcpy


def stack_rasters(
    raster_datasets: list[PathLike], 
    out_directory: PathLike = arcpy.env.workspace, 
    out_name: str = "stacked_raster.tif",
    write_metadata: bool = True
    ) -> PathLike:
    """
    combines multiple single-band rasters into one multi-band raster.
    optionally write a metadata txt file that described the bands within the output raster.
    """

    out_raster = path.join(out_directory, out_name)
    arcpy.CompositeBands_management(";".join(raster_datasets), out_raster)

    if write_metadata:
        with open(path.join(out_directory, f"{out_name.split('.')[0]}_metadata.txt"), "w") as txt:
            txt.write(f"this raster contains {len(raster_datasets)} bands made up of the following datasets in the order in which they are named:\n")
            [txt.write(rd+"\n") for rd in raster_datasets]

    return out_raster