#!/usr/bin/env python

import sys
import os
import getopt
from pathlib import Path
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import metpy.calc as mpcalc
from metpy.cbook import get_test_data
from metpy.interpolate import cross_section
from metpy.units import units

"""Extract data from a CVM netCDF file. Allow user to inspect the metadata, slice the data, plot, and 
    save the sliced data. The slicing can be performed along the existing coordinate planes or as an interpolated cross-sectional slice through gridded data.

    Call arguments:
        -h, --help: this message.
        -i, --input: [required] the input nefiletCDF filename.
        -v, --verbose [optional] provide informative information during the run.
"""

script = os.path.basename(sys.argv[0])
script = Path(script).stem
# Get the directory paths.
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

PROP_DIR = os.path.join(ROOT_DIR, "prop")
sys.path.append(PROP_DIR)
import shared_prop as prop
import slicer_prop as slicer_prop

LIB_DIR = os.path.join(ROOT_DIR, prop.LIB_DIR)
sys.path.append(LIB_DIR)
import shared_lib as lib
import slicer_lib as slicer_lib

dash = 10 * "-"

# Set up the logger.
logger = lib.get_logger()


def usage():
    logger.info(
        f"""
    Extract data from a CVM netCDF file. Interactively, user can inspect the metadata, slice the data, plot, and 
    save the sliced data. Currently slicing is only performed along the existing coordinate planes.

    Call arguments:
        -h, --help: this message.
        -v, --verbose: [optional] run in verbose mode.
        -i, --input: [required] the input nefiletCDF filename.
"""
    )


def main():
    cmap = prop.cmap
    interpolation_method = slicer_prop.interpolation_method
    xsection_steps = slicer_prop.steps

    # Capture the input parameters.
    try:
        argv = sys.argv[1:]
        opts, args = getopt.getopt(
            argv, "hvmi:o:", ["help", "verbose", "meta", "input=", "output="]
        )
    except getopt.GetoptError as err:
        # Print the error, and help information and exit:
        logger.error(err)
        usage()
        sys.exit(2)
    # Initialize the variables.
    input_file = None
    output_file = None
    verbose = False
    meta = None
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-v", "--verbose"):
            verbose = True
        elif o in ("-m", "--meta"):
            meta = True
        elif o in ("-i", "--input"):
            input_file = a
            if not os.path.isfile(input_file):
                logger.error(
                    f"[ERR] Invalid input netCDF file: [{input_file}]. File not found!"
                )
                sys.exit(2)
        elif o in ("-o", "--output"):
            output_file = a.strip()
        else:
            assert False, "unhandled option"

    # The input filename is required.
    if input_file is None:
        usage()
        logger.error("[ERR] missing -i or --input.")
        sys.exit(1)

    logger.info(f"[INFO] Working on input {input_file}")
    # Set the output filename the same as the input file.
    if output_file is None:
        output_file = f"{os.path.splitext(input_file)[0]}{prop.extension['geocsv']}"
    logger.info(f"[INFO] Will write the output to {output_file}")

    # Figure out the input's file type.
    file_type = lib.check_file_type(input_file)

    option = "start"
    coordinate_values = dict()

    # Convert netCDF to GeoCSV.
    logger.info("\n")
    if file_type["engine"] == "netcdf" and file_type["valid"] == True:
        if verbose:
            logger.info(
                f"\n\n{script}\n{dash}"
                f"\nTool for interactively extracting data from a CVM netCDF file."
                f"\nUser can inspect the metadata, slice the data, plot, and save"
                f"\nthe sliced data. Slicing can be performed along the existing"
                f"\ncoordinate planes or as an interpolated cross-sectional slice "
                f"through gridded data."
            )
        logger.info(f"\nLoaded {input_file} and it is a netCDF file")
        with xr.open_dataset(input_file, engine="netcdf4") as ds:
            data_var = list(ds.data_vars)
            dimensions = list(ds.dims)
            main_coordinates = list()
            aux_coordinates = list()
            coordinates = list(ds.coords)
            for var in coordinates:
                if var not in ds:
                    logger.error(f"[ERR] Missing data for variable {var}")
                    sys.exit(3)
                elif ds[var].data.size <= 1:
                    logger.error(f"[ERR] Missing data for variable {var}")
                    sys.exit(3)
                coordinate_values[var] = list(ds[var].data)

            for var in coordinates:
                if var in dimensions:
                    main_coordinates.append(var)
                else:
                    aux_coordinates.append(var)

            while option != "exit":
                if verbose:
                    logger.info(
                        f"\n\nThe available options are:\n\tmeta - to view file's metadata\n\tranges - to display ranges for variables\n\tsubset - to subset the data\n\t{dash}\n\thelp\n\texit "
                    )
                option = input("select option [meta, ranges, subset, help, exit]? ")

                # Done, back!
                if (
                    option.strip() == "exit"
                    or option.strip() == "back"
                    or not option.strip()
                ):
                    sys.exit()

                # Need help.
                elif option == "help":
                    usage()

                # Variables value ranges.
                elif option == "ranges":
                    logger.info("\nRanges:")
                    for var in ds.variables:
                        logger.info(
                            f"\t{var}: {np.nanmin(ds[var].data):0.2f} to  {np.nanmax(ds[var].data):0.2f} {ds[var].attrs['units']}"
                        )

                # Display the metadata.
                elif option == "meta":
                    logger.info(f"\n[Metadata] Global attributes:\n")
                    for row in ds.attrs:
                        if "geospatial" in row:
                            logger.info(f"\t{row}: {ds.attrs[row]}")

                    logger.info(f"\n[Metadata] Coordinate Variables:")
                    indent = 14
                    slicer_lib.display_var_meta(ds, dimensions, indent)
                    slicer_lib.display_var_meta(ds, data_var, indent, values=False)

                # Subset the model.
                elif option == "subset":
                    subset_type = "volume"
                    # Actions.
                    while subset_type != "back":
                        # f"\nSubset type [volume, slice, xsection, back, exit]? "
                        if verbose:
                            logger.info(
                                f"\n\nYou can subset the data as:\n\tvolume - a subvolume of data\n\tslice - a slice along a coordinate axis\n\txsection - a vertical slice in an arbitrary direction\n\t{dash}\n\tback - takes you to the previous step\n\texit "
                            )
                        subset_type = input(
                            f"\nselect [volume, slice, xsection, back, exit]? "
                        )
                        # Done, back!
                        if subset_type.strip() == "exit":
                            sys.exit()
                        elif subset_type.strip() == "back" or not option.strip():
                            break
                        if subset_type == "volume":
                            # Get the volume limits.
                            subset = ds.copy()
                            subset_limits = dict()
                            for dim in subset.dims:
                                subset_limits[dim] = (
                                    np.nanmin(subset[dim].data),
                                    np.nanmax(subset[dim].data),
                                )
                                if verbose:
                                    logger.info(
                                        f"\n\nEnter the {dim} range for the volume as minimum,maximum.\n\tPress return to accept the default values (full range) \n\t{dash}\n\tback - takes you to the previous step\n\texit "
                                    )
                                _limits = input(
                                    f"\n{dim} range [default values { subset_limits[dim]}, back, exit]? "
                                )

                                # Done, back!
                                if _limits.strip() == "exit":
                                    sys.exit()
                                elif _limits.strip() == "back":
                                    break
                                elif _limits and "," not in _limits:
                                    logger.warning(
                                        f"[WARN] Invalid limits, using the full range of {subset_limits[dim]}."
                                    )
                                    _limits = ""
                                if _limits:
                                    values = _limits.split(",")
                                    subset_limits[dim] = (
                                        float(values[0]),
                                        float(values[1]),
                                    )
                            subset_volume = slicer_lib.subsetter(ds, subset_limits)
                            # Subset the data.
                            logger.info(
                                f"\n\nThe selected volume information:\n\n{subset_volume}"
                            )
                            subset_action = "continue"
                            # Actions.
                            while subset_action != "back":
                                # Save the subset data.
                                subset_action == "save"
                                if verbose:
                                    logger.info(
                                        f"\n\nSave the data\n\tfilename -- filename for saving the file. Its extension determines the file format.\n\tfilename.csv to save as GeoCSV, filename.nc to output as netCDF\n\t{dash}\n\tback - takes you to the previous step\n\texit "
                                    )
                                filename = input(f"Output filename? [back, exit]: ")
                                # Done, back!
                                if filename.strip() == "exit":
                                    sys.exit()
                                elif filename.strip() == "back":
                                    break
                                if filename.endswith(".csv"):
                                    meta = lib.get_geocsv_metadata(subset_volume)
                                    with open(filename, "w") as fp:
                                        fp.write(meta)
                                    subset_volume.to_dataframe().to_csv(
                                        filename, mode="a"
                                    )
                                    logger.info(f"Saved as {filename}")
                                elif filename.endswith(".nc"):
                                    subset_volume.to_netcdf(filename, mode="w")
                                    logger.info(f"Saved as {filename}")
                                else:
                                    logger.error(f"[ERR] invalid file type: {filename}")
                        elif subset_type == "xsection":
                            # A cross-section of the model.
                            if verbose:
                                logger.info(
                                    f"Plotting an interpolated cross-sectional slice through gridded data\nPlease provide the cross-section limits"
                                )
                            logger.info("\nThe model coordinate ranges are:")

                            for var in ds.coords:
                                # Cross-sections use geographic coordinates.
                                if var in ["latitude", "longitude", "depth"]:
                                    logger.info(
                                        f"\t{var}: {np.nanmin(ds[var].data):0.2f} to  {np.nanmax(ds[var].data):0.2f} {ds[var].attrs['units']}"
                                    )
                            start = slicer_lib.get_point("start")
                            if start == "back":
                                break
                            end = slicer_lib.get_point("end")
                            if end == "back":
                                break
                            depth = slicer_lib.get_point("depth")
                            if depth == "break":
                                break
                            plot_data = ds.copy()
                            utm_zone = None
                            meta = ds.attrs
                            if "grid_mapping_name" not in meta:
                                logger.warning(
                                    f"[WARN] The 'grid_mapping_name' attribute not found. Assuming geographic coordinate system"
                                )
                                grid_mapping_name = "latitude_longitude"
                            else:
                                grid_mapping_name = meta["grid_mapping_name"]

                                logger.info(f"grid_mapping_name is {grid_mapping_name}")
                            if grid_mapping_name == "transverse_mercator":
                                if "utm_zone" in meta:
                                    utm_zone = meta["utm_zone"]
                                    logger.info(f"[INFO] UTM zone: {utm_zone}")
                                    projection = ccrs.UTM(utm_zone)
                                    # We use MetPy’s CF parsing to get the data ready for use, and squeeze down the size-one time dimension.
                                    plot_data = plot_data.metpy.assign_crs(
                                        projection.to_cf()
                                    )
                                    plot_data = plot_data.metpy.parse_cf().squeeze()
                                    plot_data = (
                                        plot_data.metpy.assign_latitude_longitude(
                                            force=True
                                        )
                                    )
                                else:
                                    logger.error(
                                        f"[ERR] The required attribute 'utm_zone' is missing for the grid_mapping_name of {grid_mapping_name}"
                                    )
                                    break
                            elif grid_mapping_name == "latitude_longitude":
                                projection = ccrs.Geodetic()
                                plot_data = plot_data.metpy.assign_crs(
                                    projection.to_cf()
                                )
                                plot_data = plot_data.metpy.parse_cf().squeeze()
                            else:
                                logger.warning(
                                    f"[ERR] The grid_mapping_name of {grid_mapping_name} is not supported!"
                                )
                                break

                            # Cross-section interpolation type.
                            interp_type = interpolation_method[0]
                            if verbose:
                                logger.info(
                                    f"\n\nSelect the interpolation method for creating the cross-section. The {', '.join(interpolation_method)} methods are available.\n\tPress return to select the default {interp_type} method\n\t{dash}\n\tback - takes you to the previous step\n\texit "
                                )
                            _interp_type = input(
                                f"\nInterpolation Method [{', '.join(interpolation_method+['back', 'exit'])}, default: {interp_type}]? "
                            )
                            if _interp_type.strip() == "exit":
                                sys.exit()
                            elif _interp_type.strip() == "back":
                                break
                            elif _interp_type.strip():
                                interp_type = _interp_type
                            elif interp_type not in (interpolation_method):
                                logger.warning(
                                    f"[WARN] Invalid interpolation method of {interp_type}. Will use {interpolation_method[0]}"
                                )
                                interp_type = interpolation_method[0]

                            # Steps in the cross-section.
                            steps = xsection_steps
                            if verbose:
                                logger.info(
                                    f"\n\nThe number of points along the geodesic between the start and the end point (including the end points) to use in the cross section.\n\tPress enter to select the default of 100.\n\t{dash}\n\tback - takes you to the previous step\n\texit "
                                )
                            steps = input(
                                f"\nNumber of points ['back', 'exit', default: {steps}]? "
                            )
                            if not steps.strip():
                                steps = xsection_steps
                            elif steps.strip() == "exit":
                                sys.exit()
                            elif steps.strip() == "back":
                                break
                            elif not steps.isnumeric():
                                logger.warning(
                                    f"[WARN] Invalid number of steps, expected an integer. Will use {xsection_steps}"
                                )
                                steps = xsection_steps
                            else:
                                steps = int(steps)

                                # Extract the cross-section.
                            xsection_data = cross_section(
                                plot_data,
                                start,
                                end,
                                steps=steps,
                                interp_type=interp_type,
                            )

                            # Iterate through the model variables and plot each cross-section.
                            plot_var = data_var[0]
                            while plot_var:

                                if len(data_var) > 1:
                                    if verbose:
                                        logger.info(
                                            f"\n\nThe model variable to plot.\n\tback - takes you to the previous step\n\texit "
                                        )
                                    plot_var = input(
                                        f"\nVariable to plot {data_var}, back, exit]: "
                                    )
                                # Done, back!
                                if plot_var.strip() == "exit":
                                    sys.exit()
                                elif plot_var.strip() == "back" or not option.strip():
                                    break
                                elif plot_var not in data_var:
                                    logger.error(f"[ERR] Invalid input {plot_var}")
                                    continue
                                pdata = xsection_data.copy()
                                pdata["depth"] = -pdata["depth"]

                                pdata[plot_var].plot.contourf(cmap=cmap)
                                plt.show()
                                break

                        elif subset_type == "slice":
                            # Slice the model.
                            slice_dir = "depth"
                            while slice_dir:
                                if verbose:
                                    logger.info(
                                        f"\nA slice cuts the model along one of the coordinate axis.\n\tdirection - direction of the slice, the coordinate to cut the model along\n\t{dash}\n\tback - takes you to the previous step\n\texit  "
                                    )
                                slice_dir = input(
                                    f"\ndirection [{', '.join(main_coordinates+['back', 'exit'])}]? "
                                )
                                # Done, back!
                                if slice_dir.strip() == "exit":
                                    sys.exit()
                                elif (
                                    slice_dir.strip() == "back"
                                    or slice_dir not in main_coordinates
                                    or not slice_dir.strip()
                                ):
                                    logger.error(f"[ERR] invalid variable {slice_dir}")
                                    break
                                # Explore what to do with the slice?
                                else:
                                    slice_value = input(
                                        f"\nslice {slice_dir} [{np.nanmin(ds[slice_dir].data)} to {np.nanmax(ds[slice_dir].data)}, back, exit]? "
                                    )
                                    # Exit.
                                    if slice_value.strip() == "exit":
                                        sys.exit()
                                    # A step back.
                                    elif (
                                        slice_value.strip() == "back"
                                        or not slice_value.strip()
                                    ):
                                        break
                                    # Actions.
                                    # while slice_value:
                                    closest_slice_value = lib.closest(
                                        coordinate_values[slice_dir],
                                        float(slice_value),
                                    )

                                    slice_dims = main_coordinates.copy()
                                    slice_dims.remove(slice_dir)
                                    slice_limits = dict()
                                    gmap_limits = dict()
                                    gmap_option = ""
                                    # Geographic coordinates?
                                    if (
                                        "latitude" in coordinates
                                        and "longitude" in coordinates
                                    ):
                                        # Make sure the geographical coordinates are independent of the slice direction.
                                        if (
                                            slice_dir not in ds["latitude"].dims
                                            and slice_dir not in ds["longitude"].dims
                                        ):
                                            gmap_option = ", gmap"
                                            for dim in ["longitude", "latitude"]:
                                                gmap_limits[dim] = (
                                                    np.nanmin(ds[dim].data),
                                                    np.nanmax(ds[dim].data),
                                                )
                                    # Get the slice limits.
                                    for dim in slice_dims:
                                        slice_limits[dim] = (
                                            np.nanmin(ds[dim].data),
                                            np.nanmax(ds[dim].data),
                                        )
                                        logger.info(
                                            f"\n\nSlicing at the closest {slice_dir} of {closest_slice_value}"
                                        )

                                        if verbose:
                                            logger.info(
                                                f"\nSelect slice limits in the {dim} direction.\n\tlimits - provide minimum,maximum or press enter to accept the full range default \n\t{dash}\n\tback - takes you to the previous step\n\texit  "
                                            )
                                        _limits = input(
                                            f"\n{dim} limits [default values { slice_limits[dim]}, back, exit]: "
                                        )
                                        # Done, back!
                                        if _limits.strip() == "exit":
                                            sys.exit()
                                        elif _limits.strip() == "back":
                                            break
                                        elif _limits.strip() and "," not in _limits:
                                            logger.warning(
                                                f"[WARN] Invalid limits, using the full range of {slice_limits[dim]}."
                                            )
                                            _limits = ""

                                        if _limits:
                                            values = _limits.split(",")
                                            slice_limits[dim] = (
                                                float(values[0]),
                                                float(values[1]),
                                            )
                                    # Slice the data.
                                    sliced_data = slicer_lib.slicer(
                                        ds,
                                        slice_dir,
                                        closest_slice_value,
                                        slice_limits,
                                    )
                                    logger.info(
                                        f"\n\nThe sliced data summary: \n{sliced_data}"
                                    )
                                    slice_action = "continue"
                                    # Actions.
                                    while slice_action != "back":

                                        if verbose:
                                            if gmap_option:
                                                logger.info(
                                                    f"\nWhat to do with the slice.\n\tplot2d - a 2D plot of {slice_dims}\n\tplot3d - a 3D plot of {slice_dims} and the model variable on the 3rd axis.\n\t\tThe plot is interactive and you can rotate it.\n\tgmap - a 2D plot of {slice_dims} in geographical coordinate system\n\tcmap - change the color map for the plots\n\tsave - save the slice data\n\t{dash} \n\tback - takes you to the previous step\n\texit  "
                                                )
                                            else:
                                                logger.info(
                                                    f"\nWhat to do with the slice.\n\tplot2d - a 2D plot of {slice_dims}\n\tplot3d - a 3D plot of {slice_dims} and the model variable on the 3rd axis.\n\tThe plot is interactive and you can rotate it.\n\tcmap - change the color map for the plots \n\tsave - save the slice data\n\t{dash}\n\tback - takes you to the previous step\n\texit  "
                                                )
                                        slice_action = input(
                                            f"\nAction [plot2d, plot3d{gmap_option}, cmap, save, back, exit]: "
                                        )
                                        # Done, back!
                                        if slice_action.strip() == "exit":
                                            sys.exit()
                                        elif (
                                            slice_action.strip() == "back"
                                            or not option.strip()
                                        ):
                                            break
                                        # 2D plot.
                                        if slice_action == "plot2d":
                                            # Plot each model variable.
                                            plot_var = data_var[0]
                                            while plot_var:
                                                if len(data_var) > 1:
                                                    plot_var = input(
                                                        f"\nVariable to plot {data_var}, back, exit]: "
                                                    )
                                                # Done, back!
                                                if plot_var.strip() == "exit":
                                                    sys.exit()
                                                elif (
                                                    plot_var.strip() == "back"
                                                    or not option.strip()
                                                ):
                                                    break
                                                elif plot_var not in data_var:
                                                    logger.error(
                                                        f"[ERR] Invalid input {plot_var}"
                                                    )
                                                    continue
                                                plot_data = sliced_data[plot_var].copy()
                                                # Set the depth axis (if exists) direction.
                                                if "depth" in plot_data.dims:
                                                    if (
                                                        "geospatial_vertical_positive"
                                                        in sliced_data.attrs
                                                    ):
                                                        if (
                                                            sliced_data.attrs[
                                                                "geospatial_vertical_positive"
                                                            ]
                                                            == "down"
                                                        ):
                                                            plot_data["depth"] = (
                                                                -plot_data["depth"]
                                                            )

                                                plot_data.plot(cmap=cmap)
                                                plt.show()
                                                if len(data_var) <= 1:
                                                    plot_var = "back"
                                        # 3D plot.
                                        elif slice_action == "plot3d":
                                            plot_var = data_var[0]
                                            while plot_var:
                                                if len(data_var) > 1:
                                                    plot_var = input(
                                                        f"\nVariable to plot {data_var}, back, exit]: "
                                                    )
                                                # Done, back!
                                                if plot_var.strip() == "exit":
                                                    sys.exit()
                                                elif (
                                                    plot_var.strip() == "back"
                                                    or not option.strip()
                                                ):
                                                    break
                                                elif plot_var not in data_var:
                                                    logger.error(
                                                        f"[ERR] Invalid input {plot_var}"
                                                    )
                                                    continue

                                                sliced_data[plot_var].plot.surface(
                                                    cmap=cmap,
                                                )
                                                plt.show()
                                                if len(data_var) <= 1:
                                                    plot_var = "back"
                                        # Geographic map.
                                        elif slice_action == "gmap" and gmap_option:
                                            plot_var = data_var[0]
                                            while plot_var:
                                                if len(data_var) > 1:
                                                    plot_var = input(
                                                        f"\nVariable to plot {data_var}, back, exit]: "
                                                    )
                                                # Done, back!
                                                if plot_var.strip() == "exit":
                                                    sys.exit()
                                                elif (
                                                    plot_var.strip() == "back"
                                                    or not option.strip()
                                                ):
                                                    break
                                                elif plot_var not in data_var:
                                                    logger.error(
                                                        f"[ERR] Invalid input {plot_var}"
                                                    )
                                                    continue
                                                slicer_lib.gmap(
                                                    plot_var,
                                                    cmap,
                                                    gmap_limits,
                                                    sliced_data,
                                                )
                                                if len(data_var) <= 1:
                                                    plot_var = "back"
                                        # Colormap.
                                        elif slice_action == "cmap":
                                            logger.info(f"\n[cmaps] {plt.colormaps()}")
                                            cmap = input(
                                                f"\nSelect a color map for the plot [default cmap {cmap}, back, exit]: "
                                            )
                                        # Save the slice data.
                                        elif slice_action == "save":
                                            if verbose:
                                                logger.info(
                                                    f"\n\nSave the data\n\tfilename -- filename for saving the file. Its extension determines the file format.\n\t\tfilename.csv to save as GeoCSV, filename.nc to output as netCDF\n\t{dash}\n\tback - takes you to the previous step\n\texit "
                                                )
                                            filename = input(
                                                f"Output filename? [back, exit]: "
                                            )

                                            # Done, back!
                                            if filename.strip() == "exit":
                                                sys.exit()
                                            elif (
                                                filename.strip() == "back"
                                                or not option.strip()
                                            ):
                                                break
                                            elif filename.endswith(".csv"):
                                                meta = lib.get_geocsv_metadata(
                                                    sliced_data
                                                )
                                                with open(filename, "w") as fp:
                                                    fp.write(meta)
                                                sliced_data.to_dataframe().to_csv(
                                                    filename, mode="a"
                                                )
                                            elif filename.endswith(".nc"):
                                                sliced_data.to_netcdf(
                                                    filename, mode="w"
                                                )
                                            else:
                                                logger.error(
                                                    f"[ERR] invalid file type: {filename}"
                                                )
                                # What does the user want to do next?
                                slice_value = input(
                                    f"\nSlice at {slice_dir} of (between {np.nanmin(ds[slice_dir].data)} and {np.nanmax(ds[slice_dir].data)}) [back, exit]: "
                                )
                                # Done, back!
                                if slice_value.strip() == "exit":
                                    sys.exit()
                                elif (
                                    slice_value.strip() == "back" or not option.strip()
                                ):
                                    break

    else:
        logger.error(f"[ERR] {input_file} is not a valid netCDF file")
        sys.exit(1)


if __name__ == "__main__":
    main()
