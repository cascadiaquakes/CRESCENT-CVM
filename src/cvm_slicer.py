#!/usr/bin/env python

import sys
import os
import getopt
import traceback
from pathlib import Path
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pyproj import Proj, Geod

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
    Extract data from a CVM netCDF file. Interactively, users can inspect the metadata, slice the data, plot, 
    and save the sliced data. Currently, slicing is performed only along the existing coordinate planes.

    Call arguments:
        -h, --help: this message.
        -v, --verbose: [optional] run in verbose mode.
        -i, --input: [required] the input nefiletCDF filename.
"""
    )

def display_range(ds):
    """Output ranges for the given dataset."""
    logger.info("\nRanges:")
    for var in ds.variables:
        logger.info(
            f"\t{var}: {np.nanmin(ds[var].data):0.2f} to  {np.nanmax(ds[var].data):0.2f} {ds[var].attrs['units']}"
            )
    logger.info("\n")

def display_metadata(ds):
    """Output metadata for the given dataset."""
    logger.info(f"\n[Metadata] Global attributes:\n")
    for row in ds.attrs:
        if "geospatial" in row:
            logger.info(f"\t{row}: {ds.attrs[row]}")

    logger.info(f"\n[Metadata] Coordinate Variables:")
    indent = 14
    slicer_lib.display_var_meta(ds, list(ds.dims), indent)
    slicer_lib.display_var_meta(ds, list(ds.data_vars), indent, values=False)
    logger.info("\n")

def output_messages(messages):
    """Output messages stored in a message list."""
    if messages:
        new_line = "\n"
        logger.info(f"{new_line}{new_line}NOTES:{new_line}{dash}{new_line}{new_line.join(messages)}")
    logger.info("\n")
    return list()     

def output_data(subset_volume, filename, messages):
    valid_file_extensions = [".csv", ".gcsv", ".nc"]
    
    if filename.endswith(".gcsv"):
        meta = lib.get_geocsv_metadata_from_ds(
            subset_volume
        )
        meta = f"# dataset: GeoCSV2.0\n# delimiter: ,\n{meta}"
        with open(filename, "w") as fp:
            fp.write(meta)
        subset_volume.to_dataframe().to_csv(
            filename, mode="a"
        )
        messages.append(f"[INFO] Saved as {filename}")
    elif filename.endswith(".csv"):
        meta = lib.get_geocsv_metadata_from_ds(
            subset_volume
        )
        with open(filename, "w") as fp:
            subset_volume.to_dataframe().to_csv(
            filename, mode="w"
        )
        messages.append(f"[INFO] Saved as {filename}")
    elif filename.endswith(".nc"):
        subset_volume.to_netcdf(filename, mode="w")
        messages.append(f"[INFO] Saved as {filename}")
    else:
        messages.append(f"[ERR] invalid file type: {filename} (must have one of the extensions: {valid_file_extensions})")
    return messages

def main():
    vmin = None
    vmax = None
    messages = list()
    cmap = prop.cmap
    interpolation_method = slicer_prop.interpolation_method
    xsection_steps = slicer_prop.steps
    vertical_exaggeration = 0


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
        else:
            assert False, "unhandled option"

    # The input filename is required.
    if input_file is None:
        usage()
        logger.error("[ERR] missing -i or --input.")
        sys.exit(1)

    logger.info(f"[INFO] Working on input {input_file}")

    # Figure out the input's file type.
    file_type = lib.check_file_type(input_file)

    option = "start"
    coordinate_values = dict()

    # Convert netCDF to GeoCSV.
    logger.info("\n")
    if file_type["engine"] == "netcdf" and file_type["valid"] == True:
        if verbose:
            logger.info(
                f"\n\n{script}\n{dash}\n"
                f"Tool for interactively extracting data from a CVM netCDF file. Users can inspect the metadata,\n"
                f"slice the data, plot, and save the sliced data. Slicing can be performed along the existing coordinate\n"
                f"planes or as an interpolated cross-sectional slice through gridded data.\n\n"
            )
        messages.append(f"[INFO] Loaded {input_file} and it is a netCDF file")
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
                        f"\nThe available options are:\n\tmeta - to view file's metadata\n\trange - to display value ranges for variables\n\tsubset - to subset the data\n\t{dash}\n\thelp\n\texit "
                    )
                messages = output_messages(messages)
                option = input("select option [meta, range, subset, help, exit]? ")

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

                # Variables value range.
                elif option == "range":
                    display_range(ds)

                # Display the metadata.
                elif option == "meta":
                    display_metadata(ds)

                # Subset the model.
                elif option == "subset":
                    subset_type = "volume"
                    # Actions.
                    while subset_type != "back":
                        # f"\nSubset type [volume, slice, xsection, back, exit]? "
                        if verbose:
                            logger.info(
                                f"\nYou can subset the data as a:\n\tvolume - a subvolume of data\n\tslice - a slice along a coordinate axes\n\txsection - a vertical slice in an arbitrary direction\n\t{dash}\n\tback - takes you to the previous step\n\texit "
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
                                        f"\nEnter the {dim} range for the volume as minimum,maximum.\n\tPress return to accept the default values (full range) \n\t{dash}\n\tback - takes you to the previous step\n\texit "
                                    )
                                messages = output_messages(messages)
                                _limits = input(
                                    f"{dim} range [default values { subset_limits[dim]}, back, exit]? "
                                )

                                # Done, back!
                                if _limits.strip() == "exit":
                                    sys.exit()
                                elif _limits.strip() == "back":
                                    break
                                elif _limits and "," not in _limits:
                                    messages.append(f"[WARN] Invalid limits for {dim} ({_limits}), using the full range of {subset_limits[dim]}.")
                                    logger.warning(messages[-1])
                                    _limits = ""
                                if _limits:
                                    values = _limits.split(",")
                                    if lib.is_numeric(values[0]) and lib.is_numeric(values[1]):
                                        values[0] = float(values[0].strip())
                                        values[1] = float(values[1].strip())
                                        if lib.is_in_range(values[0], dim)  and lib.is_in_range(values[1], dim): 
                                            subset_limits[dim] = (
                                                min(values), max(values)
                                            )
                                            messages.append(f"[INFO] {dim} limits: {values}.")
                                        else:
                                            messages.append(f"[WARN] Invalid limits of {_limits} for {dim}, using the full range of {subset_limits[dim]}.")
                                            logger.warning(messages[-1])
                                            _limits = ""
                                    else:
                                        messages.append(f"[WARN] Invalid limits of {_limits} for {dim} (both limits must be provided as numbers). Will be using the full range of {subset_limits[dim]}.")
                                        logger.warning(messages[-1])
                                        _limits = ""
                                else:
                                    messages.append(f"[Info] {dim} will be set to full range: {subset_limits[dim]}.")
                                    logger.warning(messages[-1])
                                    _limits = ""
                            if _limits.strip() == "back":
                                break
                            subset_volume = slicer_lib.subsetter(ds, subset_limits)
                            # Subset the data.
                            messages.append(
                                f"\nThe selected volume information:\n{subset_volume}"
                            )
                            subset_action = "continue"
                            # Actions.
                            messages = output_messages(messages)
                            while subset_action != "back":
                                # Save the subset data.
                                subset_action == "save"
                                if verbose:
                                    logger.info(
                                        f"\nSave the data\n\tFilename -- The name of the file to save to. The extension determines the file format. Use to save as CSV, .csv .gcsv to save as GeoCSV, and .nc to save as netCDF.\n\t{dash}\n\tback - takes you to the previous step\n\texit "
                                    )
                                
                                messages = output_messages(messages)
                                filename = input(f"Output filename? [back, exit]: ")
                                # Done, back!
                                if filename.strip() == "exit":
                                    sys.exit()
                                elif filename.strip() == "back":
                                    break
                                messages = output_data(subset_volume, filename, messages)
                        elif subset_type == "xsection":
                            # A cross-section of the model.
                            if verbose:
                                logger.info(
                                    f"[INFO] Plotting an interpolated cross-sectional slice through gridded data\nPlease provide the cross-section limits"
                                )
                            logger.info(
                                "\nThe model coordinate ranges for the cross-section are:"
                            )

                            for var in ds.coords:
                                # Cross-sections use geographic coordinates.
                                if var in ["latitude", "longitude", "depth"]:
                                    logger.info(
                                        f"\t{var}: {np.nanmin(ds[var].data):0.2f} to  {np.nanmax(ds[var].data):0.2f} {ds[var].attrs['units']}"
                                    )
                            
                            # Get the unit of the depth
                            z_unit = ds["depth"].attrs.get('units', 'No unit attribute found')
                            depth_unit = lib.standard_units(z_unit)

                            units = "cgs"
                            if depth_unit == "km":
                                units = "mks"

                            logger.info(
                                    f"[INFO] The depth unit is detected as: {z_unit}, will use: {depth_unit}. Units set distance units based on: {units}"
                                )
                            start = slicer_lib.get_point("start")
                            if start == "back":
                                subset_type = "back"
                                break
                            end = slicer_lib.get_point("end")
                            if end == "back":
                                subset_type = "back"
                                break
                            depth = slicer_lib.get_point("depth")
                            if depth == "back":
                                subset_type = "back"
                                break

                            try:
                                depth = [
                                    min(float(depth[0]), float(depth[1])),
                                    max(float(depth[0]), float(depth[1])),
                                ]
                            except Exception as ex:
                                logger.error(f"[ERR] Bad depth values {depth}\n{ex}")
                                subset_type = "back"
                                break

                            plot_data = ds.copy()
                            plot_data = plot_data.where(
                                (plot_data.depth >= float(depth[0]))
                                & (plot_data.depth <= float(depth[1])),
                                drop=True,
                            )
                            utm_zone = None
                            meta = ds.attrs
                            if "grid_mapping_name" not in meta:
                                messages.append(f"[WARN] The 'grid_mapping_name' attribute not found. Assuming geographic coordinate system")
                                grid_mapping_name = "latitude_longitude"
                            else:
                                grid_mapping_name = meta["grid_mapping_name"]
                                messages.append(f"[INFO] grid_mapping_name is {grid_mapping_name}")

                            # Cross-section interpolation type.
                            interp_type = interpolation_method[0]
                            if verbose:
                                logger.info(
                                    f"\nSelect the interpolation method for creating the cross-section. The {', '.join(interpolation_method)} methods are available.\n\tPress return to select the default {interp_type} method\n\t{dash}\n\tback - takes you to the previous step\n\texit "
                                )
                            messages = output_messages(messages)
                            _interp_type = input(
                                f"Interpolation Method [{', '.join(interpolation_method+['back', 'exit'])}, default: {interp_type}]? "
                            )
                            if _interp_type.strip() == "exit":
                                sys.exit()
                            elif _interp_type.strip() == "back":
                                break
                            elif _interp_type.strip():
                                interp_type = _interp_type
                            elif interp_type not in (interpolation_method):
                                messages.append(f"[WARN] Invalid interpolation method of {interp_type}. Will use {interpolation_method[0]}")
                                interp_type = interpolation_method[0]

                            # Steps in the cross-section.
                            steps = xsection_steps
                            if verbose:
                                logger.info(
                                    f"\nThe number of points along the geodesic between the start and the end point (including the end points) to use in the cross section.\n\tPress enter to select the default of 100.\n\t{dash}\n\tback - takes you to the previous step\n\texit "
                                )
                            messages = output_messages(messages)
                            steps = input(
                                f"Number of points ['back', 'exit', default: {steps}]? "
                            )
                            if not steps.strip():
                                steps = xsection_steps
                            elif steps.strip() == "exit":
                                sys.exit()
                            elif steps.strip() == "back":
                                break
                            elif not steps.isnumeric():
                                messages.append(f"[WARN] Invalid number of steps, expected an integer. Will use {xsection_steps}")
                                steps = xsection_steps
                            else:
                                steps = int(steps)

                            
                            if verbose:
                                logger.info(
                                    f"\nVertical exaggeration for the cross section.\n\tPress enter to select the default of {vertical_exaggeration}.\n\t{dash}\n\tback - takes you to the previous step\n\texit "
                                )
                            messages = output_messages(messages)
                            exaggeration = input(
                                f"The section vertical exaggeration (0: dynamic aspect ratio)['back', 'exit', default: {vertical_exaggeration}]? "
                            )
                            relabel_y = False
                            if not exaggeration.strip():
                                exaggeration = vertical_exaggeration
                            elif exaggeration.strip() == "exit":
                                sys.exit()
                            elif exaggeration.strip() == "back":
                                break
                            elif not exaggeration.replace(".", "").isnumeric():
                                messages.append(f"[WARN] Invalid vertical exaggeration, expected a number. Will use {vertical_exaggeration}")
                                exaggeration = vertical_exaggeration
                            else:
                                exaggeration = float(exaggeration)
                                relabel_y = True
                            # Extract the cross-section.
                            try:
                                xsection_data, latitudes, longitudes = slicer_lib.interpolate_path(
                                        plot_data,
                                        start,
                                        end,
                                        num_points=steps,
                                        method=interp_type,
                                        grid_mapping=grid_mapping_name,
                                        utm_zone=meta["utm_zone"],
                                        ellipsoid=meta["ellipsoid"],
                                    )
                            except Exception as ex:
                                message = f"[ERR] interpolate_path failed: {ex}\n{traceback.print_exc()}"
                                logger.error(message)
                                break

                            # Calculate distances between consecutive points
                            geod = Geod(ellps=meta["ellipsoid"])
                            _, _, distances = geod.inv(
                                longitudes[:-1], latitudes[:-1], longitudes[1:], latitudes[1:]
                            )

                            # Compute cumulative distance, starting from 0
                            cumulative_distances = np.concatenate(([0], np.cumsum(distances)))
                            
                            if "grid_mapping_name" not in meta:
                                messages.append(f"[WARN] The 'grid_mapping_name' attribute not found. Assuming geographic coordinate system")
                                grid_mapping_name = "latitude_longitude"
                                
                            else:
                                grid_mapping_name = meta["grid_mapping_name"]
                            
                            
                            if units == "mks":
                                cumulative_distances = cumulative_distances / 1000.0
                                distance_label = "distance (km)"
                            else:
                                distance_label = "distance (m)"

                            # Create a new coordinate 'distance' based on the cumulative distances
                            xsection_data = xsection_data.assign_coords(distance=("points", cumulative_distances))

                            # If you want to use 'distance' as a dimension instead of 'index',
                            # you can swap the dimensions (assuming 'index' is your current dimension)
                            xsection_data = xsection_data.swap_dims({"points": "distance"})

                            # Iterate through the model variables and plot each cross-section.
                            plot_var = data_var[0]
                            while plot_var:

                                if len(data_var) > 1:
                                    if verbose:
                                        logger.info(
                                            f"\nThe model variable to plot.\n\tback - takes you to the previous step\n\texit "
                                        )
                                    messages = output_messages(messages)
                                    plot_var = input(
                                        f"Variable to plot {data_var}, back, exit]: "
                                    )
                                # Done, back!
                                if plot_var.strip() == "exit":
                                    sys.exit()
                                elif plot_var.strip() == "back" or not option.strip():
                                    break
                                elif plot_var not in data_var:
                                    messages.append(f"[ERR] Invalid input {plot_var}")
                                    continue
                                pdata = xsection_data.copy()
                                pdata["depth"] = -pdata["depth"]

                                if vmin and vmax:
                                    pdata[plot_var].plot.contourf(cmap=cmap, vmin=vmin, vmax=vmax)
                                else:
                                       pdata[plot_var].plot.contourf(
                                    cmap=cmap,
                                )
                                # Get the current axes
                                ax = plt.gca()
                                ax.set_xlabel(distance_label)
                                if relabel_y:
                                    # Retrieve the current y-axis label
                                    current_label = ax.get_ylabel()

                                    # Append new text to the existing label
                                    new_label = f"{current_label}; VE x{exaggeration}"

                                    # Set the new label
                                    ax.set_ylabel(new_label)
                                    plt.gca().set_aspect(exaggeration)
                                
                                plt.xticks(rotation=45)  # Rotating the x-axis labels to 45 degrees

                                # Adding vertical text for start and end locations
                                plt.text(
                                    cumulative_distances[0],
                                    1.05,
                                    f"⟸{start}",
                                    rotation=90,
                                    transform=plt.gca().get_xaxis_transform(),
                                    verticalalignment="bottom",
                                    horizontalalignment='center',
                                    fontsize=9,
                                )
                                plt.text(
                                    cumulative_distances[-1],
                                    1.05,
                                    f"⟸{end}",
                                    rotation=90,
                                    transform=plt.gca().get_xaxis_transform(),
                                    verticalalignment="bottom",
                                    horizontalalignment='center',
                                    fontsize=9,
                                )

                                # Add title if the model name is available.
                                if "model" in meta:
                                    plt.title(meta["model"])

                                # Adjust the layout
                                plt.tight_layout()
                                plt.show()
                                break

                        elif subset_type == "slice":
                            # Slice the model.
                            slice_dir = "depth"
                            while slice_dir:
                                if verbose:
                                    logger.info(
                                        f"\nA slice cuts the model along one of the coordinate axes.\n\tdirection - direction of the slice, the coordinate to cut the model along\n\t{dash}\n\tback - takes you to the previous step\n\texit  "
                                    )
                                messages = output_messages(messages)
                                slice_dir = input(
                                    f"direction [{', '.join(main_coordinates+['back', 'exit'])}]? "
                                )
                                # Done, back!
                                if slice_dir.strip() == "exit":
                                    sys.exit()
                                elif (
                                    slice_dir.strip() == "back"
                                    or slice_dir not in main_coordinates
                                    or not slice_dir.strip()
                                ):
                                    messages.append(f"[ERR] invalid variable {slice_dir}")
                                    slice_dir = "back"
                                    break
                                # Explore what to do with the slice?
                                else:
                                    messages = output_messages(messages)
                                    slice_value = None
                                    while slice_value is None:
                                        slice_value = input(
                                            f"slice {slice_dir} [{np.nanmin(ds[slice_dir].data)} to {np.nanmax(ds[slice_dir].data)}, back, exit]? "
                                        )

                                        # Exit.
                                        if slice_value.strip() == "exit":
                                            sys.exit()
                                        # A step back.
                                        elif (
                                            slice_value.strip() == "back"
                                            or not slice_value.strip()
                                        ):
                                            slice_value = "back"
                                            break
                                        else:
                                            try:
                                                slice_value = float(slice_value)
                                            except Exception as ex:
                                                messages.append(
                                                    f"[ERR] invalid value {slice_value}\n{ex}"
                                                )
                                                slice_value = None
                                                slice_value = "back"
                                        if slice_value == "back":
                                            break
                                    # Actions.
                                    # while slice_value:
                                    if slice_value == "back":
                                        break
                                    closest_slice_value = lib.closest(
                                        coordinate_values[slice_dir],
                                        slice_value,
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
                                        messages.append(
                                            f"\nSlicing at the closest {slice_dir} of {closest_slice_value}"
                                        )

                                        if verbose:
                                            messages.append(
                                                f"\nSelect slice limits in the {dim} direction.\n\tlimits - provide minimum,maximum or press enter to accept the full range default \n\t{dash}\n\tback - takes you to the previous step\n\texit  "
                                            )
                                        messages = output_messages(messages)
                                        _limits = input(
                                            f"{dim} limits [default values { slice_limits[dim]}, back, exit]: "
                                        )
                                        # Done, back!
                                        if _limits.strip() == "exit":
                                            sys.exit()
                                        elif _limits.strip() == "back":
                                            break
                                        elif _limits.strip() and "," not in _limits:
                                            messages.append(
                                                f"[WARN] Invalid limits, using the full range of {slice_limits[dim]}."
                                            )
                                            _limits = ""

                                        if _limits:
                                            values = _limits.split(",")
                                            try:
                                                slice_limits[dim] = (
                                                    min(
                                                        float(values[0]),
                                                        float(values[1]),
                                                    ),
                                                    max(
                                                        float(values[0]),
                                                        float(values[1]),
                                                    ),
                                                )
                                            except Exception as ex:
                                                messages.append(
                                                    f"[ERR] Invalid limits {_limits}.\n{ex}"
                                                )
                                                _limits = "back"
                                                break
                                        else:
                                            _limits = slice_limits[dim]
                                    # Slice the data.
                                    if _limits == "back":
                                        break
                                    sliced_data = slicer_lib.slicer(
                                        ds,
                                        slice_dir,
                                        closest_slice_value,
                                        slice_limits,
                                    )
                                    messages.append(
                                        f"\nThe sliced data summary: \n{sliced_data}"
                                    )
                                    slice_action = "continue"
                                    # Actions.
                                    while slice_action != "back":

                                        if verbose:
                                            if gmap_option:
                                                logger.info(
                                                    f"\nWhat to do with the slice.\n\tplot2d - a 2D plot of {slice_dims}\n\tplot3d - a 3D plot of {slice_dims} and the model variable on the 3rd axis.\n\t\tThe plot is interactive and can be rotated.\n\tgmap - a 2D plot of {slice_dims} in geographical coordinate system\n\tChange the colormap for the plots. You can also provide cmap, vmin, and vmax, with vmin and vmax defining the minimum and maximum values for the colormap.\n\tsave - save the slice data\n\t{dash} \n\tback - takes you to the previous step\n\texit  "
                                                )
                                            else:
                                                logger.info(
                                                    f"\nWhat to do with the slice.\n\tplot2d - a 2D plot of {slice_dims}\n\tplot3d - a 3D plot of {slice_dims} and the model variable on the 3rd axis.\n\tThe plot is interactive and can be rotated.\n\tcmap - change the color map for the plots (NOTE: You may also provide cmap,vmin,vmax)\n\tsave - save the slice data\n\t{dash}\n\tback - takes you to the previous step\n\texit  "
                                                )
                                                messages = output_messages(messages)
                                        slice_action = input(
                                            f"Action [plot2d, plot3d{gmap_option}, cmap, save, back, exit]: "
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
                                                    messages=output_messages(messages)
                                                    plot_var = input(
                                                        f"Variable to plot {data_var}, back, exit]: "
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
                                                    messages.append(
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
                                                if vmin and vmax:
                                                    plot_data.plot(cmap=cmap, vmin=vmin, vmax=vmax)
                                                else:
                                                     plot_data.plot(cmap=cmap)
                                                # Adjust the layout
                                                plt.tight_layout()
                                                plt.show()
                                                if len(data_var) <= 1:
                                                    plot_var = "back"
                                        # 3D plot.
                                        elif slice_action == "plot3d":
                                            plot_var = data_var[0]
                                            while plot_var:
                                                if len(data_var) > 1:
                                                    messages=output_messages(messages)
                                                    plot_var = input(
                                                        f"Variable to plot {data_var}, back, exit]: "
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
                                                    messages.append(
                                                        f"[ERR] Invalid input {plot_var}"
                                                    )
                                                    continue

                                                if vmin and vmax:
                                                    sliced_data[plot_var].plot.surface(cmap=cmap, vmin=vmin, vmax=vmax)
                                                else:
                                                     sliced_data[plot_var].plot.surface(
                                                    cmap=cmap,
                                                )
                                                # Adjust the layout
                                                plt.tight_layout()
                                                plt.show()
                                                if len(data_var) <= 1:
                                                    plot_var = "back"
                                        # Geographic map.
                                        elif slice_action == "gmap" and gmap_option:
                                            plot_var = data_var[0]
                                            while plot_var:
                                                if len(data_var) > 1:
                                                    messages = output_messages(messages)
                                                    plot_var = input(
                                                        f"Variable to plot {data_var}, back, exit]: "
                                                    )
                                                # Done, back!
                                                if plot_var.strip() == "exit":
                                                    sys.exit()
                                                elif (
                                                    plot_var.strip() == "back"
                                                    or not option.strip()
                                                ):
                                                    plot_var = "back"
                                                    break
                                                elif plot_var not in data_var:
                                                    messages.append(
                                                        f"[ERR] Invalid input {plot_var}"
                                                    )
                                                    continue
                                                slicer_lib.gmap(
                                                    plot_var,
                                                    cmap,
                                                    gmap_limits,
                                                    sliced_data,
                                                    vmin=vmin, vmax=vmax,
                                                )
                                                if len(data_var) <= 1:
                                                    plot_var = "back"
                                        # Colormap.
                                        elif slice_action == "cmap":
                                            messages.append(f"\n[cmaps] {plt.colormaps()}")
                                            messages = output_messages(messages)
                                            cmap = input(
                                                f"Select a color map for the plot [default cmap {cmap} (NOTE: You may also provide cmap,vmin,vmax), back, exit]: "
                                            )
                                            items = cmap.split(",")
                                            if len(items) == 1:
                                                cmap = items[0]
                                                if cmap not in plt.colormaps():
                                                    messages.append(
                                                        f"[ERR] invalid cmap: {cmap}"
                                                    )
                                                    break
                                                    
                                            elif len(items) == 3:
                                                cmap = items[0]
                                                if cmap not in plt.colormaps():
                                                    messages.append(
                                                        f"[ERR] invalid cmap: {cmap}"
                                                    )
                                                    break
                                                vmin = float(items[1])
                                                vmax = float(items[2])
                                                messages.append(f"[INFO] CMAP:{cmap}, vmin:{vmin}, vmax:{vmax}")
                                            else:
                                                cmap = prop.cmap
                                                messages.append(
                                                    f"[ERR] invalid cmap: {",".join(items)}"
                                                )
                                                break
                                        # Save the slice data.
                                        elif slice_action == "save":
                                            if verbose:
                                                logger.info(
                                                    f"\nSave the data\n\tfilename -- filename for saving the file. Its extension determines the file format.\n\t\tfilename.csv to save as GeoCSV, filename.nc to output as netCDF\n\t{dash}\n\tback - takes you to the previous step\n\texit "
                                                )
                                            messages=output_messages(messages)
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
                                                filename = "back"
                                                break
                                            else:
                                                messages = output_data(sliced_data, filename, messages)
                                            

                                # What does the user want to do next?
                                messages = output_messages(messages)
                                slice_value = input(
                                    f"Slice at {slice_dir} of (between {np.nanmin(ds[slice_dir].data)} and {np.nanmax(ds[slice_dir].data)}) [back, exit]: "
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
