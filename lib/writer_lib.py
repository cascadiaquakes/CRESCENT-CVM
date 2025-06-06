import os
import sys
import logging
import json
import numpy as np
import xarray as xr
import pandas as pd
from tqdm import tqdm
from pyproj import Proj, transform


# Get the root _directory.
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

PROP_DIR = os.path.join(ROOT_DIR, "prop")
sys.path.append(PROP_DIR)
import shared_prop as prop
import writer_prop as writer_prop

LIB_DIR = os.path.join(ROOT_DIR, prop.LIB_DIR)
sys.path.append(LIB_DIR)
import shared_lib as lib
from shared_lib import get_param, get_module_param

logger = lib.get_logger()


def get_repository_info(params, writer_prop):
    """Collect the repository information.

    Keyword arguments:
    params -- [required] the parameter module
    writer_prop-- [required] the metadata properties module
    """

    repository_pid = get_param(params, "repository")
    repository = get_module_param(writer_prop, "repository")
    repository = lib.merge_dictionaries(repository, [repository_pid])

    return repository


def get_model_dictionary(params):
    """Collect the model information.

    Keyword arguments:
    params -- [required] the parameter module
    """
    model_dict = get_param(params, "model")
    if "model" not in model_dict:
        logger.error(
            f"[ERR] parameter 'model' is required but it is missing from the model section of the parameter file."
        )
        sys.exit(2)

    model_dict = lib.insert_after(model_dict, ("id", model_dict["model"]), "model")
    return model_dict


def get_h5_metadata(params, global_meta=False):
    """Assemble the metadata for HDF5 files.

    Keyword arguments:
    params -- [required] the parameter module
    global_meta -- [default False] is this a global parameter file
    """

    # Retrieve individual metadata dictionaries.
    metadata_g = dict()
    metadata_v = dict()
    if global_meta:
        metadata_g["repository"] = get_repository_info(params, writer_prop)
        metadata_g["conventions"] = get_module_param(writer_prop, "conventions")
        metadata_g["model"] = get_model_dictionary(params)
        metadata_g["geospatial"] = get_param(params, "geospatial")
        # Make sure geospatial min, max, and resolution are float.
        for _key in metadata_g["geospatial"]:
            if "_min" in _key or "_max" in _key or "_resolution" in _key:
                metadata_g["geospatial"][_key] = float(metadata_g["geospatial"][_key])

        metadata_g["corresponding_author"] = get_param(params, "corresponding_author")
        metadata_g["global_attrs"] = get_param(params, "global_attrs")

        grid_ref = get_param(params, "grid_ref")
        metadata_g["global_attrs"]["grid_ref"] = grid_ref

        utm_zone = get_param(params, "utm_zone")
        metadata_g["global_attrs"]["utm_zone"] = str(utm_zone)

        ellipsoid = get_param(params, "ellipsoid")
        metadata_g["global_attrs"]["ellipsoid"] = ellipsoid

        metadata_dict = lib.merge_dictionaries(
            dict(),
            [
                metadata_g["model"],
                metadata_g["conventions"],
                metadata_g["corresponding_author"],
                metadata_g["geospatial"],
                metadata_g["repository"],
                metadata_g["global_attrs"],
            ],
        )
        return metadata_dict, metadata_g

    metadata_v["x"] = get_param(params, "x")
    metadata_v["y"] = get_param(params, "y")
    metadata_v["z"] = get_param(params, "z", required=False)
    metadata_v["x2"] = get_param(params, "x2", required=False)
    metadata_v["y2"] = get_param(params, "y2", required=False)
    variables = get_param(params, "variables")

    data_variables = dict()
    for var in list(variables.keys()):
        data_variables[var] = dict()
        for key in variables[var]:
            data_variables[var][f"{key}"] = variables[var][key]

        # Only if there are auxiliary coordinates.
        if metadata_v["x2"] is not None:
            _vars = lib.cf_coordinate_names(
                metadata_v["x2"]["variable"], metadata_v["y2"]["variable"], aux=True
            )
            data_variables[var][f"coordinates"] = f'{_vars["x"]} {_vars["y"]}'

        if metadata_v["x"]["variable"] == "longitude":
            data_variables[var]["grid_ref"] = "latitude_longitude"
        else:
            data_variables[var]["grid_ref"] = "transverse_mercator"

    var_dict = dict()
    if metadata_v["x"] is not None:
        _vars = lib.cf_coordinate_names(
            metadata_v["x"]["variable"], metadata_v["y"]["variable"]
        )
        var_dict[_vars["x"]] = metadata_v["x"]
        var_dict[_vars["y"]] = metadata_v["y"]

    if metadata_v["z"] is not None:
        var_dict[metadata_v["z"]["variable"]] = metadata_v["z"]

    if metadata_v["x2"] is not None:
        _vars = lib.cf_coordinate_names(
            metadata_v["x2"]["variable"], metadata_v["y2"]["variable"], aux=True
        )
        var_dict[_vars["x"]] = metadata_v["x2"]
        var_dict[_vars["y"]] = metadata_v["y2"]

    var_dict = lib.merge_dictionaries(var_dict, [data_variables])

    return (
        metadata_v,
        var_dict,
        data_variables,
    )


def get_metadata(params):
    """Assemble the metadata.

    Keyword arguments:
    params -- [required] the parameter module
    """

    # Retrieve individual metadata dictionaries.
    metadata_g = dict()
    metadata_v = dict()
    metadata_g["repository"] = get_repository_info(params, writer_prop)
    metadata_g["conventions"] = get_module_param(writer_prop, "conventions")
    metadata_g["model"] = get_model_dictionary(params)
    metadata_g["geospatial"] = get_param(params, "geospatial")
    # Make sure geospatial min, max, and resolution are float.
    for _key in metadata_g["geospatial"]:
        if "_min" in _key or "_max" in _key or "_resolution" in _key:
            metadata_g["geospatial"][_key] = float(metadata_g["geospatial"][_key])

    metadata_g["corresponding_author"] = get_param(params, "corresponding_author")
    metadata_g["global_attrs"] = get_param(params, "global_attrs")

    grid_ref = get_param(params, "grid_ref")
    metadata_g["global_attrs"]["grid_ref"] = grid_ref

    utm_zone = get_param(params, "utm_zone")
    metadata_g["global_attrs"]["utm_zone"] = str(utm_zone)

    ellipsoid = get_param(params, "ellipsoid")
    metadata_g["global_attrs"]["ellipsoid"] = ellipsoid

    metadata_v["x"] = get_param(params, "x")
    metadata_v["y"] = get_param(params, "y")
    metadata_v["z"] = get_param(params, "z", required=False)
    metadata_v["x2"] = get_param(params, "x2", required=False)
    metadata_v["y2"] = get_param(params, "y2", required=False)
    variables = get_param(params, "variables")

    data_variables = dict()
    for var in list(variables.keys()):
        data_variables[var] = dict()
        for key in variables[var]:
            data_variables[var][f"{key}"] = variables[var][key]

        # Only if there are auxiliary coordinates.
        if metadata_v["x2"] is not None:
            _vars = lib.cf_coordinate_names(
                metadata_v["x2"]["variable"], metadata_v["y2"]["variable"], aux=True
            )
            data_variables[var][f"coordinates"] = f'{_vars["x"]} {_vars["y"]}'

        if metadata_v["x"]["variable"] == "longitude":
            data_variables[var]["grid_ref"] = "latitude_longitude"
        else:
            data_variables[var]["grid_ref"] = "transverse_mercator"

    var_dict = dict()
    if metadata_v["x"] is not None:
        _vars = lib.cf_coordinate_names(
            metadata_v["x"]["variable"], metadata_v["y"]["variable"]
        )
        var_dict[_vars["x"]] = metadata_v["x"]
        var_dict[_vars["y"]] = metadata_v["y"]

    if metadata_v["z"] is not None:
        var_dict[metadata_v["z"]["variable"]] = metadata_v["z"]

    if metadata_v["x2"] is not None:
        _vars = lib.cf_coordinate_names(
            metadata_v["x2"]["variable"], metadata_v["y2"]["variable"], aux=True
        )
        var_dict[_vars["x"]] = metadata_v["x2"]
        var_dict[_vars["y"]] = metadata_v["y2"]

    var_dict = lib.merge_dictionaries(var_dict, [data_variables])

    metadata_dict = lib.merge_dictionaries(
        dict(),
        [
            metadata_g["model"],
            metadata_g["conventions"],
            metadata_g["corresponding_author"],
            metadata_g["geospatial"],
            metadata_g["repository"],
            metadata_g["global_attrs"],
        ],
    )
    return (
        metadata_dict,
        lib.merge_dictionaries(metadata_v, [metadata_g]),
        var_dict,
        data_variables,
    )


def write_json_metadata(
    output, metadata_dict, var_dict, data_variables, z_var, size_kb
):
    """Write out the metadata JSON file.

    Keyword arguments:
    output -- [required] the output filename
    metadata_dict -- [required] metadata dictionary
    variables_dict -- [required] variables dictionary
    data_variables -- [required] data variables dictionary
    """
    json_file = f"{output}.json"
    metadata_dict["vars"] = list(var_dict.keys())
    metadata_dict["z_var"] = z_var
    if size_kb:
        metadata_dict["size_kb"] = round(size_kb, 2)
    metadata_dict["variables"] = var_dict
    metadata_dict["data_vars"] = list(data_variables.keys())

    with open(json_file, "w") as fp:
        json.dump(metadata_dict, fp, indent=writer_prop.json_indent)

    return f"Metadata JSON file: {json_file}"


def get_geocsv_metadata(params, metadata_dict, var_dict):
    """Get the GeoCSV metadata.

    Keyword arguments:
    params -- [required] the parameter module
    output -- [required] the output filename
    metadata_dict -- [required] metadata dictionary
    var_dic -- [Required] dictionary of variables
    """

    line_break = "\n"
    delimiter = get_param(params, "delimiter")["geocsv"]
    delimiter = delimiter.strip()
    if len(delimiter) <= 0:
        delimiter = " "
    geocsv_metadata = (
        f"# dataset: {get_module_param(writer_prop, 'dataset_version')}\n"
        f"# delimiter: {delimiter}"
    )

    for _key in metadata_dict:
        if _key.strip() not in ("variables", "data_vars"):
            geocsv_metadata = f"{geocsv_metadata}{line_break}# {writer_prop.global_prefix}{_key}: {metadata_dict[_key]}"
    for _var in var_dict:
        _dict = var_dict[_var]
        for _key in _dict:

            # The output geoCSV will have column the same as variable
            if _key == "column":
                geocsv_metadata = (
                    f"{geocsv_metadata}{line_break}# {_var}_{_key}: {_var}"
                )
            else:
                geocsv_metadata = (
                    f"{geocsv_metadata}{line_break}# {_var}_{_key}: {_dict[_key]}"
                )
    return geocsv_metadata


def write_geocsv_metadata(params, output, metadata_dict, var_dict):
    """Write out the metadata JSON file.

    Keyword arguments:
    params -- [required] the parameter module
    output -- [required] the output filename
    metadata_dict -- [required] metadata dictionary
    var_dic -- [Required] dictionary of variables
    """

    metadata = get_geocsv_metadata(params, metadata_dict, var_dict)
    geocsv_file = f"{output}.csv"
    with open(geocsv_file, "w") as fp:
        fp.write(metadata)

    return f"Metadata GeoCSV file: {geocsv_file}"


def write_geocsv_file(df, params, output, metadata_dict, var_dict):
    """Write out the the model in GeoCSV.

    Keyword arguments:
    df -- [required] the data frame
    params -- [required] the parameter module
    output -- [required] the output filename
    metadata_dict -- [required] metadata dictionary
    var_dic -- [Required] dictionary of variables
    """
    delimiter = get_param(params, "delimiter")["geocsv"]
    metadata = get_geocsv_metadata(params, metadata_dict, var_dict)
    var_columns = list()
    for var in var_dict:
        if "variable" not in var_dict[var]:
            logger.error(
                f"[ERR] Could not find 'variable' definition for variable {var}."
            )
            return "Failed to write the GeoCSV file"

        var_columns.append(var_dict[var]["variable"])

    geocsv_file = f"{output}.csv"
    with open(geocsv_file, "w") as fp:
        fp.write(f"{metadata}\n")
    df.to_csv(
        geocsv_file,
        mode="a",
        na_rep=prop.na_rep,
        sep=delimiter,
        columns=var_columns,
        index=False,
    )

    return f"Model GeoCSV file: {geocsv_file}"


def get_var_from_df(df, metadata, var):
    """Extract a coordinate variable from DataFrame.

    Keyword arguments:
    df -- [required] the DataFrame
    metadata -- [required] metadata
    var  -- [required] the variable to extract
    """
    _var = metadata[var]["variable"]
    if _var not in df:
        logger.error(f"[ERR] DataFrame does not contain a '{_var}' variable!")
        sys.exit(3)
    x = sorted(df[_var].unique())
    return {"var": metadata[var]["variable"], "data": x}


def add_aux_coord_columns(df, coords):
    """Add the auxiliary coordinate columns.

    Keyword arguments:
    df -- [required] the DataFrame
    coords -- [required] a dictionary of coordinate variables.
    """

    # Aux coordinates are already in the DataFrame, skip.
    if coords["x2"]["var"] in df:
        return df

    # Populate the variable grids.
    x2 = list()
    y2 = list()

    # Number of decimal places for the aux (set the same as main).
    fix = f'{coords["x"]["data"][0]}'[::-1].find(".")
    if fix < 0:
        fix = 0
    for idf in df.index:
        row = df.loc[idf]
        _xvar = coords["x"]["var"]
        _x = row[_xvar]
        _ix = coords["x"]["data"].index(_x)

        _yvar = coords["y"]["var"]
        _y = row[_yvar]
        _iy = coords["y"]["data"].index(_y)

        # As a default, variables with float types are attributed a _FillValue of NaN in the output file,
        # unless explicitly disabled with an encoding {'_FillValue': None}.
        x2.append(round(coords["x2"]["data"][_iy][_ix], fix))
        y2.append(round(coords["y2"]["data"][_iy][_ix], fix))
    df[coords["x2"]["var"]] = x2
    df[coords["y2"]["var"]] = y2

    return df


def get_coords(df, metadata, flat="variable", verbose=False):
    """Build a dictionary of coordinates.

    Keyword arguments:
    df -- [required] the DataFrame
    metadata -- [required] metadata
    flat -- [default variable] return the auxiliary coordinates as multi dimensional arrays.
            if set to "2d", it will return it as 2D array of latitude, longitude
            if set to "flat" as a flat array
    verbose -- [default False] Run in verbose mode
    """
    coords = dict()
    if verbose:
        logger.info(
            f"\n\n[INFO] Extracting coordinate information from DataFrame: {df}"
        )
    # X and Y coordinates.
    coords["x"] = get_var_from_df(df, metadata, "x")
    if verbose:
        logger.info(
            f"[INFO] collected X coordinate information: {coords['x']['var']} with length: {len(coords['x']['data'])}"
        )

    coords["y"] = get_var_from_df(df, metadata, "y")
    if verbose:
        logger.info(
            f"[INFO] collected Y coordinate information: {coords['y']['var']} with length: {len(coords['y']['data'])}"
        )
    # Is this a 3D model?
    if metadata["z"] is not None:
        coords["z"] = get_var_from_df(df, metadata, "z")
        if verbose:
            logger.info(
                f"[INFO] collected Z coordinate information: {coords['z']['var']} with length: {len(coords['z']['data'])}"
            )

    # Auxiliary coordinates?
    if metadata["x2"] is not None:
        # Aux data info provided in the data file.
        aux_data = True

        if (
            metadata["x2"]["column"].lower() == "none"
            or len(metadata["x2"]["column"].strip()) <= 0
        ):
            aux_data = False

        if aux_data:
            coords["x2"] = get_var_from_df(df, metadata, "x2")
            if verbose:
                logger.info(
                    f"[INFO] collected auxiliary X coordinate information: {coords['x2']['var']}"
                )

            coords["y2"] = get_var_from_df(df, metadata, "y2")
            if verbose:
                logger.info(
                    f"[INFO] collected auxiliary Y coordinate information: {coords['y2']['var']}"
                )

        # Compute the aux coords.
        else:
            grid_ref = metadata["global_attrs"]["grid_ref"]
            if grid_ref == "latitude_longitude":
                xy_to_latlon = False
                if verbose:
                    logger.info(
                        f"[INFO] grid_ref: {grid_ref}: compute auxiliary coordinates using longitude and latitude"
                    )
            else:
                xy_to_latlon = True
                if verbose:
                    logger.info(
                        f"[INFO] grid_ref: {grid_ref}: compute longitude and latitude coordinates"
                    )

            # Set the column for x2 the same as the variable name as it was either None or blank.
            metadata["x2"]["column"] = metadata["x2"]["variable"]
            if flat.lower() == "2d":
                unique_y = np.unique(coords["y"]["data"])
                unique_x = np.unique(coords["x"]["data"])
                x2 = np.empty((len(unique_y), len(unique_x)))
                x2[:] = np.nan
                y2 = np.empty((len(unique_y), len(unique_x)))
                y2[:] = np.nan
                for _yind, _y in enumerate(unique_y):
                    for _xind, _x in enumerate(unique_x):
                        x2[_yind][_xind], y2[_yind][_xind] = lib.project_lonlat_utm(
                            _x,
                            _y,
                            metadata["global_attrs"]["utm_zone"],
                            metadata["global_attrs"]["ellipsoid"],
                            xy_to_latlon=xy_to_latlon,
                        )
            elif flat.lower() == "flat":

                # Define the UTM to Lat/Lon transformer
                utm_proj = Proj(
                    proj="utm",
                    zone=metadata["global_attrs"]["utm_zone"],
                    ellps=metadata["global_attrs"]["ellipsoid"],
                )
                latlon_proj = Proj(
                    proj="latlong", datum=metadata["global_attrs"]["ellipsoid"]
                )

                # Convert UTM to Latitude and Longitude
                if xy_to_latlon:
                    y2, x2 = transform(
                        utm_proj,
                        latlon_proj,
                        df[metadata["x"]["variable"]].values,
                        df[metadata["y"]["variable"]].values,
                    )
                else:
                    x2, y2 = utm_proj(df["longitude"].values, df["latitude"].values)

            else:
                x2 = np.empty((len(coords["y"]["data"]), len(coords["x"]["data"])))
                x2[:] = np.nan
                y2 = np.empty((len(coords["y"]["data"]), len(coords["x"]["data"])))
                y2[:] = np.nan
                for _yind, _y in enumerate(coords["y"]["data"]):
                    for _xind, _x in enumerate(coords["x"]["data"]):
                        x2[_yind][_xind], y2[_yind][_xind] = lib.project_lonlat_utm(
                            _x,
                            _y,
                            metadata["global_attrs"]["utm_zone"],
                            metadata["global_attrs"]["ellipsoid"],
                            xy_to_latlon=xy_to_latlon,
                        )

            variables = lib.cf_coordinate_names(
                metadata["x2"]["variable"], metadata["y2"]["variable"], aux=True
            )
            coords["x2"] = {"var": variables["x"], "data": x2}
            if verbose:
                logger.info(
                    f"[INFO] collected auxiliary X coordinate information: {coords['x2']['var']}"
                )

            coords["y2"] = {"var": variables["y"], "data": y2}
            if verbose:
                logger.info(
                    f"[INFO] collected auxiliary Y coordinate information: {coords['y2']['var']}"
                )
    if verbose:
        logger.info("\n\n")
    return coords


def init_grid(coords):
    """Initialize a 2D or 3D grid.

    Keyword arguments:
    coords -- [required] a dictionary of coordinate variables.
    """
    if "z" not in coords:
        grid = np.empty(
            (
                len(coords["y"]["data"]),
                len(coords["x"]["data"]),
            )
        )
        grid[:] = np.nan
    else:
        grid = np.empty(
            (
                len(coords["z"]["data"]),
                len(coords["y"]["data"]),
                len(coords["x"]["data"]),
            )
        )
        grid[::] = np.nan

    return grid


def netcdf_attrs(data_variables, var):
    """Extract variable attributes and make sure they conform to the netCDF style.

    Keyword arguments:
    data_variables -- [required] a dictionary of data variables.
    car -- [required] the variable to extract attributes for.
    """
    attrs = dict()
    for key in data_variables[var]:
        _key = key
        if _key.startswith(f"{var}_"):
            _key = _key[len(f"{var}_") :]
        if _key in ("column", "dimensions"):
            continue
        value = data_variables[var][key]
        if _key == "keywords":
            value = ",".join(data_variables[var][_key])
        attrs[_key] = value
    return attrs


def create_var_grid(df, data_variables, coords, grid):
    """Create a grid for each variable.

    Keyword arguments:
    df -- [required] the DataFrame
    data_variables -- [required] a dictionary of data variables.
    coords -- [required] a dictionary of coordinate variables.
    grid -- [required] a 2D or 3D grid to use as a template.
    """
    var_grid = dict()
    # Initialize variable grids.
    for _var in data_variables.keys():
        var_grid[_var] = np.copy(grid)

    # Populate the variable grids.
    for ix in df.index:
        row = df.loc[ix]
        _xvar = coords["x"]["var"]
        _x = row[_xvar]
        _ix = coords["x"]["data"].index(_x)

        _yvar = coords["y"]["var"]
        _y = row[_yvar]
        _iy = coords["y"]["data"].index(_y)

        if "z" in coords:
            _zvar = coords["z"]["var"]
            _z = row[_zvar]
            _iz = coords["z"]["data"].index(_z)

            for _var in list(data_variables.keys()):
                var_grid[_var][_iz][_iy][_ix] = row[data_variables[_var]["variable"]]
        else:
            for _var in list(data_variables.keys()):
                var_grid[_var][_iy][_ix] = row[data_variables[_var]["variable"]]
    return var_grid


def group_data(df, coords):
    """Group the data by the primary coordinates.

    Keyword arguments:
    df -- [required] the DataFrame
    coords -- [required] a dictionary of coordinate variables.
    """

    grouped = df.groupby([coords["x"]["var"], coords["y"]["var"]])
    return grouped


def build_dataframe(data_variables, coords, var_grid, metadata):
    """Create a DataFrame.

    Keyword arguments:
    data_variables -- [required] a dictionary of data variables.
    coords -- [required] a dictionary of coordinate variables.
    var_grid -- [required] grid of variables.
    metadata -- [required] metadata
    """

    xr_df = dict()
    for _var in data_variables.keys():
        names1 = lib.cf_coordinate_names(coords["x"]["var"], coords["y"]["var"])
        if "x2" in coords:
            names2 = lib.cf_coordinate_names(
                coords["x2"]["var"], coords["y2"]["var"], aux=True
            )

        if "x2" in coords and "z" in coords:
            xr_df[_var] = xr.DataArray(
                var_grid[_var],
                coords={
                    coords["z"]["var"]: coords["z"]["data"],
                    names1["y"]: coords["y"]["data"],
                    names1["x"]: coords["x"]["data"],
                    names2["x"]: (
                        (names1["y"], names1["x"]),
                        coords["x2"]["data"],
                    ),
                    names2["y"]: (
                        (names1["y"], names1["x"]),
                        coords["y2"]["data"],
                    ),
                },
                dims=[coords["z"]["var"], names1["y"], names1["x"]],
                attrs=netcdf_attrs(data_variables, _var),
            )
        elif "x2" in coords:
            xr_df[_var] = xr.DataArray(
                var_grid[_var],
                coords={
                    names1["y"]: coords["y"]["data"],
                    names1["x"]: coords["x"]["data"],
                    names2["x"]: (
                        (names1["y"], names1["x"]),
                        coords["x2"]["data"],
                    ),
                    names2["y"]: (
                        (names1["y"], names1["x"]),
                        coords["y2"]["data"],
                    ),
                },
                dims=[names1["y"], names1["x"]],
                attrs=netcdf_attrs(data_variables, _var),
            )
        else:
            if "z" in coords:
                xr_df[_var] = xr.DataArray(
                    var_grid[_var],
                    coords=(
                        coords["z"]["data"],
                        coords["y"]["data"],
                        coords["x"]["data"],
                    ),
                    dims=(
                        coords["z"]["var"],
                        names1["y"],
                        names1["x"],
                    ),
                    attrs=netcdf_attrs(data_variables, _var),
                )
            else:
                xr_df[_var] = xr.DataArray(
                    var_grid[_var],
                    coords=(
                        coords["y"]["data"],
                        coords["x"]["data"],
                    ),
                    dims=(
                        names1["y"],
                        names1["x"],
                    ),
                    attrs=netcdf_attrs(data_variables, _var),
                )
    return xr_df


def build_dataset(xr_df, coords, metadata):
    """For a given a DataFrame, build the corresponding Dataset.

    Keyword arguments:
    xr_df -- [required] the DataFrame to use.
    metadata -- [required] metadata
    """
    ds = xr.Dataset(
        data_vars=xr_df,
        attrs=lib.merge_dictionaries(
            metadata["model"],
            [
                metadata["conventions"],
                metadata["corresponding_author"],
                metadata["repository"],
                metadata["geospatial"],
                metadata["global_attrs"],
            ],
        ),
    )

    for _var in coords:

        """
        var = _var
        # For x & y dimensions we use x and y, unless they are lat/lon.
        if _var not in ("x", "y", "x2", "y2", "z"):
            var = _var
        else:
            var = coords[_var]["var"]
        """
        if _var in ds:
            ds[_var].attrs = netcdf_attrs(metadata, _var)
        else:
            ds[coords[_var]["var"]].attrs = netcdf_attrs(metadata, _var)
    return ds


def write_netcdf_file(output, xr_dset, nc_format):
    """Output the dataset as netCDF file..

    Keyword arguments:
    output -- [required] the output filename
    mparams -- [required] the Xarray Dataset to output
    nc_format -- [required] the output netCDF type
    """
    netcdf_file = f"{output}.nc"
    if nc_format["encoding"] is None:
        xr_dset.to_netcdf(netcdf_file)  # , format=nc_format["format"])
    else:
        xr_dset.to_netcdf(
            netcdf_file
        )  # , format=nc_format["format"], engine="netcdf4", encoding=enc)

    return f"Model netCDF file: {netcdf_file}"


def read_csv(
    data_file,
    params,
):
    """Load the CSV data file and show the progress.

    Keyword arguments:
    data_file -- [required] input CSV data file
    mparams -- [required] the parameter module
    """
    delimiter = params["delimiter"]["data"].strip()
    if len(delimiter) <= 0:
        df = pd.concat(
            [
                chunk
                for chunk in tqdm(
                    pd.read_csv(
                        data_file, chunksize=writer_prop.read_chunk_size, sep=r"\s+"
                    ),
                    desc="...loading data",
                )
            ]
        )
    else:
        # df = pd.read_csv(data_file, delimiter=delimiter)
        df = pd.concat(
            [
                chunk
                for chunk in tqdm(
                    pd.read_csv(
                        data_file,
                        chunksize=writer_prop.read_chunk_size,
                        delimiter=delimiter,
                    ),
                    desc="...loading data",
                )
            ]
        )

    # Rename the DataFrame columns and set them to variable name.
    renamed_columns = dict()
    # For the coordinate variables.
    for column in df.columns:
        for par in params:
            if "column" in params[par]:
                if params[par]["column"] == column:
                    renamed_columns[column] = params[par]["variable"]
                    params[par]["column"] = params[par]["variable"]

    # For model variables
    for column in df.columns:
        for var in params["variables"]:
            if "column" in params["variables"][var]:
                if column == params["variables"][var]["column"]:
                    renamed_columns[column] = params["variables"][var]["variable"]
                    params["variables"][var]["column"] = params["variables"][var][
                        "variable"
                    ]

    if renamed_columns:
        df.rename(columns=renamed_columns, inplace=True)

    return df, params
