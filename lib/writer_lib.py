import os
import sys
import logging
import json
import numpy as np
import xarray as xr
import pandas
from tqdm import tqdm


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

    grid_mapping_name = get_param(params, "grid_mapping_name")
    metadata_g["global_attrs"]["grid_mapping_name"] = grid_mapping_name

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
            data_variables[var]["grid_mapping"] = "latitude_longitude"
        else:
            data_variables[var]["grid_mapping"] = "transverse_mercator"

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


def write_json_metadata(output, metadata_dict, var_dict, data_variables):
    """Write out the metadata JSON file.

    Keyword arguments:
    output -- [required] the output filename
    metadata_dict -- [required] metadata dictionary
    """
    json_file = f"{output}.json"
    metadata_dict["variables"] = var_dict
    metadata_dict["data_vars"] = data_variables
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
    df -- [required] the xarray data frame
    params -- [required] the parameter module
    output -- [required] the output filename
    metadata_dict -- [required] metadata dictionary
    var_dic -- [Required] dictionary of variables
    """
    delimiter = get_param(params, "delimiter")["geocsv"]
    metadata = get_geocsv_metadata(params, metadata_dict, var_dict)
    var_columns = list()
    for var in var_dict:
        if "column" not in var_dict[var]:
            logger.error(
                f"[ERR] Could not find 'column' definition for variable {var}."
            )
            return "Failed to write the GeoCSV file"

        var_columns.append(var_dict[var]["column"])

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
    _var = metadata[var]["column"]
    if _var not in df:
        logger.error(f"[ERR] DataFrame does not contain a [{_var}] column!")
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


def get_coords(df, metadata):
    """Build a dictionary of coordinates.

    Keyword arguments:
    df -- [required] the DataFrame
    metadata -- [required] metadata
    """
    coords = dict()

    # X and Y coordinates.
    coords["x"] = get_var_from_df(df, metadata, "x")
    logger.info(f"[INFO] collected X coordinate information: {coords['x']['var']}")

    coords["y"] = get_var_from_df(df, metadata, "y")
    logger.info(f"[INFO] collected Y coordinate information: {coords['y']['var']}")
    # Is this a 3D model?
    if metadata["z"] is not None:
        coords["z"] = get_var_from_df(df, metadata, "z")
        logger.info(f"[INFO] collected Z coordinate information: {coords['z']['var']}")

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
            logger.info(
                f"[INFO] collected auxiliary X coordinate information: {coords['x2']['var']}"
            )

            coords["y2"] = get_var_from_df(df, metadata, "y2")
            logger.info(
                f"[INFO] collected auxiliary Y coordinate information: {coords['y2']['var']}"
            )

        # Compute the aux coords.
        else:
            grid_mapping_name = metadata["global_attrs"]["grid_mapping_name"]
            if grid_mapping_name == "latitude_longitude":
                xy_to_latlon = False
                logger.info(
                    f"[INFO] grid_mapping_name: {grid_mapping_name}: compute auxiliary coordinates using longitude and latitude"
                )
            else:
                xy_to_latlon = True
                logger.info(
                    f"[INFO] grid_mapping_name: {grid_mapping_name}: compute longitude and latitude coordinates"
                )

            # Set the column for x2 the same as the variable name as it was either None or blank.
            metadata["x2"]["column"] = metadata["x2"]["variable"]

            x2 = np.empty((len(coords["y"]["data"]), len(coords["x"]["data"])))
            x2[:] = np.nan
            y2 = np.empty((len(coords["y"]["data"]), len(coords["x"]["data"])))
            y2[:] = np.nan
            for _yind, _y in enumerate(coords["y"]["data"]):
                for _xind, _x in enumerate(coords["x"]["data"]):
                    x2[_yind][_xind], y2[_yind, _xind] = lib.project_lonlat_utm(
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
            logger.info(
                f"[INFO] collected auxiliary X coordinate information: {coords['x2']['var']}"
            )

            coords["y2"] = {"var": variables["y"], "data": y2}
            logger.info(
                f"[INFO] collected auxiliary Y coordinate information: {coords['y2']['var']}"
            )

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
                var_grid[_var][_iz][_iy][_ix] = row[_var]
        else:
            for _var in list(data_variables.keys()):
                var_grid[_var][_iy][_ix] = row[_var]
    return var_grid


def build_xarray_dataframe(data_variables, coords, var_grid, metadata):
    """Create the Xarray DataFrame.

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


def build_xarray_dataset(xr_df, coords, metadata):
    """For a given Xarray DataFrame, build the corresponding Xarray Dataset.

    Keyword arguments:
    xr_df -- [required] the XArry DataFrame to use.
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
        var = _var
        # For x & y dimensions we use x and y, unless they are lat/lon.
        if _var not in ("x", "y", "x2", "y2", "z"):
            var = _var
        else:
            var = coords[_var]["var"]
        ds[var].attrs = netcdf_attrs(metadata, _var)
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
        df = pandas.concat(
            [
                chunk
                for chunk in tqdm(
                    pandas.read_csv(
                        data_file, chunksize=writer_prop.read_chunk_size, sep="\s+"
                    ),
                    desc="...loading data",
                )
            ]
        )
    else:
        # df = pandas.read_csv(data_file, delimiter=delimiter)
        df = pandas.concat(
            [
                chunk
                for chunk in tqdm(
                    pandas.read_csv(
                        data_file,
                        chunksize=writer_prop.read_chunk_size,
                        delimiter=delimiter,
                    ),
                    desc="...loading data",
                )
            ]
        )
    return df
