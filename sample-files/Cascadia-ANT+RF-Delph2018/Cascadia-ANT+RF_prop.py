"""The following parameters define the model parameters for constructing the metadata. 
All parameters should follow the CF Metadata Conventions (https://cfconventions.org/)."""

""" The delimiter used to separate data columns.

data: delimiter for the input CSV data file. 
geocsv: delimiter for the output GeoCSV file(S).
"""
delimiter = {"data": " ", "geocsv": "|"}

# The netCDF output format (should be one of those defined in shared_prop.py)
netcdf_format = "NETCDF3C"

"""Information on the model:

name : The model name (no spaces).
title: A short title for the model.
summary: A very short description for the model.
reference: A short reference including author names and date.
reference_pid: Persistent identifier for the model. Usually the reference's DOI.
data_revision: Data revision. This refers to the changes in the initial data. 
        Only changes to data creates a revision.
version: Data version. This goes beyond changes to data. For example, if the revised
        model uses an expanded data set.
"""
model = {
    "model": "Cascadia_ANT+RF_Delph2018",
    "title": "3D vertical shear-wave velocity model of the Cascadian forearc from the joint inversion of ambient noise dispersion and receiver functions",
    "summary": "Cascadia_ANT+RF_Delph2018 was created from the joint inversion of ambient noise Rayleigh waves dispersion measurements (8-50 seconds) and adaptive CCP-derived receiver functions (see Delph et al., 2015, 2017 for details of methodology; Delph et al., 2018 for details of this model).",
    "reference": "Delph, Levander, and Niu (2018)",
    "reference_pid": "doi:10.1029/2018gl079518",
    "data_revision": "r0.1",
    "version": "v0.0",
}


"""Repository information:

repository_pid:  Persistent identifier for the model in the repository. 
        Usually the repository provides a DOI.
"""

repository = {
    "repository_pid": "doi:10.17611/dp/cascadiaantrfd2018",
}


"""Information on the corresponding author:

author_name: Name of the author
author_email: Author's email.
author_institution: Author's institution.
author_url: "Author's website
"""
corresponding_author = {
    "author_name": "Jonathan R. Delph",
    "author_email": "jdelph@purdue.edu",
    "author_institution": "Purdue University",
    "author_url": "https://www.eaps.purdue.edu/delph/",
}


"""Geospatial limits specifies the geographical area covered by the dataset.

NOTE: For 2D models, set geospatial_vertical_min and geospatial_vertical_max to None

geospatial_lon_min: lower longitude limit
geospatial_lon_max: upper longitude limit
geospatial_lon_units: units
geospatial_lon_resolution: spacing resolution

geospatial_lat_min: lower latitude limit
geospatial_lat_max: upper latitude limit
geospatial_lat_units: units
geospatial_lat_resolution: spacing resolution

geospatial_vertical_min: minimum z
geospatial_vertical_max: maximum z
geospatial_vertical_units: units
geospatial_vertical_positive: positive direction
geospatial_vertical_units: Variable's units
"""
geospatial = {
    "geospatial_lon_min": -124.8,
    "geospatial_lon_max": -120,
    "geospatial_lon_units": "degrees_east",
    "geospatial_lon_resolution": 0.2,
    "geospatial_lat_min": 40,
    "geospatial_lat_max": 49,
    "geospatial_lat_units": "degrees_north",
    "geospatial_lat_resolution": 1,
    "geospatial_vertical_min": -3,
    "geospatial_vertical_max": 80,
    "geospatial_vertical_units": "km",
    "geospatial_vertical_positive": "down",
}


"""Miscellaneous information. Any additional parameters beyond those 
identified below, will be included as additional global attribute:

keywords: Keywords for the model.
acknowledgment: Acknowledgment
history: History of changes.
comments: Any comments on the model.
"""
global_attrs = {
    "keywords": "seismic,Rayleigh waves dispersion,shear wave,s wave,Cascadian forearc,ambient noise dispersion,receiver functions",
    "acknowledgment": "Model was provided by Jonathan R. Delph of Department of Earth, Environmental and Planetary Sciences, Rice University",
    "history": "[2024-01-03] Converted to netCDF via CRESCENT data tools.",
    "comments": "CRESCENT CVM tools development project.",
}


"""A recognized grid mapping name, for more information see:
https://cfconventions.org/cf-conventions/cf-conventions.html#appendix-grid-mappings

Currently supporting:
latitude_longitude -- This grid mapping defines the canonical 2D geographical coordinate 
        system based upon latitude and longitude coordinates. 
mercator -- Mercator

The grid_mapping_name parameter should be consistent with the x description below.
"""
grid_mapping_name = "latitude_longitude"


"""The Universal Transverse Mercator zone:

zone: The UTM zone to use when converting to/from latitude longitude
ellipsoid: Ellipsoid reference
"""
utm_zone = 10


"""The Ellipsoid reference:

ellipsoid: Ellipsoid reference
"""
ellipsoid = "WGS84"


"""
The primary X coordinate (could be longitude for the 
geographic coordinate system, or easting for UTM projection, etc.:

column: The column heading in the data file:
variable: the dimension variable name.
standard_name: a CF-compliant standard name. Please consult the CF Standard Name Table
(https://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html).
long_name: Text that describes the variable more fully.
units: Variable's units
"""
x = {
    "column": "longitude",
    "variable": "longitude",
    "standard_name": "longitude",
    "long_name": "Longitude; positive east",
    "units": "degrees_east",
}

"""
The primary Y coordinate (could be latitude for the 
geographic coordinate system, or northing for UTM projection, etc.:

column: The column heading in the data file:
variable: the dimension variable name.
standard_name: a CF-compliant standard name. Please consult the CF Standard Name Table
(https://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html).
long_name: Text that describes the variable more fully.
"""
y = {
    "column": "latitude",
    "variable": "latitude",
    "long_name": "Latitude; positive north",
    "units": "degrees_north",
    "standard_name": "latitude",
}


"""
The vertical coordinate (could be depth, elevation, etc.:
NOTE: Include only if the model is a 3D model.

column: The column heading in the data file:
variable: the dimension variable name.
standard_name: a CF-compliant standard name. Please consult the CF Standard Name Table
(https://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html).
long_name: Text that describes the variable more fully.
units: Variable's units
"""
z = {
    "column": "depth",
    "variable": "depth",
    "standard_name": "depth",
    "long_name": "depth below sea level",
    "units": "km",
}


"""
The auxiliary X coordinate.
NOTE: Set the values only if a secondary coordinate system is desired. For example easting and 
northing of the UTM projection for the primary geographic coordinate system, of longitude 
and latitude of the geographic coordinates when the primary coordinates are in UTM.

NOTE: set x2 = None to exclude addition of the auxiliary coordinates. The auxiliary X coordinate 
could be longitude for the geographic coordinate system, or easting for UTM projection, etc..

NOTE: LEAVE the column field EMPTY if you want the code to compute the auxiliary coordinate.

column: The column heading in the data file when the auxiliary coordinates are provided. Set column to None to
compute the auxiliary coordinates .
variable: the dimension variable name.
standard_name: a CF-compliant standard name. Please consult the CF Standard Name Table
(https://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html).
long_name: Text that describes the variable more fully.
units: Variable's units
"""
x2 = {
    "column": None,
    "variable": "easting",
    "standard_name": "easting",
    "long_name": "easting; UTM",
    "units": "m",
}

"""
NOTE: Include only if x2 above is included.

NOTE: LEAVE the column field EMPTY if you want the code to compute the auxiliary coordinate.

The auxiliary Y coordinate (could be longitude for the geographic coordinate system, 
or easting for UTM projection, etc.:

column: The column heading in the data file:
variable: the dimension variable name.
standard_name: a CF-compliant standard name. Please consult the CF Standard Name Table
(https://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html).
long_name: Text that describes the variable more fully.
units: Variable's units
"""
y2 = {
    "column": "northing",
    "variable": "northing",
    "long_name": "northing; UTM",
    "units": "m",
    "standard_name": "northing",
}

"""
A dictionary of model variable dictionaries. The variable names are the keys for the dictionary.

Each model variable is a dictionary that includes:
column: The column heading in the data file:
variable: the model variable name. It is recommended to use a CF-compliant standard name. 
Please consult the Guidelines for Construction of CF Standard Names
(https://cfconventions.org/Data/cf-standard-names/docs/guidelines.html).
long_name: Text that describes the variable more fully.
display_name: Text to use for variable display. 
"""
variables = {
    "vs": {
        "column": "vs",
        "variable": "vs",
        "dimensions": 3,
        "long_name": "Shear-wave velocity",
        "display_name": "S Velocity (km/s)",
        "units": "km.s-1",
    },
}


"""The following line must include all _column values defined above. The order depends on
the column order of the data block:
"""
data_columns = ["latitude", "longitude", "easting", "northing", "depth", "vs"]
