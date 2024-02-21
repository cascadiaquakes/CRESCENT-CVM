"""The following parameters define the model parameters for constructing the metadata. 
All parameters should follow the CF Metadata Conventions (https://cfconventions.org/)."""

# The delimiter used to separate data columns.
delimiter = {"data": " ", "geocsv": "|"}

# The netCDF output format (should be one of those defined in shared_prop.py).
netcdf_format = "NETCDF3C"

# Information on the model.
model = {
    "name": "MODEL NAME",
    "title": "MODEL TITLE",
    "summary": "MODEL SUMMARY",
    "reference": "SHORT REFERENCE and (YEAR)",
    "reference_pid": "TYPE:ID, for example doi:1234.567890",
    "revision": "DATA REVISION, start with r0.0",
    "Version": "DATA VERSION, start with v0.0",
}

# Repository information.
repository = {
    "repository_pid": "OBTAIN.FROM.EMC. Format is TYPE:ID, for example doi:1234.567890",
}

# Information on the corresponding author.
corresponding_author = {
    "author_name": "NAME OF THE Corresponding author",
    "author_email": "EMAIL OF THE CORRESPONDING AUTHOR",
    "author_institution": "CORRESPONDING AUTHOR'S INSTITUTION",
    "author_url": "CORRESPONDING AUTHOR'S WEBSITE",
}

# Geospatial limits specifies the geographical area covered by the dataset.
# NOTE: For 2D models, set geospatial_z_min and geospatial_z_max to None
geospatial = {
    "geospatial_lon_min": -180.0,
    "geospatial_lon_max": 180.0,
    "geospatial_on_units": "degrees_east",
    "geospatial_on_resolution": 1,
    "geospatial_lat_min": -90,
    "geospatial_lat_max": 90,
    "geospatial_lat_units": "degrees_north",
    "geospatial_lat_resolution": 1,
    "geospatial_vertical_min": 0,
    "geospatial_vertical_max": 600.0,
    "geospatial_vertical_units": "km",
    "geospatial_vertical_positive": "down",
}

# Miscellaneous information. Any additional parameters beyond those
# identified below, will be included as additional global attribute.
global_attrs = {
    "keywords": "seismic,Rayleigh waves,shear wave",
    "acknowledgment": "",
    "history": "[Date] History 1; [Date] History 2",
    "comments": "COMMENTS",
}

# A recognized grid mapping name, for more information see:
# https://cfconventions.org/cf-conventions/cf-conventions.html#appendix-grid-mappings
grid_mapping_name = "latitude_longitude"

# The Universal Transverse Mercator zone.
utm_zone = 10

# The Ellipsoid reference.
ellipsoid = "WGS84"

# The primary X coordinate (could be longitude for the
# geographic coordinate system, or easting for UTM projection, etc.:
x = {
    "column": "longitude",
    "variable": "longitude",
    "standard_name": "longitude",
    "long_name": "Longitude; positive east",
    "units": "degrees_east",
}

# The primary Y coordinate (could be latitude for the
# geographic coordinate system, or northing for UTM projection, etc.:
y = {
    "column": "latitude",
    "variable": "latitude",
    "long_name": "Latitude; positive north",
    "units": "degrees_north",
    "standard_name": "latitude",
}

# The vertical coordinate (could be depth, elevation, etc.:
# NOTE: Include only if the model is a 3D model.
z = {
    "column": "depth",
    "variable": "depth",
    "standard_name": "depth",
    "long_name": "depth below sea level",
    "units": "km",
}

# The auxiliary X coordinate.
# NOTE: Set the values only if a secondary coordinate system is desired. For example easting and
# northing of the UTM projection for the primary geographic coordinate system, of longitude
# and latitude of the geographic coordinates when the primary coordinates are in UTM.
#
# NOTE: set x2 = None to exclude addition of the auxiliary coordinates. The auxiliary X coordinate
# could be longitude for the geographic coordinate system, or easting for UTM projection, etc..
#
# NOTE: LEAVE the column field EMPTY if you want the code to compute the auxiliary coordinate.
x2 = {
    "column": "easting",
    "variable": "easting",
    "standard_name": "easting",
    "long_name": "easting; UTM",
    "units": "m",
}

# The auxiliary Y coordinate.
# NOTE: Include only if x2 above is included.
#
# NOTE: LEAVE the column field EMPTY if you want the code to compute the auxiliary coordinate.

y2 = {
    "column": "northing",
    "variable": "northing",
    "long_name": "northing; UTM",
    "units": "m",
    "standard_name": "northing",
}

# A list of model variable dictionaries. Each model variable is a dictionary.
variables = [
    {
        "column": "vs",
        "variable": "vs",
        "dimensions": 3,
        "long_name": "Shear-wave velocity",
        "display_name": "S Velocity (km/s)",
        "units": "km.s-1",
    },
]

# The following line must include all _column values defined above. The order depends on
# the column order of the data block:
data_columns = ["latitude", "longitude", "easting", "northing", "depth", "vs"]
