"""A parameter file for the simple plotter tool. 
The following values will be used as the default 
values by the simple_plotter tool."""

# The default netCDF input file.
filename = "Cascadia-ANT+RF-Delph2018.r0.1.nc"

# The default model depth for generating a horizontal slice.
depth = 20

# The default model variable.
variable = "vs"

# Valu range for the above variable.
vmin = 2.8
vmax = 4.7

# Figure information.
figure_size = (6, 10)
figure_size_s = (4, 4)
cmap = "jet_r"

# The coordinate system information.
x = "longitude"
y = "latitude"
auxiliary_x = "easting"
auxiliary_y = "northing"

# Minimum for the primary coordinate system.
y_min = 46
x_min = -122
