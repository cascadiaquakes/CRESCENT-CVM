# Cascadia Community Velocity Model Tools (CVM Tools)

Last updated: 2024-05-21 \
Release: r0.7

This repository contains Python 3 tools developed to support the Cascadia Community Velocity Models ([CVM](https://cascadiaquakes.org/cvm/)) project. These tools facilitate the storage, extraction, and visualization of the CVMs.

## Notes

Starting with release r0.6, the metadata template files have been changed from Python data files to text files. The **cvm_write** tool no longer supports the Python metadata files.

## File Formats

The supported file formats are GeoCSV and netCDF, with support for an HDF5-based format to be added next.

| Format | Type             | Description                                                                                                                                                                                                                                |
| ------ | ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| ASCII  | GeoCSV           | An extension of the human-readable CSV (Comma-Separated Values) format with extensive metadata ([more](http://geows.ds.iris.edu/files/geocsv/GeoCSV.pdf))                                                                                  |
| netCDF | netCDF 4 Classic | The netCDF 4 Classic provides the netCDF-4 performance benefits while using the classic data model. This format supports compression for smaller file sizes ([more](https://www.loc.gov/preservation/digital/formats/fdd/fdd000339.shtml)) |

### Metadata

All CVMs, regardless of format, include extensive metadata. The netCDF metadata follow the [CF Metadata Conventions](https://cfconventions.org/).

### Coordinate Systems

Both formats support projected coordinate variables (x and y) and the geographic latitude and longitude variables.

## Tools

The final package includes tools to facilitate creation, extraction, and visualization of the model data and metadata. Currently available tools are:

| Tool              | Application                  | Description                                                                                                                                                                                                                          |
| ----------------- | ---------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| cvm_writer        | Metadata and data storage    | Reads a text model metadata file and outputs the corresponding model metadata in GeoCSV and JSON formats. If a CSV model data file is also provided, then the output could also include the model in GeoCSV and/or netCDF formats.   |
| netcdf_to_geocsv  | Format conversion            | Convert a CVM netCDF file to GeoCSV and/or output its metadata in JSON                                                                                                                                                               |
| cvm_slicer        | Data extraction and plotting | To extract data from a CVM netCDF file. Users can interactively inspect the metadata, slice the data, plot, and save the sliced data along the existing coordinate planes and for vertical cross-sections in an arbitrary direction. |
| simple_plotter.py | Plotting                     | A simple Python code to extract and plot depth slices from the CVM files in netCDF format. This tool is mainly intended for checking the generated netCDF files.                                                                     |

## Installation

### Download

Clone the repository or download the files and unzip ([more information](https://docs.google.com/document/d/1NtTUNGl8H5UOGttmFx5cBKEd3Nlc3GbqYIZAxUkrqjs/edit?usp=sharing)).

### Install the required packages

- [Python](https://www.python.org/) 3
- Python modules listed in `requirements.txt`
  - Install these modules with `pip install -r requirements.txt`

This package has been tested under Python 3.12.0 on macOS 14.2.1. It may work with older Python 3 versions.

### Package testing ([more information](https://docs.google.com/presentation/d/1_pXp4xzALeo5iwJQjpp-jX8dkAWAWD3kr9EY0icYOdk/edit?usp=sharing))

The following velocity sample models are provided under the _sample-files_ directory. Please note that the **casc1.6_velmdl** model is very large, and its data and the corresponding netCDF file are available from the [Google Drive](https://drive.google.com/drive/folders/1JTN0GAf25IIFBqkTMmZCTL50VnTjiiFd?usp=sharing). The **Cascadia-ANT+RF-Delph2018** model is available from this repository.

| Model                     | Author(s)                                      | Description                                                                                                                                | Coordinates | Link                                                                      |
| ------------------------- | ---------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ | ----------- | ------------------------------------------------------------------------- |
| casc1.6_velmdl            | William Stephenson                             | P- and S-wave velocity models incorporating the Cascadia subduction zone for 3D earthquake ground motion simulations, version 1.6          | UTM         | [USGS](https://www.sciencebase.gov/catalog/item/59f1e68be4b0220bbd9dd4b4) |
| Cascadia-ANT+RF-Delph2018 | Jonathan Delph, Alan Levander, and Fenglin Niu | 3D vertical shear-wave velocity model of the Cascadian forearc from the joint inversion of ambient noise dispersion and receiver functions | Geographic  | [EMC](https://ds.iris.edu/ds/products/emc-cascadia_antrf_delph2018/)      |

**cvm_writer**: Use the cvm_writer to convert the Cascadia_ANT+RF_Delph2018 raw data to CVM model files

1. cd to the model directory _sample-files/Cascadia-ANT+RF-Delph2018_
2. The _Cascadia_ANT+RF_Delph2018.txt.gz_ file is the raw data file for the model.
3. Unzip it and copy it to a file (for example _Cascadia-ANT+RF_data.txt_) for testing. This will be our model's data file.
   All model data files need to have a header as their first line that identifies the columns in the file using the same delimiter as data between column names. Place the following
   header as the first line in _Cascadia-ANT+RF_data.txt_:
   <br />**longitude latitude depth vs**
4. Now, we need a model metadata file. The metadata files contain the model metadata as text. Templates for the metadata
   files are available under the _template_ folder. You can use either _metadata_template_detailed.txt_ or _metadata_template.txt_ for
   model's metadata. The only difference is that the _metadata_template_detailed.txt_ has a full description of the model parameters.
   For this test, you do not need to copy a template. We already have done this and the _Cascadia-ANT+RF_meta.txt_ under this directory contains the metadata
   for this model. Please look at the file's content and familiarize yourself with the model metadata variables and how they are defined.
5. Convert the metadata to GeoCSV by running the _../../src/cvm_writer.py_ code (use the _-h_ flag for usage information):
   <br />**../../src/cvm_writer.py -m Cascadia-ANT+RF_meta.txt -o Cascadia-ANT+RF-test -t metadata**<br />
   _-m Cascadia-ANT+RF_meta.txt_ tells cvm*writer that the metadata file is \_Cascadia-ANT+RF_meta.txt*<br />
   _-o Cascadia-ANT+RF-test_ tells it what the output filename should be. For metadata, a postfix of _\_metadata_ will be added to the filename.<br />
   _-t metadata_ tells the code that we only want the metadata files. After running this commend, you should see two new files in your model directory:
   <br />_ Cascadia-ANT+RF-test_metadata.csv -- the metadata file in GeoCSV format
   <br />_ Cascadia-ANT+RF-test_metadata.json -- the metadata file in JSON format
6. Now we create the model files by running:
   <br />**../../src/cvm_writer.py -m Cascadia-ANT+RF_meta.txt -o Cascadia-ANT+RF-test -t metadata -d Cascadia-ANT+RF_data.txt -t netcdf**<br />
   _-d Cascadia-ANT+RF_data.txt_ tells the code where it can find the model data (from step 3 above)
   <br />_-t netcdf_ tells it that we want the output to be in netCDF
   After running this command, we should see a new file _Cascadia-ANT+RF-test.nc_ which is our model in netCDF.<br />We could have used _-t geocsv_ to have the output as
   GeoCSV file or we could have combined them all and have _-t metadata,netcdf,geocsv_ to get all the outputs.
7. Look at the CSV and JSON files to familiarize yourself with these files.
8. To look at the netCDF file content, you will find a simple plotting script under the _src_ directory (_src/simple_plotter.py_)
   Again, use _-h_ flag to get the usage message. This simple plotter code reads the netCDF file, prints information on the data, extracts a depth slice, and
   displays it in geographic and UTM coordinates. It also extracts a small piece of the depth slice and plots it again.

   **NOTE**: The **simple_plotter** tool expects a python configuration file containing the model default plot values in the same directory where you are running the simple plotter tool. You can download a template from the template directory.

   To run the plotter script against the netCDF file you created, run:<br />

   **../../src/simple_plotter.py -i Cascadia-ANT+RF-test.nc**
   <br />As you run this command, the model information is outputted to the terminal, and a slice at the depth of 20 km is plotted. By closing this plot, another plot shows the same slice using the UTM coordinates. Finally, by closing the plot, you get a plot of the upper right section of the area.

9. Look at the model information printed to your terminal. It includes depth information. Select another depth (for example 78km).
   Now run the same command but tell it to plot the slice for the depth of 78 km:
   <br />**../../src/simple_plotter.py -i Cascadia-ANT+RF-test.nc -z 78**
   <br />which produces the same plots as before but for the depth of 78 km

**cvm_slicer**:
To interactively slice and plot sections of the netCDF file we generated in the last step, we can use the **cvm_slicer** tool (to use the cvm_slicer tool on the casc1.6_velmdl, download the corresponding netCDF file from the [Google Drive](https://drive.google.com/drive/folders/1JTN0GAf25IIFBqkTMmZCTL50VnTjiiFd) and follow the same steps outlined below).

1. cd to the model directory _sample-files/Cascadia-ANT+RF-Delph2018_
2. To run the slicer tool, execute the following command.<br />
   **NOTE:** The -v argument tells the slicer tool to run in verbose mode and provide more information on the required parameters. You can remove this tag to run the tool in a less verbose mode.<br />
   **../../src/cvm_slicer.py -v -i Cascadia-ANT+RF-test.nc**

**netcdf_to_geocsv**:
To convert a netCDF CVM file to the corresponding GeoCSV format and optionally output its JSON metadata.

In the next example, we use the **netcdf_to_geocsv** to convert the Cascadia-ANT+RF-test.nc file (from previous examples) to GeoCSV and output its metadata in JSON:

1. cd to the model directory _sample-files/Cascadia-ANT+RF-Delph2018_
2. **../../src/netcdf_to_geocsv.py -i Cascadia-ANT+RF-test.nc -g true -m true**
   <br />
   _-i Cascadia-ANT+RF-test.nc_ provides the name of the input file to convert<br />
   _-g true_ sets the GeoCSV conversion flag to True<br />
   _-m true_ sets the GeoCSV metadata output flag to True<br />
   _-o conversion_test_ sets the output filename <br />

**NOTE:** Output files will have the same name as the input file

## Package Contents

This package contains the following files:

- **README.md**: This file
- **requirements.txt**: List of the Python modules to install

- **documents/**

  - **2024-01-08-CRESCENT-CVM.pptx**: Slide presentation of 2024-01-08

- **lib/**

  - **netcdf_to_geocsv_lib.py**: Python library for netcdf_to_geocsv
  - **shared_lib.py**: Python library file shared between tools
  - **slicer_lib.py**: Python library file for cvm_slicer
  - **writer_lib.py**: Python library file for cvm_writer

- **prop/**

  - **netcdf_to_geocsv_prop.py**: Configuration file for netcdf_to_geocsv
  - **shared_prop.py**: Configuration file shared between tools
  - **slicer_prop.py**: Configuration file for cvm_slicer
  - **writer_prop.py**: Configuration file for cvm_writer

- **sample-files/**

  - **Cascadia-ANT+RF-Delph2018/**

    - **Cascadia-ANT+RF-Delph2018.r0.1.nc**: Model in netCDF format
    - **Cascadia-ANT+RF_meta.txt**: Model metadata file
    - **Cascadia_ANT+RF_Delph2018.txt.gz**: Model raw data
    - **simple_plotter_prop.py**: The simple_plotter configuration file

  - **casc1.6_velmdl/**
    - **README**: A note on where to get the model raw data
    - **casc1.6-velmdl.r0.1.nc**: Model in netCDF format
    - **casc16-velmdl_meta.txt**: Model metadata file
    - **simple_plotter_prop.py**: The simple_plotter configuration file

- **src/**

  - **cvm_writer.py**: The CVM writer tool
  - **cvm_slicer.py**: The CVM model slicer tool
  - **netcdf_to_geocsv.py**: The netCDF to GeoCSV conversion tool
  - **simple_plotter.py**: The simple_plotter tool

- **template/**
  - **metadata_template.txt**: Model metadata template
  - **metadata_template_detailed.txt**: Model metadata template with variables described
  - **simple_plotter_prop.py**: The simple_plotter configuration file

## Citation

### Comments, Questions, or Bug Reports

Please open an issue ticket or contact manochehr.bahavar@earthscope.org
