
 # Cascadia Community Velocity Model Tools (CVM Tools)

21Last update: 2024-02-14

Python 3 repository of tools developed to support the Cascadia Community Velocity Models ([CVM](https://cascadiaquakes.org/cvm/)) project by facilitating storage, extraction, and visualization of the CVMs. 

##  File Formats

The supported file formats are GeoCSV and netCDF with support for an HDF5-based format to be added next.
|Format| Type| Description|
|------|-----|------------|
|ASCII|GeoCSV| An extension of the human readable CSV (Comma-Separated Values) format with extensive metadata ([more](http://geows.ds.iris.edu/files/geocsv/GeoCSV.pdf))|
|netCDF|netCDF 4 Classic|The netCDF 4 Classic provides the netCDF-4 performance benefits while using the classic data model. This format supports compression for smaller file sizes ([more](https://www.loc.gov/preservation/digital/formats/fdd/fdd000339.shtml))|

### Metadata
All CVMs, regardless of format, include extensive metadata that follow the [CF Metadata Conventions](https://cfconventions.org/).

### Coordinate Systems
Both formats support projected coordinate variables (x and y) and the geographic latitude and longitude variables.

## Tools
The final package will include tools to facilitate creation, extraction and visualization of the model data and metadata. Currently the available tools are:

|Tool| Application| Description|
|------|-----|------------|
|cvm_writer|Metadata and data storage| Reads an input Python-based configuration file that holds model metadata and outputs the metadata in GeoCSV and JSON formats. If a CSV model data file is also provided, then the output would also include the model in GeoCSV and/or netCDF formats.|
|netcdf_to_geocsv|Format conversion| Convert a CVM netCDF file to GeoCSV|
|cvm_slicer|Data extraction and plotting|Extracts data from a CVM netCDF file. Interactively, allows user to inspect the metadata, slice the data, plot, and save the sliced data. Currently slicing is only performed along the existing coordinate planes.|
|simple_plotter.py|Plotting| A simple Python code to extract and plot depth slices from the CVM files in netCDF format. This tool is mainly intended for checking the generated netCDF files.|

## Installation

### Download

Clone the repository or download the files and unzip.

### Install the required packages

* [Python](https://www.python.org/) 3
* Python modules listed in `requirements.txt`
  * Install these modules with `pip install -r requirements.txt`

This package has been tested under Python 3.12.0 on macOS 14.2.1, it may work with older Python 3 versions.

### Package testing

The following velocity sample models are provided under the *sample-files* directory. Please note that The **casc1.6_velmdl** model is very large and should not be used as a test model. Use the **Cascadia-ANT+RF-Delph2018** model for your initial tests.

| Model | Author(s) | Description | Coordinates| link  |    
|------|----------------|----------------------------|-------------|-------------|
| casc1.6_velmdl | William Stephenson | P- and S-wave velocity models incorporating the Cascadia subduction zone for 3D earthquake ground motion simulations, version 1.6 | UTM| [USGS](https://www.sciencebase.gov/catalog/item/59f1e68be4b0220bbd9dd4b4) |
| Cascadia-ANT+RF-Delph2018 | Jonathan Delph, Alan Levander, and Fenglin Niu | 3D vertical shear-wave velocity model of the Cascadian forearc from the joint inversion of ambient noise dispersion and receiver functions | Geographic| [EMC](https://ds.iris.edu/ds/products/emc-cascadia_antrf_delph2018/) |

**cvm_writer**: Use the cvm_writer to convert the Cascadia_ANT+RF_Delph2018 raw data to CVM model files 
1. cd to the model directory *sample-files/Cascadia-ANT+RF-Delph2018*
2. The *Cascadia_ANT+RF_Delph2018.txt.gz* file is the raw data file for the model. 
3. Unzip it and copy it to a file (for example *Cascadia-ANT+RF_data.txt*) for testing. This will be our model's data file.
   All model data files need to have a header as their first line that identifies the columns in the file using the same delimiter as data between column names. Place the following
   header as the first line in *Cascadia-ANT+RF_data.txt*:
   <br />**longitude latitude depth vs**
4. Now, we need a model metadata file. The metadata files store model metadata as Python variables. Templates for the metadata 
   files are provided under the *template* folder. You could use either *metadata_template_detailed.py* or *metadata_template.py* for
   model's metadata. The only difference is that the  *metadata_template_detailed.py* as a full description of the model parameters.
   For this test, you do not need to copy a template. We already have done this and the *Cascadia-ANT+RF_prop.py* under this directory contains the metadata
   for this model. Please look at the file's content and familiarize yourself with the model metadata variables.
5. Let us convert the metadata to GeoCSV by running the *../../src/cvm_writer.py* code (use the *-h* flag for usage information):
   <br />**../../src/cvm_writer.py -m Cascadia-ANT+RF_prop.py -o Cascadia-ANT+RF-test -t metadata**<br />
   *-m Cascadia-ANT+RF_prop.py* tells cvm_writer that the metadata file is *Cascadia-ANT+RF_prop.py*<br />
   *-o Cascadia-ANT+RF-test* tells it what the output filename should be. For metadata, a postfix of *_metada* will be added to the filename.<br />
   *-t metadata* tells the code that we only want the metadata files. After running this commend, you should see two new files in your model directory: 
        <br />* Cascadia-ANT+RF-test_metadata.csv -- the metadata file in GeoCSV format
        <br />* Cascadia-ANT+RF-test_metadata.json -- the metadata file in JSON format
6. Now we create the model files by running: 
   <br />**../../src/cvm_writer.py -m Cascadia-ANT+RF_prop.py -o Cascadia-ANT+RF-test -t metadata -d Cascadia-ANT+RF_data.txt -t netcdf**<br />
   *-d Cascadia-ANT+RF_data.txt* tells the code where it can find the model data (from step 3 above) 
   <br />*-t netcdf* tells it that we want the output to be in netCDF
   After running this command, we should see a new file *Cascadia-ANT+RF-test.nc* which is our model in netCDF.<br />We could have used *-t geocsv* to have the output as
   GeoCSV file or we could have combined them all and have *-t metadata,netcdf,geocsv* to get all the outputs.
7. Look at the CSV and JSON files to familiarize yourself with these files.
8. To look ate the netCDF file content, you will find a simple plotting script under the *src* directory (*src/simple_plotter.py*)
   Again, use *-h* flag to get the usage message. This simple plotter code reads the netCDF file, prints information on the data, extracts a depth slice and 
   displays it in geographic and UTM coordinates. It also extracts a small piece of the depth slice and plots it again. To run the plotter script against the netCDF file you created, run *../../src/simple_plotter.py -m Cascadia-ANT+RF-test.nc*. As you run this command, the model information is outputted to the terminal and a slice at the depth of 20 km is plotted. By closing this plot, another plot shows the same slice using the UTM coordinates. Finally, by closing the plot you get a plot of the upper right section of the area.
9. Look ate the model information printed to your terminal. It includes depth information. Select another depth (for example 76km). 
   Now run the same command but tell it to plot the slice for the depth of 78 km:
   <br />*../../src/simple_plotter.py -m Cascadia-ANT+RF-test.nc -d 78*
   <br />which produces the same plots as before but for the depth of 78 k

<br /><br />**cvm_subset**: Subset the model netCDF files:
1. cd to the model directory *sample-files/Cascadia-ANT+RF-Delph2018*
2. *../../src/cvm_slicer.py -i Cascadia-ANT+RF-Delph2018.r0.1.nc*
3. Download the netCDF file for casc1.6_velmdl from the [GoogleDrive](https://drive.google.com/drive/folders/1JTN0GAf25IIFBqkTMmZCTL50VnTjiiFd)
4. cd to the  model directory 
5. *../../src/cvm_slicer.py -i casc1.6-velmdl.r0.1.nc*

<br /><br />**netcdf_to_geocsv**:
 Use the netcdf_to_geocsv.py tool to convert the netCDF file to convert the *Cascadia-ANT+RF-test.nc* file from step 6 above, to GeoCSV and name the GeoCSV file conversion_test.csv::
 1. cd to the model directory *sample-files/Cascadia-ANT+RF-Delph2018*
 2. *../../src/netcdf_to_geocsv.py -i Cascadia-ANT+RF-test.nc -o conversion_test.csv*

## Package contents

This package contains the following files:

    README.md           -- This file
    requirements.txt    -- List of the Python modules to install

    documents/
        2024-01-08-CRESCENT-CVM.pptx -- Slide presentation of 2024-01-08

    lib/
        netcdf_to_geocsv_lib.py   -- Python library for netcdf_to_geocsv
        shared_lib.py   -- Python library file shared between tools
        slicer_lib.py   -- Python library file for cvm_slicer
        writer_lib.py   -- Python library file for cvm_writer

    prop/
        netcdf_to_geocsv_prop.py  -- Configuration file for netcdf_to_geocsv
        shared_prop.py  -- Configuration file shared between tools
        writer_prop.py  -- Configuration file for cvm_writer

    sample-files/
        Cascadia-ANT+RF-Delph2018/              -- The Cascadia-ANT+RF-Delph2018 model folder
            Cascadia-ANT+RF-Delph2018.r0.1.nc   -- Model in netCDF format
            Cascadia-ANT+RF_prop.py             -- Model metadata file
            Cascadia_ANT+RF_Delph2018.txt.gz    -- Model raw data
            simple_plotter_prop.py              -- The simple_plotter configuration file

        casc1.6_velmdl/             -- The casc1.6_velmdl model folder
            README                  -- A note on where to get the model raw data
            casc1.6-velmdl.r0.1.nc  -- Model in netCDF format
            casc16-velmdl_prop.py   -- Model metadata file
            simple_plotter_prop.py  -- The simple_plotter configuration file

    src/
        cvm_writer.py       -- The CVM writer tool
        cvm_slicer.py       -- The CVM model slicer tool
        netcdf_to_geocsv.py -- The netCDF to GeoCSV conversion tool
        simple_plotter.py   -- The simple_plotter tool

    template/
        metadata_template.py            -- Model metadata template
        metadata_template_detailed.py   -- Model metadata template with variables described

## Citation



### Comments, questions, or bug reports

  Please open an issue ticket or contact manochehr.bahavar@earthscope.org

