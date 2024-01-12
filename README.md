*Author: Ashish Kondal (ashish.kondal@wsu.edu)
*Date: September 10, 2022
*Description: Workflow to Generate Meteorological Input Forcings for VIC-CropSyst.
*Detailed Description: This workflow is for generating input forcings for VIC-CropSyst using gridMET observations and the North American Multi-Model Ensemble (NMME) forecast product.
  for a given year and initialization month, the NMME forecast goes for the next 8 months only and thus, lacks a whole year time series. This repository provides scripts for appending spin-up data to the NMME data, 
  re-gridding it to the VIC-CropSyst model resolution, and finally, generating input forcings for the VIC-CropSyst model.

*Note: The raw NMME data is at a coarse resolution and thus, requires spatiotemporal downscaling before using any of the following scripts.
*Note: for some scripts, you might have to write batch scripts with required input variables. 

1. Converting .mat format to ASCII
  Short Description: BiasCorrected NMME data was in .mat format (MATLAB native storage format) and needs to be converted into another format (such as ASCII) for subsequent processing.
  Details: To convert it into ASCII, use the "MAT_TO_ASCII.m" script.
	THIS SCRIPT CONVERTS MAT FORMAT TO ASCII AND GENERATES A SINGLE TEXT FILE FOR EACH INITIALIZATION WITH ALL YEARS FOR A BASIN.
 
2. Extract GridMET data (observation dataset) for a particular basin from the Kamiak Server. If you have downloaded the gridMET data for a particular basin only, skip this step.
	Details: Use "Extract_GridMET_fromKamiak_Parallel.py" script.
	THIS SCRIPT EXTRACT GRIDMET DATA FROM KAMIAK SERVER (or Folder where you've downloaded/stored the GridMET data). 
	Make sure to have the required packages installed in your Python library beforehand.

3. Convert ASCII files of NMME data from Step 1. to NetCDF and divide data into small grid blocks.
  Short Description: Dividing NMME and GridMET into small grid blocks (1*1degree) for parallel computation purposes.
   Details: Use "Subsetting_NMME_GridMET.py"

4. Merge NMME and GridMET data to make a continuous 4-5 Year time series.
  Description: Here we are appending 3-4 years of spin-up data (gridMET, in this case) depending upon the NMME's initialization month and storing the files in NetCDF format for subsequent processing.
  Details: Use "Merge_NMME_gridMET_NETCDF.py"  
	
5. Convert 1/24th resolution to 1/16th resolution. 
	Short Description: HERE, MERGED NMME AND GRIDMET DATA FROM STEP 4 IS REGRIDDED TO VIC-CROPSYST'S MODEL RESOLUTION.
  Details: Use "Regridding_NetCDF_To_VIC_CropSyst_res_parallel.py".
	
6. Create VIC-CropSyst Forcings.
	Short Description: THIS PYTHON SCRIPT TAKES DATA FROM STEP 4 AND  CREATE GRID-WISE FORCINGS THAT DIRECTLY BEING USED IN VIC-CROPSYST RUNS.
	Details: Use "Create_VIC_CropSyst_Forcings_Parallel.py"

NOTE: For better performance, you can convert the final forcings files into binary too.
	
	
