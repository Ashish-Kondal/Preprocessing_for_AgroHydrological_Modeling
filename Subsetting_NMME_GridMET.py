#=================================================================
# Creating ASCII Files of NMME and GridMET data at 1 * 1 Resolution
# Author: Ashish Kondal
# Date: October 14, 2022
# Description: This scripts divide GridMET and NMME data into smaller blocks
#             
#=================================================================


#================================================================
#   IMPORTING LIBRARIES/PACKAGES
#================================================================

from __future__ import print_function
import sys
import os, string
import datetime as dt
import netCDF4 as nc
import pandas as pd
import numpy as np
import itertools as itt
from tqdm import tqdm
import concurrent.futures
import time
from multiprocessing import Pool


#================================================================
#   DEFINE VARIABLES AND INPUT FILES
#================================================================


input_directory = 'ADD_PATH_TO_YOUR_FILES_GENERATED_IN_STEP1';
variables = ['_pr','_rhsmax','_rhsmin','_srad','_tasmax','_tasmin','_vs','_sph'];           # Suffix of variables required to run VIC-Cropsyst
fullvariable_name = ['precipitation_amount','max_air_temperature','min_air_temperature','wind_speed','specific_humidity','surface_downwelling_shortwave_flux_in_air','max_relative_humidity','min_relative_humidity'];
basins = ['okanogon','wallawalla','yakima','pnw']                                       

init_month = sys.argv[1];                                       # Retrieving NMME's initialization month from sbatch script
init_month = int(init_month);
f'initialization Month = {init_month}'

varnum = sys.argv[2]                                            # Retrieving "variables" number (goes from 1-8) from sbatch script
varnum = int(varnum)-1
f'Variable Number = {varnum}'

basin_id = sys.argv[3]                                          # Retrieving "basins" number (goes from 1-4; see "basins") from batch script
basin_id = int(basin_id)-1
f'basin = {basins[basin_id]}'


nmme_header = ['Lat','Lon','Initialization_Year','Initialization_Month','Forecast_Year','Forecast_Month','Forecast_Day','Precipitation_mm', 'Max_Temperature_K', 'Min_Temperature_K', 'Wind_Speed_mps','Specific_Humidity_Kg_per_Kg','ShortWave_Radiation_W_per_m2','Max_Relative_Humidity','Min_Relative_Humidity'];
gridmet_header = ['Lat','Lon','Year','DOY','Precipitation_mm', 'Max_Temperature_K', 'Min_Temperature_K', 'Wind_Speed_mps','Specific_Humidity_Kg_per_Kg','ShortWave_Radiation_W_per_m2','Max_Relative_Humidity','Min_Relative_Humidity'];
req_cols_nmme = ['Lat','Lon','Initialization_Year','Initialization_Month','Forecast_Year','Forecast_Month','Forecast_Day',nmme_header[varnum+7]]
req_cols_gmet = ['Lat','Lon','Year','DOY',gridmet_header[varnum+4]]



# Latitudinal and Longitudinal Extent of Study Area
lat_blocks = np.arange(39.5,50.5,1)
lon_blocks = np.arange(-125,-108,1)

# Reading Input File
nmme_filename = input_directory + 'BCSD_NMME_data_' + basins[basin_id] + '_ExtractedForInitializationMonth_' + str(init_month) +'.txt' 
nmme_data = pd.read_csv(nmme_filename, header=0, index_col=None, usecols = req_cols_nmme)  
outpath = input_directory + "VIC_CropSyst_Forcings/Subsetting_SingleGrids_24thResNMME/"  
filepath = outpath +'/InitMonth_0' + str(init_month) + '/' + fullvariable_name[varnum]+ '/'

# Check if output directory exists or not. If yes, then delete old files to remove any chance of overwritting and If output directory doesn't exists, then create one.
try:
    if os.path.exists(outpath) == False:
        os.mkdir(outpath)
    if os.path.isdir(outpath + '/InitMonth_0' + str(init_month)) == False:
        os.mkdir(outpath + '/InitMonth_0' + str(init_month)+"/")
    if os.path.exists(filepath) == False:    
        os.mkdir(filepath)
    else:
        try:
            for f in os.listdir(filepath):
                os.remove(os.path.join(filepath,f))
        except IOError:
            raise IOError("Cannot Delete Previously Generated Files")
except:
    pass
    
#===============================================================================================
#  NMME: Creating ASCII Files for each grid block
#===============================================================================================


for i,j in tqdm(itt.product(range(len(lat_blocks)-1), range(len(lon_blocks)-1))):
    outfilename = filepath + "NMME_Data_{0:.5f}_{1:.5f}.txt".format(lat_blocks[i],lon_blocks[j])
    nmme_grid = nmme_data[(((np.isclose(nmme_data.Lat,lat_blocks[i])) | (np.isclose(nmme_data.Lat,lat_blocks[i+1])) | ((nmme_data.Lat > lat_blocks[i]) & (nmme_data.Lat < lat_blocks[i+1]))) & ((np.isclose(nmme_data.Lon,lon_blocks[j])) | (np.isclose(nmme_data.Lon,lon_blocks[j+1])) | ((nmme_data.Lon > lon_blocks[j]) & (nmme_data.Lon < lon_blocks[j+1]))))]
    if nmme_grid.empty:
        pass
    else:
        nmme_grid.to_csv(outfilename,index=False,header = req_cols_nmme)
        del nmme_grid


#===============================================================================================
#   GridMET: Creating ASCII Files for each grid block
#  Hardwire to run only when init_month == 1, change if you like
#===============================================================================================

if init_month == 1:
    gridMET_filename = input_directory + 'GridMET_Files/' + 'GridMET_data_' + basins[basin_id] +'.txt'
    gmet_data = pd.read_csv(gridMET_filename, header=0, index_col=None, usecols = req_cols_gmet)

    #gmet_data.info(verbose=False,memory_usage="deep")  #check memory_usage
    outpath = input_directory + "VIC_CropSyst_Forcings/Subsetting_SingleGrids_24thRes/GridMET/"
    filepath = outpath + fullvariable_name[varnum]+ '/'
    
    # Check if the output folder exists or not.
    try:
        if os.path.exists(outpath) == False:
            os.mkdir(outpath)
        if os.path.exists(filepath) == False:    
            os.mkdir(filepath)
    except:
        pass
    
    # function to write files
    def write_to_file(i,j,gmet_data):
        outfilename = filepath + "GridMET_Data_{0:.5f}_{1:.5f}.txt".format(lat_blocks[i],lon_blocks[j])
        gmet_grid = gmet_data[((np.isclose(gmet_data.Lat,lat_blocks[i])) | (np.isclose(gmet_data.Lat,lat_blocks[i+1])) | ((gmet_data.Lat > lat_blocks[i]) & (gmet_data.Lat < lat_blocks[i+1]))) & ((np.isclose(gmet_data.Lon,lon_blocks[j])) | (np.isclose(gmet_data.Lon,lon_blocks[j+1])) | ((gmet_data.Lon > lon_blocks[j]) & (gmet_data.Lon < lon_blocks[j+1])))]
        if gmet_grid.empty:
            pass
        else:
            gmet_grid.to_csv(outfilename,index=False,header = req_cols_gmet)

    # parallel computation of subsetting and saving gridMET data
    tic = time.perf_counter()
    with concurrent.futures.ThreadPoolExecutor() as executor:
        [executor.submit(write_to_file,i,j,gmet_data) for i,j in tqdm(itt.product(range(len(lat_blocks)-1), range(len(lon_blocks)-1)))];

    toc = time.perf_counter()
    print(f"GridMET: {toc - tic:0.4f} seconds")   
    del gmet_data
   