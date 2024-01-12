#=================================================
# Title: Extracting Daily GridMET data from KAMIAK 
# Date: August 31, 2022
# Author: Ashish Kondal (ashish.kondal@wsu.edu)
# Description: This script does the following: 
# Extracts gridmet daily data for any basin from HPCC based on Lat-Long provided via text file.


#=================================================
# IMPORTING REQUIRED PACKAGES
#=================================================

import os
import sys
import csv
import netCDF4 as nc
import pandas as pd
import numpy as np
import calendar
import datetime as dt
import itertools as itt
from tqdm import tqdm
import time


#=================================================
# DEFINE VARIABLES & INPUT SETTINGS
#=================================================

basin_id = sys.argv[1];		 # Retrieve basin_id from bash script. Used same numbering as mentioned in "MAT_TO_ASCII.m"
basin_id = int(basin_id)-1; 

varnum = sys.argv[2]		# Variable number; see "dispname" for variable list. 
varnum = int(varnum)-1

gridnum = sys.argv[3]		# Number of Gridblocks. For parallel computation, data was brokendown into chunks of blocks and here we are refering the number of block.
gridnum = int(gridnum)-1

lat_blocks = np.arange(39.5,50.5,1)	# Latitudinal extent
lon_blocks = np.arange(-125,-108,1)	# Longitudinal extent

cord_ind = pd.DataFrame(itt.product(range(len(lat_blocks)-1), range(len(lon_blocks)-1))) # Creating Grid based on selected "gridnum"
ii = cord_ind.iloc[gridnum,0]
jj = cord_ind.iloc[gridnum,1]

base_dir = 'workingdirectory/extracted_GridMET/'
input_directory = 'inputpath/GridMET/';
variables = ['pr_','rmax_','rmin_','srad_','tmmx_','tmmn_','vs_','sph_'];
dispname = ['precipitation_amount','max_relative_humidity','min_relative_humidity','surface_downwelling_shortwave_flux_in_air','max_air_temperature','min_air_temperature','wind_speed','specific_humidity'];
varnames = ['Lat','Lon','Year','Month','DOY',dispname[varnum]];
basins = ['okanogon','wallawalla','yakima','pnw'];

outpath = base_dir + dispname[varnum] + '/'
try:
    if os.path.exists(base_dir)== False:
        os.mkdir(base_dir)
    if os.path.exists(outpath)== False:
        os.mkdir(outpath)
except:
    pass

output_file = outpath + 'GridMET_Data_{0:.5f}_{1:.5f}'.format(lat_blocks[ii],lon_blocks[jj])

#================================================================
# READ A SINGLE GRIDMET FILE (ANYONE) TO GET LAT AND LON INFO
#================================================================

# Step 1: Read any GridMET NetCDF file.

testfilename = input_directory + variables[varnum] + str(1983) + '.nc'
try:
    testfile.close()
except:
    pass
    
testfile = nc.Dataset(testfilename,'r')

# Step 2: Save Latitude and Longitude in separate variables.
lat_np = np.array(testfile['lat'][:])
lon_np = np.array(testfile['lon'][:])
lat_df = pd.DataFrame(data = lat_np, columns = ["Lat"])
lon_df = pd.DataFrame(data = lon_np, columns = ["Lon"])

testfile.close()
del testfile
del lat_np
del lon_np


# Step 3: Extracting GridMET data for grid blocks.
lat_box = lat_df[((np.isclose(lat_df,lat_blocks[ii])) | (np.isclose(lat_df,lat_blocks[ii+1])) | ((lat_df > lat_blocks[ii]) & (lat_df < lat_blocks[ii+1])))]
lon_box = lon_df[((np.isclose(lon_df,lon_blocks[jj])) | (np.isclose(lon_df,lon_blocks[jj+1])) | ((lon_df > lon_blocks[jj]) & (lon_df < lon_blocks[jj+1])))]
lat_box = lat_box.dropna()
lon_box = lon_box.dropna()
lat_ind = lat_box.index
lon_ind = lon_box.index
  

for year in tqdm(range(1979,2022,1)):
    #================================================================
    # READ GridMET NETCDF FILES 
    #================================================================    
    var_filename = input_directory + variables[varnum] + str(year) + '.nc';
    
    if calendar.isleap(year)==True:
        mnth_filename = 'input_directory/Leap_Year_Month.csv' 		# Since the raw gridMET data doesn't have month column in it, we generated a column file with months refering to particular year.
        mnth = pd.read_csv(mnth_filename, header=None, index_col=None) 
        mnth = np.squeeze(mnth)
    else:
        mnth_filename = 'input_directory/NON_Leap_Year_Month.csv'
        mnth = pd.read_csv(mnth_filename, header=None, index_col=None) 
        mnth = np.squeeze(mnth)
 
    try:
        varfile.close()
    except:
        print("No NetCDF File Was Previously Opened")
    
       
    varfile = nc.Dataset(var_filename,mode = 'r')		# Read netcdf file of a variable
    print(var_filename)
    
    for i in lat_ind:
        clat = lat_df.loc[i]
        for j in lon_ind:
            clon = lon_df.loc[j]
            var = np.squeeze(np.array(varfile[Fullvariable_name[varnum]][:,i,j]))
            var[var > 30000] = np.nan				# missing values are replaced with NaN
            YR_c = np.squeeze(np.ones((len(var),1))*year)
            YR_c = YR_c.astype(int)
            Days = np.squeeze(np.arange(1,(len(var)+1)).reshape(len(var),1))
            Lat_c = np.squeeze(np.ones((len(var),1))*clat[0])
            Lon_c = np.squeeze(np.ones((len(var),1))*clon[0])
            data_final = pd.DataFrame({'Lat':Lat_c,'Lon':Lon_c,'Year':YR_c,'Month':mnth,'DOY':Days,Fullvariable_name[varnum]:var})
            if os.path.exists(output_file) == False:
                data_final.to_csv(output_file,index=False,header = varnames,na_rep = None)
            else:
                data_final.to_csv(output_file,mode ='a',index=False,header = False,na_rep = None)
            
            del data_final
            del var
            
    varfile.close()


