#=================================================================
# GENERATE VIC-CROPSYST FORCINGS FROM REGRIDDED NMME NETCDF DATA
# Author: Ashish Kondal
# Date: September 09, 2022
# Description: This script create grid-wise forcings that are directly used in VIC-CROPSYST Simulations.
#=================================================================

#================================================================
#   IMPORTING LIBRARIES/PACKAGES
#================================================================

from __future__ import print_function
import os
import shutil
import sys
import netCDF4
import numpy as np
import datetime as dt
import pandas as pd
import itertools as itt


#================================================================
#   DEFINE VARIABLES AND INPUT FILES
#================================================================

init_month = sys.argv[1];                                      # Retrieving NMME's initialization month from sbatch script
init_month = int(init_month);

y = sys.argv[2];                                               # Retrieving year from sbatch script 
y = int(y);         # starts from 1982

basin_id = sys.argv[3]                                         # Retrieving "basins" number (goes from 1-4; see "basins") from batch script
basin_id = int(basin_id)-1

gridnum = sys.argv[4]                                          # Retrieving gridblock number (goes from 1-160) from batch script
gridnum = int(gridnum)-1


# Latitudinal and Longitudinal Extent of Study Area
lat_blocks = np.arange(39.5,50.5,1)
lon_blocks = np.arange(-125,-108,1)


# Selecting coordinating based on current "gridnum"
cord_ind = pd.DataFrame(itt.product(range(len(lat_blocks)-1), range(len(lon_blocks)-1)))
ii = cord_ind.iloc[gridnum,0]
jj = cord_ind.iloc[gridnum,1]


fullvariable_name = ['precipitation_amount','max_air_temperature','min_air_temperature','wind_speed','specific_humidity','surface_downwelling_shortwave_flux_in_air','max_relative_humidity','min_relative_humidity'];
basins = ['okanogon','wallawalla','yakima','pnw']

base_dir = 'OutputDirectoryPath/VIC_CropSyst_Forcings/';
input_directory = 'Path_to_InputDirectory/' 

# Check if input directory exists or not
try:
    os.path.isdir(input_directory)
except:
    print("Regridded File Directory Doesn't Exist!")

# Input files directory
precip_dir = input_directory + fullvariable_name[0] +"/" + "Data_{0:.5f}_{1:.5f}".format(lat_blocks[ii],lon_blocks[jj]) + "/"
max_temp_dir = input_directory + fullvariable_name[1] +"/" + "Data_{0:.5f}_{1:.5f}".format(lat_blocks[ii],lon_blocks[jj]) + "/"
min_temp_dir = input_directory + fullvariable_name[2] +"/" + "Data_{0:.5f}_{1:.5f}".format(lat_blocks[ii],lon_blocks[jj]) + "/"
ws_dir = input_directory + fullvariable_name[3] +"/" + "Data_{0:.5f}_{1:.5f}".format(lat_blocks[ii],lon_blocks[jj]) + "/"
sph_dir = input_directory + fullvariable_name[4] +"/" + "Data_{0:.5f}_{1:.5f}".format(lat_blocks[ii],lon_blocks[jj]) + "/"
srd_dir = input_directory + fullvariable_name[5] +"/" + "Data_{0:.5f}_{1:.5f}".format(lat_blocks[ii],lon_blocks[jj]) + "/"
rhmax_dir = input_directory + fullvariable_name[6] +"/" + "Data_{0:.5f}_{1:.5f}".format(lat_blocks[ii],lon_blocks[jj]) + "/"
rhmin_dir = input_directory + fullvariable_name[7] +"/" + "Data_{0:.5f}_{1:.5f}".format(lat_blocks[ii],lon_blocks[jj]) + "/"

print(precip_dir)   # check to see if the file is named correctly or not!

start_year = y-3
if init_month >4:
    end_year = y+1
else:
    end_year = y
 
# Input file full file name - Regridded Data from previous step
pr_file = precip_dir + 'Regridded_' + 'Data_{0:.5f}_{1:.5f}'.format(lat_blocks[ii],lon_blocks[jj]) + '_Merged_VIC_Forcings_init_month_' + str(init_month) + '_' + fullvariable_name[0] + '_' + str(start_year) + '_' + str(end_year) + '.nc'
tmax_file = max_temp_dir + 'Regridded_' + 'Data_{0:.5f}_{1:.5f}'.format(lat_blocks[ii],lon_blocks[jj]) + '_Merged_VIC_Forcings_init_month_' + str(init_month) + '_' + fullvariable_name[1] + '_' + str(start_year) + '_' + str(end_year) + '.nc'
tmin_file = min_temp_dir + 'Regridded_' + 'Data_{0:.5f}_{1:.5f}'.format(lat_blocks[ii],lon_blocks[jj]) + '_Merged_VIC_Forcings_init_month_' + str(init_month) + '_' + fullvariable_name[2] + '_' + str(start_year) + '_' + str(end_year) + '.nc'
ws_file = ws_dir + 'Regridded_' + 'Data_{0:.5f}_{1:.5f}'.format(lat_blocks[ii],lon_blocks[jj]) + '_Merged_VIC_Forcings_init_month_' + str(init_month) + '_' + fullvariable_name[3] + '_' + str(start_year) + '_' + str(end_year) + '.nc'
sph_file = sph_dir + 'Regridded_' + 'Data_{0:.5f}_{1:.5f}'.format(lat_blocks[ii],lon_blocks[jj]) + '_Merged_VIC_Forcings_init_month_' + str(init_month) + '_' + fullvariable_name[4] + '_' + str(start_year) + '_' + str(end_year) + '.nc'
srd_file = srd_dir + 'Regridded_' + 'Data_{0:.5f}_{1:.5f}'.format(lat_blocks[ii],lon_blocks[jj]) + '_Merged_VIC_Forcings_init_month_' + str(init_month) + '_' + fullvariable_name[5] + '_' + str(start_year) + '_' + str(end_year) + '.nc'
rhmax_file = rhmax_dir + 'Regridded_' + 'Data_{0:.5f}_{1:.5f}'.format(lat_blocks[ii],lon_blocks[jj]) + '_Merged_VIC_Forcings_init_month_' + str(init_month) + '_' + fullvariable_name[6] + '_' + str(start_year) + '_' + str(end_year) + '.nc'
rhmin_file = rhmin_dir + 'Regridded_' + 'Data_{0:.5f}_{1:.5f}'.format(lat_blocks[ii],lon_blocks[jj]) + '_Merged_VIC_Forcings_init_month_' + str(init_month) + '_' + fullvariable_name[7] + '_' + str(start_year) + '_' + str(end_year) + '.nc'

# Loading Input Files 
try:
    pr_nc = netCDF4.Dataset(pr_file)
    tmax_nc = netCDF4.Dataset(tmax_file)
    tmin_nc = netCDF4.Dataset(tmin_file)
    ws_nc = netCDF4.Dataset(ws_file)
    sph_nc = netCDF4.Dataset(sph_file)
    srd_nc = netCDF4.Dataset(srd_file)
    rhmax_nc = netCDF4.Dataset(rhmax_file)
    rhmin_nc = netCDF4.Dataset(rhmin_file)
except:
    print("Failed to Load NetCDF Files")

# Check to see if all input files were correctly read or not!
print(pr_file)
print(tmax_file)
print(tmin_file)
print(ws_file)
print(sph_file)
print(srd_file)
print(rhmax_file)
print(rhmin_file)
print('end of files')

#================================================================
#   DEFINE OUTPUT DIRECTORIES & FILES
#================================================================

# Check if Output directory exists or not. If not, then create parent and sub-directories.
forcing_dir = base_dir + 'FORCINGS1/' + basins[basin_id] + '/InitMonth_0' + str(init_month) + '/'
try:
    if os.path.isdir(base_dir + 'FORCINGS1/')==False:
        os.mkdir(base_dir + 'FORCINGS1/')
    if os.path.isdir(base_dir + 'FORCINGS1/' + basins[basin_id] + '/')==False:
        os.mkdir(base_dir + 'FORCINGS1/' + basins[basin_id] + '/' )
    if os.path.isdir(forcing_dir) == False:
        os.mkdir(forcing_dir)
except:
    pass

# Creating/Retreiving full output file names
forcing_yr_dir1 =  forcing_dir + str(start_year) + '_' + str(end_year) + '/'
forcing_yr_dir = forcing_dir + str(start_year) + '_' + str(end_year) + '/' + 'Data_{0:.5f}_{1:.5f}'.format(lat_blocks[ii],lon_blocks[jj]) +'/'

# check if the output directory exists. If yes, then delete old files and create new directory, if not.
if os.path.isdir(forcing_yr_dir1)==False:
    os.mkdir(forcing_yr_dir1)
    
if os.path.isdir(forcing_yr_dir)==True:
    shutil.rmtree(forcing_yr_dir)
    
if os.path.isdir(forcing_yr_dir)==False:
    os.mkdir(forcing_yr_dir)
   
#================================================================
#   PREALLOCATING  MEMORY TO THE OUTPUT VARIABLES
#================================================================
   
ctr = 0
R = pr_nc['latitude'][:]                # Retrieve latitudes from one of the input files (here i chose precipitation file). You can chose any file as this infomation is same in all files. 
C = pr_nc['longitude'][:]               # Retrieve longitudes from one of the input files 
time = pr_nc['time'].shape[0]           # Retrieve time information from one of the input files 

# Creating empty files and index for later writing
file_ids = []
for i in range(len(R)):
    for j in range(len(C)):
        meteofile = forcing_yr_dir + 'data_{0:.5f}_{1:.5f}'.format(R[i],C[j])
        file_ids.append(open(meteofile,'w'))

#================================================================
#   GENERATING VIC-CROPSYST INPUT FORCINGS
#================================================================
   

for i in range(len(R)):
    for j in range(len(C)):
        cur_fid = file_ids[ctr]
        for t in range(time):
            #cur_fid = file_ids[ctr]
            #print(cur_fid)
            cur_pr = pr_nc.variables['precipitation_amount'][t,i,j]
            cur_tmax = tmax_nc.variables['max_air_temperature'][t,i,j] - 273.15
            cur_tmin = tmin_nc.variables['min_air_temperature'][t,i,j] - 273.15
            cur_ws = ws_nc.variables['wind_speed'][t,i,j]
            cur_sph = sph_nc.variables['specific_humidity'][t,i,j]
            cur_srd = srd_nc.variables['surface_downwelling_shortwave_flux_in_air'][t,i,j]
            cur_rhmax = rhmax_nc.variables['max_relative_humidity'][t,i,j]
            cur_rhmin = rhmin_nc.variables['min_relative_humidity'][t,i,j]
            cur_fid.write('{0} {1} {2} {3} {4} {5} {6} {7} \n'.format(cur_pr,cur_tmax,cur_tmin,cur_ws,cur_sph,cur_srd,cur_rhmax,cur_rhmin))
            #cur_fid.dropna()
        ctr = ctr + 1
        
try:
    pr_nc.close()
    tmax_nc.close()
    tmin_nc.close()
    ws_nc.close()
    sph_nc.close()
    srd_nc.close()
    rhmax_nc.close()
    rhmin_nc.close()
except:
    pass
    
# Saving and closing the output files
ctr = 0
for i in range(len(R)):
    for j in range(len(C)):
        cur_fid = file_ids[ctr]
        cur_fid.close()
        ctr = ctr + 1

print("Finished")