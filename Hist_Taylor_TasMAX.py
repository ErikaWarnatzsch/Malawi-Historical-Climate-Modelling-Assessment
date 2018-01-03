"""
Created on Wednesday August 2 2017

@author: s0899345

Taylor Diagram (Taylor, 2001) code from "Yannick Copin <yannick.copin@laposte.net>". This code was put in to the public domain and accessed from https://gist.github.com/ycopin/3342888
with "Time-stamp: <2017-11-20 17:36:41 ycopin>" 

"""

import iris
import iris.coord_categorisation as iriscc
import iris.analysis.cartography
import numpy as np
import matplotlib.pyplot as PLT
from Taylor_Diagram_Template_YC import TaylorDiagram
from cf_units import Unit
import math

#this file is split into parts as follows:
    #PART A: Define Necessary Maths
    #PART B: Main Script
    #PART 1: load and format all models 
    #PART 2: load and format observed data
    #PART 3: format files to synthesize data temporaly and spatially
    #PART 4: Calculate Standard Deviation
    #PART 5: Calculate Pearson's Correlation Coefficient
    #PART 6: Plot Taylor Diagram  (based on Copin, 2017)
    
#PART A : Define Necessary Maths
def mean(x):
    sum = 0.0
    for i in x:
         sum += i
    return sum / len(x) 

# calculates the sample standard deviation
def sampleStandardDeviation(x):
    sumv = 0.0
    for i in x:
         sumv += (i - mean(x))**2
    return math.sqrt(sumv/(len(x)-1))

# calculates the PCC using both the 2 functions above
def pearson(x,y):
    scorex = []
    scorey = []

    for i in x: 
        scorex.append((i - mean(x))/sampleStandardDeviation(x)) 

    for j in y:
        scorey.append((j - mean(y))/sampleStandardDeviation(y))

# multiplies both lists together into 1 list (hence zip) and sums the whole list   
    return (sum([i*j for i,j in zip(scorex,scorey)]))/(len(x)-1)
            
#PART B: Main Script
def main():
    #-------------------------------------------------------------------------
    #PART 1: LOAD and FORMAT ALL MODELS   
    #bring in all the ERAINT models we need and give them a name
    CCCma = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/ERAINT/1979-2012/tasmax_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_CCCma-CanRCM4_r2_mon_198901-200912.nc'
    CLMcom = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/ERAINT/1979-2012/tasmax_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_198901-200812.nc'
    DMI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/ERAINT/1979-2012/tasmax_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_DMI-HIRHAM5_v2_mon_198901-201012.nc'
    KNMI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/ERAINT/1979-2012/tasmax_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_KNMI-RACMO22T_v1_mon_197901-201212.nc'
    MPIE = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/ERAINT/1979-2012/tasmax_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_MPI-CSC-REMO2009_v1_mon_198902-200812.nc'
    SMHI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/ERAINT/1979-2012/tasmax_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_SMHI-RCA4_v1_mon_198001-201012.nc'
    
    #bring in all the CORDEX RCM models we need and give them a name
    CCCmaCanRCM= '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/Historical/1950-2005/tasmax_AFR-44_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_mon_195001-200512.nc'
    CCCmaSMHI= '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/Historical/1950-2005/tasmax_AFR-44_CCCma-CanESM2_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'   
    CNRM= '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/Historical/1950-2005/tasmax_AFR-44_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_195001-200512.nc'
    CNRMSMHI= '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/Historical/1950-2005/tasmax_AFR-44_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    CSIRO = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/Historical/1950-2005/tasmax_AFR-44_CSIRO-QCCCE-CSIRO-Mk3-6-0_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    ICHECDMI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/Historical/1950-2005/tasmax_AFR-44_ICHEC-EC-EARTH_historical_r3i1p1_DMI-HIRHAM5_v2_mon_195101-200512.nc'   
    ICHECCCLM = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/Historical/1950-2005/tasmax_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_CLMcom-CCLM4-8-17_v1_mon_194912-200512.nc'    
    ICHECKNMI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/Historical/1950-2005/tasmax_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_KNMI-RACMO22T_v1_mon_195001-200512.nc'
    ICHECMPI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/Historical/1950-2005/tasmax_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_MPI-CSC-REMO2009_v1_mon_195001-200512.nc'
    ICHECSMHI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/Historical/1950-2005/tasmax_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    IPSL = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/Historical/1950-2005/tasmax_AFR-44_IPSL-IPSL-CM5A-MR_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    MIROC =  '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/Historical/1950-2005/tasmax_AFR-44_MIROC-MIROC5_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc' 
    MOHCCCLM = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/Historical/1950-2005/tasmax_AFR-44_MOHC-HadGEM2-ES_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_194912-200512.nc' 
    MOHCKNMI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/Historical/1950-2005/tasmax_AFR-44_MOHC-HadGEM2-ES_historical_r1i1p1_KNMI-RACMO22T_v2_mon_195001-200512.nc'
    MOHCSMHI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/Historical/1950-2005/tasmax_AFR-44_MOHC-HadGEM2-ES_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    MPICCLM = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/Historical/1950-2005/tasmax_AFR-44_MPI-M-MPI-ESM-LR_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_194912-200512.nc'     
    MPIREMO = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/Historical/1950-2005/tasmax_AFR-44_MPI-M-MPI-ESM-LR_historical_r1i1p1_MPI-CSC-REMO2009_v1_mon_195001-200512.nc'    
    MPISMHI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/Historical/1950-2005/tasmax_AFR-44_MPI-M-MPI-ESM-LR_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    NCCDMI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/Historical/1950-2005/tasmax_AFR-44_NCC-NorESM1-M_historical_r1i1p1_DMI-HIRHAM5_v1_mon_195101-200512.nc'    
    NCCSMHI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/Historical/1950-2005/tasmax_AFR-44_NCC-NorESM1-M_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'   
    NOAA = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tasmax/Historical/1950-2005/tasmax_AFR-44_NOAA-GFDL-GFDL-ESM2M_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'      
    
    #bring in all the GCM models we need and give them a name
    CanESM2= '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/GCM_data/tasmax_Amon_CanESM2_historical_r1i1p1_185001-200512.nc'
    CNRMG= '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/GCM_data/tasmax_Amon_CNRM-CM5-2_historical_r1i1p1_195001-200512.nc'
    MK3 = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/GCM_data/tasmax_Amon_CSIRO-Mk3-6-0_historical_r1i1p1_185001-200512.nc'    
    EARTH= '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/GCM_data/tasmax_Amon_EC-EARTH_historical_r12i1p1_195001-201212.nc'
    EARTH3= '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/GCM_data/tasmax_Amon_EC-EARTH_historical+r3i1p1_1961-2005.nc'
    GFDL='/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/GCM_data/tasmax_Amon_GFDL-ESM2M_historical_r1i1p1_194601-200512.nc'
    HadGEM2='/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/GCM_data/tasmax_Amon_HadGEM2-ES_historical_r1i1p1_193412-200511.nc'
    IPSLG='/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/GCM_data/tasmax_Amon_IPSL-CM5B-LR_historical_r1i1p1_185001-200512.nc'
    MIROCG='/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/GCM_data/tasmax_Amon_MIROC5_historical_r1i1p1_185001-201212.nc'
    MPI='/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/GCM_data/tasmax_Amon_MPI-ESM-MR_historical_r1i1p1_185001-200512.nc'
    NorESM1='/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/GCM_data/tasmax_Amon_NorESM1-M_historical_r1i1p1_185001-200512.nc'
    
    #Load exactly one cube from given file
    CCCma = iris.load_cube(CCCma)
    CLMcom = iris.load_cube(CLMcom)
    DMI = iris.load_cube(DMI, 'air_temperature')
    KNMI = iris.load_cube(KNMI)
    MPIE = iris.load_cube(MPIE)
    SMHI = iris.load_cube(SMHI)
    
    CCCmaCanRCM = iris.load_cube(CCCmaCanRCM)
    CCCmaSMHI = iris.load_cube(CCCmaSMHI)
    CNRM = iris.load_cube(CNRM)
    CNRMSMHI = iris.load_cube(CNRMSMHI)
    CSIRO = iris.load_cube(CSIRO)
    ICHECDMI = iris.load_cube(ICHECDMI, 'air_temperature')
    ICHECCCLM = iris.load_cube(ICHECCCLM)
    ICHECKNMI = iris.load_cube(ICHECKNMI)
    ICHECMPI = iris.load_cube(ICHECMPI)
    ICHECSMHI = iris.load_cube(ICHECSMHI)
    IPSL = iris.load_cube(IPSL)
    MIROC = iris.load_cube(MIROC)
    MOHCCCLM = iris.load_cube(MOHCCCLM)
    MOHCKNMI = iris.load_cube(MOHCKNMI)
    MOHCSMHI = iris.load_cube(MOHCSMHI)
    MPICCLM = iris.load_cube(MPICCLM)
    MPIREMO = iris.load_cube(MPIREMO)
    MPISMHI = iris.load_cube(MPISMHI)
    NCCDMI = iris.load_cube(NCCDMI, 'air_temperature')
    NCCSMHI = iris.load_cube(NCCSMHI)
    NOAA = iris.load_cube(NOAA)
    
    CanESM2 = iris.load_cube(CanESM2)
    CNRMG = iris.load_cube(CNRMG)
    MK3 = iris.load_cube(MK3)
    EARTH = iris.load_cube(EARTH)
    EARTH3 = iris.load_cube(EARTH3)
    GFDL = iris.load_cube(GFDL, 'air_temperature')
    HadGEM2 = iris.load_cube(HadGEM2)
    IPSLG = iris.load_cube(IPSLG)
    MIROCG = iris.load_cube(MIROCG)
    MPI = iris.load_cube(MPI)
    NorESM1 = iris.load_cube(NorESM1)
    
    #remove flat latitude and longitude and only use grid latitude and grid longitude to make consistent with the observed data, also make sure all of the longitudes are monotonic. This is only applicable to the ERAINT and CORDEX Models as the Observed data and GCMs are already in standard linear format.    
    lats = iris.coords.DimCoord(CCCma.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CCCma.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
                                
    CCCma.remove_coord('latitude')
    CCCma.remove_coord('longitude')
    CCCma.remove_coord('grid_latitude')
    CCCma.remove_coord('grid_longitude')
    CCCma.add_dim_coord(lats, 1)
    CCCma.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(CLMcom.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CLMcom.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    CLMcom.remove_coord('latitude')
    CLMcom.remove_coord('longitude')
    CLMcom.remove_coord('grid_latitude')
    CLMcom.remove_coord('grid_longitude')
    CLMcom.add_dim_coord(lats, 1)
    CLMcom.add_dim_coord(lons, 2)  
    
    lats = iris.coords.DimCoord(DMI.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = DMI.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    DMI.remove_coord('latitude')
    DMI.remove_coord('longitude')
    DMI.remove_coord('grid_latitude')
    DMI.remove_coord('grid_longitude')
    DMI.add_dim_coord(lats, 1)
    DMI.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(KNMI.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = KNMI.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
        
    KNMI.remove_coord('latitude')
    KNMI.remove_coord('longitude')
    KNMI.remove_coord('grid_latitude')
    KNMI.remove_coord('grid_longitude')
    KNMI.add_dim_coord(lats, 1)
    KNMI.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(MPIE.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MPIE.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
   
    MPIE.remove_coord('latitude')
    MPIE.remove_coord('longitude')
    MPIE.remove_coord('grid_latitude')
    MPIE.remove_coord('grid_longitude')
    MPIE.add_dim_coord(lats, 1)
    MPIE.add_dim_coord(lons, 2)

    lats = iris.coords.DimCoord(SMHI.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = SMHI.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    SMHI.remove_coord('latitude')
    SMHI.remove_coord('longitude')
    SMHI.remove_coord('grid_latitude')
    SMHI.remove_coord('grid_longitude')
    SMHI.add_dim_coord(lats, 1)
    SMHI.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(CCCmaCanRCM.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CCCmaCanRCM.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
                              
    CCCmaCanRCM.remove_coord('latitude')
    CCCmaCanRCM.remove_coord('longitude')
    CCCmaCanRCM.remove_coord('grid_latitude')
    CCCmaCanRCM.remove_coord('grid_longitude')
    CCCmaCanRCM.add_dim_coord(lats, 1)
    CCCmaCanRCM.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(CCCmaSMHI.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CCCmaSMHI.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    CCCmaSMHI.remove_coord('latitude')
    CCCmaSMHI.remove_coord('longitude')
    CCCmaSMHI.remove_coord('grid_latitude')
    CCCmaSMHI.remove_coord('grid_longitude')
    CCCmaSMHI.add_dim_coord(lats, 1)
    CCCmaSMHI.add_dim_coord(lons, 2)  
    
    lats = iris.coords.DimCoord(CNRM.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CNRM.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
                                
    CNRM.remove_coord('latitude')
    CNRM.remove_coord('longitude')
    CNRM.remove_coord('grid_latitude')
    CNRM.remove_coord('grid_longitude')
    CNRM.add_dim_coord(lats, 1)
    CNRM.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(CNRMSMHI.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CNRMSMHI.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    CNRMSMHI.remove_coord('latitude')
    CNRMSMHI.remove_coord('longitude')
    CNRMSMHI.remove_coord('grid_latitude')
    CNRMSMHI.remove_coord('grid_longitude')
    CNRMSMHI.add_dim_coord(lats, 1)
    CNRMSMHI.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(CSIRO.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CSIRO.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    CSIRO.remove_coord('latitude')
    CSIRO.remove_coord('longitude')
    CSIRO.remove_coord('grid_latitude')
    CSIRO.remove_coord('grid_longitude')
    CSIRO.add_dim_coord(lats, 1)
    CSIRO.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(ICHECDMI.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECDMI.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')    
    
    ICHECDMI.remove_coord('latitude')
    ICHECDMI.remove_coord('longitude')
    ICHECDMI.remove_coord('grid_latitude')
    ICHECDMI.remove_coord('grid_longitude')
    ICHECDMI.add_dim_coord(lats, 1)
    ICHECDMI.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(ICHECCCLM.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECCCLM.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees') 
    
    ICHECCCLM.remove_coord('latitude')
    ICHECCCLM.remove_coord('longitude')
    ICHECCCLM.remove_coord('grid_latitude')
    ICHECCCLM.remove_coord('grid_longitude')
    ICHECCCLM.add_dim_coord(lats, 1)
    ICHECCCLM.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(ICHECKNMI.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECKNMI.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees') 
    
    ICHECKNMI.remove_coord('latitude')
    ICHECKNMI.remove_coord('longitude')
    ICHECKNMI.remove_coord('grid_latitude')
    ICHECKNMI.remove_coord('grid_longitude')
    ICHECKNMI.add_dim_coord(lats, 1)
    ICHECKNMI.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(ICHECMPI.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECMPI.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees') 
    
    ICHECMPI.remove_coord('latitude')
    ICHECMPI.remove_coord('longitude')
    ICHECMPI.remove_coord('grid_latitude')
    ICHECMPI.remove_coord('grid_longitude')
    ICHECMPI.add_dim_coord(lats, 1)
    ICHECMPI.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(ICHECSMHI.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECSMHI.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    ICHECSMHI.remove_coord('latitude')
    ICHECSMHI.remove_coord('longitude')
    ICHECSMHI.remove_coord('grid_latitude')
    ICHECSMHI.remove_coord('grid_longitude')
    ICHECSMHI.add_dim_coord(lats, 1)
    ICHECSMHI.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(IPSL.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = IPSL.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    IPSL.remove_coord('latitude')
    IPSL.remove_coord('longitude')
    IPSL.remove_coord('grid_latitude')
    IPSL.remove_coord('grid_longitude')
    IPSL.add_dim_coord(lats, 1)
    IPSL.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MIROC.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MIROC.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
                                
    MIROC.remove_coord('latitude')
    MIROC.remove_coord('longitude')
    MIROC.remove_coord('grid_latitude')
    MIROC.remove_coord('grid_longitude')
    MIROC.add_dim_coord(lats, 1)
    MIROC.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MOHCCCLM.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MOHCCCLM.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MOHCCCLM.remove_coord('latitude')
    MOHCCCLM.remove_coord('longitude')
    MOHCCCLM.remove_coord('grid_latitude')
    MOHCCCLM.remove_coord('grid_longitude')
    MOHCCCLM.add_dim_coord(lats, 1)
    MOHCCCLM.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MOHCKNMI.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MOHCKNMI.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MOHCKNMI.remove_coord('latitude')
    MOHCKNMI.remove_coord('longitude')
    MOHCKNMI.remove_coord('grid_latitude')
    MOHCKNMI.remove_coord('grid_longitude')
    MOHCKNMI.add_dim_coord(lats, 1)
    MOHCKNMI.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MOHCSMHI.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MOHCSMHI.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MOHCSMHI.remove_coord('latitude')
    MOHCSMHI.remove_coord('longitude')
    MOHCSMHI.remove_coord('grid_latitude')
    MOHCSMHI.remove_coord('grid_longitude')
    MOHCSMHI.add_dim_coord(lats, 1)
    MOHCSMHI.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MPICCLM.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MPICCLM.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MPICCLM.remove_coord('latitude')
    MPICCLM.remove_coord('longitude')
    MPICCLM.remove_coord('grid_latitude')
    MPICCLM.remove_coord('grid_longitude')
    MPICCLM.add_dim_coord(lats, 1)
    MPICCLM.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MPIREMO.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MPIREMO.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MPIREMO.remove_coord('latitude')
    MPIREMO.remove_coord('longitude')
    MPIREMO.remove_coord('grid_latitude')
    MPIREMO.remove_coord('grid_longitude')
    MPIREMO.add_dim_coord(lats, 1)
    MPIREMO.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MPISMHI.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MPISMHI.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MPISMHI.remove_coord('latitude')
    MPISMHI.remove_coord('longitude')
    MPISMHI.remove_coord('grid_latitude')
    MPISMHI.remove_coord('grid_longitude')
    MPISMHI.add_dim_coord(lats, 1)
    MPISMHI.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(NCCDMI.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = NCCDMI.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    NCCDMI.remove_coord('latitude')
    NCCDMI.remove_coord('longitude')
    NCCDMI.remove_coord('grid_latitude')
    NCCDMI.remove_coord('grid_longitude')
    NCCDMI.add_dim_coord(lats, 1)
    NCCDMI.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(NCCSMHI.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = NCCSMHI.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    NCCSMHI.remove_coord('latitude')
    NCCSMHI.remove_coord('longitude')
    NCCSMHI.remove_coord('grid_latitude')
    NCCSMHI.remove_coord('grid_longitude')
    NCCSMHI.add_dim_coord(lats, 1)
    NCCSMHI.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(NOAA.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = NOAA.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    NOAA.remove_coord('latitude')
    NOAA.remove_coord('longitude')
    NOAA.remove_coord('grid_latitude')
    NOAA.remove_coord('grid_longitude')
    NOAA.add_dim_coord(lats, 1)
    NOAA.add_dim_coord(lons, 2)
    
    #-------------------------------------------------------------------------
    #PART 2: OBSERVED DATA
    #bring in all the files we need and give them a name
    CRU= '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/Actual_Data/cru_ts4.01.1901.2016.tmx.dat.nc'
        
    #Load exactly one cube from given file
    CRU = iris.load_cube(CRU, 'near-surface temperature maximum')
    
    #-------------------------------------------------------------------------
    #PART 3: FORMAT DATA TO BE PLOT SPECIFIC 
    #regrid all models to have same latitude and longitude system, all regridded to model with lowest resolution
    CCCma = CCCma.regrid(CanESM2, iris.analysis.Linear())
    CLMcom =CLMcom.regrid(CanESM2, iris.analysis.Linear())
    DMI=DMI.regrid(CanESM2, iris.analysis.Linear())
    KNMI=KNMI.regrid(CanESM2, iris.analysis.Linear())
    MPIE=MPIE.regrid(CanESM2, iris.analysis.Linear())
    SMHI=SMHI.regrid(CanESM2, iris.analysis.Linear())
    
    CCCmaCanRCM = CCCmaCanRCM.regrid(CanESM2, iris.analysis.Linear())
    CCCmaSMHI = CCCmaSMHI.regrid(CanESM2, iris.analysis.Linear())
    CNRM =CNRM.regrid(CanESM2, iris.analysis.Linear())
    CNRMSMHI =CNRMSMHI.regrid(CanESM2, iris.analysis.Linear())
    CSIRO=CSIRO.regrid(CanESM2, iris.analysis.Linear())
    ICHECDMI=ICHECDMI.regrid(CanESM2, iris.analysis.Linear())
    ICHECCCLM=ICHECCCLM.regrid(CanESM2, iris.analysis.Linear())
    ICHECKNMI=ICHECKNMI.regrid(CanESM2, iris.analysis.Linear())
    ICHECMPI=ICHECMPI.regrid(CanESM2, iris.analysis.Linear())
    ICHECSMHI=ICHECSMHI.regrid(CanESM2, iris.analysis.Linear())
    IPSL=IPSL.regrid(CanESM2, iris.analysis.Linear())
    MIROC=MIROC.regrid(CanESM2, iris.analysis.Linear())
    MOHCCCLM=MOHCCCLM.regrid(CanESM2, iris.analysis.Linear())
    MOHCKNMI=MOHCKNMI.regrid(CanESM2, iris.analysis.Linear())
    MOHCSMHI=MOHCSMHI.regrid(CanESM2, iris.analysis.Linear())
    MPICCLM=MPICCLM.regrid(CanESM2, iris.analysis.Linear())
    MPIREMO=MPIREMO.regrid(CanESM2, iris.analysis.Linear())
    MPISMHI=MPISMHI.regrid(CanESM2, iris.analysis.Linear())
    NCCDMI=NCCDMI.regrid(CanESM2, iris.analysis.Linear())
    NCCSMHI=NCCSMHI.regrid(CanESM2, iris.analysis.Linear())
    NOAA=NOAA.regrid(CanESM2, iris.analysis.Linear())
    
    #CanESM2 = CanESM2.regrid(CanESM2, iris.analysis.Linear())
    CNRMG =CNRMG.regrid(CanESM2, iris.analysis.Linear())
    MK3 =MK3.regrid(CanESM2, iris.analysis.Linear())
    EARTH =EARTH.regrid(CanESM2, iris.analysis.Linear())
    EARTH3 =EARTH3.regrid(CanESM2, iris.analysis.Linear())
    GFDL=GFDL.regrid(CanESM2, iris.analysis.Linear())
    HadGEM2=HadGEM2.regrid(CanESM2, iris.analysis.Linear())
    IPSLG=IPSLG.regrid(CanESM2, iris.analysis.Linear())
    MIROCG=MIROCG.regrid(CanESM2, iris.analysis.Linear())
    MPI=MPI.regrid(CanESM2, iris.analysis.Linear())
    NorESM1=NorESM1.regrid(CanESM2, iris.analysis.Linear())
    
    CRU = CRU.regrid(CanESM2, iris.analysis.Linear())
    
    #we are only interested in the latitude and longitude relevant to Malawi (has to be slightly larger than country boundary to take into account resolution of GCMs)
    Malawi = iris.Constraint(longitude=lambda v: 32.0 <= v <= 36., latitude=lambda v: -17. <= v <= -8.)   
    
    CCCma = CCCma.extract(Malawi)
    CLMcom =CLMcom.extract(Malawi)
    DMI=DMI.extract(Malawi)
    KNMI=KNMI.extract(Malawi)
    MPIE=MPIE.extract(Malawi)
    SMHI=SMHI.extract(Malawi)
        
    CCCmaCanRCM = CCCmaCanRCM.extract(Malawi)
    CCCmaSMHI = CCCmaSMHI.extract(Malawi)
    CNRM =CNRM.extract(Malawi)
    CNRMSMHI =CNRMSMHI.extract(Malawi)
    CSIRO=CSIRO.extract(Malawi)
    ICHECDMI=ICHECDMI.extract(Malawi)
    ICHECCCLM=ICHECCCLM.extract(Malawi)
    ICHECKNMI=ICHECKNMI.extract(Malawi)
    ICHECMPI=ICHECMPI.extract(Malawi)
    ICHECSMHI=ICHECSMHI.extract(Malawi)
    IPSL=IPSL.extract(Malawi)
    MIROC=MIROC.extract(Malawi)
    MOHCCCLM=MOHCCCLM.extract(Malawi)
    MOHCKNMI=MOHCKNMI.extract(Malawi)
    MOHCSMHI=MOHCSMHI.extract(Malawi)
    MPICCLM=MPICCLM.extract(Malawi)
    MPIREMO=MPIREMO.extract(Malawi)
    MPISMHI=MPISMHI.extract(Malawi)
    NCCDMI=NCCDMI.extract(Malawi)
    NCCSMHI=NCCSMHI.extract(Malawi)
    NOAA=NOAA.extract(Malawi)
    
    CanESM2 = CanESM2.extract(Malawi)
    CNRMG =CNRMG.extract(Malawi)
    MK3 =MK3.extract(Malawi)
    EARTH =EARTH.extract(Malawi)
    EARTH3 =EARTH3.extract(Malawi)
    GFDL=GFDL.extract(Malawi)
    HadGEM2=HadGEM2.extract(Malawi)
    IPSLG=IPSLG.extract(Malawi)
    MIROCG=MIROCG.extract(Malawi)
    MPI=MPI.extract(Malawi)
    NorESM1=NorESM1.extract(Malawi)
    
    CRU = CRU.extract(Malawi)
       
    #time constraignt to make all series the same, for ERAINT this is 1990-2008 and for RCMs and GCMs this is 1961-2005
    iris.FUTURE.cell_datetime_objects = True
    t_constraint_ERAINT = iris.Constraint(time=lambda cell: 1990 <= cell.point.year <= 2008)
   
    CCCma = CCCma.extract(t_constraint_ERAINT)
    CLMcom =CLMcom.extract(t_constraint_ERAINT)
    DMI=DMI.extract(t_constraint_ERAINT)
    KNMI=KNMI.extract(t_constraint_ERAINT)
    MPIE=MPIE.extract(t_constraint_ERAINT)
    SMHI=SMHI.extract(t_constraint_ERAINT)
   
    CRUE = CRU.extract(t_constraint_ERAINT)
        
    t_constraint = iris.Constraint(time=lambda cell: 1961 <= cell.point.year <= 2005)
    
    CCCmaCanRCM = CCCmaCanRCM.extract(t_constraint)
    CCCmaSMHI = CCCmaSMHI.extract(t_constraint)
    CNRM =CNRM.extract(t_constraint)
    CNRMSMHI =CNRMSMHI.extract(t_constraint)
    CSIRO=CSIRO.extract(t_constraint)
    ICHECDMI=ICHECDMI.extract(t_constraint)
    ICHECCCLM=ICHECCCLM.extract(t_constraint)
    ICHECKNMI=ICHECKNMI.extract(t_constraint)
    ICHECMPI=ICHECMPI.extract(t_constraint)
    ICHECSMHI=ICHECSMHI.extract(t_constraint)
    IPSL=IPSL.extract(t_constraint)
    MIROC=MIROC.extract(t_constraint)
    MOHCCCLM=MOHCCCLM.extract(t_constraint)
    MOHCKNMI=MOHCKNMI.extract(t_constraint)
    MOHCSMHI=MOHCSMHI.extract(t_constraint)
    MPICCLM=MPICCLM.extract(t_constraint)
    MPIREMO=MPIREMO.extract(t_constraint)
    MPISMHI=MPISMHI.extract(t_constraint)
    NCCDMI=NCCDMI.extract(t_constraint)
    NCCSMHI=NCCSMHI.extract(t_constraint)
    NOAA=NOAA.extract(t_constraint) 
    
    CanESM2 = CanESM2.extract(t_constraint)
    CNRMG =CNRMG.extract(t_constraint)
    MK3 =MK3.extract(t_constraint)
    EARTH =EARTH.extract(t_constraint)
    EARTH3 =EARTH3.extract(t_constraint)
    GFDL=GFDL.extract(t_constraint)
    HadGEM2=HadGEM2.extract(t_constraint)
    IPSLG=IPSLG.extract(t_constraint)
    MIROCG=MIROCG.extract(t_constraint)
    MPI=MPI.extract(t_constraint)
    NorESM1=NorESM1.extract(t_constraint)
    
    CRU = CRU.extract(t_constraint)
        
    #Convert units to match, CORDEX data is in Kelvin but Observed data in Celsius, we would like to show all data in Celsius
    CCCma.convert_units('Celsius')
    CLMcom.convert_units('Celsius')
    DMI.convert_units('Celsius')
    KNMI.convert_units('Celsius')
    MPIE.convert_units('Celsius')
    SMHI.convert_units('Celsius')
    
    CCCmaCanRCM.convert_units('Celsius')
    CCCmaSMHI.convert_units('Celsius')
    CNRM.convert_units('Celsius')
    CNRMSMHI.convert_units('Celsius')
    CSIRO.convert_units('Celsius')
    ICHECDMI.convert_units('Celsius')
    ICHECCCLM.convert_units('Celsius') 
    ICHECKNMI.convert_units('Celsius')
    ICHECMPI.convert_units('Celsius')
    ICHECSMHI.convert_units('Celsius')
    IPSL.convert_units('Celsius')
    MIROC.convert_units('Celsius')
    MOHCCCLM.convert_units('Celsius')
    MOHCKNMI.convert_units('Celsius')
    MOHCSMHI.convert_units('Celsius')
    MPICCLM.convert_units('Celsius')
    MPIREMO.convert_units('Celsius')
    MPISMHI.convert_units('Celsius')
    NCCDMI.convert_units('Celsius')
    NCCSMHI.convert_units('Celsius')
    NOAA.convert_units('Celsius')   
    
    CanESM2.convert_units('Celsius')
    CNRMG.convert_units('Celsius')
    MK3.convert_units('Celsius')
    EARTH.convert_units('Celsius')
    EARTH3.units = Unit('Celsius') #this fixes EARTH3 which has no units defined
    EARTH3=EARTH3-273 #this converts the data manually from Kelvin to Celsius
    GFDL.convert_units('Celsius')
    HadGEM2.convert_units('Celsius')
    IPSLG.convert_units('Celsius')
    MIROCG.convert_units('Celsius')
    MPI.convert_units('Celsius')
    NorESM1.convert_units('Celsius') 
    
    #limit to a season
    SON = iris.Constraint(season='son')
    DJF = iris.Constraint(season='djf')
    MAM = iris.Constraint(season='mam')
    JJA = iris.Constraint(season='jja')
    
    #add season data to files
    iriscc.add_season(CCCma, 'time')
    iriscc.add_season(CLMcom, 'time')
    iriscc.add_season(DMI, 'time')
    iriscc.add_season(KNMI, 'time')
    iriscc.add_season(MPIE, 'time')
    iriscc.add_season(SMHI, 'time')
        
    iriscc.add_season(CCCmaCanRCM, 'time')
    iriscc.add_season(CCCmaSMHI, 'time')
    iriscc.add_season(CNRM, 'time')
    iriscc.add_season(CNRMSMHI, 'time')
    iriscc.add_season(CSIRO, 'time')
    iriscc.add_season(ICHECDMI, 'time')
    iriscc.add_season(ICHECCCLM, 'time')
    iriscc.add_season(ICHECKNMI, 'time')
    iriscc.add_season(ICHECMPI, 'time')
    iriscc.add_season(ICHECSMHI, 'time')
    iriscc.add_season(IPSL, 'time')
    iriscc.add_season(MIROC, 'time')
    iriscc.add_season(MOHCCCLM, 'time')
    iriscc.add_season(MOHCKNMI, 'time')
    iriscc.add_season(MOHCSMHI, 'time')
    iriscc.add_season(MPICCLM, 'time')
    iriscc.add_season(MPIREMO, 'time')
    iriscc.add_season(MPISMHI, 'time')
    iriscc.add_season(NCCDMI, 'time')
    iriscc.add_season(NCCSMHI, 'time')
    iriscc.add_season(NOAA, 'time')
    
    iriscc.add_season(CanESM2, 'time')
    iriscc.add_season(CNRMG, 'time')
    iriscc.add_season(MK3, 'time')
    iriscc.add_season(EARTH, 'time')
    iriscc.add_season(EARTH3, 'time')
    iriscc.add_season(GFDL, 'time')
    iriscc.add_season(HadGEM2, 'time')
    iriscc.add_season(IPSLG, 'time')
    iriscc.add_season(MIROCG, 'time')
    iriscc.add_season(MPI, 'time')
    iriscc.add_season(NorESM1, 'time')
    
    iriscc.add_season(CRUE, 'time')
    iriscc.add_season(CRU, 'time')
        
    #add year data to files
    iriscc.add_year(CCCma, 'time')
    iriscc.add_year(CLMcom, 'time')
    iriscc.add_year(DMI, 'time')
    iriscc.add_year(KNMI, 'time')
    iriscc.add_year(MPIE, 'time')
    iriscc.add_year(SMHI, 'time')
        
    iriscc.add_year(CCCmaCanRCM, 'time')
    iriscc.add_year(CCCmaSMHI, 'time')
    iriscc.add_year(CNRM, 'time')
    iriscc.add_year(CNRMSMHI, 'time')
    iriscc.add_year(CSIRO, 'time')
    iriscc.add_year(ICHECDMI, 'time')
    iriscc.add_year(ICHECCCLM, 'time')
    iriscc.add_year(ICHECKNMI, 'time')
    iriscc.add_year(ICHECMPI, 'time')
    iriscc.add_year(ICHECSMHI, 'time')
    iriscc.add_year(IPSL, 'time')
    iriscc.add_year(MIROC, 'time')
    iriscc.add_year(MOHCCCLM, 'time')
    iriscc.add_year(MOHCKNMI, 'time')
    iriscc.add_year(MOHCSMHI, 'time')
    iriscc.add_year(MPICCLM, 'time')
    iriscc.add_year(MPIREMO, 'time')
    iriscc.add_year(MPISMHI, 'time')
    iriscc.add_year(NCCDMI, 'time')
    iriscc.add_year(NCCSMHI, 'time')
    iriscc.add_year(NOAA, 'time')
    
    iriscc.add_year(CanESM2, 'time')
    iriscc.add_year(CNRMG, 'time')
    iriscc.add_year(MK3, 'time')
    iriscc.add_year(EARTH, 'time')
    iriscc.add_year(EARTH3, 'time')
    iriscc.add_year(GFDL, 'time')
    iriscc.add_year(HadGEM2, 'time')
    iriscc.add_year(IPSLG, 'time')
    iriscc.add_year(MIROCG, 'time')
    iriscc.add_year(MPI, 'time')
    iriscc.add_year(NorESM1, 'time')
    
    iriscc.add_year(CRUE, 'time')
    iriscc.add_year(CRU, 'time')
        
    #extract only the season we are interested in
    CCCmaSON = CCCma.extract(SON) 
    CLMcomSON = CLMcom.extract(SON) 
    DMISON=DMI.extract(SON) 
    KNMISON=KNMI.extract(SON) 
    MPIESON=MPIE.extract(SON) 
    SMHISON = SMHI.extract(SON) 
    
    CCCmaCanRCMSON = CCCmaCanRCM.extract(SON) 
    CCCmaSMHISON = CCCmaSMHI.extract(SON) 
    CNRMSON = CNRM.extract(SON) 
    CNRMSMHISON = CNRMSMHI.extract(SON) 
    CSIROSON = CSIRO.extract(SON) 
    ICHECDMISON = ICHECDMI.extract(SON) 
    ICHECCCLMSON = ICHECCCLM.extract(SON) 
    ICHECKNMISON = ICHECKNMI.extract(SON) 
    ICHECMPISON = ICHECMPI.extract(SON) 
    ICHECSMHISON = ICHECSMHI.extract(SON) 
    IPSLSON = IPSL.extract(SON) 
    MIROCSON = MIROC.extract(SON) 
    MOHCCCLMSON = MOHCCCLM.extract(SON) 
    MOHCKNMISON = MOHCKNMI.extract(SON) 
    MOHCSMHISON = MOHCSMHI.extract(SON) 
    MPICCLMSON = MPICCLM.extract(SON) 
    MPIREMOSON = MPIREMO.extract(SON) 
    MPISMHISON = MPISMHI.extract(SON) 
    NCCDMISON = NCCDMI.extract(SON) 
    NCCSMHISON = NCCSMHI.extract(SON) 
    NOAASON = NOAA.extract(SON) 
    
    CanESM2SON = CanESM2.extract(SON) 
    CNRMGSON = CNRMG.extract(SON) 
    MK3SON = MK3.extract(SON) 
    EARTHSON = EARTH.extract(SON) 
    EARTH3SON = EARTH3.extract(SON) 
    GFDLSON = GFDL.extract(SON) 
    HadGEM2SON = HadGEM2.extract(SON)  
    IPSLGSON = IPSLG.extract(SON) 
    MIROCGSON = MIROCG.extract(SON) 
    MPISON = MPI.extract(SON) 
    NorESM1SON = NorESM1.extract(SON) 
    
    CRUESON = CRUE.extract(SON)
    CRUSON = CRU.extract(SON)
           
    CCCmaDJF = CCCma.extract(DJF)
    CLMcomDJF=CLMcom.extract(DJF) 
    DMIDJF=DMI.extract(DJF) 
    KNMIDJF=KNMI.extract(DJF) 
    MPIEDJF=MPIE.extract(DJF) 
    SMHIDJF=SMHI.extract(DJF) 
    
    CCCmaCanRCMDJF = CCCmaCanRCM.extract(DJF) 
    CCCmaSMHIDJF = CCCmaSMHI.extract(DJF) 
    CNRMDJF = CNRM.extract(DJF) 
    CNRMSMHIDJF = CNRMSMHI.extract(DJF) 
    CSIRODJF = CSIRO.extract(DJF) 
    ICHECDMIDJF = ICHECDMI.extract(DJF) 
    ICHECCCLMDJF = ICHECCCLM.extract(DJF) 
    ICHECKNMIDJF = ICHECKNMI.extract(DJF) 
    ICHECMPIDJF = ICHECMPI.extract(DJF) 
    ICHECSMHIDJF = ICHECSMHI.extract(DJF) 
    IPSLDJF = IPSL.extract(DJF) 
    MIROCDJF = MIROC.extract(DJF) 
    MOHCCCLMDJF = MOHCCCLM.extract(DJF) 
    MOHCKNMIDJF = MOHCKNMI.extract(DJF) 
    MOHCSMHIDJF = MOHCSMHI.extract(DJF) 
    MPICCLMDJF = MPICCLM.extract(DJF) 
    MPIREMODJF = MPIREMO.extract(DJF) 
    MPISMHIDJF = MPISMHI.extract(DJF) 
    NCCDMIDJF = NCCDMI.extract(DJF) 
    NCCSMHIDJF = NCCSMHI.extract(DJF) 
    NOAADJF = NOAA.extract(DJF) 
    
    CanESM2DJF = CanESM2.extract(DJF) 
    CNRMGDJF = CNRMG.extract(DJF) 
    MK3DJF = MK3.extract(DJF) 
    EARTHDJF = EARTH.extract(DJF) 
    EARTH3DJF = EARTH3.extract(DJF) 
    GFDLDJF = GFDL.extract(DJF) 
    HadGEM2DJF = HadGEM2.extract(DJF)  
    IPSLGDJF = IPSLG.extract(DJF) 
    MIROCGDJF = MIROCG.extract(DJF)  
    MPIDJF = MPI.extract(DJF) 
    NorESM1DJF = NorESM1.extract(DJF)
    
    CRUEDJF = CRUE.extract(DJF)
    CRUDJF = CRU.extract(DJF)
         
    CCCmaMAM = CCCma.extract(MAM) 
    CLMcomMAM=CLMcom.extract(MAM) 
    DMIMAM=DMI.extract(MAM) 
    KNMIMAM=KNMI.extract(MAM) 
    MPIEMAM=MPIE.extract(MAM) 
    SMHIMAM=SMHI.extract(MAM) 
    
    CCCmaCanRCMMAM = CCCmaCanRCM.extract(MAM) 
    CCCmaSMHIMAM = CCCmaSMHI.extract(MAM) 
    CNRMMAM = CNRM.extract(MAM) 
    CNRMSMHIMAM = CNRMSMHI.extract(MAM) 
    CSIROMAM = CSIRO.extract(MAM) 
    ICHECDMIMAM = ICHECDMI.extract(MAM) 
    ICHECCCLMMAM = ICHECCCLM.extract(MAM) 
    ICHECKNMIMAM = ICHECKNMI.extract(MAM) 
    ICHECMPIMAM = ICHECMPI.extract(MAM) 
    ICHECSMHIMAM = ICHECSMHI.extract(MAM) 
    IPSLMAM = IPSL.extract(MAM) 
    MIROCMAM = MIROC.extract(MAM) 
    MOHCCCLMMAM = MOHCCCLM.extract(MAM) 
    MOHCKNMIMAM = MOHCKNMI.extract(MAM) 
    MOHCSMHIMAM = MOHCSMHI.extract(MAM) 
    MPICCLMMAM = MPICCLM.extract(MAM) 
    MPIREMOMAM = MPIREMO.extract(MAM) 
    MPISMHIMAM = MPISMHI.extract(MAM) 
    NCCDMIMAM = NCCDMI.extract(MAM) 
    NCCSMHIMAM = NCCSMHI.extract(MAM) 
    NOAAMAM = NOAA.extract(MAM) 
    
    CanESM2MAM = CanESM2.extract(MAM) 
    CNRMGMAM = CNRMG.extract(MAM) 
    MK3MAM = MK3.extract(MAM) 
    EARTHMAM = EARTH.extract(MAM)
    EARTH3MAM = EARTH3.extract(MAM)
    GFDLMAM = GFDL.extract(MAM) 
    HadGEM2MAM = HadGEM2.extract(MAM) 
    IPSLGMAM = IPSLG.extract(MAM) 
    MIROCGMAM = MIROCG.extract(MAM)  
    MPIMAM = MPI.extract(MAM) 
    NorESM1MAM = NorESM1.extract(MAM)
    
    CRUEMAM = CRUE.extract(MAM)
    CRUMAM = CRU.extract(MAM)
             
    CCCmaJJA = CCCma.extract(JJA)
    CLMcomJJA=CLMcom.extract(JJA) 
    DMIJJA=DMI.extract(JJA) 
    KNMIJJA=KNMI.extract(JJA) 
    MPIEJJA=MPIE.extract(JJA) 
    SMHIJJA=SMHI.extract(JJA) 
    
    CCCmaCanRCMJJA = CCCmaCanRCM.extract(JJA) 
    CCCmaSMHIJJA = CCCmaSMHI.extract(JJA) 
    CNRMJJA = CNRM.extract(JJA) 
    CNRMSMHIJJA = CNRMSMHI.extract(JJA) 
    CSIROJJA = CSIRO.extract(JJA) 
    ICHECDMIJJA = ICHECDMI.extract(JJA) 
    ICHECCCLMJJA = ICHECCCLM.extract(JJA) 
    ICHECKNMIJJA = ICHECKNMI.extract(JJA) 
    ICHECMPIJJA = ICHECMPI.extract(JJA) 
    ICHECSMHIJJA = ICHECSMHI.extract(JJA) 
    IPSLJJA = IPSL.extract(JJA) 
    MIROCJJA = MIROC.extract(JJA) 
    MOHCCCLMJJA = MOHCCCLM.extract(JJA) 
    MOHCKNMIJJA = MOHCKNMI.extract(JJA) 
    MOHCSMHIJJA = MOHCSMHI.extract(JJA) 
    MPICCLMJJA = MPICCLM.extract(JJA) 
    MPIREMOJJA = MPIREMO.extract(JJA) 
    MPISMHIJJA = MPISMHI.extract(JJA) 
    NCCDMIJJA = NCCDMI.extract(JJA) 
    NCCSMHIJJA = NCCSMHI.extract(JJA) 
    NOAAJJA = NOAA.extract(JJA) 
    
    CanESM2JJA = CanESM2.extract(JJA) 
    CNRMGJJA = CNRMG.extract(JJA) 
    MK3JJA = MK3.extract(JJA) 
    EARTHJJA = EARTH.extract(JJA) 
    EARTH3JJA = EARTH3.extract(JJA) 
    GFDLJJA = GFDL.extract(JJA) 
    HadGEM2JJA = HadGEM2.extract(JJA) 
    IPSLGJJA = IPSLG.extract(JJA) 
    MIROCGJJA = MIROCG.extract(JJA)  
    MPIJJA = MPI.extract(JJA) 
    NorESM1JJA = NorESM1.extract(JJA)
    
    CRUEJJA = CRUE.extract(JJA)
    CRUJJA = CRU.extract(JJA)
        
    #We are interested in plotting the data by year, so we need to take a mean of all the data by year
    CCCmaSON = CCCmaSON.aggregated_by('year', iris.analysis.MEAN)
    CLMcomSON = CLMcomSON.aggregated_by('year', iris.analysis.MEAN)    
    DMISON = DMISON.aggregated_by('year', iris.analysis.MEAN)    
    KNMISON = KNMISON.aggregated_by('year', iris.analysis.MEAN)
    MPIESON = MPIESON.aggregated_by('year', iris.analysis.MEAN)    
    SMHISON = SMHISON.aggregated_by('year', iris.analysis.MEAN)    
    
    CCCmaCanRCMSON = CCCmaCanRCMSON.aggregated_by('year', iris.analysis.MEAN)
    CCCmaSMHISON = CCCmaSMHISON.aggregated_by('year', iris.analysis.MEAN)
    CNRMSON = CNRMSON.aggregated_by('year', iris.analysis.MEAN)
    CNRMSMHISON = CNRMSMHISON.aggregated_by('year', iris.analysis.MEAN)
    CSIROSON = CSIROSON.aggregated_by('year', iris.analysis.MEAN)
    ICHECDMISON = ICHECDMISON.aggregated_by('year', iris.analysis.MEAN)
    ICHECCCLMSON = ICHECCCLMSON.aggregated_by('year', iris.analysis.MEAN)
    ICHECKNMISON = ICHECKNMISON.aggregated_by('year', iris.analysis.MEAN)
    ICHECMPISON = ICHECMPISON.aggregated_by('year', iris.analysis.MEAN)
    ICHECSMHISON = ICHECSMHISON.aggregated_by('year', iris.analysis.MEAN)
    IPSLSON = IPSLSON.aggregated_by('year', iris.analysis.MEAN)
    MIROCSON = MIROCSON.aggregated_by('year', iris.analysis.MEAN)
    MOHCCCLMSON = MOHCCCLMSON.aggregated_by('year', iris.analysis.MEAN)
    MOHCKNMISON = MOHCKNMISON.aggregated_by('year', iris.analysis.MEAN)
    MOHCSMHISON = MOHCSMHISON.aggregated_by('year', iris.analysis.MEAN)
    MPICCLMSON = MPICCLMSON.aggregated_by('year', iris.analysis.MEAN)
    MPIREMOSON = MPIREMOSON.aggregated_by('year', iris.analysis.MEAN)
    MPISMHISON = MPISMHISON.aggregated_by('year', iris.analysis.MEAN)
    NCCDMISON = NCCDMISON.aggregated_by('year', iris.analysis.MEAN)
    NCCSMHISON = NCCSMHISON.aggregated_by('year', iris.analysis.MEAN)
    NOAASON = NOAASON.aggregated_by('year', iris.analysis.MEAN)
    
    CanESM2SON = CanESM2SON.aggregated_by('year', iris.analysis.MEAN)
    CNRMGSON = CNRMGSON.aggregated_by('year', iris.analysis.MEAN)
    MK3SON = MK3SON.aggregated_by('year', iris.analysis.MEAN)
    EARTHSON = EARTHSON.aggregated_by('year', iris.analysis.MEAN)
    EARTH3SON = EARTH3SON.aggregated_by('year', iris.analysis.MEAN)
    GFDLSON = GFDLSON.aggregated_by('year', iris.analysis.MEAN)
    HadGEM2SON = HadGEM2SON.aggregated_by('year', iris.analysis.MEAN)
    IPSLGSON = IPSLGSON.aggregated_by('year', iris.analysis.MEAN)
    MIROCGSON = MIROCGSON.aggregated_by('year', iris.analysis.MEAN)
    MPISON = MPISON.aggregated_by('year', iris.analysis.MEAN)
    NorESM1SON = NorESM1SON.aggregated_by('year', iris.analysis.MEAN)   
    
    CRUESON = CRUESON.aggregated_by('year', iris.analysis.MEAN)  
    CRUSON = CRUSON.aggregated_by('year', iris.analysis.MEAN)  
           
    CCCmaDJF = CCCmaDJF.aggregated_by('year', iris.analysis.MEAN)
    CLMcomDJF = CLMcomDJF.aggregated_by('year', iris.analysis.MEAN)    
    DMIDJF = DMIDJF.aggregated_by('year', iris.analysis.MEAN)    
    KNMIDJF = KNMIDJF.aggregated_by('year', iris.analysis.MEAN)
    MPIEDJF = MPIEDJF.aggregated_by('year', iris.analysis.MEAN)    
    SMHIDJF = SMHIDJF.aggregated_by('year', iris.analysis.MEAN)    
    
    CCCmaCanRCMDJF = CCCmaCanRCMDJF.aggregated_by('year', iris.analysis.MEAN)
    CCCmaSMHIDJF = CCCmaSMHIDJF.aggregated_by('year', iris.analysis.MEAN)
    CNRMDJF = CNRMDJF.aggregated_by('year', iris.analysis.MEAN)
    CNRMSMHIDJF = CNRMSMHIDJF.aggregated_by('year', iris.analysis.MEAN)
    CSIRODJF = CSIRODJF.aggregated_by('year', iris.analysis.MEAN)
    ICHECDMIDJF = ICHECDMIDJF.aggregated_by('year', iris.analysis.MEAN)
    ICHECCCLMDJF = ICHECCCLMDJF.aggregated_by('year', iris.analysis.MEAN)
    ICHECKNMIDJF = ICHECKNMIDJF.aggregated_by('year', iris.analysis.MEAN)
    ICHECMPIDJF = ICHECMPIDJF.aggregated_by('year', iris.analysis.MEAN)
    ICHECSMHIDJF = ICHECSMHIDJF.aggregated_by('year', iris.analysis.MEAN)
    IPSLDJF = IPSLDJF.aggregated_by('year', iris.analysis.MEAN)
    MIROCDJF = MIROCDJF.aggregated_by('year', iris.analysis.MEAN)
    MOHCCCLMDJF = MOHCCCLMDJF.aggregated_by('year', iris.analysis.MEAN)
    MOHCKNMIDJF = MOHCKNMIDJF.aggregated_by('year', iris.analysis.MEAN)
    MOHCSMHIDJF = MOHCSMHIDJF.aggregated_by('year', iris.analysis.MEAN)
    MPICCLMDJF = MPICCLMDJF.aggregated_by('year', iris.analysis.MEAN)
    MPIREMODJF = MPIREMODJF.aggregated_by('year', iris.analysis.MEAN)
    MPISMHIDJF = MPISMHIDJF.aggregated_by('year', iris.analysis.MEAN)
    NCCDMIDJF = NCCDMIDJF.aggregated_by('year', iris.analysis.MEAN)
    NCCSMHIDJF = NCCSMHIDJF.aggregated_by('year', iris.analysis.MEAN)
    NOAADJF = NOAADJF.aggregated_by('year', iris.analysis.MEAN)
    
    CanESM2DJF = CanESM2DJF.aggregated_by('year', iris.analysis.MEAN)
    CNRMGDJF = CNRMGDJF.aggregated_by('year', iris.analysis.MEAN)
    MK3DJF = MK3DJF.aggregated_by('year', iris.analysis.MEAN)
    EARTHDJF = EARTHDJF.aggregated_by('year', iris.analysis.MEAN)
    EARTH3DJF = EARTH3DJF.aggregated_by('year', iris.analysis.MEAN)
    GFDLDJF = GFDLDJF.aggregated_by('year', iris.analysis.MEAN)
    HadGEM2DJF = HadGEM2DJF.aggregated_by('year', iris.analysis.MEAN)
    IPSLGDJF = IPSLGDJF.aggregated_by('year', iris.analysis.MEAN)
    MIROCGDJF = MIROCGDJF.aggregated_by('year', iris.analysis.MEAN)
    MPIDJF = MPIDJF.aggregated_by('year', iris.analysis.MEAN)
    NorESM1DJF = NorESM1DJF.aggregated_by('year', iris.analysis.MEAN)      
    
    CRUEDJF = CRUEDJF.aggregated_by('year', iris.analysis.MEAN)  
    CRUDJF = CRUDJF.aggregated_by('year', iris.analysis.MEAN)  
          
    CCCmaMAM = CCCmaMAM.aggregated_by('year', iris.analysis.MEAN)
    CLMcomMAM = CLMcomMAM.aggregated_by('year', iris.analysis.MEAN)    
    DMIMAM = DMIMAM.aggregated_by('year', iris.analysis.MEAN)    
    KNMIMAM = KNMIMAM.aggregated_by('year', iris.analysis.MEAN)
    MPIEMAM = MPIEMAM.aggregated_by('year', iris.analysis.MEAN)    
    SMHIMAM = SMHIMAM.aggregated_by('year', iris.analysis.MEAN)    
    
    CCCmaCanRCMMAM = CCCmaCanRCMMAM.aggregated_by('year', iris.analysis.MEAN)
    CCCmaSMHIMAM = CCCmaSMHIMAM.aggregated_by('year', iris.analysis.MEAN)
    CNRMMAM = CNRMMAM.aggregated_by('year', iris.analysis.MEAN)
    CNRMSMHIMAM = CNRMSMHIMAM.aggregated_by('year', iris.analysis.MEAN)
    CSIROMAM = CSIROMAM.aggregated_by('year', iris.analysis.MEAN)
    ICHECDMIMAM = ICHECDMIMAM.aggregated_by('year', iris.analysis.MEAN)
    ICHECCCLMMAM = ICHECCCLMMAM.aggregated_by('year', iris.analysis.MEAN)
    ICHECKNMIMAM = ICHECKNMIMAM.aggregated_by('year', iris.analysis.MEAN)
    ICHECMPIMAM = ICHECMPIMAM.aggregated_by('year', iris.analysis.MEAN)
    ICHECSMHIMAM = ICHECSMHIMAM.aggregated_by('year', iris.analysis.MEAN)
    IPSLMAM = IPSLMAM.aggregated_by('year', iris.analysis.MEAN)
    MIROCMAM = MIROCMAM.aggregated_by('year', iris.analysis.MEAN)
    MOHCCCLMMAM = MOHCCCLMMAM.aggregated_by('year', iris.analysis.MEAN)
    MOHCKNMIMAM = MOHCKNMIMAM.aggregated_by('year', iris.analysis.MEAN)
    MOHCSMHIMAM = MOHCSMHIMAM.aggregated_by('year', iris.analysis.MEAN)
    MPICCLMMAM = MPICCLMMAM.aggregated_by('year', iris.analysis.MEAN)
    MPIREMOMAM = MPIREMOMAM.aggregated_by('year', iris.analysis.MEAN)
    MPISMHIMAM = MPISMHIMAM.aggregated_by('year', iris.analysis.MEAN)
    NCCDMIMAM = NCCDMIMAM.aggregated_by('year', iris.analysis.MEAN)
    NCCSMHIMAM = NCCSMHIMAM.aggregated_by('year', iris.analysis.MEAN)
    NOAAMAM = NOAAMAM.aggregated_by('year', iris.analysis.MEAN)
    
    CanESM2MAM = CanESM2MAM.aggregated_by('year', iris.analysis.MEAN)
    CNRMGMAM = CNRMGMAM.aggregated_by('year', iris.analysis.MEAN)
    MK3MAM = MK3MAM.aggregated_by('year', iris.analysis.MEAN)
    EARTHMAM = EARTHMAM.aggregated_by('year', iris.analysis.MEAN)
    EARTH3MAM = EARTH3MAM.aggregated_by('year', iris.analysis.MEAN)
    GFDLMAM = GFDLMAM.aggregated_by('year', iris.analysis.MEAN)
    HadGEM2MAM = HadGEM2MAM.aggregated_by('year', iris.analysis.MEAN)
    IPSLGMAM = IPSLGMAM.aggregated_by('year', iris.analysis.MEAN)
    MIROCGMAM = MIROCGMAM.aggregated_by('year', iris.analysis.MEAN)
    MPIMAM = MPIMAM.aggregated_by('year', iris.analysis.MEAN)
    NorESM1MAM = NorESM1MAM.aggregated_by('year', iris.analysis.MEAN)   
    
    CRUEMAM = CRUEMAM.aggregated_by('year', iris.analysis.MEAN)  
    CRUMAM = CRUMAM.aggregated_by('year', iris.analysis.MEAN)  
            
    CCCmaJJA = CCCmaJJA.aggregated_by('year', iris.analysis.MEAN)
    CLMcomJJA = CLMcomJJA.aggregated_by('year', iris.analysis.MEAN)    
    DMIJJA = DMIJJA.aggregated_by('year', iris.analysis.MEAN)    
    KNMIJJA = KNMIJJA.aggregated_by('year', iris.analysis.MEAN)
    MPIEJJA = MPIEJJA.aggregated_by('year', iris.analysis.MEAN)    
    SMHIJJA = SMHIJJA.aggregated_by('year', iris.analysis.MEAN)    
    
    CCCmaCanRCMJJA = CCCmaCanRCMJJA.aggregated_by('year', iris.analysis.MEAN)
    CCCmaSMHIJJA = CCCmaSMHIJJA.aggregated_by('year', iris.analysis.MEAN)
    CNRMJJA = CNRMJJA.aggregated_by('year', iris.analysis.MEAN)
    CNRMSMHIJJA = CNRMSMHIJJA.aggregated_by('year', iris.analysis.MEAN)
    CSIROJJA = CSIROJJA.aggregated_by('year', iris.analysis.MEAN)
    ICHECDMIJJA = ICHECDMIJJA.aggregated_by('year', iris.analysis.MEAN)
    ICHECCCLMJJA = ICHECCCLMJJA.aggregated_by('year', iris.analysis.MEAN)
    ICHECKNMIJJA = ICHECKNMIJJA.aggregated_by('year', iris.analysis.MEAN)
    ICHECMPIJJA = ICHECMPIJJA.aggregated_by('year', iris.analysis.MEAN)
    ICHECSMHIJJA = ICHECSMHIJJA.aggregated_by('year', iris.analysis.MEAN)
    IPSLJJA = IPSLJJA.aggregated_by('year', iris.analysis.MEAN)
    MIROCJJA = MIROCJJA.aggregated_by('year', iris.analysis.MEAN)
    MOHCCCLMJJA = MOHCCCLMJJA.aggregated_by('year', iris.analysis.MEAN)
    MOHCKNMIJJA = MOHCKNMIJJA.aggregated_by('year', iris.analysis.MEAN)
    MOHCSMHIJJA = MOHCSMHIJJA.aggregated_by('year', iris.analysis.MEAN)
    MPICCLMJJA = MPICCLMJJA.aggregated_by('year', iris.analysis.MEAN)
    MPIREMOJJA = MPIREMOJJA.aggregated_by('year', iris.analysis.MEAN)
    MPISMHIJJA = MPISMHIJJA.aggregated_by('year', iris.analysis.MEAN)
    NCCDMIJJA = NCCDMIJJA.aggregated_by('year', iris.analysis.MEAN)
    NCCSMHIJJA = NCCSMHIJJA.aggregated_by('year', iris.analysis.MEAN)
    NOAAJJA = NOAAJJA.aggregated_by('year', iris.analysis.MEAN)
    
    CanESM2JJA = CanESM2JJA.aggregated_by('year', iris.analysis.MEAN)
    CNRMGJJA = CNRMGJJA.aggregated_by('year', iris.analysis.MEAN)
    MK3JJA = MK3JJA.aggregated_by('year', iris.analysis.MEAN)
    EARTHJJA = EARTHJJA.aggregated_by('year', iris.analysis.MEAN)
    EARTH3JJA = EARTH3JJA.aggregated_by('year', iris.analysis.MEAN)
    GFDLJJA = GFDLJJA.aggregated_by('year', iris.analysis.MEAN)
    HadGEM2JJA = HadGEM2JJA.aggregated_by('year', iris.analysis.MEAN)
    IPSLGJJA = IPSLGJJA.aggregated_by('year', iris.analysis.MEAN)
    MIROCGJJA = MIROCGJJA.aggregated_by('year', iris.analysis.MEAN)
    MPIJJA = MPIJJA.aggregated_by('year', iris.analysis.MEAN)
    NorESM1JJA = NorESM1JJA.aggregated_by('year', iris.analysis.MEAN)  
    
    CRUEJJA = CRUEJJA.aggregated_by('year', iris.analysis.MEAN)  
    CRUJJA = CRUJJA.aggregated_by('year', iris.analysis.MEAN)  
        
    CCCmaYR = CCCma.aggregated_by('year', iris.analysis.MEAN)
    CLMcomYR = CLMcom.aggregated_by('year', iris.analysis.MEAN)    
    DMIYR = DMI.aggregated_by('year', iris.analysis.MEAN)    
    KNMIYR = KNMI.aggregated_by('year', iris.analysis.MEAN)
    MPIEYR = MPIE.aggregated_by('year', iris.analysis.MEAN)    
    SMHIYR = SMHI.aggregated_by('year', iris.analysis.MEAN)    
    
    CCCmaCanRCMYR = CCCmaCanRCM.aggregated_by('year', iris.analysis.MEAN)
    CCCmaSMHIYR = CCCmaSMHI.aggregated_by('year', iris.analysis.MEAN)
    CNRMYR = CNRM.aggregated_by('year', iris.analysis.MEAN)
    CNRMSMHIYR = CNRMSMHI.aggregated_by('year', iris.analysis.MEAN)
    CSIROYR = CSIRO.aggregated_by('year', iris.analysis.MEAN)
    ICHECDMIYR = ICHECDMI.aggregated_by('year', iris.analysis.MEAN)
    ICHECCCLMYR = ICHECCCLM.aggregated_by('year', iris.analysis.MEAN)
    ICHECKNMIYR = ICHECKNMI.aggregated_by('year', iris.analysis.MEAN)
    ICHECMPIYR = ICHECMPI.aggregated_by('year', iris.analysis.MEAN)
    ICHECSMHIYR = ICHECSMHI.aggregated_by('year', iris.analysis.MEAN)
    IPSLYR = IPSL.aggregated_by('year', iris.analysis.MEAN)
    MIROCYR = MIROC.aggregated_by('year', iris.analysis.MEAN)
    MOHCCCLMYR = MOHCCCLM.aggregated_by('year', iris.analysis.MEAN)
    MOHCKNMIYR = MOHCKNMI.aggregated_by('year', iris.analysis.MEAN)
    MOHCSMHIYR = MOHCSMHI.aggregated_by('year', iris.analysis.MEAN)
    MPICCLMYR = MPICCLM.aggregated_by('year', iris.analysis.MEAN)
    MPIREMOYR = MPIREMO.aggregated_by('year', iris.analysis.MEAN)
    MPISMHIYR = MPISMHI.aggregated_by('year', iris.analysis.MEAN)
    NCCDMIYR = NCCDMI.aggregated_by('year', iris.analysis.MEAN)
    NCCSMHIYR = NCCSMHI.aggregated_by('year', iris.analysis.MEAN)
    NOAAYR = NOAA.aggregated_by('year', iris.analysis.MEAN)
    
    CanESM2YR=CanESM2.aggregated_by('year', iris.analysis.MEAN)
    CNRMGYR = CNRMG.aggregated_by('year', iris.analysis.MEAN)
    MK3YR = MK3.aggregated_by('year', iris.analysis.MEAN)
    EARTHYR = EARTH.aggregated_by('year', iris.analysis.MEAN)
    EARTH3YR = EARTH3.aggregated_by('year', iris.analysis.MEAN)
    GFDLYR = GFDL.aggregated_by('year', iris.analysis.MEAN)
    HadGEM2YR = HadGEM2.aggregated_by('year', iris.analysis.MEAN)
    IPSLGYR = IPSLG.aggregated_by('year', iris.analysis.MEAN)
    MIROCGYR = MIROCG.aggregated_by('year', iris.analysis.MEAN)
    MPIYR = MPI.aggregated_by('year', iris.analysis.MEAN)
    NorESM1YR = NorESM1.aggregated_by('year', iris.analysis.MEAN)  
    
    CRUEYR = CRUE.aggregated_by('year', iris.analysis.MEAN) 
    CRUYR = CRU.aggregated_by('year', iris.analysis.MEAN)  
       
    #Returns an array of area weights, with the same dimensions as the cube
    CCCmaSON_grid_areas = iris.analysis.cartography.area_weights(CCCmaSON)
    CLMcomSON_grid_areas = iris.analysis.cartography.area_weights(CLMcomSON)
    DMISON_grid_areas = iris.analysis.cartography.area_weights(DMISON)
    KNMISON_grid_areas = iris.analysis.cartography.area_weights(KNMISON)
    MPIESON_grid_areas = iris.analysis.cartography.area_weights(MPIESON)
    SMHISON_grid_areas = iris.analysis.cartography.area_weights(SMHISON)
    
    CCCmaCanRCMSON_grid_areas = iris.analysis.cartography.area_weights(CCCmaCanRCMSON)
    CCCmaSMHISON_grid_areas = iris.analysis.cartography.area_weights(CCCmaSMHISON)
    CNRMSON_grid_areas = iris.analysis.cartography.area_weights(CNRMSON)
    CNRMSMHISON_grid_areas = iris.analysis.cartography.area_weights(CNRMSMHISON)
    CSIROSON_grid_areas = iris.analysis.cartography.area_weights(CSIROSON)
    ICHECDMISON_grid_areas = iris.analysis.cartography.area_weights(ICHECDMISON)
    ICHECCCLMSON_grid_areas = iris.analysis.cartography.area_weights(ICHECCCLMSON)
    ICHECKNMISON_grid_areas = iris.analysis.cartography.area_weights(ICHECKNMISON)
    ICHECMPISON_grid_areas = iris.analysis.cartography.area_weights(ICHECMPISON)
    ICHECSMHISON_grid_areas = iris.analysis.cartography.area_weights(ICHECSMHISON)
    IPSLSON_grid_areas = iris.analysis.cartography.area_weights(IPSLSON)
    MIROCSON_grid_areas = iris.analysis.cartography.area_weights(MIROCSON)
    MOHCCCLMSON_grid_areas = iris.analysis.cartography.area_weights(MOHCCCLMSON)
    MOHCKNMISON_grid_areas = iris.analysis.cartography.area_weights(MOHCKNMISON)
    MOHCSMHISON_grid_areas = iris.analysis.cartography.area_weights(MOHCSMHISON)
    MPICCLMSON_grid_areas = iris.analysis.cartography.area_weights(MPICCLMSON)
    MPIREMOSON_grid_areas = iris.analysis.cartography.area_weights(MPIREMOSON)
    MPISMHISON_grid_areas = iris.analysis.cartography.area_weights(MPISMHISON)
    NCCDMISON_grid_areas = iris.analysis.cartography.area_weights(NCCDMISON)
    NCCSMHISON_grid_areas = iris.analysis.cartography.area_weights(NCCSMHISON)
    NOAASON_grid_areas = iris.analysis.cartography.area_weights(NOAASON)
    
    CanESM2SON_grid_areas = iris.analysis.cartography.area_weights(CanESM2SON)
    CNRMGSON_grid_areas = iris.analysis.cartography.area_weights(CNRMGSON)
    MK3SON_grid_areas = iris.analysis.cartography.area_weights(MK3SON)
    EARTHSON_grid_areas = iris.analysis.cartography.area_weights(EARTHSON)
    EARTH3SON_grid_areas = iris.analysis.cartography.area_weights(EARTH3SON)
    GFDLSON_grid_areas = iris.analysis.cartography.area_weights(GFDLSON)
    HadGEM2SON_grid_areas = iris.analysis.cartography.area_weights(HadGEM2SON)
    IPSLGSON_grid_areas = iris.analysis.cartography.area_weights(IPSLGSON)
    MIROCGSON_grid_areas = iris.analysis.cartography.area_weights(MIROCGSON)
    MPISON_grid_areas = iris.analysis.cartography.area_weights(MPISON)
    NorESM1SON_grid_areas = iris.analysis.cartography.area_weights(NorESM1SON)
    
    CRUESON_grid_areas = iris.analysis.cartography.area_weights(CRUESON)
    CRUSON_grid_areas = iris.analysis.cartography.area_weights(CRUSON)
        
    CCCmaDJF_grid_areas = iris.analysis.cartography.area_weights(CCCmaDJF)
    CLMcomDJF_grid_areas = iris.analysis.cartography.area_weights(CLMcomDJF)
    DMIDJF_grid_areas = iris.analysis.cartography.area_weights(DMIDJF)
    KNMIDJF_grid_areas = iris.analysis.cartography.area_weights(KNMIDJF)
    MPIEDJF_grid_areas = iris.analysis.cartography.area_weights(MPIEDJF)
    SMHIDJF_grid_areas = iris.analysis.cartography.area_weights(SMHIDJF)
    
    CCCmaCanRCMDJF_grid_areas = iris.analysis.cartography.area_weights(CCCmaCanRCMDJF)
    CCCmaSMHIDJF_grid_areas = iris.analysis.cartography.area_weights(CCCmaSMHIDJF)
    CNRMDJF_grid_areas = iris.analysis.cartography.area_weights(CNRMDJF)
    CNRMSMHIDJF_grid_areas = iris.analysis.cartography.area_weights(CNRMSMHIDJF)
    CSIRODJF_grid_areas = iris.analysis.cartography.area_weights(CSIRODJF)
    ICHECDMIDJF_grid_areas = iris.analysis.cartography.area_weights(ICHECDMIDJF)
    ICHECCCLMDJF_grid_areas = iris.analysis.cartography.area_weights(ICHECCCLMDJF)
    ICHECKNMIDJF_grid_areas = iris.analysis.cartography.area_weights(ICHECKNMIDJF)
    ICHECMPIDJF_grid_areas = iris.analysis.cartography.area_weights(ICHECMPIDJF)
    ICHECSMHIDJF_grid_areas = iris.analysis.cartography.area_weights(ICHECSMHIDJF)
    IPSLDJF_grid_areas = iris.analysis.cartography.area_weights(IPSLDJF)
    MIROCDJF_grid_areas = iris.analysis.cartography.area_weights(MIROCDJF)
    MOHCCCLMDJF_grid_areas = iris.analysis.cartography.area_weights(MOHCCCLMDJF)
    MOHCKNMIDJF_grid_areas = iris.analysis.cartography.area_weights(MOHCKNMIDJF)
    MOHCSMHIDJF_grid_areas = iris.analysis.cartography.area_weights(MOHCSMHIDJF)
    MPICCLMDJF_grid_areas = iris.analysis.cartography.area_weights(MPICCLMDJF)
    MPIREMODJF_grid_areas = iris.analysis.cartography.area_weights(MPIREMODJF)
    MPISMHIDJF_grid_areas = iris.analysis.cartography.area_weights(MPISMHIDJF)
    NCCDMIDJF_grid_areas = iris.analysis.cartography.area_weights(NCCDMIDJF)
    NCCSMHIDJF_grid_areas = iris.analysis.cartography.area_weights(NCCSMHIDJF)
    NOAADJF_grid_areas = iris.analysis.cartography.area_weights(NOAADJF)
    
    CanESM2DJF_grid_areas = iris.analysis.cartography.area_weights(CanESM2DJF)
    CNRMGDJF_grid_areas = iris.analysis.cartography.area_weights(CNRMGDJF)
    MK3DJF_grid_areas = iris.analysis.cartography.area_weights(MK3DJF)
    EARTHDJF_grid_areas = iris.analysis.cartography.area_weights(EARTHDJF)
    EARTH3DJF_grid_areas = iris.analysis.cartography.area_weights(EARTH3DJF)
    GFDLDJF_grid_areas = iris.analysis.cartography.area_weights(GFDLDJF)
    HadGEM2DJF_grid_areas = iris.analysis.cartography.area_weights(HadGEM2DJF)
    IPSLGDJF_grid_areas = iris.analysis.cartography.area_weights(IPSLGDJF)
    MIROCGDJF_grid_areas = iris.analysis.cartography.area_weights(MIROCGDJF)
    MPIDJF_grid_areas = iris.analysis.cartography.area_weights(MPIDJF)
    NorESM1DJF_grid_areas = iris.analysis.cartography.area_weights(NorESM1DJF)
    
    CRUEDJF_grid_areas = iris.analysis.cartography.area_weights(CRUEDJF)
    CRUDJF_grid_areas = iris.analysis.cartography.area_weights(CRUDJF)
        
    CCCmaMAM_grid_areas = iris.analysis.cartography.area_weights(CCCmaMAM)
    CLMcomMAM_grid_areas = iris.analysis.cartography.area_weights(CLMcomMAM)
    DMIMAM_grid_areas = iris.analysis.cartography.area_weights(DMIMAM)
    KNMIMAM_grid_areas = iris.analysis.cartography.area_weights(KNMIMAM)
    MPIEMAM_grid_areas = iris.analysis.cartography.area_weights(MPIEMAM)
    SMHIMAM_grid_areas = iris.analysis.cartography.area_weights(SMHIMAM)
    
    CCCmaCanRCMMAM_grid_areas = iris.analysis.cartography.area_weights(CCCmaCanRCMMAM)
    CCCmaSMHIMAM_grid_areas = iris.analysis.cartography.area_weights(CCCmaSMHIMAM)
    CNRMMAM_grid_areas = iris.analysis.cartography.area_weights(CNRMMAM)
    CNRMSMHIMAM_grid_areas = iris.analysis.cartography.area_weights(CNRMSMHIMAM)
    CSIROMAM_grid_areas = iris.analysis.cartography.area_weights(CSIROMAM)
    ICHECDMIMAM_grid_areas = iris.analysis.cartography.area_weights(ICHECDMIMAM)
    ICHECCCLMMAM_grid_areas = iris.analysis.cartography.area_weights(ICHECCCLMMAM)
    ICHECKNMIMAM_grid_areas = iris.analysis.cartography.area_weights(ICHECKNMIMAM)
    ICHECMPIMAM_grid_areas = iris.analysis.cartography.area_weights(ICHECMPIMAM)
    ICHECSMHIMAM_grid_areas = iris.analysis.cartography.area_weights(ICHECSMHIMAM)
    IPSLMAM_grid_areas = iris.analysis.cartography.area_weights(IPSLMAM)
    MIROCMAM_grid_areas = iris.analysis.cartography.area_weights(MIROCMAM)
    MOHCCCLMMAM_grid_areas = iris.analysis.cartography.area_weights(MOHCCCLMMAM)
    MOHCKNMIMAM_grid_areas = iris.analysis.cartography.area_weights(MOHCKNMIMAM)
    MOHCSMHIMAM_grid_areas = iris.analysis.cartography.area_weights(MOHCSMHIMAM)
    MPICCLMMAM_grid_areas = iris.analysis.cartography.area_weights(MPICCLMMAM)
    MPIREMOMAM_grid_areas = iris.analysis.cartography.area_weights(MPIREMOMAM)
    MPISMHIMAM_grid_areas = iris.analysis.cartography.area_weights(MPISMHIMAM)
    NCCDMIMAM_grid_areas = iris.analysis.cartography.area_weights(NCCDMIMAM)
    NCCSMHIMAM_grid_areas = iris.analysis.cartography.area_weights(NCCSMHIMAM)
    NOAAMAM_grid_areas = iris.analysis.cartography.area_weights(NOAAMAM)
    
    CanESM2MAM_grid_areas = iris.analysis.cartography.area_weights(CanESM2MAM)
    CNRMGMAM_grid_areas = iris.analysis.cartography.area_weights(CNRMGMAM)
    MK3MAM_grid_areas = iris.analysis.cartography.area_weights(MK3MAM)
    EARTHMAM_grid_areas = iris.analysis.cartography.area_weights(EARTHMAM)
    EARTH3MAM_grid_areas = iris.analysis.cartography.area_weights(EARTH3MAM)
    GFDLMAM_grid_areas = iris.analysis.cartography.area_weights(GFDLMAM)
    HadGEM2MAM_grid_areas = iris.analysis.cartography.area_weights(HadGEM2MAM)
    IPSLGMAM_grid_areas = iris.analysis.cartography.area_weights(IPSLGMAM)
    MIROCGMAM_grid_areas = iris.analysis.cartography.area_weights(MIROCGMAM)
    MPIMAM_grid_areas = iris.analysis.cartography.area_weights(MPIMAM)
    NorESM1MAM_grid_areas = iris.analysis.cartography.area_weights(NorESM1MAM)
    
    CRUEMAM_grid_areas = iris.analysis.cartography.area_weights(CRUEMAM)
    CRUMAM_grid_areas = iris.analysis.cartography.area_weights(CRUMAM)
        
    CCCmaJJA_grid_areas = iris.analysis.cartography.area_weights(CCCmaJJA)
    CLMcomJJA_grid_areas = iris.analysis.cartography.area_weights(CLMcomJJA)
    DMIJJA_grid_areas = iris.analysis.cartography.area_weights(DMIJJA)
    KNMIJJA_grid_areas = iris.analysis.cartography.area_weights(KNMIJJA)
    MPIEJJA_grid_areas = iris.analysis.cartography.area_weights(MPIEJJA)
    SMHIJJA_grid_areas = iris.analysis.cartography.area_weights(SMHIJJA)
    
    CCCmaCanRCMJJA_grid_areas = iris.analysis.cartography.area_weights(CCCmaCanRCMJJA)
    CCCmaSMHIJJA_grid_areas = iris.analysis.cartography.area_weights(CCCmaSMHIJJA)
    CNRMJJA_grid_areas = iris.analysis.cartography.area_weights(CNRMJJA)
    CNRMSMHIJJA_grid_areas = iris.analysis.cartography.area_weights(CNRMSMHIJJA)
    CSIROJJA_grid_areas = iris.analysis.cartography.area_weights(CSIROJJA)
    ICHECDMIJJA_grid_areas = iris.analysis.cartography.area_weights(ICHECDMIJJA)
    ICHECCCLMJJA_grid_areas = iris.analysis.cartography.area_weights(ICHECCCLMJJA)
    ICHECKNMIJJA_grid_areas = iris.analysis.cartography.area_weights(ICHECKNMIJJA)
    ICHECMPIJJA_grid_areas = iris.analysis.cartography.area_weights(ICHECMPIJJA)
    ICHECSMHIJJA_grid_areas = iris.analysis.cartography.area_weights(ICHECSMHIJJA)
    IPSLJJA_grid_areas = iris.analysis.cartography.area_weights(IPSLJJA)
    MIROCJJA_grid_areas = iris.analysis.cartography.area_weights(MIROCJJA)
    MOHCCCLMJJA_grid_areas = iris.analysis.cartography.area_weights(MOHCCCLMJJA)
    MOHCKNMIJJA_grid_areas = iris.analysis.cartography.area_weights(MOHCKNMIJJA)
    MOHCSMHIJJA_grid_areas = iris.analysis.cartography.area_weights(MOHCSMHIJJA)
    MPICCLMJJA_grid_areas = iris.analysis.cartography.area_weights(MPICCLMJJA)
    MPIREMOJJA_grid_areas = iris.analysis.cartography.area_weights(MPIREMOJJA)
    MPISMHIJJA_grid_areas = iris.analysis.cartography.area_weights(MPISMHIJJA)
    NCCDMIJJA_grid_areas = iris.analysis.cartography.area_weights(NCCDMIJJA)
    NCCSMHIJJA_grid_areas = iris.analysis.cartography.area_weights(NCCSMHIJJA)
    NOAAJJA_grid_areas = iris.analysis.cartography.area_weights(NOAAJJA)
    
    CanESM2JJA_grid_areas = iris.analysis.cartography.area_weights(CanESM2JJA)
    CNRMGJJA_grid_areas = iris.analysis.cartography.area_weights(CNRMGJJA)
    MK3JJA_grid_areas = iris.analysis.cartography.area_weights(MK3JJA)
    EARTHJJA_grid_areas = iris.analysis.cartography.area_weights(EARTHJJA)
    EARTH3JJA_grid_areas = iris.analysis.cartography.area_weights(EARTH3JJA)
    GFDLJJA_grid_areas = iris.analysis.cartography.area_weights(GFDLJJA)
    HadGEM2JJA_grid_areas = iris.analysis.cartography.area_weights(HadGEM2JJA)
    IPSLGJJA_grid_areas = iris.analysis.cartography.area_weights(IPSLGJJA)
    MIROCGJJA_grid_areas = iris.analysis.cartography.area_weights(MIROCGJJA)
    MPIJJA_grid_areas = iris.analysis.cartography.area_weights(MPIJJA)
    NorESM1JJA_grid_areas = iris.analysis.cartography.area_weights(NorESM1JJA)
    
    CRUEJJA_grid_areas = iris.analysis.cartography.area_weights(CRUEJJA)
    CRUJJA_grid_areas = iris.analysis.cartography.area_weights(CRUJJA)
        
    CCCmaYR_grid_areas = iris.analysis.cartography.area_weights(CCCmaYR)
    CLMcomYR_grid_areas = iris.analysis.cartography.area_weights(CLMcomYR)
    DMIYR_grid_areas = iris.analysis.cartography.area_weights(DMIYR)
    KNMIYR_grid_areas = iris.analysis.cartography.area_weights(KNMIYR)
    MPIEYR_grid_areas = iris.analysis.cartography.area_weights(MPIEYR)
    SMHIYR_grid_areas = iris.analysis.cartography.area_weights(SMHIYR)
    
    CCCmaCanRCMYR_grid_areas = iris.analysis.cartography.area_weights(CCCmaCanRCMYR)
    CCCmaSMHIYR_grid_areas = iris.analysis.cartography.area_weights(CCCmaSMHIYR)
    CNRMYR_grid_areas = iris.analysis.cartography.area_weights(CNRMYR)
    CNRMSMHIYR_grid_areas = iris.analysis.cartography.area_weights(CNRMSMHIYR)
    CSIROYR_grid_areas = iris.analysis.cartography.area_weights(CSIROYR)
    ICHECDMIYR_grid_areas = iris.analysis.cartography.area_weights(ICHECDMIYR)
    ICHECCCLMYR_grid_areas = iris.analysis.cartography.area_weights(ICHECCCLMYR)
    ICHECKNMIYR_grid_areas = iris.analysis.cartography.area_weights(ICHECKNMIYR)
    ICHECMPIYR_grid_areas = iris.analysis.cartography.area_weights(ICHECMPIYR)
    ICHECSMHIYR_grid_areas = iris.analysis.cartography.area_weights(ICHECSMHIYR)
    IPSLYR_grid_areas = iris.analysis.cartography.area_weights(IPSLYR)
    MIROCYR_grid_areas = iris.analysis.cartography.area_weights(MIROCYR)
    MOHCCCLMYR_grid_areas = iris.analysis.cartography.area_weights(MOHCCCLMYR)
    MOHCKNMIYR_grid_areas = iris.analysis.cartography.area_weights(MOHCKNMIYR)
    MOHCSMHIYR_grid_areas = iris.analysis.cartography.area_weights(MOHCSMHIYR)
    MPICCLMYR_grid_areas = iris.analysis.cartography.area_weights(MPICCLMYR)
    MPIREMOYR_grid_areas = iris.analysis.cartography.area_weights(MPIREMOYR)
    MPISMHIYR_grid_areas = iris.analysis.cartography.area_weights(MPISMHIYR)
    NCCDMIYR_grid_areas = iris.analysis.cartography.area_weights(NCCDMIYR)
    NCCSMHIYR_grid_areas = iris.analysis.cartography.area_weights(NCCSMHIYR)
    NOAAYR_grid_areas = iris.analysis.cartography.area_weights(NOAAYR)
    
    CanESM2YR_grid_areas = iris.analysis.cartography.area_weights(CanESM2YR)
    CNRMGYR_grid_areas = iris.analysis.cartography.area_weights(CNRMGYR)
    MK3YR_grid_areas = iris.analysis.cartography.area_weights(MK3YR)
    EARTHYR_grid_areas = iris.analysis.cartography.area_weights(EARTHYR)
    EARTH3YR_grid_areas = iris.analysis.cartography.area_weights(EARTH3YR)
    GFDLYR_grid_areas = iris.analysis.cartography.area_weights(GFDLYR)
    HadGEM2YR_grid_areas = iris.analysis.cartography.area_weights(HadGEM2YR)
    IPSLGYR_grid_areas = iris.analysis.cartography.area_weights(IPSLGYR)
    MIROCGYR_grid_areas = iris.analysis.cartography.area_weights(MIROCGYR)
    MPIYR_grid_areas = iris.analysis.cartography.area_weights(MPIYR)
    NorESM1YR_grid_areas = iris.analysis.cartography.area_weights(NorESM1YR)
    
    CRUEYR_grid_areas = iris.analysis.cartography.area_weights(CRUEYR)
    CRUYR_grid_areas = iris.analysis.cartography.area_weights(CRUYR)
    
    #We want to plot the mean for the whole region so we need a mean of all the lats and lons
    CCCmaSON_mean = CCCmaSON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CCCmaSON_grid_areas)                                        
    CLMcomSON_mean = CLMcomSON.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=CLMcomSON_grid_areas)
    DMISON_mean = DMISON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=DMISON_grid_areas)     
    KNMISON_mean = KNMISON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=KNMISON_grid_areas)
    MPIESON_mean = MPIESON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIESON_grid_areas)
    SMHISON_mean = SMHISON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=SMHISON_grid_areas)
    
    CCCmaCanRCMSON_mean = CCCmaCanRCMSON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CCCmaCanRCMSON_grid_areas) 
    CCCmaSMHISON_mean = CCCmaSMHISON.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=CCCmaSMHISON_grid_areas)
    CNRMSON_mean = CNRMSON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMSON_grid_areas)                                               
    CNRMSMHISON_mean = CNRMSMHISON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMSMHISON_grid_areas)  
    CSIROSON_mean = CSIROSON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CSIROSON_grid_areas)
    ICHECDMISON_mean = ICHECDMISON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECDMISON_grid_areas) 
    ICHECCCLMSON_mean = ICHECCCLMSON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECCCLMSON_grid_areas)
    ICHECKNMISON_mean = ICHECKNMISON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECKNMISON_grid_areas)
    ICHECMPISON_mean = ICHECMPISON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECMPISON_grid_areas)
    ICHECSMHISON_mean = ICHECSMHISON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECSMHISON_grid_areas)
    IPSLSON_mean = IPSLSON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=IPSLSON_grid_areas)
    MIROCSON_mean = MIROCSON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MIROCSON_grid_areas)
    MOHCCCLMSON_mean = MOHCCCLMSON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCCCLMSON_grid_areas)
    MOHCKNMISON_mean = MOHCKNMISON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCKNMISON_grid_areas)
    MOHCSMHISON_mean = MOHCSMHISON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCSMHISON_grid_areas)
    MPICCLMSON_mean = MPICCLMSON.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=MPICCLMSON_grid_areas)                                              
    MPIREMOSON_mean = MPIREMOSON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIREMOSON_grid_areas)                                              
    MPISMHISON_mean = MPISMHISON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPISMHISON_grid_areas)
    NCCDMISON_mean = NCCDMISON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NCCDMISON_grid_areas)
    NCCSMHISON_mean = NCCSMHISON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NCCSMHISON_grid_areas) 
    NOAASON_mean = NOAASON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NOAASON_grid_areas)
    
    CanESM2SON_mean = CanESM2SON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CanESM2SON_grid_areas)                        
    CNRMGSON_mean = CNRMGSON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMGSON_grid_areas) 
    MK3SON_mean = MK3SON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MK3SON_grid_areas) 
    EARTHSON_mean = EARTHSON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=EARTHSON_grid_areas)  
    EARTH3SON_mean = EARTH3SON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=EARTH3SON_grid_areas)
    GFDLSON_mean = GFDLSON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=GFDLSON_grid_areas)
    HadGEM2SON_mean = HadGEM2SON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=HadGEM2SON_grid_areas)
    IPSLGSON_mean = IPSLGSON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=IPSLGSON_grid_areas)
    MIROCGSON_mean = MIROCGSON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MIROCGSON_grid_areas)         
    MPISON_mean = MPISON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPISON_grid_areas)
    NorESM1SON_mean = NorESM1SON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NorESM1SON_grid_areas)
    
    CRUESON_mean = CRUESON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CRUESON_grid_areas) 
    CRUSON_mean = CRUSON.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CRUSON_grid_areas) 
        
    CCCmaDJF_mean = CCCmaDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CCCmaDJF_grid_areas)                                        
    CLMcomDJF_mean = CLMcomDJF.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=CLMcomDJF_grid_areas)
    DMIDJF_mean = DMIDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=DMIDJF_grid_areas)     
    KNMIDJF_mean = KNMIDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=KNMIDJF_grid_areas)
    MPIEDJF_mean = MPIEDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIEDJF_grid_areas)
    SMHIDJF_mean = SMHIDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=SMHIDJF_grid_areas)
    
    CCCmaCanRCMDJF_mean = CCCmaCanRCMDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CCCmaCanRCMDJF_grid_areas)  
    CCCmaSMHIDJF_mean = CCCmaSMHIDJF.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=CCCmaSMHIDJF_grid_areas)
    CNRMDJF_mean = CNRMDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMDJF_grid_areas)                                               
    CNRMSMHIDJF_mean = CNRMSMHIDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMSMHIDJF_grid_areas)  
    CSIRODJF_mean = CSIRODJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CSIRODJF_grid_areas)
    ICHECDMIDJF_mean = ICHECDMIDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECDMIDJF_grid_areas) 
    ICHECCCLMDJF_mean = ICHECCCLMDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECCCLMDJF_grid_areas)
    ICHECKNMIDJF_mean = ICHECKNMIDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECKNMIDJF_grid_areas)
    ICHECMPIDJF_mean = ICHECMPIDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECMPIDJF_grid_areas)
    ICHECSMHIDJF_mean = ICHECSMHIDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECSMHIDJF_grid_areas)
    IPSLDJF_mean = IPSLDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=IPSLDJF_grid_areas)
    MIROCDJF_mean = MIROCDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MIROCDJF_grid_areas)
    MOHCCCLMDJF_mean = MOHCCCLMDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCCCLMDJF_grid_areas)
    MOHCKNMIDJF_mean = MOHCKNMIDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCKNMIDJF_grid_areas)
    MOHCSMHIDJF_mean = MOHCSMHIDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCSMHIDJF_grid_areas)
    MPICCLMDJF_mean = MPICCLMDJF.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=MPICCLMDJF_grid_areas)                                              
    MPIREMODJF_mean = MPIREMODJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIREMODJF_grid_areas)                                              
    MPISMHIDJF_mean = MPISMHIDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPISMHIDJF_grid_areas)
    NCCDMIDJF_mean = NCCDMIDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NCCDMIDJF_grid_areas)
    NCCSMHIDJF_mean = NCCSMHIDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NCCSMHIDJF_grid_areas) 
    NOAADJF_mean = NOAADJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NOAADJF_grid_areas)
    
    CanESM2DJF_mean = CanESM2DJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CanESM2DJF_grid_areas)     
    CNRMGDJF_mean = CNRMGDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMGDJF_grid_areas)           
    MK3DJF_mean = MK3DJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MK3DJF_grid_areas)  
    EARTHDJF_mean = EARTHDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=EARTHDJF_grid_areas)  
    EARTH3DJF_mean = EARTH3DJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=EARTH3DJF_grid_areas)  
    GFDLDJF_mean = GFDLDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=GFDLDJF_grid_areas)
    HadGEM2DJF_mean = HadGEM2DJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=HadGEM2DJF_grid_areas)
    IPSLGDJF_mean = IPSLGDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=IPSLGDJF_grid_areas)
    MIROCGDJF_mean = MIROCGDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MIROCGDJF_grid_areas) 
    MPIDJF_mean = MPIDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIDJF_grid_areas)
    NorESM1DJF_mean = NorESM1DJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NorESM1DJF_grid_areas)
    
    CRUEDJF_mean = CRUEDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CRUEDJF_grid_areas) 
    CRUDJF_mean = CRUDJF.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CRUDJF_grid_areas) 
        
    CCCmaMAM_mean = CCCmaMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CCCmaMAM_grid_areas)                                        
    CLMcomMAM_mean = CLMcomMAM.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=CLMcomMAM_grid_areas)
    DMIMAM_mean = DMIMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=DMIMAM_grid_areas)     
    KNMIMAM_mean = KNMIMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=KNMIMAM_grid_areas)
    MPIEMAM_mean = MPIEMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIEMAM_grid_areas)
    SMHIMAM_mean = SMHIMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=SMHIMAM_grid_areas)
    
    CCCmaCanRCMMAM_mean = CCCmaCanRCMMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CCCmaCanRCMMAM_grid_areas)  
    CCCmaSMHIMAM_mean = CCCmaSMHIMAM.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=CCCmaSMHIMAM_grid_areas)
    CNRMMAM_mean = CNRMMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMMAM_grid_areas)                                               
    CNRMSMHIMAM_mean = CNRMSMHIMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMSMHIMAM_grid_areas)  
    CSIROMAM_mean = CSIROMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CSIROMAM_grid_areas)
    ICHECDMIMAM_mean = ICHECDMIMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECDMIMAM_grid_areas) 
    ICHECCCLMMAM_mean = ICHECCCLMMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECCCLMMAM_grid_areas)
    ICHECKNMIMAM_mean = ICHECKNMIMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECKNMIMAM_grid_areas)
    ICHECMPIMAM_mean = ICHECMPIMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECMPIMAM_grid_areas)
    ICHECSMHIMAM_mean = ICHECSMHIMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECSMHIMAM_grid_areas)
    IPSLMAM_mean = IPSLMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=IPSLMAM_grid_areas)
    MIROCMAM_mean = MIROCMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MIROCMAM_grid_areas)
    MOHCCCLMMAM_mean = MOHCCCLMMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCCCLMMAM_grid_areas)
    MOHCKNMIMAM_mean = MOHCKNMIMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCKNMIMAM_grid_areas)
    MOHCSMHIMAM_mean = MOHCSMHIMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCSMHIMAM_grid_areas)
    MPICCLMMAM_mean = MPICCLMMAM.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=MPICCLMMAM_grid_areas)                                              
    MPIREMOMAM_mean = MPIREMOMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIREMOMAM_grid_areas)                                              
    MPISMHIMAM_mean = MPISMHIMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPISMHIMAM_grid_areas)
    NCCDMIMAM_mean = NCCDMIMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NCCDMIMAM_grid_areas)
    NCCSMHIMAM_mean = NCCSMHIMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NCCSMHIMAM_grid_areas) 
    NOAAMAM_mean = NOAAMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NOAAMAM_grid_areas)
    
    CanESM2MAM_mean = CanESM2MAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CanESM2MAM_grid_areas) 
    CNRMGMAM_mean = CNRMGMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMGMAM_grid_areas)     
    MK3MAM_mean = MK3MAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MK3MAM_grid_areas) 
    EARTHMAM_mean = EARTHMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=EARTHMAM_grid_areas) 
    EARTH3MAM_mean = EARTH3MAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=EARTH3MAM_grid_areas)
    GFDLMAM_mean = GFDLMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=GFDLMAM_grid_areas)
    HadGEM2MAM_mean = HadGEM2MAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=HadGEM2MAM_grid_areas)
    IPSLGMAM_mean = IPSLGMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=IPSLGMAM_grid_areas)
    MIROCGMAM_mean = MIROCGMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MIROCGMAM_grid_areas)            
    MPIMAM_mean = MPIMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIMAM_grid_areas)
    NorESM1MAM_mean = NorESM1MAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NorESM1MAM_grid_areas)
       
    CRUEMAM_mean = CRUEMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CRUEMAM_grid_areas) 
    CRUMAM_mean = CRUMAM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CRUMAM_grid_areas) 
    
    CCCmaJJA_mean = CCCmaJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CCCmaJJA_grid_areas)                                        
    CLMcomJJA_mean = CLMcomJJA.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=CLMcomJJA_grid_areas)
    DMIJJA_mean = DMIJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=DMIJJA_grid_areas)     
    KNMIJJA_mean = KNMIJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=KNMIJJA_grid_areas)
    MPIEJJA_mean = MPIEJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIEJJA_grid_areas)
    SMHIJJA_mean = SMHIJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=SMHIJJA_grid_areas)
    
    CCCmaCanRCMJJA_mean = CCCmaCanRCMJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CCCmaCanRCMJJA_grid_areas)  
    CCCmaSMHIJJA_mean = CCCmaSMHIJJA.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=CCCmaSMHIJJA_grid_areas)
    CNRMJJA_mean = CNRMJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMJJA_grid_areas)                                               
    CNRMSMHIJJA_mean = CNRMSMHIJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMSMHIJJA_grid_areas)  
    CSIROJJA_mean = CSIROJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CSIROJJA_grid_areas)
    ICHECDMIJJA_mean = ICHECDMIJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECDMIJJA_grid_areas) 
    ICHECCCLMJJA_mean = ICHECCCLMJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECCCLMJJA_grid_areas)
    ICHECKNMIJJA_mean = ICHECKNMIJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECKNMIJJA_grid_areas)
    ICHECMPIJJA_mean = ICHECMPIJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECMPIJJA_grid_areas)
    ICHECSMHIJJA_mean = ICHECSMHIJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECSMHIJJA_grid_areas)
    IPSLJJA_mean = IPSLJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=IPSLJJA_grid_areas)
    MIROCJJA_mean = MIROCJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MIROCJJA_grid_areas)
    MOHCCCLMJJA_mean = MOHCCCLMJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCCCLMJJA_grid_areas)
    MOHCKNMIJJA_mean = MOHCKNMIJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCKNMIJJA_grid_areas)
    MOHCSMHIJJA_mean = MOHCSMHIJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCSMHIJJA_grid_areas)
    MPICCLMJJA_mean = MPICCLMJJA.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=MPICCLMJJA_grid_areas)                                              
    MPIREMOJJA_mean = MPIREMOJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIREMOJJA_grid_areas)                                              
    MPISMHIJJA_mean = MPISMHIJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPISMHIJJA_grid_areas)
    NCCDMIJJA_mean = NCCDMIJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NCCDMIJJA_grid_areas)
    NCCSMHIJJA_mean = NCCSMHIJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NCCSMHIJJA_grid_areas) 
    NOAAJJA_mean = NOAAJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NOAAJJA_grid_areas)
    
    CanESM2JJA_mean = CanESM2JJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CanESM2JJA_grid_areas)  
    CNRMGJJA_mean = CNRMGJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMGJJA_grid_areas)    
    MK3JJA_mean = MK3JJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MK3JJA_grid_areas)
    EARTHJJA_mean = EARTHJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=EARTHJJA_grid_areas)
    EARTH3JJA_mean = EARTH3JJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=EARTH3JJA_grid_areas)
    GFDLJJA_mean = GFDLJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=GFDLJJA_grid_areas)
    HadGEM2JJA_mean = HadGEM2JJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=HadGEM2JJA_grid_areas)
    IPSLGJJA_mean = IPSLGJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=IPSLGJJA_grid_areas)
    MIROCGJJA_mean = MIROCGJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MIROCGJJA_grid_areas)         
    MPIJJA_mean = MPIJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIJJA_grid_areas)
    NorESM1JJA_mean = NorESM1JJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NorESM1JJA_grid_areas)
       
    CRUEJJA_mean = CRUEJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CRUEJJA_grid_areas)   
    CRUJJA_mean = CRUJJA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CRUJJA_grid_areas) 
        
    CCCmaYR_mean = CCCmaYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CCCmaYR_grid_areas) 
    CLMcomYR_mean = CLMcomYR.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=CLMcomYR_grid_areas)
    DMIYR_mean = DMIYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=DMIYR_grid_areas)     
    KNMIYR_mean = KNMIYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=KNMIYR_grid_areas)
    MPIEYR_mean = MPIEYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIEYR_grid_areas)
    SMHIYR_mean = SMHIYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=SMHIYR_grid_areas)
    
    CCCmaCanRCMYR_mean = CCCmaCanRCMYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CCCmaCanRCMYR_grid_areas)   
    CCCmaSMHIYR_mean = CCCmaSMHIYR.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=CCCmaSMHIYR_grid_areas)
    CNRMYR_mean = CNRMYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMYR_grid_areas)                                               
    CNRMSMHIYR_mean = CNRMSMHIYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMSMHIYR_grid_areas)  
    CSIROYR_mean = CSIROYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CSIROYR_grid_areas)
    ICHECDMIYR_mean = ICHECDMIYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECDMIYR_grid_areas) 
    ICHECCCLMYR_mean = ICHECCCLMYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECCCLMYR_grid_areas)
    ICHECKNMIYR_mean = ICHECKNMIYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECKNMIYR_grid_areas)
    ICHECMPIYR_mean = ICHECMPIYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECMPIYR_grid_areas)
    ICHECSMHIYR_mean = ICHECSMHIYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECSMHIYR_grid_areas)
    IPSLYR_mean = IPSLYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=IPSLYR_grid_areas)
    MIROCYR_mean = MIROCYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MIROCYR_grid_areas)
    MOHCCCLMYR_mean = MOHCCCLMYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCCCLMYR_grid_areas)
    MOHCKNMIYR_mean = MOHCKNMIYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCKNMIYR_grid_areas)
    MOHCSMHIYR_mean = MOHCSMHIYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCSMHIYR_grid_areas)
    MPICCLMYR_mean = MPICCLMYR.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=MPICCLMYR_grid_areas)                                              
    MPIREMOYR_mean = MPIREMOYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIREMOYR_grid_areas)                                              
    MPISMHIYR_mean = MPISMHIYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPISMHIYR_grid_areas)
    NCCDMIYR_mean = NCCDMIYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NCCDMIYR_grid_areas)
    NCCSMHIYR_mean = NCCSMHIYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NCCSMHIYR_grid_areas) 
    NOAAYR_mean = NOAAYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NOAAYR_grid_areas)
    
    CanESM2YR_mean = CanESM2YR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CanESM2YR_grid_areas)   
    CNRMGYR_mean = CNRMGYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMGYR_grid_areas)
    MK3YR_mean = MK3YR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MK3YR_grid_areas)
    EARTHYR_mean = EARTHYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=EARTHYR_grid_areas) 
    EARTH3YR_mean = EARTH3YR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=EARTH3YR_grid_areas) 
    GFDLYR_mean = GFDLYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=GFDLYR_grid_areas)
    HadGEM2YR_mean = HadGEM2YR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=HadGEM2YR_grid_areas)
    IPSLGYR_mean = IPSLGYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=IPSLGYR_grid_areas)
    MIROCGYR_mean = MIROCGYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MIROCGYR_grid_areas) 
    MPIYR_mean = MPIYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIYR_grid_areas)
    NorESM1YR_mean = NorESM1YR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NorESM1YR_grid_areas)
   
    CRUEYR_mean = CRUEYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CRUEYR_grid_areas) 
    CRUYR_mean = CRUYR.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CRUYR_grid_areas) 
    
    #create averages
    AverageSONEY = (CCCmaSON_mean.data + CLMcomSON_mean.data + DMISON_mean.data + KNMISON_mean.data + MPIESON_mean.data + SMHISON_mean.data)/6.
    AverageDJFEY = (CCCmaDJF_mean.data + CLMcomDJF_mean.data + DMIDJF_mean.data + KNMIDJF_mean.data + MPIEDJF_mean.data + SMHIDJF_mean.data)/6.
    AverageMAMEY = (CCCmaMAM_mean.data + CLMcomMAM_mean.data + DMIMAM_mean.data + KNMIMAM_mean.data + MPIEMAM_mean.data + SMHIMAM_mean.data)/6.
    AverageJJAEY = (CCCmaJJA_mean.data + CLMcomJJA_mean.data + DMIJJA_mean.data + KNMIJJA_mean.data + MPIEJJA_mean.data + SMHIJJA_mean.data)/6.
    AverageEY = (CCCmaYR_mean.data + CLMcomYR_mean.data + DMIYR_mean.data + KNMIYR_mean.data + MPIEYR_mean.data + SMHIYR_mean.data)/6.
    
    AverageSONRY = (CCCmaCanRCMSON_mean.data + CCCmaSMHISON_mean.data + CNRMSON_mean.data + CNRMSMHISON_mean.data + CSIROSON_mean.data + ICHECDMISON_mean.data + ICHECCCLMSON_mean.data + ICHECKNMISON_mean.data + ICHECMPISON_mean.data + ICHECSMHISON_mean.data + IPSLSON_mean.data + MIROCSON_mean.data + MOHCCCLMSON_mean.data + MOHCKNMISON_mean.data + MOHCSMHISON_mean.data + MPICCLMSON_mean.data + MPIREMOSON_mean.data + MPISMHISON_mean.data + NCCDMISON_mean.data + NCCSMHISON_mean.data + NOAASON_mean.data)/21.
    AverageDJFRY = (CCCmaCanRCMDJF_mean.data + CCCmaSMHIDJF_mean.data + CNRMDJF_mean.data + CNRMSMHIDJF_mean.data + CSIRODJF_mean.data + ICHECDMIDJF_mean.data + ICHECCCLMDJF_mean.data + ICHECKNMIDJF_mean.data + ICHECMPIDJF_mean.data + ICHECSMHIDJF_mean.data + IPSLDJF_mean.data + MIROCDJF_mean.data + MOHCCCLMDJF_mean.data + MOHCKNMIDJF_mean.data + MOHCSMHIDJF_mean.data + MPICCLMDJF_mean.data + MPIREMODJF_mean.data + MPISMHIDJF_mean.data + NCCDMIDJF_mean.data + NCCSMHIDJF_mean.data + NOAADJF_mean.data)/21.
    AverageMAMRY = (CCCmaCanRCMMAM_mean.data + CCCmaSMHIMAM_mean.data + CNRMMAM_mean.data + CNRMSMHIMAM_mean.data + CSIROMAM_mean.data + ICHECDMIMAM_mean.data + ICHECCCLMMAM_mean.data + ICHECKNMIMAM_mean.data + ICHECMPIMAM_mean.data + ICHECSMHIMAM_mean.data + IPSLMAM_mean.data + MIROCMAM_mean.data + MOHCCCLMMAM_mean.data + MOHCKNMIMAM_mean.data + MOHCSMHIMAM_mean.data + MPICCLMMAM_mean.data + MPIREMOMAM_mean.data + MPISMHIMAM_mean.data + NCCDMIMAM_mean.data + NCCSMHIMAM_mean.data + NOAAMAM_mean.data)/21.
    AverageJJARY = (CCCmaCanRCMJJA_mean.data + CCCmaSMHIJJA_mean.data + CNRMJJA_mean.data + CNRMSMHIJJA_mean.data + CSIROJJA_mean.data + ICHECDMIJJA_mean.data + ICHECCCLMJJA_mean.data + ICHECKNMIJJA_mean.data + ICHECMPIJJA_mean.data + ICHECSMHIJJA_mean.data + IPSLJJA_mean.data + MIROCJJA_mean.data + MOHCCCLMJJA_mean.data + MOHCKNMIJJA_mean.data + MOHCSMHIJJA_mean.data + MPICCLMJJA_mean.data + MPIREMOJJA_mean.data + MPISMHIJJA_mean.data + NCCDMIJJA_mean.data + NCCSMHIJJA_mean.data + NOAAJJA_mean.data)/21.
    AverageRY = (CCCmaCanRCMYR_mean.data + CCCmaSMHIYR_mean.data + CNRMYR_mean.data + CNRMSMHIYR_mean.data + CSIROYR_mean.data + ICHECDMIYR_mean.data + ICHECCCLMYR_mean.data + ICHECKNMIYR_mean.data + ICHECMPIYR_mean.data + ICHECSMHIYR_mean.data + IPSLYR_mean.data + MIROCYR_mean.data + MOHCCCLMYR_mean.data + MOHCKNMIYR_mean.data + MOHCSMHIYR_mean.data + MPICCLMYR_mean.data + MPIREMOYR_mean.data + MPISMHIYR_mean.data + NCCDMIYR_mean.data + NCCSMHIYR_mean.data + NOAAYR_mean.data)/21.
    
    AverageSONGY = (CanESM2SON_mean.data + CNRMGSON_mean.data + MK3SON_mean.data + EARTHSON_mean.data + EARTH3SON_mean.data + GFDLSON_mean.data + HadGEM2SON_mean.data + IPSLGSON_mean.data + MIROCGSON_mean.data + MPISON_mean.data + NorESM1SON_mean.data)/11.
    AverageDJFGY = (CanESM2DJF_mean.data + CNRMGDJF_mean.data + MK3DJF_mean.data + EARTHDJF_mean.data + EARTH3DJF_mean.data + GFDLDJF_mean.data + HadGEM2DJF_mean.data + IPSLGDJF_mean.data + MIROCGDJF_mean.data + MPIDJF_mean.data + NorESM1DJF_mean.data)/11.
    AverageMAMGY = (CanESM2MAM_mean.data + CNRMGMAM_mean.data + MK3MAM_mean.data + EARTHMAM_mean.data + EARTH3MAM_mean.data + GFDLMAM_mean.data + HadGEM2MAM_mean.data + IPSLGMAM_mean.data + MIROCGMAM_mean.data + MPIMAM_mean.data + NorESM1MAM_mean.data)/11.
    AverageJJAGY = (CanESM2JJA_mean.data + CNRMGJJA_mean.data + MK3JJA_mean.data + EARTHJJA_mean.data + EARTH3JJA_mean.data + GFDLJJA_mean.data + HadGEM2JJA_mean.data + IPSLGJJA_mean.data + MIROCGJJA_mean.data + MPIJJA_mean.data + NorESM1JJA_mean.data)/11.
    AverageGY = (CanESM2YR_mean.data + CNRMGYR_mean.data + MK3YR_mean.data + EARTHYR_mean.data + EARTH3YR_mean.data + GFDLYR_mean.data + HadGEM2YR_mean.data + IPSLGYR_mean.data + MIROCGYR_mean.data + MPIYR_mean.data + NorESM1YR_mean.data)/11.
       
    ObsESONY = (CRUESON_mean.data)
    ObsEDJFY = (CRUEDJF_mean.data)
    ObsEMAMY = (CRUEMAM_mean.data)
    ObsEJJAY = (CRUEJJA_mean.data)
    ObsEY = (CRUEYR_mean.data)
    
    ObsSONY = (CRUSON_mean.data)
    ObsDJFY = (CRUDJF_mean.data)
    ObsMAMY = (CRUMAM_mean.data)
    ObsJJAY = (CRUJJA_mean.data)
    ObsY = (CRUYR_mean.data)
    
       
    #-------------------------------------------------------------------------
    #PART 4: Calculate Standard Deviation
    
    AverageSONEY_std=np.std(AverageSONEY)
    CCCmaSON_std=np.std(CCCmaSON_mean.data)
    CLMcomSON_std=np.std(CLMcomSON_mean.data)
    DMISON_std=np.std(DMISON_mean.data)
    KNMISON_std=np.std(KNMISON_mean.data)
    MPIESON_std=np.std(MPIESON_mean.data)
    SMHISON_std=np.std(SMHISON_mean.data)
    
    AverageDJFEY_std=np.std(AverageDJFEY)
    CCCmaDJF_std=np.std(CCCmaDJF_mean.data) 
    CLMcomDJF_std=np.std(CLMcomDJF_mean.data)
    DMIDJF_std=np.std(DMIDJF_mean.data)
    KNMIDJF_std=np.std(KNMIDJF_mean.data)
    MPIEDJF_std=np.std(MPIEDJF_mean.data)
    SMHIDJF_std=np.std(SMHIDJF_mean.data)
    
    AverageMAMEY_std=np.std(AverageMAMEY)
    CCCmaMAM_std=np.std(CCCmaMAM_mean.data) 
    CLMcomMAM_std=np.std(CLMcomMAM_mean.data)
    DMIMAM_std=np.std(DMIMAM_mean.data)
    KNMIMAM_std=np.std(KNMIMAM_mean.data) 
    MPIEMAM_std=np.std(MPIEMAM_mean.data)  
    SMHIMAM_std=np.std(SMHIMAM_mean.data) 
    
    AverageJJAEY_std=np.std(AverageJJAEY)  
    CCCmaJJA_std=np.std(CCCmaJJA_mean.data)  
    CLMcomJJA_std=np.std(CLMcomJJA_mean.data)  
    DMIJJA_std=np.std(DMIJJA_mean.data)  
    KNMIJJA_std=np.std(KNMIJJA_mean.data)  
    MPIEJJA_std=np.std(MPIEJJA_mean.data)  
    SMHIJJA_std=np.std(SMHIJJA_mean.data) 
    
    AverageEY_std=np.std(AverageEY) 
    CCCmaYR_std=np.std(CCCmaYR_mean.data)
    CLMcomYR_std=np.std(CLMcomYR_mean.data)  
    DMIYR_std=np.std(DMIYR_mean.data)  
    KNMIYR_std=np.std(KNMIYR_mean.data)  
    MPIEYR_std=np.std(MPIEYR_mean.data)  
    SMHIYR_std=np.std(SMHIYR_mean.data)
     
    AverageSONRY_std=np.std(AverageSONRY)
    CCCmaCanRCMSON_std=np.std(CCCmaCanRCMSON_mean.data)
    CCCmaSMHISON_std=np.std(CCCmaSMHISON_mean.data) 
    CNRMSON_std=np.std(CNRMSON_mean.data) 
    CNRMSMHISON_std=np.std(CNRMSMHISON_mean.data) 
    CSIROSON_std=np.std(CSIROSON_mean.data) 
    ICHECDMISON_std=np.std(ICHECDMISON_mean.data) 
    ICHECCCLMSON_std=np.std(ICHECCCLMSON_mean.data) 
    ICHECKNMISON_std=np.std(ICHECKNMISON_mean.data) 
    ICHECMPISON_std=np.std(ICHECMPISON_mean.data)
    ICHECSMHISON_std=np.std(ICHECSMHISON_mean.data) 
    IPSLSON_std=np.std(IPSLSON_mean.data) 
    MIROCSON_std=np.std(MIROCSON_mean.data) 
    MOHCCCLMSON_std=np.std(MOHCCCLMSON_mean.data) 
    MOHCKNMISON_std=np.std(MOHCKNMISON_mean.data) 
    MOHCSMHISON_std=np.std(MOHCSMHISON_mean.data) 
    MPICCLMSON_std=np.std(MPICCLMSON_mean.data) 
    MPIREMOSON_std=np.std(MPIREMOSON_mean.data) 
    MPISMHISON_std=np.std(MPISMHISON_mean.data) 
    NCCDMISON_std=np.std(NCCDMISON_mean.data) 
    NCCSMHISON_std=np.std(NCCSMHISON_mean.data) 
    NOAASON_std=np.std(NOAASON_mean.data)
     
    AverageDJFRY_std=np.std(AverageDJFRY)
    CCCmaCanRCMDJF_std=np.std(CCCmaCanRCMDJF_mean.data)
    CCCmaSMHIDJF_std=np.std(CCCmaSMHIDJF_mean.data) 
    CNRMDJF_std=np.std(CNRMDJF_mean.data) 
    CNRMSMHIDJF_std=np.std(CNRMSMHIDJF_mean.data) 
    CSIRODJF_std=np.std(CSIRODJF_mean.data) 
    ICHECDMIDJF_std=np.std(ICHECDMIDJF_mean.data) 
    ICHECCCLMDJF_std=np.std(ICHECCCLMDJF_mean.data) 
    ICHECKNMIDJF_std=np.std(ICHECKNMIDJF_mean.data) 
    ICHECMPIDJF_std=np.std(ICHECMPIDJF_mean.data)
    ICHECSMHIDJF_std=np.std(ICHECSMHIDJF_mean.data) 
    IPSLDJF_std=np.std(IPSLDJF_mean.data) 
    MIROCDJF_std=np.std(MIROCDJF_mean.data) 
    MOHCCCLMDJF_std=np.std(MOHCCCLMDJF_mean.data) 
    MOHCKNMIDJF_std=np.std(MOHCKNMIDJF_mean.data) 
    MOHCSMHIDJF_std=np.std(MOHCSMHIDJF_mean.data) 
    MPICCLMDJF_std=np.std(MPICCLMDJF_mean.data) 
    MPIREMODJF_std=np.std(MPIREMODJF_mean.data) 
    MPISMHIDJF_std=np.std(MPISMHIDJF_mean.data) 
    NCCDMIDJF_std=np.std(NCCDMIDJF_mean.data) 
    NCCSMHIDJF_std=np.std(NCCSMHIDJF_mean.data) 
    NOAADJF_std=np.std(NOAADJF_mean.data)
     
    AverageMAMRY_std=np.std(AverageMAMRY)
    CCCmaCanRCMMAM_std=np.std(CCCmaCanRCMMAM_mean.data)
    CCCmaSMHIMAM_std=np.std(CCCmaSMHIMAM_mean.data) 
    CNRMMAM_std=np.std(CNRMMAM_mean.data) 
    CNRMSMHIMAM_std=np.std(CNRMSMHIMAM_mean.data) 
    CSIROMAM_std=np.std(CSIROMAM_mean.data) 
    ICHECDMIMAM_std=np.std(ICHECDMIMAM_mean.data) 
    ICHECCCLMMAM_std=np.std(ICHECCCLMMAM_mean.data) 
    ICHECKNMIMAM_std=np.std(ICHECKNMIMAM_mean.data) 
    ICHECMPIMAM_std=np.std(ICHECMPIMAM_mean.data)
    ICHECSMHIMAM_std=np.std(ICHECSMHIMAM_mean.data) 
    IPSLMAM_std=np.std(IPSLMAM_mean.data) 
    MIROCMAM_std=np.std(MIROCMAM_mean.data) 
    MOHCCCLMMAM_std=np.std(MOHCCCLMMAM_mean.data) 
    MOHCKNMIMAM_std=np.std(MOHCKNMIMAM_mean.data) 
    MOHCSMHIMAM_std=np.std(MOHCSMHIMAM_mean.data) 
    MPICCLMMAM_std=np.std(MPICCLMMAM_mean.data) 
    MPIREMOMAM_std=np.std(MPIREMOMAM_mean.data) 
    MPISMHIMAM_std=np.std(MPISMHIMAM_mean.data) 
    NCCDMIMAM_std=np.std(NCCDMIMAM_mean.data) 
    NCCSMHIMAM_std=np.std(NCCSMHIMAM_mean.data) 
    NOAAMAM_std=np.std(NOAAMAM_mean.data)
    
    AverageJJARY_std=np.std(AverageJJARY)
    CCCmaCanRCMJJA_std=np.std(CCCmaCanRCMJJA_mean.data)
    CCCmaSMHIJJA_std=np.std(CCCmaSMHIJJA_mean.data) 
    CNRMJJA_std=np.std(CNRMJJA_mean.data) 
    CNRMSMHIJJA_std=np.std(CNRMSMHIJJA_mean.data) 
    CSIROJJA_std=np.std(CSIROJJA_mean.data) 
    ICHECDMIJJA_std=np.std(ICHECDMIJJA_mean.data) 
    ICHECCCLMJJA_std=np.std(ICHECCCLMJJA_mean.data) 
    ICHECKNMIJJA_std=np.std(ICHECKNMIJJA_mean.data) 
    ICHECMPIJJA_std=np.std(ICHECMPIJJA_mean.data)
    ICHECSMHIJJA_std=np.std(ICHECSMHIJJA_mean.data) 
    IPSLJJA_std=np.std(IPSLJJA_mean.data) 
    MIROCJJA_std=np.std(MIROCJJA_mean.data) 
    MOHCCCLMJJA_std=np.std(MOHCCCLMJJA_mean.data) 
    MOHCKNMIJJA_std=np.std(MOHCKNMIJJA_mean.data) 
    MOHCSMHIJJA_std=np.std(MOHCSMHIJJA_mean.data) 
    MPICCLMJJA_std=np.std(MPICCLMJJA_mean.data) 
    MPIREMOJJA_std=np.std(MPIREMOJJA_mean.data) 
    MPISMHIJJA_std=np.std(MPISMHIJJA_mean.data) 
    NCCDMIJJA_std=np.std(NCCDMIJJA_mean.data) 
    NCCSMHIJJA_std=np.std(NCCSMHIJJA_mean.data) 
    NOAAJJA_std=np.std(NOAAJJA_mean.data)
    
    AverageRY_std=np.std(AverageRY)
    CCCmaCanRCMYR_std=np.std(CCCmaCanRCMYR_mean.data)
    CCCmaSMHIYR_std=np.std(CCCmaSMHIYR_mean.data) 
    CNRMYR_std=np.std(CNRMYR_mean.data) 
    CNRMSMHIYR_std=np.std(CNRMSMHIYR_mean.data) 
    CSIROYR_std=np.std(CSIROYR_mean.data) 
    ICHECDMIYR_std=np.std(ICHECDMIYR_mean.data) 
    ICHECCCLMYR_std=np.std(ICHECCCLMYR_mean.data) 
    ICHECKNMIYR_std=np.std(ICHECKNMIYR_mean.data) 
    ICHECMPIYR_std=np.std(ICHECMPIYR_mean.data)
    ICHECSMHIYR_std=np.std(ICHECSMHIYR_mean.data) 
    IPSLYR_std=np.std(IPSLYR_mean.data) 
    MIROCYR_std=np.std(MIROCYR_mean.data) 
    MOHCCCLMYR_std=np.std(MOHCCCLMYR_mean.data) 
    MOHCKNMIYR_std=np.std(MOHCKNMIYR_mean.data) 
    MOHCSMHIYR_std=np.std(MOHCSMHIYR_mean.data) 
    MPICCLMYR_std=np.std(MPICCLMYR_mean.data) 
    MPIREMOYR_std=np.std(MPIREMOYR_mean.data) 
    MPISMHIYR_std=np.std(MPISMHIYR_mean.data) 
    NCCDMIYR_std=np.std(NCCDMIYR_mean.data) 
    NCCSMHIYR_std=np.std(NCCSMHIYR_mean.data) 
    NOAAYR_std=np.std(NOAAYR_mean.data)
    
    AverageSONGY_std=np.std(AverageSONGY)
    CanESM2SON_std=np.std(CanESM2SON_mean.data)
    CNRMGSON_std=np.std(CNRMGSON_mean.data)
    MK3SON_std=np.std(MK3SON_mean.data)
    EARTHSON_std=np.std(EARTHSON_mean.data)
    EARTH3SON_std=np.std(EARTH3SON_mean.data)
    GFDLSON_std=np.std(GFDLSON_mean.data) 
    HadGEM2SON_std=np.std(HadGEM2SON_mean.data) 
    IPSLGSON_std=np.std(IPSLGSON_mean.data) 
    MIROCGSON_std=np.std(MIROCGSON_mean.data) 
    MPISON_std=np.std(MPISON_mean.data) 
    NorESM1SON_std=np.std(NorESM1SON_mean.data)
 
    AverageDJFGY_std=np.std(AverageDJFGY) 
    CanESM2DJF_std=np.std(CanESM2DJF_mean.data) 
    CNRMGDJF_std=np.std(CNRMGDJF_mean.data) 
    MK3DJF_std=np.std(MK3DJF_mean.data) 
    EARTHDJF_std=np.std(EARTHDJF_mean.data) 
    EARTH3DJF_std=np.std(EARTH3DJF_mean.data) 
    GFDLDJF_std=np.std(GFDLDJF_mean.data) 
    HadGEM2DJF_std=np.std(HadGEM2DJF_mean.data) 
    IPSLGDJF_std=np.std(IPSLGDJF_mean.data) 
    MIROCGDJF_std=np.std(MIROCGDJF_mean.data) 
    MPIDJF_std=np.std(MPIDJF_mean.data) 
    NorESM1DJF_std=np.std(NorESM1DJF_mean.data)
     
    AverageMAMGY_std=np.std(AverageMAMGY) 
    CanESM2MAM_std=np.std(CanESM2MAM_mean.data) 
    CNRMGMAM_std=np.std(CNRMGMAM_mean.data) 
    MK3MAM_std=np.std(MK3MAM_mean.data) 
    EARTHMAM_std=np.std(EARTHMAM_mean.data) 
    EARTH3MAM_std=np.std(EARTH3MAM_mean.data)
    GFDLMAM_std=np.std(GFDLMAM_mean.data) 
    HadGEM2MAM_std=np.std(HadGEM2MAM_mean.data) 
    IPSLGMAM_std=np.std(IPSLGMAM_mean.data) 
    MIROCGMAM_std=np.std(MIROCGMAM_mean.data) 
    MPIMAM_std=np.std(MPIMAM_mean.data) 
    NorESM1MAM_std=np.std(NorESM1MAM_mean.data)
 
    AverageJJAGY_std=np.std(AverageJJAGY) 
    CanESM2JJA_std=np.std(CanESM2JJA_mean.data) 
    CNRMGJJA_std=np.std(CNRMGJJA_mean.data) 
    MK3JJA_std=np.std(MK3JJA_mean.data) 
    EARTHJJA_std=np.std(EARTHJJA_mean.data) 
    EARTH3JJA_std=np.std(EARTH3JJA_mean.data) 
    GFDLJJA_std=np.std(GFDLJJA_mean.data) 
    HadGEM2JJA_std=np.std(HadGEM2JJA_mean.data) 
    IPSLGJJA_std=np.std(IPSLGJJA_mean.data) 
    MIROCGJJA_std=np.std(MIROCGJJA_mean.data) 
    MPIJJA_std=np.std(MPIJJA_mean.data) 
    NorESM1JJA_std=np.std(NorESM1JJA_mean.data)
   
    AverageGY_std=np.std(AverageGY) 
    CanESM2YR_std=np.std(CanESM2YR_mean.data) 
    CNRMGYR_std=np.std(CNRMGYR_mean.data) 
    MK3YR_std=np.std(MK3YR_mean.data) 
    EARTHYR_std=np.std(EARTHYR_mean.data) 
    EARTH3YR_std=np.std(EARTH3YR_mean.data) 
    GFDLYR_std=np.std(GFDLYR_mean.data) 
    HadGEM2YR_std=np.std(HadGEM2YR_mean.data) 
    IPSLGYR_std=np.std(IPSLGYR_mean.data) 
    MIROCGYR_std=np.std(MIROCGYR_mean.data) 
    MPIYR_std=np.std(MPIYR_mean.data) 
    NorESM1YR_std=np.std(NorESM1YR_mean.data)
    
    ObsESONY_std=np.std(ObsESONY) 
    ObsEDJFY_std=np.std(ObsEDJFY) 
    ObsEMAMY_std=np.std(ObsEMAMY) 
    ObsEJJAY_std=np.std(ObsEJJAY) 
    ObsEY_std=np.std(ObsEY) 
    
    ObsSONY_std=np.std(ObsSONY) 
    ObsDJFY_std=np.std(ObsDJFY) 
    ObsMAMY_std=np.std(ObsMAMY) 
    ObsJJAY_std=np.std(ObsJJAY) 
    ObsY_std=np.std(ObsY)
    
    #-------------------------------------------------------------------------
    #PART 5: Calculate Pearson's Correlation Coefficient
    #use calculation in def Pearson
    
    AverageSONEY_PCC=pearson(AverageSONEY,ObsESONY)
    CCCmaSON_PCC=pearson(CCCmaSON_mean.data,ObsESONY)
    CLMcomSON_PCC=pearson(CLMcomSON_mean.data,ObsESONY)
    DMISON_PCC=pearson(DMISON_mean.data,ObsESONY)
    KNMISON_PCC=pearson(KNMISON_mean.data,ObsESONY)
    MPIESON_PCC=pearson(MPIESON_mean.data,ObsESONY)
    SMHISON_PCC=pearson(SMHISON_mean.data,ObsESONY)
    
    AverageDJFEY_PCC=pearson(AverageDJFEY,ObsEDJFY)
    CCCmaDJF_PCC=pearson(CCCmaDJF_mean.data,ObsEDJFY) 
    CLMcomDJF_PCC=pearson(CLMcomDJF_mean.data,ObsEDJFY)  
    DMIDJF_PCC=pearson(DMIDJF_mean.data,ObsEDJFY)
    KNMIDJF_PCC=pearson(KNMIDJF_mean.data,ObsEDJFY) 
    MPIEDJF_PCC=pearson(MPIEDJF_mean.data,ObsEDJFY)  
    SMHIDJF_PCC=pearson(SMHIDJF_mean.data,ObsEDJFY)  
    
    AverageMAMEY_PCC=pearson(AverageMAMEY,ObsEMAMY)
    CCCmaMAM_PCC=pearson(CCCmaMAM_mean.data,ObsEMAMY)
    CLMcomMAM_PCC=pearson(CLMcomMAM_mean.data,ObsEMAMY)
    DMIMAM_PCC=pearson(DMIMAM_mean.data,ObsEMAMY)
    KNMIMAM_PCC=pearson(KNMIMAM_mean.data,ObsEMAMY)
    MPIEMAM_PCC=pearson(MPIEMAM_mean.data,ObsEMAMY)
    SMHIMAM_PCC=pearson(SMHIMAM_mean.data,ObsEMAMY)  
    
    AverageJJAEY_PCC=pearson(AverageJJAEY,ObsEJJAY)
    CCCmaJJA_PCC=pearson(CCCmaJJA_mean.data,ObsEJJAY)
    CLMcomJJA_PCC=pearson(CLMcomJJA_mean.data,ObsEJJAY)
    DMIJJA_PCC=pearson(DMIJJA_mean.data,ObsEJJAY) 
    KNMIJJA_PCC=pearson(KNMIJJA_mean.data,ObsEJJAY)  
    MPIEJJA_PCC=pearson(MPIEJJA_mean.data,ObsEJJAY) 
    SMHIJJA_PCC=pearson(SMHIJJA_mean.data,ObsEJJAY)  
    
    AverageEY_PCC=pearson(AverageEY,ObsEY)
    CCCmaYR_PCC=pearson(CCCmaYR_mean.data,ObsEY)
    CLMcomYR_PCC=pearson(CLMcomYR_mean.data,ObsEY)
    DMIYR_PCC=pearson(DMIYR_mean.data,ObsEY) 
    KNMIYR_PCC=pearson(KNMIYR_mean.data,ObsEY)  
    MPIEYR_PCC=pearson(MPIEYR_mean.data,ObsEY) 
    SMHIYR_PCC=pearson(SMHIYR_mean.data,ObsEY)  
    
    AverageSONRY_PCC=pearson(AverageSONRY,ObsSONY)
    CCCmaCanRCMSON_PCC=pearson(CCCmaCanRCMSON_mean.data,ObsSONY)
    CCCmaSMHISON_PCC=pearson(CCCmaSMHISON_mean.data,ObsSONY)
    CNRMSON_PCC=pearson(CNRMSON_mean.data,ObsSONY)
    CNRMSMHISON_PCC=pearson(CNRMSMHISON_mean.data,ObsSONY)
    CSIROSON_PCC=pearson(CSIROSON_mean.data,ObsSONY) 
    ICHECDMISON_PCC=pearson(ICHECDMISON_mean.data,ObsSONY)  
    ICHECCCLMSON_PCC=pearson(ICHECCCLMSON_mean.data,ObsSONY) 
    ICHECKNMISON_PCC=pearson(ICHECKNMISON_mean.data,ObsSONY)
    ICHECMPISON_PCC=pearson(ICHECMPISON_mean.data,ObsSONY)
    ICHECSMHISON_PCC=pearson(ICHECSMHISON_mean.data,ObsSONY)
    IPSLSON_PCC=pearson(IPSLSON_mean.data,ObsSONY) 
    MIROCSON_PCC=pearson(MIROCSON_mean.data,ObsSONY)
    MOHCCCLMSON_PCC=pearson(MOHCCCLMSON_mean.data,ObsSONY) 
    MOHCKNMISON_PCC=pearson(MOHCKNMISON_mean.data,ObsSONY)  
    MOHCSMHISON_PCC=pearson(MOHCSMHISON_mean.data,ObsSONY) 
    MPICCLMSON_PCC=pearson(MPICCLMSON_mean.data,ObsSONY)
    MPIREMOSON_PCC=pearson(MPIREMOSON_mean.data,ObsSONY)
    MPISMHISON_PCC=pearson(MPISMHISON_mean.data,ObsSONY)  
    NCCDMISON_PCC=pearson(NCCDMISON_mean.data,ObsSONY) 
    NCCSMHISON_PCC=pearson(NCCSMHISON_mean.data,ObsSONY)  
    NOAASON_PCC=pearson(NOAASON_mean.data,ObsSONY)
    
    AverageDJFRY_PCC=pearson(AverageDJFRY,ObsDJFY)
    CCCmaCanRCMDJF_PCC=pearson(CCCmaCanRCMDJF_mean.data,ObsDJFY)
    CCCmaSMHIDJF_PCC=pearson(CCCmaSMHIDJF_mean.data,ObsDJFY)
    CNRMDJF_PCC=pearson(CNRMDJF_mean.data,ObsDJFY)
    CNRMSMHIDJF_PCC=pearson(CNRMSMHIDJF_mean.data,ObsDJFY)
    CSIRODJF_PCC=pearson(CSIRODJF_mean.data,ObsDJFY) 
    ICHECDMIDJF_PCC=pearson(ICHECDMIDJF_mean.data,ObsDJFY)  
    ICHECCCLMDJF_PCC=pearson(ICHECCCLMDJF_mean.data,ObsDJFY) 
    ICHECKNMIDJF_PCC=pearson(ICHECKNMIDJF_mean.data,ObsDJFY)
    ICHECMPIDJF_PCC=pearson(ICHECMPIDJF_mean.data,ObsDJFY)
    ICHECSMHIDJF_PCC=pearson(ICHECSMHIDJF_mean.data,ObsDJFY)
    IPSLDJF_PCC=pearson(IPSLDJF_mean.data,ObsDJFY) 
    MIROCDJF_PCC=pearson(MIROCDJF_mean.data,ObsDJFY)
    MOHCCCLMDJF_PCC=pearson(MOHCCCLMDJF_mean.data,ObsDJFY) 
    MOHCKNMIDJF_PCC=pearson(MOHCKNMIDJF_mean.data,ObsDJFY)  
    MOHCSMHIDJF_PCC=pearson(MOHCSMHIDJF_mean.data,ObsDJFY) 
    MPICCLMDJF_PCC=pearson(MPICCLMDJF_mean.data,ObsDJFY)
    MPIREMODJF_PCC=pearson(MPIREMODJF_mean.data,ObsDJFY)
    MPISMHIDJF_PCC=pearson(MPISMHIDJF_mean.data,ObsDJFY)  
    NCCDMIDJF_PCC=pearson(NCCDMIDJF_mean.data,ObsDJFY) 
    NCCSMHIDJF_PCC=pearson(NCCSMHIDJF_mean.data,ObsDJFY)  
    NOAADJF_PCC=pearson(NOAADJF_mean.data,ObsDJFY)
    
    AverageMAMRY_PCC=pearson(AverageMAMRY,ObsMAMY)
    CCCmaCanRCMMAM_PCC=pearson(CCCmaCanRCMMAM_mean.data,ObsMAMY)
    CCCmaSMHIMAM_PCC=pearson(CCCmaSMHIMAM_mean.data,ObsMAMY)
    CNRMMAM_PCC=pearson(CNRMMAM_mean.data,ObsMAMY)
    CNRMSMHIMAM_PCC=pearson(CNRMSMHIMAM_mean.data,ObsMAMY)
    CSIROMAM_PCC=pearson(CSIROMAM_mean.data,ObsMAMY) 
    ICHECDMIMAM_PCC=pearson(ICHECDMIMAM_mean.data,ObsMAMY)  
    ICHECCCLMMAM_PCC=pearson(ICHECCCLMMAM_mean.data,ObsMAMY) 
    ICHECKNMIMAM_PCC=pearson(ICHECKNMIMAM_mean.data,ObsMAMY)
    ICHECMPIMAM_PCC=pearson(ICHECMPIMAM_mean.data,ObsMAMY)
    ICHECSMHIMAM_PCC=pearson(ICHECSMHIMAM_mean.data,ObsMAMY)
    IPSLMAM_PCC=pearson(IPSLMAM_mean.data,ObsMAMY) 
    MIROCMAM_PCC=pearson(MIROCMAM_mean.data,ObsMAMY)
    MOHCCCLMMAM_PCC=pearson(MOHCCCLMMAM_mean.data,ObsMAMY) 
    MOHCKNMIMAM_PCC=pearson(MOHCKNMIMAM_mean.data,ObsMAMY)  
    MOHCSMHIMAM_PCC=pearson(MOHCSMHIMAM_mean.data,ObsMAMY) 
    MPICCLMMAM_PCC=pearson(MPICCLMMAM_mean.data,ObsMAMY)
    MPIREMOMAM_PCC=pearson(MPIREMOMAM_mean.data,ObsMAMY)
    MPISMHIMAM_PCC=pearson(MPISMHIMAM_mean.data,ObsMAMY)  
    NCCDMIMAM_PCC=pearson(NCCDMIMAM_mean.data,ObsMAMY) 
    NCCSMHIMAM_PCC=pearson(NCCSMHIMAM_mean.data,ObsMAMY)  
    NOAAMAM_PCC=pearson(NOAAMAM_mean.data,ObsMAMY)
    
    AverageJJARY_PCC=pearson(AverageJJARY,ObsJJAY)
    CCCmaCanRCMJJA_PCC=pearson(CCCmaCanRCMJJA_mean.data,ObsJJAY)
    CCCmaSMHIJJA_PCC=pearson(CCCmaSMHIJJA_mean.data,ObsJJAY)
    CNRMJJA_PCC=pearson(CNRMJJA_mean.data,ObsJJAY)
    CNRMSMHIJJA_PCC=pearson(CNRMSMHIJJA_mean.data,ObsJJAY)
    CSIROJJA_PCC=pearson(CSIROJJA_mean.data,ObsJJAY) 
    ICHECDMIJJA_PCC=pearson(ICHECDMIJJA_mean.data,ObsJJAY)  
    ICHECCCLMJJA_PCC=pearson(ICHECCCLMJJA_mean.data,ObsJJAY) 
    ICHECKNMIJJA_PCC=pearson(ICHECKNMIJJA_mean.data,ObsJJAY)
    ICHECMPIJJA_PCC=pearson(ICHECMPIJJA_mean.data,ObsJJAY)
    ICHECSMHIJJA_PCC=pearson(ICHECSMHIJJA_mean.data,ObsJJAY)
    IPSLJJA_PCC=pearson(IPSLJJA_mean.data,ObsJJAY) 
    MIROCJJA_PCC=pearson(MIROCJJA_mean.data,ObsJJAY)
    MOHCCCLMJJA_PCC=pearson(MOHCCCLMJJA_mean.data,ObsJJAY) 
    MOHCKNMIJJA_PCC=pearson(MOHCKNMIJJA_mean.data,ObsJJAY)  
    MOHCSMHIJJA_PCC=pearson(MOHCSMHIJJA_mean.data,ObsJJAY) 
    MPICCLMJJA_PCC=pearson(MPICCLMJJA_mean.data,ObsJJAY)
    MPIREMOJJA_PCC=pearson(MPIREMOJJA_mean.data,ObsJJAY)
    MPISMHIJJA_PCC=pearson(MPISMHIJJA_mean.data,ObsJJAY)  
    NCCDMIJJA_PCC=pearson(NCCDMIJJA_mean.data,ObsJJAY) 
    NCCSMHIJJA_PCC=pearson(NCCSMHIJJA_mean.data,ObsJJAY)  
    NOAAJJA_PCC=pearson(NOAAJJA_mean.data,ObsJJAY)
    
    AverageRY_PCC=pearson(AverageRY,ObsY)
    CCCmaCanRCMYR_PCC=pearson(CCCmaCanRCMYR_mean.data,ObsY)
    CCCmaSMHIYR_PCC=pearson(CCCmaSMHIYR_mean.data,ObsY)
    CNRMYR_PCC=pearson(CNRMYR_mean.data,ObsY)
    CNRMSMHIYR_PCC=pearson(CNRMSMHIYR_mean.data,ObsY)
    CSIROYR_PCC=pearson(CSIROYR_mean.data,ObsY) 
    ICHECDMIYR_PCC=pearson(ICHECDMIYR_mean.data,ObsY)  
    ICHECCCLMYR_PCC=pearson(ICHECCCLMYR_mean.data,ObsY) 
    ICHECKNMIYR_PCC=pearson(ICHECKNMIYR_mean.data,ObsY)
    ICHECMPIYR_PCC=pearson(ICHECMPIYR_mean.data,ObsY)
    ICHECSMHIYR_PCC=pearson(ICHECSMHIYR_mean.data,ObsY)
    IPSLYR_PCC=pearson(IPSLYR_mean.data,ObsY) 
    MIROCYR_PCC=pearson(MIROCYR_mean.data,ObsY)
    MOHCCCLMYR_PCC=pearson(MOHCCCLMYR_mean.data,ObsY) 
    MOHCKNMIYR_PCC=pearson(MOHCKNMIYR_mean.data,ObsY)  
    MOHCSMHIYR_PCC=pearson(MOHCSMHIYR_mean.data,ObsY) 
    MPICCLMYR_PCC=pearson(MPICCLMYR_mean.data,ObsY)
    MPIREMOYR_PCC=pearson(MPIREMOYR_mean.data,ObsY)
    MPISMHIYR_PCC=pearson(MPISMHIYR_mean.data,ObsY)  
    NCCDMIYR_PCC=pearson(NCCDMIYR_mean.data,ObsY) 
    NCCSMHIYR_PCC=pearson(NCCSMHIYR_mean.data,ObsY)  
    NOAAYR_PCC=pearson(NOAAYR_mean.data,ObsY)  
    
    AverageSONGY_PCC=pearson(AverageSONGY,ObsSONY)
    CanESM2SON_PCC=pearson(CanESM2SON_mean.data,ObsSONY)
    CNRMGSON_PCC=pearson(CNRMGSON_mean.data,ObsSONY) 
    MK3SON_PCC=pearson(MK3SON_mean.data,ObsSONY) 
    EARTHSON_PCC=pearson(EARTHSON_mean.data,ObsSONY) 
    EARTH3SON_PCC=pearson(EARTH3SON_mean.data,ObsSONY)
    GFDLSON_PCC=pearson(GFDLSON_mean.data,ObsSONY)  
    HadGEM2SON_PCC=pearson(HadGEM2SON_mean.data,ObsSONY) 
    IPSLGSON_PCC=pearson(IPSLGSON_mean.data,ObsSONY)
    MIROCGSON_PCC=pearson(MIROCGSON_mean.data,ObsSONY)  
    MPISON_PCC=pearson(MPISON_mean.data,ObsSONY)  
    NorESM1SON_PCC=pearson(NorESM1SON_mean.data,ObsSONY)
    
    AverageDJFGY_PCC=pearson(AverageDJFGY,ObsDJFY)
    CanESM2DJF_PCC=pearson(CanESM2DJF_mean.data,ObsDJFY)
    CNRMGDJF_PCC=pearson(CNRMGDJF_mean.data,ObsDJFY) 
    MK3DJF_PCC=pearson(MK3DJF_mean.data,ObsDJFY) 
    EARTHDJF_PCC=pearson(EARTHDJF_mean.data,ObsDJFY) 
    EARTH3DJF_PCC=pearson(EARTH3DJF_mean.data,ObsDJFY)
    GFDLDJF_PCC=pearson(GFDLDJF_mean.data,ObsDJFY)  
    HadGEM2DJF_PCC=pearson(HadGEM2DJF_mean.data,ObsDJFY) 
    IPSLGDJF_PCC=pearson(IPSLGDJF_mean.data,ObsDJFY)
    MIROCGDJF_PCC=pearson(MIROCGDJF_mean.data,ObsDJFY)  
    MPIDJF_PCC=pearson(MPIDJF_mean.data,ObsDJFY)  
    NorESM1DJF_PCC=pearson(NorESM1DJF_mean.data,ObsDJFY)
    
    AverageMAMGY_PCC=pearson(AverageMAMGY,ObsMAMY)
    CanESM2MAM_PCC=pearson(CanESM2MAM_mean.data,ObsMAMY)
    CNRMGMAM_PCC=pearson(CNRMGMAM_mean.data,ObsMAMY) 
    MK3MAM_PCC=pearson(MK3MAM_mean.data,ObsMAMY) 
    EARTHMAM_PCC=pearson(EARTHMAM_mean.data,ObsMAMY) 
    EARTH3MAM_PCC=pearson(EARTH3MAM_mean.data,ObsMAMY)
    GFDLMAM_PCC=pearson(GFDLMAM_mean.data,ObsMAMY)  
    HadGEM2MAM_PCC=pearson(HadGEM2MAM_mean.data,ObsMAMY) 
    IPSLGMAM_PCC=pearson(IPSLGMAM_mean.data,ObsMAMY)
    MIROCGMAM_PCC=pearson(MIROCGMAM_mean.data,ObsMAMY)  
    MPIMAM_PCC=pearson(MPIMAM_mean.data,ObsMAMY)  
    NorESM1MAM_PCC=pearson(NorESM1MAM_mean.data,ObsMAMY)
    
    AverageJJAGY_PCC=pearson(AverageJJAGY,ObsJJAY)
    CanESM2JJA_PCC=pearson(CanESM2JJA_mean.data,ObsJJAY)
    CNRMGJJA_PCC=pearson(CNRMGJJA_mean.data,ObsJJAY) 
    MK3JJA_PCC=pearson(MK3JJA_mean.data,ObsJJAY) 
    EARTHJJA_PCC=pearson(EARTHJJA_mean.data,ObsJJAY) 
    EARTH3JJA_PCC=pearson(EARTH3JJA_mean.data,ObsJJAY)
    GFDLJJA_PCC=pearson(GFDLJJA_mean.data,ObsJJAY)  
    HadGEM2JJA_PCC=pearson(HadGEM2JJA_mean.data,ObsJJAY) 
    IPSLGJJA_PCC=pearson(IPSLGJJA_mean.data,ObsJJAY)
    MIROCGJJA_PCC=pearson(MIROCGJJA_mean.data,ObsJJAY)  
    MPIJJA_PCC=pearson(MPIJJA_mean.data,ObsJJAY)  
    NorESM1JJA_PCC=pearson(NorESM1JJA_mean.data,ObsJJAY)
    
    AverageGY_PCC=pearson(AverageGY,ObsY)
    CanESM2YR_PCC=pearson(CanESM2YR_mean.data,ObsY)
    CNRMGYR_PCC=pearson(CNRMGYR_mean.data,ObsY)
    MK3YR_PCC=pearson(MK3YR_mean.data,ObsY)
    EARTHYR_PCC=pearson(EARTHYR_mean.data,ObsY)
    EARTH3YR_PCC=pearson(EARTH3YR_mean.data,ObsY)
    GFDLYR_PCC=pearson(GFDLYR_mean.data,ObsY)
    HadGEM2YR_PCC=pearson(HadGEM2YR_mean.data,ObsY)
    IPSLGYR_PCC=pearson(IPSLGYR_mean.data,ObsY)
    MIROCGYR_PCC=pearson(MIROCGYR_mean.data,ObsY)
    MPIYR_PCC=pearson(MPIYR_mean.data,ObsY)
    NorESM1YR_PCC=pearson(NorESM1YR_mean.data,ObsY)
    
    #-------------------------------------------------------------------------
    #PART 6: Plot Taylor Diagram 

    # Plot Observed data
    stdref = ObsEY_std
    # Plot data: ("Name of model", Standard deviation, Pearsons Correlation Coefficient)
    ERAINT = [
        ('CanRCM4_r2_ERAINT', CCCmaYR_std, CCCmaYR_PCC),
        ('CCLM4-8-17_v1_ERAINT', CLMcomYR_std, CLMcomYR_PCC),
        ('HIRHAM5_v2_ERAINT', DMIYR_std, DMIYR_PCC),
        ('RACMO22T_v1_ERAINT', KNMIYR_std, KNMIYR_PCC),
        ('REMO2009_v1_ERAINT', MPIEYR_std, MPIEYR_PCC),
        ('RCA4_v1_ERAINT', SMHIYR_std, SMHIYR_PCC)
        ]
    Ave = [('Average_ERAINT', AverageEY_std, AverageEY_PCC)]
    fig = PLT.figure()
    dia = TaylorDiagram(stdref, fig=fig, label='Reference')
    dia.samplePoints[0].set_color('r')  # Mark reference point as a red star
    # Add models to Taylor diagram
    for i, (name, stddev, corrcoef) in enumerate(ERAINT):
        dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='', mfc='k', mec='k', label=name) 
    for i, (name, stddev, corrcoef) in enumerate(Ave):
        dia.add_sample(stddev, corrcoef,marker="o", ms=10, ls='', mfc='k', mec='k', label=name)                
    #Add RMS contours, and label them
    contours = dia.add_contours(levels=5, colors='0.5')  # 5 levels in grey
    PLT.clabel(contours, inline=1, fontsize=10, fmt='%.1f')
    # Add a figure legend and title
    fig.legend(dia.samplePoints,
               [ p.get_label() for p in dia.samplePoints ], numpoints=1, prop=dict(size='small'), bbox_to_anchor=(0.,1.02,1.,.102), loc=8, ncol=3, mode='expand', borderaxespad=0.)
    fig.suptitle("TasMAX for Malawi 1990-2008", size='x-large')
    PLT.savefig('TasMAX_Taylor_ERAINT_Annual', bbox_inches='tight')
    PLT.show()
    
    # Plot Observed data
    stdref = ObsESONY_std
    # Plot data: ("Name of model", Standard deviation, Pearsons Correlation Coefficient)
    ERAINT = [
        ('CanRCM4_r2_ERAINT', CCCmaSON_std, CCCmaSON_PCC),
        ('CCLM4-8-17_v1_ERAINT', CLMcomSON_std, CLMcomSON_PCC),
        ('HIRHAM5_v2_ERAINT', DMISON_std, DMISON_PCC),
        ('RACMO22T_v1_ERAINT', KNMISON_std, KNMISON_PCC),
        ('REMO2009_v1_ERAINT', MPIESON_std, MPIESON_PCC),
        ('RCA4_v1_ERAINT', SMHISON_std, SMHISON_PCC),
        ]
    Ave = [('Average_ERAINT', AverageSONEY_std, AverageSONEY_PCC)]
    fig = PLT.figure()
    dia = TaylorDiagram(stdref, fig=fig, label='Reference')
    dia.samplePoints[0].set_color('r')  # Mark reference point as a red star
    # Add models to Taylor diagram
    for i, (name, stddev, corrcoef) in enumerate(ERAINT):
        dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='', mfc='k', mec='k', label=name)                   
    for i, (name, stddev, corrcoef) in enumerate(Ave):
        dia.add_sample(stddev, corrcoef,marker="o", ms=10, ls='', mfc='k', mec='k', label=name)
    #Add RMS contours, and label them
    contours = dia.add_contours(levels=5, colors='0.5')  # 5 levels in grey
    PLT.clabel(contours, inline=1, fontsize=10, fmt='%.1f')
    # Add a figure legend and title
    fig.legend(dia.samplePoints,
               [ p.get_label() for p in dia.samplePoints ], numpoints=1, prop=dict(size='small'), bbox_to_anchor=(0.,1.02,1.,.102), loc=8, ncol=3, mode='expand', borderaxespad=0.)
    fig.suptitle("TasMAX for Malawi SON 1990-2008", size='x-large')
    PLT.savefig('TasMAX_Taylor_ERAINT_SON', bbox_inches='tight')
    PLT.show()
    
    # Plot Observed data
    stdref = ObsEDJFY_std
    # Plot data: ("Name of model", Standard deviation, Pearsons Correlation Coefficient)
    ERAINT = [
        ('CanRCM4_r2_ERAINT', CCCmaDJF_std, CCCmaDJF_PCC),
        ('CCLM4-8-17_v1_ERAINT', CLMcomDJF_std, CLMcomDJF_PCC),
        ('HIRHAM5_v2_ERAINT', DMIDJF_std, DMIDJF_PCC),
        ('RACMO22T_v1_ERAINT', KNMIDJF_std, KNMIDJF_PCC),
        ('REMO2009_v1_ERAINT', MPIEDJF_std, MPIEDJF_PCC),
        ('RCA4_v1_ERAINT', SMHIDJF_std, SMHIDJF_PCC),
        ]
    Ave = [('Average_ERAINT', AverageDJFEY_std, AverageDJFEY_PCC)]
    fig = PLT.figure()
    dia = TaylorDiagram(stdref, fig=fig, label='Reference')
    dia.samplePoints[0].set_color('r')  # Mark reference point as a red star
    # Add models to Taylor diagram
    for i, (name, stddev, corrcoef) in enumerate(ERAINT):
        dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='', mfc='k', mec='k', label=name)                   
    for i, (name, stddev, corrcoef) in enumerate(Ave):
        dia.add_sample(stddev, corrcoef,marker="o", ms=10, ls='', mfc='k', mec='k', label=name)
    #Add RMS contours, and label them
    contours = dia.add_contours(levels=5, colors='0.5')  # 5 levels in grey
    PLT.clabel(contours, inline=1, fontsize=10, fmt='%.1f')
    # Add a figure legend and title
    fig.legend(dia.samplePoints,
               [ p.get_label() for p in dia.samplePoints ], numpoints=1, prop=dict(size='small'), bbox_to_anchor=(0.,1.02,1.,.102), loc=8, ncol=3, mode='expand', borderaxespad=0.)
    fig.suptitle("TasMAX for Malawi DJF 1990-2008", size='x-large')
    PLT.savefig('TasMAX_Taylor_ERAINT_DJF', bbox_inches='tight')
    PLT.show()
    
    # Plot Observed data
    stdref = ObsEMAMY_std
    # Plot data: ("Name of model", Standard deviation, Pearsons Correlation Coefficient)
    ERAINT = [
        ('CanRCM4_r2_ERAINT', CCCmaMAM_std, CCCmaMAM_PCC),
        ('CCLM4-8-17_v1_ERAINT', CLMcomMAM_std, CLMcomMAM_PCC),
        ('HIRHAM5_v2_ERAINT', DMIMAM_std, DMIMAM_PCC),
        ('RACMO22T_v1_ERAINT', KNMIMAM_std, KNMIMAM_PCC),
        ('REMO2009_v1_ERAINT', MPIEMAM_std, MPIEMAM_PCC),
        ('RCA4_v1_ERAINT', SMHIMAM_std, SMHIMAM_PCC),
        ]
    Ave = [('Average_ERAINT', AverageMAMEY_std, AverageMAMEY_PCC)]
    fig = PLT.figure()
    dia = TaylorDiagram(stdref, fig=fig, label='Reference')
    dia.samplePoints[0].set_color('r')  # Mark reference point as a red star
    # Add models to Taylor diagram
    for i, (name, stddev, corrcoef) in enumerate(ERAINT):
        dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='', mfc='k', mec='k', label=name)                   
    for i, (name, stddev, corrcoef) in enumerate(Ave):
        dia.add_sample(stddev, corrcoef,marker="o", ms=10, ls='', mfc='k', mec='k', label=name)
    #Add RMS contours, and label them
    contours = dia.add_contours(levels=5, colors='0.5')  # 5 levels in grey
    PLT.clabel(contours, inline=1, fontsize=10, fmt='%.1f')
    # Add a figure legend and title
    fig.legend(dia.samplePoints,
               [ p.get_label() for p in dia.samplePoints ], numpoints=1, prop=dict(size='small'), bbox_to_anchor=(0.,1.02,1.,.102), loc=8, ncol=3, mode='expand', borderaxespad=0.)
    fig.suptitle("TasMAX for Malawi MAM 1990-2008", size='x-large')
    PLT.savefig('TasMAX_Taylor_ERAINT_MAM', bbox_inches='tight')
    PLT.show()
    
    # Plot Observed data
    stdref = ObsEJJAY_std
    # Plot data: ("Name of model", Standard deviation, Pearsons Correlation Coefficient)
    ERAINT = [
        ('CanRCM4_r2_ERAINT', CCCmaJJA_std, CCCmaJJA_PCC),
        ('CCLM4-8-17_v1_ERAINT', CLMcomJJA_std, CLMcomJJA_PCC),
        ('HIRHAM5_v2_ERAINT', DMIJJA_std, DMIJJA_PCC),
        ('RACMO22T_v1_ERAINT', KNMIJJA_std, KNMIJJA_PCC),
        ('REMO2009_v1_ERAINT', MPIEJJA_std, MPIEJJA_PCC),
        ('RCA4_v1_ERAINT', SMHIJJA_std, SMHIJJA_PCC),
        ]
    Ave = [('Average_ERAINT', AverageJJAEY_std, AverageJJAEY_PCC)]
    fig = PLT.figure()
    dia = TaylorDiagram(stdref, fig=fig, label='Reference')
    dia.samplePoints[0].set_color('r')  # Mark reference point as a red star
    # Add models to Taylor diagram-
    for i, (name, stddev, corrcoef) in enumerate(ERAINT):
        dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='', mfc='k', mec='k', label=name)                   
    for i, (name, stddev, corrcoef) in enumerate(Ave):
        dia.add_sample(stddev, corrcoef,marker="o", ms=10, ls='', mfc='k', mec='k', label=name)
    #Add RMS contours, and label them
    contours = dia.add_contours(levels=5, colors='0.5')  # 5 levels in grey
    PLT.clabel(contours, inline=1, fontsize=10, fmt='%.1f')
    # Add a figure legend and title
    fig.legend(dia.samplePoints,
               [ p.get_label() for p in dia.samplePoints ], numpoints=1, prop=dict(size='small'), bbox_to_anchor=(0.,1.02,1.,.102), loc=8, ncol=3, mode='expand', borderaxespad=0.)
    fig.suptitle("TasMAX for Malawi JJA 1990-2008", size='x-large')
    PLT.savefig('TasMAX_Taylor_ERAINT_JJA', bbox_inches='tight')
    PLT.show()

    # Plot Observed data
    stdref = ObsY_std
    # Plot data: ("Name of model", Standard deviation, Pearsons Correlation Coefficient)
    RCM = [
        ('CanRCM4_CanESM2', CCCmaCanRCMYR_std, CCCmaCanRCMYR_PCC),
        ('RCA4_CanESM2', CCCmaSMHIYR_std, CCCmaSMHIYR_PCC),
        ('CCLM4-8-17_CNRM-CM5', CNRMYR_std, CNRMYR_PCC),
        ('RCA4_CNRM-CM5', CNRMSMHIYR_std, CNRMSMHIYR_PCC),
        ('RC4_CSIRO-MK3', CSIROYR_std, CSIROYR_PCC),
        ('HIRHAM5_EC-EARTH', ICHECDMIYR_std, ICHECDMIYR_PCC),
        ('CCLM4-8-17_EC-EARTH', ICHECCCLMYR_std, ICHECCCLMYR_PCC),
        ('RACMO22T_EC-EARTH', ICHECKNMIYR_std, ICHECKNMIYR_PCC),
        ('REMO2009_EC-EARTH', ICHECMPIYR_std, ICHECMPIYR_PCC),
        ('RCA4_EC-EARTH', ICHECSMHIYR_std, ICHECSMHIYR_PCC),
        ('RCA4_IPSL-CM5A-MR', IPSLYR_std, IPSLYR_PCC),
        ('RCA4_MIROC5', MIROCYR_std, MIROCYR_PCC),
        ('CCLM4-8-17_HadGEM2-ES', MOHCCCLMYR_std, MOHCCCLMYR_PCC),
        ('RACMO22T_HadGEM2-ES', MOHCKNMIYR_std, MOHCKNMIYR_PCC),
        ('RCA4_HadGEM2-ES', MOHCSMHIYR_std, MOHCSMHIYR_PCC),
        ('CCLM4-8-17_MPI-ESM-LR', MPICCLMYR_std, MPICCLMYR_PCC),
        ('REMO2009_MPI-ESM-LR', MPIREMOYR_std, MPIREMOYR_PCC),
        ('RCA4_MPI-ESM-LR', MPISMHIYR_std, MPISMHIYR_PCC),
        ('HIRHAM5_NorESM1-M', NCCDMIYR_std, NCCDMIYR_PCC),
        ('RCA4_NorESM1-M', NCCSMHIYR_std, NCCSMHIYR_PCC),
        ('RCA4_GFDL-ESM2M', NOAAYR_std, NOAAYR_PCC)
        ]
    RCM_Ave = [('Average_RCM', AverageRY_std, AverageRY_PCC)]
    GCM = [   
        ('CanESM2', CanESM2YR_std, CanESM2YR_PCC),
        ('CNRM_CM5', CNRMGYR_std, CNRMGYR_PCC),
        ('CSIRO MK3-6-0', MK3YR_std, MK3YR_PCC),
        ('EC-EARTH (r12i1p1)', EARTHYR_std, EARTHYR_PCC),
        ('EC-EARTH3 (r3i1p1)', EARTH3YR_std, EARTH3YR_PCC),
        ('GFDL-ESM2M', GFDLYR_std, GFDLYR_PCC),
        ('HadGEM2-ES', HadGEM2YR_std, HadGEM2YR_PCC),
        ('IPSL-CM5A-MR', IPSLGYR_std, IPSLGYR_PCC),
        ('MIROC5', MIROCGYR_std, MIROCGYR_PCC),
        ('MPI-ESM-LR', MPIYR_std, MPIYR_PCC),
        ('NorESM1-M', NorESM1YR_std, NorESM1YR_PCC)
        ] 
    GCM_Ave = [('Average_GCM', AverageGY_std, AverageGY_PCC)]
    fig = PLT.figure()
    dia = TaylorDiagram(stdref, fig=fig, label='Reference')
    dia.samplePoints[0].set_color('r')  # Mark reference point as a red star
    # Add models to Taylor diagram
    for i, (name, stddev, corrcoef) in enumerate(RCM):
        dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='', mfc='cyan', mec='cyan', label=name) 
    for i, (name, stddev, corrcoef) in enumerate(RCM_Ave):
        dia.add_sample(stddev, corrcoef,marker="o", ms=10, ls='', mfc='cyan', mec='cyan', label=name) 
    for i, (name, stddev, corrcoef) in enumerate(GCM):
        dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='', mfc='magenta', mec='magenta', label=name)                
    for i, (name, stddev, corrcoef) in enumerate(GCM_Ave):
        dia.add_sample(stddev, corrcoef,marker="o", ms=10, ls='', mfc='magenta', mec='magenta', label=name)                   
    #Add RMS contours, and label them
    contours = dia.add_contours(levels=5, colors='0.5')  # 5 levels in grey
    PLT.clabel(contours, inline=1, fontsize=10, fmt='%.1f')
    # Add a figure legend and title
    fig.legend(dia.samplePoints,
               [ p.get_label() for p in dia.samplePoints ], numpoints=1, prop=dict(size='small'), bbox_to_anchor=(0.,1.02,1.,.102), loc=8, ncol=3, mode='expand', borderaxespad=0.)
    fig.suptitle("TasMAX for Malawi 1961-2005", size='x-large')
    PLT.savefig('TasMAX_Taylor_Annual', bbox_inches='tight')
    PLT.show()
    
    # Plot Observed data
    stdref = ObsSONY_std
    # Plot data: ("Name of model", Standard deviation, Pearsons Correlation Coefficient)
    RCM = [
        ('CanRCM4_CanESM2', CCCmaCanRCMSON_std, CCCmaCanRCMSON_PCC),
        ('RCA4_CanESM2', CCCmaSMHISON_std, CCCmaSMHISON_PCC),
        ('CCLM4-8-17_CNRM-CM5', CNRMSON_std, CNRMSON_PCC),
        ('RCA4_CNRM-CM5', CNRMSMHISON_std, CNRMSMHISON_PCC),
        ('RC4_CSIRO-MK3', CSIROSON_std, CSIROSON_PCC),
        ('HIRHAM5_EC-EARTH', ICHECDMISON_std, ICHECDMISON_PCC),
        ('CCLM4-8-17_EC-EARTH', ICHECCCLMSON_std, ICHECCCLMSON_PCC),
        ('RACMO22T_EC-EARTH', ICHECKNMISON_std, ICHECKNMISON_PCC),
        ('REMO2009_EC-EARTH', ICHECMPISON_std, ICHECMPISON_PCC),
        ('RCA4_EC-EARTH', ICHECSMHISON_std, ICHECSMHISON_PCC),
        ('RCA4_IPSL-CM5A-MR', IPSLSON_std, IPSLSON_PCC),
        ('RCA4_MIROC5', MIROCSON_std, MIROCSON_PCC),
        ('CCLM4-8-17_HadGEM2-ES', MOHCCCLMSON_std, MOHCCCLMSON_PCC),
        ('RACMO22T_HadGEM2-ES', MOHCKNMISON_std, MOHCKNMISON_PCC),
        ('RCA4_HadGEM2-ES', MOHCSMHISON_std, MOHCSMHISON_PCC),
        ('CCLM4-8-17_MPI-ESM-LR', MPICCLMSON_std, MPICCLMSON_PCC),
        ('REMO2009_MPI-ESM-LR', MPIREMOSON_std, MPIREMOSON_PCC),
        ('RCA4_MPI-ESM-LR', MPISMHISON_std, MPISMHISON_PCC),
        ('HIRHAM5_NorESM1-M', NCCDMISON_std, NCCDMISON_PCC),
        ('RCA4_NorESM1-M', NCCSMHISON_std, NCCSMHISON_PCC),
        ('RCA4_GFDL-ESM2M', NOAASON_std, NOAASON_PCC)
        ]
    RCM_Ave = [('Average_RCM', AverageSONRY_std, AverageSONRY_PCC)]
    GCM = [   
        ('CanESM2', CanESM2SON_std, CanESM2SON_PCC),
        ('CNRM_CM5', CNRMGSON_std, CNRMGSON_PCC),
        ('CSIRO MK3-6-0', MK3SON_std, MK3SON_PCC),
        ('EC-EARTH (r12i1p1)', EARTHSON_std, EARTHSON_PCC),
        ('EC-EARTH3 (r3i1p1)', EARTH3SON_std, EARTH3SON_PCC),
        ('GFDL-ESM2M', GFDLSON_std, GFDLSON_PCC),
        ('HadGEM2-ES', HadGEM2SON_std, HadGEM2SON_PCC),
        ('IPSL-CM5A-MR', IPSLGSON_std, IPSLGSON_PCC),
        ('MIROC5', MIROCGSON_std, MIROCGSON_PCC),
        ('MPI-ESM-LR', MPISON_std, MPISON_PCC),
        ('NorESM1-M', NorESM1SON_std, NorESM1SON_PCC)
        ] 
    GCM_Ave = [('Average_GCM', AverageSONGY_std, AverageSONGY_PCC)]
    fig = PLT.figure()
    dia = TaylorDiagram(stdref, fig=fig, label='Reference')
    dia.samplePoints[0].set_color('r')  # Mark reference point as a red star
    # Add models to Taylor diagram
    for i, (name, stddev, corrcoef) in enumerate(RCM):
        dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='', mfc='cyan', mec='cyan', label=name) 
    for i, (name, stddev, corrcoef) in enumerate(RCM_Ave):
        dia.add_sample(stddev, corrcoef,marker="o", ms=10, ls='', mfc='cyan', mec='cyan', label=name) 
    for i, (name, stddev, corrcoef) in enumerate(GCM):
        dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='', mfc='magenta', mec='magenta', label=name)                
    for i, (name, stddev, corrcoef) in enumerate(GCM_Ave):
        dia.add_sample(stddev, corrcoef,marker="o", ms=10, ls='', mfc='magenta', mec='magenta', label=name)                   
    #Add RMS contours, and label them
    contours = dia.add_contours(levels=5, colors='0.5')  # 5 levels in grey
    PLT.clabel(contours, inline=1, fontsize=10, fmt='%.1f')
    # Add a figure legend and title
    fig.legend(dia.samplePoints,
               [ p.get_label() for p in dia.samplePoints ], numpoints=1, prop=dict(size='small'), bbox_to_anchor=(0.,1.02,1.,.102), loc=8, ncol=3, mode='expand', borderaxespad=0.)
    fig.suptitle("TasMAX for Malawi SON 1961-2005", size='x-large')
    PLT.savefig('TasMAX_Taylor_SON', bbox_inches='tight')
    PLT.show()
    
    # Plot Observed data
    stdref = ObsDJFY_std
    # Plot data: ("Name of model", Standard deviation, Pearsons Correlation Coefficient)
    RCM = [
        ('CanRCM4_CanESM2', CCCmaCanRCMDJF_std, CCCmaCanRCMDJF_PCC),
        ('RCA4_CanESM2', CCCmaSMHIDJF_std, CCCmaSMHIDJF_PCC),
        ('CCLM4-8-17_CNRM-CM5', CNRMDJF_std, CNRMDJF_PCC),
        ('RCA4_CNRM-CM5', CNRMSMHIDJF_std, CNRMSMHIDJF_PCC),
        ('RC4_CSIRO-MK3', CSIRODJF_std, CSIRODJF_PCC),
        ('HIRHAM5_EC-EARTH', ICHECDMIDJF_std, ICHECDMIDJF_PCC),
        ('CCLM4-8-17_EC-EARTH', ICHECCCLMDJF_std, ICHECCCLMDJF_PCC),
        ('RACMO22T_EC-EARTH', ICHECKNMIDJF_std, ICHECKNMIDJF_PCC),
        ('REMO2009_EC-EARTH', ICHECMPIDJF_std, ICHECMPIDJF_PCC),
        ('RCA4_EC-EARTH', ICHECSMHIDJF_std, ICHECSMHIDJF_PCC),
        ('RCA4_IPSL-CM5A-MR', IPSLDJF_std, IPSLDJF_PCC),
        ('RCA4_MIROC5', MIROCDJF_std, MIROCDJF_PCC),
        ('CCLM4-8-17_HadGEM2-ES', MOHCCCLMDJF_std, MOHCCCLMDJF_PCC),
        ('RACMO22T_HadGEM2-ES', MOHCKNMIDJF_std, MOHCKNMIDJF_PCC),
        ('RCA4_HadGEM2-ES', MOHCSMHIDJF_std, MOHCSMHIDJF_PCC),
        ('CCLM4-8-17_MPI-ESM-LR', MPICCLMDJF_std, MPICCLMDJF_PCC),
        ('REMO2009_MPI-ESM-LR', MPIREMODJF_std, MPIREMODJF_PCC),
        ('RCA4_MPI-ESM-LR', MPISMHIDJF_std, MPISMHIDJF_PCC),
        ('HIRHAM5_NorESM1-M', NCCDMIDJF_std, NCCDMIDJF_PCC),
        ('RCA4_NorESM1-M', NCCSMHIDJF_std, NCCSMHIDJF_PCC),
        ('RCA4_GFDL-ESM2M', NOAADJF_std, NOAADJF_PCC)
        ]
    RCM_Ave = [('Average_RCM', AverageDJFRY_std, AverageDJFRY_PCC)]
    GCM = [   
        ('CanESM2', CanESM2DJF_std, CanESM2DJF_PCC),
        ('CNRM_CM5', CNRMGDJF_std, CNRMGDJF_PCC),
        ('CSIRO MK3-6-0', MK3DJF_std, MK3DJF_PCC),
        ('EC-EARTH (r12i1p1)', EARTHDJF_std, EARTHDJF_PCC),
        ('EC-EARTH3 (r3i1p1)', EARTH3DJF_std, EARTH3DJF_PCC),
        ('GFDL-ESM2M', GFDLDJF_std, GFDLDJF_PCC),
        ('HadGEM2-ES', HadGEM2DJF_std, HadGEM2DJF_PCC),
        ('IPSL-CM5A-MR', IPSLGDJF_std, IPSLGDJF_PCC),
        ('MIROC5', MIROCGDJF_std, MIROCGDJF_PCC),
        ('MPI-ESM-LR', MPIDJF_std, MPIDJF_PCC),
        ('NorESM1-M', NorESM1DJF_std, NorESM1DJF_PCC)
        ] 
    GCM_Ave = [('Average_GCM', AverageDJFGY_std, AverageDJFGY_PCC)]
    fig = PLT.figure()
    dia = TaylorDiagram(stdref, fig=fig, label='Reference')
    dia.samplePoints[0].set_color('r')  # Mark reference point as a red star
    # Add models to Taylor diagram
    for i, (name, stddev, corrcoef) in enumerate(RCM):
        dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='', mfc='cyan', mec='cyan', label=name) 
    for i, (name, stddev, corrcoef) in enumerate(RCM_Ave):
        dia.add_sample(stddev, corrcoef,marker="o", ms=10, ls='', mfc='cyan', mec='cyan', label=name) 
    for i, (name, stddev, corrcoef) in enumerate(GCM):
        dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='', mfc='magenta', mec='magenta', label=name)                
    for i, (name, stddev, corrcoef) in enumerate(GCM_Ave):
        dia.add_sample(stddev, corrcoef,marker="o", ms=10, ls='', mfc='magenta', mec='magenta', label=name)                   
    #Add RMS contours, and label them
    contours = dia.add_contours(levels=5, colors='0.5')  # 5 levels in grey
    PLT.clabel(contours, inline=1, fontsize=10, fmt='%.1f')
    # Add a figure legend and title
    fig.legend(dia.samplePoints,
               [ p.get_label() for p in dia.samplePoints ], numpoints=1, prop=dict(size='small'), bbox_to_anchor=(0.,1.02,1.,.102), loc=8, ncol=3, mode='expand', borderaxespad=0.)
    fig.suptitle("TasMAX for Malawi DJF 1961-2005", size='x-large')
    PLT.savefig('TasMAX_Taylor_DJF', bbox_inches='tight')
    PLT.show()
    
    # Plot Observed data
    stdref = ObsMAMY_std
    # Plot data: ("Name of model", Standard deviation, Pearsons Correlation Coefficient)
    RCM = [
        ('CanRCM4_CanESM2', CCCmaCanRCMMAM_std, CCCmaCanRCMMAM_PCC),
        ('RCA4_CanESM2', CCCmaSMHIMAM_std, CCCmaSMHIMAM_PCC),
        ('CCLM4-8-17_CNRM-CM5', CNRMMAM_std, CNRMMAM_PCC),
        ('RCA4_CNRM-CM5', CNRMSMHIMAM_std, CNRMSMHIMAM_PCC),
        ('RC4_CSIRO-MK3', CSIROMAM_std, CSIROMAM_PCC),
        ('HIRHAM5_EC-EARTH', ICHECDMIMAM_std, ICHECDMIMAM_PCC),
        ('CCLM4-8-17_EC-EARTH', ICHECCCLMMAM_std, ICHECCCLMMAM_PCC),
        ('RACMO22T_EC-EARTH', ICHECKNMIMAM_std, ICHECKNMIMAM_PCC),
        ('REMO2009_EC-EARTH', ICHECMPIMAM_std, ICHECMPIMAM_PCC),
        ('RCA4_EC-EARTH', ICHECSMHIMAM_std, ICHECSMHIMAM_PCC),
        ('RCA4_IPSL-CM5A-MR', IPSLMAM_std, IPSLMAM_PCC),
        ('RCA4_MIROC5', MIROCMAM_std, MIROCMAM_PCC),
        ('CCLM4-8-17_HadGEM2-ES', MOHCCCLMMAM_std, MOHCCCLMMAM_PCC),
        ('RACMO22T_HadGEM2-ES', MOHCKNMIMAM_std, MOHCKNMIMAM_PCC),
        ('RCA4_HadGEM2-ES', MOHCSMHIMAM_std, MOHCSMHIMAM_PCC),
        ('CCLM4-8-17_MPI-ESM-LR', MPICCLMMAM_std, MPICCLMMAM_PCC),
        ('REMO2009_MPI-ESM-LR', MPIREMOMAM_std, MPIREMOMAM_PCC),
        ('RCA4_MPI-ESM-LR', MPISMHIMAM_std, MPISMHIMAM_PCC),
        ('HIRHAM5_NorESM1-M', NCCDMIMAM_std, NCCDMIMAM_PCC),
        ('RCA4_NorESM1-M', NCCSMHIMAM_std, NCCSMHIMAM_PCC),
        ('RCA4_GFDL-ESM2M', NOAAMAM_std, NOAAMAM_PCC)
        ]
    RCM_Ave = [('Average_RCM', AverageMAMRY_std, AverageMAMRY_PCC)]
    GCM = [   
        ('CanESM2', CanESM2MAM_std, CanESM2MAM_PCC),
        ('CNRM_CM5', CNRMGMAM_std, CNRMGMAM_PCC),
        ('CSIRO MK3-6-0', MK3MAM_std, MK3MAM_PCC),
        ('EC-EARTH (r12i1p1)', EARTHMAM_std, EARTHMAM_PCC),
        ('EC-EARTH3 (r3i1p1)', EARTH3MAM_std, EARTH3MAM_PCC),
        ('GFDL-ESM2M', GFDLMAM_std, GFDLMAM_PCC),
        ('HadGEM2-ES', HadGEM2MAM_std, HadGEM2MAM_PCC),
        ('IPSL-CM5A-MR', IPSLGMAM_std, IPSLGMAM_PCC),
        ('MIROC5', MIROCGMAM_std, MIROCGMAM_PCC),
        ('MPI-ESM-LR', MPIMAM_std, MPIMAM_PCC),
        ('NorESM1-M', NorESM1MAM_std, NorESM1MAM_PCC)
        ] 
    GCM_Ave = [('Average_GCM', AverageMAMGY_std, AverageMAMGY_PCC)]
    fig = PLT.figure()
    dia = TaylorDiagram(stdref, fig=fig, label='Reference')
    dia.samplePoints[0].set_color('r')  # Mark reference point as a red star
    # Add models to Taylor diagram
    for i, (name, stddev, corrcoef) in enumerate(RCM):
        dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='', mfc='cyan', mec='cyan', label=name) 
    for i, (name, stddev, corrcoef) in enumerate(RCM_Ave):
        dia.add_sample(stddev, corrcoef,marker="o", ms=10, ls='', mfc='cyan', mec='cyan', label=name) 
    for i, (name, stddev, corrcoef) in enumerate(GCM):
        dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='', mfc='magenta', mec='magenta', label=name)                
    for i, (name, stddev, corrcoef) in enumerate(GCM_Ave):
        dia.add_sample(stddev, corrcoef,marker="o", ms=10, ls='', mfc='magenta', mec='magenta', label=name)                   
    #Add RMS contours, and label them
    contours = dia.add_contours(levels=5, colors='0.5')  # 5 levels in grey
    PLT.clabel(contours, inline=1, fontsize=10, fmt='%.1f')
    # Add a figure legend and title
    fig.legend(dia.samplePoints,
               [ p.get_label() for p in dia.samplePoints ], numpoints=1, prop=dict(size='small'), bbox_to_anchor=(0.,1.02,1.,.102), loc=8, ncol=3, mode='expand', borderaxespad=0.)
    fig.suptitle("TasMAX for Malawi MAM 1961-2005", size='x-large')
    PLT.savefig('TasMAX_Taylor_MAM', bbox_inches='tight')
    PLT.show()
    
    # Plot Observed data
    stdref = ObsJJAY_std
    # Plot data: ("Name of model", Standard deviation, Pearsons Correlation Coefficient)
    RCM = [
        ('CanRCM4_CanESM2', CCCmaCanRCMJJA_std, CCCmaCanRCMJJA_PCC),
        ('RCA4_CanESM2', CCCmaSMHIJJA_std, CCCmaSMHIJJA_PCC),
        ('CCLM4-8-17_CNRM-CM5', CNRMJJA_std, CNRMJJA_PCC),
        ('RCA4_CNRM-CM5', CNRMSMHIJJA_std, CNRMSMHIJJA_PCC),
        ('RC4_CSIRO-MK3', CSIROJJA_std, CSIROJJA_PCC),
        ('HIRHAM5_EC-EARTH', ICHECDMIJJA_std, ICHECDMIJJA_PCC),
        ('CCLM4-8-17_EC-EARTH', ICHECCCLMJJA_std, ICHECCCLMJJA_PCC),
        ('RACMO22T_EC-EARTH', ICHECKNMIJJA_std, ICHECKNMIJJA_PCC),
        ('REMO2009_EC-EARTH', ICHECMPIJJA_std, ICHECMPIJJA_PCC),
        ('RCA4_EC-EARTH', ICHECSMHIJJA_std, ICHECSMHIJJA_PCC),
        ('RCA4_IPSL-CM5A-MR', IPSLJJA_std, IPSLJJA_PCC),
        ('RCA4_MIROC5', MIROCJJA_std, MIROCJJA_PCC),
        ('CCLM4-8-17_HadGEM2-ES', MOHCCCLMJJA_std, MOHCCCLMJJA_PCC),
        ('RACMO22T_HadGEM2-ES', MOHCKNMIJJA_std, MOHCKNMIJJA_PCC),
        ('RCA4_HadGEM2-ES', MOHCSMHIJJA_std, MOHCSMHIJJA_PCC),
        ('CCLM4-8-17_MPI-ESM-LR', MPICCLMJJA_std, MPICCLMJJA_PCC),
        ('REMO2009_MPI-ESM-LR', MPIREMOJJA_std, MPIREMOJJA_PCC),
        ('RCA4_MPI-ESM-LR', MPISMHIJJA_std, MPISMHIJJA_PCC),
        ('HIRHAM5_NorESM1-M', NCCDMIJJA_std, NCCDMIJJA_PCC),
        ('RCA4_NorESM1-M', NCCSMHIJJA_std, NCCSMHIJJA_PCC),
        ('RCA4_GFDL-ESM2M', NOAAJJA_std, NOAAJJA_PCC)
        ]
    RCM_Ave = [('Average_RCM', AverageJJARY_std, AverageJJARY_PCC)]
    GCM = [   
        ('CanESM2', CanESM2JJA_std, CanESM2JJA_PCC),
        ('CNRM_CM5', CNRMGJJA_std, CNRMGJJA_PCC),
        ('CSIRO MK3-6-0', MK3JJA_std, MK3JJA_PCC),
        ('EC-EARTH (r12i1p1)', EARTHJJA_std, EARTHJJA_PCC),
        ('EC-EARTH3 (r3i1p1)', EARTH3JJA_std, EARTH3JJA_PCC),
        ('GFDL-ESM2M', GFDLJJA_std, GFDLJJA_PCC),
        ('HadGEM2-ES', HadGEM2JJA_std, HadGEM2JJA_PCC),
        ('IPSL-CM5A-MR', IPSLGJJA_std, IPSLGJJA_PCC),
        ('MIROC5', MIROCGJJA_std, MIROCGJJA_PCC),
        ('MPI-ESM-LR', MPIJJA_std, MPIJJA_PCC),
        ('NorESM1-M', NorESM1JJA_std, NorESM1JJA_PCC)
        ] 
    GCM_Ave = [('Average_GCM', AverageJJAGY_std, AverageJJAGY_PCC)]
    fig = PLT.figure()
    dia = TaylorDiagram(stdref, fig=fig, label='Reference')
    dia.samplePoints[0].set_color('r')  # Mark reference point as a red star
    # Add models to Taylor diagram
    for i, (name, stddev, corrcoef) in enumerate(RCM):
        dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='', mfc='cyan', mec='cyan', label=name) 
    for i, (name, stddev, corrcoef) in enumerate(RCM_Ave):
        dia.add_sample(stddev, corrcoef,marker="o", ms=10, ls='', mfc='cyan', mec='cyan', label=name) 
    for i, (name, stddev, corrcoef) in enumerate(GCM):
        dia.add_sample(stddev, corrcoef,marker='$%d$' % (i+1), ms=10, ls='', mfc='magenta', mec='magenta', label=name)                
    for i, (name, stddev, corrcoef) in enumerate(GCM_Ave):
        dia.add_sample(stddev, corrcoef,marker="o", ms=10, ls='', mfc='magenta', mec='magenta', label=name)                   
    #Add RMS contours, and label them
    contours = dia.add_contours(levels=5, colors='0.5')  # 5 levels in grey
    PLT.clabel(contours, inline=1, fontsize=10, fmt='%.1f')
    # Add a figure legend and title
    fig.legend(dia.samplePoints,
               [ p.get_label() for p in dia.samplePoints ], numpoints=1, prop=dict(size='small'), bbox_to_anchor=(0.,1.02,1.,.102), loc=8, ncol=3, mode='expand', borderaxespad=0.)
    fig.suptitle("TasMAX for Malawi JJA 1961-2005", size='x-large')
    PLT.savefig('TasMAX_Taylor_JJA', bbox_inches='tight')
    PLT.show()
    
if __name__ == '__main__':
    main()
