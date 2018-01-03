# -*- coding: utf-8 -*-
"""
Created on Wednesday August 2 2017

@author: s0899345
"""

import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import numpy as np
from cf_units import Unit
import iris
import iris.plot as iplt
import cartopy
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import iris.analysis.cartography
import iris.coord_categorisation as iriscc

#this file is split into parts as follows:
    #PART 1: load and format all models 
    #PART 2: load and format observed data
    #PART 3: format files to be plot specific
    #PART 4: plot data

def main():
    iris.FUTURE.netcdf_promote=True    
    
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
      
    #---------------------------------------------------------------------------------------------------------------------
     #PART 2: OBSERVED DATA
    #bring in all the files we need and give them a name
    CRU= '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/Actual_Data/cru_ts4.01.1901.2016.tmx.dat.nc'
        
    #Load exactly one cube from given file
    CRU = iris.load_cube(CRU, 'near-surface temperature maximum')
    
    #---------------------------------------------------------------------------------------------------------------------
    #PART 3: FORMAT DATA TO BE PLOT SPECIFIC 
    #regrid all models to have same latitude and longitude system, all regridded to model with lowest resolution
    #CCCma = CCCma.regrid(CCCma, iris.analysis.Linear())
    CLMcom =CLMcom.regrid(CCCma, iris.analysis.Linear())
    DMI=DMI.regrid(CCCma, iris.analysis.Linear())
    KNMI=KNMI.regrid(CCCma, iris.analysis.Linear())
    MPIE=MPIE.regrid(CCCma, iris.analysis.Linear())
    SMHI=SMHI.regrid(CCCma, iris.analysis.Linear())
    
    #CCCmaCanRCM = CCCmaCanRCM.regrid(CCCmaCanRCM, iris.analysis.Linear())
    CCCmaSMHI = CCCmaSMHI.regrid(CCCmaCanRCM, iris.analysis.Linear())
    CNRM =CNRM.regrid(CCCmaCanRCM, iris.analysis.Linear())
    CNRMSMHI =CNRMSMHI.regrid(CCCmaCanRCM, iris.analysis.Linear())
    CSIRO=CSIRO.regrid(CCCmaCanRCM, iris.analysis.Linear())
    ICHECDMI=ICHECDMI.regrid(CCCmaCanRCM, iris.analysis.Linear())
    ICHECCCLM=ICHECCCLM.regrid(CCCmaCanRCM, iris.analysis.Linear())
    ICHECKNMI=ICHECKNMI.regrid(CCCmaCanRCM, iris.analysis.Linear())
    ICHECMPI=ICHECMPI.regrid(CCCmaCanRCM, iris.analysis.Linear())
    ICHECSMHI=ICHECSMHI.regrid(CCCmaCanRCM, iris.analysis.Linear())
    IPSL=IPSL.regrid(CCCmaCanRCM, iris.analysis.Linear())
    MIROC=MIROC.regrid(CCCmaCanRCM, iris.analysis.Linear())
    MOHCCCLM=MOHCCCLM.regrid(CCCmaCanRCM, iris.analysis.Linear())
    MOHCKNMI=MOHCKNMI.regrid(CCCmaCanRCM, iris.analysis.Linear())
    MOHCSMHI=MOHCSMHI.regrid(CCCmaCanRCM, iris.analysis.Linear())
    MPICCLM=MPICCLM.regrid(CCCmaCanRCM, iris.analysis.Linear())
    MPIREMO=MPIREMO.regrid(CCCmaCanRCM, iris.analysis.Linear())
    MPISMHI=MPISMHI.regrid(CCCmaCanRCM, iris.analysis.Linear())
    NCCDMI=NCCDMI.regrid(CCCmaCanRCM, iris.analysis.Linear())
    NCCSMHI=NCCSMHI.regrid(CCCmaCanRCM, iris.analysis.Linear())
    NOAA=NOAA.regrid(CCCmaCanRCM, iris.analysis.Linear())
    
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
    
    CRU = CRU.regrid(CCCma, iris.analysis.Linear())
    CRUG = CRU.regrid(CanESM2, iris.analysis.Linear())
            
    #we are only interested in the latitude and longitude relevant to Malawi (has to be slightly larger than country boundary to take into account resolution of GCMs)
    Malawi = iris.Constraint(longitude=lambda v: 32.0 <= v <= 36.5, latitude=lambda v: -17.5 <= v <= -8.)   
    
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
    CRUG = CRUG.extract(Malawi)
                                                  
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
    CRUG = CRUG.extract(t_constraint)
        
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
    
    CRU.units = Unit('Celsius') # This fixes CRU which is in 'Degrees Celsius' to read 'Celsius'
    CRUE.units = Unit('Celsius') # This fixes CRU which is in 'Degrees Celsius' to read 'Celsius'
    CRUG.units = Unit('Celsius') # This fixes CRU which is in 'Degrees Celsius' to read 'Celsius'
        
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
    iriscc.add_season(CRUG, 'time')
        
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
    iriscc.add_year(CRUG, 'time')
        
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
    CRUGSON = CRUG.extract(SON)
            
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
    CRUGDJF = CRUG.extract(DJF)
          
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
    CRUGMAM = CRUG.extract(MAM)
              
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
    CRUGJJA = CRUG.extract(JJA)
   
    #We are interested in plotting the data for the average of the time period. 
    CCCmaSON = CCCmaSON.collapsed('year', iris.analysis.MEAN)
    CLMcomSON = CLMcomSON.collapsed('year', iris.analysis.MEAN)    
    DMISON = DMISON.collapsed('year', iris.analysis.MEAN)    
    KNMISON = KNMISON.collapsed('year', iris.analysis.MEAN)
    MPIESON = MPIESON.collapsed('year', iris.analysis.MEAN)    
    SMHISON = SMHISON.collapsed('year', iris.analysis.MEAN)    
    
    CCCmaCanRCMSON = CCCmaCanRCMSON.collapsed('year', iris.analysis.MEAN)
    CCCmaSMHISON = CCCmaSMHISON.collapsed('year', iris.analysis.MEAN)
    CNRMSON = CNRMSON.collapsed('year', iris.analysis.MEAN)
    CNRMSMHISON = CNRMSMHISON.collapsed('year', iris.analysis.MEAN)
    CSIROSON = CSIROSON.collapsed('year', iris.analysis.MEAN)
    ICHECDMISON = ICHECDMISON.collapsed('year', iris.analysis.MEAN)
    ICHECCCLMSON = ICHECCCLMSON.collapsed('year', iris.analysis.MEAN)
    ICHECKNMISON = ICHECKNMISON.collapsed('year', iris.analysis.MEAN)
    ICHECMPISON = ICHECMPISON.collapsed('year', iris.analysis.MEAN)
    ICHECSMHISON = ICHECSMHISON.collapsed('year', iris.analysis.MEAN)
    IPSLSON = IPSLSON.collapsed('year', iris.analysis.MEAN)
    MIROCSON = MIROCSON.collapsed('year', iris.analysis.MEAN)
    MOHCCCLMSON = MOHCCCLMSON.collapsed('year', iris.analysis.MEAN)
    MOHCKNMISON = MOHCKNMISON.collapsed('year', iris.analysis.MEAN)
    MOHCSMHISON = MOHCSMHISON.collapsed('year', iris.analysis.MEAN)
    MPICCLMSON = MPICCLMSON.collapsed('year', iris.analysis.MEAN)
    MPIREMOSON = MPIREMOSON.collapsed('year', iris.analysis.MEAN)
    MPISMHISON = MPISMHISON.collapsed('year', iris.analysis.MEAN)
    NCCDMISON = NCCDMISON.collapsed('year', iris.analysis.MEAN)
    NCCSMHISON = NCCSMHISON.collapsed('year', iris.analysis.MEAN)
    NOAASON = NOAASON.collapsed('year', iris.analysis.MEAN)
    
    CanESM2SON = CanESM2SON.collapsed('year', iris.analysis.MEAN)
    CNRMGSON = CNRMGSON.collapsed('year', iris.analysis.MEAN)
    MK3SON = MK3SON.collapsed('year', iris.analysis.MEAN)
    EARTHSON = EARTHSON.collapsed('year', iris.analysis.MEAN)
    EARTH3SON = EARTH3SON.collapsed('year', iris.analysis.MEAN)                       
    GFDLSON = GFDLSON.collapsed('year', iris.analysis.MEAN)
    HadGEM2SON = HadGEM2SON.collapsed('year', iris.analysis.MEAN)
    IPSLGSON = IPSLGSON.collapsed('year', iris.analysis.MEAN)
    MIROCGSON = MIROCGSON.collapsed('year', iris.analysis.MEAN)
    MPISON = MPISON.collapsed('year', iris.analysis.MEAN)
    NorESM1SON = NorESM1SON.collapsed('year', iris.analysis.MEAN)   
    
    CRUESON = CRUESON.collapsed('year', iris.analysis.MEAN)  
    CRUSON = CRUSON.collapsed('year', iris.analysis.MEAN)  
    CRUGSON = CRUGSON.collapsed('year', iris.analysis.MEAN)  
           
    CCCmaDJF = CCCmaDJF.collapsed('year', iris.analysis.MEAN)
    CLMcomDJF = CLMcomDJF.collapsed('year', iris.analysis.MEAN)    
    DMIDJF = DMIDJF.collapsed('year', iris.analysis.MEAN)    
    KNMIDJF = KNMIDJF.collapsed('year', iris.analysis.MEAN)
    MPIEDJF = MPIEDJF.collapsed('year', iris.analysis.MEAN)    
    SMHIDJF = SMHIDJF.collapsed('year', iris.analysis.MEAN)    
    
    CCCmaCanRCMDJF = CCCmaCanRCMDJF.collapsed('year', iris.analysis.MEAN)
    CCCmaSMHIDJF = CCCmaSMHIDJF.collapsed('year', iris.analysis.MEAN)
    CNRMDJF = CNRMDJF.collapsed('year', iris.analysis.MEAN)
    CNRMSMHIDJF = CNRMSMHIDJF.collapsed('year', iris.analysis.MEAN)
    CSIRODJF = CSIRODJF.collapsed('year', iris.analysis.MEAN)
    ICHECDMIDJF = ICHECDMIDJF.collapsed('year', iris.analysis.MEAN)
    ICHECCCLMDJF = ICHECCCLMDJF.collapsed('year', iris.analysis.MEAN)
    ICHECKNMIDJF = ICHECKNMIDJF.collapsed('year', iris.analysis.MEAN)
    ICHECMPIDJF = ICHECMPIDJF.collapsed('year', iris.analysis.MEAN)
    ICHECSMHIDJF = ICHECSMHIDJF.collapsed('year', iris.analysis.MEAN)
    IPSLDJF = IPSLDJF.collapsed('year', iris.analysis.MEAN)
    MIROCDJF = MIROCDJF.collapsed('year', iris.analysis.MEAN)
    MOHCCCLMDJF = MOHCCCLMDJF.collapsed('year', iris.analysis.MEAN)
    MOHCKNMIDJF = MOHCKNMIDJF.collapsed('year', iris.analysis.MEAN)
    MOHCSMHIDJF = MOHCSMHIDJF.collapsed('year', iris.analysis.MEAN)
    MPICCLMDJF = MPICCLMDJF.collapsed('year', iris.analysis.MEAN)
    MPIREMODJF = MPIREMODJF.collapsed('year', iris.analysis.MEAN)
    MPISMHIDJF = MPISMHIDJF.collapsed('year', iris.analysis.MEAN)
    NCCDMIDJF = NCCDMIDJF.collapsed('year', iris.analysis.MEAN)
    NCCSMHIDJF = NCCSMHIDJF.collapsed('year', iris.analysis.MEAN)
    NOAADJF = NOAADJF.collapsed('year', iris.analysis.MEAN)
    
    CanESM2DJF = CanESM2DJF.collapsed('year', iris.analysis.MEAN)
    CNRMGDJF = CNRMGDJF.collapsed('year', iris.analysis.MEAN)
    MK3DJF = MK3DJF.collapsed('year', iris.analysis.MEAN)
    EARTHDJF = EARTHDJF.collapsed('year', iris.analysis.MEAN)
    EARTH3DJF = EARTH3DJF.collapsed('year', iris.analysis.MEAN)                       
    GFDLDJF = GFDLDJF.collapsed('year', iris.analysis.MEAN)
    HadGEM2DJF = HadGEM2DJF.collapsed('year', iris.analysis.MEAN)
    IPSLGDJF = IPSLGDJF.collapsed('year', iris.analysis.MEAN)
    MIROCGDJF = MIROCGDJF.collapsed('year', iris.analysis.MEAN)
    MPIDJF = MPIDJF.collapsed('year', iris.analysis.MEAN)
    NorESM1DJF = NorESM1DJF.collapsed('year', iris.analysis.MEAN)      
    
    CRUEDJF = CRUEDJF.collapsed('year', iris.analysis.MEAN)  
    CRUDJF = CRUDJF.collapsed('year', iris.analysis.MEAN)  
    CRUGDJF = CRUGDJF.collapsed('year', iris.analysis.MEAN)  
            
    CCCmaMAM = CCCmaMAM.collapsed('year', iris.analysis.MEAN)
    CLMcomMAM = CLMcomMAM.collapsed('year', iris.analysis.MEAN)    
    DMIMAM = DMIMAM.collapsed('year', iris.analysis.MEAN)    
    KNMIMAM = KNMIMAM.collapsed('year', iris.analysis.MEAN)
    MPIEMAM = MPIEMAM.collapsed('year', iris.analysis.MEAN)    
    SMHIMAM = SMHIMAM.collapsed('year', iris.analysis.MEAN)    
    
    CCCmaCanRCMMAM = CCCmaCanRCMMAM.collapsed('year', iris.analysis.MEAN)
    CCCmaSMHIMAM = CCCmaSMHIMAM.collapsed('year', iris.analysis.MEAN)
    CNRMMAM = CNRMMAM.collapsed('year', iris.analysis.MEAN)
    CNRMSMHIMAM = CNRMSMHIMAM.collapsed('year', iris.analysis.MEAN)
    CSIROMAM = CSIROMAM.collapsed('year', iris.analysis.MEAN)
    ICHECDMIMAM = ICHECDMIMAM.collapsed('year', iris.analysis.MEAN)
    ICHECCCLMMAM = ICHECCCLMMAM.collapsed('year', iris.analysis.MEAN)
    ICHECKNMIMAM = ICHECKNMIMAM.collapsed('year', iris.analysis.MEAN)
    ICHECMPIMAM = ICHECMPIMAM.collapsed('year', iris.analysis.MEAN)
    ICHECSMHIMAM = ICHECSMHIMAM.collapsed('year', iris.analysis.MEAN)
    IPSLMAM = IPSLMAM.collapsed('year', iris.analysis.MEAN)
    MIROCMAM = MIROCMAM.collapsed('year', iris.analysis.MEAN)
    MOHCCCLMMAM = MOHCCCLMMAM.collapsed('year', iris.analysis.MEAN)
    MOHCKNMIMAM = MOHCKNMIMAM.collapsed('year', iris.analysis.MEAN)
    MOHCSMHIMAM = MOHCSMHIMAM.collapsed('year', iris.analysis.MEAN)
    MPICCLMMAM = MPICCLMMAM.collapsed('year', iris.analysis.MEAN)
    MPIREMOMAM = MPIREMOMAM.collapsed('year', iris.analysis.MEAN)
    MPISMHIMAM = MPISMHIMAM.collapsed('year', iris.analysis.MEAN)
    NCCDMIMAM = NCCDMIMAM.collapsed('year', iris.analysis.MEAN)
    NCCSMHIMAM = NCCSMHIMAM.collapsed('year', iris.analysis.MEAN)
    NOAAMAM = NOAAMAM.collapsed('year', iris.analysis.MEAN)
    
    CanESM2MAM = CanESM2MAM.collapsed('year', iris.analysis.MEAN)
    CNRMGMAM = CNRMGMAM.collapsed('year', iris.analysis.MEAN)
    MK3MAM = MK3MAM.collapsed('year', iris.analysis.MEAN)
    EARTHMAM = EARTHMAM.collapsed('year', iris.analysis.MEAN)
    EARTH3MAM = EARTH3MAM.collapsed('year', iris.analysis.MEAN)                       
    GFDLMAM = GFDLMAM.collapsed('year', iris.analysis.MEAN)
    HadGEM2MAM = HadGEM2MAM.collapsed('year', iris.analysis.MEAN)
    IPSLGMAM = IPSLGMAM.collapsed('year', iris.analysis.MEAN)
    MIROCGMAM = MIROCGMAM.collapsed('year', iris.analysis.MEAN)
    MPIMAM = MPIMAM.collapsed('year', iris.analysis.MEAN)
    NorESM1MAM = NorESM1MAM.collapsed('year', iris.analysis.MEAN)   
    
    CRUEMAM = CRUEMAM.collapsed('year', iris.analysis.MEAN)  
    CRUMAM = CRUMAM.collapsed('year', iris.analysis.MEAN)  
    CRUGMAM = CRUGMAM.collapsed('year', iris.analysis.MEAN)  
           
    CCCmaJJA = CCCmaJJA.collapsed('year', iris.analysis.MEAN)
    CLMcomJJA = CLMcomJJA.collapsed('year', iris.analysis.MEAN)    
    DMIJJA = DMIJJA.collapsed('year', iris.analysis.MEAN)    
    KNMIJJA = KNMIJJA.collapsed('year', iris.analysis.MEAN)
    MPIEJJA = MPIEJJA.collapsed('year', iris.analysis.MEAN)    
    SMHIJJA = SMHIJJA.collapsed('year', iris.analysis.MEAN)    
    
    CCCmaCanRCMJJA = CCCmaCanRCMJJA.collapsed('year', iris.analysis.MEAN)
    CCCmaSMHIJJA = CCCmaSMHIJJA.collapsed('year', iris.analysis.MEAN)
    CNRMJJA = CNRMJJA.collapsed('year', iris.analysis.MEAN)
    CNRMSMHIJJA = CNRMSMHIJJA.collapsed('year', iris.analysis.MEAN)
    CSIROJJA = CSIROJJA.collapsed('year', iris.analysis.MEAN)
    ICHECDMIJJA = ICHECDMIJJA.collapsed('year', iris.analysis.MEAN)
    ICHECCCLMJJA = ICHECCCLMJJA.collapsed('year', iris.analysis.MEAN)
    ICHECKNMIJJA = ICHECKNMIJJA.collapsed('year', iris.analysis.MEAN)
    ICHECMPIJJA = ICHECMPIJJA.collapsed('year', iris.analysis.MEAN)
    ICHECSMHIJJA = ICHECSMHIJJA.collapsed('year', iris.analysis.MEAN)
    IPSLJJA = IPSLJJA.collapsed('year', iris.analysis.MEAN)
    MIROCJJA = MIROCJJA.collapsed('year', iris.analysis.MEAN)
    MOHCCCLMJJA = MOHCCCLMJJA.collapsed('year', iris.analysis.MEAN)
    MOHCKNMIJJA = MOHCKNMIJJA.collapsed('year', iris.analysis.MEAN)
    MOHCSMHIJJA = MOHCSMHIJJA.collapsed('year', iris.analysis.MEAN)
    MPICCLMJJA = MPICCLMJJA.collapsed('year', iris.analysis.MEAN)
    MPIREMOJJA = MPIREMOJJA.collapsed('year', iris.analysis.MEAN)
    MPISMHIJJA = MPISMHIJJA.collapsed('year', iris.analysis.MEAN)
    NCCDMIJJA = NCCDMIJJA.collapsed('year', iris.analysis.MEAN)
    NCCSMHIJJA = NCCSMHIJJA.collapsed('year', iris.analysis.MEAN)
    NOAAJJA = NOAAJJA.collapsed('year', iris.analysis.MEAN)
    
    CanESM2JJA = CanESM2JJA.collapsed('year', iris.analysis.MEAN)
    CNRMGJJA = CNRMGJJA.collapsed('year', iris.analysis.MEAN)
    MK3JJA = MK3JJA.collapsed('year', iris.analysis.MEAN)
    EARTHJJA = EARTHJJA.collapsed('year', iris.analysis.MEAN)
    EARTH3JJA = EARTH3JJA.collapsed('year', iris.analysis.MEAN)                       
    GFDLJJA = GFDLJJA.collapsed('year', iris.analysis.MEAN)
    HadGEM2JJA = HadGEM2JJA.collapsed('year', iris.analysis.MEAN)
    IPSLGJJA = IPSLGJJA.collapsed('year', iris.analysis.MEAN)
    MIROCGJJA = MIROCGJJA.collapsed('year', iris.analysis.MEAN)
    MPIJJA = MPIJJA.collapsed('year', iris.analysis.MEAN)
    NorESM1JJA = NorESM1JJA.collapsed('year', iris.analysis.MEAN)  
    
    CRUEJJA = CRUEJJA.collapsed('year', iris.analysis.MEAN)  
    CRUJJA = CRUJJA.collapsed('year', iris.analysis.MEAN)  
    CRUGJJA = CRUGJJA.collapsed('year', iris.analysis.MEAN)  
           
    CCCmaYR = CCCma.collapsed('year', iris.analysis.MEAN)
    CLMcomYR = CLMcom.collapsed('year', iris.analysis.MEAN)    
    DMIYR = DMI.collapsed('year', iris.analysis.MEAN)    
    KNMIYR = KNMI.collapsed('year', iris.analysis.MEAN)
    MPIEYR = MPIE.collapsed('year', iris.analysis.MEAN)    
    SMHIYR = SMHI.collapsed('year', iris.analysis.MEAN)    
    
    CCCmaCanRCMYR = CCCmaCanRCM.collapsed('year', iris.analysis.MEAN)
    CCCmaSMHIYR = CCCmaSMHI.collapsed('year', iris.analysis.MEAN)
    CNRMYR = CNRM.collapsed('year', iris.analysis.MEAN)
    CNRMSMHIYR = CNRMSMHI.collapsed('year', iris.analysis.MEAN)
    CSIROYR = CSIRO.collapsed('year', iris.analysis.MEAN)
    ICHECDMIYR = ICHECDMI.collapsed('year', iris.analysis.MEAN)
    ICHECCCLMYR = ICHECCCLM.collapsed('year', iris.analysis.MEAN)
    ICHECKNMIYR = ICHECKNMI.collapsed('year', iris.analysis.MEAN)
    ICHECMPIYR = ICHECMPI.collapsed('year', iris.analysis.MEAN)
    ICHECSMHIYR = ICHECSMHI.collapsed('year', iris.analysis.MEAN)
    IPSLYR = IPSL.collapsed('year', iris.analysis.MEAN)
    MIROCYR = MIROC.collapsed('year', iris.analysis.MEAN)
    MOHCCCLMYR = MOHCCCLM.collapsed('year', iris.analysis.MEAN)
    MOHCKNMIYR = MOHCKNMI.collapsed('year', iris.analysis.MEAN)
    MOHCSMHIYR = MOHCSMHI.collapsed('year', iris.analysis.MEAN)
    MPICCLMYR = MPICCLM.collapsed('year', iris.analysis.MEAN)
    MPIREMOYR = MPIREMO.collapsed('year', iris.analysis.MEAN)
    MPISMHIYR = MPISMHI.collapsed('year', iris.analysis.MEAN)
    NCCDMIYR = NCCDMI.collapsed('year', iris.analysis.MEAN)
    NCCSMHIYR = NCCSMHI.collapsed('year', iris.analysis.MEAN)
    NOAAYR = NOAA.collapsed('year', iris.analysis.MEAN)
    
    CanESM2YR=CanESM2.collapsed('year', iris.analysis.MEAN)
    CNRMGYR = CNRMG.collapsed('year', iris.analysis.MEAN)
    MK3YR = MK3.collapsed('year', iris.analysis.MEAN)
    EARTHYR = EARTH.collapsed('year', iris.analysis.MEAN)
    EARTH3YR = EARTH3.collapsed('year', iris.analysis.MEAN)                       
    GFDLYR = GFDL.collapsed('year', iris.analysis.MEAN)
    HadGEM2YR = HadGEM2.collapsed('year', iris.analysis.MEAN)
    IPSLGYR = IPSLG.collapsed('year', iris.analysis.MEAN)
    MIROCGYR = MIROCG.collapsed('year', iris.analysis.MEAN)
    MPIYR = MPI.collapsed('year', iris.analysis.MEAN)
    NorESM1YR = NorESM1.collapsed('year', iris.analysis.MEAN)
    
    CRUEYR = CRUE.collapsed('year', iris.analysis.MEAN) 
    CRUYR = CRU.collapsed('year', iris.analysis.MEAN) 
    CRUGYR = CRUG.collapsed('year', iris.analysis.MEAN) 
          
    #make units match, they already do but the name of the units doesn't. 
    CCCmaSON.units = Unit('Celsius')
    CLMcomSON.units = Unit('Celsius')
    DMISON.units = Unit('Celsius')
    KNMISON.units = Unit('Celsius')
    MPIESON.units = Unit('Celsius')
    SMHISON.units = Unit('Celsius')
    
    CCCmaCanRCMSON.units = Unit('Celsius')
    CCCmaSMHISON.units = Unit('Celsius')
    CNRMSON.units = Unit('Celsius')
    CNRMSMHISON.units = Unit('Celsius')
    CSIROSON.units = Unit('Celsius')
    ICHECDMISON.units = Unit('Celsius')
    ICHECCCLMSON.units = Unit('Celsius')
    ICHECKNMISON.units = Unit('Celsius')
    ICHECMPISON.units = Unit('Celsius')
    ICHECSMHISON.units = Unit('Celsius')
    IPSLSON.units = Unit('Celsius')
    MIROCSON.units = Unit('Celsius')
    MOHCCCLMSON.units = Unit('Celsius')
    MOHCKNMISON.units = Unit('Celsius')
    MOHCSMHISON.units = Unit('Celsius')
    MPICCLMSON.units = Unit('Celsius')
    MPIREMOSON.units = Unit('Celsius')
    MPISMHISON.units = Unit('Celsius')
    NCCDMISON.units = Unit('Celsius')
    NCCSMHISON.units = Unit('Celsius')
    NOAASON.units = Unit('Celsius')
    
    CanESM2SON.units = Unit('Celsius')
    CNRMGSON.units = Unit('Celsius')
    MK3SON.units = Unit('Celsius')
    EARTHSON.units = Unit('Celsius')
    EARTH3SON.units = Unit('Celsius')                       
    GFDLSON.units = Unit('Celsius')
    HadGEM2SON.units = Unit('Celsius')
    IPSLGSON.units = Unit('Celsius')
    MIROCGSON.units = Unit('Celsius')
    MPISON.units = Unit('Celsius')
    NorESM1SON.units = Unit('Celsius')
                           
    CRUESON.units = Unit('Celsius')  
    CRUSON.units = Unit('Celsius')  
    CRUGSON.units = Unit('Celsius')  
    
    CCCmaDJF.units = Unit('Celsius')
    CLMcomDJF.units = Unit('Celsius')
    DMIDJF.units = Unit('Celsius')
    KNMIDJF.units = Unit('Celsius')
    MPIEDJF.units = Unit('Celsius')
    SMHIDJF.units = Unit('Celsius')
    
    CCCmaCanRCMDJF.units = Unit('Celsius')
    CCCmaSMHIDJF.units = Unit('Celsius')
    CNRMDJF.units = Unit('Celsius')
    CNRMSMHIDJF.units = Unit('Celsius')
    CSIRODJF.units = Unit('Celsius')
    ICHECDMIDJF.units = Unit('Celsius')
    ICHECCCLMDJF.units = Unit('Celsius')
    ICHECKNMIDJF.units = Unit('Celsius')
    ICHECMPIDJF.units = Unit('Celsius')
    ICHECSMHIDJF.units = Unit('Celsius')
    IPSLDJF.units = Unit('Celsius')
    MIROCDJF.units = Unit('Celsius')
    MOHCCCLMDJF.units = Unit('Celsius')
    MOHCKNMIDJF.units = Unit('Celsius')
    MOHCSMHIDJF.units = Unit('Celsius')
    MPICCLMDJF.units = Unit('Celsius')
    MPIREMODJF.units = Unit('Celsius')
    MPISMHIDJF.units = Unit('Celsius')
    NCCDMIDJF.units = Unit('Celsius')
    NCCSMHIDJF.units = Unit('Celsius')
    NOAADJF.units = Unit('Celsius')
    
    CanESM2DJF.units = Unit('Celsius')
    CNRMGDJF.units = Unit('Celsius')
    MK3DJF.units = Unit('Celsius')
    EARTHDJF.units = Unit('Celsius')
    EARTH3DJF.units = Unit('Celsius')
    GFDLDJF.units = Unit('Celsius')
    HadGEM2DJF.units = Unit('Celsius')
    IPSLGDJF.units = Unit('Celsius')
    MIROCGDJF.units = Unit('Celsius')
    MPIDJF.units = Unit('Celsius')
    NorESM1DJF.units = Unit('Celsius')
                           
    CRUEDJF.units = Unit('Celsius') 
    CRUDJF.units = Unit('Celsius')
    CRUGDJF.units = Unit('Celsius')
    
    CCCmaMAM.Aunits = Unit('Celsius')
    CLMcomMAM.units = Unit('Celsius')
    DMIMAM.units = Unit('Celsius')
    KNMIMAM.units = Unit('Celsius')
    MPIEMAM.units = Unit('Celsius')
    SMHIMAM.units = Unit('Celsius')
    
    CCCmaCanRCMMAM.units = Unit('Celsius')
    CCCmaSMHIMAM.units = Unit('Celsius')
    CNRMMAM.units = Unit('Celsius')
    CNRMSMHIMAM.units = Unit('Celsius')
    CSIROMAM.units = Unit('Celsius')
    ICHECDMIMAM.units = Unit('Celsius')
    ICHECCCLMMAM.units = Unit('Celsius')
    ICHECKNMIMAM.units = Unit('Celsius')
    ICHECMPIMAM.units = Unit('Celsius')
    ICHECSMHIMAM.units = Unit('Celsius')
    IPSLMAM.units = Unit('Celsius')
    MIROCMAM.units = Unit('Celsius')
    MOHCCCLMMAM.units = Unit('Celsius')
    MOHCKNMIMAM.units = Unit('Celsius')
    MOHCSMHIMAM.units = Unit('Celsius')
    MPICCLMMAM.units = Unit('Celsius')
    MPIREMOMAM.units = Unit('Celsius')
    MPISMHIMAM.units = Unit('Celsius')
    NCCDMIMAM.units = Unit('Celsius')
    NCCSMHIMAM.units = Unit('Celsius')
    NOAAMAM.units = Unit('Celsius')
    
    CanESM2MAM.units = Unit('Celsius')
    CNRMGMAM.units = Unit('Celsius')
    MK3MAM.units = Unit('Celsius')
    EARTHMAM.units = Unit('Celsius')
    EARTH3MAM.units = Unit('Celsius')
    GFDLMAM.units = Unit('Celsius')
    HadGEM2MAM.units = Unit('Celsius')
    IPSLGMAM.units = Unit('Celsius')
    MIROCGMAM.units = Unit('Celsius')
    MPIMAM.units = Unit('Celsius')
    NorESM1MAM.units = Unit('Celsius')
                           
    CRUEMAM.units = Unit('Celsius')  
    CRUMAM.units = Unit('Celsius')
    CRUGMAM.units = Unit('Celsius')
    
    CCCmaJJA.units = Unit('Celsius')
    CLMcomJJA.units = Unit('Celsius')
    DMIJJA.units = Unit('Celsius')
    KNMIJJA.units = Unit('Celsius')
    MPIEJJA.units = Unit('Celsius')
    SMHIJJA.units = Unit('Celsius')
    
    CCCmaCanRCMJJA.units = Unit('Celsius')
    CCCmaSMHIJJA.units = Unit('Celsius')
    CNRMJJA.units = Unit('Celsius')
    CNRMSMHIJJA.units = Unit('Celsius')
    CSIROJJA.units = Unit('Celsius')
    ICHECDMIJJA.units = Unit('Celsius')
    ICHECCCLMJJA.units = Unit('Celsius')
    ICHECKNMIJJA.units = Unit('Celsius')
    ICHECMPIJJA.units = Unit('Celsius')
    ICHECSMHIJJA.units = Unit('Celsius')
    IPSLJJA.units = Unit('Celsius')
    MIROCJJA.units = Unit('Celsius')
    MOHCCCLMJJA.units = Unit('Celsius')
    MOHCKNMIJJA.units = Unit('Celsius')
    MOHCSMHIJJA.units = Unit('Celsius')
    MPICCLMJJA.units = Unit('Celsius')
    MPIREMOJJA.units = Unit('Celsius')
    MPISMHIJJA.units = Unit('Celsius')
    NCCDMIJJA.units = Unit('Celsius')
    NCCSMHIJJA.units = Unit('Celsius')
    NOAAJJA.units = Unit('Celsius')
    
    CanESM2JJA.units = Unit('Celsius')
    CNRMGJJA.units = Unit('Celsius')
    MK3JJA.units = Unit('Celsius')
    EARTHJJA.units = Unit('Celsius')
    EARTH3JJA.units = Unit('Celsius')
    GFDLJJA.units = Unit('Celsius')
    HadGEM2JJA.units = Unit('Celsius')
    IPSLGJJA.units = Unit('Celsius')
    MIROCGJJA.units = Unit('Celsius')
    MPIJJA.units = Unit('Celsius')
    NorESM1JJA.units = Unit('Celsius')       
                           
    CRUEJJA.units = Unit('Celsius')
    CRUJJA.units = Unit('Celsius')  
    CRUGJJA.units = Unit('Celsius')  
    
    CCCmaYR.units = Unit('Celsius')
    CLMcomYR.units = Unit('Celsius')
    DMIYR.units = Unit('Celsius')
    KNMIYR.units = Unit('Celsius')
    MPIEYR.units = Unit('Celsius')
    SMHIYR.units = Unit('Celsius')
    
    CCCmaCanRCMYR.units = Unit('Celsius')
    CCCmaSMHIYR.units = Unit('Celsius')
    CNRMYR.units = Unit('Celsius')
    CNRMSMHIYR.units = Unit('Celsius')
    CSIROYR.units = Unit('Celsius')
    ICHECDMIYR.units = Unit('Celsius')
    ICHECCCLMYR.units = Unit('Celsius')
    ICHECKNMIYR.units = Unit('Celsius')
    ICHECMPIYR.units = Unit('Celsius')
    ICHECSMHIYR.units = Unit('Celsius')
    IPSLYR.units = Unit('Celsius')
    MIROCYR.units = Unit('Celsius')
    MOHCCCLMYR.units = Unit('Celsius')
    MOHCKNMIYR.units = Unit('Celsius')
    MOHCSMHIYR.units = Unit('Celsius')
    MPICCLMYR.units = Unit('Celsius')
    MPIREMOYR.units = Unit('Celsius')
    MPISMHIYR.units = Unit('Celsius')
    NCCDMIYR.units = Unit('Celsius')
    NCCSMHIYR.units = Unit('Celsius')
    NOAAYR.units = Unit('Celsius')
    
    CanESM2YR.units = Unit('Celsius')
    CNRMGYR.units = Unit('Celsius')
    MK3YR.units = Unit('Celsius')
    EARTHYR.units = Unit('Celsius')
    EARTH3YR.units = Unit('Celsius')
    GFDLYR.units = Unit('Celsius')
    HadGEM2YR.units = Unit('Celsius')
    IPSLGYR.units = Unit('Celsius')
    MIROCGYR.units = Unit('Celsius')
    MPIYR.units = Unit('Celsius')
    NorESM1YR.units = Unit('Celsius')
    
    CRUEYR.units = Unit('Celsius')
    CRUYR.units = Unit('Celsius')  
    CRUGYR.units = Unit('Celsius')  
    
    #Create averages
    AverageSONE = (CCCmaSON + CLMcomSON + DMISON + KNMISON + MPIESON + SMHISON)/6.
    AverageDJFE = (CCCmaDJF + CLMcomDJF + DMIDJF + KNMIDJF + MPIEDJF + SMHIDJF)/6.
    AverageMAME = (CCCmaMAM + CLMcomMAM + DMIMAM + KNMIMAM + MPIEMAM + SMHIMAM)/6.
    AverageJJAE = (CCCmaJJA + CLMcomJJA + DMIJJA + KNMIJJA + MPIEJJA + SMHIJJA)/6.
    AverageE = (CCCmaYR + CLMcomYR + DMIYR + KNMIYR + MPIEYR + SMHIYR)/6.
    
    AverageSONR = (CCCmaCanRCMSON + CCCmaSMHISON + CNRMSON + CNRMSMHISON + CSIROSON + ICHECDMISON + ICHECCCLMSON + ICHECKNMISON + ICHECMPISON + ICHECSMHISON + IPSLSON + MIROCSON + MOHCCCLMSON + MOHCKNMISON + MOHCSMHISON + MPICCLMSON + MPIREMOSON + MPISMHISON + NCCDMISON + NCCSMHISON + NOAASON)/21.
    AverageDJFR = (CCCmaCanRCMDJF + CCCmaSMHIDJF + CNRMDJF + CNRMSMHIDJF + CSIRODJF + ICHECDMIDJF + ICHECCCLMDJF + ICHECKNMIDJF + ICHECMPIDJF + ICHECSMHIDJF + IPSLDJF + MIROCDJF + MOHCCCLMDJF + MOHCKNMIDJF + MOHCSMHIDJF + MPICCLMDJF + MPIREMODJF + MPISMHIDJF + NCCDMIDJF + NCCSMHIDJF + NOAADJF)/21.
    AverageMAMR = (CCCmaCanRCMMAM + CCCmaSMHIMAM + CNRMMAM + CNRMSMHIMAM + CSIROMAM + ICHECDMIMAM + ICHECCCLMMAM + ICHECKNMIMAM + ICHECMPIMAM + ICHECSMHIMAM + IPSLMAM + MIROCMAM + MOHCCCLMMAM + MOHCKNMIMAM + MOHCSMHIMAM + MPICCLMMAM + MPIREMOMAM + MPISMHIMAM + NCCDMIMAM + NCCSMHIMAM + NOAAMAM)/21.
    AverageJJAR = (CCCmaCanRCMJJA + CCCmaSMHIJJA + CNRMJJA + CNRMSMHIJJA + CSIROJJA + ICHECDMIJJA + ICHECCCLMJJA + ICHECKNMIJJA + ICHECMPIJJA + ICHECSMHIJJA + IPSLJJA + MIROCJJA + MOHCCCLMJJA + MOHCKNMIJJA + MOHCSMHIJJA + MPICCLMJJA + MPIREMOJJA + MPISMHIJJA + NCCDMIJJA + NCCSMHIJJA + NOAAJJA)/21.
    AverageR = (CCCmaCanRCMYR + CCCmaSMHIYR + CNRMYR + CNRMSMHIYR + CSIROYR + ICHECDMIYR + ICHECCCLMYR + ICHECKNMIYR + ICHECMPIYR + ICHECSMHIYR + IPSLYR + MIROCYR + MOHCCCLMYR + MOHCKNMIYR + MOHCSMHIYR + MPICCLMYR + MPIREMOYR + MPISMHIYR + NCCDMIYR + NCCSMHIYR + NOAAYR)/21.
     
    AverageSONG = (CanESM2SON + CNRMGSON + MK3SON + EARTHSON + EARTH3SON + GFDLSON + HadGEM2SON + IPSLGSON + MIROCGSON + MPISON + NorESM1SON)/11.
    AverageDJFG = (CanESM2DJF + CNRMGDJF + MK3DJF + EARTHDJF + EARTH3DJF + GFDLDJF + HadGEM2DJF + IPSLGDJF + MIROCGDJF + MPIDJF + NorESM1DJF)/11.
    AverageMAMG = (CanESM2MAM + CNRMGMAM + MK3MAM + EARTHMAM + EARTH3MAM + GFDLMAM + HadGEM2MAM + IPSLGMAM + MIROCGMAM + MPIMAM + NorESM1MAM)/11.
    AverageJJAG = (CanESM2JJA + CNRMGJJA + MK3JJA + EARTHJJA + EARTH3JJA + GFDLJJA + HadGEM2JJA + IPSLGJJA + MIROCGJJA + MPIJJA + NorESM1JJA)/11.
    AverageG = (CanESM2YR + CNRMGYR + MK3YR + EARTHYR + EARTH3YR + GFDLYR + HadGEM2YR + IPSLGYR + MIROCGYR + MPIYR + NorESM1YR)/11.
    
    ObsSONE = (CRUESON)
    ObsDJFE = (CRUEDJF)
    ObsMAME = (CRUEMAM)
    ObsJJAE = (CRUEJJA)
    ObsE = (CRUEYR)
    
    ObsSON = (CRUSON)
    ObsDJF = (CRUDJF)
    ObsMAM = (CRUMAM)
    ObsJJA = (CRUJJA)   
    Obs = (CRUYR)
    
    ObsSONG = (CRUGSON)
    ObsDJFG = (CRUGDJF)
    ObsMAMG = (CRUGMAM)
    ObsJJAG = (CRUGJJA)   
    ObsG = (CRUGYR)
        
    #fix unit names, data in Celsius, but called 'Kelvin'...     
    AverageSONE.units = Unit('Celsius')
    AverageDJFE.units = Unit('Celsius')
    AverageMAME.units = Unit('Celsius')
    AverageJJAE.units = Unit('Celsius')
    AverageE.units = Unit('Celsius')
    
    AverageSONR.units = Unit('Celsius')
    AverageDJFR.units = Unit('Celsius')
    AverageMAMR.units = Unit('Celsius')
    AverageJJAR.units = Unit('Celsius')
    AverageR.units = Unit('Celsius')
    
    AverageSONG.units = Unit('Celsius')
    AverageDJFG.units = Unit('Celsius')
    AverageMAMG.units = Unit('Celsius')
    AverageJJAG.units = Unit('Celsius')
    AverageG.units = Unit('Celsius')
    
    ObsSON.units = Unit('Celsius')
    ObsDJF.units = Unit('Celsius')
    ObsMAM.units = Unit('Celsius')
    ObsJJA.units = Unit('Celsius')
    Obs.units = Unit('Celsius')
    
    ObsSONG.units = Unit('Celsius')
    ObsDJFG.units = Unit('Celsius')
    ObsMAMG.units = Unit('Celsius')
    ObsJJAG.units = Unit('Celsius')
    ObsG.units = Unit('Celsius')
    
    ObsSONE.units = Unit('Celsius')
    ObsDJFE.units = Unit('Celsius')
    ObsMAME.units = Unit('Celsius')
    ObsJJAE.units = Unit('Celsius')
    ObsE.units = Unit('Celsius')
    
    #Take difference between two datasets
    BiasSONE = AverageSONE-ObsSONE
    BiasDJFE = AverageDJFE-ObsDJFE
    BiasMAME = AverageMAME-ObsMAME
    BiasJJAE = AverageJJAE-ObsJJAE
    BiasE = AverageE-ObsE
    
    BiasSONR = AverageSONR-ObsSON
    BiasDJFR = AverageDJFR-ObsDJF
    BiasMAMR = AverageMAMR-ObsMAM
    BiasJJAR = AverageJJAR-ObsJJA
    BiasR = AverageR-Obs
    
    BiasSONG = AverageSONG-ObsSONG
    BiasDJFG = AverageDJFG-ObsDJFG
    BiasMAMG = AverageMAMG-ObsMAMG
    BiasJJAG = AverageJJAG-ObsJJAG
    BiasG = AverageG-ObsG
    
    print "ERAINT"
    print np.amax(AverageE.data)
    print np.amin(AverageE.data)
    
    print "ERAINTSON"
    print np.amax(AverageSONE.data)
    print np.amin(AverageSONE.data)
    
    print "ERAINTDJF"
    print np.amax(AverageDJFE.data)
    print np.amin(AverageDJFE.data)
    
    print "ERAINTMAM"
    print np.amax(AverageMAME.data)
    print np.amin(AverageMAME.data)
    
    print "ERAINTJJA"
    print np.amax(AverageJJAE.data)
    print np.amin(AverageJJAE.data)
    
    print "RCM"
    print np.amax(AverageR.data)
    print np.amin(AverageR.data)
    
    print "RCMSON"
    print np.amax(AverageSONR.data)
    print np.amin(AverageSONR.data)
    
    print "RCMDJF"
    print np.amax(AverageDJFR.data)
    print np.amin(AverageDJFR.data)
    
    print "CORDEXMAM"
    print np.amax(AverageMAMR.data)
    print np.amin(AverageMAMR.data)
    
    print "RCMJJA"
    print np.amax(AverageJJAR.data)
    print np.amin(AverageJJAR.data)
    
    print "GCM"
    print np.amax(AverageG.data)
    print np.amin(AverageG.data)
    
    print "GCMSON"
    print np.amax(AverageSONG.data)
    print np.amin(AverageSONG.data)
    
    print "GCMDJF"
    print np.amax(AverageDJFG.data)
    print np.amin(AverageDJFG.data)
    
    print "GCMMAM"
    print np.amax(AverageMAMG.data)
    print np.amin(AverageMAMG.data)
    
    print "GCMJJA"
    print np.amax(AverageJJAG.data)
    print np.amin(AverageJJAG.data)
    
    print "ObsE"
    print np.amax(ObsE.data)
    print np.amin(ObsE.data)
    
    print "ObsSONE"
    print np.amax(ObsSONE.data)
    print np.amin(ObsSONE.data)
    
    print "ObsDJFE"
    print np.amax(ObsDJFE.data)
    print np.amin(ObsDJFE.data)
    
    print "ObsMAME"
    print np.amax(ObsMAME.data)
    print np.amin(ObsMAME.data)

    print "ObsJJAE"
    print np.amax(ObsJJAE.data)
    print np.amin(ObsJJAE.data)
    
    print "Obs"
    print np.amax(Obs.data)
    print np.amin(Obs.data)
    
    print "ObsSON"
    print np.amax(ObsSON.data)
    print np.amin(ObsSON.data)
    
    print "ObsDJF"
    print np.amax(ObsDJF.data)
    print np.amin(ObsDJF.data)
    
    print "ObsMAM"
    print np.amax(ObsMAM.data)
    print np.amin(ObsMAM.data)

    print "ObsJJA"
    print np.amax(ObsJJA.data)
    print np.amin(ObsJJA.data)
    
    print "ObsG"
    print np.amax(ObsG.data)
    print np.amin(ObsG.data)
    
    print "ObsSONG"
    print np.amax(ObsSONG.data)
    print np.amin(ObsSONG.data)
    
    print "ObsDJFG"
    print np.amax(ObsDJFG.data)
    print np.amin(ObsDJFG.data)
    
    print "ObsMAM"
    print np.amax(ObsMAMG.data)
    print np.amin(ObsMAMG.data)

    print "ObsJJA"
    print np.amax(ObsJJAG.data)
    print np.amin(ObsJJAG.data)
    
    print "BiasE"
    print np.amax(BiasE.data)
    print np.amin(BiasE.data)

    print "BiasSONE"
    print np.amax(BiasSONE.data)
    print np.amin(BiasSONE.data)

    print "BiasDJFE"
    print np.amax(BiasDJFE.data)
    print np.amin(BiasDJFE.data)

    print "BiasMAME"
    print np.amax(BiasMAME.data)
    print np.amin(BiasMAME.data) 
   
    print "BiasJJAE"
    print np.amax(BiasJJAE.data)
    print np.amin(BiasJJAE.data)
    
    print "BiasR"
    print np.amax(BiasR.data)
    print np.amin(BiasR.data)

    print "BiasSONR"
    print np.amax(BiasSONR.data)
    print np.amin(BiasSONR.data)

    print "BiasDJFR"
    print np.amax(BiasDJFR.data)
    print np.amin(BiasDJFR.data)

    print "BiasMAMR"
    print np.amax(BiasMAMR.data)
    print np.amin(BiasMAMR.data) 
   
    print "BiasJJAR"
    print np.amax(BiasJJAR.data)
    print np.amin(BiasJJAR.data)
    
    print "BiasG"
    print np.amax(BiasG.data)
    print np.amin(BiasG.data)

    print "BiasSONG"
    print np.amax(BiasSONG.data)
    print np.amin(BiasSONG.data)

    print "BiasDJFG"
    print np.amax(BiasDJFG.data)
    print np.amin(BiasDJFG.data)

    print "BiasMAMG"
    print np.amax(BiasMAMG.data)
    print np.amin(BiasMAMG.data) 
   
    print "BiasJJAG"
    print np.amax(BiasJJAG.data)
    print np.amin(BiasJJAG.data)
    
    #---------------------------------------------------------------------------------------------------------------------
    #PART 4: PLOT MAP
    #load color palettes
    colourA = mpl_cm.get_cmap('YlOrRd')
    colourB = mpl_cm.get_cmap('bwr')
    
    #4Ai: ABSOLUTE TEMPERATURES - Annual
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(AverageE, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1990-2008 - Average_ERAINT', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Average_ERAINT_MAP_Annual_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(ObsE, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1990-2008 - Observed', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Observed_ERAINT_MAP_Annual_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(AverageR, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1961-2005 - Average_RCM', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Average_RCM_MAP_Annual_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(AverageG, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1961-2005 - Average_GCM', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Average_GCM_MAP_Annual', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(Obs, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1961-2005 - Observed', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Observed_MAP_Annual_H', bbox_inches='tight')
    plt.show()
    
    
    #-------------------
    #4Aii: ABSOLUTE TEMPERATURES -Seasonal
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(AverageSONE, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1990-2008 - Average_ERAINT SON', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Average_ERAINT_MAP_SON_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(AverageDJFE, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1990-2008 - Average_ERAINT DJF', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Average_ERAINT_MAP_DJF_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(AverageMAME, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1990-2008 - Average_ERAINT MAM', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Average_ERAINT_MAP_MAM_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(AverageJJAE, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1990-2008 - Average_ERAINT JJA', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Average_ERAINT_MAP_JJA_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(ObsSONE, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1990-2008 - Observed SON', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Observed_ERAINT_MAP_SON_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(ObsDJFE, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1990-2008 - Observed DJF', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Observed_ERAINT_MAP_DJF_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(ObsMAME, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1990-2008 - Observed MAM', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Observed_ERAINT_MAP_MAM_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(ObsJJAE, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1990-2008 - Observed JJA', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Observed_ERAINT_MAP_JJA_H', bbox_inches='tight')
    plt.show()
              
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(AverageSONR, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1961-2005 - Average_RCM_SON', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Average_RCM_MAP_SON_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(AverageDJFR, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1961-2005 - Average_RCM DJF', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Average_RCM_MAP_DJF_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(AverageMAMR, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1961-2005 - Average_RCM MAM', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Average_RCM_MAP_MAM_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(AverageJJAR, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1961-2005 - Average_RCM JJA', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Average_RCM_MAP_JJA_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(AverageSONG, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1961-2005 - Average_GCM_SON', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Average_GCM_MAP_SON', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(AverageDJFG, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1961-2005 - Average_GCM DJF', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Average_GCM_MAP_DJF', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(AverageMAMG, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1961-2005 - Average_GCM MAM', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Average_GCM_MAP_MAM', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(AverageJJAG, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1961-2005 - Average_GCM JJA', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Average_GCM_MAP_JJA', bbox_inches='tight')
    plt.show()
              
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(ObsSON, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1961-2005 - Observed SON', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Observed_MAP_SON_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(ObsDJF, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1961-2005 - Observed DJF', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Observed_MAP_DJF_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(ObsMAM, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1961-2005 - Observed MAM', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Observed_MAP_MAM_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(ObsJJA, cmap=colourA, levels=np.arange(20,36,0.5), extend='both')
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX 1961-2005 - Observed JJA', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Observed_MAP_JJA_H', bbox_inches='tight')
    plt.show()
    
    
    #-------------------
    #4B: DIFFERENCE FROM OBSERVED

    #4Bi: DIFFERENCE FROM OBSERVED - Annual
    #load color palette
        
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(BiasE, cmap=colourB, levels=np.arange(-5,5,0.25), extend='both') 
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX Bias 1990-2008 - Average_ERAINT', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAXBIAS_Average_ERAINT_MAP_Annual_H', bbox_inches='tight')
    plt.show()
              
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(BiasR, cmap=colourB, levels=np.arange(-5,5,0.25), extend='both') 
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX Bias 1961-2005 - Average_RCM', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAXBIAS_Average_RCM_MAP_Annual_H', bbox_inches='tight')
    plt.show()
              
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(BiasG, cmap=colourB, levels=np.arange(-5,5,0.25), extend='both') 
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX Bias 1961-2005 - Average_GCM', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAXBIAS_Average_GCM_MAP_Annual', bbox_inches='tight')
    plt.show()
    
    #-------------------
    #4Bii: DIFFERENCE FROM OBSERVED-Seasonal
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(BiasSONE, cmap=colourB, levels=np.arange(-5,5,0.25), extend='both') 
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX Bias SON 1990-2008 - Average_ERAINT', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAXBIAS_Average_ERAINT_MAP_SON_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(BiasDJFE, cmap=colourB, levels=np.arange(-5,5,0.25), extend='both') 
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX Bias DJF 1990-2008 - Average_ERAINT', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAXBIAS_Average_ERAINT_MAP_DJF_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(BiasMAME, cmap=colourB, levels=np.arange(-5,5,0.25), extend='both') 
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX Bias MAM 1990-2008 - Average_ERAINT', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAXBIAS_Average_ERAINT_MAP_MAM_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(BiasJJAE, cmap=colourB, levels=np.arange(-5,5,0.25), extend='both') 
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX Bias JJA 1990-2008 - Average_ERAINT', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAXBIAS_Average_ERAINT_MAP_JJA_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(BiasSONR, cmap=colourB, levels=np.arange(-5,5,0.25), extend='both') 
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX Bias SON 1961-2005 - Average_RCM', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAXBIAS_Average_RCM_MAP_SON_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(BiasDJFR, cmap=colourB, levels=np.arange(-5,5,0.25), extend='both') 
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX Bias DJF 1961-2005 - Average_RCM', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAXBIAS_Average_RCM_MAP_DJF_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(BiasMAMR, cmap=colourB, levels=np.arange(-5,5,0.25), extend='both') 
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX Bias MAM 1961-2005 - Average_RCM', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAXBIAS_Average_RCM_MAP_MAM_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(BiasJJAR, cmap=colourB, levels=np.arange(-5,5,0.25), extend='both') 
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX Bias JJA 1961-2005 - Average_RCM', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAXBIAS_Average_RCM_MAP_JJA_H', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(BiasSONG, cmap=colourB, levels=np.arange(-5,5,0.25), extend='both') 
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX Bias SON 1961-2005 - Average_GCM', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAXBIAS_Average_GCM_MAP_SON', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(BiasDJFG, cmap=colourB, levels=np.arange(-5,5,0.25), extend='both') 
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX Bias DJF 1961-2005 - Average_GCM', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAXBIAS_Average_GCM_MAP_DJF', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(BiasMAMG, cmap=colourB, levels=np.arange(-5,5,0.25), extend='both') 
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX Bias MAM 1961-2005 - Average_GCM', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAXBIAS_Average_GCM_MAP_MAM', bbox_inches='tight')
    plt.show()
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32.5, 36., -9, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #plot data and set colour range
    plot = iplt.contourf(BiasJJAG, cmap=colourB, levels=np.arange(-5,5,0.25), extend='both') 
    #add colour bar index and a label
    plt.colorbar(plot, label='Celsius')
    #give map a title
    plt.title('TasMAX Bias JJA 1961-2005 - Average_GCM', fontsize=10)
    #save the image of the graph and include full legend
    plt.savefig('TasMAXBIAS_Average_GCM_MAP_JJA', bbox_inches='tight')
    plt.show()
    
if __name__ == '__main__':
    main()