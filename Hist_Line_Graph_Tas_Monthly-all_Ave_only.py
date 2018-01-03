"""
Created on Wednesday August 2 2017

@author: s0899345
"""

import matplotlib.pyplot as plt
import iris
import iris.coord_categorisation as iriscc
import iris.plot as iplt
import iris.analysis.cartography
import numpy as np
import cf_units
from cf_units import Unit
import datetime
import calendar

#this file is split into parts as follows:
    #PART 1: load and format all models 
    #PART 2: load and format observed data
    #PART 3: format files to be plot specific
    #PART 4: plot data
    
def main():
    #-------------------------------------------------------------------------
    #PART 1: LOAD and FORMAT ALL MODELS   
    #bring in all the ERAINT models we need and give them a name
    CCCma = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/ERAINT/1979-2012/tas_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_CCCma-CanRCM4_r2_mon_198901-200912.nc'
    CLMcom = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/ERAINT/1979-2012/tas_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_198901-200812.nc'
    DMI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/ERAINT/1979-2012/tas_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_DMI-HIRHAM5_v2_mon_198901-201012.nc'
    KNMI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/ERAINT/1979-2012/tas_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_KNMI-RACMO22T_v1_mon_197901-201212.nc'
    MPIE = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/ERAINT/1979-2012/tas_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_MPI-CSC-REMO2009_v1_mon_198902-200812.nc'
    SMHI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/ERAINT/1979-2012/tas_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_SMHI-RCA4_v1_mon_198001-201012.nc'
    
    #bring in all the CORDEX RCM models we need and give them a name
    CCCmaCanRCM= '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/Historical/1950-2005/tas_AFR-44_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_mon_195001-200512.nc'
    CCCmaSMHI= '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/Historical/1950-2005/tas_AFR-44_CCCma-CanESM2_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'   
    CNRM= '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/Historical/1950-2005/tas_AFR-44_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_195101-200512.nc'
    CNRMSMHI= '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/Historical/1950-2005/tas_AFR-44_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    CSIRO = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/Historical/1950-2005/tas_AFR-44_CSIRO-QCCCE-CSIRO-Mk3-6-0_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    ICHECDMI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/Historical/1950-2005/tas_AFR-44_ICHEC-EC-EARTH_historical_r3i1p1_DMI-HIRHAM5_v2_mon_195101-200512.nc'   
    ICHECCCLM = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/Historical/1950-2005/tas_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_CLMcom-CCLM4-8-17_v1_mon_194912-200512.nc'    
    ICHECKNMI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/Historical/1950-2005/tas_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_KNMI-RACMO22T_v1_mon_195001-200512.nc'
    ICHECMPI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/Historical/1950-2005/tas_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_MPI-CSC-REMO2009_v1_mon_195001-200512.nc'
    ICHECSMHI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/Historical/1950-2005/tas_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    IPSL = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/Historical/1950-2005/tas_AFR-44_IPSL-IPSL-CM5A-MR_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    MIROC =  '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/Historical/1950-2005/tas_AFR-44_MIROC-MIROC5_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc' 
    MOHCCCLM = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/Historical/1950-2005/tas_AFR-44_MOHC-HadGEM2-ES_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_194912-200512.nc' 
    MOHCKNMI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/Historical/1950-2005/tas_AFR-44_MOHC-HadGEM2-ES_historical_r1i1p1_KNMI-RACMO22T_v2_mon_195001-200512.nc'
    MOHCSMHI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/Historical/1950-2005/tas_AFR-44_MOHC-HadGEM2-ES_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    MPICCLM = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/Historical/1950-2005/tas_AFR-44_MPI-M-MPI-ESM-LR_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_194912-200512.nc'     
    MPIREMO = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/Historical/1950-2005/tas_AFR-44_MPI-M-MPI-ESM-LR_historical_r1i1p1_MPI-CSC-REMO2009_v1_mon_195001-200512.nc'    
    MPISMHI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/Historical/1950-2005/tas_AFR-44_MPI-M-MPI-ESM-LR_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    NCCDMI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/Historical/1950-2005/tas_AFR-44_NCC-NorESM1-M_historical_r1i1p1_DMI-HIRHAM5_v1_mon_195101-200512.nc'    
    NCCSMHI = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/Historical/1950-2005/tas_AFR-44_NCC-NorESM1-M_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'   
    NOAA = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/AFR_44_tas/Historical/1950-2005/tas_AFR-44_NOAA-GFDL-GFDL-ESM2M_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'   
    
    #bring in all the GCM models we need and give them a name
    CanESM2= '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/GCM_data/tas_Amon_CanESM2_historical_r1i1p1_185001-200512.nc'
    CNRMG= '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/GCM_data/tas_Amon_CNRM-CM5-2_historical_r1i1p1_195001-200512.nc'
    MK3 = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/GCM_data/tas_Amon_CSIRO-Mk3-6-0_historical_r1i1p1_185001-200512.nc'    
    EARTH= '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/GCM_data/tas_Amon_EC-EARTH_historical_r12i1p1_195001-201212.nc'
    EARTH3= '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/GCM_data/tas_Amon_EC-EARTH_historical_r3i1p1_1850-2009.nc'
    GFDL='/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/GCM_data/tas_Amon_GFDL-ESM2M_historical_r1i1p1_194601-200512.nc'
    HadGEM2='/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/GCM_data/tas_Amon_HadGEM2-ES_historical_r1i1p1_193412-200511.nc'
    IPSLG='/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/GCM_data/tas_Amon_IPSL-CM5B-LR_historical_r1i1p1_185001-200512.nc'
    MIROCG='/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/GCM_data/tas_Amon_MIROC5_historical_r1i1p1_185001-201212.nc'
    MPI='/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/GCM_data/tas_Amon_MPI-ESM-MR_historical_r1i1p1_185001-200512.nc'
    NorESM1='/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/GCM_data/tas_Amon_NorESM1-M_historical_r1i1p1_185001-200512.nc'
    
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
    CRU= '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/Actual_Data/cru_ts4.00.1901.2015.tmp.dat.nc'
    UDel= '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/Actual_Data/UDel_air.mon.mean.v401.nc'
        
    #Load exactly one cube from given file
    CRU = iris.load_cube(CRU, 'near-surface temperature')
    UDel = iris.load_cube(UDel)
    
    #-------------------------------------------------------------------------
    #PART 3: FORMAT DATA TO BE PLOT SPECIFIC 
    #fix EARTH3 time units as they differ from all other models
    t_coord=EARTH3.coord('time')
    
    t_unit = t_coord.attributes['invalid_units']
    timestep, _, t_fmt_str = t_unit.split(' ')
    new_t_unit_str= '{} since 1850-01-01 00:00:00'.format(timestep) 
    new_t_unit = cf_units.Unit(new_t_unit_str, calendar=cf_units.CALENDAR_STANDARD)
    
    new_datetimes = [datetime.datetime.strptime(str(dt), t_fmt_str) for dt in t_coord.points]
    new_dt_points = [new_t_unit.date2num(new_dt) for new_dt in new_datetimes]
    new_t_coord = iris.coords.DimCoord(new_dt_points, standard_name='time', units=new_t_unit)

    t_coord_dim = EARTH3.coord_dims('time')
    EARTH3.remove_coord('time')
    EARTH3.add_dim_coord(new_t_coord,t_coord_dim)
    
    print "EARTH3 new time"
    print (EARTH3.coord('time'))
    
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
    UDel = UDel.regrid(CanESM2, iris.analysis.Linear())
    
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
    UDel = UDel.extract(Malawi)
    
    #time constraignt to make all series the same, for ERAINT this is 1990-2008 and for RCMs and GCMs this is 1951-2005
    iris.FUTURE.cell_datetime_objects = True
    t_constraint_ERAINT = iris.Constraint(time=lambda cell: 1990 <= cell.point.year <= 2005)
   
    CCCma = CCCma.extract(t_constraint_ERAINT)
    CLMcom =CLMcom.extract(t_constraint_ERAINT)
    DMI=DMI.extract(t_constraint_ERAINT)
    KNMI=KNMI.extract(t_constraint_ERAINT)
    MPIE=MPIE.extract(t_constraint_ERAINT)
    SMHI=SMHI.extract(t_constraint_ERAINT)
            
    CCCmaCanRCM = CCCmaCanRCM.extract(t_constraint_ERAINT)
    CCCmaSMHI = CCCmaSMHI.extract(t_constraint_ERAINT)
    CNRM =CNRM.extract(t_constraint_ERAINT)
    CNRMSMHI =CNRMSMHI.extract(t_constraint_ERAINT)
    CSIRO=CSIRO.extract(t_constraint_ERAINT)
    ICHECDMI=ICHECDMI.extract(t_constraint_ERAINT)
    ICHECCCLM=ICHECCCLM.extract(t_constraint_ERAINT)
    ICHECKNMI=ICHECKNMI.extract(t_constraint_ERAINT)
    ICHECMPI=ICHECMPI.extract(t_constraint_ERAINT)
    ICHECSMHI=ICHECSMHI.extract(t_constraint_ERAINT)
    IPSL=IPSL.extract(t_constraint_ERAINT)
    MIROC=MIROC.extract(t_constraint_ERAINT)
    MOHCCCLM=MOHCCCLM.extract(t_constraint_ERAINT)
    MOHCKNMI=MOHCKNMI.extract(t_constraint_ERAINT)
    MOHCSMHI=MOHCSMHI.extract(t_constraint_ERAINT)
    MPICCLM=MPICCLM.extract(t_constraint_ERAINT)
    MPIREMO=MPIREMO.extract(t_constraint_ERAINT)
    MPISMHI=MPISMHI.extract(t_constraint_ERAINT)
    NCCDMI=NCCDMI.extract(t_constraint_ERAINT)
    NCCSMHI=NCCSMHI.extract(t_constraint_ERAINT)
    NOAA=NOAA.extract(t_constraint_ERAINT) 
    
    CanESM2 = CanESM2.extract(t_constraint_ERAINT)
    CNRMG =CNRMG.extract(t_constraint_ERAINT)
    MK3 =MK3.extract(t_constraint_ERAINT)
    EARTH =EARTH.extract(t_constraint_ERAINT)
    EARTH3 =EARTH3.extract(t_constraint_ERAINT)
    GFDL=GFDL.extract(t_constraint_ERAINT)
    HadGEM2=HadGEM2.extract(t_constraint_ERAINT)
    IPSLG=IPSLG.extract(t_constraint_ERAINT)
    MIROCG=MIROCG.extract(t_constraint_ERAINT)
    MPI=MPI.extract(t_constraint_ERAINT)
    NorESM1=NorESM1.extract(t_constraint_ERAINT)
    
    CRU = CRU.extract(t_constraint_ERAINT)
    UDel = UDel.extract(t_constraint_ERAINT)
    
    #Convert units to match, CORDEX data is in Kelvin but Observed data in Celsius, we would like to show all data in Celsius
    CCCma.convert_units('Celsius')
    CLMcom.convert_units('Celsius')
    DMI.convert_units('Celsius')
    KNMI.convert_units('Celsius')
    MPI.convert_units('Celsius')
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
    
    #add month numbers to data
    iriscc.add_month_number(CCCma, 'time')
    iriscc.add_month_number(CLMcom, 'time')
    iriscc.add_month_number(DMI, 'time')
    iriscc.add_month_number(KNMI, 'time')
    iriscc.add_month_number(MPIE, 'time')
    iriscc.add_month_number(SMHI, 'time')
        
    iriscc.add_month_number(CCCmaCanRCM, 'time')
    iriscc.add_month_number(CCCmaSMHI, 'time')
    iriscc.add_month_number(CNRM, 'time')
    iriscc.add_month_number(CNRMSMHI, 'time')
    iriscc.add_month_number(CSIRO, 'time')
    iriscc.add_month_number(ICHECDMI, 'time')
    iriscc.add_month_number(ICHECCCLM, 'time')
    iriscc.add_month_number(ICHECKNMI, 'time')
    iriscc.add_month_number(ICHECMPI, 'time')
    iriscc.add_month_number(ICHECSMHI, 'time')
    iriscc.add_month_number(IPSL, 'time')
    iriscc.add_month_number(MIROC, 'time')
    iriscc.add_month_number(MOHCCCLM, 'time')
    iriscc.add_month_number(MOHCKNMI, 'time')
    iriscc.add_month_number(MOHCSMHI, 'time')
    iriscc.add_month_number(MPICCLM, 'time')
    iriscc.add_month_number(MPIREMO, 'time')
    iriscc.add_month_number(MPISMHI, 'time')
    iriscc.add_month_number(NCCDMI, 'time')
    iriscc.add_month_number(NCCSMHI, 'time')
    iriscc.add_month_number(NOAA, 'time')
    
    iriscc.add_month_number(CanESM2, 'time')
    iriscc.add_month_number(CNRMG, 'time')
    iriscc.add_month_number(MK3, 'time')
    iriscc.add_month_number(EARTH, 'time')
    iriscc.add_month_number(EARTH3, 'time')
    iriscc.add_month_number(GFDL, 'time')
    iriscc.add_month_number(HadGEM2, 'time')
    iriscc.add_month_number(IPSLG, 'time')
    iriscc.add_month_number(MIROCG, 'time')
    iriscc.add_month_number(MPI, 'time')
    iriscc.add_month_number(NorESM1, 'time')
    
    iriscc.add_month_number(CRU, 'time')
    iriscc.add_month_number(UDel, 'time')
    
    #We are interested in plotting the data by month, so we need to take a mean of all the data by month    
    CCCma = CCCma.aggregated_by('month_number', iris.analysis.MEAN)
    CLMcom = CLMcom.aggregated_by('month_number', iris.analysis.MEAN)    
    DMI = DMI.aggregated_by('month_number', iris.analysis.MEAN)    
    KNMI = KNMI.aggregated_by('month_number', iris.analysis.MEAN)
    MPIE = MPI.aggregated_by('month_number', iris.analysis.MEAN)    
    SMHI = SMHI.aggregated_by('month_number', iris.analysis.MEAN)    
    
    CCCmaCanRCM = CCCmaCanRCM.aggregated_by('month_number', iris.analysis.MEAN)
    CCCmaSMHI = CCCmaSMHI.aggregated_by('month_number', iris.analysis.MEAN)
    CNRM = CNRM.aggregated_by('month_number', iris.analysis.MEAN)
    CNRMSMHI = CNRMSMHI.aggregated_by('month_number', iris.analysis.MEAN)
    CSIRO = CSIRO.aggregated_by('month_number', iris.analysis.MEAN)
    ICHECDMI = ICHECDMI.aggregated_by('month_number', iris.analysis.MEAN)
    ICHECCCLM = ICHECCCLM.aggregated_by('month_number', iris.analysis.MEAN)
    ICHECKNMI = ICHECKNMI.aggregated_by('month_number', iris.analysis.MEAN)
    ICHECMPI = ICHECMPI.aggregated_by('month_number', iris.analysis.MEAN)
    ICHECSMHI = ICHECSMHI.aggregated_by('month_number', iris.analysis.MEAN)
    IPSL = IPSL.aggregated_by('month_number', iris.analysis.MEAN)
    MIROC = MIROC.aggregated_by('month_number', iris.analysis.MEAN)
    MOHCCCLM = MOHCCCLM.aggregated_by('month_number', iris.analysis.MEAN)
    MOHCKNMI = MOHCKNMI.aggregated_by('month_number', iris.analysis.MEAN)
    MOHCSMHI = MOHCSMHI.aggregated_by('month_number', iris.analysis.MEAN)
    MPICCLM = MPICCLM.aggregated_by('month_number', iris.analysis.MEAN)
    MPIREMO = MPIREMO.aggregated_by('month_number', iris.analysis.MEAN)
    MPISMHI = MPISMHI.aggregated_by('month_number', iris.analysis.MEAN)
    NCCDMI = NCCDMI.aggregated_by('month_number', iris.analysis.MEAN)
    NCCSMHI = NCCSMHI.aggregated_by('month_number', iris.analysis.MEAN)
    NOAA = NOAA.aggregated_by('month_number', iris.analysis.MEAN)
    
    CanESM2 = CanESM2.aggregated_by('month_number', iris.analysis.MEAN)
    CNRMG = CNRMG.aggregated_by('month_number', iris.analysis.MEAN)
    MK3 = MK3.aggregated_by('month_number', iris.analysis.MEAN)
    EARTH = EARTH.aggregated_by('month_number', iris.analysis.MEAN)
    EARTH3 = EARTH3.aggregated_by('month_number', iris.analysis.MEAN)
    GFDL = GFDL.aggregated_by('month_number', iris.analysis.MEAN)
    HadGEM2 = HadGEM2.aggregated_by('month_number', iris.analysis.MEAN)
    IPSLG = IPSLG.aggregated_by('month_number', iris.analysis.MEAN)
    MIROCG = MIROCG.aggregated_by('month_number', iris.analysis.MEAN)
    MPI = MPI.aggregated_by('month_number', iris.analysis.MEAN)
    NorESM1 = NorESM1.aggregated_by('month_number', iris.analysis.MEAN)   
     
    CRU = CRU.aggregated_by('month_number', iris.analysis.MEAN) 
    UDel = UDel.aggregated_by('month_number', iris.analysis.MEAN) 
    
    #Returns an array of area weights, with the same dimensions as the cube.
    CCCma_grid_areas = iris.analysis.cartography.area_weights(CCCma)
    CLMcom_grid_areas = iris.analysis.cartography.area_weights(CLMcom)
    DMI_grid_areas = iris.analysis.cartography.area_weights(DMI)
    KNMI_grid_areas = iris.analysis.cartography.area_weights(KNMI)
    MPIE_grid_areas = iris.analysis.cartography.area_weights(MPIE)
    SMHI_grid_areas = iris.analysis.cartography.area_weights(SMHI)
    
    CCCmaCanRCM_grid_areas = iris.analysis.cartography.area_weights(CCCmaCanRCM)
    CCCmaSMHI_grid_areas = iris.analysis.cartography.area_weights(CCCmaSMHI)
    CNRM_grid_areas = iris.analysis.cartography.area_weights(CNRM)
    CNRMSMHI_grid_areas = iris.analysis.cartography.area_weights(CNRMSMHI)
    CSIRO_grid_areas = iris.analysis.cartography.area_weights(CSIRO)
    ICHECDMI_grid_areas = iris.analysis.cartography.area_weights(ICHECDMI)
    ICHECCCLM_grid_areas = iris.analysis.cartography.area_weights(ICHECCCLM)
    ICHECKNMI_grid_areas = iris.analysis.cartography.area_weights(ICHECKNMI)
    ICHECMPI_grid_areas = iris.analysis.cartography.area_weights(ICHECMPI)
    ICHECSMHI_grid_areas = iris.analysis.cartography.area_weights(ICHECSMHI)
    IPSL_grid_areas = iris.analysis.cartography.area_weights(IPSL)
    MIROC_grid_areas = iris.analysis.cartography.area_weights(MIROC)
    MOHCCCLM_grid_areas = iris.analysis.cartography.area_weights(MOHCCCLM)
    MOHCKNMI_grid_areas = iris.analysis.cartography.area_weights(MOHCKNMI)
    MOHCSMHI_grid_areas = iris.analysis.cartography.area_weights(MOHCSMHI)
    MPICCLM_grid_areas = iris.analysis.cartography.area_weights(MPICCLM)
    MPIREMO_grid_areas = iris.analysis.cartography.area_weights(MPIREMO)
    MPISMHI_grid_areas = iris.analysis.cartography.area_weights(MPISMHI)
    NCCDMI_grid_areas = iris.analysis.cartography.area_weights(NCCDMI)
    NCCSMHI_grid_areas = iris.analysis.cartography.area_weights(NCCSMHI)
    NOAA_grid_areas = iris.analysis.cartography.area_weights(NOAA)
    
    CanESM2_grid_areas = iris.analysis.cartography.area_weights(CanESM2)
    CNRMG_grid_areas = iris.analysis.cartography.area_weights(CNRMG)
    MK3_grid_areas = iris.analysis.cartography.area_weights(MK3)
    EARTH_grid_areas = iris.analysis.cartography.area_weights(EARTH)
    EARTH3_grid_areas = iris.analysis.cartography.area_weights(EARTH3)
    GFDL_grid_areas = iris.analysis.cartography.area_weights(GFDL)
    HadGEM2_grid_areas = iris.analysis.cartography.area_weights(HadGEM2)
    IPSLG_grid_areas = iris.analysis.cartography.area_weights(IPSLG)
    MIROCG_grid_areas = iris.analysis.cartography.area_weights(MIROCG)
    MPI_grid_areas = iris.analysis.cartography.area_weights(MPI)
    NorESM1_grid_areas = iris.analysis.cartography.area_weights(NorESM1)
    
    CRU_grid_areas = iris.analysis.cartography.area_weights(CRU)
    UDel_grid_areas = iris.analysis.cartography.area_weights(UDel)
    
    #We want to plot the mean for the whole region so we need a mean of all the lats and lons
    CCCma_mean = CCCma.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CCCma_grid_areas) 
    CLMcom_mean = CLMcom.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=CLMcom_grid_areas)
    DMI_mean = DMI.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=DMI_grid_areas)     
    KNMI_mean = KNMI.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=KNMI_grid_areas)
    MPIE_mean = MPIE.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIE_grid_areas)
    SMHI_mean = SMHI.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=SMHI_grid_areas)
    
    CCCmaCanRCM_mean = CCCmaCanRCM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CCCmaCanRCM_grid_areas)      
    CCCmaSMHI_mean = CCCmaSMHI.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=CCCmaSMHI_grid_areas)
    CNRM_mean = CNRM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRM_grid_areas)                                               
    CNRMSMHI_mean = CNRMSMHI.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMSMHI_grid_areas)  
    CSIRO_mean = CSIRO.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CSIRO_grid_areas)
    ICHECDMI_mean = ICHECDMI.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECDMI_grid_areas) 
    ICHECCCLM_mean = ICHECCCLM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECCCLM_grid_areas)
    ICHECKNMI_mean = ICHECKNMI.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECKNMI_grid_areas)
    ICHECMPI_mean = ICHECMPI.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECMPI_grid_areas)
    ICHECSMHI_mean = ICHECSMHI.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECSMHI_grid_areas)
    IPSL_mean = IPSL.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=IPSL_grid_areas)
    MIROC_mean = MIROC.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MIROC_grid_areas)
    MOHCCCLM_mean = MOHCCCLM.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCCCLM_grid_areas)
    MOHCKNMI_mean = MOHCKNMI.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCKNMI_grid_areas)
    MOHCSMHI_mean = MOHCSMHI.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCSMHI_grid_areas)
    MPICCLM_mean = MPICCLM.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=MPICCLM_grid_areas)                                              
    MPIREMO_mean = MPIREMO.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIREMO_grid_areas)                                              
    MPISMHI_mean = MPISMHI.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPISMHI_grid_areas)
    NCCDMI_mean = NCCDMI.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NCCDMI_grid_areas)
    NCCSMHI_mean = NCCSMHI.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NCCSMHI_grid_areas) 
    NOAA_mean = NOAA.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NOAA_grid_areas)
   
    CanESM2_mean = CanESM2.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CanESM2_grid_areas)                                               
    CNRMG_mean = CNRMG.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=CNRMG_grid_areas)
    MK3_mean = MK3.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MK3_grid_areas) 
    EARTH_mean = EARTH.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=EARTH_grid_areas)   
    EARTH3_mean = EARTH3.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=EARTH3_grid_areas)   
    GFDL_mean = GFDL.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=GFDL_grid_areas)  
    HadGEM2_mean = HadGEM2.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=HadGEM2_grid_areas)
    IPSLG_mean = IPSLG.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=IPSLG_grid_areas) 
    MIROCG_mean = MIROCG.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MIROCG_grid_areas)
    MPI_mean = MPI.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPI_grid_areas)
    NorESM1_mean = NorESM1.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NorESM1_grid_areas)
    
    CRU_mean = CRU.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CRU_grid_areas) 
    UDel_mean = UDel.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=UDel_grid_areas) 
    
    #Create averages      
    AverageEY = (CCCma_mean.data + CLMcom_mean.data + DMI_mean.data + KNMI_mean.data + MPIE_mean.data + SMHI_mean.data)/6.
    
    AverageRY = (CCCmaCanRCM_mean.data + CCCmaSMHI_mean.data + CNRM_mean.data + CNRMSMHI_mean.data + CSIRO_mean.data + ICHECDMI_mean.data + ICHECCCLM_mean.data + ICHECKNMI_mean.data + ICHECMPI_mean.data + ICHECSMHI_mean.data + IPSL_mean.data + MIROC_mean.data + MOHCCCLM_mean.data + MOHCKNMI_mean.data + MOHCSMHI_mean.data + MPICCLM_mean.data + MPIREMO_mean.data + MPISMHI_mean.data + NCCDMI_mean.data + NCCSMHI_mean.data + NOAA_mean.data)/21.
    
    AverageGY = (CanESM2_mean.data + CNRMG_mean.data + MK3_mean.data + EARTH_mean.data + EARTH3_mean.data + GFDL_mean.data + HadGEM2_mean.data + IPSLG_mean.data + MIROCG_mean.data + MPI_mean.data + NorESM1_mean.data )/11. 
    
    ObsY = (CRU_mean.data + UDel_mean.data)/2.
    
    X = np.arange(1,13,1) 
    
    #-------------------------------------------------------------------------
    
    #PART 4: PLOT LINE GRAPH 
    #set x-axis ticks                                                                                            
    plt.xticks(range(12), calendar.month_abbr[0:12])
    
    #assign the line colours and set x axis to 'month' rather than 'time'        
    plt.plot(X, ObsY, label='Observed', lw=3, color='black')
    plt.plot(X, AverageRY, label='Average RCM', lw=1.5, color='cyan')
    plt.plot(X, AverageGY, label='Average GCM', lw=1.5, color='magenta') 
    plt.plot(X, AverageEY, label='Average ERAINT', lw=1.5, color='blue') 
        
    #set a title for the y axis
    plt.ylabel('Near-Surface Temperature (degrees Celsius)')
    
    #create a legend and set its location to under the graph
    plt.legend(loc="upper center", bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=2)
    
    #create a title
    plt.title('Tas for Malawi by Month 1990-2005', fontsize=11)   
    
    #add grid lines
    plt.grid()
    
    #save the image of the graph and include full legend
    plt.savefig('Tas_Ave_LineGraph_Monthly3', bbox_inches='tight')

    #show the graph in the console
    iplt.show()
    
    print ObsY
    print AverageRY
    print AverageGY
    print AverageEY
if __name__ == '__main__':
    main() 
    