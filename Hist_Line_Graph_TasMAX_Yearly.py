"""
Created on Wednesday August 2 2017

@author: s0899345
"""

import matplotlib.pyplot as plt
import iris
import iris.coord_categorisation as iriscc
import iris.plot as iplt
import iris.quickplot as qplt
import iris.analysis.cartography
import numpy as np
from cf_units import Unit

#this file is split into parts as follows:
    #PART 1: load and format all models 
    #PART 2: load and format observed data
    #PART 3: format files to be plot specific
    #PART 4: plot data
    
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
    
    XE = np.array([1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008])
    X = np.array([1961,1962,1963,1964,1965,1966,1967,1968,1969,1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005])  
    
    #-------------------------------------------------------------------------
    #PART 4: PLOT LINE GRAPH
    #PART 4A: ERAINT Models Line Graph
    #limit x axis    
    plt.xlim((1990,2008)) 
    
    #assign the line colours and set x axis to 'year' rather than 'time'
    qplt.plot(CCCmaSON_mean.coord('year'), CCCmaSON_mean, label='CanRCM4_ERAINT', lw=1.5, color='blue')
    qplt.plot(CLMcomSON_mean.coord('year'), CLMcomSON_mean, label='CCLM4-8-17_ERAINT', lw=1.5, color='green')
    qplt.plot(DMISON_mean.coord('year'), DMISON_mean, label='HIRHAM5_ERAINT', lw=1.5, color='red')
    qplt.plot(KNMISON_mean.coord('year'), KNMISON_mean, label='RACMO22T_ERAINT', lw=1.5, color='cyan')
    qplt.plot(MPIESON_mean.coord('year'), MPIESON_mean, label='REMO2009_ERAINT', lw=1.5, color='magenta')
    qplt.plot(SMHISON_mean.coord('year'), SMHISON_mean, label='RCA4_ERAINT', lw=1.5, color='yellow')
    plt.plot(XE, ObsESONY, label='Observed', lw=3, color='black')
    plt.plot(XE, AverageSONEY, label='Average ERAINT', lw=3, color='grey')
        
    #set a title for the y axis
    plt.ylabel('Near-Surface Temperature (degrees Celsius)')
               
    #create a legend and set its location to under the graph
    plt.legend(loc="upper center", bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=2)
    
    #create a title
    plt.title('TasMax for Malawi SON 1990-2008', fontsize=11)   
    
    #add grid lines
    plt.grid()
    
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_ERAINT_LineGraph_SON', bbox_inches='tight')

    #show the graph in the console
    iplt.show()
    
    #limit x axis    
    plt.xlim((1990,2008)) 
    
    #assign the line colours and set x axis to 'year' rather than 'time'
    qplt.plot(CCCmaDJF_mean.coord('year'), CCCmaDJF_mean, label='CanRCM4_ERAINT', lw=1.5, color='blue')
    qplt.plot(CLMcomDJF_mean.coord('year'), CLMcomDJF_mean, label='CCLM4-8-17_ERAINT', lw=1.5, color='green')
    qplt.plot(DMIDJF_mean.coord('year'), DMIDJF_mean, label='HIRHAM5_ERAINT', lw=1.5, color='red')
    qplt.plot(KNMIDJF_mean.coord('year'), KNMIDJF_mean, label='RACMO22T_ERAINT', lw=1.5, color='cyan')
    qplt.plot(MPIEDJF_mean.coord('year'), MPIEDJF_mean, label='REMO2009_ERAINT', lw=1.5, color='magenta')
    qplt.plot(SMHIDJF_mean.coord('year'), SMHIDJF_mean, label='RCA4_ERAINT', lw=1.5, color='yellow')
    plt.plot(XE, ObsEDJFY, label='Observed', lw=3, color='black')
    plt.plot(XE, AverageDJFEY, label='Average ERAINT', lw=3, color='grey')
       
    #set a title for the y axis
    plt.ylabel('Near-Surface Temperature (degrees Celsius)')
               
    #create a legend and set its location to under the graph
    plt.legend(loc="upper center", bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=2)
    
    #create a title
    plt.title('TasMax for Malaw DJF 1990-2008', fontsize=11)   
    
    #add grid lines
    plt.grid()
    
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_ERAINT_LineGraph_DJF', bbox_inches='tight')

    #show the graph in the console
    iplt.show()
    
    #limit x axis    
    plt.xlim((1990,2008)) 
    
    #assign the line colours and set x axis to 'year' rather than 'time'
    qplt.plot(CCCmaMAM_mean.coord('year'), CCCmaMAM_mean, label='CanRCM4_ERAINT', lw=1.5, color='blue')
    qplt.plot(CLMcomMAM_mean.coord('year'), CLMcomMAM_mean, label='CCLM4-8-17_ERAINT', lw=1.5, color='green')
    qplt.plot(DMIMAM_mean.coord('year'), DMIMAM_mean, label='HIRHAM5_ERAINT', lw=1.5, color='red')
    qplt.plot(KNMIMAM_mean.coord('year'), KNMIMAM_mean, label='RACMO22T_ERAINT', lw=1.5, color='cyan')
    qplt.plot(MPIEMAM_mean.coord('year'), MPIEMAM_mean, label='REMO2009_ERAINT', lw=1.5, color='magenta')
    qplt.plot(SMHIMAM_mean.coord('year'), SMHIMAM_mean, label='RCA4_ERAINT', lw=1.5, color='yellow')
    plt.plot(XE, ObsEMAMY, label='Observed', lw=3, color='black')
    plt.plot(XE, AverageMAMEY, label='Average ERAINT', lw=3, color='grey')
        
    #set a title for the y axis
    plt.ylabel('Near-Surface Temperature (degrees Celsius)')   
               
    #create a legend and set its location to under the graph
    plt.legend(loc="upper center", bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=2)
    
    #create a title
    plt.title('TasMax for Malaw MAM 1990-2008', fontsize=11)   
    
    #add grid lines
    plt.grid()
    
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_ERAINT_LineGraph_MAM', bbox_inches='tight')

    #show the graph in the console
    iplt.show()
    
    #limit x axis    
    plt.xlim((1990,2008)) 
    
    #assign the line colours and set x axis to 'year' rather than 'time'
    qplt.plot(CCCmaJJA_mean.coord('year'), CCCmaJJA_mean, label='CanRCM4_ERAINT', lw=1.5, color='blue')
    qplt.plot(CLMcomJJA_mean.coord('year'), CLMcomJJA_mean, label='CCLM4-8-17_ERAINT', lw=1.5, color='green')
    qplt.plot(DMIJJA_mean.coord('year'), DMIJJA_mean, label='HIRHAM5_ERAINT', lw=1.5, color='red')
    qplt.plot(KNMIJJA_mean.coord('year'), KNMIJJA_mean, label='RACMO22T_ERAINT', lw=1.5, color='cyan')
    qplt.plot(MPIEJJA_mean.coord('year'), MPIEJJA_mean, label='REMO2009_ERAINT', lw=1.5, color='magenta')
    qplt.plot(SMHIJJA_mean.coord('year'), SMHIJJA_mean, label='RCA4_ERAINT', lw=1.5, color='yellow')
    plt.plot(XE, ObsEJJAY, label='Observed', lw=3, color='black')
    plt.plot(XE, AverageJJAEY, label='Average ERAINT', lw=3, color='grey')
        
    #set a title for the y axis
    plt.ylabel('Near-Surface Temperature (degrees Celsius)')
               
    #create a legend and set its location to under the graph
    plt.legend(loc="upper center", bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=2)
    
    #create a title
    plt.title('TasMax for Malaw JJA 1990-2008', fontsize=11)   
    
    #add grid lines
    plt.grid()
    
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_ERAINT_LineGraph_JJA', bbox_inches='tight')

    #show the graph in the console
    iplt.show()
    
    #limit x axis    
    plt.xlim((1990,2008)) 
    
    #assign the line colours and set x axis to 'year' rather than 'time'
    qplt.plot(CCCmaYR_mean.coord('year'),CCCmaYR_mean, label='CanRCM4_ERAINT', lw=1.5, color='blue')
    qplt.plot(CLMcomYR_mean.coord('year'), CLMcomYR_mean, label='CCLM4-8-17_ERAINT', lw=1.5, color='green')
    qplt.plot(DMIYR_mean.coord('year'), DMIYR_mean, label='HIRHAM5_ERAINT', lw=1.5, color='red')
    qplt.plot(KNMIYR_mean.coord('year'), KNMIYR_mean, label='RACMO22T_ERAINT', lw=1.5, color='cyan')
    qplt.plot(MPIEYR_mean.coord('year'), MPIEYR_mean, label='REMO2009_ERAINT', lw=1.5, color='magenta')
    qplt.plot(SMHIYR_mean.coord('year'), SMHIYR_mean, label='RCA4_ERAINT', lw=1.5, color='yellow')
    plt.plot(XE, ObsEY, label='Observed', lw=3, color='black')
    plt.plot(XE, AverageEY, label='Average ERAINT', lw=3, color='grey')
        
    #set a title for the y axis
    plt.ylabel('Near-Surface Temperature (degrees Celsius)')
    
    #create a legend and set its location to under the graph
    plt.legend(loc="upper center", bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=2)
    
    #create a title
    plt.title('TasMax for Malawi 1990-2008', fontsize=11)   
    
    #add grid lines
    plt.grid()
    
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_ERAINT_LineGraph_Annual', bbox_inches='tight')

    #show the graph in the console
    iplt.show()   
    
    #PART 4B: Regional Climate Models Line Graph
    #limit x axis    
    plt.xlim((1961,2005)) 
    
    #assign the line colours and set x axis to 'year' rather than 'time'
    qplt.plot(CCCmaCanRCMSON_mean.coord('year'), CCCmaCanRCMSON_mean, label='CanRCM4_CanESM2', lw=1.5, color='blue')
    qplt.plot(CCCmaSMHISON_mean.coord('year'), CCCmaSMHISON_mean, label='RCA4_CanESM2', lw=1.5, color='green')
    qplt.plot(CNRMSON_mean.coord('year'), CNRMSON_mean, label='CCLM4-8-17_CNRM-CM5', lw=1.5, color='red')
    qplt.plot(CNRMSMHISON_mean.coord('year'), CNRMSMHISON_mean, label='RCA4_CNRM-CM5', lw=1.5, color='cyan')
    qplt.plot(CSIROSON_mean.coord('year'), CSIROSON_mean, label='RCA4_CSIRO-MK3', lw=1.5, color='magenta')
    qplt.plot(ICHECDMISON_mean.coord('year'), ICHECDMISON_mean, label='HIRHAM5_EC-EARTH', lw=1.5, color='yellow')
    qplt.plot(ICHECCCLMSON_mean.coord('year'), ICHECCCLMSON_mean, label='CCLM4-8-17_EC-EARTH ', lw=1.5, color='blue', linestyle='--')
    qplt.plot(ICHECKNMISON_mean.coord('year'), ICHECKNMISON_mean, label='RACMO22T_EC-EARTH', lw=1.5, color='green', linestyle='--')
    qplt.plot(ICHECMPISON_mean.coord('year'), ICHECMPISON_mean, label='REMO2009_EC-EARTH ', lw=1.5, color='red', linestyle='--') 
    qplt.plot(ICHECSMHISON_mean.coord('year'), ICHECSMHISON_mean, label='RCA4_EC-EARTH', lw=1.5, color='cyan', linestyle='--') 
    qplt.plot(IPSLSON_mean.coord('year'), IPSLSON_mean, label='RCA4_IPSL-CM5A-MR', lw=1.5, color='magenta', linestyle='--') 
    qplt.plot(MIROCSON_mean.coord('year'), MIROCSON_mean, label='RCA4_MIROC5', lw=1.5, color='yellow', linestyle='--')
    qplt.plot(MOHCCCLMSON_mean.coord('year'), MOHCCCLMSON_mean, label='CCLM4-8-17_HadGEM2-ES', lw=1.5, color='blue', linestyle='-.')
    qplt.plot(MOHCKNMISON_mean.coord('year'), MOHCKNMISON_mean, label='RACMO22T_HadGEM2-ES', lw=1.5, color='green', linestyle='-.')
    qplt.plot(MOHCSMHISON_mean.coord('year'), MOHCSMHISON_mean, label='RCA4_HadGEM2-ES', lw=1.5, color='red', linestyle='-.')
    qplt.plot(MPICCLMSON_mean.coord('year'), MPICCLMSON_mean, label='CCLM4-8-17_MPI-ESM-LR', lw=1.5, color='cyan', linestyle='-.')
    qplt.plot(MPIREMOSON_mean.coord('year'), MPIREMOSON_mean, label='REMO2009_MPI-ESM-LR ', lw=1.5, color='magenta', linestyle='-.')
    qplt.plot(MPISMHISON_mean.coord('year'), MPISMHISON_mean, label='RCA4_MPI-ESM-LR ', lw=1.5, color='yellow', linestyle='-.')
    qplt.plot(NCCDMISON_mean.coord('year'), NCCDMISON_mean, label='HIRHAM5_NorESM1-M', lw=1.5, color='blue', linestyle=':')
    qplt.plot(NCCSMHISON_mean.coord('year'), NCCSMHISON_mean, label='RCA4_NorESM1-M', lw=1.5, color='green', linestyle=':')
    qplt.plot(NOAASON_mean.coord('year'), NOAASON_mean, label='RCA4_GFDL-ESM2M', lw=1.5, color='red', linestyle=':') 
    plt.plot(X, ObsSONY, label='Observed', lw=3, color='black')
    plt.plot(X, AverageSONRY, label='Average RCM', lw=3, color='grey')
            
    #set a title for the y axis
    plt.ylabel('Near-Surface Temperature (degrees Celsius)')
    
    #create a legend and set its location to under the graph
    plt.legend(loc="upper center", bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=2)
    
    #create a title
    plt.title('TasMax for Malawi SON 1961-2005', fontsize=11)   
    
    #add grid lines
    plt.grid()
    
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_RCM_LineGraph_SON', bbox_inches='tight')

    #show the graph in the console
    iplt.show() 
    
    #limit x axis    
    plt.xlim((1961,2005)) 
    
    #assign the line colours and set x axis to 'year' rather than 'time'
    qplt.plot(CCCmaCanRCMDJF_mean.coord('year'), CCCmaCanRCMDJF_mean, label='CanRCM4_CanESM2', lw=1.5, color='blue')
    qplt.plot(CCCmaSMHIDJF_mean.coord('year'), CCCmaSMHIDJF_mean, label='RCA4_CanESM2', lw=1.5, color='green')
    qplt.plot(CNRMDJF_mean.coord('year'), CNRMDJF_mean, label='CCLM4-8-17_CNRM-CM5', lw=1.5, color='red')
    qplt.plot(CNRMSMHIDJF_mean.coord('year'), CNRMSMHIDJF_mean, label='RCA4_CNRM-CM5', lw=1.5, color='cyan')
    qplt.plot(CSIRODJF_mean.coord('year'), CSIRODJF_mean, label='RCA4_CSIRO-MK3', lw=1.5, color='magenta')
    qplt.plot(ICHECDMIDJF_mean.coord('year'), ICHECDMIDJF_mean, label='HIRHAM5_EC-EARTH', lw=1.5, color='yellow')
    qplt.plot(ICHECCCLMDJF_mean.coord('year'), ICHECCCLMDJF_mean, label='CCLM4-8-17_EC-EARTH ', lw=1.5, color='blue', linestyle='--')
    qplt.plot(ICHECKNMIDJF_mean.coord('year'), ICHECKNMIDJF_mean, label='RACMO22T_EC-EARTH', lw=1.5, color='green', linestyle='--')
    qplt.plot(ICHECMPIDJF_mean.coord('year'), ICHECMPIDJF_mean, label='REMO2009_EC-EARTH ', lw=1.5, color='red', linestyle='--') 
    qplt.plot(ICHECSMHIDJF_mean.coord('year'), ICHECSMHIDJF_mean, label='RCA4_EC-EARTH', lw=1.5, color='cyan', linestyle='--') 
    qplt.plot(IPSLDJF_mean.coord('year'), IPSLDJF_mean, label='RCA4_IPSL-CM5A-MR', lw=1.5, color='magenta', linestyle='--') 
    qplt.plot(MIROCDJF_mean.coord('year'), MIROCDJF_mean, label='RCA4_MIROC5', lw=1.5, color='yellow', linestyle='--')
    qplt.plot(MOHCCCLMDJF_mean.coord('year'), MOHCCCLMDJF_mean, label='CCLM4-8-17_HadGEM2-ES', lw=1.5, color='blue', linestyle='-.')
    qplt.plot(MOHCKNMIDJF_mean.coord('year'), MOHCKNMIDJF_mean, label='RACMO22T_HadGEM2-ES', lw=1.5, color='green', linestyle='-.')
    qplt.plot(MOHCSMHIDJF_mean.coord('year'), MOHCSMHIDJF_mean, label='RCA4_HadGEM2-ES', lw=1.5, color='red', linestyle='-.')
    qplt.plot(MPICCLMDJF_mean.coord('year'), MPICCLMDJF_mean, label='CCLM4-8-17_MPI-ESM-LR', lw=1.5, color='cyan', linestyle='-.')
    qplt.plot(MPIREMODJF_mean.coord('year'), MPIREMODJF_mean, label='REMO2009_MPI-ESM-LR ', lw=1.5, color='magenta', linestyle='-.')
    qplt.plot(MPISMHIDJF_mean.coord('year'), MPISMHIDJF_mean, label='RCA4_MPI-ESM-LR ', lw=1.5, color='yellow', linestyle='-.')
    qplt.plot(NCCDMIDJF_mean.coord('year'), NCCDMIDJF_mean, label='HIRHAM5_NorESM1-M', lw=1.5, color='blue', linestyle=':')
    qplt.plot(NCCSMHIDJF_mean.coord('year'), NCCSMHIDJF_mean, label='RCA4_NorESM1-M', lw=1.5, color='green', linestyle=':')
    qplt.plot(NOAADJF_mean.coord('year'), NOAADJF_mean, label='RCA4_GFDL-ESM2M', lw=1.5, color='red', linestyle=':') 
    plt.plot(X, ObsDJFY, label='Observed', lw=3, color='black')
    plt.plot(X, AverageDJFRY, label='Average RCM', lw=3, color='grey')
            
    #set a title for the y axis
    plt.ylabel('Near-Surface Temperature (degrees Celsius)')
    
    #create a legend and set its location to under the graph
    plt.legend(loc="upper center", bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=2)
    
    #create a title
    plt.title('TasMax for Malawi DJF 1961-2005', fontsize=11)   
    
    #add grid lines
    plt.grid()
    
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_RCM_LineGraph_DJF', bbox_inches='tight')

    #show the graph in the console
    iplt.show() 
    
    #limit x axis    
    plt.xlim((1961,2005)) 
    
    #assign the line colours and set x axis to 'year' rather than 'time'
    qplt.plot(CCCmaCanRCMMAM_mean.coord('year'), CCCmaCanRCMMAM_mean, label='CanRCM4_CanESM2', lw=1.5, color='blue')
    qplt.plot(CCCmaSMHIMAM_mean.coord('year'), CCCmaSMHIMAM_mean, label='RCA4_CanESM2', lw=1.5, color='green')
    qplt.plot(CNRMMAM_mean.coord('year'), CNRMMAM_mean, label='CCLM4-8-17_CNRM-CM5', lw=1.5, color='red')
    qplt.plot(CNRMSMHIMAM_mean.coord('year'), CNRMSMHIMAM_mean, label='RCA4_CNRM-CM5', lw=1.5, color='cyan')
    qplt.plot(CSIROMAM_mean.coord('year'), CSIROMAM_mean, label='RCA4_CSIRO-MK3', lw=1.5, color='magenta')
    qplt.plot(ICHECDMIMAM_mean.coord('year'), ICHECDMIMAM_mean, label='HIRHAM5_EC-EARTH', lw=1.5, color='yellow')
    qplt.plot(ICHECCCLMMAM_mean.coord('year'), ICHECCCLMMAM_mean, label='CCLM4-8-17_EC-EARTH ', lw=1.5, color='blue', linestyle='--')
    qplt.plot(ICHECKNMIMAM_mean.coord('year'), ICHECKNMIMAM_mean, label='RACMO22T_EC-EARTH', lw=1.5, color='green', linestyle='--')
    qplt.plot(ICHECMPIMAM_mean.coord('year'), ICHECMPIMAM_mean, label='REMO2009_EC-EARTH ', lw=1.5, color='red', linestyle='--') 
    qplt.plot(ICHECSMHIMAM_mean.coord('year'), ICHECSMHIMAM_mean, label='RCA4_EC-EARTH', lw=1.5, color='cyan', linestyle='--') 
    qplt.plot(IPSLMAM_mean.coord('year'), IPSLMAM_mean, label='RCA4_IPSL-CM5A-MR', lw=1.5, color='magenta', linestyle='--') 
    qplt.plot(MIROCMAM_mean.coord('year'), MIROCMAM_mean, label='RCA4_MIROC5', lw=1.5, color='yellow', linestyle='--')
    qplt.plot(MOHCCCLMMAM_mean.coord('year'), MOHCCCLMMAM_mean, label='CCLM4-8-17_HadGEM2-ES', lw=1.5, color='blue', linestyle='-.')
    qplt.plot(MOHCKNMIMAM_mean.coord('year'), MOHCKNMIMAM_mean, label='RACMO22T_HadGEM2-ES', lw=1.5, color='green', linestyle='-.')
    qplt.plot(MOHCSMHIMAM_mean.coord('year'), MOHCSMHIMAM_mean, label='RCA4_HadGEM2-ES', lw=1.5, color='red', linestyle='-.')
    qplt.plot(MPICCLMMAM_mean.coord('year'), MPICCLMMAM_mean, label='CCLM4-8-17_MPI-ESM-LR', lw=1.5, color='cyan', linestyle='-.')
    qplt.plot(MPIREMOMAM_mean.coord('year'), MPIREMOMAM_mean, label='REMO2009_MPI-ESM-LR ', lw=1.5, color='magenta', linestyle='-.')
    qplt.plot(MPISMHIMAM_mean.coord('year'), MPISMHIMAM_mean, label='RCA4_MPI-ESM-LR ', lw=1.5, color='yellow', linestyle='-.')
    qplt.plot(NCCDMIMAM_mean.coord('year'), NCCDMIMAM_mean, label='HIRHAM5_NorESM1-M', lw=1.5, color='blue', linestyle=':')
    qplt.plot(NCCSMHIMAM_mean.coord('year'), NCCSMHIMAM_mean, label='RCA4_NorESM1-M', lw=1.5, color='green', linestyle=':')
    qplt.plot(NOAAMAM_mean.coord('year'), NOAAMAM_mean, label='RCA4_GFDL-ESM2M', lw=1.5, color='red', linestyle=':') 
    plt.plot(X, ObsMAMY, label='Observed', lw=3, color='black')
    plt.plot(X, AverageMAMRY, label='Average RCM', lw=3, color='grey')
            
    #set a title for the y axis
    plt.ylabel('Near-Surface Temperature (degrees Celsius)')
    
    #create a legend and set its location to under the graph
    plt.legend(loc="upper center", bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=2)
    
    #create a title
    plt.title('TasMax for Malawi MAM 1961-2005', fontsize=11)   
    
    #add grid lines
    plt.grid()
    
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_RCM_LineGraph_MAM', bbox_inches='tight')

    #show the graph in the console
    iplt.show()
   
    #limit x axis    
    plt.xlim((1961,2005)) 
    
    #assign the line colours and set x axis to 'year' rather than 'time'
    qplt.plot(CCCmaCanRCMJJA_mean.coord('year'), CCCmaCanRCMJJA_mean, label='CanRCM4_CanESM2', lw=1.5, color='blue')
    qplt.plot(CCCmaSMHIJJA_mean.coord('year'), CCCmaSMHIJJA_mean, label='RCA4_CanESM2', lw=1.5, color='green')
    qplt.plot(CNRMJJA_mean.coord('year'), CNRMJJA_mean, label='CCLM4-8-17_CNRM-CM5', lw=1.5, color='red')
    qplt.plot(CNRMSMHIJJA_mean.coord('year'), CNRMSMHIJJA_mean, label='RCA4_CNRM-CM5', lw=1.5, color='cyan')
    qplt.plot(CSIROJJA_mean.coord('year'), CSIROJJA_mean, label='RCA4_CSIRO-MK3', lw=1.5, color='magenta')
    qplt.plot(ICHECDMIJJA_mean.coord('year'), ICHECDMIJJA_mean, label='HIRHAM5_EC-EARTH', lw=1.5, color='yellow')
    qplt.plot(ICHECCCLMJJA_mean.coord('year'), ICHECCCLMJJA_mean, label='CCLM4-8-17_EC-EARTH ', lw=1.5, color='blue', linestyle='--')
    qplt.plot(ICHECKNMIJJA_mean.coord('year'), ICHECKNMIJJA_mean, label='RACMO22T_EC-EARTH', lw=1.5, color='green', linestyle='--')
    qplt.plot(ICHECMPIJJA_mean.coord('year'), ICHECMPIJJA_mean, label='REMO2009_EC-EARTH ', lw=1.5, color='red', linestyle='--') 
    qplt.plot(ICHECSMHIJJA_mean.coord('year'), ICHECSMHIJJA_mean, label='RCA4_EC-EARTH', lw=1.5, color='cyan', linestyle='--') 
    qplt.plot(IPSLJJA_mean.coord('year'), IPSLJJA_mean, label='RCA4_IPSL-CM5A-MR', lw=1.5, color='magenta', linestyle='--') 
    qplt.plot(MIROCJJA_mean.coord('year'), MIROCJJA_mean, label='RCA4_MIROC5', lw=1.5, color='yellow', linestyle='--')
    qplt.plot(MOHCCCLMJJA_mean.coord('year'), MOHCCCLMJJA_mean, label='CCLM4-8-17_HadGEM2-ES', lw=1.5, color='blue', linestyle='-.')
    qplt.plot(MOHCKNMIJJA_mean.coord('year'), MOHCKNMIJJA_mean, label='RACMO22T_HadGEM2-ES', lw=1.5, color='green', linestyle='-.')
    qplt.plot(MOHCSMHIJJA_mean.coord('year'), MOHCSMHIJJA_mean, label='RCA4_HadGEM2-ES', lw=1.5, color='red', linestyle='-.')
    qplt.plot(MPICCLMJJA_mean.coord('year'), MPICCLMJJA_mean, label='CCLM4-8-17_MPI-ESM-LR', lw=1.5, color='cyan', linestyle='-.')
    qplt.plot(MPIREMOJJA_mean.coord('year'), MPIREMOJJA_mean, label='REMO2009_MPI-ESM-LR ', lw=1.5, color='magenta', linestyle='-.')
    qplt.plot(MPISMHIJJA_mean.coord('year'), MPISMHIJJA_mean, label='RCA4_MPI-ESM-LR ', lw=1.5, color='yellow', linestyle='-.')
    qplt.plot(NCCDMIJJA_mean.coord('year'), NCCDMIJJA_mean, label='HIRHAM5_NorESM1-M', lw=1.5, color='blue', linestyle=':')
    qplt.plot(NCCSMHIJJA_mean.coord('year'), NCCSMHIJJA_mean, label='RCA4_NorESM1-M', lw=1.5, color='green', linestyle=':')
    qplt.plot(NOAAJJA_mean.coord('year'), NOAAJJA_mean, label='RCA4_GFDL-ESM2M', lw=1.5, color='red', linestyle=':') 
    plt.plot(X, ObsJJAY, label='Observed', lw=3, color='black')
    plt.plot(X, AverageJJARY, label='Average RCM', lw=3, color='grey')
            
    #set a title for the y axis
    plt.ylabel('Near-Surface Temperature (degrees Celsius)')
    
    #create a legend and set its location to under the graph
    plt.legend(loc="upper center", bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=2)
    
    #create a title
    plt.title('TasMax for Malawi JJA 1961-2005', fontsize=11)   
    
    #add grid lines
    plt.grid()
    
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_RCM_LineGraph_JJA', bbox_inches='tight')

    #show the graph in the console
    iplt.show()
    
    #limit x axis    
    plt.xlim((1961,2005)) 
    
    #assign the line colours and set x axis to 'year' rather than 'time'
    qplt.plot(CCCmaCanRCMYR_mean.coord('year'), CCCmaCanRCMYR_mean, label='CanRCM4_CanESM2', lw=1.5, color='blue')
    qplt.plot(CCCmaSMHIYR_mean.coord('year'), CCCmaSMHIYR_mean, label='RCA4_CanESM2', lw=1.5, color='green')
    qplt.plot(CNRMYR_mean.coord('year'), CNRMYR_mean, label='CCLM4-8-17_CNRM-CM5', lw=1.5, color='red')
    qplt.plot(CNRMSMHIYR_mean.coord('year'), CNRMSMHIYR_mean, label='RCA4_CNRM-CM5', lw=1.5, color='cyan')
    qplt.plot(CSIROYR_mean.coord('year'), CSIROYR_mean, label='RCA4_CSIRO-MK3', lw=1.5, color='magenta')
    qplt.plot(ICHECDMIYR_mean.coord('year'), ICHECDMIYR_mean, label='HIRHAM5_EC-EARTH', lw=1.5, color='yellow')
    qplt.plot(ICHECCCLMYR_mean.coord('year'), ICHECCCLMYR_mean, label='CCLM4-8-17_EC-EARTH ', lw=1.5, color='blue', linestyle='--')
    qplt.plot(ICHECKNMIYR_mean.coord('year'), ICHECKNMIYR_mean, label='RACMO22T_EC-EARTH', lw=1.5, color='green', linestyle='--')
    qplt.plot(ICHECMPIYR_mean.coord('year'), ICHECMPIYR_mean, label='REMO2009_EC-EARTH ', lw=1.5, color='red', linestyle='--') 
    qplt.plot(ICHECSMHIYR_mean.coord('year'), ICHECSMHIYR_mean, label='RCA4_EC-EARTH', lw=1.5, color='cyan', linestyle='--') 
    qplt.plot(IPSLYR_mean.coord('year'), IPSLYR_mean, label='RCA4_IPSL-CM5A-MR', lw=1.5, color='magenta', linestyle='--') 
    qplt.plot(MIROCYR_mean.coord('year'), MIROCYR_mean, label='RCA4_MIROC5', lw=1.5, color='yellow', linestyle='--')
    qplt.plot(MOHCCCLMYR_mean.coord('year'), MOHCCCLMYR_mean, label='CCLM4-8-17_HadGEM2-ES', lw=1.5, color='blue', linestyle='-.')
    qplt.plot(MOHCKNMIYR_mean.coord('year'), MOHCKNMIYR_mean, label='RACMO22T_HadGEM2-ES', lw=1.5, color='green', linestyle='-.')
    qplt.plot(MOHCSMHIYR_mean.coord('year'), MOHCSMHIYR_mean, label='RCA4_HadGEM2-ES', lw=1.5, color='red', linestyle='-.')
    qplt.plot(MPICCLMYR_mean.coord('year'), MPICCLMYR_mean, label='CCLM4-8-17_MPI-ESM-LR', lw=1.5, color='cyan', linestyle='-.')
    qplt.plot(MPIREMOYR_mean.coord('year'), MPIREMOYR_mean, label='REMO2009_MPI-ESM-LR ', lw=1.5, color='magenta', linestyle='-.')
    qplt.plot(MPISMHIYR_mean.coord('year'), MPISMHIYR_mean, label='RCA4_MPI-ESM-LR ', lw=1.5, color='yellow', linestyle='-.')
    qplt.plot(NCCDMIYR_mean.coord('year'), NCCDMIYR_mean, label='HIRHAM5_NorESM1-M', lw=1.5, color='blue', linestyle=':')
    qplt.plot(NCCSMHIYR_mean.coord('year'), NCCSMHIYR_mean, label='RCA4_NorESM1-M', lw=1.5, color='green', linestyle=':')
    qplt.plot(NOAAYR_mean.coord('year'), NOAAYR_mean, label='RCA4_GFDL-ESM2M', lw=1.5, color='red', linestyle=':') 
    plt.plot(X, ObsY, label='Observed', lw=3, color='black')
    plt.plot(X, AverageRY, label='Average RCM', lw=3, color='grey')
            
    #set a title for the y axis
    plt.ylabel('Near-Surface Temperature (degrees Celsius)')
    
    #create a legend and set its location to under the graph
    plt.legend(loc="upper center", bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=2)
    
    #create a title
    plt.title('TasMax for Malawi 1961-2005', fontsize=11)   
    
    #add grid lines
    plt.grid()
    
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_RCM_LineGraph_Annual', bbox_inches='tight')

    #show the graph in the console
    iplt.show()
    
    #PART 4C: Global Climate Models Line Graph
    #limit x axis    
    plt.xlim((1961,2005)) 
    
    #assign the line colours and set x axis to 'year' rather than 'time'
    qplt.plot(CanESM2SON_mean.coord('year'), CanESM2SON_mean, label='CanESM2', lw=1.5, color='blue')
    qplt.plot(CNRMGSON_mean.coord('year'), CNRMGSON_mean, label='CNRM_CM5', lw=1.5, color='green')
    qplt.plot(MK3SON_mean.coord('year'), MK3SON_mean, label='CSIRO MK3-6-0', lw=1.5, color='red')
    qplt.plot(EARTHSON_mean.coord('year'), EARTHSON_mean, label='EC-EARTH (r12i1p1)', lw=1.5, color='cyan')
    qplt.plot(EARTH3SON_mean.coord('year'), EARTH3SON_mean, label='EC-EARTH (r3i1p1', lw=1.5, color='magenta')
    qplt.plot(GFDLSON_mean.coord('year'), GFDLSON_mean, label='GFDL-ESM2M', lw=1.5, color='yellow')
    qplt.plot(HadGEM2SON_mean.coord('year'), HadGEM2SON_mean, label='HadGEM2-ES', lw=1.5, color='blue', linestyle='--')
    qplt.plot(IPSLGSON_mean.coord('year'), IPSLGSON_mean, label='IPSL-CM5A-MR', lw=1.5, color='green', linestyle='--')
    qplt.plot(MIROCGSON_mean.coord('year'), MIROCGSON_mean, label='MIROC5 ', lw=1.5, color='red', linestyle='--')
    qplt.plot(MPISON_mean.coord('year'), MPISON_mean, label='MPI-ESM-LR', lw=1.5, color='cyan', linestyle='--')
    qplt.plot(NorESM1SON_mean.coord('year'), NorESM1SON_mean, label='NorESM1-M', lw=1.5, color='magenta', linestyle='--')
    plt.plot(X, ObsSONY, label='Observed', lw=3, color='black')
    plt.plot(X, AverageSONGY, label='Average GCM', lw=3, color='grey')
    
    #set a title for the y axis
    plt.ylabel('Near-Surface Temperature (degrees Celsius)')
    
    #create a legend and set its location to under the graph
    plt.legend(loc="upper center", bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=2)
    
    #create a title
    plt.title('TasMax for Malawi SON 1961-2005', fontsize=11)   
    
    #add grid lines
    plt.grid()
    
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_GCM_LineGraph_SON', bbox_inches='tight')

    #show the graph in the console
    iplt.show() 
    
    #limit x axis    
    plt.xlim((1961,2005)) 
    
    #assign the line colours and set x axis to 'year' rather than 'time'
    qplt.plot(CanESM2DJF_mean.coord('year'), CanESM2DJF_mean, label='CanESM2', lw=1.5, color='blue')
    qplt.plot(CNRMGDJF_mean.coord('year'), CNRMGDJF_mean, label='CNRM_CM5', lw=1.5, color='green')
    qplt.plot(MK3DJF_mean.coord('year'), MK3DJF_mean, label='CSIRO MK3-6-0', lw=1.5, color='red')
    qplt.plot(EARTHDJF_mean.coord('year'), EARTHDJF_mean, label='EC-EARTH (r12i1p1)', lw=1.5, color='cyan')
    qplt.plot(EARTH3DJF_mean.coord('year'), EARTH3DJF_mean, label='EC-EARTH (r3i1p1', lw=1.5, color='magenta')
    qplt.plot(GFDLDJF_mean.coord('year'), GFDLDJF_mean, label='GFDL-ESM2M', lw=1.5, color='yellow')
    qplt.plot(HadGEM2DJF_mean.coord('year'), HadGEM2DJF_mean, label='HadGEM2-ES', lw=1.5, color='blue', linestyle='--')
    qplt.plot(IPSLGDJF_mean.coord('year'), IPSLGDJF_mean, label='IPSL-CM5A-MR', lw=1.5, color='green', linestyle='--')
    qplt.plot(MIROCGDJF_mean.coord('year'), MIROCGDJF_mean, label='MIROC5 ', lw=1.5, color='red', linestyle='--')
    qplt.plot(MPIDJF_mean.coord('year'), MPIDJF_mean, label='MPI-ESM-LR', lw=1.5, color='cyan', linestyle='--')
    qplt.plot(NorESM1DJF_mean.coord('year'), NorESM1DJF_mean, label='NorESM1-M', lw=1.5, color='magenta', linestyle='--')
    plt.plot(X, ObsDJFY, label='Observed', lw=3, color='black')
    plt.plot(X, AverageDJFGY, label='Average GCM', lw=3, color='grey')
    
    #set a title for the y axis
    plt.ylabel('Near-Surface Temperature (degrees Celsius)')
    
    #create a legend and set its location to under the graph
    plt.legend(loc="upper center", bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=2)
    
    #create a title
    plt.title('TasMax for Malawi DJF 1961-2005', fontsize=11)   
    
    #add grid lines
    plt.grid()
    
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_GCM_LineGraph_DJF', bbox_inches='tight')

    #show the graph in the console
    iplt.show() 
    
    #limit x axis    
    plt.xlim((1961,2005)) 
    
    #assign the line colours and set x axis to 'year' rather than 'time'
    qplt.plot(CanESM2MAM_mean.coord('year'), CanESM2MAM_mean, label='CanESM2', lw=1.5, color='blue')
    qplt.plot(CNRMGMAM_mean.coord('year'), CNRMGMAM_mean, label='CNRM_CM5', lw=1.5, color='green')
    qplt.plot(MK3MAM_mean.coord('year'), MK3MAM_mean, label='CSIRO MK3-6-0', lw=1.5, color='red')
    qplt.plot(EARTHMAM_mean.coord('year'), EARTHMAM_mean, label='EC-EARTH (r12i1p1)', lw=1.5, color='cyan')
    qplt.plot(EARTH3MAM_mean.coord('year'), EARTH3MAM_mean, label='EC-EARTH (r3i1p1', lw=1.5, color='magenta')
    qplt.plot(GFDLMAM_mean.coord('year'), GFDLMAM_mean, label='GFDL-ESM2M', lw=1.5, color='yellow')
    qplt.plot(HadGEM2MAM_mean.coord('year'), HadGEM2MAM_mean, label='HadGEM2-ES', lw=1.5, color='blue', linestyle='--')
    qplt.plot(IPSLGMAM_mean.coord('year'), IPSLGMAM_mean, label='IPSL-CM5A-MR', lw=1.5, color='green', linestyle='--')
    qplt.plot(MIROCGMAM_mean.coord('year'), MIROCGMAM_mean, label='MIROC5 ', lw=1.5, color='red', linestyle='--')
    qplt.plot(MPIMAM_mean.coord('year'), MPIMAM_mean, label='MPI-ESM-LR', lw=1.5, color='cyan', linestyle='--')
    qplt.plot(NorESM1MAM_mean.coord('year'), NorESM1MAM_mean, label='NorESM1-M', lw=1.5, color='magenta', linestyle='--')
    plt.plot(X, ObsMAMY, label='Observed', lw=3, color='black')
    plt.plot(X, AverageMAMGY, label='Average GCM', lw=3, color='grey')
    
    #set a title for the y axis
    plt.ylabel('Near-Surface Temperature (degrees Celsius)')
    
    #create a legend and set its location to under the graph
    plt.legend(loc="upper center", bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=2)
    
    #create a title
    plt.title('TasMax for Malawi MAM 1961-2005', fontsize=11)   
    
    #add grid lines
    plt.grid()
    
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_GCM_LineGraph_MAM', bbox_inches='tight')

    #show the graph in the console
    iplt.show()
   
    #limit x axis    
    plt.xlim((1961,2005)) 
    
    #assign the line colours and set x axis to 'year' rather than 'time'
    qplt.plot(CanESM2JJA_mean.coord('year'), CanESM2JJA_mean, label='CanESM2', lw=1.5, color='blue')
    qplt.plot(CNRMGJJA_mean.coord('year'), CNRMGJJA_mean, label='CNRM_CM5', lw=1.5, color='green')
    qplt.plot(MK3JJA_mean.coord('year'), MK3JJA_mean, label='CSIRO MK3-6-0', lw=1.5, color='red')
    qplt.plot(EARTHJJA_mean.coord('year'), EARTHJJA_mean, label='EC-EARTH (r12i1p1)', lw=1.5, color='cyan')
    qplt.plot(EARTH3JJA_mean.coord('year'), EARTH3JJA_mean, label='EC-EARTH (r3i1p1', lw=1.5, color='magenta')
    qplt.plot(GFDLJJA_mean.coord('year'), GFDLJJA_mean, label='GFDL-ESM2M', lw=1.5, color='yellow')
    qplt.plot(HadGEM2JJA_mean.coord('year'), HadGEM2JJA_mean, label='HadGEM2-ES', lw=1.5, color='blue', linestyle='--')
    qplt.plot(IPSLGJJA_mean.coord('year'), IPSLGJJA_mean, label='IPSL-CM5A-MR', lw=1.5, color='green', linestyle='--')
    qplt.plot(MIROCGJJA_mean.coord('year'), MIROCGJJA_mean, label='MIROC5 ', lw=1.5, color='red', linestyle='--')
    qplt.plot(MPIJJA_mean.coord('year'), MPIJJA_mean, label='MPI-ESM-LR', lw=1.5, color='cyan', linestyle='--')
    qplt.plot(NorESM1JJA_mean.coord('year'), NorESM1JJA_mean, label='NorESM1-M', lw=1.5, color='magenta', linestyle='--')
    plt.plot(X, ObsJJAY, label='Observed', lw=3, color='black')
    plt.plot(X, AverageJJAGY, label='Average GCM', lw=3, color='grey')
    
    #set a title for the y axis
    plt.ylabel('Near-Surface Temperature (degrees Celsius)')
    
    #create a legend and set its location to under the graph
    plt.legend(loc="upper center", bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=2)
    
    #create a title
    plt.title('TasMax for Malawi JJA 1961-2005', fontsize=11)   
    
    #add grid lines
    plt.grid()
    
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_GCM_LineGraph_JJA', bbox_inches='tight')

    #show the graph in the console
    iplt.show()
    
    #limit x axis    
    plt.xlim((1961,2005)) 
    
    #assign the line colours and set x axis to 'year' rather than 'time'
    qplt.plot(CanESM2YR_mean.coord('year'), CanESM2YR_mean, label='CanESM2', lw=1.5, color='blue')
    qplt.plot(CNRMGYR_mean.coord('year'), CNRMGYR_mean, label='CNRM_CM5', lw=1.5, color='green')
    qplt.plot(MK3YR_mean.coord('year'), MK3YR_mean, label='CSIRO MK3-6-0', lw=1.5, color='red')
    qplt.plot(EARTHYR_mean.coord('year'), EARTHYR_mean, label='EC-EARTH (r12i1p1)', lw=1.5, color='cyan')
    qplt.plot(EARTH3YR_mean.coord('year'), EARTH3YR_mean, label='EC-EARTH (r3i1p1', lw=1.5, color='magenta')
    qplt.plot(GFDLYR_mean.coord('year'), GFDLYR_mean, label='GFDL-ESM2M', lw=1.5, color='yellow')
    qplt.plot(HadGEM2YR_mean.coord('year'), HadGEM2YR_mean, label='HadGEM2-ES', lw=1.5, color='blue', linestyle='--')
    qplt.plot(IPSLGYR_mean.coord('year'), IPSLGYR_mean, label='IPSL-CM5A-MR', lw=1.5, color='green', linestyle='--')
    qplt.plot(MIROCGYR_mean.coord('year'), MIROCGYR_mean, label='MIROC5 ', lw=1.5, color='red', linestyle='--')
    qplt.plot(MPIYR_mean.coord('year'), MPIYR_mean, label='MPI-ESM-LR', lw=1.5, color='cyan', linestyle='--')
    qplt.plot(NorESM1YR_mean.coord('year'), NorESM1YR_mean, label='NorESM1-M', lw=1.5, color='magenta', linestyle='--')
    plt.plot(X, ObsY, label='Observed', lw=3, color='black')
    plt.plot(X, AverageGY, label='Average GCM', lw=3, color='grey')
    
    #set a title for the y axis
    plt.ylabel('Near-Surface Temperature (degrees Celsius)')
    
    #create a legend and set its location to under the graph
    plt.legend(loc="upper center", bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=2)
    
    #create a title
    plt.title('TasMax for Malawi 1961-2005', fontsize=11)   
    
    #add grid lines
    plt.grid()
    
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_GCM_LineGraph_Annual', bbox_inches='tight')

    #show the graph in the console
    iplt.show()
    
    #PART D: Average RCM and GCM Line Graph
    #limit x axis    
    plt.xlim((1961,2005)) 
    
    #assign the line colours and set x axis to 'year' rather than 'time'
    plt.plot(X, ObsSONY, label='Observed', lw=1.5, color='black')
    plt.plot(X, AverageSONRY, label='Average RCM', lw=1.5, color='cyan')
    plt.plot(X, AverageSONGY, label='Average GCM', lw=1.5, color='magenta')
            
    #set a title for the y axis
    plt.ylabel('Near-Surface Temperature (degrees Celsius)')
    
    #create a legend and set its location to under the graph
    plt.legend(loc="upper center", bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=3)
    
    #create a title
    plt.title('TasMax for Malawi SON 1961-2005', fontsize=11)   
    
    #add grid lines
    plt.grid()
    
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Ave_LineGraph_SON', bbox_inches='tight')

    #show the graph in the console
    iplt.show() 
    
    #limit x axis    
    plt.xlim((1961,2005)) 
    
    #assign the line colours and set x axis to 'year' rather than 'time'
    plt.plot(X, ObsDJFY, label='Observed', lw=1.5, color='black')
    plt.plot(X, AverageDJFRY, label='Average RCM', lw=1.5, color='cyan')
    plt.plot(X, AverageDJFGY, label='Average GCM', lw=1.5, color='magenta')
            
    #set a title for the y axis
    plt.ylabel('Near-Surface Temperature (degrees Celsius)')
    
    #create a legend and set its location to under the graph
    plt.legend(loc="upper center", bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=3)
    
    #create a title
    plt.title('TasMax for Malawi DJF 1961-2005', fontsize=11)   
    
    #add grid lines
    plt.grid()
    
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Ave_LineGraph_DJF', bbox_inches='tight')

    #show the graph in the console
    iplt.show() 
    
    #limit x axis    
    plt.xlim((1961,2005)) 
    
    #assign the line colours and set x axis to 'year' rather than 'time'
    plt.plot(X, ObsMAMY, label='Observed', lw=1.5, color='black')
    plt.plot(X, AverageMAMRY, label='Average RCM', lw=1.5, color='cyan')
    plt.plot(X, AverageMAMGY, label='Average GCM', lw=1.5, color='magenta')
            
    #set a title for the y axis
    plt.ylabel('Near-Surface Temperature (degrees Celsius)')
    
    #create a legend and set its location to under the graph
    plt.legend(loc="upper center", bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=3)
    
    #create a title
    plt.title('TasMax for Malawi MAM 1961-2005', fontsize=11)   
    
    #add grid lines
    plt.grid()
    
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Ave_LineGraph_MAM', bbox_inches='tight')

    #show the graph in the console
    iplt.show()
   
    #limit x axis    
    plt.xlim((1961,2005)) 
    
    #assign the line colours and set x axis to 'year' rather than 'time'
    plt.plot(X, ObsJJAY, label='Observed', lw=1.5, color='black')
    plt.plot(X, AverageJJARY, label='AverageRCM', lw=1.5, color='cyan')
    plt.plot(X, AverageJJAGY, label='Average GCM', lw=1.5, color='magenta')
            
    #set a title for the y axis
    plt.ylabel('Near-Surface Temperature (degrees Celsius)')
    
    #create a legend and set its location to under the graph
    plt.legend(loc="upper center", bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=3)
    
    #create a title
    plt.title('TasMax for Malawi JJA 1961-2005', fontsize=11)   
    
    #add grid lines
    plt.grid()
    
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Ave_LineGraph_JJA', bbox_inches='tight')

    #show the graph in the console
    iplt.show()
    
    #limit x axis    
    plt.xlim((1961,2005)) 
    
    #assign the line colours and set x axis to 'year' rather than 'time'
    plt.plot(X, ObsY, label='Observed', lw=1.5, color='black')
    plt.plot(X, AverageRY, label='Average RCM', lw=1.5, color='cyan')
    plt.plot(X, AverageGY, label='Average GCM', lw=1.5, color='magenta')
            
    #set a title for the y axis
    plt.ylabel('Near-Surface Temperature (degrees Celsius)')
    
    #create a legend and set its location to under the graph
    plt.legend(loc="upper center", bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=3)
    
    #create a title
    plt.title('TasMax for Malawi 1961-2005', fontsize=11)   
    
    #add grid lines
    plt.grid()
    
    #save the image of the graph and include full legend
    plt.savefig('TasMAX_Ave_LineGraph_Annual', bbox_inches='tight')

    #show the graph in the console
    iplt.show()
    
    
if __name__ == '__main__':
    main()
    
           