"""
Created on Thursday October 25 2018

@author: s0899345
"""

import iris
import iris.analysis.cartography
import iris.coord_categorisation as iriscc
import numpy as np
import cf_units
from cf_units import Unit
import datetime


#this file is split into parts as follows:
    #PART 1: load and format all models 
    #PART 2: load and format observed data
    #PART 3: format files to be plot specific
    #PART 4: print data
    
    
def main():
    #-------------------------------------------------------------------------
    #PART 1: LOAD and FORMAT ALL MODELS   
    #bring in all the ERAINT models we need and give them a name
    CCCma_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/ERAINT/pr_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_CCCma-CanRCM4_r2_mon_198901-200912.nc'
    CLMcom_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/ERAINT/pr_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_198901-200812.nc'
    DMI_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/ERAINT/pr_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_DMI-HIRHAM5_v2_mon_198901-201012.nc'
    KNMI_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/ERAINT/pr_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_KNMI-RACMO22T_v1_mon_197901-201212.nc'
    MPIE_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/ERAINT/pr_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_MPI-CSC-REMO2009_v1_mon_198902-200812.nc'
    SMHI_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/ERAINT/pr_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_SMHI-RCA4_v1_mon_198001-201012.nc'
    
    CCCma_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/ERAINT/tas_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_CCCma-CanRCM4_r2_mon_198901-200912.nc'
    CLMcom_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/ERAINT/tas_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_198901-200812.nc'
    DMI_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/ERAINT/tas_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_DMI-HIRHAM5_v2_mon_198901-201012.nc'
    KNMI_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/ERAINT/tas_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_KNMI-RACMO22T_v1_mon_197901-201212.nc'
    MPIE_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/ERAINT/tas_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_MPI-CSC-REMO2009_v1_mon_198902-200812.nc'
    SMHI_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/ERAINT/tas_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_SMHI-RCA4_v1_mon_198001-201012.nc'
    
    CCCma_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/ERAINT/tasmin_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_CCCma-CanRCM4_r2_mon_198901-200912.nc'
    CLMcom_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/ERAINT/tasmin_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_198901-200812.nc'
    DMI_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/ERAINT/tasmin_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_DMI-HIRHAM5_v2_mon_198901-201012.nc'
    KNMI_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/ERAINT/tasmin_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_KNMI-RACMO22T_v1_mon_197901-201212.nc'
    MPIE_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/ERAINT/tasmin_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_MPI-CSC-REMO2009_v1_mon_198902-200812.nc'
    SMHI_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/ERAINT/tasmin_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_SMHI-RCA4_v1_mon_198001-201012.nc'
    
    CCCma_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/ERAINT/tasmax_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_CCCma-CanRCM4_r2_mon_198901-200912.nc'
    CLMcom_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/ERAINT/tasmax_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_198901-200812.nc'
    DMI_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/ERAINT/tasmax_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_DMI-HIRHAM5_v2_mon_198901-201012.nc'
    KNMI_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/ERAINT/tasmax_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_KNMI-RACMO22T_v1_mon_197901-201212.nc'
    MPIE_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/ERAINT/tasmax_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_MPI-CSC-REMO2009_v1_mon_198902-200812.nc'
    SMHI_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/ERAINT/tasmax_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_SMHI-RCA4_v1_mon_198001-201012.nc'
    
    #bring in all the CORDEX RCM models we need and give them a name
    CCCmaCanRCM_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/Historical_monthly/pr_AFR-44_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_mon_195101-200512.nc'
    CCCmaSMHI_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/Historical_monthly/pr_AFR-44_CCCma-CanESM2_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    CNRM_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/Historical_monthly/pr_AFR-44_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_195001-200512.nc'
    CNRMSMHI_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/Historical_monthly/pr_AFR-44_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    CSIRO_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/Historical_monthly/pr_AFR-44_CSIRO-QCCCE-CSIRO-Mk3-6-0_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    ICHECDMI_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/Historical_monthly/pr_AFR-44_ICHEC-EC-EARTH_historical_r3i1p1_DMI-HIRHAM5_v2_mon_195101-200512.nc'   
    ICHECCCLM_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/Historical_monthly/pr_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_CLMcom-CCLM4-8-17_v1_mon_194912-200512.nc'  
    ICHECKNMI_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/Historical_monthly/pr_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_KNMI-RACMO22T_v1_mon_195001-200512.nc'
    ICHECMPI_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/Historical_monthly/pr_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_MPI-CSC-REMO2009_v1_mon_195001-200512.nc'
    ICHECSMHI_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/Historical_monthly/pr_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    IPSL_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/Historical_monthly/pr_AFR-44_IPSL-IPSL-CM5A-MR_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    MIROC_pr =  '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/Historical_monthly/pr_AFR-44_MIROC-MIROC5_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc' 
    MOHCCCLM_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/Historical_monthly/pr_AFR-44_MOHC-HadGEM2-ES_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_194912-200512.nc' 
    MOHCKNMI_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/Historical_monthly/pr_AFR-44_MOHC-HadGEM2-ES_historical_r1i1p1_KNMI-RACMO22T_v2_mon_195001-200512.nc'
    MOHCSMHI_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/Historical_monthly/pr_AFR-44_MOHC-HadGEM2-ES_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    MPICCLM_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/Historical_monthly/pr_AFR-44_MPI-M-MPI-ESM-LR_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_194912-200512.nc' 
    MPIREMO_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/Historical_monthly/pr_AFR-44_MPI-M-MPI-ESM-LR_historical_r1i1p1_MPI-CSC-REMO2009_v1_mon_195001-200512.nc' 
    MPISMHI_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/Historical_monthly/pr_AFR-44_MPI-M-MPI-ESM-LR_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    NCCDMI_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/Historical_monthly/pr_AFR-44_NCC-NorESM1-M_historical_r1i1p1_DMI-HIRHAM5_v1_mon_195101-200512.nc'    
    NCCSMHI_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/Historical_monthly/pr_AFR-44_NCC-NorESM1-M_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'   
    NOAA_pr = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_pr/Historical_monthly/pr_AFR-44_NOAA-GFDL-GFDL-ESM2M_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'     
    
    CCCmaCanRCM_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/Historical_monthly/tas_AFR-44_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_mon_195001-200512.nc'
    CCCmaSMHI_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/Historical_monthly/tas_AFR-44_CCCma-CanESM2_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'   
    CNRM_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/Historical_monthly/tas_AFR-44_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_195101-200512.nc'
    CNRMSMHI_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/Historical_monthly/tas_AFR-44_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    CSIRO_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/Historical_monthly/tas_AFR-44_CSIRO-QCCCE-CSIRO-Mk3-6-0_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    ICHECDMI_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/Historical_monthly/tas_AFR-44_ICHEC-EC-EARTH_historical_r3i1p1_DMI-HIRHAM5_v2_mon_195101-200512.nc'   
    ICHECCCLM_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/Historical_monthly/tas_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_CLMcom-CCLM4-8-17_v1_mon_194912-200512.nc' 
    ICHECKNMI_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/Historical_monthly/tas_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_KNMI-RACMO22T_v1_mon_195001-200512.nc'
    ICHECMPI_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/Historical_monthly/tas_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_MPI-CSC-REMO2009_v1_mon_195001-200512.nc'
    ICHECSMHI_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/Historical_monthly/tas_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    IPSL_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/Historical_monthly/tas_AFR-44_IPSL-IPSL-CM5A-MR_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    MIROC_tas =  '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/Historical_monthly/tas_AFR-44_MIROC-MIROC5_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc' 
    MOHCCCLM_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/Historical_monthly/tas_AFR-44_MOHC-HadGEM2-ES_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_194912-200512.nc' 
    MOHCKNMI_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/Historical_monthly/tas_AFR-44_MOHC-HadGEM2-ES_historical_r1i1p1_KNMI-RACMO22T_v2_mon_195001-200512.nc'
    MOHCSMHI_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/Historical_monthly/tas_AFR-44_MOHC-HadGEM2-ES_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    MPICCLM_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/Historical_monthly/tas_AFR-44_MPI-M-MPI-ESM-LR_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_194912-200512.nc'
    MPIREMO_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/Historical_monthly/tas_AFR-44_MPI-M-MPI-ESM-LR_historical_r1i1p1_MPI-CSC-REMO2009_v1_mon_195001-200512.nc' 
    MPISMHI_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/Historical_monthly/tas_AFR-44_MPI-M-MPI-ESM-LR_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    NCCDMI_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/Historical_monthly/tas_AFR-44_NCC-NorESM1-M_historical_r1i1p1_DMI-HIRHAM5_v1_mon_195101-200512.nc'    
    NCCSMHI_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/Historical_monthly/tas_AFR-44_NCC-NorESM1-M_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'   
    NOAA_tas = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tas/Historical_monthly/tas_AFR-44_NOAA-GFDL-GFDL-ESM2M_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'  
    
    CCCmaCanRCM_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/Historical_monthly/tasmin_AFR-44_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_mon_195001-200512.nc'
    CCCmaSMHI_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/Historical_monthly/tasmin_AFR-44_CCCma-CanESM2_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'   
    CNRM_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/Historical_monthly/tasmin_AFR-44_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_195001-200512.nc'
    CNRMSMHI_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/Historical_monthly/tasmin_AFR-44_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    CSIRO_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/Historical_monthly/tasmin_AFR-44_CSIRO-QCCCE-CSIRO-Mk3-6-0_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    ICHECDMI_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/Historical_monthly/tasmin_AFR-44_ICHEC-EC-EARTH_historical_r3i1p1_DMI-HIRHAM5_v2_mon_195101-200512.nc'   
    ICHECCCLM_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/Historical_monthly/tasmin_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_CLMcom-CCLM4-8-17_v1_mon_194912-200512.nc' 
    ICHECKNMI_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/Historical_monthly/tasmin_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_KNMI-RACMO22T_v1_mon_195001-200512.nc'
    ICHECMPI_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/Historical_monthly/tasmin_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_MPI-CSC-REMO2009_v1_mon_195001-200512.nc'
    ICHECSMHI_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/Historical_monthly/tasmin_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    IPSL_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/Historical_monthly/tasmin_AFR-44_IPSL-IPSL-CM5A-MR_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    MIROC_tasmin =  '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/Historical_monthly/tasmin_AFR-44_MIROC-MIROC5_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc' 
    MOHCCCLM_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/Historical_monthly/tasmin_AFR-44_MOHC-HadGEM2-ES_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_194912-200512.nc' 
    MOHCKNMI_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/Historical_monthly/tasmin_AFR-44_MOHC-HadGEM2-ES_historical_r1i1p1_KNMI-RACMO22T_v2_mon_195001-200512.nc'
    MOHCSMHI_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/Historical_monthly/tasmin_AFR-44_MOHC-HadGEM2-ES_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    MPICCLM_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/Historical_monthly/tasmin_AFR-44_MPI-M-MPI-ESM-LR_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_194912-200512.nc'
    MPIREMO_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/Historical_monthly/tasmin_AFR-44_MPI-M-MPI-ESM-LR_historical_r1i1p1_MPI-CSC-REMO2009_v1_mon_195001-200512.nc' 
    MPISMHI_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/Historical_monthly/tasmin_AFR-44_MPI-M-MPI-ESM-LR_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    NCCDMI_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/Historical_monthly/tasmin_AFR-44_NCC-NorESM1-M_historical_r1i1p1_DMI-HIRHAM5_v1_mon_195101-200512.nc'    
    NCCSMHI_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/Historical_monthly/tasmin_AFR-44_NCC-NorESM1-M_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'   
    NOAA_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmin/Historical_monthly/tasmin_AFR-44_NOAA-GFDL-GFDL-ESM2M_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    
    CCCmaCanRCM_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/Historical_monthly/tasmax_AFR-44_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_mon_195001-200512.nc'
    CCCmaSMHI_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/Historical_monthly/tasmax_AFR-44_CCCma-CanESM2_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'   
    CNRM_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/Historical_monthly/tasmax_AFR-44_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_195001-200512.nc'
    CNRMSMHI_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/Historical_monthly/tasmax_AFR-44_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    CSIRO_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/Historical_monthly/tasmax_AFR-44_CSIRO-QCCCE-CSIRO-Mk3-6-0_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    ICHECDMI_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/Historical_monthly/tasmax_AFR-44_ICHEC-EC-EARTH_historical_r3i1p1_DMI-HIRHAM5_v2_mon_195101-200512.nc'   
    ICHECCCLM_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/Historical_monthly/tasmax_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_CLMcom-CCLM4-8-17_v1_mon_194912-200512.nc' 
    ICHECKNMI_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/Historical_monthly/tasmax_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_KNMI-RACMO22T_v1_mon_195001-200512.nc'
    ICHECMPI_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/Historical_monthly/tasmax_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_MPI-CSC-REMO2009_v1_mon_195001-200512.nc'
    ICHECSMHI_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/Historical_monthly/tasmax_AFR-44_ICHEC-EC-EARTH_historical_r12i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    IPSL_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/Historical_monthly/tasmax_AFR-44_IPSL-IPSL-CM5A-MR_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    MIROC_tasmax =  '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/Historical_monthly/tasmax_AFR-44_MIROC-MIROC5_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc' 
    MOHCCCLM_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/Historical_monthly/tasmax_AFR-44_MOHC-HadGEM2-ES_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_194912-200512.nc' 
    MOHCKNMI_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/Historical_monthly/tasmax_AFR-44_MOHC-HadGEM2-ES_historical_r1i1p1_KNMI-RACMO22T_v2_mon_195001-200512.nc'
    MOHCSMHI_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/Historical_monthly/tasmax_AFR-44_MOHC-HadGEM2-ES_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    MPICCLM_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/Historical_monthly/tasmax_AFR-44_MPI-M-MPI-ESM-LR_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_194912-200512.nc'
    MPIREMO_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/Historical_monthly/tasmax_AFR-44_MPI-M-MPI-ESM-LR_historical_r1i1p1_MPI-CSC-REMO2009_v1_mon_195001-200512.nc' 
    MPISMHI_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/Historical_monthly/tasmax_AFR-44_MPI-M-MPI-ESM-LR_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    NCCDMI_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/Historical_monthly/tasmax_AFR-44_NCC-NorESM1-M_historical_r1i1p1_DMI-HIRHAM5_v1_mon_195101-200512.nc'    
    NCCSMHI_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/Historical_monthly/tasmax_AFR-44_NCC-NorESM1-M_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'   
    NOAA_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/AFR_44_tasmax/Historical_monthly/tasmax_AFR-44_NOAA-GFDL-GFDL-ESM2M_historical_r1i1p1_SMHI-RCA4_v1_mon_195101-200512.nc'
    
    #bring in all the GCM models we need and give them a name
    CanESM2_pr = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/pr_Amon_CanESM2_historical_r1i1p1_185001-200512.nc'
    CNRMG_pr = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/pr_Amon_CNRM-CM5-2_historical_r1i1p1_195001-200512.nc'
    MK3_pr = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/pr_Amon_CSIRO-Mk3-6-0_historical_r1i1p1_185001-200512.nc'
    EARTH_pr = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/pr_Amon_EC-EARTH_historical_r12i1p1_195001-201212.nc'
    EARTH3_pr = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/pr_Amon_EC-EARTH_historical_r3i1p1_1850-2009.nc'
    GFDL_pr ='/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/pr_Amon_GFDL-ESM2M_historical_r1i1p1_194601-200512.nc'
    HadGEM2_pr ='/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/pr_Amon_HadGEM2-ES_historical_r1i1p1_193412-200511.nc'
    IPSLG_pr ='/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/pr_Amon_IPSL-CM5B-LR_historical_r1i1p1_185001-200512.nc'
    MIROCG_pr ='/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/pr_Amon_MIROC5_historical_r1i1p1_185001-201212.nc'
    MPI_pr ='/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/pr_Amon_MPI-ESM-MR_historical_r1i1p1_185001-200512.nc'
    NorESM1_pr ='/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/pr_Amon_NorESM1-M_historical_r1i1p1_185001-200512.nc'
    
    CanESM2_tas = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tas_Amon_CanESM2_historical_r1i1p1_185001-200512.nc'
    CNRMG_tas = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tas_Amon_CNRM-CM5-2_historical_r1i1p1_195001-200512.nc'
    MK3_tas = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tas_Amon_CSIRO-Mk3-6-0_historical_r1i1p1_185001-200512.nc'    
    EARTH_tas = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tas_Amon_EC-EARTH_historical_r12i1p1_195001-201212.nc'
    EARTH3_tas = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tas_Amon_EC-EARTH_historical_r3i1p1_1850-2009.nc'
    GFDL_tas = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tas_Amon_GFDL-ESM2M_historical_r1i1p1_194601-200512.nc'
    HadGEM2_tas = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tas_Amon_HadGEM2-ES_historical_r1i1p1_193412-200511.nc'
    IPSLG_tas = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tas_Amon_IPSL-CM5B-LR_historical_r1i1p1_185001-200512.nc'
    MIROCG_tas = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tas_Amon_MIROC5_historical_r1i1p1_185001-201212.nc'
    MPI_tas = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tas_Amon_MPI-ESM-MR_historical_r1i1p1_185001-200512.nc'
    NorESM1_tas = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tas_Amon_NorESM1-M_historical_r1i1p1_185001-200512.nc'
    
    CanESM2_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tasmin_Amon_CanESM2_historical_r1i1p1_185001-200512.nc'
    CNRMG_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tasmin_Amon_CNRM-CM5-2_historical_r1i1p1_195001-200512.nc'
    MK3_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tasmin_Amon_CSIRO-Mk3-6-0_historical_r1i1p1_185001-200512.nc'    
    EARTH_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tasmin_Amon_EC-EARTH_historical_r12i1p1_195001-201212.nc'
    EARTH3_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tasmax_Amon_EC-EARTH_historical+r3i1p1_1961-2005.nc'
    GFDL_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tasmin_Amon_GFDL-ESM2M_historical_r1i1p1_194601-200512.nc'
    HadGEM2_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tasmin_Amon_HadGEM2-ES_historical_r1i1p1_193412-200511.nc'
    IPSLG_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tasmin_Amon_IPSL-CM5B-LR_historical_r1i1p1_185001-200512.nc'
    MIROCG_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tasmin_Amon_MIROC5_historical_r1i1p1_185001-201212.nc'
    MPI_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tasmin_Amon_MPI-ESM-MR_historical_r1i1p1_185001-200512.nc'
    NorESM1_tasmin = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tasmin_Amon_NorESM1-M_historical_r1i1p1_185001-200512.nc'
    
    CanESM2_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tasmax_Amon_CanESM2_historical_r1i1p1_185001-200512.nc'
    CNRMG_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tasmax_Amon_CNRM-CM5-2_historical_r1i1p1_195001-200512.nc'
    MK3_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tasmax_Amon_CSIRO-Mk3-6-0_historical_r1i1p1_185001-200512.nc'    
    EARTH_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tasmax_Amon_EC-EARTH_historical_r12i1p1_195001-201212.nc'
    EARTH3_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tasmax_Amon_EC-EARTH_historical+r3i1p1_1961-2005.nc'
    GFDL_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tasmax_Amon_GFDL-ESM2M_historical_r1i1p1_194601-200512.nc'
    HadGEM2_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tasmax_Amon_HadGEM2-ES_historical_r1i1p1_193412-200511.nc'
    IPSLG_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tasmax_Amon_IPSL-CM5B-LR_historical_r1i1p1_185001-200512.nc'
    MIROCG_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tasmax_Amon_MIROC5_historical_r1i1p1_185001-201212.nc'
    MPI_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tasmax_Amon_MPI-ESM-MR_historical_r1i1p1_185001-200512.nc'
    NorESM1_tasmax = '/exports/csce/datastore/geos/groups/s0899345_2/GCM_data/tasmax_Amon_NorESM1-M_historical_r1i1p1_185001-200512.nc'
    
    #Load exactly one cube from given file
    CCCma_pr = iris.load_cube(CCCma_pr)
    CLMcom_pr = iris.load_cube(CLMcom_pr)
    DMI_pr = iris.load_cube(DMI_pr, 'precipitation_flux')
    KNMI_pr = iris.load_cube(KNMI_pr)
    MPIE_pr = iris.load_cube(MPIE_pr)
    SMHI_pr = iris.load_cube(SMHI_pr)
    
    CCCmaCanRCM_pr = iris.load_cube(CCCmaCanRCM_pr)
    CCCmaSMHI_pr = iris.load_cube(CCCmaSMHI_pr)
    CNRM_pr = iris.load_cube(CNRM_pr)
    CNRMSMHI_pr = iris.load_cube(CNRMSMHI_pr)
    CSIRO_pr = iris.load_cube(CSIRO_pr)
    ICHECDMI_pr = iris.load_cube(ICHECDMI_pr, 'precipitation_flux')
    ICHECCCLM_pr = iris.load_cube(ICHECCCLM_pr)
    ICHECKNMI_pr = iris.load_cube(ICHECKNMI_pr)
    ICHECMPI_pr = iris.load_cube(ICHECMPI_pr)
    ICHECSMHI_pr = iris.load_cube(ICHECSMHI_pr)
    IPSL_pr = iris.load_cube(IPSL_pr)
    MIROC_pr = iris.load_cube(MIROC_pr)
    MOHCCCLM_pr = iris.load_cube(MOHCCCLM_pr)
    MOHCKNMI_pr = iris.load_cube(MOHCKNMI_pr)
    MOHCSMHI_pr = iris.load_cube(MOHCSMHI_pr)
    MPICCLM_pr = iris.load_cube(MPICCLM_pr)
    MPIREMO_pr = iris.load_cube(MPIREMO_pr)
    MPISMHI_pr = iris.load_cube(MPISMHI_pr)
    NCCDMI_pr = iris.load_cube(NCCDMI_pr, 'precipitation_flux')
    NCCSMHI_pr = iris.load_cube(NCCSMHI_pr)
    NOAA_pr = iris.load_cube(NOAA_pr)
    
    CanESM2_pr = iris.load_cube(CanESM2_pr)
    CNRMG_pr = iris.load_cube(CNRMG_pr)
    MK3_pr = iris.load_cube(MK3_pr)
    EARTH_pr = iris.load_cube(EARTH_pr)
    EARTH3_pr = iris.load_cube(EARTH3_pr)
    GFDL_pr = iris.load_cube(GFDL_pr, 'precipitation_flux')
    HadGEM2_pr = iris.load_cube(HadGEM2_pr)
    IPSLG_pr = iris.load_cube(IPSLG_pr)
    MIROCG_pr = iris.load_cube(MIROCG_pr)
    MPI_pr = iris.load_cube(MPI_pr)
    NorESM1_pr = iris.load_cube(NorESM1_pr)
    
    CCCma_tas = iris.load_cube(CCCma_tas)
    CLMcom_tas = iris.load_cube(CLMcom_tas)
    DMI_tas = iris.load_cube(DMI_tas, 'air_temperature')
    KNMI_tas = iris.load_cube(KNMI_tas)
    MPIE_tas = iris.load_cube(MPIE_tas)
    SMHI_tas = iris.load_cube(SMHI_tas)
    
    CCCmaCanRCM_tas = iris.load_cube(CCCmaCanRCM_tas)
    CCCmaSMHI_tas = iris.load_cube(CCCmaSMHI_tas)
    CNRM_tas = iris.load_cube(CNRM_tas)
    CNRMSMHI_tas = iris.load_cube(CNRMSMHI_tas)
    CSIRO_tas = iris.load_cube(CSIRO_tas)
    ICHECDMI_tas = iris.load_cube(ICHECDMI_tas, 'air_temperature')
    ICHECCCLM_tas = iris.load_cube(ICHECCCLM_tas)
    ICHECKNMI_tas = iris.load_cube(ICHECKNMI_tas)
    ICHECMPI_tas = iris.load_cube(ICHECMPI_tas)
    ICHECSMHI_tas = iris.load_cube(ICHECSMHI_tas)
    IPSL_tas = iris.load_cube(IPSL_tas)
    MIROC_tas = iris.load_cube(MIROC_tas)
    MOHCCCLM_tas = iris.load_cube(MOHCCCLM_tas)
    MOHCKNMI_tas = iris.load_cube(MOHCKNMI_tas)
    MOHCSMHI_tas = iris.load_cube(MOHCSMHI_tas)
    MPICCLM_tas = iris.load_cube(MPICCLM_tas)
    MPIREMO_tas = iris.load_cube(MPIREMO_tas)
    MPISMHI_tas = iris.load_cube(MPISMHI_tas)
    NCCDMI_tas = iris.load_cube(NCCDMI_tas, 'air_temperature')
    NCCSMHI_tas = iris.load_cube(NCCSMHI_tas)
    NOAA_tas = iris.load_cube(NOAA_tas)
    
    CanESM2_tas = iris.load_cube(CanESM2_tas)
    CNRMG_tas = iris.load_cube(CNRMG_tas)
    MK3_tas = iris.load_cube(MK3_tas)
    EARTH_tas = iris.load_cube(EARTH_tas)
    EARTH3_tas = iris.load_cube(EARTH3_tas)
    GFDL_tas = iris.load_cube(GFDL_tas, 'air_temperature')
    HadGEM2_tas = iris.load_cube(HadGEM2_tas)
    IPSLG_tas = iris.load_cube(IPSLG_tas)
    MIROCG_tas = iris.load_cube(MIROCG_tas)
    MPI_tas = iris.load_cube(MPI_tas)
    NorESM1_tas = iris.load_cube(NorESM1_tas)
    
    CCCma_tasmin = iris.load_cube(CCCma_tasmin)
    CLMcom_tasmin = iris.load_cube(CLMcom_tasmin)
    DMI_tasmin = iris.load_cube(DMI_tasmin, 'air_temperature')
    KNMI_tasmin = iris.load_cube(KNMI_tasmin)
    MPIE_tasmin = iris.load_cube(MPIE_tasmin)
    SMHI_tasmin = iris.load_cube(SMHI_tasmin)
    
    CCCmaCanRCM_tasmin = iris.load_cube(CCCmaCanRCM_tasmin)
    CCCmaSMHI_tasmin = iris.load_cube(CCCmaSMHI_tasmin)
    CNRM_tasmin = iris.load_cube(CNRM_tasmin)
    CNRMSMHI_tasmin = iris.load_cube(CNRMSMHI_tasmin)
    CSIRO_tasmin = iris.load_cube(CSIRO_tasmin)
    ICHECDMI_tasmin = iris.load_cube(ICHECDMI_tasmin, 'air_temperature')
    ICHECCCLM_tasmin = iris.load_cube(ICHECCCLM_tasmin)
    ICHECKNMI_tasmin = iris.load_cube(ICHECKNMI_tasmin)
    ICHECMPI_tasmin = iris.load_cube(ICHECMPI_tasmin)
    ICHECSMHI_tasmin = iris.load_cube(ICHECSMHI_tasmin)
    IPSL_tasmin = iris.load_cube(IPSL_tasmin)
    MIROC_tasmin = iris.load_cube(MIROC_tasmin)
    MOHCCCLM_tasmin = iris.load_cube(MOHCCCLM_tasmin)
    MOHCKNMI_tasmin = iris.load_cube(MOHCKNMI_tasmin)
    MOHCSMHI_tasmin = iris.load_cube(MOHCSMHI_tasmin)
    MPICCLM_tasmin = iris.load_cube(MPICCLM_tasmin)
    MPIREMO_tasmin = iris.load_cube(MPIREMO_tasmin)
    MPISMHI_tasmin = iris.load_cube(MPISMHI_tasmin)
    NCCDMI_tasmin = iris.load_cube(NCCDMI_tasmin, 'air_temperature')
    NCCSMHI_tasmin = iris.load_cube(NCCSMHI_tasmin)
    NOAA_tasmin = iris.load_cube(NOAA_tasmin)
    
    CanESM2_tasmin = iris.load_cube(CanESM2_tasmin)
    CNRMG_tasmin = iris.load_cube(CNRMG_tasmin)
    MK3_tasmin = iris.load_cube(MK3_tasmin)
    EARTH_tasmin = iris.load_cube(EARTH_tasmin)
    EARTH3_tasmin = iris.load_cube(EARTH3_tasmin)
    GFDL_tasmin = iris.load_cube(GFDL_tasmin, 'air_temperature')
    HadGEM2_tasmin = iris.load_cube(HadGEM2_tasmin)
    IPSLG_tasmin = iris.load_cube(IPSLG_tasmin)
    MIROCG_tasmin = iris.load_cube(MIROCG_tasmin)
    MPI_tasmin = iris.load_cube(MPI_tasmin)
    NorESM1_tasmin = iris.load_cube(NorESM1_tasmin)
    
    CCCma_tasmax = iris.load_cube(CCCma_tasmax)
    CLMcom_tasmax = iris.load_cube(CLMcom_tasmax)
    DMI_tasmax = iris.load_cube(DMI_tasmax, 'air_temperature')
    KNMI_tasmax = iris.load_cube(KNMI_tasmax)
    MPIE_tasmax = iris.load_cube(MPIE_tasmax)
    SMHI_tasmax = iris.load_cube(SMHI_tasmax)
    
    CCCmaCanRCM_tasmax = iris.load_cube(CCCmaCanRCM_tasmax)
    CCCmaSMHI_tasmax = iris.load_cube(CCCmaSMHI_tasmax)
    CNRM_tasmax = iris.load_cube(CNRM_tasmax)
    CNRMSMHI_tasmax = iris.load_cube(CNRMSMHI_tasmax)
    CSIRO_tasmax = iris.load_cube(CSIRO_tasmax)
    ICHECDMI_tasmax = iris.load_cube(ICHECDMI_tasmax, 'air_temperature')
    ICHECCCLM_tasmax = iris.load_cube(ICHECCCLM_tasmax)
    ICHECKNMI_tasmax = iris.load_cube(ICHECKNMI_tasmax)
    ICHECMPI_tasmax = iris.load_cube(ICHECMPI_tasmax)
    ICHECSMHI_tasmax = iris.load_cube(ICHECSMHI_tasmax)
    IPSL_tasmax = iris.load_cube(IPSL_tasmax)
    MIROC_tasmax = iris.load_cube(MIROC_tasmax)
    MOHCCCLM_tasmax = iris.load_cube(MOHCCCLM_tasmax)
    MOHCKNMI_tasmax = iris.load_cube(MOHCKNMI_tasmax)
    MOHCSMHI_tasmax = iris.load_cube(MOHCSMHI_tasmax)
    MPICCLM_tasmax = iris.load_cube(MPICCLM_tasmax)
    MPIREMO_tasmax = iris.load_cube(MPIREMO_tasmax)
    MPISMHI_tasmax = iris.load_cube(MPISMHI_tasmax)
    NCCDMI_tasmax = iris.load_cube(NCCDMI_tasmax, 'air_temperature')
    NCCSMHI_tasmax = iris.load_cube(NCCSMHI_tasmax)
    NOAA_tasmax = iris.load_cube(NOAA_tasmax)
    
    CanESM2_tasmax = iris.load_cube(CanESM2_tasmax)
    CNRMG_tasmax = iris.load_cube(CNRMG_tasmax)
    MK3_tasmax = iris.load_cube(MK3_tasmax)
    EARTH_tasmax = iris.load_cube(EARTH_tasmax)
    EARTH3_tasmax = iris.load_cube(EARTH3_tasmax)
    GFDL_tasmax = iris.load_cube(GFDL_tasmax, 'air_temperature')
    HadGEM2_tasmax = iris.load_cube(HadGEM2_tasmax)
    IPSLG_tasmax = iris.load_cube(IPSLG_tasmax)
    MIROCG_tasmax = iris.load_cube(MIROCG_tasmax)
    MPI_tasmax = iris.load_cube(MPI_tasmax)
    NorESM1_tasmax = iris.load_cube(NorESM1_tasmax)
    
    #remove flat latitude and longitude and only use grid latitude and grid longitude to make consistent with the observed data, also make sure all of the longitudes are monotonic. This is only applicable to the ERAINT and CORDEX Models as the Observed data and GCMs are already in standard linear format.    
    lats = iris.coords.DimCoord(CCCma_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CCCma_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
                                
    CCCma_pr.remove_coord('latitude')
    CCCma_pr.remove_coord('longitude')
    CCCma_pr.remove_coord('grid_latitude')
    CCCma_pr.remove_coord('grid_longitude')
    CCCma_pr.add_dim_coord(lats, 1)
    CCCma_pr.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(CLMcom_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CLMcom_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    CLMcom_pr.remove_coord('latitude')
    CLMcom_pr.remove_coord('longitude')
    CLMcom_pr.remove_coord('grid_latitude')
    CLMcom_pr.remove_coord('grid_longitude')
    CLMcom_pr.add_dim_coord(lats, 1)
    CLMcom_pr.add_dim_coord(lons, 2)  
    
    lats = iris.coords.DimCoord(DMI_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = DMI_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    DMI_pr.remove_coord('latitude')
    DMI_pr.remove_coord('longitude')
    DMI_pr.remove_coord('grid_latitude')
    DMI_pr.remove_coord('grid_longitude')
    DMI_pr.add_dim_coord(lats, 1)
    DMI_pr.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(KNMI_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = KNMI_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
        
    KNMI_pr.remove_coord('latitude')
    KNMI_pr.remove_coord('longitude')
    KNMI_pr.remove_coord('grid_latitude')
    KNMI_pr.remove_coord('grid_longitude')
    KNMI_pr.add_dim_coord(lats, 1)
    KNMI_pr.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(MPIE_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MPIE_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
   
    MPIE_pr.remove_coord('latitude')
    MPIE_pr.remove_coord('longitude')
    MPIE_pr.remove_coord('grid_latitude')
    MPIE_pr.remove_coord('grid_longitude')
    MPIE_pr.add_dim_coord(lats, 1)
    MPIE_pr.add_dim_coord(lons, 2)

    lats = iris.coords.DimCoord(SMHI_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = SMHI_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    SMHI_pr.remove_coord('latitude')
    SMHI_pr.remove_coord('longitude')
    SMHI_pr.remove_coord('grid_latitude')
    SMHI_pr.remove_coord('grid_longitude')
    SMHI_pr.add_dim_coord(lats, 1)
    SMHI_pr.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(CCCmaCanRCM_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CCCmaCanRCM_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
                              
    CCCmaCanRCM_pr.remove_coord('latitude')
    CCCmaCanRCM_pr.remove_coord('longitude')
    CCCmaCanRCM_pr.remove_coord('grid_latitude')
    CCCmaCanRCM_pr.remove_coord('grid_longitude')
    CCCmaCanRCM_pr.add_dim_coord(lats, 1)
    CCCmaCanRCM_pr.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(CCCmaSMHI_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CCCmaSMHI_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    CCCmaSMHI_pr.remove_coord('latitude')
    CCCmaSMHI_pr.remove_coord('longitude')
    CCCmaSMHI_pr.remove_coord('grid_latitude')
    CCCmaSMHI_pr.remove_coord('grid_longitude')
    CCCmaSMHI_pr.add_dim_coord(lats, 1)
    CCCmaSMHI_pr.add_dim_coord(lons, 2)  
    
    lats = iris.coords.DimCoord(CNRM_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CNRM_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
                                
    CNRM_pr.remove_coord('latitude')
    CNRM_pr.remove_coord('longitude')
    CNRM_pr.remove_coord('grid_latitude')
    CNRM_pr.remove_coord('grid_longitude')
    CNRM_pr.add_dim_coord(lats, 1)
    CNRM_pr.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(CNRMSMHI_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CNRMSMHI_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    CNRMSMHI_pr.remove_coord('latitude')
    CNRMSMHI_pr.remove_coord('longitude')
    CNRMSMHI_pr.remove_coord('grid_latitude')
    CNRMSMHI_pr.remove_coord('grid_longitude')
    CNRMSMHI_pr.add_dim_coord(lats, 1)
    CNRMSMHI_pr.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(CSIRO_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CSIRO_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    CSIRO_pr.remove_coord('latitude')
    CSIRO_pr.remove_coord('longitude')
    CSIRO_pr.remove_coord('grid_latitude')
    CSIRO_pr.remove_coord('grid_longitude')
    CSIRO_pr.add_dim_coord(lats, 1)
    CSIRO_pr.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(ICHECDMI_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECDMI_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')    
    
    ICHECDMI_pr.remove_coord('latitude')
    ICHECDMI_pr.remove_coord('longitude')
    ICHECDMI_pr.remove_coord('grid_latitude')
    ICHECDMI_pr.remove_coord('grid_longitude')
    ICHECDMI_pr.add_dim_coord(lats, 1)
    ICHECDMI_pr.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(ICHECCCLM_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECCCLM_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees') 
    
    ICHECCCLM_pr.remove_coord('latitude')
    ICHECCCLM_pr.remove_coord('longitude')
    ICHECCCLM_pr.remove_coord('grid_latitude')
    ICHECCCLM_pr.remove_coord('grid_longitude')
    ICHECCCLM_pr.add_dim_coord(lats, 1)
    ICHECCCLM_pr.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(ICHECKNMI_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECKNMI_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees') 
    
    ICHECKNMI_pr.remove_coord('latitude')
    ICHECKNMI_pr.remove_coord('longitude')
    ICHECKNMI_pr.remove_coord('grid_latitude')
    ICHECKNMI_pr.remove_coord('grid_longitude')
    ICHECKNMI_pr.add_dim_coord(lats, 1)
    ICHECKNMI_pr.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(ICHECMPI_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECMPI_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees') 
    
    ICHECMPI_pr.remove_coord('latitude')
    ICHECMPI_pr.remove_coord('longitude')
    ICHECMPI_pr.remove_coord('grid_latitude')
    ICHECMPI_pr.remove_coord('grid_longitude')
    ICHECMPI_pr.add_dim_coord(lats, 1)
    ICHECMPI_pr.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(ICHECSMHI_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECSMHI_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    ICHECSMHI_pr.remove_coord('latitude')
    ICHECSMHI_pr.remove_coord('longitude')
    ICHECSMHI_pr.remove_coord('grid_latitude')
    ICHECSMHI_pr.remove_coord('grid_longitude')
    ICHECSMHI_pr.add_dim_coord(lats, 1)
    ICHECSMHI_pr.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(IPSL_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = IPSL_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    IPSL_pr.remove_coord('latitude')
    IPSL_pr.remove_coord('longitude')
    IPSL_pr.remove_coord('grid_latitude')
    IPSL_pr.remove_coord('grid_longitude')
    IPSL_pr.add_dim_coord(lats, 1)
    IPSL_pr.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MIROC_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MIROC_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
                                
    MIROC_pr.remove_coord('latitude')
    MIROC_pr.remove_coord('longitude')
    MIROC_pr.remove_coord('grid_latitude')
    MIROC_pr.remove_coord('grid_longitude')
    MIROC_pr.add_dim_coord(lats, 1)
    MIROC_pr.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MOHCCCLM_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MOHCCCLM_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MOHCCCLM_pr.remove_coord('latitude')
    MOHCCCLM_pr.remove_coord('longitude')
    MOHCCCLM_pr.remove_coord('grid_latitude')
    MOHCCCLM_pr.remove_coord('grid_longitude')
    MOHCCCLM_pr.add_dim_coord(lats, 1)
    MOHCCCLM_pr.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MOHCKNMI_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MOHCKNMI_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MOHCKNMI_pr.remove_coord('latitude')
    MOHCKNMI_pr.remove_coord('longitude')
    MOHCKNMI_pr.remove_coord('grid_latitude')
    MOHCKNMI_pr.remove_coord('grid_longitude')
    MOHCKNMI_pr.add_dim_coord(lats, 1)
    MOHCKNMI_pr.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MOHCSMHI_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MOHCSMHI_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MOHCSMHI_pr.remove_coord('latitude')
    MOHCSMHI_pr.remove_coord('longitude')
    MOHCSMHI_pr.remove_coord('grid_latitude')
    MOHCSMHI_pr.remove_coord('grid_longitude')
    MOHCSMHI_pr.add_dim_coord(lats, 1)
    MOHCSMHI_pr.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MPICCLM_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MPICCLM_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MPICCLM_pr.remove_coord('latitude')
    MPICCLM_pr.remove_coord('longitude')
    MPICCLM_pr.remove_coord('grid_latitude')
    MPICCLM_pr.remove_coord('grid_longitude')
    MPICCLM_pr.add_dim_coord(lats, 1)
    MPICCLM_pr.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MPIREMO_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MPIREMO_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MPIREMO_pr.remove_coord('latitude')
    MPIREMO_pr.remove_coord('longitude')
    MPIREMO_pr.remove_coord('grid_latitude')
    MPIREMO_pr.remove_coord('grid_longitude')
    MPIREMO_pr.add_dim_coord(lats, 1)
    MPIREMO_pr.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MPISMHI_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MPISMHI_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MPISMHI_pr.remove_coord('latitude')
    MPISMHI_pr.remove_coord('longitude')
    MPISMHI_pr.remove_coord('grid_latitude')
    MPISMHI_pr.remove_coord('grid_longitude')
    MPISMHI_pr.add_dim_coord(lats, 1)
    MPISMHI_pr.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(NCCDMI_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = NCCDMI_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    NCCDMI_pr.remove_coord('latitude')
    NCCDMI_pr.remove_coord('longitude')
    NCCDMI_pr.remove_coord('grid_latitude')
    NCCDMI_pr.remove_coord('grid_longitude')
    NCCDMI_pr.add_dim_coord(lats, 1)
    NCCDMI_pr.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(NCCSMHI_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = NCCSMHI_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    NCCSMHI_pr.remove_coord('latitude')
    NCCSMHI_pr.remove_coord('longitude')
    NCCSMHI_pr.remove_coord('grid_latitude')
    NCCSMHI_pr.remove_coord('grid_longitude')
    NCCSMHI_pr.add_dim_coord(lats, 1)
    NCCSMHI_pr.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(NOAA_pr.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = NOAA_pr.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    NOAA_pr.remove_coord('latitude')
    NOAA_pr.remove_coord('longitude')
    NOAA_pr.remove_coord('grid_latitude')
    NOAA_pr.remove_coord('grid_longitude')
    NOAA_pr.add_dim_coord(lats, 1)
    NOAA_pr.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(CCCma_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CCCma_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
                                
    CCCma_tas.remove_coord('latitude')
    CCCma_tas.remove_coord('longitude')
    CCCma_tas.remove_coord('grid_latitude')
    CCCma_tas.remove_coord('grid_longitude')
    CCCma_tas.add_dim_coord(lats, 1)
    CCCma_tas.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(CLMcom_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CLMcom_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    CLMcom_tas.remove_coord('latitude')
    CLMcom_tas.remove_coord('longitude')
    CLMcom_tas.remove_coord('grid_latitude')
    CLMcom_tas.remove_coord('grid_longitude')
    CLMcom_tas.add_dim_coord(lats, 1)
    CLMcom_tas.add_dim_coord(lons, 2)  
    
    lats = iris.coords.DimCoord(DMI_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = DMI_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    DMI_tas.remove_coord('latitude')
    DMI_tas.remove_coord('longitude')
    DMI_tas.remove_coord('grid_latitude')
    DMI_tas.remove_coord('grid_longitude')
    DMI_tas.add_dim_coord(lats, 1)
    DMI_tas.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(KNMI_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = KNMI_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
        
    KNMI_tas.remove_coord('latitude')
    KNMI_tas.remove_coord('longitude')
    KNMI_tas.remove_coord('grid_latitude')
    KNMI_tas.remove_coord('grid_longitude')
    KNMI_tas.add_dim_coord(lats, 1)
    KNMI_tas.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(MPIE_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MPIE_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
   
    MPIE_tas.remove_coord('latitude')
    MPIE_tas.remove_coord('longitude')
    MPIE_tas.remove_coord('grid_latitude')
    MPIE_tas.remove_coord('grid_longitude')
    MPIE_tas.add_dim_coord(lats, 1)
    MPIE_tas.add_dim_coord(lons, 2)

    lats = iris.coords.DimCoord(SMHI_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = SMHI_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    SMHI_tas.remove_coord('latitude')
    SMHI_tas.remove_coord('longitude')
    SMHI_tas.remove_coord('grid_latitude')
    SMHI_tas.remove_coord('grid_longitude')
    SMHI_tas.add_dim_coord(lats, 1)
    SMHI_tas.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(CCCmaCanRCM_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CCCmaCanRCM_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
                              
    CCCmaCanRCM_tas.remove_coord('latitude')
    CCCmaCanRCM_tas.remove_coord('longitude')
    CCCmaCanRCM_tas.remove_coord('grid_latitude')
    CCCmaCanRCM_tas.remove_coord('grid_longitude')
    CCCmaCanRCM_tas.add_dim_coord(lats, 1)
    CCCmaCanRCM_tas.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(CCCmaSMHI_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CCCmaSMHI_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    CCCmaSMHI_tas.remove_coord('latitude')
    CCCmaSMHI_tas.remove_coord('longitude')
    CCCmaSMHI_tas.remove_coord('grid_latitude')
    CCCmaSMHI_tas.remove_coord('grid_longitude')
    CCCmaSMHI_tas.add_dim_coord(lats, 1)
    CCCmaSMHI_tas.add_dim_coord(lons, 2)  
    
    lats = iris.coords.DimCoord(CNRM_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CNRM_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
                                
    CNRM_tas.remove_coord('latitude')
    CNRM_tas.remove_coord('longitude')
    CNRM_tas.remove_coord('grid_latitude')
    CNRM_tas.remove_coord('grid_longitude')
    CNRM_tas.add_dim_coord(lats, 1)
    CNRM_tas.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(CNRMSMHI_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CNRMSMHI_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    CNRMSMHI_tas.remove_coord('latitude')
    CNRMSMHI_tas.remove_coord('longitude')
    CNRMSMHI_tas.remove_coord('grid_latitude')
    CNRMSMHI_tas.remove_coord('grid_longitude')
    CNRMSMHI_tas.add_dim_coord(lats, 1)
    CNRMSMHI_tas.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(CSIRO_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CSIRO_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    CSIRO_tas.remove_coord('latitude')
    CSIRO_tas.remove_coord('longitude')
    CSIRO_tas.remove_coord('grid_latitude')
    CSIRO_tas.remove_coord('grid_longitude')
    CSIRO_tas.add_dim_coord(lats, 1)
    CSIRO_tas.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(ICHECDMI_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECDMI_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')    
    
    ICHECDMI_tas.remove_coord('latitude')
    ICHECDMI_tas.remove_coord('longitude')
    ICHECDMI_tas.remove_coord('grid_latitude')
    ICHECDMI_tas.remove_coord('grid_longitude')
    ICHECDMI_tas.add_dim_coord(lats, 1)
    ICHECDMI_tas.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(ICHECCCLM_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECCCLM_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees') 
    
    ICHECCCLM_tas.remove_coord('latitude')
    ICHECCCLM_tas.remove_coord('longitude')
    ICHECCCLM_tas.remove_coord('grid_latitude')
    ICHECCCLM_tas.remove_coord('grid_longitude')
    ICHECCCLM_tas.add_dim_coord(lats, 1)
    ICHECCCLM_tas.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(ICHECKNMI_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECKNMI_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees') 
    
    ICHECKNMI_tas.remove_coord('latitude')
    ICHECKNMI_tas.remove_coord('longitude')
    ICHECKNMI_tas.remove_coord('grid_latitude')
    ICHECKNMI_tas.remove_coord('grid_longitude')
    ICHECKNMI_tas.add_dim_coord(lats, 1)
    ICHECKNMI_tas.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(ICHECMPI_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECMPI_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees') 
    
    ICHECMPI_tas.remove_coord('latitude')
    ICHECMPI_tas.remove_coord('longitude')
    ICHECMPI_tas.remove_coord('grid_latitude')
    ICHECMPI_tas.remove_coord('grid_longitude')
    ICHECMPI_tas.add_dim_coord(lats, 1)
    ICHECMPI_tas.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(ICHECSMHI_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECSMHI_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    ICHECSMHI_tas.remove_coord('latitude')
    ICHECSMHI_tas.remove_coord('longitude')
    ICHECSMHI_tas.remove_coord('grid_latitude')
    ICHECSMHI_tas.remove_coord('grid_longitude')
    ICHECSMHI_tas.add_dim_coord(lats, 1)
    ICHECSMHI_tas.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(IPSL_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = IPSL_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    IPSL_tas.remove_coord('latitude')
    IPSL_tas.remove_coord('longitude')
    IPSL_tas.remove_coord('grid_latitude')
    IPSL_tas.remove_coord('grid_longitude')
    IPSL_tas.add_dim_coord(lats, 1)
    IPSL_tas.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MIROC_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MIROC_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
                                
    MIROC_tas.remove_coord('latitude')
    MIROC_tas.remove_coord('longitude')
    MIROC_tas.remove_coord('grid_latitude')
    MIROC_tas.remove_coord('grid_longitude')
    MIROC_tas.add_dim_coord(lats, 1)
    MIROC_tas.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MOHCCCLM_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MOHCCCLM_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MOHCCCLM_tas.remove_coord('latitude')
    MOHCCCLM_tas.remove_coord('longitude')
    MOHCCCLM_tas.remove_coord('grid_latitude')
    MOHCCCLM_tas.remove_coord('grid_longitude')
    MOHCCCLM_tas.add_dim_coord(lats, 1)
    MOHCCCLM_tas.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MOHCKNMI_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MOHCKNMI_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MOHCKNMI_tas.remove_coord('latitude')
    MOHCKNMI_tas.remove_coord('longitude')
    MOHCKNMI_tas.remove_coord('grid_latitude')
    MOHCKNMI_tas.remove_coord('grid_longitude')
    MOHCKNMI_tas.add_dim_coord(lats, 1)
    MOHCKNMI_tas.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MOHCSMHI_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MOHCSMHI_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MOHCSMHI_tas.remove_coord('latitude')
    MOHCSMHI_tas.remove_coord('longitude')
    MOHCSMHI_tas.remove_coord('grid_latitude')
    MOHCSMHI_tas.remove_coord('grid_longitude')
    MOHCSMHI_tas.add_dim_coord(lats, 1)
    MOHCSMHI_tas.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MPICCLM_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MPICCLM_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MPICCLM_tas.remove_coord('latitude')
    MPICCLM_tas.remove_coord('longitude')
    MPICCLM_tas.remove_coord('grid_latitude')
    MPICCLM_tas.remove_coord('grid_longitude')
    MPICCLM_tas.add_dim_coord(lats, 1)
    MPICCLM_tas.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MPIREMO_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MPIREMO_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MPIREMO_tas.remove_coord('latitude')
    MPIREMO_tas.remove_coord('longitude')
    MPIREMO_tas.remove_coord('grid_latitude')
    MPIREMO_tas.remove_coord('grid_longitude')
    MPIREMO_tas.add_dim_coord(lats, 1)
    MPIREMO_tas.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MPISMHI_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MPISMHI_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MPISMHI_tas.remove_coord('latitude')
    MPISMHI_tas.remove_coord('longitude')
    MPISMHI_tas.remove_coord('grid_latitude')
    MPISMHI_tas.remove_coord('grid_longitude')
    MPISMHI_tas.add_dim_coord(lats, 1)
    MPISMHI_tas.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(NCCDMI_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = NCCDMI_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    NCCDMI_tas.remove_coord('latitude')
    NCCDMI_tas.remove_coord('longitude')
    NCCDMI_tas.remove_coord('grid_latitude')
    NCCDMI_tas.remove_coord('grid_longitude')
    NCCDMI_tas.add_dim_coord(lats, 1)
    NCCDMI_tas.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(NCCSMHI_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = NCCSMHI_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    NCCSMHI_tas.remove_coord('latitude')
    NCCSMHI_tas.remove_coord('longitude')
    NCCSMHI_tas.remove_coord('grid_latitude')
    NCCSMHI_tas.remove_coord('grid_longitude')
    NCCSMHI_tas.add_dim_coord(lats, 1)
    NCCSMHI_tas.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(NOAA_tas.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = NOAA_tas.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    NOAA_tas.remove_coord('latitude')
    NOAA_tas.remove_coord('longitude')
    NOAA_tas.remove_coord('grid_latitude')
    NOAA_tas.remove_coord('grid_longitude')
    NOAA_tas.add_dim_coord(lats, 1)
    NOAA_tas.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(CCCma_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CCCma_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
                                
    CCCma_tasmin.remove_coord('latitude')
    CCCma_tasmin.remove_coord('longitude')
    CCCma_tasmin.remove_coord('grid_latitude')
    CCCma_tasmin.remove_coord('grid_longitude')
    CCCma_tasmin.add_dim_coord(lats, 1)
    CCCma_tasmin.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(CLMcom_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CLMcom_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    CLMcom_tasmin.remove_coord('latitude')
    CLMcom_tasmin.remove_coord('longitude')
    CLMcom_tasmin.remove_coord('grid_latitude')
    CLMcom_tasmin.remove_coord('grid_longitude')
    CLMcom_tasmin.add_dim_coord(lats, 1)
    CLMcom_tasmin.add_dim_coord(lons, 2)  
    
    lats = iris.coords.DimCoord(DMI_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = DMI_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    DMI_tasmin.remove_coord('latitude')
    DMI_tasmin.remove_coord('longitude')
    DMI_tasmin.remove_coord('grid_latitude')
    DMI_tasmin.remove_coord('grid_longitude')
    DMI_tasmin.add_dim_coord(lats, 1)
    DMI_tasmin.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(KNMI_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = KNMI_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
        
    KNMI_tasmin.remove_coord('latitude')
    KNMI_tasmin.remove_coord('longitude')
    KNMI_tasmin.remove_coord('grid_latitude')
    KNMI_tasmin.remove_coord('grid_longitude')
    KNMI_tasmin.add_dim_coord(lats, 1)
    KNMI_tasmin.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(MPIE_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MPIE_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
   
    MPIE_tasmin.remove_coord('latitude')
    MPIE_tasmin.remove_coord('longitude')
    MPIE_tasmin.remove_coord('grid_latitude')
    MPIE_tasmin.remove_coord('grid_longitude')
    MPIE_tasmin.add_dim_coord(lats, 1)
    MPIE_tasmin.add_dim_coord(lons, 2)

    lats = iris.coords.DimCoord(SMHI_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = SMHI_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    SMHI_tasmin.remove_coord('latitude')
    SMHI_tasmin.remove_coord('longitude')
    SMHI_tasmin.remove_coord('grid_latitude')
    SMHI_tasmin.remove_coord('grid_longitude')
    SMHI_tasmin.add_dim_coord(lats, 1)
    SMHI_tasmin.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(CCCmaCanRCM_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CCCmaCanRCM_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
                              
    CCCmaCanRCM_tasmin.remove_coord('latitude')
    CCCmaCanRCM_tasmin.remove_coord('longitude')
    CCCmaCanRCM_tasmin.remove_coord('grid_latitude')
    CCCmaCanRCM_tasmin.remove_coord('grid_longitude')
    CCCmaCanRCM_tasmin.add_dim_coord(lats, 1)
    CCCmaCanRCM_tasmin.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(CCCmaSMHI_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CCCmaSMHI_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    CCCmaSMHI_tasmin.remove_coord('latitude')
    CCCmaSMHI_tasmin.remove_coord('longitude')
    CCCmaSMHI_tasmin.remove_coord('grid_latitude')
    CCCmaSMHI_tasmin.remove_coord('grid_longitude')
    CCCmaSMHI_tasmin.add_dim_coord(lats, 1)
    CCCmaSMHI_tasmin.add_dim_coord(lons, 2)  
    
    lats = iris.coords.DimCoord(CNRM_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CNRM_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
                                
    CNRM_tasmin.remove_coord('latitude')
    CNRM_tasmin.remove_coord('longitude')
    CNRM_tasmin.remove_coord('grid_latitude')
    CNRM_tasmin.remove_coord('grid_longitude')
    CNRM_tasmin.add_dim_coord(lats, 1)
    CNRM_tasmin.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(CNRMSMHI_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CNRMSMHI_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    CNRMSMHI_tasmin.remove_coord('latitude')
    CNRMSMHI_tasmin.remove_coord('longitude')
    CNRMSMHI_tasmin.remove_coord('grid_latitude')
    CNRMSMHI_tasmin.remove_coord('grid_longitude')
    CNRMSMHI_tasmin.add_dim_coord(lats, 1)
    CNRMSMHI_tasmin.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(CSIRO_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CSIRO_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    CSIRO_tasmin.remove_coord('latitude')
    CSIRO_tasmin.remove_coord('longitude')
    CSIRO_tasmin.remove_coord('grid_latitude')
    CSIRO_tasmin.remove_coord('grid_longitude')
    CSIRO_tasmin.add_dim_coord(lats, 1)
    CSIRO_tasmin.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(ICHECDMI_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECDMI_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')    
    
    ICHECDMI_tasmin.remove_coord('latitude')
    ICHECDMI_tasmin.remove_coord('longitude')
    ICHECDMI_tasmin.remove_coord('grid_latitude')
    ICHECDMI_tasmin.remove_coord('grid_longitude')
    ICHECDMI_tasmin.add_dim_coord(lats, 1)
    ICHECDMI_tasmin.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(ICHECCCLM_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECCCLM_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees') 
    
    ICHECCCLM_tasmin.remove_coord('latitude')
    ICHECCCLM_tasmin.remove_coord('longitude')
    ICHECCCLM_tasmin.remove_coord('grid_latitude')
    ICHECCCLM_tasmin.remove_coord('grid_longitude')
    ICHECCCLM_tasmin.add_dim_coord(lats, 1)
    ICHECCCLM_tasmin.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(ICHECKNMI_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECKNMI_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees') 
    
    ICHECKNMI_tasmin.remove_coord('latitude')
    ICHECKNMI_tasmin.remove_coord('longitude')
    ICHECKNMI_tasmin.remove_coord('grid_latitude')
    ICHECKNMI_tasmin.remove_coord('grid_longitude')
    ICHECKNMI_tasmin.add_dim_coord(lats, 1)
    ICHECKNMI_tasmin.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(ICHECMPI_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECMPI_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees') 
    
    ICHECMPI_tasmin.remove_coord('latitude')
    ICHECMPI_tasmin.remove_coord('longitude')
    ICHECMPI_tasmin.remove_coord('grid_latitude')
    ICHECMPI_tasmin.remove_coord('grid_longitude')
    ICHECMPI_tasmin.add_dim_coord(lats, 1)
    ICHECMPI_tasmin.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(ICHECSMHI_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECSMHI_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    ICHECSMHI_tasmin.remove_coord('latitude')
    ICHECSMHI_tasmin.remove_coord('longitude')
    ICHECSMHI_tasmin.remove_coord('grid_latitude')
    ICHECSMHI_tasmin.remove_coord('grid_longitude')
    ICHECSMHI_tasmin.add_dim_coord(lats, 1)
    ICHECSMHI_tasmin.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(IPSL_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = IPSL_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    IPSL_tasmin.remove_coord('latitude')
    IPSL_tasmin.remove_coord('longitude')
    IPSL_tasmin.remove_coord('grid_latitude')
    IPSL_tasmin.remove_coord('grid_longitude')
    IPSL_tasmin.add_dim_coord(lats, 1)
    IPSL_tasmin.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MIROC_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MIROC_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
                                
    MIROC_tasmin.remove_coord('latitude')
    MIROC_tasmin.remove_coord('longitude')
    MIROC_tasmin.remove_coord('grid_latitude')
    MIROC_tasmin.remove_coord('grid_longitude')
    MIROC_tasmin.add_dim_coord(lats, 1)
    MIROC_tasmin.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MOHCCCLM_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MOHCCCLM_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MOHCCCLM_tasmin.remove_coord('latitude')
    MOHCCCLM_tasmin.remove_coord('longitude')
    MOHCCCLM_tasmin.remove_coord('grid_latitude')
    MOHCCCLM_tasmin.remove_coord('grid_longitude')
    MOHCCCLM_tasmin.add_dim_coord(lats, 1)
    MOHCCCLM_tasmin.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MOHCKNMI_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MOHCKNMI_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MOHCKNMI_tasmin.remove_coord('latitude')
    MOHCKNMI_tasmin.remove_coord('longitude')
    MOHCKNMI_tasmin.remove_coord('grid_latitude')
    MOHCKNMI_tasmin.remove_coord('grid_longitude')
    MOHCKNMI_tasmin.add_dim_coord(lats, 1)
    MOHCKNMI_tasmin.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MOHCSMHI_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MOHCSMHI_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MOHCSMHI_tasmin.remove_coord('latitude')
    MOHCSMHI_tasmin.remove_coord('longitude')
    MOHCSMHI_tasmin.remove_coord('grid_latitude')
    MOHCSMHI_tasmin.remove_coord('grid_longitude')
    MOHCSMHI_tasmin.add_dim_coord(lats, 1)
    MOHCSMHI_tasmin.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MPICCLM_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MPICCLM_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MPICCLM_tasmin.remove_coord('latitude')
    MPICCLM_tasmin.remove_coord('longitude')
    MPICCLM_tasmin.remove_coord('grid_latitude')
    MPICCLM_tasmin.remove_coord('grid_longitude')
    MPICCLM_tasmin.add_dim_coord(lats, 1)
    MPICCLM_tasmin.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MPIREMO_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MPIREMO_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MPIREMO_tasmin.remove_coord('latitude')
    MPIREMO_tasmin.remove_coord('longitude')
    MPIREMO_tasmin.remove_coord('grid_latitude')
    MPIREMO_tasmin.remove_coord('grid_longitude')
    MPIREMO_tasmin.add_dim_coord(lats, 1)
    MPIREMO_tasmin.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MPISMHI_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MPISMHI_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MPISMHI_tasmin.remove_coord('latitude')
    MPISMHI_tasmin.remove_coord('longitude')
    MPISMHI_tasmin.remove_coord('grid_latitude')
    MPISMHI_tasmin.remove_coord('grid_longitude')
    MPISMHI_tasmin.add_dim_coord(lats, 1)
    MPISMHI_tasmin.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(NCCDMI_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = NCCDMI_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    NCCDMI_tasmin.remove_coord('latitude')
    NCCDMI_tasmin.remove_coord('longitude')
    NCCDMI_tasmin.remove_coord('grid_latitude')
    NCCDMI_tasmin.remove_coord('grid_longitude')
    NCCDMI_tasmin.add_dim_coord(lats, 1)
    NCCDMI_tasmin.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(NCCSMHI_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = NCCSMHI_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    NCCSMHI_tasmin.remove_coord('latitude')
    NCCSMHI_tasmin.remove_coord('longitude')
    NCCSMHI_tasmin.remove_coord('grid_latitude')
    NCCSMHI_tasmin.remove_coord('grid_longitude')
    NCCSMHI_tasmin.add_dim_coord(lats, 1)
    NCCSMHI_tasmin.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(NOAA_tasmin.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = NOAA_tasmin.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    NOAA_tasmin.remove_coord('latitude')
    NOAA_tasmin.remove_coord('longitude')
    NOAA_tasmin.remove_coord('grid_latitude')
    NOAA_tasmin.remove_coord('grid_longitude')
    NOAA_tasmin.add_dim_coord(lats, 1)
    NOAA_tasmin.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(CCCma_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CCCma_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
                                
    CCCma_tasmax.remove_coord('latitude')
    CCCma_tasmax.remove_coord('longitude')
    CCCma_tasmax.remove_coord('grid_latitude')
    CCCma_tasmax.remove_coord('grid_longitude')
    CCCma_tasmax.add_dim_coord(lats, 1)
    CCCma_tasmax.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(CLMcom_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CLMcom_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    CLMcom_tasmax.remove_coord('latitude')
    CLMcom_tasmax.remove_coord('longitude')
    CLMcom_tasmax.remove_coord('grid_latitude')
    CLMcom_tasmax.remove_coord('grid_longitude')
    CLMcom_tasmax.add_dim_coord(lats, 1)
    CLMcom_tasmax.add_dim_coord(lons, 2)  
    
    lats = iris.coords.DimCoord(DMI_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = DMI_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    DMI_tasmax.remove_coord('latitude')
    DMI_tasmax.remove_coord('longitude')
    DMI_tasmax.remove_coord('grid_latitude')
    DMI_tasmax.remove_coord('grid_longitude')
    DMI_tasmax.add_dim_coord(lats, 1)
    DMI_tasmax.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(KNMI_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = KNMI_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
        
    KNMI_tasmax.remove_coord('latitude')
    KNMI_tasmax.remove_coord('longitude')
    KNMI_tasmax.remove_coord('grid_latitude')
    KNMI_tasmax.remove_coord('grid_longitude')
    KNMI_tasmax.add_dim_coord(lats, 1)
    KNMI_tasmax.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(MPIE_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MPIE_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
   
    MPIE_tasmax.remove_coord('latitude')
    MPIE_tasmax.remove_coord('longitude')
    MPIE_tasmax.remove_coord('grid_latitude')
    MPIE_tasmax.remove_coord('grid_longitude')
    MPIE_tasmax.add_dim_coord(lats, 1)
    MPIE_tasmax.add_dim_coord(lons, 2)

    lats = iris.coords.DimCoord(SMHI_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = SMHI_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    SMHI_tasmax.remove_coord('latitude')
    SMHI_tasmax.remove_coord('longitude')
    SMHI_tasmax.remove_coord('grid_latitude')
    SMHI_tasmax.remove_coord('grid_longitude')
    SMHI_tasmax.add_dim_coord(lats, 1)
    SMHI_tasmax.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(CCCmaCanRCM_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CCCmaCanRCM_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
                              
    CCCmaCanRCM_tasmax.remove_coord('latitude')
    CCCmaCanRCM_tasmax.remove_coord('longitude')
    CCCmaCanRCM_tasmax.remove_coord('grid_latitude')
    CCCmaCanRCM_tasmax.remove_coord('grid_longitude')
    CCCmaCanRCM_tasmax.add_dim_coord(lats, 1)
    CCCmaCanRCM_tasmax.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(CCCmaSMHI_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CCCmaSMHI_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    CCCmaSMHI_tasmax.remove_coord('latitude')
    CCCmaSMHI_tasmax.remove_coord('longitude')
    CCCmaSMHI_tasmax.remove_coord('grid_latitude')
    CCCmaSMHI_tasmax.remove_coord('grid_longitude')
    CCCmaSMHI_tasmax.add_dim_coord(lats, 1)
    CCCmaSMHI_tasmax.add_dim_coord(lons, 2)  
    
    lats = iris.coords.DimCoord(CNRM_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CNRM_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
                                
    CNRM_tasmax.remove_coord('latitude')
    CNRM_tasmax.remove_coord('longitude')
    CNRM_tasmax.remove_coord('grid_latitude')
    CNRM_tasmax.remove_coord('grid_longitude')
    CNRM_tasmax.add_dim_coord(lats, 1)
    CNRM_tasmax.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(CNRMSMHI_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CNRMSMHI_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    CNRMSMHI_tasmax.remove_coord('latitude')
    CNRMSMHI_tasmax.remove_coord('longitude')
    CNRMSMHI_tasmax.remove_coord('grid_latitude')
    CNRMSMHI_tasmax.remove_coord('grid_longitude')
    CNRMSMHI_tasmax.add_dim_coord(lats, 1)
    CNRMSMHI_tasmax.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(CSIRO_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = CSIRO_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    CSIRO_tasmax.remove_coord('latitude')
    CSIRO_tasmax.remove_coord('longitude')
    CSIRO_tasmax.remove_coord('grid_latitude')
    CSIRO_tasmax.remove_coord('grid_longitude')
    CSIRO_tasmax.add_dim_coord(lats, 1)
    CSIRO_tasmax.add_dim_coord(lons, 2) 
    
    lats = iris.coords.DimCoord(ICHECDMI_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECDMI_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')    
    
    ICHECDMI_tasmax.remove_coord('latitude')
    ICHECDMI_tasmax.remove_coord('longitude')
    ICHECDMI_tasmax.remove_coord('grid_latitude')
    ICHECDMI_tasmax.remove_coord('grid_longitude')
    ICHECDMI_tasmax.add_dim_coord(lats, 1)
    ICHECDMI_tasmax.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(ICHECCCLM_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECCCLM_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees') 
    
    ICHECCCLM_tasmax.remove_coord('latitude')
    ICHECCCLM_tasmax.remove_coord('longitude')
    ICHECCCLM_tasmax.remove_coord('grid_latitude')
    ICHECCCLM_tasmax.remove_coord('grid_longitude')
    ICHECCCLM_tasmax.add_dim_coord(lats, 1)
    ICHECCCLM_tasmax.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(ICHECKNMI_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECKNMI_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees') 
    
    ICHECKNMI_tasmax.remove_coord('latitude')
    ICHECKNMI_tasmax.remove_coord('longitude')
    ICHECKNMI_tasmax.remove_coord('grid_latitude')
    ICHECKNMI_tasmax.remove_coord('grid_longitude')
    ICHECKNMI_tasmax.add_dim_coord(lats, 1)
    ICHECKNMI_tasmax.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(ICHECMPI_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECMPI_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees') 
    
    ICHECMPI_tasmax.remove_coord('latitude')
    ICHECMPI_tasmax.remove_coord('longitude')
    ICHECMPI_tasmax.remove_coord('grid_latitude')
    ICHECMPI_tasmax.remove_coord('grid_longitude')
    ICHECMPI_tasmax.add_dim_coord(lats, 1)
    ICHECMPI_tasmax.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(ICHECSMHI_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = ICHECSMHI_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    ICHECSMHI_tasmax.remove_coord('latitude')
    ICHECSMHI_tasmax.remove_coord('longitude')
    ICHECSMHI_tasmax.remove_coord('grid_latitude')
    ICHECSMHI_tasmax.remove_coord('grid_longitude')
    ICHECSMHI_tasmax.add_dim_coord(lats, 1)
    ICHECSMHI_tasmax.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(IPSL_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = IPSL_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    IPSL_tasmax.remove_coord('latitude')
    IPSL_tasmax.remove_coord('longitude')
    IPSL_tasmax.remove_coord('grid_latitude')
    IPSL_tasmax.remove_coord('grid_longitude')
    IPSL_tasmax.add_dim_coord(lats, 1)
    IPSL_tasmax.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MIROC_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MIROC_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
                                
    MIROC_tasmax.remove_coord('latitude')
    MIROC_tasmax.remove_coord('longitude')
    MIROC_tasmax.remove_coord('grid_latitude')
    MIROC_tasmax.remove_coord('grid_longitude')
    MIROC_tasmax.add_dim_coord(lats, 1)
    MIROC_tasmax.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MOHCCCLM_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MOHCCCLM_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MOHCCCLM_tasmax.remove_coord('latitude')
    MOHCCCLM_tasmax.remove_coord('longitude')
    MOHCCCLM_tasmax.remove_coord('grid_latitude')
    MOHCCCLM_tasmax.remove_coord('grid_longitude')
    MOHCCCLM_tasmax.add_dim_coord(lats, 1)
    MOHCCCLM_tasmax.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MOHCKNMI_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MOHCKNMI_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MOHCKNMI_tasmax.remove_coord('latitude')
    MOHCKNMI_tasmax.remove_coord('longitude')
    MOHCKNMI_tasmax.remove_coord('grid_latitude')
    MOHCKNMI_tasmax.remove_coord('grid_longitude')
    MOHCKNMI_tasmax.add_dim_coord(lats, 1)
    MOHCKNMI_tasmax.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MOHCSMHI_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MOHCSMHI_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MOHCSMHI_tasmax.remove_coord('latitude')
    MOHCSMHI_tasmax.remove_coord('longitude')
    MOHCSMHI_tasmax.remove_coord('grid_latitude')
    MOHCSMHI_tasmax.remove_coord('grid_longitude')
    MOHCSMHI_tasmax.add_dim_coord(lats, 1)
    MOHCSMHI_tasmax.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MPICCLM_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MPICCLM_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MPICCLM_tasmax.remove_coord('latitude')
    MPICCLM_tasmax.remove_coord('longitude')
    MPICCLM_tasmax.remove_coord('grid_latitude')
    MPICCLM_tasmax.remove_coord('grid_longitude')
    MPICCLM_tasmax.add_dim_coord(lats, 1)
    MPICCLM_tasmax.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MPIREMO_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MPIREMO_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MPIREMO_tasmax.remove_coord('latitude')
    MPIREMO_tasmax.remove_coord('longitude')
    MPIREMO_tasmax.remove_coord('grid_latitude')
    MPIREMO_tasmax.remove_coord('grid_longitude')
    MPIREMO_tasmax.add_dim_coord(lats, 1)
    MPIREMO_tasmax.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(MPISMHI_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = MPISMHI_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    MPISMHI_tasmax.remove_coord('latitude')
    MPISMHI_tasmax.remove_coord('longitude')
    MPISMHI_tasmax.remove_coord('grid_latitude')
    MPISMHI_tasmax.remove_coord('grid_longitude')
    MPISMHI_tasmax.add_dim_coord(lats, 1)
    MPISMHI_tasmax.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(NCCDMI_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = NCCDMI_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    NCCDMI_tasmax.remove_coord('latitude')
    NCCDMI_tasmax.remove_coord('longitude')
    NCCDMI_tasmax.remove_coord('grid_latitude')
    NCCDMI_tasmax.remove_coord('grid_longitude')
    NCCDMI_tasmax.add_dim_coord(lats, 1)
    NCCDMI_tasmax.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(NCCSMHI_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = NCCSMHI_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    NCCSMHI_tasmax.remove_coord('latitude')
    NCCSMHI_tasmax.remove_coord('longitude')
    NCCSMHI_tasmax.remove_coord('grid_latitude')
    NCCSMHI_tasmax.remove_coord('grid_longitude')
    NCCSMHI_tasmax.add_dim_coord(lats, 1)
    NCCSMHI_tasmax.add_dim_coord(lons, 2)
    
    lats = iris.coords.DimCoord(NOAA_tasmax.coord('latitude').points[:,0], standard_name='latitude', units='degrees')
    lons = NOAA_tasmax.coord('longitude').points[0]
    for i in range(len(lons)):
        if lons[i]>100.:
            lons[i] = lons[i]-360.
    lons = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    NOAA_tasmax.remove_coord('latitude')
    NOAA_tasmax.remove_coord('longitude')
    NOAA_tasmax.remove_coord('grid_latitude')
    NOAA_tasmax.remove_coord('grid_longitude')
    NOAA_tasmax.add_dim_coord(lats, 1)
    NOAA_tasmax.add_dim_coord(lons, 2)
    
    
    #-------------------------------------------------------------------------
    #PART 2: LOAD AND FORMAT OBSERVED DATA
    #bring in all the files we need and give them a name
    CRU_pr= '/exports/csce/datastore/geos/users/s0899345/Actual_Data/cru_ts4.00.1901.2015.pre.dat.nc'
    UDel_pr= '/exports/csce/datastore/geos/users/s0899345/Actual_Data/UDel_precip.mon.total.v401.nc'
    GPCC_pr= '/exports/csce/datastore/geos/users/s0899345/Actual_Data/ESRL_PSD_GPCC_precip.mon.combined.total.v7.nc'
    
    CRU_tas= '/exports/csce/datastore/geos/users/s0899345/Actual_Data/cru_ts4.00.1901.2015.tmp.dat.nc'
    UDel_tas= '/exports/csce/datastore/geos/users/s0899345/Actual_Data/UDel_air.mon.mean.v401.nc'
    
    CRU_tasmin= '/exports/csce/datastore/geos/users/s0899345/Actual_Data/cru_ts4.01.1901.2016.tmn.dat.nc'
    
    CRU_tasmax= '/exports/csce/datastore/geos/users/s0899345/Actual_Data/cru_ts4.01.1901.2016.tmx.dat.nc'
    
    #Load exactly one cube from given file
    CRU_pr = iris.load_cube(CRU_pr, 'precipitation')
    UDel_pr = iris.load_cube(UDel_pr)
    GPCC_pr = iris.load_cube(GPCC_pr)
    
    CRU_tas = iris.load_cube(CRU_tas, 'near-surface temperature')
    UDel_tas = iris.load_cube(UDel_tas)
    
    CRU_tasmin = iris.load_cube(CRU_tasmin, 'near-surface temperature minimum')
    
    CRU_tasmax = iris.load_cube(CRU_tasmax, 'near-surface temperature maximum')
    
    
    #-------------------------------------------------------------------------
    #PART 3: FORMAT FILES TO BE PLOT SPECIFIC
    #fix EARTH3 time units as they differ from all other models
    t_coord=EARTH3_pr.coord('time')
    t_unit = t_coord.attributes['invalid_units']
    timestep, _, t_fmt_str = t_unit.split(' ')
    new_t_unit_str= '{} since 1850-01-01 00:00:00'.format(timestep) 
    new_t_unit = cf_units.Unit(new_t_unit_str, calendar=cf_units.CALENDAR_STANDARD)
    new_datetimes = [datetime.datetime.strptime(str(dt), t_fmt_str) for dt in t_coord.points]
    new_dt_points = [new_t_unit.date2num(new_dt) for new_dt in new_datetimes]
    new_t_coord = iris.coords.DimCoord(new_dt_points, standard_name='time', units=new_t_unit)
    t_coord_dim = EARTH3_pr.coord_dims('time')
    EARTH3_pr.remove_coord('time')
    EARTH3_pr.add_dim_coord(new_t_coord,t_coord_dim)
    
    t_coord=EARTH3_tas.coord('time')
#    t_unit = t_coord.attributes['invalid_units']
    timestep, _, t_fmt_str = t_unit.split(' ')
    new_t_unit_str= '{} since 1850-01-01 00:00:00'.format(timestep) 
    new_t_unit = cf_units.Unit(new_t_unit_str, calendar=cf_units.CALENDAR_STANDARD)
    new_datetimes = [datetime.datetime.strptime(str(dt), t_fmt_str) for dt in t_coord.points]
    new_dt_points = [new_t_unit.date2num(new_dt) for new_dt in new_datetimes]
    new_t_coord = iris.coords.DimCoord(new_dt_points, standard_name='time', units=new_t_unit)
    t_coord_dim = EARTH3_tas.coord_dims('time')
    EARTH3_tas.remove_coord('time')
    EARTH3_tas.add_dim_coord(new_t_coord,t_coord_dim)
                                       
    t_coord=EARTH3_tasmin.coord('time')
#    t_unit = t_coord.attributes['invalid_units']
    timestep, _, t_fmt_str = t_unit.split(' ')
    new_t_unit_str= '{} since 1850-01-01 00:00:00'.format(timestep) 
    new_t_unit = cf_units.Unit(new_t_unit_str, calendar=cf_units.CALENDAR_STANDARD)
    new_datetimes = [datetime.datetime.strptime(str(dt), t_fmt_str) for dt in t_coord.points]
    new_dt_points = [new_t_unit.date2num(new_dt) for new_dt in new_datetimes]
    new_t_coord = iris.coords.DimCoord(new_dt_points, standard_name='time', units=new_t_unit)
    t_coord_dim = EARTH3_tasmin.coord_dims('time')
    EARTH3_tasmin.remove_coord('time')
    EARTH3_tasmin.add_dim_coord(new_t_coord,t_coord_dim)
                                       
    t_coord=EARTH3_tasmax.coord('time')
#    t_unit = t_coord.attributes['invalid_units']
    timestep, _, t_fmt_str = t_unit.split(' ')
    new_t_unit_str= '{} since 1850-01-01 00:00:00'.format(timestep) 
    new_t_unit = cf_units.Unit(new_t_unit_str, calendar=cf_units.CALENDAR_STANDARD)
    new_datetimes = [datetime.datetime.strptime(str(dt), t_fmt_str) for dt in t_coord.points]
    new_dt_points = [new_t_unit.date2num(new_dt) for new_dt in new_datetimes]
    new_t_coord = iris.coords.DimCoord(new_dt_points, standard_name='time', units=new_t_unit)
    t_coord_dim = EARTH3_tasmax.coord_dims('time')
    EARTH3_tasmax.remove_coord('time')
    EARTH3_tasmax.add_dim_coord(new_t_coord,t_coord_dim)

    
    #regrid all models to have same latitude and longitude system, all regridded to model with lowest resolution
    CCCma_pr = CCCma_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    CLMcom_pr = CLMcom_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    DMI_pr = DMI_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    KNMI_pr = KNMI_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    MPIE_pr = MPIE_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    SMHI_pr = SMHI_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    
    CCCmaCanRCM_pr = CCCmaCanRCM_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    CCCmaSMHI_pr = CCCmaSMHI_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    CNRM_pr = CNRM_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    CNRMSMHI_pr = CNRMSMHI_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    CSIRO_pr = CSIRO_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    ICHECDMI_pr = ICHECDMI_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    ICHECCCLM_pr = ICHECCCLM_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    ICHECKNMI_pr = ICHECKNMI_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    ICHECMPI_pr = ICHECMPI_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    ICHECSMHI_pr = ICHECSMHI_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    IPSL_pr = IPSL_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    MIROC_pr = MIROC_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    MOHCCCLM_pr = MOHCCCLM_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    MOHCKNMI_pr = MOHCKNMI_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    MOHCSMHI_pr = MOHCSMHI_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    MPICCLM_pr = MPICCLM_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    MPIREMO_pr = MPIREMO_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    MPISMHI_pr = MPISMHI_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    NCCDMI_pr = NCCDMI_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    NCCSMHI_pr = NCCSMHI_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    NOAA_pr = NOAA_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    
    #CanESM2 = CanESM2_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    CNRMG_pr = CNRMG_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    MK3_pr = MK3_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    EARTH_pr = EARTH_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    EARTH3_pr = EARTH3_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    GFDL_pr = GFDL_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    HadGEM2_pr = HadGEM2_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    IPSLG_pr = IPSLG_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    MIROCG_pr = MIROCG_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    MPI_pr = MPI_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    NorESM1_pr = NorESM1_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    
    CRU_pr = CRU_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    UDel_pr = UDel_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    GPCC_pr = GPCC_pr.regrid(CanESM2_pr, iris.analysis.Linear())
    
    CCCma_tas = CCCma_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    CLMcom_tas = CLMcom_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    DMI_tas = DMI_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    KNMI_tas = KNMI_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    MPIE_tas = MPIE_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    SMHI_tas = SMHI_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    
    CCCmaCanRCM_tas = CCCmaCanRCM_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    CCCmaSMHI_tas = CCCmaSMHI_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    CNRM_tas = CNRM_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    CNRMSMHI_tas = CNRMSMHI_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    CSIRO_tas = CSIRO_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    ICHECDMI_tas = ICHECDMI_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    ICHECCCLM_tas = ICHECCCLM_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    ICHECKNMI_tas = ICHECKNMI_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    ICHECMPI_tas = ICHECMPI_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    ICHECSMHI_tas = ICHECSMHI_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    IPSL_tas = IPSL_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    MIROC_tas = MIROC_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    MOHCCCLM_tas = MOHCCCLM_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    MOHCKNMI_tas = MOHCKNMI_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    MOHCSMHI_tas = MOHCSMHI_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    MPICCLM_tas = MPICCLM_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    MPIREMO_tas = MPIREMO_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    MPISMHI_tas = MPISMHI_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    NCCDMI_tas = NCCDMI_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    NCCSMHI_tas = NCCSMHI_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    NOAA_tas = NOAA_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    
    #CanESM2 = CanESM2_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    CNRMG_tas = CNRMG_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    MK3_tas = MK3_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    EARTH_tas = EARTH_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    EARTH3_tas = EARTH3_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    GFDL_tas = GFDL_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    HadGEM2_tas = HadGEM2_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    IPSLG_tas = IPSLG_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    MIROCG_tas = MIROCG_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    MPI_tas = MPI_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    NorESM1_tas = NorESM1_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    
    CRU_tas = CRU_tas.regrid(CanESM2_tas, iris.analysis.Linear())
    UDel_tas = UDel_tas.regrid(CanESM2_tas, iris.analysis.Linear())
   
    CCCma_tasmin = CCCma_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    CLMcom_tasmin = CLMcom_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    DMI_tasmin = DMI_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    KNMI_tasmin = KNMI_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    MPIE_tasmin = MPIE_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    SMHI_tasmin = SMHI_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    
    CCCmaCanRCM_tasmin = CCCmaCanRCM_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    CCCmaSMHI_tasmin = CCCmaSMHI_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    CNRM_tasmin = CNRM_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    CNRMSMHI_tasmin = CNRMSMHI_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    CSIRO_tasmin = CSIRO_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    ICHECDMI_tasmin = ICHECDMI_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    ICHECCCLM_tasmin = ICHECCCLM_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    ICHECKNMI_tasmin = ICHECKNMI_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    ICHECMPI_tasmin = ICHECMPI_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    ICHECSMHI_tasmin = ICHECSMHI_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    IPSL_tasmin = IPSL_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    MIROC_tasmin = MIROC_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    MOHCCCLM_tasmin = MOHCCCLM_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    MOHCKNMI_tasmin = MOHCKNMI_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    MOHCSMHI_tasmin = MOHCSMHI_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    MPICCLM_tasmin = MPICCLM_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    MPIREMO_tasmin = MPIREMO_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    MPISMHI_tasmin = MPISMHI_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    NCCDMI_tasmin = NCCDMI_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    NCCSMHI_tasmin = NCCSMHI_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    NOAA_tasmin = NOAA_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    
    #CanESM2 = CanESM2_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    CNRMG_tasmin = CNRMG_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    MK3_tasmin = MK3_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    EARTH_tasmin = EARTH_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    EARTH3_tasmin = EARTH3_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    GFDL_tasmin = GFDL_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    HadGEM2_tasmin = HadGEM2_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    IPSLG_tasmin = IPSLG_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    MIROCG_tasmin = MIROCG_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    MPI_tasmin = MPI_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    NorESM1_tasmin = NorESM1_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
    
    CRU_tasmin = CRU_tasmin.regrid(CanESM2_tasmin, iris.analysis.Linear())
  
    CCCma_tasmax = CCCma_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    CLMcom_tasmax = CLMcom_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    DMI_tasmax = DMI_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    KNMI_tasmax = KNMI_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    MPIE_tasmax = MPIE_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    SMHI_tasmax = SMHI_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    
    CCCmaCanRCM_tasmax = CCCmaCanRCM_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    CCCmaSMHI_tasmax = CCCmaSMHI_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    CNRM_tasmax = CNRM_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    CNRMSMHI_tasmax = CNRMSMHI_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    CSIRO_tasmax = CSIRO_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    ICHECDMI_tasmax = ICHECDMI_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    ICHECCCLM_tasmax = ICHECCCLM_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    ICHECKNMI_tasmax = ICHECKNMI_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    ICHECMPI_tasmax = ICHECMPI_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    ICHECSMHI_tasmax = ICHECSMHI_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    IPSL_tasmax = IPSL_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    MIROC_tasmax = MIROC_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    MOHCCCLM_tasmax = MOHCCCLM_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    MOHCKNMI_tasmax = MOHCKNMI_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    MOHCSMHI_tasmax = MOHCSMHI_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    MPICCLM_tasmax = MPICCLM_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    MPIREMO_tasmax = MPIREMO_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    MPISMHI_tasmax = MPISMHI_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    NCCDMI_tasmax = NCCDMI_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    NCCSMHI_tasmax = NCCSMHI_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    NOAA_tasmax = NOAA_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    
    #CanESM2 = CanESM2_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    CNRMG_tasmax = CNRMG_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    MK3_tasmax = MK3_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    EARTH_tasmax = EARTH_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    EARTH3_tasmax = EARTH3_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    GFDL_tasmax = GFDL_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    HadGEM2_tasmax = HadGEM2_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    IPSLG_tasmax = IPSLG_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    MIROCG_tasmax = MIROCG_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    MPI_tasmax = MPI_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    NorESM1_tasmax = NorESM1_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    
    CRU_tasmax = CRU_tasmax.regrid(CanESM2_tasmax, iris.analysis.Linear())
    
    #we are only interested in the latitude and longitude relevant to Malawi 
    Malawi = iris.Constraint(longitude=lambda v: 32.0 <= v <= 36., latitude=lambda v: -17. <= v <= -8.)    
    
    CCCma_pr = CCCma_pr.extract(Malawi)
    CLMcom_pr = CLMcom_pr.extract(Malawi)
    DMI_pr  =DMI_pr.extract(Malawi)
    KNMI_pr =KNMI_pr.extract(Malawi)
    MPIE_pr = MPIE_pr.extract(Malawi)
    SMHI_pr = SMHI_pr.extract(Malawi)
        
    CCCmaCanRCM_pr = CCCmaCanRCM_pr.extract(Malawi)
    CCCmaSMHI_pr = CCCmaSMHI_pr.extract(Malawi)
    CNRM_pr = CNRM_pr.extract(Malawi)
    CNRMSMHI_pr = CNRMSMHI_pr.extract(Malawi)
    CSIRO_pr = CSIRO_pr.extract(Malawi)
    ICHECDMI_pr = ICHECDMI_pr.extract(Malawi)
    ICHECCCLM_pr = ICHECCCLM_pr.extract(Malawi)
    ICHECKNMI_pr = ICHECKNMI_pr.extract(Malawi)
    ICHECMPI_pr = ICHECMPI_pr.extract(Malawi)
    ICHECSMHI_pr = ICHECSMHI_pr.extract(Malawi)
    IPSL_pr = IPSL_pr.extract(Malawi)
    MIROC_pr = MIROC_pr.extract(Malawi)
    MOHCCCLM_pr = MOHCCCLM_pr.extract(Malawi)
    MOHCKNMI_pr = MOHCKNMI_pr.extract(Malawi)
    MOHCSMHI_pr = MOHCSMHI_pr.extract(Malawi)
    MPICCLM_pr = MPICCLM_pr.extract(Malawi)
    MPIREMO_pr = MPIREMO_pr.extract(Malawi)
    MPISMHI_pr = MPISMHI_pr.extract(Malawi)
    NCCDMI_pr = NCCDMI_pr.extract(Malawi)
    NCCSMHI_pr = NCCSMHI_pr.extract(Malawi)
    NOAA_pr = NOAA_pr.extract(Malawi)
    
    CanESM2_pr = CanESM2_pr.extract(Malawi)
    CNRMG_pr = CNRMG_pr.extract(Malawi)
    MK3_pr = MK3_pr.extract(Malawi)
    EARTH_pr = EARTH_pr.extract(Malawi)
    EARTH3_pr = EARTH3_pr.extract(Malawi)
    GFDL_pr = GFDL_pr.extract(Malawi)
    HadGEM2_pr = HadGEM2_pr.extract(Malawi)
    IPSLG_pr = IPSLG_pr.extract(Malawi)
    MIROCG_pr = MIROCG_pr.extract(Malawi)
    MPI_pr = MPI_pr.extract(Malawi)
    NorESM1_pr = NorESM1_pr.extract(Malawi)
    
    CRU_pr = CRU_pr.extract(Malawi)
    UDel_pr = UDel_pr.extract(Malawi)
    GPCC_pr = GPCC_pr.extract(Malawi)
    
    CCCma_tas = CCCma_tas.extract(Malawi)
    CLMcom_tas = CLMcom_tas.extract(Malawi)
    DMI_tas  =DMI_tas.extract(Malawi)
    KNMI_tas =KNMI_tas.extract(Malawi)
    MPIE_tas = MPIE_tas.extract(Malawi)
    SMHI_tas = SMHI_tas.extract(Malawi)
        
    CCCmaCanRCM_tas = CCCmaCanRCM_tas.extract(Malawi)
    CCCmaSMHI_tas = CCCmaSMHI_tas.extract(Malawi)
    CNRM_tas = CNRM_tas.extract(Malawi)
    CNRMSMHI_tas = CNRMSMHI_tas.extract(Malawi)
    CSIRO_tas = CSIRO_tas.extract(Malawi)
    ICHECDMI_tas = ICHECDMI_tas.extract(Malawi)
    ICHECCCLM_tas = ICHECCCLM_tas.extract(Malawi)
    ICHECKNMI_tas = ICHECKNMI_tas.extract(Malawi)
    ICHECMPI_tas = ICHECMPI_tas.extract(Malawi)
    ICHECSMHI_tas = ICHECSMHI_tas.extract(Malawi)
    IPSL_tas = IPSL_tas.extract(Malawi)
    MIROC_tas = MIROC_tas.extract(Malawi)
    MOHCCCLM_tas = MOHCCCLM_tas.extract(Malawi)
    MOHCKNMI_tas = MOHCKNMI_tas.extract(Malawi)
    MOHCSMHI_tas = MOHCSMHI_tas.extract(Malawi)
    MPICCLM_tas = MPICCLM_tas.extract(Malawi)
    MPIREMO_tas = MPIREMO_tas.extract(Malawi)
    MPISMHI_tas = MPISMHI_tas.extract(Malawi)
    NCCDMI_tas = NCCDMI_tas.extract(Malawi)
    NCCSMHI_tas = NCCSMHI_tas.extract(Malawi)
    NOAA_tas = NOAA_tas.extract(Malawi)
    
    CanESM2_tas = CanESM2_tas.extract(Malawi)
    CNRMG_tas = CNRMG_tas.extract(Malawi)
    MK3_tas = MK3_tas.extract(Malawi)
    EARTH_tas = EARTH_tas.extract(Malawi)
    EARTH3_tas = EARTH3_tas.extract(Malawi)
    GFDL_tas = GFDL_tas.extract(Malawi)
    HadGEM2_tas = HadGEM2_tas.extract(Malawi)
    IPSLG_tas = IPSLG_tas.extract(Malawi)
    MIROCG_tas = MIROCG_tas.extract(Malawi)
    MPI_tas = MPI_tas.extract(Malawi)
    NorESM1_tas = NorESM1_tas.extract(Malawi)
    
    CRU_tas = CRU_tas.extract(Malawi)
    UDel_tas = UDel_tas.extract(Malawi)
    
    CCCma_tasmin = CCCma_tasmin.extract(Malawi)
    CLMcom_tasmin = CLMcom_tasmin.extract(Malawi)
    DMI_tasmin  =DMI_tasmin.extract(Malawi)
    KNMI_tasmin =KNMI_tasmin.extract(Malawi)
    MPIE_tasmin = MPIE_tasmin.extract(Malawi)
    SMHI_tasmin = SMHI_tasmin.extract(Malawi)
        
    CCCmaCanRCM_tasmin = CCCmaCanRCM_tasmin.extract(Malawi)
    CCCmaSMHI_tasmin = CCCmaSMHI_tasmin.extract(Malawi)
    CNRM_tasmin = CNRM_tasmin.extract(Malawi)
    CNRMSMHI_tasmin = CNRMSMHI_tasmin.extract(Malawi)
    CSIRO_tasmin = CSIRO_tasmin.extract(Malawi)
    ICHECDMI_tasmin = ICHECDMI_tasmin.extract(Malawi)
    ICHECCCLM_tasmin = ICHECCCLM_tasmin.extract(Malawi)
    ICHECKNMI_tasmin = ICHECKNMI_tasmin.extract(Malawi)
    ICHECMPI_tasmin = ICHECMPI_tasmin.extract(Malawi)
    ICHECSMHI_tasmin = ICHECSMHI_tasmin.extract(Malawi)
    IPSL_tasmin = IPSL_tasmin.extract(Malawi)
    MIROC_tasmin = MIROC_tasmin.extract(Malawi)
    MOHCCCLM_tasmin = MOHCCCLM_tasmin.extract(Malawi)
    MOHCKNMI_tasmin = MOHCKNMI_tasmin.extract(Malawi)
    MOHCSMHI_tasmin = MOHCSMHI_tasmin.extract(Malawi)
    MPICCLM_tasmin = MPICCLM_tasmin.extract(Malawi)
    MPIREMO_tasmin = MPIREMO_tasmin.extract(Malawi)
    MPISMHI_tasmin = MPISMHI_tasmin.extract(Malawi)
    NCCDMI_tasmin = NCCDMI_tasmin.extract(Malawi)
    NCCSMHI_tasmin = NCCSMHI_tasmin.extract(Malawi)
    NOAA_tasmin = NOAA_tasmin.extract(Malawi)
    
    CanESM2_tasmin = CanESM2_tasmin.extract(Malawi)
    CNRMG_tasmin = CNRMG_tasmin.extract(Malawi)
    MK3_tasmin = MK3_tasmin.extract(Malawi)
    EARTH_tasmin = EARTH_tasmin.extract(Malawi)
    EARTH3_tasmin = EARTH3_tasmin.extract(Malawi)
    GFDL_tasmin = GFDL_tasmin.extract(Malawi)
    HadGEM2_tasmin = HadGEM2_tasmin.extract(Malawi)
    IPSLG_tasmin = IPSLG_tasmin.extract(Malawi)
    MIROCG_tasmin = MIROCG_tasmin.extract(Malawi)
    MPI_tasmin = MPI_tasmin.extract(Malawi)
    NorESM1_tasmin = NorESM1_tasmin.extract(Malawi)
    
    CRU_tasmin = CRU_tasmin.extract(Malawi)
  
    CCCma_tasmax = CCCma_tasmax.extract(Malawi)
    CLMcom_tasmax = CLMcom_tasmax.extract(Malawi)
    DMI_tasmax  =DMI_tasmax.extract(Malawi)
    KNMI_tasmax =KNMI_tasmax.extract(Malawi)
    MPIE_tasmax = MPIE_tasmax.extract(Malawi)
    SMHI_tasmax = SMHI_tasmax.extract(Malawi)
        
    CCCmaCanRCM_tasmax = CCCmaCanRCM_tasmax.extract(Malawi)
    CCCmaSMHI_tasmax = CCCmaSMHI_tasmax.extract(Malawi)
    CNRM_tasmax = CNRM_tasmax.extract(Malawi)
    CNRMSMHI_tasmax = CNRMSMHI_tasmax.extract(Malawi)
    CSIRO_tasmax = CSIRO_tasmax.extract(Malawi)
    ICHECDMI_tasmax = ICHECDMI_tasmax.extract(Malawi)
    ICHECCCLM_tasmax = ICHECCCLM_tasmax.extract(Malawi)
    ICHECKNMI_tasmax = ICHECKNMI_tasmax.extract(Malawi)
    ICHECMPI_tasmax = ICHECMPI_tasmax.extract(Malawi)
    ICHECSMHI_tasmax = ICHECSMHI_tasmax.extract(Malawi)
    IPSL_tasmax = IPSL_tasmax.extract(Malawi)
    MIROC_tasmax = MIROC_tasmax.extract(Malawi)
    MOHCCCLM_tasmax = MOHCCCLM_tasmax.extract(Malawi)
    MOHCKNMI_tasmax = MOHCKNMI_tasmax.extract(Malawi)
    MOHCSMHI_tasmax = MOHCSMHI_tasmax.extract(Malawi)
    MPICCLM_tasmax = MPICCLM_tasmax.extract(Malawi)
    MPIREMO_tasmax = MPIREMO_tasmax.extract(Malawi)
    MPISMHI_tasmax = MPISMHI_tasmax.extract(Malawi)
    NCCDMI_tasmax = NCCDMI_tasmax.extract(Malawi)
    NCCSMHI_tasmax = NCCSMHI_tasmax.extract(Malawi)
    NOAA_tasmax = NOAA_tasmax.extract(Malawi)
    
    CanESM2_tasmax = CanESM2_tasmax.extract(Malawi)
    CNRMG_tasmax = CNRMG_tasmax.extract(Malawi)
    MK3_tasmax = MK3_tasmax.extract(Malawi)
    EARTH_tasmax = EARTH_tasmax.extract(Malawi)
    EARTH3_tasmax = EARTH3_tasmax.extract(Malawi)
    GFDL_tasmax = GFDL_tasmax.extract(Malawi)
    HadGEM2_tasmax = HadGEM2_tasmax.extract(Malawi)
    IPSLG_tasmax = IPSLG_tasmax.extract(Malawi)
    MIROCG_tasmax = MIROCG_tasmax.extract(Malawi)
    MPI_tasmax = MPI_tasmax.extract(Malawi)
    NorESM1_tasmax = NorESM1_tasmax.extract(Malawi)
    
    CRU_tasmax = CRU_tasmax.extract(Malawi)
    
    #add year to all files
    iriscc.add_year(CCCma_pr, 'time')
    iriscc.add_year(CLMcom_pr, 'time')
    iriscc.add_year(DMI_pr, 'time')
    iriscc.add_year(KNMI_pr, 'time')
    iriscc.add_year(MPIE_pr, 'time')
    iriscc.add_year(SMHI_pr, 'time')
        
    iriscc.add_year(CCCmaCanRCM_pr, 'time')
    iriscc.add_year(CCCmaSMHI_pr, 'time')
    iriscc.add_year(CNRM_pr, 'time')
    iriscc.add_year(CNRMSMHI_pr, 'time')
    iriscc.add_year(CSIRO_pr, 'time')
    iriscc.add_year(ICHECDMI_pr, 'time')
    iriscc.add_year(ICHECCCLM_pr, 'time')
    iriscc.add_year(ICHECKNMI_pr, 'time')
    iriscc.add_year(ICHECMPI_pr, 'time')
    iriscc.add_year(ICHECSMHI_pr, 'time')
    iriscc.add_year(IPSL_pr, 'time')
    iriscc.add_year(MIROC_pr, 'time')
    iriscc.add_year(MOHCCCLM_pr, 'time')
    iriscc.add_year(MOHCKNMI_pr, 'time')
    iriscc.add_year(MOHCSMHI_pr, 'time')
    iriscc.add_year(MPICCLM_pr, 'time')
    iriscc.add_year(MPIREMO_pr, 'time')
    iriscc.add_year(MPISMHI_pr, 'time')
    iriscc.add_year(NCCDMI_pr, 'time')
    iriscc.add_year(NCCSMHI_pr, 'time')
    iriscc.add_year(NOAA_pr, 'time')
    
    iriscc.add_year(CanESM2_pr, 'time')
    iriscc.add_year(CNRMG_pr, 'time')
    iriscc.add_year(MK3_pr, 'time')
    iriscc.add_year(EARTH_pr, 'time')
    iriscc.add_year(EARTH3_pr, 'time')
    iriscc.add_year(GFDL_pr, 'time')
    iriscc.add_year(HadGEM2_pr, 'time')
    iriscc.add_year(IPSLG_pr, 'time')
    iriscc.add_year(MIROCG_pr, 'time')
    iriscc.add_year(MPI_pr, 'time')
    iriscc.add_year(NorESM1_pr, 'time')
    
    iriscc.add_year(CRU_pr, 'time')
    iriscc.add_year(UDel_pr, 'time')
    iriscc.add_year(GPCC_pr, 'time')
    
    iriscc.add_year(CCCma_tas, 'time')
    iriscc.add_year(CLMcom_tas, 'time')
    iriscc.add_year(DMI_tas, 'time')
    iriscc.add_year(KNMI_tas, 'time')
    iriscc.add_year(MPIE_tas, 'time')
    iriscc.add_year(SMHI_tas, 'time')
        
    iriscc.add_year(CCCmaCanRCM_tas, 'time')
    iriscc.add_year(CCCmaSMHI_tas, 'time')
    iriscc.add_year(CNRM_tas, 'time')
    iriscc.add_year(CNRMSMHI_tas, 'time')
    iriscc.add_year(CSIRO_tas, 'time')
    iriscc.add_year(ICHECDMI_tas, 'time')
    iriscc.add_year(ICHECCCLM_tas, 'time')
    iriscc.add_year(ICHECKNMI_tas, 'time')
    iriscc.add_year(ICHECMPI_tas, 'time')
    iriscc.add_year(ICHECSMHI_tas, 'time')
    iriscc.add_year(IPSL_tas, 'time')
    iriscc.add_year(MIROC_tas, 'time')
    iriscc.add_year(MOHCCCLM_tas, 'time')
    iriscc.add_year(MOHCKNMI_tas, 'time')
    iriscc.add_year(MOHCSMHI_tas, 'time')
    iriscc.add_year(MPICCLM_tas, 'time')
    iriscc.add_year(MPIREMO_tas, 'time')
    iriscc.add_year(MPISMHI_tas, 'time')
    iriscc.add_year(NCCDMI_tas, 'time')
    iriscc.add_year(NCCSMHI_tas, 'time')
    iriscc.add_year(NOAA_tas, 'time')
    
    iriscc.add_year(CanESM2_tas, 'time')
    iriscc.add_year(CNRMG_tas, 'time')
    iriscc.add_year(MK3_tas, 'time')
    iriscc.add_year(EARTH_tas, 'time')
    iriscc.add_year(EARTH3_tas, 'time')
    iriscc.add_year(GFDL_tas, 'time')
    iriscc.add_year(HadGEM2_tas, 'time')
    iriscc.add_year(IPSLG_tas, 'time')
    iriscc.add_year(MIROCG_tas, 'time')
    iriscc.add_year(MPI_tas, 'time')
    iriscc.add_year(NorESM1_tas, 'time')
    
    iriscc.add_year(CRU_tas, 'time')
    iriscc.add_year(UDel_tas, 'time')
    
    iriscc.add_year(CCCma_tasmin, 'time')
    iriscc.add_year(CLMcom_tasmin, 'time')
    iriscc.add_year(DMI_tasmin, 'time')
    iriscc.add_year(KNMI_tasmin, 'time')
    iriscc.add_year(MPIE_tasmin, 'time')
    iriscc.add_year(SMHI_tasmin, 'time')
        
    iriscc.add_year(CCCmaCanRCM_tasmin, 'time')
    iriscc.add_year(CCCmaSMHI_tasmin, 'time')
    iriscc.add_year(CNRM_tasmin, 'time')
    iriscc.add_year(CNRMSMHI_tasmin, 'time')
    iriscc.add_year(CSIRO_tasmin, 'time')
    iriscc.add_year(ICHECDMI_tasmin, 'time')
    iriscc.add_year(ICHECCCLM_tasmin, 'time')
    iriscc.add_year(ICHECKNMI_tasmin, 'time')
    iriscc.add_year(ICHECMPI_tasmin, 'time')
    iriscc.add_year(ICHECSMHI_tasmin, 'time')
    iriscc.add_year(IPSL_tasmin, 'time')
    iriscc.add_year(MIROC_tasmin, 'time')
    iriscc.add_year(MOHCCCLM_tasmin, 'time')
    iriscc.add_year(MOHCKNMI_tasmin, 'time')
    iriscc.add_year(MOHCSMHI_tasmin, 'time')
    iriscc.add_year(MPICCLM_tasmin, 'time')
    iriscc.add_year(MPIREMO_tasmin, 'time')
    iriscc.add_year(MPISMHI_tasmin, 'time')
    iriscc.add_year(NCCDMI_tasmin, 'time')
    iriscc.add_year(NCCSMHI_tasmin, 'time')
    iriscc.add_year(NOAA_tasmin, 'time')
    
    iriscc.add_year(CanESM2_tasmin, 'time')
    iriscc.add_year(CNRMG_tasmin, 'time')
    iriscc.add_year(MK3_tasmin, 'time')
    iriscc.add_year(EARTH_tasmin, 'time')
    iriscc.add_year(EARTH3_tasmin, 'time')
    iriscc.add_year(GFDL_tasmin, 'time')
    iriscc.add_year(HadGEM2_tasmin, 'time')
    iriscc.add_year(IPSLG_tasmin, 'time')
    iriscc.add_year(MIROCG_tasmin, 'time')
    iriscc.add_year(MPI_tasmin, 'time')
    iriscc.add_year(NorESM1_tasmin, 'time')
    
    iriscc.add_year(CRU_tasmin, 'time')
    
    iriscc.add_year(CCCma_tasmax, 'time')
    iriscc.add_year(CLMcom_tasmax, 'time')
    iriscc.add_year(DMI_tasmax, 'time')
    iriscc.add_year(KNMI_tasmax, 'time')
    iriscc.add_year(MPIE_tasmax, 'time')
    iriscc.add_year(SMHI_tasmax, 'time')
        
    iriscc.add_year(CCCmaCanRCM_tasmax, 'time')
    iriscc.add_year(CCCmaSMHI_tasmax, 'time')
    iriscc.add_year(CNRM_tasmax, 'time')
    iriscc.add_year(CNRMSMHI_tasmax, 'time')
    iriscc.add_year(CSIRO_tasmax, 'time')
    iriscc.add_year(ICHECDMI_tasmax, 'time')
    iriscc.add_year(ICHECCCLM_tasmax, 'time')
    iriscc.add_year(ICHECKNMI_tasmax, 'time')
    iriscc.add_year(ICHECMPI_tasmax, 'time')
    iriscc.add_year(ICHECSMHI_tasmax, 'time')
    iriscc.add_year(IPSL_tasmax, 'time')
    iriscc.add_year(MIROC_tasmax, 'time')
    iriscc.add_year(MOHCCCLM_tasmax, 'time')
    iriscc.add_year(MOHCKNMI_tasmax, 'time')
    iriscc.add_year(MOHCSMHI_tasmax, 'time')
    iriscc.add_year(MPICCLM_tasmax, 'time')
    iriscc.add_year(MPIREMO_tasmax, 'time')
    iriscc.add_year(MPISMHI_tasmax, 'time')
    iriscc.add_year(NCCDMI_tasmax, 'time')
    iriscc.add_year(NCCSMHI_tasmax, 'time')
    iriscc.add_year(NOAA_tasmax, 'time')
    
    iriscc.add_year(CanESM2_tasmax, 'time')
    iriscc.add_year(CNRMG_tasmax, 'time')
    iriscc.add_year(MK3_tasmax, 'time')
    iriscc.add_year(EARTH_tasmax, 'time')
    iriscc.add_year(EARTH3_tasmax, 'time')
    iriscc.add_year(GFDL_tasmax, 'time')
    iriscc.add_year(HadGEM2_tasmax, 'time')
    iriscc.add_year(IPSLG_tasmax, 'time')
    iriscc.add_year(MIROCG_tasmax, 'time')
    iriscc.add_year(MPI_tasmax, 'time')
    iriscc.add_year(NorESM1_tasmax, 'time')
    
    iriscc.add_year(CRU_tasmax, 'time')

    #add month to all files
    iriscc.add_month_number(CCCma_pr, 'time')
    iriscc.add_month_number(CLMcom_pr, 'time')
    iriscc.add_month_number(DMI_pr, 'time')
    iriscc.add_month_number(KNMI_pr, 'time')
    iriscc.add_month_number(MPIE_pr, 'time')
    iriscc.add_month_number(SMHI_pr, 'time')
        
    iriscc.add_month_number(CCCmaCanRCM_pr, 'time')
    iriscc.add_month_number(CCCmaSMHI_pr, 'time')
    iriscc.add_month_number(CNRM_pr, 'time')
    iriscc.add_month_number(CNRMSMHI_pr, 'time')
    iriscc.add_month_number(CSIRO_pr, 'time')
    iriscc.add_month_number(ICHECDMI_pr, 'time')
    iriscc.add_month_number(ICHECCCLM_pr, 'time')
    iriscc.add_month_number(ICHECKNMI_pr, 'time')
    iriscc.add_month_number(ICHECMPI_pr, 'time')
    iriscc.add_month_number(ICHECSMHI_pr, 'time')
    iriscc.add_month_number(IPSL_pr, 'time')
    iriscc.add_month_number(MIROC_pr, 'time')
    iriscc.add_month_number(MOHCCCLM_pr, 'time')
    iriscc.add_month_number(MOHCKNMI_pr, 'time')
    iriscc.add_month_number(MOHCSMHI_pr, 'time')
    iriscc.add_month_number(MPICCLM_pr, 'time')
    iriscc.add_month_number(MPIREMO_pr, 'time')
    iriscc.add_month_number(MPISMHI_pr, 'time')
    iriscc.add_month_number(NCCDMI_pr, 'time')
    iriscc.add_month_number(NCCSMHI_pr, 'time')
    iriscc.add_month_number(NOAA_pr, 'time')
    
    iriscc.add_month_number(CanESM2_pr, 'time')
    iriscc.add_month_number(CNRMG_pr, 'time')
    iriscc.add_month_number(MK3_pr, 'time')
    iriscc.add_month_number(EARTH_pr, 'time')
    iriscc.add_month_number(EARTH3_pr, 'time')
    iriscc.add_month_number(GFDL_pr, 'time')
    iriscc.add_month_number(HadGEM2_pr, 'time')
    iriscc.add_month_number(IPSLG_pr, 'time')
    iriscc.add_month_number(MIROCG_pr, 'time')
    iriscc.add_month_number(MPI_pr, 'time')
    iriscc.add_month_number(NorESM1_pr, 'time')
    
    iriscc.add_month_number(CRU_pr, 'time')
    iriscc.add_month_number(UDel_pr, 'time')
    iriscc.add_month_number(GPCC_pr, 'time')
    
    iriscc.add_month_number(CCCma_tas, 'time')
    iriscc.add_month_number(CLMcom_tas, 'time')
    iriscc.add_month_number(DMI_tas, 'time')
    iriscc.add_month_number(KNMI_tas, 'time')
    iriscc.add_month_number(MPIE_tas, 'time')
    iriscc.add_month_number(SMHI_tas, 'time')
        
    iriscc.add_month_number(CCCmaCanRCM_tas, 'time')
    iriscc.add_month_number(CCCmaSMHI_tas, 'time')
    iriscc.add_month_number(CNRM_tas, 'time')
    iriscc.add_month_number(CNRMSMHI_tas, 'time')
    iriscc.add_month_number(CSIRO_tas, 'time')
    iriscc.add_month_number(ICHECDMI_tas, 'time')
    iriscc.add_month_number(ICHECCCLM_tas, 'time')
    iriscc.add_month_number(ICHECKNMI_tas, 'time')
    iriscc.add_month_number(ICHECMPI_tas, 'time')
    iriscc.add_month_number(ICHECSMHI_tas, 'time')
    iriscc.add_month_number(IPSL_tas, 'time')
    iriscc.add_month_number(MIROC_tas, 'time')
    iriscc.add_month_number(MOHCCCLM_tas, 'time')
    iriscc.add_month_number(MOHCKNMI_tas, 'time')
    iriscc.add_month_number(MOHCSMHI_tas, 'time')
    iriscc.add_month_number(MPICCLM_tas, 'time')
    iriscc.add_month_number(MPIREMO_tas, 'time')
    iriscc.add_month_number(MPISMHI_tas, 'time')
    iriscc.add_month_number(NCCDMI_tas, 'time')
    iriscc.add_month_number(NCCSMHI_tas, 'time')
    iriscc.add_month_number(NOAA_tas, 'time')
    
    iriscc.add_month_number(CanESM2_tas, 'time')
    iriscc.add_month_number(CNRMG_tas, 'time')
    iriscc.add_month_number(MK3_tas, 'time')
    iriscc.add_month_number(EARTH_tas, 'time')
    iriscc.add_month_number(EARTH3_tas, 'time')
    iriscc.add_month_number(GFDL_tas, 'time')
    iriscc.add_month_number(HadGEM2_tas, 'time')
    iriscc.add_month_number(IPSLG_tas, 'time')
    iriscc.add_month_number(MIROCG_tas, 'time')
    iriscc.add_month_number(MPI_tas, 'time')
    iriscc.add_month_number(NorESM1_tas, 'time')
    
    iriscc.add_month_number(CRU_tas, 'time')
    iriscc.add_month_number(UDel_tas, 'time')
    
    iriscc.add_month_number(CCCma_tasmin, 'time')
    iriscc.add_month_number(CLMcom_tasmin, 'time')
    iriscc.add_month_number(DMI_tasmin, 'time')
    iriscc.add_month_number(KNMI_tasmin, 'time')
    iriscc.add_month_number(MPIE_tasmin, 'time')
    iriscc.add_month_number(SMHI_tasmin, 'time')
        
    iriscc.add_month_number(CCCmaCanRCM_tasmin, 'time')
    iriscc.add_month_number(CCCmaSMHI_tasmin, 'time')
    iriscc.add_month_number(CNRM_tasmin, 'time')
    iriscc.add_month_number(CNRMSMHI_tasmin, 'time')
    iriscc.add_month_number(CSIRO_tasmin, 'time')
    iriscc.add_month_number(ICHECDMI_tasmin, 'time')
    iriscc.add_month_number(ICHECCCLM_tasmin, 'time')
    iriscc.add_month_number(ICHECKNMI_tasmin, 'time')
    iriscc.add_month_number(ICHECMPI_tasmin, 'time')
    iriscc.add_month_number(ICHECSMHI_tasmin, 'time')
    iriscc.add_month_number(IPSL_tasmin, 'time')
    iriscc.add_month_number(MIROC_tasmin, 'time')
    iriscc.add_month_number(MOHCCCLM_tasmin, 'time')
    iriscc.add_month_number(MOHCKNMI_tasmin, 'time')
    iriscc.add_month_number(MOHCSMHI_tasmin, 'time')
    iriscc.add_month_number(MPICCLM_tasmin, 'time')
    iriscc.add_month_number(MPIREMO_tasmin, 'time')
    iriscc.add_month_number(MPISMHI_tasmin, 'time')
    iriscc.add_month_number(NCCDMI_tasmin, 'time')
    iriscc.add_month_number(NCCSMHI_tasmin, 'time')
    iriscc.add_month_number(NOAA_tasmin, 'time')
    
    iriscc.add_month_number(CanESM2_tasmin, 'time')
    iriscc.add_month_number(CNRMG_tasmin, 'time')
    iriscc.add_month_number(MK3_tasmin, 'time')
    iriscc.add_month_number(EARTH_tasmin, 'time')
    iriscc.add_month_number(EARTH3_tasmin, 'time')
    iriscc.add_month_number(GFDL_tasmin, 'time')
    iriscc.add_month_number(HadGEM2_tasmin, 'time')
    iriscc.add_month_number(IPSLG_tasmin, 'time')
    iriscc.add_month_number(MIROCG_tasmin, 'time')
    iriscc.add_month_number(MPI_tasmin, 'time')
    iriscc.add_month_number(NorESM1_tasmin, 'time')
    
    iriscc.add_month_number(CRU_tasmin, 'time')
    
    iriscc.add_month_number(CCCma_tasmax, 'time')
    iriscc.add_month_number(CLMcom_tasmax, 'time')
    iriscc.add_month_number(DMI_tasmax, 'time')
    iriscc.add_month_number(KNMI_tasmax, 'time')
    iriscc.add_month_number(MPIE_tasmax, 'time')
    iriscc.add_month_number(SMHI_tasmax, 'time')
        
    iriscc.add_month_number(CCCmaCanRCM_tasmax, 'time')
    iriscc.add_month_number(CCCmaSMHI_tasmax, 'time')
    iriscc.add_month_number(CNRM_tasmax, 'time')
    iriscc.add_month_number(CNRMSMHI_tasmax, 'time')
    iriscc.add_month_number(CSIRO_tasmax, 'time')
    iriscc.add_month_number(ICHECDMI_tasmax, 'time')
    iriscc.add_month_number(ICHECCCLM_tasmax, 'time')
    iriscc.add_month_number(ICHECKNMI_tasmax, 'time')
    iriscc.add_month_number(ICHECMPI_tasmax, 'time')
    iriscc.add_month_number(ICHECSMHI_tasmax, 'time')
    iriscc.add_month_number(IPSL_tasmax, 'time')
    iriscc.add_month_number(MIROC_tasmax, 'time')
    iriscc.add_month_number(MOHCCCLM_tasmax, 'time')
    iriscc.add_month_number(MOHCKNMI_tasmax, 'time')
    iriscc.add_month_number(MOHCSMHI_tasmax, 'time')
    iriscc.add_month_number(MPICCLM_tasmax, 'time')
    iriscc.add_month_number(MPIREMO_tasmax, 'time')
    iriscc.add_month_number(MPISMHI_tasmax, 'time')
    iriscc.add_month_number(NCCDMI_tasmax, 'time')
    iriscc.add_month_number(NCCSMHI_tasmax, 'time')
    iriscc.add_month_number(NOAA_tasmax, 'time')
    
    iriscc.add_month_number(CanESM2_tasmax, 'time')
    iriscc.add_month_number(CNRMG_tasmax, 'time')
    iriscc.add_month_number(MK3_tasmax, 'time')
    iriscc.add_month_number(EARTH_tasmax, 'time')
    iriscc.add_month_number(EARTH3_tasmax, 'time')
    iriscc.add_month_number(GFDL_tasmax, 'time')
    iriscc.add_month_number(HadGEM2_tasmax, 'time')
    iriscc.add_month_number(IPSLG_tasmax, 'time')
    iriscc.add_month_number(MIROCG_tasmax, 'time')
    iriscc.add_month_number(MPI_tasmax, 'time')
    iriscc.add_month_number(NorESM1_tasmax, 'time')
    
    iriscc.add_month_number(CRU_tasmax, 'time')

    #time constraint to make all series the same, for ERAINT this is 1990-2008 and for RCMs and GCMs this is 1961-2005
    iris.FUTURE.cell_datetime_objects = True
    t_constraint_ERAINT = iris.Constraint(time=lambda cell: 1990 <= cell.point.year <= 2008)
   
    CCCma_pr = CCCma_pr.extract(t_constraint_ERAINT)
    CLMcom_pr = CLMcom_pr.extract(t_constraint_ERAINT)
    DMI_pr = DMI_pr.extract(t_constraint_ERAINT)
    KNMI_pr = KNMI_pr.extract(t_constraint_ERAINT)
    MPIE_pr = MPIE_pr.extract(t_constraint_ERAINT)
    SMHI_pr = SMHI_pr.extract(t_constraint_ERAINT)
    
    CCCma_tas = CCCma_tas.extract(t_constraint_ERAINT)
    CLMcom_tas = CLMcom_tas.extract(t_constraint_ERAINT)
    DMI_tas = DMI_tas.extract(t_constraint_ERAINT)
    KNMI_tas = KNMI_tas.extract(t_constraint_ERAINT)
    MPIE_tas = MPIE_tas.extract(t_constraint_ERAINT)
    SMHI_tas = SMHI_tas.extract(t_constraint_ERAINT)
    
    CCCma_tasmin = CCCma_tasmin.extract(t_constraint_ERAINT)
    CLMcom_tasmin = CLMcom_tasmin.extract(t_constraint_ERAINT)
    DMI_tasmin = DMI_tasmin.extract(t_constraint_ERAINT)
    KNMI_tasmin = KNMI_tasmin.extract(t_constraint_ERAINT)
    MPIE_tasmin = MPIE_tasmin.extract(t_constraint_ERAINT)
    SMHI_tasmin = SMHI_tasmin.extract(t_constraint_ERAINT)
    
    CCCma_tasmax = CCCma_tasmax.extract(t_constraint_ERAINT)
    CLMcom_tasmax = CLMcom_tasmax.extract(t_constraint_ERAINT)
    DMI_tasmax = DMI_tasmax.extract(t_constraint_ERAINT)
    KNMI_tasmax = KNMI_tasmax.extract(t_constraint_ERAINT)
    MPIE_tasmax = MPIE_tasmax.extract(t_constraint_ERAINT)
    SMHI_tasmax = SMHI_tasmax.extract(t_constraint_ERAINT)
    
    t_constraint = iris.Constraint(time=lambda cell: 1961 <= cell.point.year <= 2005)
    
    CCCmaCanRCM_pr = CCCmaCanRCM_pr.extract(t_constraint)
    CCCmaSMHI_pr = CCCmaSMHI_pr.extract(t_constraint)
    CNRM_pr = CNRM_pr.extract(t_constraint)
    CNRMSMHI_pr = CNRMSMHI_pr.extract(t_constraint)
    CSIRO_pr = CSIRO_pr.extract(t_constraint)
    ICHECDMI_pr = ICHECDMI_pr.extract(t_constraint)
    ICHECCCLM_pr = ICHECCCLM_pr.extract(t_constraint)
    ICHECKNMI_pr = ICHECKNMI_pr.extract(t_constraint)
    ICHECMPI_pr = ICHECMPI_pr.extract(t_constraint)
    ICHECSMHI_pr = ICHECSMHI_pr.extract(t_constraint)
    IPSL_pr = IPSL_pr.extract(t_constraint)
    MIROC_pr = MIROC_pr.extract(t_constraint)
    MOHCCCLM_pr = MOHCCCLM_pr.extract(t_constraint)
    MOHCKNMI_pr = MOHCKNMI_pr.extract(t_constraint)
    MOHCSMHI_pr = MOHCSMHI_pr.extract(t_constraint)
    MPICCLM_pr = MPICCLM_pr.extract(t_constraint)
    MPIREMO_pr = MPIREMO_pr.extract(t_constraint)
    MPISMHI_pr = MPISMHI_pr.extract(t_constraint)
    NCCDMI_pr = NCCDMI_pr.extract(t_constraint)
    NCCSMHI_pr = NCCSMHI_pr.extract(t_constraint)
    NOAA_pr = NOAA_pr.extract(t_constraint) 
    
    CanESM2_pr =  CanESM2_pr.extract(t_constraint)
    CNRMG_pr = CNRMG_pr.extract(t_constraint)
    MK3_pr = MK3_pr.extract(t_constraint)
    EARTH_pr = EARTH_pr.extract(t_constraint)
    EARTH3_pr = EARTH3_pr.extract(t_constraint)
    GFDL_pr = GFDL_pr.extract(t_constraint)
    HadGEM2_pr = HadGEM2_pr.extract(t_constraint)
    IPSLG_pr = IPSLG_pr.extract(t_constraint)
    MIROCG_pr = MIROCG_pr.extract(t_constraint)
    MPI_pr = MPI_pr.extract(t_constraint)
    NorESM1_pr = NorESM1_pr.extract(t_constraint)
    
    CRU_pr = CRU_pr.extract(t_constraint)
    UDel_pr = UDel_pr.extract(t_constraint)
    GPCC_pr = GPCC_pr.extract(t_constraint)
    
    CCCmaCanRCM_tas = CCCmaCanRCM_tas.extract(t_constraint)
    CCCmaSMHI_tas = CCCmaSMHI_tas.extract(t_constraint)
    CNRM_tas = CNRM_tas.extract(t_constraint)
    CNRMSMHI_tas = CNRMSMHI_tas.extract(t_constraint)
    CSIRO_tas = CSIRO_tas.extract(t_constraint)
    ICHECDMI_tas = ICHECDMI_tas.extract(t_constraint)
    ICHECCCLM_tas = ICHECCCLM_tas.extract(t_constraint)
    ICHECKNMI_tas = ICHECKNMI_tas.extract(t_constraint)
    ICHECMPI_tas = ICHECMPI_tas.extract(t_constraint)
    ICHECSMHI_tas = ICHECSMHI_tas.extract(t_constraint)
    IPSL_tas = IPSL_tas.extract(t_constraint)
    MIROC_tas = MIROC_tas.extract(t_constraint)
    MOHCCCLM_tas = MOHCCCLM_tas.extract(t_constraint)
    MOHCKNMI_tas = MOHCKNMI_tas.extract(t_constraint)
    MOHCSMHI_tas = MOHCSMHI_tas.extract(t_constraint)
    MPICCLM_tas = MPICCLM_tas.extract(t_constraint)
    MPIREMO_tas = MPIREMO_tas.extract(t_constraint)
    MPISMHI_tas = MPISMHI_tas.extract(t_constraint)
    NCCDMI_tas = NCCDMI_tas.extract(t_constraint)
    NCCSMHI_tas = NCCSMHI_tas.extract(t_constraint)
    NOAA_tas = NOAA_tas.extract(t_constraint) 
    
    CanESM2_tas =  CanESM2_tas.extract(t_constraint)
    CNRMG_tas = CNRMG_tas.extract(t_constraint)
    MK3_tas = MK3_tas.extract(t_constraint)
    EARTH_tas = EARTH_tas.extract(t_constraint)
    EARTH3_tas = EARTH3_tas.extract(t_constraint)
    GFDL_tas = GFDL_tas.extract(t_constraint)
    HadGEM2_tas = HadGEM2_tas.extract(t_constraint)
    IPSLG_tas = IPSLG_tas.extract(t_constraint)
    MIROCG_tas = MIROCG_tas.extract(t_constraint)
    MPI_tas = MPI_tas.extract(t_constraint)
    NorESM1_tas = NorESM1_tas.extract(t_constraint)
    
    CRU_tas = CRU_tas.extract(t_constraint)
    UDel_tas = UDel_tas.extract(t_constraint)
    
    CCCmaCanRCM_tasmin = CCCmaCanRCM_tasmin.extract(t_constraint)
    CCCmaSMHI_tasmin = CCCmaSMHI_tasmin.extract(t_constraint)
    CNRM_tasmin = CNRM_tasmin.extract(t_constraint)
    CNRMSMHI_tasmin = CNRMSMHI_tasmin.extract(t_constraint)
    CSIRO_tasmin = CSIRO_tasmin.extract(t_constraint)
    ICHECDMI_tasmin = ICHECDMI_tasmin.extract(t_constraint)
    ICHECCCLM_tasmin = ICHECCCLM_tasmin.extract(t_constraint)
    ICHECKNMI_tasmin = ICHECKNMI_tasmin.extract(t_constraint)
    ICHECMPI_tasmin = ICHECMPI_tasmin.extract(t_constraint)
    ICHECSMHI_tasmin = ICHECSMHI_tasmin.extract(t_constraint)
    IPSL_tasmin = IPSL_tasmin.extract(t_constraint)
    MIROC_tasmin = MIROC_tasmin.extract(t_constraint)
    MOHCCCLM_tasmin = MOHCCCLM_tasmin.extract(t_constraint)
    MOHCKNMI_tasmin = MOHCKNMI_tasmin.extract(t_constraint)
    MOHCSMHI_tasmin = MOHCSMHI_tasmin.extract(t_constraint)
    MPICCLM_tasmin = MPICCLM_tasmin.extract(t_constraint)
    MPIREMO_tasmin = MPIREMO_tasmin.extract(t_constraint)
    MPISMHI_tasmin = MPISMHI_tasmin.extract(t_constraint)
    NCCDMI_tasmin = NCCDMI_tasmin.extract(t_constraint)
    NCCSMHI_tasmin = NCCSMHI_tasmin.extract(t_constraint)
    NOAA_tasmin = NOAA_tasmin.extract(t_constraint) 
    
    CanESM2_tasmin =  CanESM2_tasmin.extract(t_constraint)
    CNRMG_tasmin = CNRMG_tasmin.extract(t_constraint)
    MK3_tasmin = MK3_tasmin.extract(t_constraint)
    EARTH_tasmin = EARTH_tasmin.extract(t_constraint)
    EARTH3_tasmin = EARTH3_tasmin.extract(t_constraint)
    GFDL_tasmin = GFDL_tasmin.extract(t_constraint)
    HadGEM2_tasmin = HadGEM2_tasmin.extract(t_constraint)
    IPSLG_tasmin = IPSLG_tasmin.extract(t_constraint)
    MIROCG_tasmin = MIROCG_tasmin.extract(t_constraint)
    MPI_tasmin = MPI_tasmin.extract(t_constraint)
    NorESM1_tasmin = NorESM1_tasmin.extract(t_constraint)
    
    CRU_tasmin = CRU_tasmin.extract(t_constraint)
    
    CCCmaCanRCM_tasmax = CCCmaCanRCM_tasmax.extract(t_constraint)
    CCCmaSMHI_tasmax = CCCmaSMHI_tasmax.extract(t_constraint)
    CNRM_tasmax = CNRM_tasmax.extract(t_constraint)
    CNRMSMHI_tasmax = CNRMSMHI_tasmax.extract(t_constraint)
    CSIRO_tasmax = CSIRO_tasmax.extract(t_constraint)
    ICHECDMI_tasmax = ICHECDMI_tasmax.extract(t_constraint)
    ICHECCCLM_tasmax = ICHECCCLM_tasmax.extract(t_constraint)
    ICHECKNMI_tasmax = ICHECKNMI_tasmax.extract(t_constraint)
    ICHECMPI_tasmax = ICHECMPI_tasmax.extract(t_constraint)
    ICHECSMHI_tasmax = ICHECSMHI_tasmax.extract(t_constraint)
    IPSL_tasmax = IPSL_tasmax.extract(t_constraint)
    MIROC_tasmax = MIROC_tasmax.extract(t_constraint)
    MOHCCCLM_tasmax = MOHCCCLM_tasmax.extract(t_constraint)
    MOHCKNMI_tasmax = MOHCKNMI_tasmax.extract(t_constraint)
    MOHCSMHI_tasmax = MOHCSMHI_tasmax.extract(t_constraint)
    MPICCLM_tasmax = MPICCLM_tasmax.extract(t_constraint)
    MPIREMO_tasmax = MPIREMO_tasmax.extract(t_constraint)
    MPISMHI_tasmax = MPISMHI_tasmax.extract(t_constraint)
    NCCDMI_tasmax = NCCDMI_tasmax.extract(t_constraint)
    NCCSMHI_tasmax = NCCSMHI_tasmax.extract(t_constraint)
    NOAA_tasmax = NOAA_tasmax.extract(t_constraint) 
    
    CanESM2_tasmax =  CanESM2_tasmax.extract(t_constraint)
    CNRMG_tasmax = CNRMG_tasmax.extract(t_constraint)
    MK3_tasmax = MK3_tasmax.extract(t_constraint)
    EARTH_tasmax = EARTH_tasmax.extract(t_constraint)
    EARTH3_tasmax = EARTH3_tasmax.extract(t_constraint)
    GFDL_tasmax = GFDL_tasmax.extract(t_constraint)
    HadGEM2_tasmax = HadGEM2_tasmax.extract(t_constraint)
    IPSLG_tasmax = IPSLG_tasmax.extract(t_constraint)
    MIROCG_tasmax = MIROCG_tasmax.extract(t_constraint)
    MPI_tasmax = MPI_tasmax.extract(t_constraint)
    NorESM1_tasmax = NorESM1_tasmax.extract(t_constraint)
    
    CRU_tasmax = CRU_tasmax.extract(t_constraint)
    
    #Convert units to match, CORDEX data in precipitation_flux (kg m-2 s-1) but want all data in precipitation rate (mm month-1).
    #Since 1 kg of rain spread over 1 m of surface is 1mm in thickness, and there are 60*60*24*365.25=31557600 seconds in a year and 12 months in the year, the conversion is:
    Convert=31557600/12
    CCCma_pr = CCCma_pr*Convert
    CLMcom_pr = CLMcom_pr*Convert
    DMI_pr = DMI_pr*Convert
    KNMI_pr = KNMI_pr*Convert
    MPIE_pr = MPIE_pr*Convert
    SMHI_pr = SMHI_pr*Convert
    
    CCCmaCanRCM_pr = CCCmaCanRCM_pr*Convert
    CCCmaSMHI_pr = CCCmaSMHI_pr*Convert
    CNRM_pr = CNRM_pr*Convert
    CNRMSMHI_pr = CNRMSMHI_pr*Convert 
    CSIRO_pr = CSIRO_pr*Convert
    ICHECDMI_pr = ICHECDMI_pr*Convert
    ICHECCCLM_pr = ICHECCCLM_pr*Convert
    ICHECKNMI_pr = ICHECKNMI_pr*Convert
    ICHECMPI_pr = ICHECMPI_pr*Convert
    ICHECSMHI_pr = ICHECSMHI_pr*Convert
    IPSL_pr = IPSL_pr*Convert
    MIROC_pr = MIROC_pr*Convert
    MOHCCCLM_pr = MOHCCCLM_pr*Convert
    MOHCKNMI_pr = MOHCKNMI_pr*Convert
    MOHCSMHI_pr = MOHCSMHI_pr*Convert
    MPICCLM_pr = MPICCLM_pr*Convert
    MPIREMO_pr = MPIREMO_pr*Convert
    MPISMHI_pr = MPISMHI_pr*Convert
    NCCDMI_pr = NCCDMI_pr*Convert
    NCCSMHI_pr = NCCSMHI_pr*Convert
    NOAA_pr = NOAA_pr*Convert
    
    CanESM2_pr = CanESM2_pr*Convert
    CNRMG_pr = CNRMG_pr*Convert
    MK3_pr = MK3_pr*Convert
    EARTH_pr = EARTH_pr*Convert
    EARTH3_pr = EARTH3_pr*Convert
    GFDL_pr = GFDL_pr*Convert
    HadGEM2_pr = HadGEM2_pr*Convert
    IPSLG_pr = IPSLG_pr*Convert
    MIROCG_pr = MIROCG_pr*Convert
    MPI_pr = MPI_pr*Convert
    NorESM1_pr = NorESM1_pr*Convert 
    
    #Convert units to match, UDel_pr data in cm per month but want precipitation rate in mm per month.
    #Since there are 10mm in a cm, the conversion is:
    Convert=10
    UDel_pr = UDel_pr*Convert
    
    #Convert units to match, CORDEX data is in Kelvin but Observed data in Celsius, we would like to show all data in Celsius
    CCCma_tas.convert_units('Celsius')
    CLMcom_tas.convert_units('Celsius')
    DMI_tas.convert_units('Celsius')
    KNMI_tas.convert_units('Celsius')
    MPIE_tas.convert_units('Celsius')
    SMHI_tas.convert_units('Celsius')
    
    CCCmaCanRCM_tas.convert_units('Celsius')
    CCCmaSMHI_tas.convert_units('Celsius')
    CNRM_tas.convert_units('Celsius')
    CNRMSMHI_tas.convert_units('Celsius')
    CSIRO_tas.convert_units('Celsius')
    ICHECDMI_tas.convert_units('Celsius')
    ICHECCCLM_tas.convert_units('Celsius') 
    ICHECKNMI_tas.convert_units('Celsius')
    ICHECMPI_tas.convert_units('Celsius')
    ICHECSMHI_tas.convert_units('Celsius')
    IPSL_tas.convert_units('Celsius')
    MIROC_tas.convert_units('Celsius')
    MOHCCCLM_tas.convert_units('Celsius')
    MOHCKNMI_tas.convert_units('Celsius')
    MOHCSMHI_tas.convert_units('Celsius')
    MPICCLM_tas.convert_units('Celsius')
    MPIREMO_tas.convert_units('Celsius')
    MPISMHI_tas.convert_units('Celsius')
    NCCDMI_tas.convert_units('Celsius')
    NCCSMHI_tas.convert_units('Celsius')
    NOAA_tas.convert_units('Celsius')   
    
    CanESM2_tas.convert_units('Celsius')
    CNRMG_tas.convert_units('Celsius')
    MK3_tas.convert_units('Celsius')
    EARTH_tas.convert_units('Celsius')
    EARTH3_tas.units = Unit('Celsius') #this fixes EARTH3 which has no units defined
    EARTH3_tas=EARTH3_tas-273 #this converts the data manually from Kelvin to Celsius
    GFDL_tas.convert_units('Celsius')
    HadGEM2_tas.convert_units('Celsius')
    IPSLG_tas.convert_units('Celsius')
    MIROCG_tas.convert_units('Celsius')
    MPI_tas.convert_units('Celsius')
    NorESM1_tas.convert_units('Celsius') 
    
    CCCma_tasmin.convert_units('Celsius')
    CLMcom_tasmin.convert_units('Celsius')
    DMI_tasmin.convert_units('Celsius')
    KNMI_tasmin.convert_units('Celsius')
    MPIE_tasmin.convert_units('Celsius')
    SMHI_tasmin.convert_units('Celsius')
    
    CCCmaCanRCM_tasmin.convert_units('Celsius')
    CCCmaSMHI_tasmin.convert_units('Celsius')
    CNRM_tasmin.convert_units('Celsius')
    CNRMSMHI_tasmin.convert_units('Celsius')
    CSIRO_tasmin.convert_units('Celsius')
    ICHECDMI_tasmin.convert_units('Celsius')
    ICHECCCLM_tasmin.convert_units('Celsius') 
    ICHECKNMI_tasmin.convert_units('Celsius')
    ICHECMPI_tasmin.convert_units('Celsius')
    ICHECSMHI_tasmin.convert_units('Celsius')
    IPSL_tasmin.convert_units('Celsius')
    MIROC_tasmin.convert_units('Celsius')
    MOHCCCLM_tasmin.convert_units('Celsius')
    MOHCKNMI_tasmin.convert_units('Celsius')
    MOHCSMHI_tasmin.convert_units('Celsius')
    MPICCLM_tasmin.convert_units('Celsius')
    MPIREMO_tasmin.convert_units('Celsius')
    MPISMHI_tasmin.convert_units('Celsius')
    NCCDMI_tasmin.convert_units('Celsius')
    NCCSMHI_tasmin.convert_units('Celsius')
    NOAA_tasmin.convert_units('Celsius')   
    
    CanESM2_tasmin.convert_units('Celsius')
    CNRMG_tasmin.convert_units('Celsius')
    MK3_tasmin.convert_units('Celsius')
    EARTH_tasmin.convert_units('Celsius')
    EARTH3_tasmin.units = Unit('Celsius') #this fixes EARTH3 which has no units defined
    EARTH3_tasmin=EARTH3_tasmin-273 #this converts the data manually from Kelvin to Celsius
    GFDL_tasmin.convert_units('Celsius')
    HadGEM2_tasmin.convert_units('Celsius')
    IPSLG_tasmin.convert_units('Celsius')
    MIROCG_tasmin.convert_units('Celsius')
    MPI_tasmin.convert_units('Celsius')
    NorESM1_tasmin.convert_units('Celsius') 
    
    CCCma_tasmax.convert_units('Celsius')
    CLMcom_tasmax.convert_units('Celsius')
    DMI_tasmax.convert_units('Celsius')
    KNMI_tasmax.convert_units('Celsius')
    MPIE_tasmax.convert_units('Celsius')
    SMHI_tasmax.convert_units('Celsius')
    
    CCCmaCanRCM_tasmax.convert_units('Celsius')
    CCCmaSMHI_tasmax.convert_units('Celsius')
    CNRM_tasmax.convert_units('Celsius')
    CNRMSMHI_tasmax.convert_units('Celsius')
    CSIRO_tasmax.convert_units('Celsius')
    ICHECDMI_tasmax.convert_units('Celsius')
    ICHECCCLM_tasmax.convert_units('Celsius') 
    ICHECKNMI_tasmax.convert_units('Celsius')
    ICHECMPI_tasmax.convert_units('Celsius')
    ICHECSMHI_tasmax.convert_units('Celsius')
    IPSL_tasmax.convert_units('Celsius')
    MIROC_tasmax.convert_units('Celsius')
    MOHCCCLM_tasmax.convert_units('Celsius')
    MOHCKNMI_tasmax.convert_units('Celsius')
    MOHCSMHI_tasmax.convert_units('Celsius')
    MPICCLM_tasmax.convert_units('Celsius')
    MPIREMO_tasmax.convert_units('Celsius')
    MPISMHI_tasmax.convert_units('Celsius')
    NCCDMI_tasmax.convert_units('Celsius')
    NCCSMHI_tasmax.convert_units('Celsius')
    NOAA_tasmax.convert_units('Celsius')   
    
    CanESM2_tasmax.convert_units('Celsius')
    CNRMG_tasmax.convert_units('Celsius')
    MK3_tasmax.convert_units('Celsius')
    EARTH_tasmax.convert_units('Celsius')
    EARTH3_tasmax.units = Unit('Celsius') #this fixes EARTH3 which has no units defined
    EARTH3_tasmax=EARTH3_tasmax-273 #this converts the data manually from Kelvin to Celsius
    GFDL_tasmax.convert_units('Celsius')
    HadGEM2_tasmax.convert_units('Celsius')
    IPSLG_tasmax.convert_units('Celsius')
    MIROCG_tasmax.convert_units('Celsius')
    MPI_tasmax.convert_units('Celsius')
    NorESM1_tasmax.convert_units('Celsius') 
    
    #Returns an array of area weights, with the same dimensions as the cube
    CCCma_pr_grid_areas = iris.analysis.cartography.area_weights(CCCma_pr)
    CLMcom_pr_grid_areas = iris.analysis.cartography.area_weights(CLMcom_pr)
    DMI_pr_grid_areas = iris.analysis.cartography.area_weights(DMI_pr)
    KNMI_pr_grid_areas = iris.analysis.cartography.area_weights(KNMI_pr)
    MPIE_pr_grid_areas = iris.analysis.cartography.area_weights(MPIE_pr)
    SMHI_pr_grid_areas = iris.analysis.cartography.area_weights(SMHI_pr)
    
    CCCmaCanRCM_pr_grid_areas = iris.analysis.cartography.area_weights(CCCmaCanRCM_pr)
    CCCmaSMHI_pr_grid_areas = iris.analysis.cartography.area_weights(CCCmaSMHI_pr)
    CNRM_pr_grid_areas = iris.analysis.cartography.area_weights(CNRM_pr)
    CNRMSMHI_pr_grid_areas = iris.analysis.cartography.area_weights(CNRMSMHI_pr)
    CSIRO_pr_grid_areas = iris.analysis.cartography.area_weights(CSIRO_pr)
    ICHECDMI_pr_grid_areas = iris.analysis.cartography.area_weights(ICHECDMI_pr)
    ICHECCCLM_pr_grid_areas = iris.analysis.cartography.area_weights(ICHECCCLM_pr)
    ICHECKNMI_pr_grid_areas = iris.analysis.cartography.area_weights(ICHECKNMI_pr)
    ICHECMPI_pr_grid_areas = iris.analysis.cartography.area_weights(ICHECMPI_pr)
    ICHECSMHI_pr_grid_areas = iris.analysis.cartography.area_weights(ICHECSMHI_pr)
    IPSL_pr_grid_areas = iris.analysis.cartography.area_weights(IPSL_pr)
    MIROC_pr_grid_areas = iris.analysis.cartography.area_weights(MIROC_pr)
    MOHCCCLM_pr_grid_areas = iris.analysis.cartography.area_weights(MOHCCCLM_pr)
    MOHCKNMI_pr_grid_areas = iris.analysis.cartography.area_weights(MOHCKNMI_pr)
    MOHCSMHI_pr_grid_areas = iris.analysis.cartography.area_weights(MOHCSMHI_pr)
    MPICCLM_pr_grid_areas = iris.analysis.cartography.area_weights(MPICCLM_pr)
    MPIREMO_pr_grid_areas = iris.analysis.cartography.area_weights(MPIREMO_pr)
    MPISMHI_pr_grid_areas = iris.analysis.cartography.area_weights(MPISMHI_pr)
    NCCDMI_pr_grid_areas = iris.analysis.cartography.area_weights(NCCDMI_pr)
    NCCSMHI_pr_grid_areas = iris.analysis.cartography.area_weights(NCCSMHI_pr)
    NOAA_pr_grid_areas = iris.analysis.cartography.area_weights(NOAA_pr)
    
    CanESM2_pr_grid_areas = iris.analysis.cartography.area_weights(CanESM2_pr)
    CNRMG_pr_grid_areas = iris.analysis.cartography.area_weights(CNRMG_pr)
    MK3_pr_grid_areas = iris.analysis.cartography.area_weights(MK3_pr)
    EARTH_pr_grid_areas = iris.analysis.cartography.area_weights(EARTH_pr)
    EARTH3_pr_grid_areas = iris.analysis.cartography.area_weights(EARTH3_pr)
    GFDL_pr_grid_areas = iris.analysis.cartography.area_weights(GFDL_pr)
    HadGEM2_pr_grid_areas = iris.analysis.cartography.area_weights(HadGEM2_pr)
    IPSLG_pr_grid_areas = iris.analysis.cartography.area_weights(IPSLG_pr)
    MIROCG_pr_grid_areas = iris.analysis.cartography.area_weights(MIROCG_pr)
    MPI_pr_grid_areas = iris.analysis.cartography.area_weights(MPI_pr)
    NorESM1_pr_grid_areas = iris.analysis.cartography.area_weights(NorESM1_pr)
    
    CRU_pr_grid_areas = iris.analysis.cartography.area_weights(CRU_pr)
    UDel_pr_grid_areas = iris.analysis.cartography.area_weights(UDel_pr)
    GPCC_pr_grid_areas = iris.analysis.cartography.area_weights(GPCC_pr)
    
    CCCma_tas_grid_areas = iris.analysis.cartography.area_weights(CCCma_tas)
    CLMcom_tas_grid_areas = iris.analysis.cartography.area_weights(CLMcom_tas)
    DMI_tas_grid_areas = iris.analysis.cartography.area_weights(DMI_tas)
    KNMI_tas_grid_areas = iris.analysis.cartography.area_weights(KNMI_tas)
    MPIE_tas_grid_areas = iris.analysis.cartography.area_weights(MPIE_tas)
    SMHI_tas_grid_areas = iris.analysis.cartography.area_weights(SMHI_tas)
    
    CCCmaCanRCM_tas_grid_areas = iris.analysis.cartography.area_weights(CCCmaCanRCM_tas)
    CCCmaSMHI_tas_grid_areas = iris.analysis.cartography.area_weights(CCCmaSMHI_tas)
    CNRM_tas_grid_areas = iris.analysis.cartography.area_weights(CNRM_tas)
    CNRMSMHI_tas_grid_areas = iris.analysis.cartography.area_weights(CNRMSMHI_tas)
    CSIRO_tas_grid_areas = iris.analysis.cartography.area_weights(CSIRO_tas)
    ICHECDMI_tas_grid_areas = iris.analysis.cartography.area_weights(ICHECDMI_tas)
    ICHECCCLM_tas_grid_areas = iris.analysis.cartography.area_weights(ICHECCCLM_tas)
    ICHECKNMI_tas_grid_areas = iris.analysis.cartography.area_weights(ICHECKNMI_tas)
    ICHECMPI_tas_grid_areas = iris.analysis.cartography.area_weights(ICHECMPI_tas)
    ICHECSMHI_tas_grid_areas = iris.analysis.cartography.area_weights(ICHECSMHI_tas)
    IPSL_tas_grid_areas = iris.analysis.cartography.area_weights(IPSL_tas)
    MIROC_tas_grid_areas = iris.analysis.cartography.area_weights(MIROC_tas)
    MOHCCCLM_tas_grid_areas = iris.analysis.cartography.area_weights(MOHCCCLM_tas)
    MOHCKNMI_tas_grid_areas = iris.analysis.cartography.area_weights(MOHCKNMI_tas)
    MOHCSMHI_tas_grid_areas = iris.analysis.cartography.area_weights(MOHCSMHI_tas)
    MPICCLM_tas_grid_areas = iris.analysis.cartography.area_weights(MPICCLM_tas)
    MPIREMO_tas_grid_areas = iris.analysis.cartography.area_weights(MPIREMO_tas)
    MPISMHI_tas_grid_areas = iris.analysis.cartography.area_weights(MPISMHI_tas)
    NCCDMI_tas_grid_areas = iris.analysis.cartography.area_weights(NCCDMI_tas)
    NCCSMHI_tas_grid_areas = iris.analysis.cartography.area_weights(NCCSMHI_tas)
    NOAA_tas_grid_areas = iris.analysis.cartography.area_weights(NOAA_tas)
    
    CanESM2_tas_grid_areas = iris.analysis.cartography.area_weights(CanESM2_tas)
    CNRMG_tas_grid_areas = iris.analysis.cartography.area_weights(CNRMG_tas)
    MK3_tas_grid_areas = iris.analysis.cartography.area_weights(MK3_tas)
    EARTH_tas_grid_areas = iris.analysis.cartography.area_weights(EARTH_tas)
    EARTH3_tas_grid_areas = iris.analysis.cartography.area_weights(EARTH3_tas)
    GFDL_tas_grid_areas = iris.analysis.cartography.area_weights(GFDL_tas)
    HadGEM2_tas_grid_areas = iris.analysis.cartography.area_weights(HadGEM2_tas)
    IPSLG_tas_grid_areas = iris.analysis.cartography.area_weights(IPSLG_tas)
    MIROCG_tas_grid_areas = iris.analysis.cartography.area_weights(MIROCG_tas)
    MPI_tas_grid_areas = iris.analysis.cartography.area_weights(MPI_tas)
    NorESM1_tas_grid_areas = iris.analysis.cartography.area_weights(NorESM1_tas)
    
    CRU_tas_grid_areas = iris.analysis.cartography.area_weights(CRU_tas)
    UDel_tas_grid_areas = iris.analysis.cartography.area_weights(UDel_tas)
    
    CCCma_tasmin_grid_areas = iris.analysis.cartography.area_weights(CCCma_tasmin)
    CLMcom_tasmin_grid_areas = iris.analysis.cartography.area_weights(CLMcom_tasmin)
    DMI_tasmin_grid_areas = iris.analysis.cartography.area_weights(DMI_tasmin)
    KNMI_tasmin_grid_areas = iris.analysis.cartography.area_weights(KNMI_tasmin)
    MPIE_tasmin_grid_areas = iris.analysis.cartography.area_weights(MPIE_tasmin)
    SMHI_tasmin_grid_areas = iris.analysis.cartography.area_weights(SMHI_tasmin)
    
    CCCmaCanRCM_tasmin_grid_areas = iris.analysis.cartography.area_weights(CCCmaCanRCM_tasmin)
    CCCmaSMHI_tasmin_grid_areas = iris.analysis.cartography.area_weights(CCCmaSMHI_tasmin)
    CNRM_tasmin_grid_areas = iris.analysis.cartography.area_weights(CNRM_tasmin)
    CNRMSMHI_tasmin_grid_areas = iris.analysis.cartography.area_weights(CNRMSMHI_tasmin)
    CSIRO_tasmin_grid_areas = iris.analysis.cartography.area_weights(CSIRO_tasmin)
    ICHECDMI_tasmin_grid_areas = iris.analysis.cartography.area_weights(ICHECDMI_tasmin)
    ICHECCCLM_tasmin_grid_areas = iris.analysis.cartography.area_weights(ICHECCCLM_tasmin)
    ICHECKNMI_tasmin_grid_areas = iris.analysis.cartography.area_weights(ICHECKNMI_tasmin)
    ICHECMPI_tasmin_grid_areas = iris.analysis.cartography.area_weights(ICHECMPI_tasmin)
    ICHECSMHI_tasmin_grid_areas = iris.analysis.cartography.area_weights(ICHECSMHI_tasmin)
    IPSL_tasmin_grid_areas = iris.analysis.cartography.area_weights(IPSL_tasmin)
    MIROC_tasmin_grid_areas = iris.analysis.cartography.area_weights(MIROC_tasmin)
    MOHCCCLM_tasmin_grid_areas = iris.analysis.cartography.area_weights(MOHCCCLM_tasmin)
    MOHCKNMI_tasmin_grid_areas = iris.analysis.cartography.area_weights(MOHCKNMI_tasmin)
    MOHCSMHI_tasmin_grid_areas = iris.analysis.cartography.area_weights(MOHCSMHI_tasmin)
    MPICCLM_tasmin_grid_areas = iris.analysis.cartography.area_weights(MPICCLM_tasmin)
    MPIREMO_tasmin_grid_areas = iris.analysis.cartography.area_weights(MPIREMO_tasmin)
    MPISMHI_tasmin_grid_areas = iris.analysis.cartography.area_weights(MPISMHI_tasmin)
    NCCDMI_tasmin_grid_areas = iris.analysis.cartography.area_weights(NCCDMI_tasmin)
    NCCSMHI_tasmin_grid_areas = iris.analysis.cartography.area_weights(NCCSMHI_tasmin)
    NOAA_tasmin_grid_areas = iris.analysis.cartography.area_weights(NOAA_tasmin)
    
    CanESM2_tasmin_grid_areas = iris.analysis.cartography.area_weights(CanESM2_tasmin)
    CNRMG_tasmin_grid_areas = iris.analysis.cartography.area_weights(CNRMG_tasmin)
    MK3_tasmin_grid_areas = iris.analysis.cartography.area_weights(MK3_tasmin)
    EARTH_tasmin_grid_areas = iris.analysis.cartography.area_weights(EARTH_tasmin)
    EARTH3_tasmin_grid_areas = iris.analysis.cartography.area_weights(EARTH3_tasmin)
    GFDL_tasmin_grid_areas = iris.analysis.cartography.area_weights(GFDL_tasmin)
    HadGEM2_tasmin_grid_areas = iris.analysis.cartography.area_weights(HadGEM2_tasmin)
    IPSLG_tasmin_grid_areas = iris.analysis.cartography.area_weights(IPSLG_tasmin)
    MIROCG_tasmin_grid_areas = iris.analysis.cartography.area_weights(MIROCG_tasmin)
    MPI_tasmin_grid_areas = iris.analysis.cartography.area_weights(MPI_tasmin)
    NorESM1_tasmin_grid_areas = iris.analysis.cartography.area_weights(NorESM1_tasmin)
    
    CRU_tasmin_grid_areas = iris.analysis.cartography.area_weights(CRU_tasmin)
    
    CCCma_tasmax_grid_areas = iris.analysis.cartography.area_weights(CCCma_tasmax)
    CLMcom_tasmax_grid_areas = iris.analysis.cartography.area_weights(CLMcom_tasmax)
    DMI_tasmax_grid_areas = iris.analysis.cartography.area_weights(DMI_tasmax)
    KNMI_tasmax_grid_areas = iris.analysis.cartography.area_weights(KNMI_tasmax)
    MPIE_tasmax_grid_areas = iris.analysis.cartography.area_weights(MPIE_tasmax)
    SMHI_tasmax_grid_areas = iris.analysis.cartography.area_weights(SMHI_tasmax)
    
    CCCmaCanRCM_tasmax_grid_areas = iris.analysis.cartography.area_weights(CCCmaCanRCM_tasmax)
    CCCmaSMHI_tasmax_grid_areas = iris.analysis.cartography.area_weights(CCCmaSMHI_tasmax)
    CNRM_tasmax_grid_areas = iris.analysis.cartography.area_weights(CNRM_tasmax)
    CNRMSMHI_tasmax_grid_areas = iris.analysis.cartography.area_weights(CNRMSMHI_tasmax)
    CSIRO_tasmax_grid_areas = iris.analysis.cartography.area_weights(CSIRO_tasmax)
    ICHECDMI_tasmax_grid_areas = iris.analysis.cartography.area_weights(ICHECDMI_tasmax)
    ICHECCCLM_tasmax_grid_areas = iris.analysis.cartography.area_weights(ICHECCCLM_tasmax)
    ICHECKNMI_tasmax_grid_areas = iris.analysis.cartography.area_weights(ICHECKNMI_tasmax)
    ICHECMPI_tasmax_grid_areas = iris.analysis.cartography.area_weights(ICHECMPI_tasmax)
    ICHECSMHI_tasmax_grid_areas = iris.analysis.cartography.area_weights(ICHECSMHI_tasmax)
    IPSL_tasmax_grid_areas = iris.analysis.cartography.area_weights(IPSL_tasmax)
    MIROC_tasmax_grid_areas = iris.analysis.cartography.area_weights(MIROC_tasmax)
    MOHCCCLM_tasmax_grid_areas = iris.analysis.cartography.area_weights(MOHCCCLM_tasmax)
    MOHCKNMI_tasmax_grid_areas = iris.analysis.cartography.area_weights(MOHCKNMI_tasmax)
    MOHCSMHI_tasmax_grid_areas = iris.analysis.cartography.area_weights(MOHCSMHI_tasmax)
    MPICCLM_tasmax_grid_areas = iris.analysis.cartography.area_weights(MPICCLM_tasmax)
    MPIREMO_tasmax_grid_areas = iris.analysis.cartography.area_weights(MPIREMO_tasmax)
    MPISMHI_tasmax_grid_areas = iris.analysis.cartography.area_weights(MPISMHI_tasmax)
    NCCDMI_tasmax_grid_areas = iris.analysis.cartography.area_weights(NCCDMI_tasmax)
    NCCSMHI_tasmax_grid_areas = iris.analysis.cartography.area_weights(NCCSMHI_tasmax)
    NOAA_tasmax_grid_areas = iris.analysis.cartography.area_weights(NOAA_tasmax)
    
    CanESM2_tasmax_grid_areas = iris.analysis.cartography.area_weights(CanESM2_tasmax)
    CNRMG_tasmax_grid_areas = iris.analysis.cartography.area_weights(CNRMG_tasmax)
    MK3_tasmax_grid_areas = iris.analysis.cartography.area_weights(MK3_tasmax)
    EARTH_tasmax_grid_areas = iris.analysis.cartography.area_weights(EARTH_tasmax)
    EARTH3_tasmax_grid_areas = iris.analysis.cartography.area_weights(EARTH3_tasmax)
    GFDL_tasmax_grid_areas = iris.analysis.cartography.area_weights(GFDL_tasmax)
    HadGEM2_tasmax_grid_areas = iris.analysis.cartography.area_weights(HadGEM2_tasmax)
    IPSLG_tasmax_grid_areas = iris.analysis.cartography.area_weights(IPSLG_tasmax)
    MIROCG_tasmax_grid_areas = iris.analysis.cartography.area_weights(MIROCG_tasmax)
    MPI_tasmax_grid_areas = iris.analysis.cartography.area_weights(MPI_tasmax)
    NorESM1_tasmax_grid_areas = iris.analysis.cartography.area_weights(NorESM1_tasmax)
    
    CRU_tasmax_grid_areas = iris.analysis.cartography.area_weights(CRU_tasmax)
    
    #We want to plot the mean for the whole region so we need a mean of all the lats and lons
    CCCma_pr_mean = CCCma_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CCCma_pr_grid_areas)                                        
    CLMcom_pr_mean = CLMcom_pr.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=CLMcom_pr_grid_areas)
    DMI_pr_mean = DMI_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=DMI_pr_grid_areas)     
    KNMI_pr_mean = KNMI_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=KNMI_pr_grid_areas)
    MPIE_pr_mean = MPIE_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIE_pr_grid_areas)
    SMHI_pr_mean = SMHI_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=SMHI_pr_grid_areas)
    
    CCCmaCanRCM_pr_mean = CCCmaCanRCM_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CCCmaCanRCM_pr_grid_areas) 
    CCCmaSMHI_pr_mean = CCCmaSMHI_pr.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=CCCmaSMHI_pr_grid_areas)
    CNRM_pr_mean = CNRM_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRM_pr_grid_areas)                                               
    CNRMSMHI_pr_mean = CNRMSMHI_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMSMHI_pr_grid_areas)  
    CSIRO_pr_mean = CSIRO_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CSIRO_pr_grid_areas)
    ICHECDMI_pr_mean = ICHECDMI_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECDMI_pr_grid_areas) 
    ICHECCCLM_pr_mean = ICHECCCLM_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECCCLM_pr_grid_areas)
    ICHECKNMI_pr_mean = ICHECKNMI_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECKNMI_pr_grid_areas)
    ICHECMPI_pr_mean = ICHECMPI_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECMPI_pr_grid_areas)
    ICHECSMHI_pr_mean = ICHECSMHI_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECSMHI_pr_grid_areas)
    IPSL_pr_mean = IPSL_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=IPSL_pr_grid_areas)
    MIROC_pr_mean = MIROC_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MIROC_pr_grid_areas)
    MOHCCCLM_pr_mean = MOHCCCLM_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCCCLM_pr_grid_areas)
    MOHCKNMI_pr_mean = MOHCKNMI_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCKNMI_pr_grid_areas)
    MOHCSMHI_pr_mean = MOHCSMHI_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCSMHI_pr_grid_areas)
    MPICCLM_pr_mean = MPICCLM_pr.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=MPICCLM_pr_grid_areas)                                              
    MPIREMO_pr_mean = MPIREMO_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIREMO_pr_grid_areas)                                              
    MPISMHI_pr_mean = MPISMHI_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPISMHI_pr_grid_areas)
    NCCDMI_pr_mean = NCCDMI_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NCCDMI_pr_grid_areas)
    NCCSMHI_pr_mean = NCCSMHI_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NCCSMHI_pr_grid_areas) 
    NOAA_pr_mean = NOAA_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NOAA_pr_grid_areas)
    
    CanESM2_pr_mean = CanESM2_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CanESM2_pr_grid_areas)                        
    CNRMG_pr_mean = CNRMG_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMG_pr_grid_areas) 
    MK3_pr_mean = MK3_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MK3_pr_grid_areas) 
    EARTH_pr_mean = EARTH_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=EARTH_pr_grid_areas)  
    EARTH3_pr_mean = EARTH3_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=EARTH3_pr_grid_areas)  
    GFDL_pr_mean = GFDL_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=GFDL_pr_grid_areas)
    HadGEM2_pr_mean = HadGEM2_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=HadGEM2_pr_grid_areas)
    IPSLG_pr_mean = IPSLG_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=IPSLG_pr_grid_areas)
    MIROCG_pr_mean = MIROCG_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MIROCG_pr_grid_areas)         
    MPI_pr_mean = MPI_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPI_pr_grid_areas)
    NorESM1_pr_mean = NorESM1_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NorESM1_pr_grid_areas)
    
    CRU_pr_mean = CRU_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CRU_pr_grid_areas) 
    UDel_pr_mean = UDel_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=UDel_pr_grid_areas)
    GPCC_pr_mean = GPCC_pr.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=GPCC_pr_grid_areas)
    
    CCCma_tas_mean = CCCma_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CCCma_tas_grid_areas)                                        
    CLMcom_tas_mean = CLMcom_tas.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=CLMcom_tas_grid_areas)
    DMI_tas_mean = DMI_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=DMI_tas_grid_areas)     
    KNMI_tas_mean = KNMI_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=KNMI_tas_grid_areas)
    MPIE_tas_mean = MPIE_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIE_tas_grid_areas)
    SMHI_tas_mean = SMHI_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=SMHI_tas_grid_areas)
    
    CCCmaCanRCM_tas_mean = CCCmaCanRCM_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CCCmaCanRCM_tas_grid_areas) 
    CCCmaSMHI_tas_mean = CCCmaSMHI_tas.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=CCCmaSMHI_tas_grid_areas)
    CNRM_tas_mean = CNRM_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRM_tas_grid_areas)                                               
    CNRMSMHI_tas_mean = CNRMSMHI_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMSMHI_tas_grid_areas)  
    CSIRO_tas_mean = CSIRO_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CSIRO_tas_grid_areas)
    ICHECDMI_tas_mean = ICHECDMI_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECDMI_tas_grid_areas) 
    ICHECCCLM_tas_mean = ICHECCCLM_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECCCLM_tas_grid_areas)
    ICHECKNMI_tas_mean = ICHECKNMI_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECKNMI_tas_grid_areas)
    ICHECMPI_tas_mean = ICHECMPI_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECMPI_tas_grid_areas)
    ICHECSMHI_tas_mean = ICHECSMHI_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECSMHI_tas_grid_areas)
    IPSL_tas_mean = IPSL_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=IPSL_tas_grid_areas)
    MIROC_tas_mean = MIROC_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MIROC_tas_grid_areas)
    MOHCCCLM_tas_mean = MOHCCCLM_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCCCLM_tas_grid_areas)
    MOHCKNMI_tas_mean = MOHCKNMI_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCKNMI_tas_grid_areas)
    MOHCSMHI_tas_mean = MOHCSMHI_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCSMHI_tas_grid_areas)
    MPICCLM_tas_mean = MPICCLM_tas.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=MPICCLM_tas_grid_areas)                                              
    MPIREMO_tas_mean = MPIREMO_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIREMO_tas_grid_areas)                                              
    MPISMHI_tas_mean = MPISMHI_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPISMHI_tas_grid_areas)
    NCCDMI_tas_mean = NCCDMI_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NCCDMI_tas_grid_areas)
    NCCSMHI_tas_mean = NCCSMHI_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NCCSMHI_tas_grid_areas) 
    NOAA_tas_mean = NOAA_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NOAA_tas_grid_areas)
    
    CanESM2_tas_mean = CanESM2_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CanESM2_tas_grid_areas)                        
    CNRMG_tas_mean = CNRMG_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMG_tas_grid_areas) 
    MK3_tas_mean = MK3_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MK3_tas_grid_areas) 
    EARTH_tas_mean = EARTH_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=EARTH_tas_grid_areas)  
    EARTH3_tas_mean = EARTH3_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=EARTH3_tas_grid_areas)  
    GFDL_tas_mean = GFDL_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=GFDL_tas_grid_areas)
    HadGEM2_tas_mean = HadGEM2_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=HadGEM2_tas_grid_areas)
    IPSLG_tas_mean = IPSLG_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=IPSLG_tas_grid_areas)
    MIROCG_tas_mean = MIROCG_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MIROCG_tas_grid_areas)         
    MPI_tas_mean = MPI_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPI_tas_grid_areas)
    NorESM1_tas_mean = NorESM1_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NorESM1_tas_grid_areas)
    
    CRU_tas_mean = CRU_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CRU_tas_grid_areas) 
    UDel_tas_mean = UDel_tas.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=UDel_tas_grid_areas)
    CCCma_tasmin_mean = CCCma_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CCCma_tasmin_grid_areas)                                        
    CLMcom_tasmin_mean = CLMcom_tasmin.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=CLMcom_tasmin_grid_areas)
    DMI_tasmin_mean = DMI_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=DMI_tasmin_grid_areas)     
    KNMI_tasmin_mean = KNMI_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=KNMI_tasmin_grid_areas)
    MPIE_tasmin_mean = MPIE_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIE_tasmin_grid_areas)
    SMHI_tasmin_mean = SMHI_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=SMHI_tasmin_grid_areas)
    
    CCCmaCanRCM_tasmin_mean = CCCmaCanRCM_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CCCmaCanRCM_tasmin_grid_areas) 
    CCCmaSMHI_tasmin_mean = CCCmaSMHI_tasmin.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=CCCmaSMHI_tasmin_grid_areas)
    CNRM_tasmin_mean = CNRM_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRM_tasmin_grid_areas)                                               
    CNRMSMHI_tasmin_mean = CNRMSMHI_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMSMHI_tasmin_grid_areas)  
    CSIRO_tasmin_mean = CSIRO_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CSIRO_tasmin_grid_areas)
    ICHECDMI_tasmin_mean = ICHECDMI_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECDMI_tasmin_grid_areas) 
    ICHECCCLM_tasmin_mean = ICHECCCLM_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECCCLM_tasmin_grid_areas)
    ICHECKNMI_tasmin_mean = ICHECKNMI_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECKNMI_tasmin_grid_areas)
    ICHECMPI_tasmin_mean = ICHECMPI_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECMPI_tasmin_grid_areas)
    ICHECSMHI_tasmin_mean = ICHECSMHI_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECSMHI_tasmin_grid_areas)
    IPSL_tasmin_mean = IPSL_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=IPSL_tasmin_grid_areas)
    MIROC_tasmin_mean = MIROC_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MIROC_tasmin_grid_areas)
    MOHCCCLM_tasmin_mean = MOHCCCLM_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCCCLM_tasmin_grid_areas)
    MOHCKNMI_tasmin_mean = MOHCKNMI_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCKNMI_tasmin_grid_areas)
    MOHCSMHI_tasmin_mean = MOHCSMHI_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCSMHI_tasmin_grid_areas)
    MPICCLM_tasmin_mean = MPICCLM_tasmin.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=MPICCLM_tasmin_grid_areas)                                              
    MPIREMO_tasmin_mean = MPIREMO_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIREMO_tasmin_grid_areas)                                              
    MPISMHI_tasmin_mean = MPISMHI_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPISMHI_tasmin_grid_areas)
    NCCDMI_tasmin_mean = NCCDMI_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NCCDMI_tasmin_grid_areas)
    NCCSMHI_tasmin_mean = NCCSMHI_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NCCSMHI_tasmin_grid_areas) 
    NOAA_tasmin_mean = NOAA_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NOAA_tasmin_grid_areas)
    
    CanESM2_tasmin_mean = CanESM2_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CanESM2_tasmin_grid_areas)                        
    CNRMG_tasmin_mean = CNRMG_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMG_tasmin_grid_areas) 
    MK3_tasmin_mean = MK3_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MK3_tasmin_grid_areas) 
    EARTH_tasmin_mean = EARTH_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=EARTH_tasmin_grid_areas)  
    EARTH3_tasmin_mean = EARTH3_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=EARTH3_tasmin_grid_areas)  
    GFDL_tasmin_mean = GFDL_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=GFDL_tasmin_grid_areas)
    HadGEM2_tasmin_mean = HadGEM2_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=HadGEM2_tasmin_grid_areas)
    IPSLG_tasmin_mean = IPSLG_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=IPSLG_tasmin_grid_areas)
    MIROCG_tasmin_mean = MIROCG_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MIROCG_tasmin_grid_areas)         
    MPI_tasmin_mean = MPI_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPI_tasmin_grid_areas)
    NorESM1_tasmin_mean = NorESM1_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NorESM1_tasmin_grid_areas)
    
    CRU_tasmin_mean = CRU_tasmin.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CRU_tasmin_grid_areas) 
   
    CCCma_tasmax_mean = CCCma_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CCCma_tasmax_grid_areas)                                        
    CLMcom_tasmax_mean = CLMcom_tasmax.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=CLMcom_tasmax_grid_areas)
    DMI_tasmax_mean = DMI_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=DMI_tasmax_grid_areas)     
    KNMI_tasmax_mean = KNMI_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=KNMI_tasmax_grid_areas)
    MPIE_tasmax_mean = MPIE_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIE_tasmax_grid_areas)
    SMHI_tasmax_mean = SMHI_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=SMHI_tasmax_grid_areas)
    
    CCCmaCanRCM_tasmax_mean = CCCmaCanRCM_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CCCmaCanRCM_tasmax_grid_areas) 
    CCCmaSMHI_tasmax_mean = CCCmaSMHI_tasmax.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=CCCmaSMHI_tasmax_grid_areas)
    CNRM_tasmax_mean = CNRM_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRM_tasmax_grid_areas)                                               
    CNRMSMHI_tasmax_mean = CNRMSMHI_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMSMHI_tasmax_grid_areas)  
    CSIRO_tasmax_mean = CSIRO_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CSIRO_tasmax_grid_areas)
    ICHECDMI_tasmax_mean = ICHECDMI_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECDMI_tasmax_grid_areas) 
    ICHECCCLM_tasmax_mean = ICHECCCLM_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECCCLM_tasmax_grid_areas)
    ICHECKNMI_tasmax_mean = ICHECKNMI_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECKNMI_tasmax_grid_areas)
    ICHECMPI_tasmax_mean = ICHECMPI_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECMPI_tasmax_grid_areas)
    ICHECSMHI_tasmax_mean = ICHECSMHI_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=ICHECSMHI_tasmax_grid_areas)
    IPSL_tasmax_mean = IPSL_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=IPSL_tasmax_grid_areas)
    MIROC_tasmax_mean = MIROC_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MIROC_tasmax_grid_areas)
    MOHCCCLM_tasmax_mean = MOHCCCLM_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCCCLM_tasmax_grid_areas)
    MOHCKNMI_tasmax_mean = MOHCKNMI_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCKNMI_tasmax_grid_areas)
    MOHCSMHI_tasmax_mean = MOHCSMHI_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MOHCSMHI_tasmax_grid_areas)
    MPICCLM_tasmax_mean = MPICCLM_tasmax.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=MPICCLM_tasmax_grid_areas)                                              
    MPIREMO_tasmax_mean = MPIREMO_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPIREMO_tasmax_grid_areas)                                              
    MPISMHI_tasmax_mean = MPISMHI_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPISMHI_tasmax_grid_areas)
    NCCDMI_tasmax_mean = NCCDMI_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NCCDMI_tasmax_grid_areas)
    NCCSMHI_tasmax_mean = NCCSMHI_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NCCSMHI_tasmax_grid_areas) 
    NOAA_tasmax_mean = NOAA_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NOAA_tasmax_grid_areas)
    
    CanESM2_tasmax_mean = CanESM2_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CanESM2_tasmax_grid_areas)                        
    CNRMG_tasmax_mean = CNRMG_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CNRMG_tasmax_grid_areas) 
    MK3_tasmax_mean = MK3_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MK3_tasmax_grid_areas) 
    EARTH_tasmax_mean = EARTH_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=EARTH_tasmax_grid_areas)  
    EARTH3_tasmax_mean = EARTH3_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=EARTH3_tasmax_grid_areas)  
    GFDL_tasmax_mean = GFDL_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=GFDL_tasmax_grid_areas)
    HadGEM2_tasmax_mean = HadGEM2_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=HadGEM2_tasmax_grid_areas)
    IPSLG_tasmax_mean = IPSLG_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=IPSLG_tasmax_grid_areas)
    MIROCG_tasmax_mean = MIROCG_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MIROCG_tasmax_grid_areas)         
    MPI_tasmax_mean = MPI_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=MPI_tasmax_grid_areas)
    NorESM1_tasmax_mean = NorESM1_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=NorESM1_tasmax_grid_areas)
    
    CRU_tasmax_mean = CRU_tasmax.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=CRU_tasmax_grid_areas) 
   
    #create averages
    Average_Eraint_pr = (CCCma_pr_mean.data + CLMcom_pr_mean.data + DMI_pr_mean.data + KNMI_pr_mean.data + MPIE_pr_mean.data + SMHI_pr_mean.data)/6.
    Average_Eraint_tas = (CCCma_tas_mean.data + CLMcom_tas_mean.data + DMI_tas_mean.data + KNMI_tas_mean.data + MPIE_tas_mean.data + SMHI_tas_mean.data)/6.
    Average_Eraint_tasmin = (CCCma_tasmin_mean.data + CLMcom_tasmin_mean.data + DMI_tasmin_mean.data + KNMI_tasmin_mean.data + MPIE_tasmin_mean.data + SMHI_tasmin_mean.data)/6.
    Average_Eraint_tasmax = (CCCma_tasmax_mean.data + CLMcom_tasmax_mean.data + DMI_tasmax_mean.data + KNMI_tasmax_mean.data + MPIE_tasmax_mean.data + SMHI_tasmax_mean.data)/6.
    
    Average_RCM_pr = (CCCmaCanRCM_pr_mean.data + CCCmaSMHI_pr_mean.data + CNRM_pr_mean.data + CNRMSMHI_pr_mean.data + CSIRO_pr_mean.data + ICHECDMI_pr_mean.data + ICHECCCLM_pr_mean.data + ICHECKNMI_pr_mean.data + ICHECMPI_pr_mean.data + ICHECSMHI_pr_mean.data + IPSL_pr_mean.data + MIROC_pr_mean.data + MOHCCCLM_pr_mean.data + MOHCKNMI_pr_mean.data + MOHCSMHI_pr_mean.data + MPICCLM_pr_mean.data + MPIREMO_pr_mean.data + MPISMHI_pr_mean.data + NCCDMI_pr_mean.data + NCCSMHI_pr_mean.data + NOAA_pr_mean.data)/21.
    Average_RCM_tas = (CCCmaCanRCM_tas_mean.data + CCCmaSMHI_tas_mean.data + CNRM_tas_mean.data + CNRMSMHI_tas_mean.data + CSIRO_tas_mean.data + ICHECDMI_tas_mean.data + ICHECCCLM_tas_mean.data + ICHECKNMI_tas_mean.data + ICHECMPI_tas_mean.data + ICHECSMHI_tas_mean.data + IPSL_tas_mean.data + MIROC_tas_mean.data + MOHCCCLM_tas_mean.data + MOHCKNMI_tas_mean.data + MOHCSMHI_tas_mean.data + MPICCLM_tas_mean.data + MPIREMO_tas_mean.data + MPISMHI_tas_mean.data + NCCDMI_tas_mean.data + NCCSMHI_tas_mean.data + NOAA_tas_mean.data)/21.
    Average_RCM_tasmin = (CCCmaCanRCM_tasmin_mean.data + CCCmaSMHI_tasmin_mean.data + CNRM_tasmin_mean.data + CNRMSMHI_tasmin_mean.data + CSIRO_tasmin_mean.data + ICHECDMI_tasmin_mean.data + ICHECCCLM_tasmin_mean.data + ICHECKNMI_tasmin_mean.data + ICHECMPI_tasmin_mean.data + ICHECSMHI_tasmin_mean.data + IPSL_tasmin_mean.data + MIROC_tasmin_mean.data + MOHCCCLM_tasmin_mean.data + MOHCKNMI_tasmin_mean.data + MOHCSMHI_tasmin_mean.data + MPICCLM_tasmin_mean.data + MPIREMO_tasmin_mean.data + MPISMHI_tasmin_mean.data + NCCDMI_tasmin_mean.data + NCCSMHI_tasmin_mean.data + NOAA_tasmin_mean.data)/21.
    Average_RCM_tasmax = (CCCmaCanRCM_tasmax_mean.data + CCCmaSMHI_tasmax_mean.data + CNRM_tasmax_mean.data + CNRMSMHI_tasmax_mean.data + CSIRO_tasmax_mean.data + ICHECDMI_tasmax_mean.data + ICHECCCLM_tasmax_mean.data + ICHECKNMI_tasmax_mean.data + ICHECMPI_tasmax_mean.data + ICHECSMHI_tasmax_mean.data + IPSL_tasmax_mean.data + MIROC_tasmax_mean.data + MOHCCCLM_tasmax_mean.data + MOHCKNMI_tasmax_mean.data + MOHCSMHI_tasmax_mean.data + MPICCLM_tasmax_mean.data + MPIREMO_tasmax_mean.data + MPISMHI_tasmax_mean.data + NCCDMI_tasmax_mean.data + NCCSMHI_tasmax_mean.data + NOAA_tasmax_mean.data)/21.
    
    Average_GCM_pr = (CanESM2_pr_mean.data + CNRMG_pr_mean.data + MK3_pr_mean.data + EARTH_pr_mean.data + EARTH3_pr_mean.data + GFDL_pr_mean.data + HadGEM2_pr_mean.data + IPSLG_pr_mean.data + MIROCG_pr_mean.data + MPI_pr_mean.data + NorESM1_pr_mean.data)/11.
    Average_GCM_tas = (CanESM2_tas_mean.data + CNRMG_tas_mean.data + MK3_tas_mean.data + EARTH_tas_mean.data + EARTH3_tas_mean.data + GFDL_tas_mean.data + HadGEM2_tas_mean.data + IPSLG_tas_mean.data + MIROCG_tas_mean.data + MPI_tas_mean.data + NorESM1_tas_mean.data)/11.
    Average_GCM_tasmin = (CanESM2_tasmin_mean.data + CNRMG_tasmin_mean.data + MK3_tasmin_mean.data + EARTH_tasmin_mean.data + EARTH3_tasmin_mean.data + GFDL_tasmin_mean.data + HadGEM2_tasmin_mean.data + IPSLG_tasmin_mean.data + MIROCG_tasmin_mean.data + MPI_tasmin_mean.data + NorESM1_tasmin_mean.data)/11.
    Average_GCM_tasmax = (CanESM2_tasmax_mean.data + CNRMG_tasmax_mean.data + MK3_tasmax_mean.data + EARTH_tasmax_mean.data + EARTH3_tasmax_mean.data + GFDL_tasmax_mean.data + HadGEM2_tasmax_mean.data + IPSLG_tasmax_mean.data + MIROCG_tasmax_mean.data + MPI_tasmax_mean.data + NorESM1_tasmax_mean.data)/11.
    
    Obs_pr = (CRU_pr_mean.data + UDel_pr_mean.data + GPCC_pr_mean.data)/3
    Obs_tas = (CRU_tas_mean.data + UDel_tas_mean.data)/2
    
    
    
    #-------------------------------------------------------------------------
    #PART 7: PRINT DATA
    print CCCma_pr_mean.coord('year')
    print CCCma_pr_mean.coord('month_number')
    print CCCmaCanRCM_pr_mean.coord('year')
    print CCCmaCanRCM_pr_mean.coord('month_number')
    
    import csv
    with open('output_print_pr.csv', 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        
        writer.writerow(['Parameter', 'Means'])
        
        writer.writerow(["CCCma_pr_mean"] + CCCma_pr_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["CLMcom_pr_mean"] + CLMcom_pr_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["DMI_pr_mean"] + DMI_pr_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["KNMI_pr_mean"] + KNMI_pr_mean.data.flatten().astype(np.str).tolist())   
        writer.writerow(["MPIE_pr_mean"] + MPIE_pr_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["SMHI_pr_mean"] + SMHI_pr_mean.data.flatten().astype(np.str).tolist()) 
        writer.writerow(["Average_Eraint_pr"] + Average_Eraint_pr.data.flatten().astype(np.str).tolist()) 
        
        writer.writerow(["CCCmaCanRCM_pr_mean"] + CCCmaCanRCM_pr_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["CCCmaSMHI_pr_mean"] + CCCmaSMHI_pr_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["CNRM_pr_mean"] + CNRM_pr_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["CNRMSMHI_pr_mean"] + CNRMSMHI_pr_mean.data.flatten().astype(np.str).tolist())   
        writer.writerow(["CSIRO_pr_mean"] + CSIRO_pr_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["ICHECDMI_pr_mean"] + ICHECDMI_pr_mean.data.flatten().astype(np.str).tolist())       
        writer.writerow(["ICHECCCLM_pr_mean"] + ICHECCCLM_pr_mean.data.flatten().astype(np.str).tolist()) 
        writer.writerow(["ICHECKNMI_pr_mean"] + ICHECKNMI_pr_mean.data.flatten().astype(np.str).tolist()) 
        writer.writerow(["ICHECMPI_pr_mean"] + ICHECMPI_pr_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["ICHECSMHI_pr_mean"] + ICHECSMHI_pr_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["IPSL_pr_mean"] + IPSL_pr_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["MIROC_pr_mean"] + MIROC_pr_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["MOHCCCLM_pr_mean"] + MOHCCCLM_pr_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["MOHCKNMI_pr_mean"] + MOHCKNMI_pr_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["MOHCSMHI_pr_mean"] + MOHCSMHI_pr_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["MPICCLM_pr_mean"] + MPICCLM_pr_mean.data.flatten().astype(np.str).tolist())    
        writer.writerow(["MPIREMO_pr_mean"] + MPIREMO_pr_mean.data.flatten().astype(np.str).tolist()) 
        writer.writerow(["MPISMHI_pr_mean"] + MPISMHI_pr_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["NCCDMI_pr_mean"] + NCCDMI_pr_mean.data.flatten().astype(np.str).tolist())   
        writer.writerow(["NCCSMHI_pr_mean"] + NCCSMHI_pr_mean.data.flatten().astype(np.str).tolist())    
        writer.writerow(["NOAA_pr_mean"] + NOAA_pr_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["Average_RCM_pr"] + Average_RCM_pr.data.flatten().astype(np.str).tolist())
        
        writer.writerow(["CanESM2_pr_mean"] + CanESM2_pr_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["CNRMG_pr_mean"] + CNRMG_pr_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["MK3_pr_mean"] + MK3_pr_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["EARTH_pr_mean"] + EARTH_pr_mean.data.flatten().astype(np.str).tolist())   
        writer.writerow(["EARTH3_pr_mean"] + EARTH3_pr_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["GFDL_pr_mean"] + GFDL_pr_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["HadGEM2_pr_mean"] + HadGEM2_pr_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["IPSLG_pr_mean"] + IPSLG_pr_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["MIROCG_pr_mean"] + MIROCG_pr_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["MPI_pr_mean"] + MPI_pr_mean.data.flatten().astype(np.str).tolist())   
        writer.writerow(["NorESM1_pr_mean"] + NorESM1_pr_mean.data.flatten().astype(np.str).tolist()) 
        #writer.writerow(["Average_GCM_pr"] + Average_GCM.data.flatten().astype(np.str).tolist())
        
        writer.writerow(["CRU_pr_mean"] + CRU_pr_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["UDel_pr_mean"] + UDel_pr_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["GPCC_pr_mean"] + GPCC_pr_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["Obs_pr"] + Obs_pr.data.flatten().astype(np.str).tolist())
        
        writer.writerow(["CCCma_tas_mean"] + CCCma_tas_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["CLMcom_tas_mean"] + CLMcom_tas_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["DMI_tas_mean"] + DMI_tas_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["KNMI_tas_mean"] + KNMI_tas_mean.data.flatten().astype(np.str).tolist())   
        writer.writerow(["MPIE_tas_mean"] + MPIE_tas_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["SMHI_tas_mean"] + SMHI_tas_mean.data.flatten().astype(np.str).tolist()) 
        writer.writerow(["Average_Eraint_tas"] + Average_Eraint_tas.data.flatten().astype(np.str).tolist()) 
        
        writer.writerow(["CCCmaCanRCM_tas_mean"] + CCCmaCanRCM_tas_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["CCCmaSMHI_tas_mean"] + CCCmaSMHI_tas_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["CNRM_tas_mean"] + CNRM_tas_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["CNRMSMHI_tas_mean"] + CNRMSMHI_tas_mean.data.flatten().astype(np.str).tolist())   
        writer.writerow(["CSIRO_tas_mean"] + CSIRO_tas_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["ICHECDMI_tas_mean"] + ICHECDMI_tas_mean.data.flatten().astype(np.str).tolist())       
        writer.writerow(["ICHECCCLM_tas_mean"] + ICHECCCLM_tas_mean.data.flatten().astype(np.str).tolist()) 
        writer.writerow(["ICHECKNMI_tas_mean"] + ICHECKNMI_tas_mean.data.flatten().astype(np.str).tolist()) 
        writer.writerow(["ICHECMPI_tas_mean"] + ICHECMPI_tas_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["ICHECSMHI_tas_mean"] + ICHECSMHI_tas_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["IPSL_tas_mean"] + IPSL_tas_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["MIROC_tas_mean"] + MIROC_tas_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["MOHCCCLM_tas_mean"] + MOHCCCLM_tas_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["MOHCKNMI_tas_mean"] + MOHCKNMI_tas_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["MOHCSMHI_tas_mean"] + MOHCSMHI_tas_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["MPICCLM_tas_mean"] + MPICCLM_tas_mean.data.flatten().astype(np.str).tolist())    
        writer.writerow(["MPIREMO_tas_mean"] + MPIREMO_tas_mean.data.flatten().astype(np.str).tolist()) 
        writer.writerow(["MPISMHI_tas_mean"] + MPISMHI_tas_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["NCCDMI_tas_mean"] + NCCDMI_tas_mean.data.flatten().astype(np.str).tolist())   
        writer.writerow(["NCCSMHI_tas_mean"] + NCCSMHI_tas_mean.data.flatten().astype(np.str).tolist())    
        writer.writerow(["NOAA_tas_mean"] + NOAA_tas_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["Average_RCM_tas"] + Average_RCM_tas.data.flatten().astype(np.str).tolist())
        
        writer.writerow(["CanESM2_tas_mean"] + CanESM2_tas_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["CNRMG_tas_mean"] + CNRMG_tas_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["MK3_tas_mean"] + MK3_tas_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["EARTH_tas_mean"] + EARTH_tas_mean.data.flatten().astype(np.str).tolist())   
        writer.writerow(["EARTH3_tas_mean"] + EARTH3_tas_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["GFDL_tas_mean"] + GFDL_tas_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["HadGEM2_tas_mean"] + HadGEM2_tas_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["IPSLG_tas_mean"] + IPSLG_tas_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["MIROCG_tas_mean"] + MIROCG_tas_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["MPI_tas_mean"] + MPI_tas_mean.data.flatten().astype(np.str).tolist())   
        writer.writerow(["NorESM1_tas_mean"] + NorESM1_tas_mean.data.flatten().astype(np.str).tolist()) 
#        #writer.writerow(["Average_GCM_tas"] + Average_GCM.data.flatten().astype(np.str).tolist())
#        
        writer.writerow(["CRU_tas_mean"] + CRU_tas_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["UDel_tas_mean"] + UDel_tas_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["Obs_tas"] + Obs_tas.data.flatten().astype(np.str).tolist())
        
        writer.writerow(["CCCma_tasmin_mean"] + CCCma_tasmin_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["CLMcom_tasmin_mean"] + CLMcom_tasmin_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["DMI_tasmin_mean"] + DMI_tasmin_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["KNMI_tasmin_mean"] + KNMI_tasmin_mean.data.flatten().astype(np.str).tolist())   
        writer.writerow(["MPIE_tasmin_mean"] + MPIE_tasmin_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["SMHI_tasmin_mean"] + SMHI_tasmin_mean.data.flatten().astype(np.str).tolist()) 
        writer.writerow(["Average_Eraint_tasmin"] + Average_Eraint_tasmin.data.flatten().astype(np.str).tolist()) 
        
        writer.writerow(["CCCmaCanRCM_tasmin_mean"] + CCCmaCanRCM_tasmin_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["CCCmaSMHI_tasmin_mean"] + CCCmaSMHI_tasmin_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["CNRM_tasmin_mean"] + CNRM_tasmin_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["CNRMSMHI_tasmin_mean"] + CNRMSMHI_tasmin_mean.data.flatten().astype(np.str).tolist())   
        writer.writerow(["CSIRO_tasmin_mean"] + CSIRO_tasmin_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["ICHECDMI_tasmin_mean"] + ICHECDMI_tasmin_mean.data.flatten().astype(np.str).tolist())       
        writer.writerow(["ICHECCCLM_tasmin_mean"] + ICHECCCLM_tasmin_mean.data.flatten().astype(np.str).tolist()) 
        writer.writerow(["ICHECKNMI_tasmin_mean"] + ICHECKNMI_tasmin_mean.data.flatten().astype(np.str).tolist()) 
        writer.writerow(["ICHECMPI_tasmin_mean"] + ICHECMPI_tasmin_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["ICHECSMHI_tasmin_mean"] + ICHECSMHI_tasmin_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["IPSL_tasmin_mean"] + IPSL_tasmin_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["MIROC_tasmin_mean"] + MIROC_tasmin_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["MOHCCCLM_tasmin_mean"] + MOHCCCLM_tasmin_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["MOHCKNMI_tasmin_mean"] + MOHCKNMI_tasmin_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["MOHCSMHI_tasmin_mean"] + MOHCSMHI_tasmin_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["MPICCLM_tasmin_mean"] + MPICCLM_tasmin_mean.data.flatten().astype(np.str).tolist())    
        writer.writerow(["MPIREMO_tasmin_mean"] + MPIREMO_tasmin_mean.data.flatten().astype(np.str).tolist()) 
        writer.writerow(["MPISMHI_tasmin_mean"] + MPISMHI_tasmin_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["NCCDMI_tasmin_mean"] + NCCDMI_tasmin_mean.data.flatten().astype(np.str).tolist())   
        writer.writerow(["NCCSMHI_tasmin_mean"] + NCCSMHI_tasmin_mean.data.flatten().astype(np.str).tolist())    
        writer.writerow(["NOAA_tasmin_mean"] + NOAA_tasmin_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["Average_RCM_tasmin"] + Average_RCM_tasmin.data.flatten().astype(np.str).tolist())
        
        writer.writerow(["CanESM2_tasmin_mean"] + CanESM2_tasmin_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["CNRMG_tasmin_mean"] + CNRMG_tasmin_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["MK3_tasmin_mean"] + MK3_tasmin_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["EARTH_tasmin_mean"] + EARTH_tasmin_mean.data.flatten().astype(np.str).tolist())   
        writer.writerow(["EARTH3_tasmin_mean"] + EARTH3_tasmin_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["GFDL_tasmin_mean"] + GFDL_tasmin_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["HadGEM2_tasmin_mean"] + HadGEM2_tasmin_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["IPSLG_tasmin_mean"] + IPSLG_tasmin_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["MIROCG_tasmin_mean"] + MIROCG_tasmin_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["MPI_tasmin_mean"] + MPI_tasmin_mean.data.flatten().astype(np.str).tolist())   
        writer.writerow(["NorESM1_tasmin_mean"] + NorESM1_tasmin_mean.data.flatten().astype(np.str).tolist()) 
#        #writer.writerow(["Average_GCM_tasmin"] + Average_GCM.data.flatten().astype(np.str).tolist())
        
        writer.writerow(["CRU_tasmin_mean"] + CRU_tasmin_mean.data.flatten().astype(np.str).tolist())
       
        writer.writerow(["CCCma_tasmax_mean"] + CCCma_tasmax_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["CLMcom_tasmax_mean"] + CLMcom_tasmax_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["DMI_tasmax_mean"] + DMI_tasmax_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["KNMI_tasmax_mean"] + KNMI_tasmax_mean.data.flatten().astype(np.str).tolist())   
        writer.writerow(["MPIE_tasmax_mean"] + MPIE_tasmax_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["SMHI_tasmax_mean"] + SMHI_tasmax_mean.data.flatten().astype(np.str).tolist()) 
        writer.writerow(["Average_Eraint_tasmax"] + Average_Eraint_tasmax.data.flatten().astype(np.str).tolist()) 
        
        writer.writerow(["CCCmaCanRCM_tasmax_mean"] + CCCmaCanRCM_tasmax_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["CCCmaSMHI_tasmax_mean"] + CCCmaSMHI_tasmax_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["CNRM_tasmax_mean"] + CNRM_tasmax_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["CNRMSMHI_tasmax_mean"] + CNRMSMHI_tasmax_mean.data.flatten().astype(np.str).tolist())   
        writer.writerow(["CSIRO_tasmax_mean"] + CSIRO_tasmax_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["ICHECDMI_tasmax_mean"] + ICHECDMI_tasmax_mean.data.flatten().astype(np.str).tolist())       
        writer.writerow(["ICHECCCLM_tasmax_mean"] + ICHECCCLM_tasmax_mean.data.flatten().astype(np.str).tolist()) 
        writer.writerow(["ICHECKNMI_tasmax_mean"] + ICHECKNMI_tasmax_mean.data.flatten().astype(np.str).tolist()) 
        writer.writerow(["ICHECMPI_tasmax_mean"] + ICHECMPI_tasmax_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["ICHECSMHI_tasmax_mean"] + ICHECSMHI_tasmax_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["IPSL_tasmax_mean"] + IPSL_tasmax_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["MIROC_tasmax_mean"] + MIROC_tasmax_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["MOHCCCLM_tasmax_mean"] + MOHCCCLM_tasmax_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["MOHCKNMI_tasmax_mean"] + MOHCKNMI_tasmax_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["MOHCSMHI_tasmax_mean"] + MOHCSMHI_tasmax_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["MPICCLM_tasmax_mean"] + MPICCLM_tasmax_mean.data.flatten().astype(np.str).tolist())    
        writer.writerow(["MPIREMO_tasmax_mean"] + MPIREMO_tasmax_mean.data.flatten().astype(np.str).tolist()) 
        writer.writerow(["MPISMHI_tasmax_mean"] + MPISMHI_tasmax_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["NCCDMI_tasmax_mean"] + NCCDMI_tasmax_mean.data.flatten().astype(np.str).tolist())   
        writer.writerow(["NCCSMHI_tasmax_mean"] + NCCSMHI_tasmax_mean.data.flatten().astype(np.str).tolist())    
        writer.writerow(["NOAA_tasmax_mean"] + NOAA_tasmax_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["Average_RCM_tasmax"] + Average_RCM_tasmax.data.flatten().astype(np.str).tolist())
        
        writer.writerow(["CanESM2_tasmax_mean"] + CanESM2_tasmax_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["CNRMG_tasmax_mean"] + CNRMG_tasmax_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["MK3_tasmax_mean"] + MK3_tasmax_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["EARTH_tasmax_mean"] + EARTH_tasmax_mean.data.flatten().astype(np.str).tolist())   
        writer.writerow(["EARTH3_tasmax_mean"] + EARTH3_tasmax_mean.data.flatten().astype(np.str).tolist())      
        writer.writerow(["GFDL_tasmax_mean"] + GFDL_tasmax_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["HadGEM2_tasmax_mean"] + HadGEM2_tasmax_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["IPSLG_tasmax_mean"] + IPSLG_tasmax_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["MIROCG_tasmax_mean"] + MIROCG_tasmax_mean.data.flatten().astype(np.str).tolist())
        writer.writerow(["MPI_tasmax_mean"] + MPI_tasmax_mean.data.flatten().astype(np.str).tolist())   
        writer.writerow(["NorESM1_tasmax_mean"] + NorESM1_tasmax_mean.data.flatten().astype(np.str).tolist()) 
#        #writer.writerow(["Average_GCM_tasmax"] + Average_GCM.data.flatten().astype(np.str).tolist())
        
        writer.writerow(["CRU_tasmax_mean"] + CRU_tasmax_mean.data.flatten().astype(np.str).tolist())
       
        
if __name__ == '__main__':
    main()
    