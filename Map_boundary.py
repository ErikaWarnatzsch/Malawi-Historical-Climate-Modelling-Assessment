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
    
    #plot map with physical features 
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.COASTLINE)   
    ax.add_feature(cartopy.feature.BORDERS)
    ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
    ax.add_feature(cartopy.feature.RIVERS)
    #set map boundary
    ax.set_extent([32., 36., -8, -17]) 
    #set axis tick marks
    ax.set_xticks([33, 34, 35]) 
    ax.set_yticks([-10, -12, -14, -16]) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #save the image of the graph and include full legend
    plt.savefig('Map_data_boundary', bbox_inches='tight')
    plt.show()
    
    
    
if __name__ == '__main__':
    main()