# -*- coding: utf-8 -*-
"""
Created on Wednesday August 2 2017

@author: s0899345
"""

import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import numpy as np
import iris
import cartopy
import iris.plot as iplt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter


def main():
    
    #bring in altitude data
    Elev = '/exports/csce/datastore/geos/users/s0899345/Climate_Modelling/Actual_Data/elev.0.25-deg.nc'
    Elev= iris.load_cube(Elev)
    
    #remove variable for time
    Elev = iris.util.squeeze(Elev)
    
    #define colours for contour map    
    cmap = mpl_cm.get_cmap('YlGn') 
    
    #set levels and extent
    levels = np.arange(0,2000,150)
    extend = 'max'
            
    #plot data 
    Elev = iplt.contourf(Elev, cmap=cmap, levels=levels, extend=extend) 
    ax = plt.gca()
    
    #add colour bar index and a label
    plt.colorbar(Elev, label='Meters above sea-level')
    
    #plot map with physical features 
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
    plt.savefig('Map_Elev', bbox_inches='tight')
    plt.show() 
 
if __name__ == '__main__':
    main()