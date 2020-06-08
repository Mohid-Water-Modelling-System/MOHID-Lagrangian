#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 15:48:27 2020

@author: gfnl143
"""

from src.XMLReader import getPlotTypeFromRecipe, getPlotMeasuresFromRecipe, getPlotTimeFromRecipe

import matplotlib.pyplot as plt
from math import floor
from matplotlib.cm import get_cmap
from matplotlib import patheffects
import xarray as xr
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


def utm_from_lon(lon):
    """
    utm_from_lon - UTM zone for a longitude

    Not right for some polar regions (Norway, Svalbard, Antartica)

    :param float lon: longitude
    :return: UTM zone number
    :rtype: int
    """
    return floor((lon + 180) / 6) + 1


def scale_bar(ax, proj, length, location=(0.5, 0.05), linewidth=3,
              units='km', m_per_unit=1000):
    """

    http://stackoverflow.com/a/35705477/1072212
    ax is the axes to draw the scalebar on.
    proj is the projection the axes are in
    location is center of the scalebar in axis coordinates ie.
    0.5 is the middle of the plot
    length is the length of the scalebar in km.
    linewidth is the thickness of the scalebar.
    units is the name of the unit
    m_per_unit is the number of meters in a unit
    """
    # find lat/lon center to find best UTM zone
    x0, x1, y0, y1 = ax.get_extent(proj.as_geodetic())
    # Projection in metres
    utm = ccrs.UTM(utm_from_lon((x0+x1)/2))
    # Get the extent of the plotted area in coordinates in metres
    x0, x1, y0, y1 = ax.get_extent(utm)
    # Turn the specified scalebar location into coordinates in metres
    sbcx, sbcy = x0 + (x1 - x0) * location[0], y0 + (y1 - y0) * location[1]
    # Generate the x coordinate for the ends of the scalebar
    bar_xs = [sbcx - length * m_per_unit/2, sbcx + length * m_per_unit/2]
    # buffer for scalebar
    buffer = [patheffects.withStroke(linewidth=5, foreground="w")]
    # Plot the scalebar with buffer
    ax.plot(bar_xs, [sbcy, sbcy], transform=utm, color='k',
            linewidth=linewidth, path_effects=buffer)
    # buffer for text
    buffer = [patheffects.withStroke(linewidth=3, foreground="w")]
    # Plot the scalebar label
    t0 = ax.text(sbcx, sbcy, str(length) + ' ' + units, transform=utm,
                 horizontalalignment='center', verticalalignment='bottom',
                 path_effects=buffer, zorder=2)

    left = x0+(x1-x0)*0.05
    # Plot the N arrow
    t1 = ax.text(left, sbcy, u'\u25B2\nN', transform=utm,
                 horizontalalignment='center', verticalalignment='bottom',
                 path_effects=buffer, zorder=2)
    # Plot the scalebar without buffer, in case covered by text buffer
    ax.plot(bar_xs, [sbcy, sbcy], transform=utm, color='k',
            linewidth=linewidth, zorder=3)


def toPlot(dataArray, output_filename, plot_type='contourf'):
    fig = plt.figure(figsize=(15, 10))
    request = cimgt.Stamen('terrain-background')
    ax = plt.subplot(projection=ccrs.PlateCarree())

    extent = [dataArray.longitude.min(),
              dataArray.longitude.max(),
              dataArray.latitude.min(),
              dataArray.latitude.max()]

    ax.set_extent(extent)

    getattr(dataArray.plot, plot_type)('longitude', 'latitude',
                                       ax=ax,
                                       cmap=get_cmap("jet"),
                                       zorder=10
                                       )

    ax.add_image(request, 11)

    gl = ax.gridlines(draw_labels=True, color='black')
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    scale_bar(ax, ccrs.PlateCarree(), 10)

    fig.savefig(output_filename, dpi=300)
    plt.close()


def plotResultsFromRecipe(outDir, xml_recipe):
    time_group = getPlotTimeFromRecipe(xml_recipe)
    methods = getPlotMeasuresFromRecipe(xml_recipe)
    plot_type = getPlotTypeFromRecipe(xml_recipe)

    print('-> Plotting results:')
    print('-> Grouping time steps:', time_group)
    print('-> Methods:', methods)
    print('-> Plot_type', plot_type)

    if 'all' in time_group:
        ds = xr.open_mfdataset(outDir+'*.nc')
    # else:
    #     ds = xr.open_mfdataset(outDir+'*.nc').groupby(time_group)

    if ds['depth'].size == 1:
        ds = ds.drop_dims('depth')
    else:
        ds = ds.isel(depth=-1)

    variables = list(ds.keys())

    for variable in variables:
        _da = ds[variable]
        if 'depth' in ds:
            _da = ds[variable].isel(depth=-1)
        else:
            _da = ds[variable]

        _method = ''
        for method in methods:
            _da = getattr(_da, method)(dim='time')
            _method = _method + method
        _da = _da.where(_da != 0)
        toPlot(_da, _method + '_' + variable + '.png', plot_type[0])
