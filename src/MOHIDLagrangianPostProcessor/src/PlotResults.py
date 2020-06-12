#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 15:48:27 2020

@author: gfnl143
"""

import pandas as pd
import xarray as xr
from math import floor

import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib import patheffects

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from src.XMLReader import getPlotTypeFromRecipe, getPlotMeasuresFromRecipe
from src.XMLReader import getPlotWeightFromRecipe, getPlotTimeFromRecipe


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


def toPlotTimeSteps(dataArray, output_filename, plot_type='contourf'):
    """


    Args:
        dataArray (TYPE): DESCRIPTION.
        output_filename (TYPE): DESCRIPTION.
        plot_type (TYPE, optional): DESCRIPTION. Defaults to 'contourf'.

    Returns:
        None.

    """

    time_slice = slice(0, -1, int(dataArray.time.size/4))
    fig, axarr = plt.subplots(nrows=2, ncols=2, figsize=(15, 10),
                              subplot_kw={'projection': ccrs.PlateCarree()})

    request = cimgt.Stamen('terrain-background')

    extent = [dataArray.longitude.min(),
              dataArray.longitude.max(),
              dataArray.latitude.min(),
              dataArray.latitude.max()]

    time_step = 0
    for ax in axarr.flat:
        dataArray_step = dataArray.isel(time=time_step)
        getattr(dataArray_step.plot, plot_type)('longitude', 'latitude',
                                                cmap=get_cmap("jet"),
                                                robust=True,
                                                ax=ax,
                                                projection=ccrs.PlateCarree(),
                                                cbar_kwargs={'shrink': 0.6})

        ax.set_extent(extent)
        ax.add_image(request, 11)
        gl = ax.gridlines(draw_labels=True, color='black')
        gl.xlabels_top = gl.ylabels_right = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        scale_bar(ax, ccrs.PlateCarree(), 10)

        time_step += 1

    # Creating the title from the filename
    # If you do a cumsum(mean), the title should look in the same way.
    title = getTitleFromFilename(output_filename)
    fig.suptitle(title, fontsize = 'x-large')

    fig.savefig(output_filename, dpi=300)
    plt.close()


def toPlotTimeStep(dataArray, output_filename, plot_type='contourf'):
    """


    Args:
        dataArray (TYPE): DESCRIPTION.
        output_filename (TYPE): DESCRIPTION.
        plot_type (TYPE, optional): DESCRIPTION. Defaults to 'contourf'.

    Returns:
        None.

    """
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
                                       zorder=10,
                                       robust=True
                                       )

    ax.add_image(request, 11)

    gl = ax.gridlines(draw_labels=True, color='black')
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    scale_bar(ax, ccrs.PlateCarree(), 10)

    # Creating the title from the filename
    # If you do a cumsum(mean), the title should look in the same way.
    title = getTitleFromFilename(output_filename)
    fig.suptitle(title, fontsize='x-large')

    fig.savefig(output_filename, dpi=300)
    plt.close()


def toPlot(dataArray, output_filename, plot_type='contourf'):
    """
    

    Args:
        dataArray (TYPE): DESCRIPTION.
        output_filename (TYPE): DESCRIPTION.
        plot_type (TYPE, optional): DESCRIPTION. Defaults to 'contourf'.

    Returns:
        None.

    """
    if dataArray.time.size == 1:
        toPlotTimeStep(dataArray, output_filename, plot_type)
    else:
        toPlotTimeSteps(dataArray, output_filename, plot_type)


def getTitleFromFilename(filename, debug_title=False):
    """
    

    Args:
        filename (TYPE): DESCRIPTION.
        debug_title (TYPE, optional): DESCRIPTION. Defaults to False.

    Returns:
        TYPE: DESCRIPTION.

    """

    methods_str = filename.split('.')[0].split('-')[:-1]  # Turn the results in operation order.
    var_str = filename.split('.')[0].split('-')[-1]       # Get the var name
    methods_funcstyle = [methods_str[0]] + ['('+methods_str[k] for k in range(1, len(methods_str))]
    title_funcstyle = methods_funcstyle + ['(']+[var_str] + ['(t)'] + [')']*(len(methods_str)-1)

    if debug_title is True:
        print('-> method_str', methods_str)
        print('-> title_funcsyle: ', title_funcstyle)
        print('-> var_str: ', var_str)
        print('-> methods_funcstyle: ', methods_funcstyle)

    return ''.join(title_funcstyle)


def isResamplable(ds, time_group, measure):
    """


    Args:
        ds (TYPE): DESCRIPTION.
        time_group (TYPE): DESCRIPTION.
        measure (TYPE): DESCRIPTION.

    Returns:
        flag (TYPE): DESCRIPTION.

    """
    flag = measure in dir(ds.resample(time=time_group))
    return flag


def weightDataset(dataset, weight_file):
    df = pd.read_csv(weight_file)
    df = df.set_index('name')

    for datasetVar in list(dataset.keys()):
        for indexVar in df.index:
            if indexVar in datasetVar:
                weight = df.loc[indexVar, 'weight']
                dataset[datasetVar] = dataset[datasetVar]*df.loc[indexVar, 'weight']
                print('-> Weighting', datasetVar, ' by ', weight)

    return dataset


def plotResultsFromRecipe(outDir, xml_recipe):
    """


    Args:
        outDir (TYPE): DESCRIPTION.
        xml_recipe (TYPE): DESCRIPTION.

    Returns:
        None.

    """
    # they returns a list.
    time_group = getPlotTimeFromRecipe(xml_recipe)[0]
    methods = getPlotMeasuresFromRecipe(xml_recipe)
    plot_type = getPlotTypeFromRecipe(xml_recipe)[0]
    weight_file = getPlotWeightFromRecipe(xml_recipe)[0]

    print('-> Plotting results:')
    print('-> Grouping time steps:', time_group)
    print('-> Methods:', methods)
    print('-> Plot_type', plot_type)

    ds = xr.open_mfdataset(outDir+'*.nc')
    variables = list(ds.keys())

    if weight_file:
        ds = weightDataset(ds, weight_file)

    if ds['depth'].size == 1:
        ds = ds.drop_dims('depth')
    else:
        ds = ds.isel(depth=-1)

#    if time_group != 'all':
#        ds = xr.open_mfdataset(outDir+'*.nc').resample(time=time_group[0])

    for variable in variables:
        _da = ds[variable]
        if 'depth' in ds:
            _da = ds[variable].isel(depth=-1)
        else:
            _da = ds[variable]

        # if time_group != 'all':
        #     _da = _da.resample(time=time_group[0])

        _method = ''
        for method in methods:
            if time_group != 'all':
                if isResamplable(_da, time_group, method):
                    _da = _da.resample(time=time_group)
            _da = getattr(_da, method)(dim='time')
            _method = method + '-' + _method
        _da = _da.where(_da != 0)
        toPlot(_da, _method + variable + '.png', plot_type)
