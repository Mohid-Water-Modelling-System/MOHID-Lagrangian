#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 15:48:27 2020

@author: gfnl143
"""

import pandas as pd
import xarray as xr

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap

from src.PlotFeatures import scale_bar, get_extent, get_color_lims
from src.PlotFeatures import get_title_methods, get_horizontal_scale

import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from src.XMLReader import getPlotTypeFromRecipe, getPlotMeasuresFromRecipe
from src.XMLReader import getPlotWeightFromRecipe, getPlotTimeFromRecipe


def toPlotTimeSteps(dataArray, output_filename, title,
                    plot_type='contourf', robust=True):
    """


    Args:
        dataArray (TYPE): DESCRIPTION.
        output_filename (TYPE): DESCRIPTION.
        plot_type (TYPE, optional): DESCRIPTION. Defaults to 'contourf'.

    Returns:
        None.

    """
    time_key = dataArray.dims[0]
    time_size = dataArray.shape[0]

    if time_size == 4:
        nrows = 2
        ncols = 2
    elif time_size == 12:
        nrows = 3
        ncols = 4
    else:
        time_slice = slice(0, -1, int(time_size/4))
        dataArray = dataArray.isel({time_key: time_slice})
        nrows = 2
        ncols = 2

    fig, axarr = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 10),
                              subplot_kw={'projection': ccrs.PlateCarree()})

    #request = cimgt.Stamen('terrain-background')
    #ax.add_image(request, 2)

    extent = get_extent(dataArray)
    vmin, vmax = get_color_lims(dataArray, robust=True)
    scale_bar_lenght = get_horizontal_scale(dataArray)


    land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                            edgecolor='face',
                                            facecolor=np.array((0.75, 0.75, 0.75)))

    time_step = 0
    for ax in axarr.flat:
        dataArray_step = dataArray[time_step,:,:] # BUG with.isel({time_key:time_step})
        getattr(dataArray_step.plot, plot_type)(x='longitude', y='latitude',
                                                cmap=get_cmap("rainbow"),
                                                ax=ax,
                                                cbar_kwargs={'shrink': 0.6},
                                                vmin=vmin,
                                                vmax=vmax
                                                )

        ax.set_extent(extent)
        ax.add_feature(land_10m)
        ax.add_feature(cfeature.OCEAN, color='white')
        gl = ax.gridlines(draw_labels=True, color='gray', linestyle='--')
        gl.xlabels_top = gl.ylabels_right = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        scale_bar(ax, ccrs.PlateCarree(), scale_bar_lenght)
        time_step += 1

    # Creating the title from the filename
    # If you do a cumsum(mean), the title should look in the same way.
    fig.suptitle(title, fontsize = 'x-large')

    fig.savefig(output_filename, dpi=300)
    plt.close()


def toPlotTimeStep(dataArray, output_filename, title,
                   plot_type='contourf', robust=True):
    """


    Args:
        dataArray (TYPE): DESCRIPTION.
        output_filename (TYPE): DESCRIPTION.
        plot_type (TYPE, optional): DESCRIPTION. Defaults to 'contourf'.

    Returns:
        None.

    """
    fig = plt.figure(figsize=(15, 10))
    #request = cimgt.Stamen('terrain-background')
    #ax.add_image(request, 2)

    ax = plt.subplot(projection=ccrs.PlateCarree())

    extent = get_extent(dataArray)
    vmin, vmax = get_color_lims(dataArray, robust=True)

    ax.set_extent(extent)

    land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                            edgecolor='face',
                                            facecolor=np.array((0.75, 0.75, 0.75)))

    getattr(dataArray.plot, plot_type)(x='longitude', y='latitude',
                                       ax=ax,
                                       cmap=get_cmap("rainbow"),
                                       zorder=10,
                                       cbar_kwargs={'shrink': 0.6},
                                       vmin=vmin,
                                       vmax=vmax
                                       )

    ax.add_feature(land_10m)
    ax.add_feature(cfeature.OCEAN, color='white')
    gl = ax.gridlines(draw_labels=True, color='gray', linestyle='--')
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    scale_bar(ax, ccrs.PlateCarree(), 10)

    # Creating the title from the filename
    # If you do a cumsum(mean), the title should look in the same way.
    fig.suptitle(title, fontsize='x-large')

    fig.savefig(output_filename, dpi=300)
    plt.close()


def toPlot(dataArray, output_filename, title, plot_type='contourf'):
    """


    Args:
        dataArray (TYPE): DESCRIPTION.
        output_filename (TYPE): DESCRIPTION.
        plot_type (TYPE, optional): DESCRIPTION. Defaults to 'contourf'.

    Returns:
        None.

    """

    if (len(dataArray.dims) > 2) and (dataArray.shape[0] > 1):
        toPlotTimeSteps(dataArray, output_filename, title, plot_type)

    else:
        dataArray = dataArray[0,:,:]
        toPlotTimeStep(dataArray, output_filename, title, plot_type)


def GroupResample(da, time_group, time_freq, measure):
    """


    Args:
        ds (TYPE): DESCRIPTION.
        time_group (TYPE): DESCRIPTION.
        measure (TYPE): DESCRIPTION.

    Returns:
        flag (TYPE): DESCRIPTION.

    """
    if time_group == 'resample':
         da = da.resample(time=time_freq)
    elif time_group == 'groupby':
         da = da.groupby(time_freq)
    return da

def isGroupable(da, time_group, time_freq, measure):
    """


    Args:
        ds (TYPE): DESCRIPTION.
        time_group (TYPE): DESCRIPTION.
        measure (TYPE): DESCRIPTION.

    Returns:
        flag (TYPE): DESCRIPTION.

    """
    if time_group == 'resample':
        flag = (measure in dir(da.resample(time=time_freq)))
    elif time_group == 'groupby':
        flag = (measure in dir(da.groupby(time_freq)))
    return flag


def weightDataset(dataset, weight_file):
    """


    Args:
        dataset (TYPE): DESCRIPTION.
        weight_file (TYPE): DESCRIPTION.

    Returns:
        dataset (TYPE): DESCRIPTION.

    """
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
    group_freq, group_type = getPlotTimeFromRecipe(xml_recipe)
    methods = getPlotMeasuresFromRecipe(xml_recipe)
    plot_type = getPlotTypeFromRecipe(xml_recipe)
    weight_file = getPlotWeightFromRecipe(xml_recipe)

    print('-> Plotting results:')
    print('-> Grouping time steps:',group_freq, group_type)
    print('-> Methods:', methods)
    print('-> Plot_type', plot_type)

    ds = xr.open_mfdataset(outDir+'*.nc')
    variables = list(ds.keys())

    if weight_file:
        ds = weightDataset(ds, weight_file[0])

    for variable in variables:
        da = ds[variable].load()
        if 'depth' in da.dims:
            da = da.isel(depth=-1)

        methods_list = []
        for method in methods:
            if group_freq != 'all':
                if isGroupable(da, group_type, group_freq, method):
                    da = GroupResample(da, group_type, group_freq, method)

            da = getattr(da, method)(dim='time')
            methods_list.append(method)

        if ('contour' in plot_type[0]) == False:
            da = da.where(da != 0)

        output_filename = outDir + '-'.join(methods_list) + '-' + variable + '.png'
        title = get_title_methods(methods_list, variable)

        toPlot(da, output_filename, title, plot_type[0])
