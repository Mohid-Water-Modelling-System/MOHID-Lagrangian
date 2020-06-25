#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 15:28:44 2020

@author: gfnl143
"""

import numpy as np
import pandas as pd
import xarray as xr
from math import floor
from matplotlib import patheffects
import matplotlib.cm, matplotlib.colors

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt


def hastime(dataArray: xr.DataArray) -> bool:
    """
    Check if dataArray has time or not

    Alternative future implementation based on key search:

    dims = len(dataArray.shape)
    time_dim_size = dataArray.shape[0]
    if (dims >= 3) and (time_dim_size >= 1):
        dim_name = ['time','year','week','month','season','day',
                    'hour','minute','second']
        flag = any([name in dataArray for name in dim_name])

    Args:
        dataArray (xr.DataArray): DESCRIPTION.


    Returns:
        bool: DESCRIPTION.

    """
    flag = False
    dims = len(dataArray.shape)
    time_dim_size = dataArray.shape[0]
    if (dims >= 3) and (time_dim_size >= 1):
        flag = True
    else:
        flag = False

    return flag


def group_resample(da, time_group, time_freq, measure):
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


def isgroupable(da, time_group, time_freq, measure):
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


def weight_dataset(dataset, weight_file):
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


def get_background_map(ax, extent):
    """
    Get the background map based on the width of extent in degrees.

    Args:
        ax (cartopy.mpl.geoaxes.GeoAxes): DESCRIPTION.
        extent (list): DESCRIPTION.

    Returns:
        ax (TYPE): DESCRIPTION.

    """
    gray_color_RGB = np.array((0.75, 0.75, 0.75))
    extent_width = np.abs(extent[1]-extent[0])
    if extent_width > 2:  # one degree ~ 111 km
        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                                edgecolor='face',
                                                facecolor=gray_color_RGB)
        ax.add_feature(land_10m)
        ax.add_feature(cfeature.OCEAN, color='white')

    else:
        stamen_terrain = cimgt.Stamen('toner-background')
        stamen_terrain.desired_tile_form = 'L'
        # tile - color limist RGB [0, 256]
        # 0 -  5 water - white
        # 5 - 64 land - white gray
        boundaries_RGB_tiles = [0, 5, 64]
        norm = matplotlib.colors.BoundaryNorm(boundaries=boundaries_RGB_tiles, ncolors=64)
        ax.add_image(stamen_terrain, 11, cmap='gray_r', norm=norm)
    return ax


def get_color_lims(dataArray: xr.DataArray, robust: bool = True,
                   min_quartile=0.05, max_quartile=0.99):
    """
    Get the vmax and vmin from the dataArray

    Args:
        dataArray (xr.DataArray): DESCRIPTION.
        robust (bool, optional): DESCRIPTION. Defaults to True.
        min_quartile (float, optional): DESCRIPTION. Defaults to 0.01.
        max_quartile (float, optional): DESCRIPTION. Defaults to 0.99.

    Returns:
        vmin (float): DESCRIPTION.
        vmax (float): DESCRIPTION.

    """
    if robust is True:
        vmin = dataArray.quantile(min_quartile).values
        vmax = dataArray.quantile(max_quartile).values
    else:
        vmin = dataArray.min().values
        vmax = dataArray.max().values
    return vmin, vmax


def get_extent(dataArray):
    extent = [dataArray.longitude.min(),
              dataArray.longitude.max(),
              dataArray.latitude.min(),
              dataArray.latitude.max()]
    return extent


def get_horizontal_scale(dataArray: xr.DataArray) -> np.int:
    """Get the horizontal scale to build the scale bar.

    It approaches 1º-100km to tenth part of horizontal extent

    Args:
        dataArray (xr.DataArray): DESCRIPTION.

    Returns:
        int: Scale in km to set the bar.

    """
    _max = dataArray.longitude.max().values
    _min = dataArray.longitude.min().values
    dlon = np.abs(_max - _min)

    extent_horizontal_fraction = 0.1
    degrees_to_km_approach = 100

    scale = np.round(dlon*extent_horizontal_fraction)*degrees_to_km_approach
    # If scale < 1 - 10 km.
    if scale < 1:
        scale = 10
    return scale


def get_title_methods(methods: list, variable: str,
                      debug_title: bool = True) -> str:
    """

    Args:
        methods (list): list with strings of methods used
        variable (str): str with the variable where the methods applied.
        debug_title (bool, optional): Debug flag with info. Defaults to True.

    Returns:
        str: Title.

    """
    source_name = get_source_from_variable(variable)
    measure_name = get_measure_from_variable(variable)

    methods_variable_list = methods[::-1] + [measure_name] + ['t']
    start_parenthesis = '('.join(methods_variable_list)
    end_parenthesis = ''.join((len(methods_variable_list)-1)*[')'])

    title = r"$\bf{Method :}$" + start_parenthesis + end_parenthesis + '\n' + \
        r"$\bf{Source :}$" + source_name

    return title


def get_source_from_variable(variable: str):
    """
    get the source name from variable

    Args:
        variable (str): variable name from postprocessor

    Returns:
        source_name (str): source name

    """

    post_measures = ['concentration_volume_',
                     'concentration_area_',
                     'n_counts_',
                     'residence_time_']

    for post_measure in post_measures:
        if post_measure in variable:
            source_name = variable.replace(post_measure, '')

    return source_name


def get_cbar_position(axarr: list):
    """
    get the source name from variable

    Args:
        axarr (list): list with axis

    Returns:
        source_name (str): source name

    """
    # limits = [x, y, size_x, size_y]
    axarr_limits = np.zeros((axarr.flatten().size, 4))
    for i, ax in enumerate(axarr.flat):
        axarr_limits[i, :] = ax.get_position().bounds

    axarr_limits[:, 2] = axarr_limits[:, 0] + axarr_limits[:, 2]
    axarr_limits[:, 3] = axarr_limits[:, 1] + axarr_limits[:, 3]

    # array with min_x, min_y, max_x, max_y
    min_x = np.min(axarr_limits[:, 0])
    min_y = np.min(axarr_limits[:, 1])
    max_x = np.max(axarr_limits[:, 2])
    max_y = np.max(axarr_limits[:, 3])

    cbar_x = max_x + (1 - max_x)*0.35
    cbar_y = min_y + (max_y - min_y)*0.15

    size_x = 0.015
    size_y = (max_y-min_y)*0.65
    return cbar_x, cbar_y, size_x, size_y


def get_measure_from_variable(variable: str):
    """
    get the measure from variable

    Args:
        variable (str): DESCRIPTION.

    Returns:
        post_measure (TYPE): DESCRIPTION.

    """

    post_measures = ['concentration_volume_', 'concentration_area_',
                     'n_counts_', 'residence_time_']

    for post_measure in post_measures:
        if post_measure in variable:
            return post_measure


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

    left = x0+(x1-x0)*0.1
    up = y0+(y1-y0)*0.9
    # Plot the N arrow
    t1 = ax.text(left, up, u'\u25B2\nN', transform=utm,
                 horizontalalignment='center', verticalalignment='bottom',
                 path_effects=buffer, zorder=2)
    # Plot the scalebar without buffer, in case covered by text buffer
    ax.plot(bar_xs, [sbcy, sbcy], transform=utm, color='k',
            linewidth=linewidth, zorder=3)