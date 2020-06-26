#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 15:48:27 2020

@author: gfnl143
"""


import xarray as xr
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cmx

import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from src.PlotFeatures import *

from src.XMLReader import getPlotTypeFromRecipe, getPlotMeasuresFromRecipe
from src.XMLReader import getPlotWeightFromRecipe, getPlotTimeFromRecipe


def plot_it(dataArray, output_filename, title, units, plot_type='contourf',):
    """


    Args:
        dataArray (TYPE): DESCRIPTION.
        output_filename (TYPE): DESCRIPTION.
        plot_type (TYPE, optional): DESCRIPTION. Defaults to 'contourf'.

    Returns:
        None.

    """

    time_plot_flag = hastime(dataArray)

    if time_plot_flag:
        time_key = dataArray.dims[0]
        time_size = dataArray.shape[0]
        if time_size < 4:
            nrows = 1
            ncols = time_size
        elif time_size == 4:
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
    else:
        nrows = ncols = 1

    if nrows == ncols:
        figsize = (17, 15)
    elif ncols > nrows:
        figsize = (20, 15)

    fig, axarr = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize,
                              subplot_kw={'projection': ccrs.PlateCarree()})

    # generalization -> make one subplot iterable.
    if not isinstance(axarr, np.ndarray):
        axarr = np.array([axarr])

    extent = get_extent(dataArray)
    vmin, vmax = get_color_lims(dataArray, robust=True)
    scale_bar_lenght = get_horizontal_scale(dataArray)
    cmap_key = get_cmap_key(vmin, vmax)
    cmap_norm = get_cmap_norm(vmin,vmax)
    time_step = 0
    for ax in axarr.flat:

        if time_plot_flag:
            dataArray_step = dataArray.isel({time_key: time_step})
        else:
            dataArray_step = dataArray

        try:
            setup_plot = {'x':'longitude',
                          'y':'latitude',
                          'cmap': cmap_key,
                          'ax': ax,
                          'vmin': vmin,
                          'vmax': vmax,
                          'zorder': 1,
                          'add_colorbar': False}

            p = getattr(dataArray_step.plot, plot_type)(**setup_plot)

            ax = get_background_map(ax, extent)
            gl = ax.gridlines(draw_labels=True, color='gray', linestyle='--')
            gl.xlabels_top = gl.ylabels_right = False
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER

            if np.abs((extent[3]-extent[2]) > ((extent[1]-extent[0]))):
                gl.xlabel_style = {'rotation': 45}

            scale_bar(ax, ccrs.PlateCarree(), scale_bar_lenght)
            time_step += 1

        except Exception as exc:
              print('-> WARNING: Could not plot', title)
              print('-> ERROR:',exc)
              return

    # Creating a common colorbar.
    cbar_x, cbar_y, size_x, size_y = get_cbar_position(axarr) # get extent
    cbar_ax = fig.add_axes([cbar_x, cbar_y, size_x, size_y])
    scalarMap = cmx.ScalarMappable(norm=cmap_norm, cmap=cmap_key)
    cbar = fig.colorbar(scalarMap, cax=cbar_ax)
    cbar.set_label(units)
    # Creating the title from the filename
    fig.suptitle(title, fontsize='x-large')
    #fig.tight_layout()
    # Save the title
    fig.savefig(output_filename, dpi=150)
    plt.close()



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
    print('-> Grouping time steps:', group_freq, group_type)
    print('-> Methods:', methods)
    print('-> Plot_type', plot_type)

    ds = xr.open_mfdataset(outDir+'*.nc')
    variables = list(ds.keys())

    # sources = map(get_source_from_variable, variables)
    # materials = map(get_material_from_variable, variables)

    # # remove duplicates elements
    # sources = list(dict.fromkeys(sources))
    # materials = list(dict.fromkeys(materials))

    # print(list(sources))
    # print(list(materials))



    for idxvar, variable in enumerate(variables):
        da = ds[variable].load()
        units = da.units
        
        if weight_file:
            if idxvar == 0:
                print("->Weighting sources with:", weight_file[0] )
                print("   %-30s | %9s" % ('Source', 'Weight'))
            da = weight_dataarray_with_csv(da, weight_file[0])

        if 'depth' in da.dims:
            da = da.isel(depth=-1)

        methods_list = []
        for method in methods:
            if group_freq != 'all':
                if isgroupable(da, group_type, group_freq, method):
                    da = group_resample(da, group_type, group_freq, method)

            da = getattr(da, method)(dim='time')
            methods_list.append(method)


        da = da.where(da != 0)

        output_filename = outDir + '-'.join(methods_list) + '-' + variable + '.png'
        title = get_title_methods(methods_list, variable)

        plot_it(da, output_filename, title, units, plot_type[0])
