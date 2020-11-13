# -*- coding: utf-8 -*-
"""
Plot results module.
"""

from tqdm import tqdm
import cartopy.crs as ccrs
import xarray as xr
import numpy as np
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from src.PlotFeatures import *
from src.XMLReader import *
from .Utils import *
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cmx


class Plot:
    """ Parent class to plot a dataArray. """

    def __init__(self, dataArray, output_filename, title, units):
        self.dataArray = dataArray
        self.fig = []
        self.axarr = []
        self.cbar = []
        self.extent = []
        self.time_key = []
        self.time_size = self.dataArray.shape[0]
        self.setup_plot = {}
        self.output_filename = output_filename
        self.title = title
        self.units = units
        self.plot_type = []
        self.n_plot_max_threshold = 12
        self.plot_every_n_step = 4

    def setTimeKey(self):
        """ set the time key from dataArray if exist"""
        if hastime(self.dataArray):
            self.time_key = get_time_key(self.dataArray)
        else:
            self.time_key = None

    def setSliceTimeDataArray(self):
        """ Slice to dataArray to fit into a graph. """

        if hastime(self.dataArray):
            # If time_size is higher than 12 (for example, 30 days, it plots
            # one of each four days.
            if self.time_size > self.n_plot_max_threshold: # Timesteps > 12 slices 
                time_slice = slice(0, -1, int(self.time_size/self.plot_every_n_step))
                self.dataArray = self.dataArray.isel({self.time_key: time_slice})

    def setFigureAxisLayout(self):
        """
        Creates the figure layout to plot.

        Returns:
            None.

        """
        time_plot_flag = hastime(self.dataArray)
        if time_plot_flag:
            if self.time_size < 4:
                nrows = 1
                ncols = self.time_size
            elif self.time_size == 4:
                nrows = 2
                ncols = 2
            elif self.time_size == 6:
                nrows = 2
                ncols = 3
            elif self.time_size == 8:
                nrows = 2
                ncols = 4
            elif self.time_size == 12:
                nrows = 3
                ncols = 4
            else:
                nrows = 2
                ncols = 2
        else:
            nrows = ncols = 1

        if nrows == ncols:
            figsize = (17, 15)
        elif ncols > nrows:
            figsize = (20, 15)

        self.fig, self.axarr = plt.subplots(nrows=nrows, ncols=ncols,
                                            figsize=figsize,
                                            subplot_kw={'projection': ccrs.PlateCarree()})

        if not isinstance(self.axarr, np.ndarray):
            self.axarr = np.array([self.axarr])

    def setColorbar(self):
        vmin, vmax = get_color_lims(self.dataArray, robust=True)
        cmap_key = get_cmap_key(vmin, vmax)
        cmap_norm = get_cmap_norm(vmin, vmax)
        cbar_x, cbar_y, size_x, size_y = get_cbar_position(self.axarr) # get extent
        cbar_ax = self.fig.add_axes([cbar_x, cbar_y, size_x, size_y])
        scalarMap = cmx.ScalarMappable(norm=cmap_norm, cmap=cmap_key)
        self.cbar = self.fig.colorbar(scalarMap, cax=cbar_ax)
        self.cbar.set_label(self.units)

    def getBackground(self, ax, extent):
        """Get the background map, gridlines and ticks for the current axis."""
        ax = get_background_map(ax, extent)
        gl = ax.gridlines(draw_labels=True, color='gray', linestyle='--')
        gl.xlabels_top = gl.ylabels_right = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

        if np.abs((extent[3]-extent[2]) > ((extent[1]-extent[0]))):
            gl.xlabel_style = {'rotation': 45}
        return ax

    def getScaleBar(self, ax):
        """Set the scale bar in km for the current axis"""
        scale_bar_lenght = get_horizontal_scale(self.dataArray)
        scale_bar(ax, ccrs.PlateCarree(), scale_bar_lenght)


class PlotPolygon(Plot):
    """ Module to plot polygon dataArrays"""

    def __init__(self, dataArray, output_filename, title, units, polygon_file):
        Plot.__init__(self, dataArray, output_filename, title, units)
        self.polygon_file = polygon_file

    def getSetupDict(self, ax):
        vmin, vmax = get_color_lims(self.dataArray, robust=True)
        cmap_key = get_cmap_key(vmin, vmax)
        cmap_norm = get_cmap_norm(vmin, vmax)

        setup_plot = {
              'cmap': cmap_key,
              'ax': ax,
              'vmin': vmin,
              'vmax': vmax,
              'zorder': 1}

        return setup_plot

    def getPlots(self):
        geoDataFrame = gpd.read_file(self.polygon_file).to_crs({"init": 'EPSG:4326'})
        geoDataFrame['index'] = np.arange(0, geoDataFrame.shape[0])
        self.setTimeKey()
        self.setFigureAxisLayout()
        self.setSliceTimeDataArray()
        self.setColorbar()

        time_step = 0
        for ax in self.axarr.flat:

            extent = get_extent(geoDataFrame)
            ax = self.getBackground(ax, extent)

            if hastime(self.dataArray):
                dataArray_step = self.dataArray.isel({self.time_key: time_step})
            else:
                dataArray_step = self.dataArray

            varName = self.dataArray.name
            geoDataFrame[varName] = dataArray_step.values

            extent = get_extent(geoDataFrame)
            ax = self.getBackground(ax, extent)
            setupPlotDict = self.getSetupDict(ax)

            geoDataFrame.plot(column=varName, **setupPlotDict)
            time_step += 1

        # Creating the suptitle from the filename
        self.fig.suptitle(self.title, fontsize='x-large')
        # Save the image
        self.fig.savefig(self.output_filename, dpi=150)
        # Save the shapefile
        shapefile = self.output_filename.split('.')[0]+'.shp'
        geoDataFrame.to_file(shapefile)
        plt.close()


class PlotGrid(Plot):
    """ Module to plot grid dataArrays"""
    def __init__(self, dataArray, output_filename, title, units, plot_type):
        Plot.__init__(self, dataArray, output_filename, title, units)
        self.plot_type = plot_type

    def getSetupDict(self, ax):
        vmin, vmax = get_color_lims(self.dataArray, robust=True)
        cmap_key = get_cmap_key(vmin, vmax)
        cmap_norm = get_cmap_norm(vmin, vmax)

        if self.plot_type == 'hist':
            setup_plot = {'ax': ax, 'zorder': 1}

        else:
            setup_plot = {'x':'longitude',
                          'y':'latitude',
                          'ax': ax,
                          'cmap': cmap_key,
                          'vmin': vmin,
                          'vmax': vmax,
                          'zorder': 1,
                          'add_colorbar': False}
        return setup_plot

    def getPlots(self):
        self.setTimeKey()
        self.setFigureAxisLayout()
        self.setSliceTimeDataArray
        self.setColorbar()

        time_step = 0
        for ax in self.axarr.flat:

            extent = get_extent(self.dataArray)
            if hastime(self.dataArray):
                dataArray_step = self.dataArray.isel({self.time_key: time_step})
            else:
                dataArray_step = self.dataArray

            setup_plot = self.getSetupDict(ax)
            _ = getattr(dataArray_step.plot, self.plot_type)(**setup_plot)

            self.getBackground(ax, extent)
            self.getScaleBar(ax)
            time_step += 1

        # Creating the title from the filename
        self.fig.suptitle(self.title, fontsize='x-large')
        # fig.tight_layout()
        # Save the title
        self.fig.savefig(self.output_filename, dpi=150)
        plt.close()



def plotResultsFromRecipe(outDir, xml_recipe):
    """


    Args:
        outDir (TYPE): DESCRIPTION.
        xml_recipe (TYPE): DESCRIPTION.

    Returns:
        None.

    """
    # Read the xml attributes.
    group_freq, group_type = getPlotTimeFromRecipe(xml_recipe)
    groups = getGroupFromRecipe(xml_recipe)
    methods = getPlotMeasuresFromRecipe(xml_recipe)
    plot_type = getPlotTypeFromRecipe(xml_recipe)
    weight_file = getPlotWeightFromRecipe(xml_recipe)
    normalize_method = getNormalizeFromRecipe(xml_recipe)
    polygon_file = getPolygonFileFromRecipe(xml_recipe)

    print('-> Plotting results:')
    print('-> Grouping time steps:', group_freq, group_type)
    print('-> Methods:', methods)
    print('-> Plot_type:', plot_type)
    print('-> Normalizing:', normalize_method)

    # Read the dataset and the variables keys
    ds = xr.open_mfdataset(outDir + '*.nc')
    variables = list(ds.keys())

    # Select just the desire group of variables
    if groups:
        variables = groupByWords(variables, groups)
        if len(variables) == 0:
            print('-> You did an empty group. Stopping')
            raise ValueError
        ds = ds[variables]

    # Normalization
    if  normalize_method:
        normalizer = Normalizer(normalize_method, ds)
        normalizer.setFactor()
        normalizer.setMethodName()

    for idxvar, variable in enumerate(tqdm(variables)):
        da = ds[variable].load()
        units = da.units

        if weight_file:
            if idxvar == 0:
                print("-> Weighting sources with:", weight_file[0])
                print("   %-40s | %9s" % ('Source', 'Weight'))
            da = weight_dataarray_with_csv(da, weight_file[0])

        if 'depth' in da.dims:
            da = da.isel(depth=-1)

        methods_list = []
        for method in methods:
            if group_freq != 'all':
                if isgroupable(da, group_type, group_freq, method):
                    da = group_resample(da, group_type, group_freq)

            da = getattr(da, method)(dim='time')
            methods_list.append(method)

        if normalize_method:
            units = normalizer.getNormalizedUnits(units)
            da = normalizer.getNormalizedDataArray(da)

        da = da.where(da != 0)  # Values with 0 consider as missing value.

        output_filename = outDir + '-'.join(methods_list) + '-' + variable + '.png'
        title = get_title_methods(methods_list, variable)

        if polygon_file:
            plotter = PlotPolygon(da, output_filename, title, units, polygon_file[0])
            plotter.getPlots()
        else:
            plotter = PlotGrid(da, output_filename, title, units, plot_type[0])
            plotter.getPlots()
