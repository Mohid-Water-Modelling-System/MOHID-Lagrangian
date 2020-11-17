# -*- coding: utf-8 -*-
"""
Module with functions related to plotting stage.
"""

import numpy as np
import pandas as pd
import xarray as xr
from math import floor
from matplotlib import patheffects
import matplotlib.colors as mcolors

import geopandas as gpd
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


def get_time_key(da: xr.DataArray) -> str:
    """
    Get the time "key" from a dataArray.

    Args:
        da (xr.DataArray): Input Datarray.

    Returns:
        str: DESCRIPTION.

    """
    dim_name = ['time', 'year', 'week', 'month', 'season', 'day',
                'hour', 'minute', 'second']
    da_dim_names = list(da.dims)
    time_key = [name for name in dim_name if name in da_dim_names]
    return time_key[0]


def group_resample(da: xr.DataArray, group_type: str,
                   time_freq: str) -> xr.DataArray:
    """
    Resample/group the DataArray based on the time frequency.

    Args:
        da (xr.DataArray): Input dataArray
        group_type (str): 'resample/group' method to group data.
        time_freq (TYPE): Frequency to apply one or other group method.

    'resample': Any pandas resample method: '3M', '2Y',...
    'group': Any xarray group method: 'time.week','time.season','time.'

     Returns:
        da (TYPE): grouped or resampled dataArray

    Example:
        Consider a timeseries of two years with dayly data.

        Resample: split the dataset into time intervals from freq.
        Ex: '1M', will split the timeseries
        split the timeseries into months.
        Group: group dataset time steps based on time_freq. Ex. 'time.month'.
        Will split the series into 12 months, and group all Januarys



    """

    if group_type == 'resample':
        da = da.resample(time=time_freq)
    elif group_type == 'groupby':
        da = da.groupby(time_freq)
    return da


def isgroupable(da: xr.DataArray, group_type: str, time_freq: str,
                method: str) -> bool:
    """
    Check if the dataset is groupable using the provided measure

    Args:
        da (xr.DataArray): Input dataArray
        group_type (str): 'resample/group' method to group data.
        time_freq (srtE): Frequency to apply one or other group method.

    Returns:
        flag (bool): True if the measure can be applied.

    """
    flag = False
    if group_type == 'resample':
        flag = (method in dir(da.resample(time=time_freq)))
    elif group_type == 'groupby':
        flag = (method in dir(da.groupby(time_freq)))
    return flag


def weight_dataset_with_csv(dataset: xr.Dataset, weight_file: str) -> xr.Dataset:
    """
    Multiply the dataset using the source weight provided in the csv file.

    Args:
        dataset (xr.Dataset): Input dataset to weight.
        weight_file (str): File with each weight per source

    Returns:
        dataset (xr.Dataset): Weighted dataset

    """

    df = pd.read_csv(weight_file)
    df = df.set_index('source')

    for datasetVar in list(dataset.keys()):
        for indexVar in df.index:
            if indexVar in datasetVar.names:
                weight = df.loc[indexVar, 'weight']
                dataset[datasetVar] = dataset[datasetVar]*df.loc[indexVar, 'weight']
                print("-> %-40s | %4.1f" % (datasetVar, weight))

    return dataset


def weight_dataarray_with_csv(dataArray: xr.DataArray,
                              weight_file: str) -> xr.DataArray:
    """
    Multiply the dataArray within the source weight provided in the csv file.

    Args:
        dataArray (xr.DataArray): Input dataArray to weight.
        weight_file (str): File with each weight per source

    Returns:
        dataset (xr.Dataset): Weighted dataset

    """

    df = pd.read_csv(weight_file)
    df = df.set_index('source')

    for indexVar in df.index:
        if indexVar in dataArray.name:
            weight = df.loc[indexVar, 'weight']
            dataArray = dataArray*df.loc[indexVar, 'weight']
            print("-> %-40s | %4.1f" % (dataArray.name, weight))

    return dataArray


def get_cmap_key(vmin: float, vmax: float) -> str:
    """
    Get the cmap key based on limits.

    If the boundary limits are positive and negative, returns a divergent cmap.

    Args:
        vmin (float): Data minimum
        vmax (float): Daata Maximum

    Returns:
        cmap_key (str): Returns a str with the color key code.

    """
    if vmin*vmax > 0:
        cmap_key = 'rainbow'
    else:
        cmap_key = 'RdBu_r'
    return cmap_key


def get_cmap_norm(vmin: float, vmax: float) -> mcolors:
    """
    Get the cmap norm  key based on limits.

    If the boundary limits are positive and negative, returns a  0 divergent
    norm.

    Args:
        vmin (float): Data minimum
        vmax (float): Daata Maximum


    Returns:
        mcolors: norm.

    """
    if vmin*vmax >= 0:
        cNorm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    else:
        color_max = np.max((np.abs(vmax), np.abs(vmin)))
        cNorm = mcolors.TwoSlopeNorm(vmin=-color_max,
                                     vmax=color_max,
                                     vcenter=0)
    return cNorm


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
        norm = mcolors.BoundaryNorm(boundaries=boundaries_RGB_tiles, ncolors=64)
        ax.add_image(stamen_terrain, 11, cmap='gray_r', norm=norm)
    return ax


def get_color_lims(dataArray: xr.DataArray, robust: bool = True,
                   min_quartile=0.05, max_quartile=0.99):
    """
    Get the vmax and vmin from the dataArray

    Args:
        dataArray (xr.DataArray): DESCRIPTION.
        robust (bool, optional): True, slice data between min and max quartile.
        extremes. Defaults to True.
        min_quartile (float, optional): DESCRIPTION. Defaults to 0.05.
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


def get_grid_extent(dataArray: xr.DataArray) -> list:
    """
    Get the dataArray domains limits [x_min, x_max, y_min, y_max]

    Args:
        dataArray (xr.DataArray): Input dataArray.

    Returns:
        list: extent with [x_min, x_max, y_min, y_max]

    """
    extent = [dataArray.longitude.min(),
              dataArray.longitude.max(),
              dataArray.latitude.min(),
              dataArray.latitude.max()]
    return extent


def get_polygon_extent(geoDataframe: gpd.GeoDataFrame) -> list:
    """
    Get the GeoDataframe domains limits x_min, x_max, y_min, y_max]
    Args:
        dataArray (gpd.GeoDataFrame): Input geodataframe.

    Returns:
        list:  extent with [x_min, x_max, y_min, y_max]

    """
    extent = geoDataframe.total_bounds
    return [extent[0], extent[2], extent[1], extent[3]]


def get_extent(inputdata) -> list:
    """
    Get extent from Datarray or GeoDataframe  [x_min, x_max, y_min, y_max]
    """
    if isinstance(inputdata, xr.DataArray):
        return get_grid_extent(inputdata)
    elif isinstance(inputdata, gpd.GeoDataFrame):
        return get_polygon_extent(inputdata)


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


def get_material_from_variable(variable: str):
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

    material_names = ['light', 'heavy', 'medium']

    for post_measure in post_measures:
        if post_measure in variable:
            source_name = variable.replace(post_measure, '')

    for material_name in material_names:
        if material_name in source_name:
            material_name = source_name.replace(material_name, '')

    return material_name


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



class Normalizer:
    
    def __init__(self, method, dataset):
        self.method = method
        self.ds = dataset
        self.factor = []
        self.normalized_name = []
        self.setFactor()
        self.setMethodName()
        
    @staticmethod
    def getTotalAmmountEmitted(dataset):
        if 'n_counts_global' in dataset:
            n = dataset['n_counts_global'].diff('time').sum().values
            return n
    
    @staticmethod
    def getMeanStdGlobal(dataset, dim=None):
        if dim is None:
            return dataset['n_counts_global'].mean(), dataset['n_counts_global'].std()
        else:
            return dataset['n_counts_global'].mean(dim=dim), dataset['n_counts_global'].std(dim=dim)

    @staticmethod
    def getMaxGlobal(dataset):
        return dataset['n_counts_global'].max()

    def setFactor(self):
        if self.method == 'total':
            self.factor = Normalizer.getTotalAmmountEmitted(self.ds)
        if self.method == 'mean':
            self.factor= Normalizer.getMeanStdGlobal(self.ds)
        if self.method == 'mean-zonal':
            self.factor = Normalizer.getMeanStdGlobal(self.ds,dim='time')
        if self.method == 'max':
            self.factor = Normalizer.getMaxGlobal(self.ds)

    def setMethodName(self):
        if self.method == 'total':
            self.normalized_name = '% over total emmited'
        if self.method == 'mean':
            self.normalized_name = 'normalized over climate mean'
        if self.method == 'mean-zonal':
            self.normalized_name = 'normalized over climate zonal-mean'
        if self.method == 'max':
            self.normalize_name = 'normalized over maximum'
    
    def getNormalizedUnits(self, units):
        if self.method == 'total':
            units = units + self.normalized_name
        elif (self.method == 'mean') or (self.method == 'mean-zonal'):
            units = units + self.normalized_name
        elif self.method == 'max':
            units = units + self.normalized_name
        return units
            
    def getNormalizedDataArray(self, da):
        if self.method == 'total':
            # units are place before.
            # Datarray lost attributes whern dataArray arithmetic broadcast.
            # is done
            da.values = (da.values/self.factor)*100.
        elif (self.method == 'mean') or (self.method == 'mean-zonal'):
            da = da/self.factor[0]
        elif self.method == 'max':
            da.values = da.values/self.factor

        return da
