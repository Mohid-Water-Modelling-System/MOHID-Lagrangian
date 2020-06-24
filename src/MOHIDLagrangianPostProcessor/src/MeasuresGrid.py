#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 10:24:41 2019

@author: gfnl143
"""
import numpy as np


def to_dict(dict_coords, dims, data, units, long_name):
    d = {}
    d['coords'] = {k: v for k, v in dict_coords if k in dims}
    d['dims'] = dims
    d['data'] = data
    d['attrs'] = {'units': units,
                  'long_name': long_name}
    return d

def is2D(array):
    return array.ndim == 2

def is2Dlayer(array):
    return (array.ndim == 3) and (array.shape[0] == 1)


def getConcentrationsVolume(gridTimeInstance):
    """
    Get the number of particles in a cell divided by the volume.

    Args:
        gridTimeInstance (RectangularGridBase): Initialized RectangularGridBase.

    Returns:
        d (dict): Dictionary with dataArray format

    """
    long_name = 'concentration_volume'
    units = 'particles/km^3'
    dims = ['time', 'depth', 'latitude', 'longitude']
    data = gridTimeInstance.countsInCell/gridTimeInstance.cellVolume
    dict_coords = gridTimeInstance.coords.items()

    # If depth has one layer, squeeze depth 'layer'
    if is2Dlayer(data):
        dims = ['time', 'latitude', 'longitude']
        data = np.squeeze(data, axis=0)

    # turn into dataArray-ready dict
    da_dict = to_dict(dict_coords, dims, data, units, long_name)

    return da_dict


def getConcentrationsArea(gridTimeInstance):
    """
    Get the number of particles in the whole column divided by the area column.

    Args:
        gridTimeInstance (RectangularGridBase): Initialized RectangularGridBase.

    Returns:
        d (dict): Dictionary with dataArray format.

    """

    long_name = 'concentration_area'
    units = 'particles/km^2'
    dims = ['time', 'depth', 'latitude', 'longitude']
    dict_coords = gridTimeInstance.coords.items()
    data = gridTimeInstance.countsInCell.sum(axis=0)/gridTimeInstance.cellArea

    if is2Dlayer(data):
        dims = ['time', 'latitude', 'longitude']
        data = np.squeeze(data, axis=0)

    # turn into dataArray-ready dict
    da_dict = to_dict(dict_coords, dims, data, units, long_name)

    return da_dict


def getResidenceTime(gridTimeInstance, dt):
    """
    Get the time that a cell contains particles.

    Args:
        gridTimeInstance (RectangularGridBase): Initialized RectangularGridBase.

    Returns:
        d (dict): Dictionary with dataArray format

    """
    long_name = 'residence_time'
    units = 's'
    dims = ['time', 'depth', 'latitude', 'longitude']
    data = (gridTimeInstance.countsInCell > 0)*dt
    dict_coords = gridTimeInstance.coords.items()

    if is2Dlayer(data):
        dims = ['time', 'latitude', 'longitude']
        data = np.squeeze(data, axis=0)

    # turn into dataArray-ready dict
    da_dict = to_dict(dict_coords, dims, data, units, long_name)

    return da_dict


def getCountsInCell(gridTimeInstance):
    """
    Get the number of particles in a cell.

    Args:
        gridTimeInstance (RectangularGridBase): Initialized RectangularGridBase.

    Returns:
        d (dict): Dictionary with dataArray format

    """
    long_name = 'counts'
    units = 'particles-per-box'
    dims = ['time', 'depth', 'latitude', 'longitude']
    data = gridTimeInstance.countsInCell
    dict_coords = gridTimeInstance.coords.items()

    if is2Dlayer(data):
        dims = ['time', 'latitude', 'longitude']
        data = np.squeeze(data, axis=0)

    # turn into dataArray-ready dict
    da_dict = to_dict(dict_coords, dims, data, units, long_name)

    return da_dict


def getVariableMeanCell(gridTimeInstance, varName, units='n.u'):
    """
    Get the mean value of variable in a cell

    Args:
        gridTimeInstance (RectangularGridBase): Initialized RectangularGridBase.
        units (string): Units of the variable

    Returns:
        d (dict): Dictionary with dataArray format

    """

    long_name = varName
    dims = ['time', 'depth', 'latitude', 'longitude']
    data = gridTimeInstance.meanDataInCell
    dict_coords = gridTimeInstance.coords.items()

    if is2Dlayer(data):
        dims = ['time', 'latitude', 'longitude']
        data = np.squeeze(data, axis=0)

    # turn into dataArray-ready dict
    da_dict = to_dict(dict_coords, dims, data, units, long_name)

    return da_dict
