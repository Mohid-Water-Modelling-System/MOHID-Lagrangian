#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 10:24:41 2019

@author: gfnl143
"""
import numpy as np


def getConcentrationsVolume(gridTimeInstance):
    """
    Get the number of particles in a cell divided by the volume.

    Args:
        gridTimeInstance (RectangularGridBase): Initialized RectangularGridBase.

    Returns:
        d (dict): Dictionary with dataArray format

    """
    baseName = 'concentration_volume'
    units = 'pp / m^3'
    dims = ['time', 'depth', 'latitude', 'longitude']
    data = gridTimeInstance.countsInCell/gridTimeInstance.cellVolume[np.newaxis]
    if data.shape[0] == 1:
        dims = ['time', 'latitude', 'longitude']
        data = np.squeeze(data, axis=0)

    d = {}
    d['coords'] = {k: v for k, v in gridTimeInstance.coords.items() if k in dims}
    d['dims'] = dims
    d['data'] = data
    d['attrs'] = {'units': units,
                  'long_name': baseName}

    return d


def getConcentrationsArea(gridTimeInstance):
    """
    Get the number of particles in the whole column divided by the area column.

    Args:
        gridTimeInstance (RectangularGridBase): Initialized RectangularGridBase.

    Returns:
        d (dict): Dictionary with dataArray format.

    """

    baseName = 'concentration_area'
    units = 'pp / m^2'
    dims = ['time', 'depth', 'latitude', 'longitude']
    data = gridTimeInstance.countsInCell.sum(axis=1)/gridTimeInstance.cellArea
    if data.shape[0] == 1:
        dims = ['time', 'latitude', 'longitude']
        data = np.squeeze(data, axis=0)

    d = {}
    d['coords'] = {k: v for k, v in gridTimeInstance.coords.items() if k in dims}
    d['dims'] = dims
    d['data'] = data
    d['attrs'] = {'units': units,
                  'long_name': baseName}

    return d


def getResidenceTime(gridTimeInstance, dt):
    """
    Get the time that a cell contains particles.

    Args:
        gridTimeInstance (RectangularGridBase): Initialized RectangularGridBase.

    Returns:
        d (dict): Dictionary with dataArray format

    """
    baseName = 'residence_time'
    units = 'pp / m^3'
    dims = ['time', 'depth', 'latitude', 'longitude']
    data = (gridTimeInstance.countsInCell > 0)*dt
    if data.shape[0] == 1:
        dims = ['time', 'latitude', 'longitude']
        data = np.squeeze(data, axis=0)

    d = {}
    d['coords'] = {k: v for k, v in gridTimeInstance.coords.items() if k in dims}
    d['dims'] = dims
    d['data'] = data
    d['attrs'] = {'units': units,
                  'long_name': baseName}

    return d


def getCountsInCell(gridTimeInstance):
    """
    Get the number of particles in a cell.

    Args:
        gridTimeInstance (RectangularGridBase): Initialized RectangularGridBase.

    Returns:
        d (dict): Dictionary with dataArray format

    """
    baseName = 'counts'
    units = 'ppb'
    dims = ['time', 'depth', 'latitude', 'longitude']
    data = gridTimeInstance.countsInCell
    if data.shape[0] == 1:
        dims = ['time', 'latitude', 'longitude']
        data = np.squeeze(data, axis=0)

    d = {}
    d['coords'] = {k: v for k, v in gridTimeInstance.coords.items() if k in dims}
    d['dims'] = dims
    d['data'] = data
    d['attrs'] = {'units': units,
                  'long_name': baseName}

    return d


def getVariableMeanCell(gridTimeInstance, varName, units=''):
    """
    Get the mean value of variable in a cell

    Args:
        gridTimeInstance (RectangularGridBase): Initialized RectangularGridBase.
        units (string): Units of the variable 

    Returns:
        d (dict): Dictionary with dataArray format

    """

    baseName = varName
    dims = ['time', 'depth', 'latitude', 'longitude']
    data = gridTimeInstance.meanDataInCell
    if data.shape[0] == 1:
        dims = ['time', 'latitude', 'longitude']
        data = np.squeeze(data, axis=0)

    d = {}
    d['coords'] = {k: v for k, v in gridTimeInstance.coords.items() if k in dims}
    d['dims'] = dims
    d['data'] = data
    d['attrs'] = {'units': units,
                  'long_name': baseName}

    return d
