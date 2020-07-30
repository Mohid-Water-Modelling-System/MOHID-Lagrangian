#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 15:06:29 2019

@author: gfnl143
"""
import xarray as xr
import numpy as np
from netCDF4 import Dataset
from src.Grid import Grid
from src.Polygon import Polygon


def getDimsAttrs(dimensionName, dimensionData=None):
    if dimensionName == 'longitude':
        attrs = {'long_name': 'longitude',
                 'standard_name': 'longitude',
                 'units': 'degrees_east',
                 '_FillValue': -9.8999995E15,
                 'valid_min': -180.0,
                 'valid_max': 180.0}
        if dimensionData is not None:
            attrs = {'valid_min': np.min(dimensionData),
                     'valid_max': np.max(dimensionData)
                     }
    elif dimensionName == 'latitude':
        attrs = {'long_name': 'latitude',
                 'standard_name': 'latitude',
                 'units': 'degrees_north',
                 '_FillValue': -9.8999995E15,
                 'valid_min': -90.0,
                 'valid_max': 90.0}
        if dimensionData is not None:
            attrs = {'valid_min': np.min(dimensionData),
                     'valid_max': np.max(dimensionData)
                     }
    elif dimensionName == 'depth':
        attrs = {'long_name': 'depth',
                 'standard_name': 'depth',
                 'units': 'meters',
                 '_FillValue': -9.8999995E15}
        if dimensionData is not None:
            attrs = {'valid_min': np.min(dimensionData),
                     'valid_max': np.max(dimensionData)
                     }
    elif dimensionName == 'time':
        attrs = {'long_name': 'time',
                 'units': 'days since 1950-01-01 00:00:00'}

    elif dimensionName == 'index':
        attrs = {'long_name': 'polygon_identified',
                 'units': 'id'}

    return attrs


def getVarsAttrs(variableName):
    if variableName == 'concentration_area':
        attrs = {'long_name': 'concentration',
                 'units': 'p/km^2'}
    elif variableName == 'concentration_volume':
        attrs = {'long_name': 'concentration',
                 'units': 'p/km^3'}
    elif variableName == 'residence_time':
        attrs = {'long_name': 'residence_time',
                 'units': 's'}
    elif variableName == 'age':
        attrs = {'long_name': 'age',
                 'units': 's'}
    elif variableName == 'velocity':
        attrs = {'long_name': 'velocity_magnitude',
                 'units': 'm/s'}
    elif variableName == 'state':
        attrs = {'long_name': 'beached'}
    else:
        attrs = {}

    return attrs


class NetcdfParser:

    def __init__(self, fileName):
        self.fileName = fileName

    def initPolygonDataset(self, polygon, timeGrid):
        coords = {'time': ('time', timeGrid.timeAxis.round(decimals=10)),
                  polygon.dim: (polygon.dim, polygon.ids)}

        dataset = xr.Dataset(None, coords=coords)

        for dimensionName in dataset.dims:
            dataset[dimensionName].attrs = getDimsAttrs(dimensionName)

        dataset.to_netcdf(self.fileName)
        print('-> Dataset initizalized in: ', self.fileName)


    def initGridDataset(self, spatialGrid, timeGrid):
        coords = {'time': ('time', timeGrid.timeAxis.round(decimals=10)),
                  spatialGrid.dims[0]: (spatialGrid.dims[0], spatialGrid.cellCenters[0]), 
                  spatialGrid.dims[1]: (spatialGrid.dims[1], spatialGrid.cellCenters[1]), 
                  spatialGrid.dims[2]: (spatialGrid.dims[2], spatialGrid.cellCenters[2]), 
                  }
        dataset = xr.Dataset(None, coords=coords)

        for dimensionName in dataset.dims:
            dataset[dimensionName].attrs = getDimsAttrs(dimensionName)

        dataset.to_netcdf(self.fileName)
        print('-> Dataset initizalized in: ', self.fileName)


    def initDataset(self, spatialData, timeData):
        if isinstance(spatialData, Polygon):
            self.initPolygonDataset(spatialData, timeData)
        elif isinstance(spatialData, Grid):
            self.initGridDataset(spatialData, timeData)

    def appendVariableTimeStepToDataset(self, variableNetcdfName, dataArray, step):
        ds = Dataset(self.fileName, 'a')
        if step == 0:
            formatData = str(dataArray['data'].dtype.kind) + str(dataArray['data'].dtype.alignment)
            ncvar = ds.createVariable(variableNetcdfName, formatData, dataArray['dims'])
            ncvar.units = dataArray['attrs']['units']
        appendvar = ds.variables[variableNetcdfName]
        appendvar[step] = dataArray['data']
        ds.close()
