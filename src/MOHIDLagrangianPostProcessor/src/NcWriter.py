# -*- coding: utf-8 -*-
""" Module to manage the store of the results in netcdf format. """

import xarray as xr
import numpy as np
from netCDF4 import Dataset
from src.Grid import Grid
from src.Polygon import Polygon

class NetcdfParser:

    def __init__(self, fileName):
        self.fileName = fileName

    def initDataset(self, spatialGrid, timeGrid):

        coords = {**timeGrid.coords, **spatialGrid.coords}

        dataset = xr.Dataset(None, coords=coords)

        for dimensionName in dataset.dims:
            dataset[dimensionName].attrs = NetcdfParser.getDimsAttrs(dimensionName)

        dataset.to_netcdf(self.fileName)

        print('-> Dataset initizalized in: ', self.fileName)

    def appendVariableTimeStepToDataset(self, variableName, dataArray, step):
        ds = Dataset(self.fileName, 'a')
        if step == 0:
            formatData = str(dataArray['data'].dtype.kind) + str(dataArray['data'].dtype.alignment)
            ncvar = ds.createVariable(variableName, formatData, dataArray['dims'])
            ncvar.units = dataArray['attrs']['units']
        appendvar = ds.variables[variableName]
        appendvar[step] = dataArray['data']
        ds.close()

    @staticmethod
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

    @staticmethod
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

