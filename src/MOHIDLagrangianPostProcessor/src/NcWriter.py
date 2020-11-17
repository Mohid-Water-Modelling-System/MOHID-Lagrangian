# -*- coding: utf-8 -*-
""" Module to manage and store the results in netcdf format. """

import xarray as xr
import numpy as np
from netCDF4 import Dataset


class NetcdfParser:
    
    """Class to manage and store the results in netcdf format. """

    def __init__(self, fileName):
        self.fileName = fileName
        self.timeIdx = 0

    def initDataset(self, spatialGrid, timeGrid):
        """ Initialize an empyt netcdf dataset using the Grid/time
        coordinates.

        Args:
            spatialGrid (Grid): Grid instance.
            timeGrid (Time): Time instance

        Returns:
            None.

        """

        coords = {**timeGrid.coords, **spatialGrid.coords}
        dataset = xr.Dataset(None, coords=coords)

        # Add attributes to dimensions (CF-compliant).
        for dimensionName in dataset.dims:
            dataset[dimensionName].attrs = NetcdfParser.getDimsAttrs(dimensionName)

        dataset.to_netcdf(self.fileName)
        print('-> Dataset initizalized in: ', self.fileName)

    def appendVariableTimeStepToDataset(self, variableName: str,
                                        dataArray: xr.DataArray):
        """
        Append variable data on time dimension

        Args:
            variableName (str): Name of the variable to append.
            dataArray (xr.DataArray): variable data to append
            step (int): time dimension index where to append.

        Returns:
            None.

        """

        ds = Dataset(self.fileName, 'a')
        # At first step -> initilize headers for each variable.
        if self.timeIdx == 0:
            # get the data type: (f4,f8, i4, i8...)
            data_kind = str(dataArray['data'].dtype.kind)
            data_alingment = str(dataArray['data'].dtype.alignment)
            formatData = data_kind + data_alingment
            ncvar = ds.createVariable(variableName, formatData, dataArray['dims'])
            ncvar.units = dataArray['attrs']['units']
        appendvar = ds.variables[variableName]
        appendvar[self.timeIdx] = dataArray['data']
        ds.close()

    def increaseTimeIdx(self):
        """Increase time index."""
        self.timeIdx += 1

    def resetTimeIdx(self):
        """Reset time index."""
        self.timeIdx = 0

    @staticmethod
    def getDimsAttrs(dimensionName: str, dimensionData=None) -> dict:
        """
         Get dimension attributtes (CF-Compliant).

        Args:
            dimensionName (str): Name of the dimension to inquire attrs
            dimensionData (np.array/xr.ataArray, optional): Data dimension
            to find min and max values.Defaults to None.

        Returns:
            dict: DESCRIPTION.

        """

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
    def getVarsAttrs(variableName: str) -> dict:
        """
            Gets the variable attributtes (CF-Compliant).

            Args:
                variableName (str): Variable name to inquire attrs
        """
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
