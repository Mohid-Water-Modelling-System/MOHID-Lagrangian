#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 15:06:29 2019

@author: gfnl143
"""
import xarray as xr
import numpy as np
from netCDF4 import Dataset

def getDimsAttrs(dimensionName,dimensionData=None):

    if dimensionName == 'longitude':
        attrs = {'long_name': 'longitude',\
        'standard_name': 'longitude',\
        'units': 'degrees_east',
        '_FillValue': -9.8999995E15,
        'valid_min': -180.0,\
        'valid_max': 180.0}
        if dimensionData != None:
            attrs['valid_min']= np.min(dimensionData)
            attrs['valid_max']: np.max(dimensionData)
    elif dimensionName == 'latitude':
        attrs = {'long_name': 'latitude',\
        'standard_name': 'latitude',\
        'units': 'degrees_north',\
        '_FillValue': -9.8999995E15,\
        'valid_min': -90.0,\
        'valid_max': 90.0}
        if dimensionData != None:
            attrs['valid_min']= np.min(dimensionData)
            attrs['valid_max']: np.max(dimensionData)
    elif dimensionName == 'depth':
        attrs = {'long_name': 'depth',\
        'standard_name': 'depth',\
        'units': 'meters',\
        '_FillValue': -9.8999995E15}
        if dimensionData != None:
            attrs['valid_min']= np.min(dimensionData)
            attrs['valid_max']= np.max(dimensionData)
    elif dimensionName == 'time':
        attrs = {'long_name':'time',
         'units':'days since 1950-01-01 00:00:00'}

    return attrs


def getVarsAttrs(variableName):
    if variableName =='concentration':
        attrs = {'long_name':'concentration','units':'pp/m^2'}
    elif variableName == 'residence_time':
        attrs = {'long_name':'residence_time', 'units':'s'}
    elif variableName == 'age':
        attrs = {'long_name':'age','units':'s'}
    elif variableName == 'velocity':
        attrs = {'long_name':'velocity_magnitude','units':'m/s'}
    elif variableName == 'state':
        attrs = {'long_name':'beached'}
    else:
        attrs = {}
    
    return attrs


class NetcdfParser:
    

    def __init__(self,fileName):
        self.dataset = []
        self.fileName = fileName
        
    def initDataset(self,spatialGrid,timeGrid):
        coords = {'time':('time',timeGrid.timeAxis.round(decimals=10)),
                  spatialGrid.dims[0]: (spatialGrid.dims[0],spatialGrid.cellCenters[0]),
                  spatialGrid.dims[1]: (spatialGrid.dims[1],spatialGrid.cellCenters[1]),
                  spatialGrid.dims[2]: (spatialGrid.dims[2],spatialGrid.cellCenters[2]),
                  }
        self.dataset = xr.Dataset(None,coords=coords)
        
        for dimensionName in self.dataset.dims:
            self.dataset[dimensionName].attrs = getDimsAttrs(dimensionName)
        
        self.dataset.to_netcdf(self.fileName)
        print('-> Dataset initizalized in: ',self.fileName)
        return

    def appendVariableTimeStepToDataset(self,variableNetcdfName,dataArray,step):
        #print('--> Writing '+ variableNetcdfName)
        ds = Dataset(self.fileName,'a')
        if step == 0:
            formatData = str(dataArray['data'].dtype.kind) + str(dataArray['data'].dtype.alignment)
            _ = ds.createVariable(variableNetcdfName,formatData,dataArray['dims'])
        appendvar = ds.variables[variableNetcdfName]
        appendvar[step] = dataArray['data']
        ds.close()


""" Xarray implementation -- deprecated
#   The Xarray do not let to append data to a existing dataset. This formulation
#   was replaced by native netcdf Dataset package                  
#    def appendVariableToDataset(self,variableNetcdfName,dataArray):
#        print('--> Writing '+ variableNetcdfName) #' for source: ',self.sources['id'][str(source)])
#        ds = xr.open_dataset(self.fileName)
#        self.dataset[variableNetcdfName] = dataArray
#        ds.close()
#        ds.to_netcdf(self.fileName,'a')
#        return
#    
#    def checkDataset(self):
#        with xr.open_dataset(self.fileName) as ds:
#            if ds.depth.size == 1:
#                ds_sq = ds.load()
#                ds.close()
#                print('->The nc has a depth degenerated dimension: Squeezing...' )
#                ds_sq = ds_sq.squeeze(dim='depth',drop=True)
#                # When you alter the dimensions of the netcdf, you cannot overwrite the
#                # netcdf. We ha create a new one without extension, remove the old one,
#                # and rename the first one.
#                #nc_squeezed = self.netcdf_output_file.replace('.nc','_sq.nc')
#                ds_sq.to_netcdf(self.fileName)
"""