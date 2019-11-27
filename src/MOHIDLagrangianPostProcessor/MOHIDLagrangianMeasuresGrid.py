#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 10:24:41 2019

@author: gfnl143
"""
import numpy as np
import xarray as xr


def getConcentrationsVolume(gridTimeInstance):
    baseName = 'concentration_volume'
    units =  'pp / m^3'
    dims = ['time','depth','latitude','longitude']

    d = {}
    d['coords'] = {k:v for k,v in gridTimeInstance.coords.items() if k in dims}
    d['dims'] = dims
    d['data'] = gridTimeInstance.countsInCell/gridTimeInstance.cellVolume[np.newaxis]
    d['attrs'] = {'units':units,'long_name':baseName}
    
    return d
    

def getConcentrationsArea(gridTimeInstance):
    baseName = 'concentration_area'
    units =  'pp / m^2'
    dims = ['time','depth','latitude','longitude']
    
    d = {}
    d['coords'] = {k:v for k,v in gridTimeInstance.coords.items() if k in dims}
    d['dims'] = dims
    d['data'] = gridTimeInstance.countsInCell.sum(axis=1)/gridTimeInstance.cellArea
    d['attrs'] = {'units':units,'long_name':baseName}
    
    return d


def getResidenceTime(gridTimeInstance,dt):
    baseName = 'residence_time'
    units =  'pp / m^3'
    dims = ['time','depth','latitude','longitude']
    
    d = {}
    d['coords'] = {k:v for k,v in gridTimeInstance.coords.items() if k in dims}
    d['dims'] = dims
    #gridTimeInstance.countsInCell[gridTimeInstance.countsInCell > 0] = dt
    d['data'] = (gridTimeInstance.countsInCell > 0)*dt
    d['attrs'] = {'units':units,'long_name':baseName}
    
    return d

    
def getCountsInCell(gridTimeInstance):
   baseName = 'counts'
   units =  'ppb'
   dims = ['time','depth','latitude','longitude']
   
   d = {}
   d['coords'] = {k:v for k,v in gridTimeInstance.coords.items() if k in dims}
   d['dims'] = dims
   d['data'] = gridTimeInstance.countsInCell
   d['attrs'] = {'units':units,'long_name':baseName}

   return d


def getVariableMeanCell(gridTimeInstance,varName,units=''):
    baseName = varName
    units =  units
    dims = ['time','depth','latitude','longitude']
    
    d = {}
    d['coords'] = {k:v for k,v in gridTimeInstance.coords.items() if k in dims}
    d['dims'] = dims
    d['data'] = gridTimeInstance.meanDataInCell
    d['attrs'] = {'units':units,'long_name':baseName}
    

    return d
     