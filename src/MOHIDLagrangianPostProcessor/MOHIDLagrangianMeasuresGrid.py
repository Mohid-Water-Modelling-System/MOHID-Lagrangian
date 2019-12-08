#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 10:24:41 2019

@author: gfnl143
"""
import numpy as np


def getConcentrationsVolume(gridTimeInstance):
    baseName = 'concentration_volume'
    units =  'pp / m^3'
    dims = ['time','depth','latitude','longitude']
    data = gridTimeInstance.countsInCell/gridTimeInstance.cellVolume[np.newaxis]
    if data.shape[0] == 1:
        dims = ['time','latitude','longitude']
        data = np.squeeze(data,axis=0)

    d = {}    
    d['coords'] = {k:v for k,v in gridTimeInstance.coords.items() if k in dims}
    d['dims'] = dims
    d['data'] = data
    d['attrs'] = {'units':units,'long_name':baseName}
    
    return d
    

def getConcentrationsArea(gridTimeInstance):
    baseName = 'concentration_area'
    units =  'pp / m^2'
    dims = ['time','depth','latitude','longitude']
    data = gridTimeInstance.countsInCell.sum(axis=1)/gridTimeInstance.cellArea
    if data.shape[0] == 1:
       dims = ['time','latitude','longitude']
       data = np.squeeze(data,axis=0)

    d = {}
    d['coords'] = {k:v for k,v in gridTimeInstance.coords.items() if k in dims}
    d['dims'] = dims
    d['data'] = data
    d['attrs'] = {'units':units,'long_name':baseName}
    
    return d


def getResidenceTime(gridTimeInstance,dt):
    baseName = 'residence_time'
    units =  'pp / m^3'
    dims = ['time','depth','latitude','longitude']
    data = (gridTimeInstance.countsInCell > 0)*dt
    if data.shape[0] == 1:
        dims = ['time','latitude','longitude']
        data = np.squeeze(data,axis=0)
    
    d = {}
    d['coords'] = {k:v for k,v in gridTimeInstance.coords.items() if k in dims}
    d['dims'] = dims
    #gridTimeInstance.countsInCell[gridTimeInstance.countsInCell > 0] = dt
    d['data'] = data
    d['attrs'] = {'units':units,'long_name':baseName}
    
    return d

    
def getCountsInCell(gridTimeInstance):
   baseName = 'counts'
   units =  'ppb'
   dims = ['time','depth','latitude','longitude']
   data = gridTimeInstance.countsInCell
   if data.shape[0] == 1:
       dims = ['time','latitude','longitude']
       data = np.squeeze(data,axis=0)
   
   d = {}
   d['coords'] = {k:v for k,v in gridTimeInstance.coords.items() if k in dims}
   d['dims'] = dims
   d['data'] = data
   d['attrs'] = {'units':units,'long_name':baseName}

   return d


def getVariableMeanCell(gridTimeInstance,varName,units=''):
    baseName = varName
    units =  units
    dims = ['time','depth','latitude','longitude']
    data = gridTimeInstance.meanDataInCell
    if data.shape[0] == 1:
       dims = ['time','latitude','longitude']
       data = np.squeeze(data,axis=0)
        
    d = {}
    d['coords'] = {k:v for k,v in gridTimeInstance.coords.items() if k in dims}
    d['dims'] = dims
    d['data'] = data
    d['attrs'] = {'units':units,'long_name':baseName}
    

    return d
     