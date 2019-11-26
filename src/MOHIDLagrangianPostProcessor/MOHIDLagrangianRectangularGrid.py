#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 14:44:32 2019

@author: gfnl143
"""
import numpy as np
import xml.etree.ElementTree as ET
import constants as cte
from numba import jit

@jit(nopython=True)
def cellCountingJIT(rIdCell,nCells):
    cellCounts = np.empty(nCells)
    for idCell in range(0,nCells):
        cellCounts[idCell] = np.sum(idCell == rIdCell)
    return cellCounts

@jit(nopython=True)
def cellMeanDataJIT(rIdCell,nCells,varData):
    cellMean = np.empty(nCells)
    for idCell in range(0,nCells):
        dataInCell = (idCell == rIdCell)*varData
        if dataInCell.size == 0:
            cellMean[idCell] = 0
        else:
           cellMean[idCell] = np.sum(dataInCell)/dataInCell.size
    return cellMean


class RectangularGridBase:
    
    def __init__(self,xml_recipe, xml_file, dims = ['depth','latitude','longitude']):
        self.xml_recipe = xml_recipe
        self.xml_file = xml_file
        self.grid = []
        self.cellCenters = []
        self.cellArea = []
        self.cellVolume = []
        self.dims = dims
        self.coords = {}
        self.countsInCell = []
        self.meanDataInCell=[]
        self.rIdCell = []
    
    
    def getGrid(self):
        root = ET.parse(self.xml_recipe).getroot()
        self.grid = len(self.dims)*[[]]
        for parameter in root.findall('EulerianMeasures/gridDefinition/'):
            if parameter.tag == 'BoundingBoxMin':
                x_min = np.float(parameter.get('x'))
                y_min = np.float(parameter.get('y'))
                z_min = np.float(parameter.get('z'))
            else:
                root_global = ET.parse(self.xml_file).getroot()
                bbox_min = root_global.find('caseDefinitions/simulation/BoundingBoxMin')
                x_min = np.float(bbox_min.get('x'))
                y_min = np.float(bbox_min.get('y'))
                z_min = np.float(bbox_min.get('z'))                
            if parameter.tag == 'BoundingBoxMax':
                x_max = np.float(parameter.get('x'))
                y_max = np.float(parameter.get('y'))
                z_max = np.float(parameter.get('z'))
            else:
                root_global = ET.parse(self.xml_file).getroot()
                bbox_max = root_global.find('caseDefinitions/simulation/BoundingBoxMax')
                x_max = np.float(bbox_max.get('x'))
                y_max = np.float(bbox_max.get('y'))
                z_max = np.float(bbox_max.get('z'))                
            if parameter.tag == 'resolution':
                x_step = np.float(parameter.get('x'))
                y_step = np.float(parameter.get('y'))
                z_step = np.float(parameter.get('z'))
            if parameter.tag == 'units':
                units_value = parameter.get('value')                
        
        print('Domain limits: ',x_min,x_max,y_min,y_max,z_min,z_max)
        
        if units_value == 'degrees':
            self.grid[2] = np.arange(x_min,x_max,x_step)
            self.grid[1] = np.arange(y_min,y_max,y_step)
            self.grid[0] = np.arange(z_min,z_max,z_step)
                      
        elif units_value == 'relative':
            
             self.grid[2] = np.linspace(x_min,x_max,np.int(x_step+1))
             self.grid[1] = np.linspace(y_min,y_max,np.int(y_step+1))
             self.grid[0] = np.linspace(z_min,z_max,np.int(z_step+1))
    
        elif units_value == 'meters':
            y_c = (y_max+y_min)/2.
            dlat = y_step/(cte.degreesToRad*cte.earthRadius)
            dlon = x_step/(cte.degreesToRad*cte.earthRadius * np.cos(cte.degreesToRad*(y_c)))
            self.grid[2] = np.arange(x_min,x_max,dlon)
            self.grid[1] = np.arange(y_min,y_max,dlat)
            self.grid[0] = np.arange(z_min,z_max,z_step)
        
        return
    
    def getCellCenters(self):
        for value in self.grid:
            self.cellCenters.append((value[:-1] + value[1:])/2.)
        return
    
    
    def getCellAreas(self):
            # check in order to reverse the axis dimensions
            #depths, lats, lons = np.meshgrid(self.grid['depth'],self.grid['latitude'],self.grid['longitude'],indexing='ij')
            #dx = (lons[1:]-lons[:-1])*(np.pi/180.)*6371837. * np.cos((np.pi/180.)*((lats[:-1] + lats[1:])/2.))
            dlon = (self.grid[2][1:]-self.grid[2][:-1])
            dlat = (self.grid[1][1:]-self.grid[1][:-1])
            y_c = (self.grid[1][1:]+self.grid[1][:-1])/2.
            dx = dlon[np.newaxis,:]*(cte.degreesToRad*cte.earthRadius * np.cos(cte.degreesToRad*(y_c[:,np.newaxis])))
            dy = dlat*cte.degreesToRad*cte.earthRadius
            self.cellArea = dx*dy[:,np.newaxis]
    
    def getCellVolumes(self):
        dz = self.grid[0][1:]-self.grid[0][:-1]
        self.cellVolume = dz[:,np.newaxis,np.newaxis]*self.cellArea[np.newaxis,:,:]
    
    
    def getCoords(self):
        self.coords = {self.dims[0]:([self.dims[0]],self.cellCenters[0]),
                  self.dims[1]:([self.dims[1]],self.cellCenters[1]),
                  self.dims[2]:([self.dims[2]],self.cellCenters[2])
                  }
        
    
    def initializeGrid(self):
        self.getGrid()
        self.getCellCenters()
        self.getCoords()
        self.getCellAreas()
        self.getCellVolumes()


    def PositionsToIdCell(self,particlePositions):
        
        nz = self.cellCenters[0].size
        ny = self.cellCenters[1].size
        nx = self.cellCenters[2].size
        if np.size(particlePositions.shape) == 1: particlePositions = particlePositions[np.newaxis,:]
        z_dig = np.digitize(particlePositions[:,0],self.grid[0],right=True)
        y_dig = np.digitize(particlePositions[:,1],self.grid[1],right=True)
        x_dig = np.digitize(particlePositions[:,2],self.grid[2],right=True)

#        z_dig = np.int32((particlePositions[:,0] - min(self.grid[0]))/abs(self.grid[0][1]-self.grid[0][0]))
#        z_dig[z_dig >= (nz-1)] = nz-1 
#        z_dig[0 > z_dig] = 0 
#        y_dig = np.int32((particlePositions[:,1] - min(self.grid[1]))/abs(self.grid[1][1]-self.grid[1][0]))
#        y_dig[y_dig >= (ny-1)] = ny-1 
#        y_dig[0 > y_dig] = 0
#        x_dig = np.int32((particlePositions[:,2] - min(self.grid[2]))/abs(self.grid[2][1]-self.grid[2][0]))
#        x_dig[x_dig >= (nx-1)] = nx-1 
#        x_dig[0 > x_dig] = 0 

        self.rIdCell = np.ravel_multi_index((z_dig,y_dig,x_dig),(nz,ny,nx),mode='clip')
        #print(self.rIdCell) 
      
    def getCountsInCell(self,particlePositions):
        nz = self.cellCenters[0].size
        ny = self.cellCenters[1].size
        nx = self.cellCenters[2].size
        self.PositionsToIdCell(particlePositions)
        nCells = nx*ny*nz
        self.countsInCell = np.reshape(cellCountingJIT(self.rIdCell,nCells),(nz,ny,nx))
          
    def getMeanDataInCell(self,varData):
        nz = self.cellCenters[0].size
        ny = self.cellCenters[1].size
        nx = self.cellCenters[2].size
        nCells = nx*ny*nz
        cellMean = cellMeanDataJIT(self.rIdCell,nCells,varData)
        self.meanDataInCell = np.reshape(cellMean,(nz,ny,nx))

    def shape(self):
        return map(len,self.grid)
        

    


class RectangularGridBaseTime(RectangularGridBase):
    
    def __init__(self,TimeGrid):
        self.cellCountsT = []
        self.activeCells = []
        self.timeAxis = []
        
    
    def shape(self):
        nr = map(len,self.grid)
        nt = len(self.timeAxis.time)
        return nt+nr 
        
    def setnCounts(self):
        self.cellCounts = np.zeros((self.shape))
    
    def getCountsPerCell(self,particlePositions,timeInstant):
        self.cellCounts[timeInstant] = self.cellCounting(particlePositions)
    
    def getMeanDataInCell(self,varData,timeInstant):
        self.cellMean[timeInstant] = self.cellMeanData(varData)
    
    def getCoords(self):
        self.getCoords()
        self.coords['time'] = (['time'],self.timeAxis.time) 


