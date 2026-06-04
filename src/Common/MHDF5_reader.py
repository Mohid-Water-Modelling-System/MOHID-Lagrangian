# -*- coding: utf-8 -*-

#    !------------------------------------------------------------------------------
#    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
#    !------------------------------------------------------------------------------
#    !
#    ! TITLE         : MHDF5_reader
#    ! PROJECT       : Mohid python tools
#    ! MODULE        : background
#    ! URL           : http://www.mohid.com
#    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
#    ! DATE          : January 2019
#    ! REVISION      : Canelas 0.1
#    !> @author
#    !> Ricardo Birjukovs Canelas
#    !
#    ! DESCRIPTION:
#    !This class provides an API to read and extract data from MOHID hdf5 outputs.
#    !------------------------------------------------------------------------------
#    
#    MIT License
#    
#    Copyright (c) 2018 RBCanelas
#    
#    Permission is hereby granted, free of charge, to any person obtaining a copy
#    of this software and associated documentation files (the "Software"), to deal
#    in the Software without restriction, including without limitation the rights
#    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#    copies of the Software, and to permit persons to whom the Software is
#    furnished to do so, subject to the following conditions:
#    
#    The above copyright notice and this permission notice shall be included in all
#    copies or substantial portions of the Software.
#    
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#    SOFTWARE.

import h5py
import MDateTime

class MHDF5Reader:
    
    #given a file name for an hdf5 file and a directory where it is located
    #this constructor: 
    #opens the file
    #checks if it is a valid MOHID output
    #checks the type of MOHID file (hydro, lagrangian, etc)
    #checks if the version of eulerian files is high enough (must have Corners3D) 
    def __init__(self, fileName, directory, mandatoryMesh = True):
        self.fileName = fileName
        self.directory = directory
        self.f = h5py.File(self.directory +'/'+ self.fileName, 'r')
        self.validFile = 0
        self.fileType = []
        self.possibleFileTypes = ['Hydrodynamic', 'Hydrodynamic2D', 'Lagrangian', 'WaterProperties', 'WaterProperties2D', 'InterfaceSedimentWater', 'InterfaceWaterAir', 'Turbulence', 'Generic', 'Generic2D']
        self.fileKeys = self.f.keys()
        self.fileTimeSteps = list(self.f['Time'].keys())        
        
        self.MOHIDkeys = ['Grid', 'Results', 'Time']
        
        #check if file is a valid MOHID output
        if all(i in self.fileKeys for i in self.MOHIDkeys):
            self.validFile = 1
        else:
            for key in self.MOHIDkeys:
                if key not in self.fileKeys:
                    print('- [MHDF5Reader::init]: file does not have', key, 'group, not a MOHID output, ignoring')
        #check for file type
        if self.validFile == 1:
            #checking for Hydrodynamic files
            if 'water level' in list(self.f['Results'].keys()):
                self.fileType = 'Hydrodynamic'
                #Because 2D fiels are mixed with 3D fields
                exclusions = ['Error','TidePotential','water column','water level']                
                if self.getGeoDims() == 2:
                    self.fileType = 'Hydrodynamic2D'
                    exclusions = [] 
                if mandatoryMesh and self.fileType == 'Hydrodynamic':
                    if 'Corners3D' not in list(self.f['Grid'].keys()):
                        print('- [MHDF5Reader::init]: old hydrodynamic file, without mesh information, ignoring')
                        self.validFile = 0
                self.fVars = list(self.f['Results'].keys())                
                for exc in exclusions:
                    if exc in self.fVars:
                        self.fVars.remove(exc)
            #checking for WaterProperies files
            if 'temperature' in list(self.f['Results'].keys()):
                self.fileType = 'WaterProperties'
                #weird props
                exclusions = ['Assimila']
                if self.getGeoDims() == 2:
                    self.fileType = 'WaterProperties2D'
                    exclusions = []
                if mandatoryMesh and self.fileType == 'WaterProperties':
                    if 'Corners3D' not in list(self.f['Grid'].keys()):
                        print('- [MHDF5Reader::init]: old WaterProperties file, without mesh information, ignoring')
                        self.validFile = 0
                self.fVars = list(self.f['Results'].keys())                
                for exc in exclusions:
                    if exc in self.fVars:
                        self.fVars.remove(exc)
            #cheking for Lagrangian files 
            if 'Group_1' in list(self.f['Results'].keys()):
                self.fileType = 'Lagrangian'
                #storing all variables in 
                self.fVars = list(self.f['Results']['Group_1']['Data_1D'].keys())
                #Because conventions are not followed (name of the variable 
                #is not the name of the field, mixing diferent dimenisionalities on the same group,...)
                exclusions = ['X Pos','Y Pos','Z Pos','Latitude average','Longitude average', 'googlemaps_x_average', 'googlemaps_y_average']
                for exc in exclusions:
                    if exc in self.fVars:
                        self.fVars.remove(exc)
            #checking for interface files
            if 'Deposition' in list(self.f['Results'].keys()):
                self.fileType = 'InterfaceSedimentWater'
                self.validFile = 0
            if 'evaporation' in list(self.f['Results'].keys()):
                self.fileType = 'InterfaceWaterAir'
                self.validFile = 0
            #checking for turbulence files
            if 'Diffusivity' in list(self.f['Results'].keys()):
                self.fileType = 'Turbulence'
                self.validFile = 0
            #maybe it can be a generic file...
            if self.fileType == []:
                self.fileType = 'Generic'
                exclusions = []                
                if self.getGeoDims() == 2:
                    self.fileType = 'Generic2D'
                    exclusions = []
                self.fVars = list(self.f['Results'].keys())                
                for exc in exclusions:
                    if exc in self.fVars:
                        self.fVars.remove(exc)
                
                
    def isValidFile(self):
        return self.validFile
                
    #returns the file type as a string
    def getFileType(self):
        if self.validFile == 1:
            return self.fileType
        else:
            print('- [MHDF5Reader::getfileType]: invalid file, no type, ignoring')
     
    #returns the number of time steps in the file
    def getNumbTimeSteps(self):
        if self.validFile == 1:
            return len(self.fileTimeSteps)
        else:
            print('- [MHDF5Reader::getNumbTimeSteps]: invalid file, ignoring')
    
    #returns the date of a time step in the file in string format
    def getDateStr(self, timeIndex):
        if self.validFile == 1:            
            return MDateTime.getDateStringFromMOHIDDate(self.getDate(timeIndex))
        else:
            print('- [MHDF5Reader::getDate]: invalid file, ignoring')
            
    #returns the date of a time step in the file in list format
    def getDate(self, timeIndex):
        if self.validFile == 1:
            date = list(self.f['Time'][self.fileTimeSteps[timeIndex-1]][:].transpose())
            return date
        else:
            print('- [MHDF5Reader::getDate]: invalid file, ignoring')
                
    #returns the geometry dimensionality
    def getGeoDims(self):
        if self.validFile == 1:
            if self.fileType != 'Lagrangian':
                if 'WaterPoints3D' in list(self.f['Grid'].keys()):
                    dims = self.f['Grid']['WaterPoints3D'].shape
                    if dims[0] == 1: #first layer, it's a 2D file
                        return 2
                    else:
                        return 3
                elif 'WaterPoints2D' in list(self.f['Grid'].keys()):
                    return 2
                else: 
                    self.validFile = 0
                    return 0
            if self.fileType == 'Lagrangian':
                return 3
        else:
            print('- [MHDF5Reader::getGeoDims]: invalid file, no geometry, ignoring')
            
    #returns an array with the mesh dimensions
    def getMeshDims(self, timeIndex):
        if self.validFile == 1:
            if self.fileType != 'Lagrangian' and self.getGeoDims() == 3:
                return self.f['Grid']['Corners3D']['Latitude'].shape
            if self.fileType != 'Lagrangian' and self.getGeoDims() == 2:
                return self.f['Grid']['Latitude'].shape
            if self.fileType == 'Lagrangian':
                timeVar = 'Latitude_' + str(timeIndex).zfill(5)
                return self.f['Results']['Group_1']['Data_1D']['Latitude'][timeVar].size
        else:
            print('- [MHDF5Reader::getMeshDims]: invalid file, no mesh geometry, ignoring')
            
    #returns a list with (name,attPath) with all variables
    def getAllAttributesPath(self, timeIndex):
        if self.validFile == 1:
            Attr = []
            if self.fileType != 'Lagrangian':
                for var in self.fVars:
                    timeVar = var + '_' + str(timeIndex).zfill(5)
                    pathVar = '/Results/'+var+'/'+timeVar
                    Attr.append([var, pathVar])
            if self.fileType == 'Lagrangian':
                for var in self.fVars:
                    timeVar = var + '_' + str(timeIndex).zfill(5)
                    pathVar = '/Results/Group_1/Data_1D/'+var+'/'+timeVar
                    Attr.append([var, pathVar])
            return Attr
        else:
            print('- [MHDF5Reader::getAllAttributesPath]: invalid file ignoring')        
    
    #returns true if the file has bathymetry field
    def hasBathymetry(self):
        if self.validFile == 1:            
            if self.fileType != 'Lagrangian':
                if 'Bathymetry' in list(self.f['Grid'].keys()):
                    return 1
                else:
                    return 0
            if self.fileType == 'Lagrangian':
                return 0
        else:
            print('- [MHDF5Reader::hasBathymetry]: invalid file ignoring')