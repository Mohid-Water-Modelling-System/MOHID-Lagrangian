# -*- coding: utf-8 -*-


import vtk
from vtk.util.numpy_support import vtk_to_numpy
import xml.etree.ElementTree as ET
import numpy as np
import xarray as xr
import os
import sys

basePath = os.path.dirname(os.path.realpath(__file__))
commonPath = os.path.abspath(os.path.join(basePath, "../Common"))
sys.path.append(commonPath)

import os_dir
import MDateTime

import vtuParser


class vtu2hdf5:
    def __init__(self,xmlCaseFile, xmlRecipe, dataDir):
        self.xmlCaseFile = xmlCaseFile
        self.xmlRecipe = xmlRecipe
        self.dataDir = dataDir
        self.vtuFiles = vtuParser.validVtuFilesList(dataDir)
        self.nVtuFiles = len(vtuParser.validVtuFilesList(dataDir))
        self.ISO_time_origin = MDateTime.getDateStringFromDateTime(MDateTime.BaseDateTime())
        self.time = []
        self.timeMask = []
        
        self.get_time_axis()
        
    def get_time_axis(self):
        root= ET.parse(self.xmlCaseFile).getroot()

        for parameter in root.findall('execution/parameters/parameter'):
            if parameter.get('key') == 'Start':
                self.start_time = parameter.get('value')
            if parameter.get('key') == 'End':
                self.end_time = parameter.get('value')
            if parameter.get('key') == 'OutputWriteTime':
                self.dt = np.float(parameter.get('value'))
        startTimeStamp = MDateTime.getTimeStampFromISODateString(self.start_time)
        self.time = np.array([startTimeStamp + i*self.dt/(3600.0*24.0) for i in range(0,self.nVtuFiles)])
        self.timeMask = np.ones(self.time.size,dtype=np.bool)
        self.get_time_xml_recipe()
        
        
    def get_time_xml_recipe(self):
         root= ET.parse(self.xmlRecipe).getroot()
         for parameter in root.findall('time/'):
            if parameter.tag == 'start':
                time_start = MDateTime.getTimeStampFromISODateString(parameter.get('value'))
                self.timeMask = self.timeMask & (self.time > time_start)
            if parameter.tag == 'end':
                time_end = MDateTime.getTimeStampFromISODateString(parameter.get('value'))
                self.timeMask = self.timeMask & (self.time < time_end)
        
        
    def convertFile(self):
        f=0
        for file in self.vtuFiles:
            if self.timeMask:
                #convert the file
                vtuFile = vtuParser.VTUParser(file)
                r = vtuFile.points('coords')
            f = f+1






def run(caseXML, recipeXML, directory):
    converter = vtu2hdf5(caseXML, recipeXML, directory)
    
    
run('bla','blabla','C:/Users/RBC_workhorse/Documents/GitHub/MOHID-Lagrangian/RUN_Cases/Tagus3D_case/Tagus3D_out')