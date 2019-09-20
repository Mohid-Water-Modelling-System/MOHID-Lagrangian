# -*- coding: utf-8 -*-

import h5py
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
    def __init__(self, xmlCaseFile, xmlRecipe, dataDir):
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
        
        
    def convertFiles(self):
        outdirLocal = self.dataDir + '/postProcess_'+os.path.basename(os_dir.filename_without_ext(self.xmlRecipe))+'_hdf5'                
        if os.path.exists(outdirLocal):
            os_dir.deleteDirForce(outdirLocal)
        os.mkdir(outdirLocal)
        f=0
        print('-> Converting .vtu to .hdf5, MOHID formated')
        for file in self.vtuFiles:
            if self.timeMask[f]:
                #convert the file
                vtuFile = vtuParser.VTUParser(file)
                r = vtuFile.points('coords')
                #print(r[1,:]) #[depth, lat, lon] of a tracer
                hdf5FileName = os_dir.filename_without_ext(os.path.basename(vtuFile.fileName))+'.hdf5'
                print('--> '+os.path.basename(vtuFile.fileName)+' -> '+hdf5FileName)
                with h5py.File(outdirLocal + '/' + hdf5FileName ,'a') as hdf5File:
                    #main groups
                    grid = hdf5File.create_group("Grid")
                    results = hdf5File.create_group("Results")
                    time = hdf5File.create_group("Time")
                    #subgroups
                    group1 = results.create_group("Group_1/Data_1D")
                    #writing data
                    lon = group1.create_dataset('Longitude/Longitude_00001', data=r[:,2], dtype='f')
                    lon.attrs['Maximum'] = max(r[:,2])
                    lon.attrs['Minimum'] = min(r[:,2])
                    lon.attrs['Units'] = 'º'
                    lat = group1.create_dataset('Latitude/Latitude_00001', data=r[:,1], dtype='f')
                    lat.attrs['Maximum'] = max(r[:,1])
                    lat.attrs['Minimum'] = min(r[:,1])
                    lat.attrs['Units'] = 'º'
                    zz = group1.create_dataset('Z Pos/Z Position_00001', data=r[:,0], dtype='f')
                    zz.attrs['Maximum'] = max(r[:,0])
                    zz.attrs['Minimum'] = min(r[:,0])
                    zz.attrs['Units'] = 'm'
                    r = vtuFile.points('source')
                    source = group1.create_dataset('Origin ID/Origin ID_00001', data=r, dtype='f')
                    source.attrs['Maximum'] = max(r)
                    source.attrs['Minimum'] = min(r)
                    source.attrs['Units'] = '-'
                    #writing time
                    dateArray = MDateTime.getMOHIDDateFromTimeStamp(self.time[f])                    
                    date = time.create_dataset('Time_00001', data=dateArray, dtype='f')
                    date.attrs['Maximum'] = max(dateArray)
                    date.attrs['Minimum'] = min(dateArray)
                    date.attrs['Units'] = 'YYYY/MM/DD HH:MM:SS'                
            f = f+1


def run(caseXML, recipeXML, directory):
    converter = vtu2hdf5(caseXML, recipeXML, directory)
    converter.convertFiles()
    
    
#run('C:/Users/RBC_workhorse/Documents/GitHub/MOHID-Lagrangian/RUN_Cases/Arousa_2D_test_case/Arousa2D_case.xml','C:/Users/RBC_workhorse/Documents/GitHub/MOHID-Lagrangian/RUN_Cases/Arousa_2D_test_case/Post_scripts/PostRecipe_Arousa.xml','C:/Users/RBC_workhorse/Documents/GitHub/MOHID-Lagrangian/RUN_Cases/Arousa_2D_test_case/Arousa2D_case_out')