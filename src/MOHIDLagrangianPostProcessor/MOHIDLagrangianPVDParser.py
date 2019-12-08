# -*- coding: utf-8 -*-


import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
import glob


class VTUParser:
    reader = vtk.vtkXMLUnstructuredGridReader()
    nvars = []
    parentFile = []
    availableVtuVars = []
    reader = []
           
    @classmethod
    def getVtuVarsFromInitialFile(cls,initialVtuFile):
        cls.parentFile = initialVtuFile
        cls.reader = vtk.vtkXMLUnstructuredGridReader()
        cls.reader.SetFileName(cls.parentFile)
        cls.reader.Update()
        cls.nvars = cls.reader.GetOutput().GetPointData().GetNumberOfArrays()               
        for i in range(0,cls.nvars):
            cls.availableVtuVars.append(cls.reader.GetOutput().GetPointData().GetArrayName(i))
    
    def __init__(self,vtuFile):
        self.fileName = vtuFile
        self.part_vars = ['coords']

                   
    def getVtuVariableData(self, variableName, source='global', beachCondition=None):
        VTUParser.reader.SetFileName(self.fileName)
        VTUParser.reader.Update()
        if source:
            if source == 'global': 
                sourceMask = np.bool(True)
            else:                
                sourceMask = vtk_to_numpy(VTUParser.reader.GetOutput().GetPointData().GetArray('source')) == np.int(source)
        if beachCondition:
            state = vtk_to_numpy(VTUParser.reader.GetOutput().GetPointData().GetArray('state'))
            if beachCondition == '0' : beachMask = True
            elif beachCondition == '1' : beachMask = state < 0.5
            elif beachCondition == '2' : beachMask = state >= 0.5
        else: 
            beachMask = np.bool(True)
            
        if variableName == 'coords':
            vtu_vars = vtk_to_numpy(VTUParser.reader.GetOutput().GetPoints().GetData())[:,::-1]
        elif variableName == 'velocity':
            vtu_vars = vtk_to_numpy(VTUParser.reader.GetOutput().GetPointData().GetArray(variableName))
            vtu_vars = np.sqrt(vtu_vars[:,0]**2 + vtu_vars[:,1]**2 + vtu_vars[:,2]**2)
        elif variableName in VTUParser.availableVtuVars:
            vtu_vars = vtk_to_numpy(VTUParser.reader.GetOutput().GetPointData().GetArray(variableName))
        
        if source or beachCondition:
            if ( (np.size(sourceMask & beachMask) == 1) and ((sourceMask & beachMask) == True)):
                return vtu_vars
            else:
                vtu_vars=vtu_vars[sourceMask*beachMask]
            
        return vtu_vars
            
            

class PVDParser:
    
    def __init__(self,pvdFile):
        self.pvdFile = pvdFile
        self.vtuFilelist = []
        self.vtuFileHandlers = []
        self.vtuParentFile = []
    
    def getVtuFileList(self,outDir):
        vtu_list = glob.glob(outDir+'/*_?????.vtu')
        vtu_list.sort()
        self.vtuParentFile = vtu_list[0]
        self.vtuFilelist = vtu_list[1:]
            
    def getVtuFileHandlers(self, outDir):
        self.getVtuFileList(outDir)
        VTUParser.getVtuVarsFromInitialFile(self.vtuParentFile)
        for vtuFile in self.vtuFilelist:
            self.vtuFileHandlers.append(VTUParser(vtuFile))
    
    def updateVtuFileHandlers(self, fileList):
        self.vtuFileHandlers = []
        self.vtuFilelist = fileList
        for vtuFile in self.vtuFilelist:
            self.vtuFileHandlers.append(VTUParser(vtuFile))



def test(outDir):
    pvd = PVDParser('Tagus3D.pvd')
    pvd.getVtuFileHandlers(outDir)
    for vtuInstance in pvd.vtuHeaders:
        print(vtuInstance.getVarData('coords'))
    

