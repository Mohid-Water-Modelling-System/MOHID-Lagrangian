"""
Created on Wed Jun  3 12:29:02 2020

@author: gfnl143
"""

import vtk
import glob
import numpy as np
from src.VTUReader import getSourceMaskFromVTU, getBeachMaskFromVTU
from src.VTUReader import getVariableArrayFromVTU


class VTUParser:

    def __init__(self, outDir):
        self.outDir = outDir
        self.fileList, self.parentFile = self.getfileList()
        self.part_vars = ['coords']
        self.vtkReader = vtk.vtkXMLUnstructuredGridReader()
        self.nvars = self.getNumberOfVars()
        self.availableVtuVars = self.getAvailableVars()

    def getfileList(self):
        vtuList = glob.glob(self.outDir+'/*_?????.vtu')
        vtuList.sort()
        parentFile = vtuList[0]
        fileList = vtuList[1:]
        return fileList, parentFile

    def updateFileList(self, fileList):
        self.fileList = fileList

    def getNumberOfVars(self):
        self.vtkReader.SetFileName(self.parentFile)
        self.vtkReader.Update()
        number_of_arrays = self.vtkReader.GetOutput().GetPointData().GetNumberOfArrays()
        return number_of_arrays

    def getAvailableVars(self) -> list:
        self.vtkReader.SetFileName(self.parentFile)
        self.vtkReader.Update()
        variableList = []
        for i in range(0, self.nvars):
            variableList.append(self.vtkReader.GetOutput().GetPointData().GetArrayName(i))
        return variableList

    def updateReaderWithFile(self, fileName):
        self.vtkReader.SetFileName(fileName)
        self.vtkReader.Update()
        self.nvars = self.getNumberOfVars()
        self.availableVtuVars = self.getAvailableVars()

    def getVariableData(self, variableName, source='global', beachCondition=None):
        sourceMask = getSourceMaskFromVTU(self.vtkReader, source)
        beachMask = getBeachMaskFromVTU(self.vtkReader, beachCondition)
        vtuVarArray = getVariableArrayFromVTU(self.vtkReader, variableName)

        if source or beachCondition:
            if ((np.size(sourceMask & beachMask) == 1) and
               ((sourceMask & beachMask) == True)):
                return vtuVarArray
            else:
                vtuVarArray = vtuVarArray[sourceMask*beachMask]

        return vtuVarArray
