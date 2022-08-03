"""
Module to read vtu files. 
"""

import vtk
import glob
import numpy as np
from src.VTUReader import getSourceMaskFromVTU, getBeachMaskFromVTU
from src.VTUReader import getVariableFromVTU


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

    def getNumberOfVars(self, file=None) -> int:
        """Get the number of variables available in a vtu file.

        Args:
            file (str, optional): name of the input file. Defaults to None.

        Returns:
            int: number of available variables.

        """

        if file is None:
            self.vtkReader.SetFileName(self.parentFile)
            self.vtkReader.Update()
        else:
            self.vtkReader.SetFileName(file)
            self.vtkReader.Update()
        number_of_arrays = self.vtkReader.GetOutput().GetPointData().GetNumberOfArrays()
        return number_of_arrays

    def getAvailableVars(self, file=None) -> list:
        """Get the names of the variables available in a vtu files.

        Args:
            file (str, optional): name of the input file. Defaults to None.

        Returns:
            int: number of available variables.

        """
        if file is None:
            self.vtkReader.SetFileName(self.parentFile)
            self.vtkReader.Update()
        else:
            self.vtkReader.SetFileName(file)
            self.vtkReader.Update()
        variableList = []
        for i in range(0, self.nvars):
            variableList.append(self.vtkReader.GetOutput().GetPointData().GetArrayName(i))
        return variableList

    def updateReaderWithFile(self, fileName):
        """Updates the VTUParser reader with the provided filename attributes."""
        self.vtkReader.SetFileName(fileName)
        self.vtkReader.Update()
        self.nvars = self.getNumberOfVars(fileName)
        self.availableVtuVars = self.getAvailableVars(fileName)

    def getVariableData(self, variableName, source='global', beachCondition=None):
        """Reads the variable from the current VTU filename.

        Args:
            variableName (str): available vtu variable.
            source (str, optional): source name to read. Defaults to 'global'.
            beachCondition (str, optional): '0,1,2'. . Defaults to None.

        Returns:
            vtuVarArray (np.array): Array with variable choosen

        """
        
        sourceMask = getSourceMaskFromVTU(self.vtkReader, source)
        beachMask = getBeachMaskFromVTU(self.vtkReader, beachCondition)
        vtuVarArray = getVariableFromVTU(self.vtkReader, variableName)

        if source or beachCondition:
            if ((np.size(sourceMask & beachMask) == 1) and
               ((sourceMask & beachMask) == True)):
                return vtuVarArray
            else:
                vtuVarArray = vtuVarArray[sourceMask*beachMask]
        return vtuVarArray
