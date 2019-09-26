# -*- coding: utf-8 -*-

import vtk
from vtk.util.numpy_support import vtk_to_numpy
import glob
import numpy as np

class VTUParser:
    def __init__(self,vtu_file):
        self.fileName = vtu_file
        self.part_vars = ['coords']
        
        
    def points(self,var):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(self.fileName)
        reader.Update()
        nvars = reader.GetOutput().GetPointData().GetNumberOfArrays()
        for i in range(0,nvars):
            self.part_vars.append(reader.GetOutput().GetPointData().GetArrayName(i))
        if var == 'coords':
            vtu_vars = vtk_to_numpy(reader.GetOutput().GetPoints().GetData())[:,::-1]
        elif var == 'velocity':
            vtu_vars = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray(var))
            vtu_vars = np.sqrt(vtu_vars[:,0]**2 + vtu_vars[:,1]**2 + vtu_vars[:,2]**2)
        elif var in self.part_vars:
            vtu_vars = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray(var))
        return vtu_vars
    
   
def validVtuFilesList(directory):
    vtu_list = glob.glob(directory+'/*_?????.vtu')
    vtu_list.sort()
    vtu_list = vtu_list[1:]
    return vtu_list