# -*- coding: utf-8 -*-

import vtk
from vtk.util.numpy_support import vtk_to_numpy
import glob

class VTUParser:
    def __init__(self,vtu_file):
        self.vtu_file = vtu_file
        self.part_coords = ['longitude','latitude','depth']
        self.part_vars = ['coords','id','source','velocity']
        
        
    def points(self,var):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(self.vtu_file)
        reader.Update()
        if var == 'coords':
            vtu_vars = vtk_to_numpy(reader.GetOutput().GetPoints().GetData())[:,::-1]
        else:
            vtu_vars = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray(var))
        return vtu_vars
    
   
def validVtuFilesList(directory):
    vtu_list = glob.glob(directory+'/*_?????.vtu')
    vtu_list.sort()
    vtu_list = vtu_list[1:]            
    return vtu_list


#class VTUParser:
#    def __init__(self,vtu_file):
#        self.vtu_file = vtu_file
#        self.part_coords = ['longitude','latitude','depth']
#        self.part_vars = ['coords','id','source','velocity']
#        
#        
#    def points(self):
#        reader = vtk.vtkXMLUnstructuredGridReader()
#        reader.SetFileName(self.vtu_file)
#        reader.Update()       
#        vtu_vars = {}
#        for var in self.part_vars:
#            if var == 'coords':
#                dim = 0
#                for coord in self.part_coords:
#                    vtu_vars[coord] = vtk_to_numpy(reader.GetOutput().GetPoints().GetData())[:,dim]
#                    dim = dim +1
#            else:
#                vtu_vars[var] = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray(var))
#        return vtu_vars