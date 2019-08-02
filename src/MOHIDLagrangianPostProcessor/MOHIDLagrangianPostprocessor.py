# -*- coding: utf-8 -*-

#    !------------------------------------------------------------------------------
#    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
#    !------------------------------------------------------------------------------
#    !
#    ! TITLE         : MOHIDLagrangianPreProcessor
#    ! PROJECT       : MOHIDLagrangian
#    ! URL           : http://www.mohid.com
#    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
#    ! DATE          : April 2019
#    ! REVISION      : Canelas 0.1
#    !> @author
#    !> Angel Daniel Garaboa Paz Canelas
#    !
#    ! DESCRIPTION:
#    ! Preprocessing script for MOHID Lagrangian. Lists input files, composes config 
#    ! files, etc 
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

import os
import vtk
import sys
from vtk.util.numpy_support import vtk_to_numpy
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr



class VTUParser:
    def __init__(self,vtu_file):
        self.vtu_file = vtu_file
        
    def points(self):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(self.vtu_file)
        reader.Update()
        my_vtk_array = reader.GetOutput().GetPoints().GetData()
        return vtk_to_numpy(my_vtk_array)
    
        
class PVDParser:
    
    def __init__(self,pvd_file):
        self.pvd_file = pvd_file
        self.vtu_list = []
        self.files = []
        self.timesteps = []
        self.vtu_data = []
    
    def get_vtu_file_info(self):
        tree = ET.parse(self.pvd_file)
        self.vtu_list = tree.getroot()[0]
        for vtu_file in self.vtu_list:
            self.files.append(vtu_file.attrib['file'])
            self.timesteps.append(vtu_file.attrib['timestep'])
            self.vtu_data.append(VTUParser(vtu_file.attrib['file']))             
    

class GridBasedMeasures:
    def __init__(self,pvd_file):
        self.pvd_file = pvd_file
        self.xml_output_file = self.pvd_file.replace('.pvd','.xml')
        self.pvd_data = PVDParser(pvd_file)
        self.grid_steps = [50.,0.0025,0.0025]
        self.ISO_time_origin = '1950-01-01 00:00:00'
        self.coords = []
        self.dims = ['time','depth','latitude','longitude']
        self.dims_2d = ['time','latitude','longitude']
        self.grid  = []
        self.centers = []
        self.ds = []
        
    def get_pvd(self):
        self.pvd_data.get_vtu_file_info()
    
        
    def get_time_axis(self):
        tree = ET.parse(self.xml_output_file)
        root = tree.getroot()
        self.start_time = root.find('execution').find('parameters')[0].attrib['value']
        self.end_time = root.find('execution').find('parameters')[1].attrib['value']
        self.time = np.array(list(map(float,self.pvd_data.timesteps)))
        self.dt = np.float(root.find('caseDefinitions').find('simulation').find('timestep').attrib['dt'])
    
    
    def get_domain_grid(self):
        tree = ET.parse(self.xml_output_file)
        root = tree.getroot()
        min_dict = root.find('caseDefinitions').find('simulation').find('BoundingBoxMin').attrib
        max_dict = root.find('caseDefinitions').find('simulation').find('BoundingBoxMax').attrib
        
        x_min, x_max = np.float(min_dict['x']),np.float(max_dict['x'])
        y_min, y_max = np.float(min_dict['y']),np.float(max_dict['y'])
        z_min, z_max = np.float(min_dict['z']),np.float(max_dict['z'])
        
        self.grid ={'longitude': np.arange(x_min,x_max,self.grid_steps[2]),
                    'latitude': np.arange(y_min,y_max,self.grid_steps[1]),
                    'depth': np.arange(z_min,z_max,self.grid_steps[0])                   
                    }
                    
        self.grid=[np.arange(z_min,z_max,self.grid_steps[0]),\
                   np.arange(y_min,y_max,self.grid_steps[1]),
                   np.arange(x_min,x_max,self.grid_steps[2])]
        
        self.centers = [(array[:-1] + array[1:])/2. for array in self.grid]
        
    
    def netcdf_header(self):
        coords = {'time':('time',self.time),
                  'depth': ('depth',self.centers[0]),
                  'latitude' : ('latitude', self.centers[1]),
                  'longitude': ('longitude',self.centers[2]),
                  }
                  
        self.ds = xr.Dataset(None,coords=coords)
    
    
    def concentrations_2d(self):
        nz,ny,nx = list(map(np.size,self.centers))
        nt = len(self.pvd_data.vtu_data)
        self.counts_t = np.zeros((nt,ny,nx))
        self.residence_time = np.zeros((ny,nx))
        
        i = 0 
        for vtu_step in self.pvd_data.vtu_data:
            counts, _, _ = np.histogram2d(vtu_step.points()[:,0], vtu_step.points()[:,1], bins=(self.grid[-1],self.grid[-2]))
            counts = counts.T
            self.counts_t[i] = counts            
            counts[counts > 0] = 1.
            self.residence_time = counts * self.dt + self.residence_time
            i += 1
            
        self.ds['concentration_2d'] = (self.dims_2d, self.counts_t)
        self.ds['residence_time_2d'] = (self.dims_2d[1:], self.residence_time)
    
    
    def concentrations_3d(self):
        nz,ny,nx = list(map(np.size,self.centers))
        nt = len(self.pvd_data.vtu_data)
        self.counts_t = np.zeros((nt,nz,ny,nx))
        self.residence_time = np.zeros((nz,ny,nx))
        
        i = 0 
        for vtu_step in self.pvd_data.vtu_data:
            counts, _, _ = np.histogramdd(vtu_step.points(), bins=(self.grid))
            counts = counts.T
            self.counts_t[i] = counts
            counts[counts > 0] = 1.
            self.residence_time = counts * self.dt + self.residence_time
            i += 1
        self.ds['concentration_3d'] = (self.dims,self.counts_t)
        self.ds['residence_time_3d'] = (self.dims[1:],self.residence_time)
    
    def to_netcdf(self):
        lon_attributtes = {'long_name': 'longitude',
        'standard_name': 'longitude',\
        'units': 'degrees_east',\
        '_FillValue': -9.8999995E15,\
        'valid_min': -180.0,\
        'valid_max': 180.0}
        
        lat_attributtes = {'long_name': 'latitude',
        'standard_name': 'latitude',\
        'units': 'degrees_north',\
        '_FillValue': -9.8999995E15,\
        'valid_min': -90.0,\
        'valid_max': 90.0}
        
        depth_attributtes = {'long_name': 'depth',
        'standard_name': 'depth',\
        'units': 'meters',\
        '_FillValue': -9.8999995E15,\
        'valid_min': np.min(self.grid[0]),\
        'valid_max': np.max(self.grid[1])}
        
        time_atts = {'long_name':'time',
             'units':'seconds since 1950-01-01 00:00:00'}
        
        concentration_atts = {'long_name':'concentration',
                              'units':'ppb'}
        
        residence_time_atts = {'long_name':'residence_time',
                      'units':'s'}
        
        
        self.ds.longitude.attrs = lon_attributtes
        self.ds.latitude.attrs = lat_attributtes
        self.ds.depth.attrs = depth_attributtes
        self.ds.time.attrs = time_atts
        if 'concentration_2d' in self.ds.keys():
            self.ds.concentration_2d.attrs = concentration_atts
        if 'concentration_3d' in self.ds.keys():
            self.ds.concentration_3d.attrs = concentration_atts
        if 'residence_time_2d' in self.ds.keys():
            self.ds.residence_time_2d.attrs = residence_time_atts
        if 'residence_time_3d' in self.ds.keys():
            self.ds.residence_time_3d.attrs = residence_time_atts

        self.ds.to_netcdf(self.pvd_file.replace('.pvd','.nc'))
    
    def run_postprocessing(self):
        self.get_pvd()
        self.get_time_axis()
        self.get_domain_grid()
        self.netcdf_header()
        self.concentrations_2d()
        self.to_netcdf()
        

def main(pvd_file):
    post = GridBasedMeasures(pvd_file)
    post.run_postprocessing()

main(sys.argv[0])
    

