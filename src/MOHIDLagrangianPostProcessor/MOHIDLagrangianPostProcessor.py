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

import vtk
from vtk.util.numpy_support import vtk_to_numpy
import xml.etree.ElementTree as ET
import numpy as np
import xarray as xr
import os
import sys
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'



class VTUParser:
    def __init__(self,vtu_file):
        self.vtu_file = vtu_file
        self.part_coords = ['longitude','latitude','depth']
        self.part_vars = ['coords','age','id','landIntMask','source','velocity']
        
    def points(self):
        print('Reading vtu_file',self.vtu_file)
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(self.vtu_file)
        reader.Update()
       
        vtu_vars = {}
        for var in self.part_vars:
            if var == 'coords':
                dim = 0
                for coord in self.part_coords:
                    vtu_vars[coord] = vtk_to_numpy(reader.GetOutput().GetPoints().GetData())[:,dim]
                    dim = dim +1
            else:
                vtu_vars[var] = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray(var))
        return vtu_vars
        
        
    
        
class PVDParser:
    
    def __init__(self,pvd_file):
        self.pvd_file = pvd_file
        self.vtu_list = []
        self.files = []
        self.timesteps = []
        self.vtu_data = []
    
    def get_vtu_files(self):
        tree = ET.parse(self.pvd_file)
        self.vtu_list = tree.find('Collection')[:]
        for vtu_file in self.vtu_list:
            self.files.append(vtu_file.attrib['file'])
            self.timesteps.append(float(vtu_file.attrib['timestep']))
            self.vtu_data.append(VTUParser(vtu_file.attrib['file'])) 
    



class GridBasedMeasures:
    def __init__(self,xml_file):
        self.xml_file = xml_file
        self.pvd_file = self.xml_file.replace('.xml','.pvd')
        self.pvd_data = PVDParser(self.pvd_file)
        self.sources = {'id':{}}
        self.grid_steps = [50.,0.005,0.005]
        self.ISO_time_origin = '1950-01-01 00:00:00'
        self.grid  = {}
        self.centers = {}
        self.ds = []
        self.area = []
        self.volume = []
        self.variables = {}
        self.dims = ['time','depth','latitude','longitude']
        self.netcdf_output_file = self.pvd_file.replace('.pvd','.nc')
        
    def get_pvd(self):
        self.pvd_data.get_vtu_files()
            
    
    def get_time_axis(self):
        tree = ET.parse(self.xml_file)
        root = tree.getroot()
        self.start_time = root.find('execution').find('parameters')[0].attrib['value']
        self.end_time = root.find('execution').find('parameters')[1].attrib['value']
        self.time = np.array(list(map(float,self.pvd_data.timesteps)))
        self.dt = np.float(root.find('caseDefinitions').find('simulation').find('timestep').attrib['dt'])
    
    
    def get_sources(self):
        tree = ET.parse(self.xml_file)
        self.sources['id'] = {}
        for source in tree.find('caseDefinitions').find('sourceDefinitions')[:]:
            id_source = source.find('setsource').attrib['id']
            self.sources['id'] = {id_source: None}
            # at the moment we don't require specific information just the ids.
            # In a future releases if we are going to compute measures using 
            # integration time, we will 
            
            
        
        
    def get_grid(self):
        tree = ET.parse(self.xml_file)
        root = tree.getroot()
        min_dict = root.find('caseDefinitions').find('simulation').find('BoundingBoxMin').attrib
        max_dict = root.find('caseDefinitions').find('simulation').find('BoundingBoxMax').attrib
        
        x_min, x_max = np.float(min_dict['x']),np.float(max_dict['x'])
        y_min, y_max = np.float(min_dict['y']),np.float(max_dict['y'])
        z_min, z_max = np.float(min_dict['z']),np.float(max_dict['z'])
        
        print('Domain limits:',x_min,x_max,y_min,y_max,z_min,z_max)
        
        self.grid ={'longitude': np.arange(x_min,x_max,self.grid_steps[2]),
                    'latitude': np.arange(y_min,y_max,self.grid_steps[1]),
                    'depth': np.arange(z_min,z_max,self.grid_steps[0])                   
                    }
                    
        for key,value in self.grid.items():
            self.centers[key] = (value[:-1] + value[1:])/2.
        return
        
    
    def get_areas(self):
        # check in order to reverse the axis dimensions
        #depths, lats, lons = np.meshgrid(self.grid['depth'],self.grid['latitude'],self.grid['longitude'],indexing='ij')
        #dx = (lons[1:]-lons[:-1])*(np.pi/180.)*6371837. * np.cos((np.pi/180.)*((lats[:-1] + lats[1:])/2.))
        
        dx = (self.grid['longitude'][1:]-self.grid['longitude'][:-1])
        dy = (self.grid['latitude'][1:]-self.grid['latitude'][:-1])
        y_c = (self.grid['latitude'][1:]+self.grid['latitude'][:-1])/2.
        dlon = dx[np.newaxis,:]*(np.pi/180.)*6371837. * np.cos((np.pi/180.)*(y_c[:,np.newaxis]))

        dlat = dy*(np.pi/180.)*6371837.
        dz = self.grid['depth'][1:]-self.grid['depth'][:-1]
        self.area = dlon*dlat[:,np.newaxis]
        self.volume = dz[:,np.newaxis,np.newaxis]*self.area[np.newaxis,:,:]
        return
    
        
    def get_netcdf_header(self):
        coords = {'time':('time',self.time),
                  'depth': ('depth',self.centers['depth']),
                  'latitude' : ('latitude', self.centers['latitude']),
                  'longitude': ('longitude',self.centers['longitude']),
                  }
        
        self.ds = xr.Dataset(None,coords=coords)
        
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
        'valid_min': np.min(self.grid['depth']),\
        'valid_max': np.max(self.grid['depth'])}
        
        time_atts = {'long_name':'time',
             'units':'seconds since 1950-01-01 00:00:00'}
        

        self.ds.longitude.attrs = lon_attributtes
        self.ds.latitude.attrs = lat_attributtes
        self.ds.depth.attrs = depth_attributtes
        self.ds.time.attrs = time_atts
        
        self.ds.to_netcdf(self.netcdf_output_file)
        
        
        return

    
    def counts(self,source = None):
        # counts 2d and 3d are splitted in two functions. 
        nz,ny,nx = [np.size(self.centers[key]) for key in ['depth','latitude','longitude']]
        nt = len(self.pvd_data.vtu_data)
        counts_t = np.zeros((nt,nz,ny,nx))
        
        if nz > 1:
            i = 0
            for vtu_step in self.pvd_data.vtu_data:
                position = vtu_step.points()
                r = np.c_[position['depth'], position['latitude'], position['longitude']]
                if source:
                    source_mask = vtu_step.points()['source'] == int(source)
                    r = r[source_mask]
                bins = [self.grid['depth'],self.grid['latitude'],self.grid['longitude']]
                counts_t[i], _ = np.histogramdd(r,bins=bins) 
                i=i+1
        else:
            i = 0 
            for vtu_step in self.pvd_data.vtu_data:
                position = vtu_step.points()
                r = np.c_[position['latitude'],['longitude']]
                if source:
                    source_mask = vtu_step.points()['source'] == int(source)
                    r = r[source_mask]
                bins = [self.grid['latitude'],self.grid['longitude']]
                print('Processing',r.shape)
                counts_t[i], _ = np.histogramdd(r,bins=bins) 
                i=i+1
        
        return counts_t
        
        
            

        
    def residence_time(self):
        
        
        # Read netcdf, compute concentrations, compute residence tiem
        # compute global 
        ds = xr.open_dataset(self.netcdf_output_file)
        
        counts_t = self.counts()
        self.residence_time = np.zeros(counts_t.shape[1:])
        for t_s in range(0,counts_t.shape[0]):
            self.residence_time = (counts_t[t_s]>0) * self.dt + self.residence_time
            
        ds['residence_time_global'] = (['depth','latitude','longitude'], self.residence_time)
        
        # Read netcdf, compute concentrations, compute residence tiem
        # compute per source
        for source in self.sources['id'].keys():
            
            
            counts_t = self.counts(source=source)
            self.residence_time = np.zeros(counts_t.shape[1:])
            for t_s in range(0,counts_t.shape[0]): 
                self.residence_time = (counts_t[t_s] > 0) * self.dt + self.residence_time   
                
            var_name ='residence_time_source_' + source.zfill(3)
            ds[var_name] = (self.dims[1:],self.residence_time)
            ds[var_name].attrs = {'long_name':'residence_time',
              'units':'s'}
            
        
        ds.to_netcdf(self.netcdf_output_file,'a')
        return

        
            
    
    def concentrations(self):
        ds = xr.open_dataset(self.netcdf_output_file)
        
        counts_t = self.counts()        
        conc_area = counts_t.sum(axis=1)/self.area
        conc_volume = counts_t/self.volume
        
        ds['concentration_area'] = (['time','latitude','longitude'], conc_area)
        ds['concentration_volume'] = (self.dims, conc_volume)
        
                
        for source in self.sources['id'].keys():
           
            counts_t = self.counts(source=source)
            conc_area = counts_t.sum(axis=1)/self.area
            conc_volume = counts_t/self.volume
            
            var_name = 'concentration_area_source_' + source.zfill(3)
            ds[var_name] = (['time','latitude','longitude'], conc_area)
            ds[var_name].attrs = {'long_name':'concentration',
                                      'units':'ppm*m'}
            var_name = 'concentration_volume_source_' + source.zfill(3)
            ds[var_name] = (self.dims, conc_volume)
            ds[var_name].attrs = {'long_name':'concentration',
                          'units':'ppm*m*m'}
        
        ds.to_netcdf(self.netcdf_output_file,'a')
        return

        
            
 
    def run_postprocessing(self,measures):
        self.get_pvd()
        self.get_time_axis()
        self.get_grid()
        self.get_sources()
        self.get_areas()
        self.get_netcdf_header()
        if 'concentrations' in measures:
            self.concentrations()
        if 'residence_time' in measures:
            self.residence_time()

        

def main(xml_file, measures):
    post = GridBasedMeasures(xml_file)
    post.run_postprocessing(measures)

#xml_file = 'Tagus3D.xml'
#measures = ['residence_time','concentrations']
#post = GridBasedMeasures(xml_file)
#post.run_postprocessing(measures)
xml = sys.argv[1]
measures = [measure for measure in sys.argv[2:]]
print('XML_FILE:',xml, 'MEASURES',measures)
main(xml, measures)
