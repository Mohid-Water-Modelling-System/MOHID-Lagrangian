# -*- coding: utf-8 -*-

#    !------------------------------------------------------------------------------
#    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
#    !------------------------------------------------------------------------------
#    !
#    ! TITLE         : MOHIDLagrangianPreProcessor
#    ! PROJECT       : MOHIDLagrangian
#    ! URL           : http://www.mohid.com
#    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
#    ! DATE          : June 2019
#    ! REVISION      : Garaboa 0.1
#    !> @author
#    !> Angel Daniel Garaboa Paz
#    !
#    ! DESCRIPTION:
#    ! Preprocessing script for MOHID Lagrangian. Lists input files, composes config 
#    ! files, etc 
#    !------------------------------------------------------------------------------
#    
#    MIT License
#    
#    Copyright (c) 2018 DGaraboa
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

from about import License
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import xml.etree.ElementTree as ET
import numpy as np
import xarray as xr
import os
import sys
import argparse
import glob
from scipy.interpolate import griddata
#from fastHistogram import histogramdd

# This environment variable avoids error on locking file when writing.
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

basePath = os.path.dirname(os.path.realpath(__file__))
commonPath = os.path.abspath(os.path.join(basePath, "../Common"))
sys.path.append(commonPath)

import os_dir
import MDateTime


class VTUParser:
    def __init__(self,vtu_file):
        self.vtu_file = vtu_file
        self.part_coords = ['longitude','latitude','depth']
        self.part_vars = ['coords','id','source','velocity']
        
        
    def points(self):
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
             

def validVtuFilesList(directory):
    vtu_list = glob.glob(directory+'/*_?????.vtu')
    vtu_list.sort()
    vtu_list = vtu_list[1:]            
    return vtu_list

        
class PVDParser:
    
    def __init__(self,pvd_file):
        self.pvd_file = pvd_file
        self.vtu_list = []
        self.files = []
        self.timesteps = []
        self.vtu_data = []
    
    def get_vtu_files(self, outDir):
        #tree = ET.parse(self.pvd_file)
        #self.vtu_list = tree.find('Collection')[:]
        self.vtu_list = validVtuFilesList(outDir)
        for vtu_file in self.vtu_list:
            self.files.append(vtu_file)
            #self.timesteps.append(float(vtu_file.attrib['timestep']))
            self.vtu_data.append(VTUParser(vtu_file))
    

class GridBasedMeasures:
    def __init__(self,xml_file, xml_recipe, outdir, outdirLocal):
        self.xml_file = xml_file
        self.xml_recipe = xml_recipe
        self.pvd_file = outdir +'/'+self.xml_file.replace('.xml','.pvd')
        self.pvd_data = PVDParser(self.pvd_file)
        self.nFiles = len(validVtuFilesList(outdir))
        self.sources = {'id':{}}
        self.grid_steps = [50.,0.005,0.005]
        self.ISO_time_origin = MDateTime.getDateStringFromDateTime(MDateTime.BaseDateTime())
        self.grid  = {}
        self.centers = {}
        self.ds = []
        self.area = []
        self.volume = []
        self.dims = ['time','depth','latitude','longitude']
        self.time = []
        self.timeMask = []
        self.netcdf_output_file = outdirLocal+'/'+os.path.basename(self.xml_file.replace('.xml','.nc'))
        if os.path.exists(outdirLocal):
            os_dir.deleteDirForce(outdirLocal)
        os.mkdir(outdirLocal)
        
        
    def get_pvd(self, outDir):
        self.pvd_data.get_vtu_files(outDir)
            
    
    def get_time_axis(self):
        root= ET.parse(self.xml_file).getroot()

        for parameter in root.findall('execution/parameters/parameter'):
            if parameter.get('key') == 'Start':
                self.start_time = parameter.get('value')
            if parameter.get('key') == 'End':
                self.end_time = parameter.get('value')
            if parameter.get('key') == 'OutputWriteTime':
                self.dt = np.float(parameter.get('value'))
        startTimeStamp = MDateTime.getTimeStampFromISODateString(self.start_time)
        self.time = np.array([startTimeStamp + i*self.dt/(3600.0*24.0) for i in range(0,self.nFiles)])
        self.timeMask = np.ones(self.time.size,dtype=np.bool)
        self.get_time_xml_recipe()
        
    
    def get_time_xml_recipe(self):
         root= ET.parse(self.xml_recipe).getroot()
         for parameter in root.findall('time/'):
            if parameter.tag == 'start':
                time_start = MDateTime.getTimeStampFromISODateString(parameter.get('value'))
                self.timeMask = self.timeMask & (self.time > time_start)
            if parameter.tag == 'end':
                time_end = MDateTime.getTimeStampFromISODateString(parameter.get('value'))
                self.timeMask = self.timeMask & (self.time < time_end)
    
    
    def get_sources(self):
        tree = ET.parse(self.xml_file)
        self.sources['id'] = {}
        dict_source = {}
        for source in tree.find('caseDefinitions/sourceDefinitions')[:]:
            id_source = source.find('setsource').attrib['id']
            dict_source[id_source] = None
        self.sources['id'] = dict_source
            # at the moment we don't require specific information just the ids.
            # In a future releases if we are going to compute measures using 
            # integration time, we will 
    
        
    def get_grid(self):
        root = ET.parse(self.xml_recipe).getroot()
        for parameter in root.findall('gridDefinition/'):
            if parameter.tag == 'BoundingBoxMin':
                x_min = np.float(parameter.get('x'))
                y_min = np.float(parameter.get('y'))
                z_min = np.float(parameter.get('z'))
            else:
                root_global = ET.parse(self.xml_file).getroot()
                bbox_min = root_global.find('caseDefinitions/simulation/BoundingBoxMin')
                x_min = np.float(bbox_min.get('x'))
                y_min = np.float(bbox_min.get('y'))
                z_min = np.float(bbox_min.get('z'))                
            if parameter.tag == 'BoundingBoxMax':
                x_max = np.float(parameter.get('x'))
                y_max = np.float(parameter.get('y'))
                z_max = np.float(parameter.get('z'))
            else:
                root_global = ET.parse(self.xml_file).getroot()
                bbox_max = root_global.find('caseDefinitions/simulation/BoundingBoxMax')
                x_max = np.float(bbox_max.get('x'))
                y_max = np.float(bbox_max.get('y'))
                z_max = np.float(bbox_max.get('z'))                
            if parameter.tag == 'resolution':
                x_step = np.float(parameter.get('x'))
                y_step = np.float(parameter.get('y'))
                z_step = np.float(parameter.get('z'))
            if parameter.tag == 'units':
                units_value = parameter.get('value')                
        
        print('Domain limits:',x_min,x_max,y_min,y_max,z_min,z_max)        
        if units_value == 'degrees':
            self.grid ={'longitude': np.arange(x_min,x_max,x_step),
                        'latitude': np.arange(y_min,y_max,y_step),
                        'depth': np.arange(z_min,z_max,z_step)                   
            }        
        elif units_value == 'relative':
            if z_step < 3: z_step = 3
            self.grid ={'longitude': np.linspace(x_min,x_max,np.int(x_step)),
                        'latitude': np.linspace(y_min,y_max,np.int(y_step)),
                        'depth': np.linspace(z_min,z_max,np.int(z_step))                   
            }        
        elif units_value == 'meters':
            y_c = (y_max+y_min)/2.
            dlat = y_step/((np.pi/180.)*6371837.)
            dlon = x_step/((np.pi/180.)*6371837. * np.cos((np.pi/180.)*(y_c)))
            self.grid['latitude'] = np.arange(y_min,y_max,dlat)
            self.grid['longitude'] = np.arange(x_min,x_max,dlon)
            self.grid['depth'] = np.arange(z_min,z_max,z_step)                                
        for key,value in self.grid.items():
            self.centers[key] = (value[:-1] + value[1:])/2.            
        return
        
    
    def get_areas(self):
        # check in order to reverse the axis dimensions
        #depths, lats, lons = np.meshgrid(self.grid['depth'],self.grid['latitude'],self.grid['longitude'],indexing='ij')
        #dx = (lons[1:]-lons[:-1])*(np.pi/180.)*6371837. * np.cos((np.pi/180.)*((lats[:-1] + lats[1:])/2.))
        
        dlon = (self.grid['longitude'][1:]-self.grid['longitude'][:-1])
        dlat = (self.grid['latitude'][1:]-self.grid['latitude'][:-1])
        y_c = (self.grid['latitude'][1:]+self.grid['latitude'][:-1])/2.
        dx = dlon[np.newaxis,:]*(np.pi/180.)*6371837. * np.cos((np.pi/180.)*(y_c[:,np.newaxis]))
        dy = dlat*(np.pi/180.)*6371837.
        dz = self.grid['depth'][1:]-self.grid['depth'][:-1]
        self.area = dx*dy[:,np.newaxis]
        self.volume = dz[:,np.newaxis,np.newaxis]*self.area[np.newaxis,:,:]
        return
    
        
    def get_netcdf_header(self):
        coords = {'time':('time',self.time[self.timeMask]),
                  'depth': ('depth',self.centers['depth']),
                  'latitude' : ('latitude', self.centers['latitude']),
                  'longitude': ('longitude',self.centers['longitude']),
                  }
        
        ds = xr.Dataset(None,coords=coords)
        
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

        ds.longitude.attrs = lon_attributtes
        ds.latitude.attrs = lat_attributtes
        ds.depth.attrs = depth_attributtes
        ds.time.attrs = time_atts
        ds.to_netcdf(self.netcdf_output_file)
        ds.close()
        return

    
    def counts(self,source = None):
        # counts 2d and 3d are splitted in two functions. 
        nz,ny,nx = [np.size(self.centers[key]) for key in ['depth','latitude','longitude']]
        nt = self.time[self.timeMask].size
        
        counts_t = np.zeros((nt,nz,ny,nx))
        i = 0
        t = 0
        for vtu_step in self.pvd_data.vtu_data:
            if self.timeMask[t] == True:
                position = vtu_step.points()
                r = np.c_[position['depth'], position['latitude'], position['longitude']]
                if source:
                    source_mask = vtu_step.points()['source'] == int(source)
                    r = r[source_mask]
                bins = [self.grid['depth'],self.grid['latitude'],self.grid['longitude']]
                counts_t[i], _ = np.histogramdd(r,bins=bins) 
                i=i+1
            t=t+1
#        else:
#            counts_t = np.zeros((nt,ny,nx))
#            i = 0
#            t = 0
#            for vtu_step in self.pvd_data.vtu_data:
#                if self.timeMask[t] == True:
#                    position = vtu_step.points()
#                    r = np.c_[position['latitude'],position['longitude']]
#                    if source:
#                        source_mask = vtu_step.points()['source'] == int(source)
#                        r = r[source_mask]
#                    bins = [self.grid['latitude'],self.grid['longitude']]
#                    counts_t[i], _ = np.histogramdd(r,bins=bins) 
#                    i=i+1
#                t=t+1        
        return counts_t
    
                
    def writeResidence_time(self):        
        print('--> Computing residence time on grid')        
        # Read netcdf, compute concentrations, compute residence tiem
        # compute global 
        ds = xr.open_dataset(self.netcdf_output_file)        
        counts_t = self.counts()
        self.residence_time = np.zeros(counts_t.shape[1:])
        for t_s in range(0,counts_t.shape[0]):
            self.residence_time = (counts_t[t_s]>0) * self.dt + self.residence_time            
        ds['residence_time_global'] = (['depth','latitude','longitude'], self.residence_time)
        ds.close()
        ds.to_netcdf(self.netcdf_output_file,'a')           
        
        # compute per source        
        for source in self.sources['id'].keys():  
            ds = xr.open_dataset(self.netcdf_output_file) 
            counts_t = self.counts(source=source)
            self.residence_time = np.zeros(counts_t.shape[1:])
            for t_s in range(0,counts_t.shape[0]):
                self.residence_time = (counts_t[t_s] > 0) * self.dt + self.residence_time                   
            var_name ='residence_time_source_' + source.zfill(3)
            ds[var_name] = (self.dims[1:],self.residence_time)
            ds[var_name].attrs = {'long_name':'residence_time', 'units':'s'}
            ds.close()
            ds.to_netcdf(self.netcdf_output_file,'a')        
        return
 
    
    def writeConcentrations(self):        
        print('--> Computing concentrations on grid')
        ds = xr.open_dataset(self.netcdf_output_file)        
        # Compute concentrations: total number of particles
        counts_t = self.counts()        
        conc_area = counts_t.sum(axis=1)/self.area
        conc_volume = counts_t/self.volume
        ds['concentration_area'] = (['time','latitude','longitude'], conc_area)
        ds['concentration_volume'] = (self.dims, conc_volume)
        ds.close()
        ds.to_netcdf(self.netcdf_output_file,'a')   
        # Compute concentrations: particles of each source        
        for source in self.sources['id'].keys():
            ds = xr.open_dataset(self.netcdf_output_file)
            counts_t = self.counts(source=source)
            conc_area = counts_t.sum(axis=1)/self.area
            conc_volume = counts_t/self.volume        
            var_name_as = 'concentration_area_source_' + source.zfill(3)
            ds[var_name_as] = (['time','latitude','longitude'], conc_area)
            ds[var_name_as].attrs = {'long_name':'concentration',
                                      'units':'ppm*m'}            
            var_name_vs = 'concentration_volume_source_' + source.zfill(3)
            ds[var_name_vs] = (self.dims, conc_volume)
            ds[var_name_vs].attrs = {'long_name':'concentration',
                                      'units':'ppm*m*m'}        
            ds.close()
            ds.to_netcdf(self.netcdf_output_file,'a')        
        return

    
    def writeCount(self):         
        print('--> Sampling tracers on grid')
        ds = xr.open_dataset(self.netcdf_output_file)
        ds['number_of_tracers'] = (self.dims, self.counts())
        ds.close()
        ds.to_netcdf(self.netcdf_output_file,'a')
        
    
    def writeVolume(self):         
        print('--> Writing grid cell volumes')
        ds = xr.open_dataset(self.netcdf_output_file)
        ds['cell_volume'] = (self.dims[1:], self.volume)
        ds.close()
        ds.to_netcdf(self.netcdf_output_file,'a')

    def write(self,variable):
        print('--> Writing '+ variable)
        ds = xr.open_dataset(self.netcdf_output_file)
        ds[variable] = (getattr(self,variable)['dims'], getattr(self,variable))
        ds.close()
        ds.to_netcdf(self.netcdf_output_file,'a')
        
#    def age(self):#        
#        print('--> Computing age on grid')       
#        ds = xr.open_dataset(self.netcdf_output_file)
#        # Compute concentrations: total number of particles
#        nz,ny,nx = [np.size(self.centers[key]) for key in ['depth','latitude','longitude']]
#        nt = len(self.pvd_data.vtu_data)
#        age_t = np.zeros((nt,nz,ny,nx))#        
#        grid = (self.grid['depth'],self.grid['latitude'],self.grid['longitude'])        
#        k=0
#        for vtu_step in self.pvd_data.vtu_data:
#                position = vtu_step.points()
#                r = np.c_[position['depth'], position['latitude'], position['longitude']]
#                age = position['age']
#                age_t[k] = griddata(r, age, grid, method='linear')
#                k = k + 1#        
#        var_name = 'age_global'
#        ds[var_name] = (self.dims, age_t)
#        ds[var_name].attrs = {'long_name':'age', 'units':'s'}
#        # COMPUTE AGE PER SOURCE
#        for source in self.sources['id'].keys():#            
#            k=0
#            for vtu_step in self.pvd_data.vtu_data:
#                    position = vtu_step.points()
#                    source_mask = vtu_step.points()['source'] == int(source)
#                    r = np.c_[position['depth'], position['latitude'], position['longitude']]
#                    r = r[source_mask]
#                    age = position['age'][source_mask]
#                    age_t[k] = griddata(r, age, grid, method='linear')
#                    k = k + 1#          
#            var_name ='age_source_' + source.zfill(3)
#            ds[var_name] = (self.dims[1:],age_t)
#            ds[var_name].attrs = {'long_name':'age', 'units':'s'}        
#        ds.close()
#        ds.to_netcdf(self.netcdf_output_file,'a')       
#        return
        

    def run_postprocessing(self, outDir, measures):
        self.get_pvd(outDir)
        self.get_time_axis()
        self.get_grid()
        self.get_sources()
        self.get_areas()
        self.get_netcdf_header()
        #writting fields
        self.writeCount()
        self.writeVolume()
        #self.write('volume')
        if 'concentrations' in measures:
            self.writeConcentrations()
        if 'residence_time' in measures:
            self.writeResidence_time()

        
def getRecipeListFromCase(xmlFile):
    recipeList =[]        
    #parsing case definition file
    root = ET.parse(xmlFile).getroot()                
    for fileName in root.findall('execution/postProcessing/file'):
        recipeList.append(fileName.get('name'))    
    return recipeList


def getFieldsFromRecipe(xmlFile):
    fieldList = []
    root = ET.parse(xmlFile).getroot()    
    for fieldName in root.findall('measures/field'):
        fieldList.append(fieldName.get('key'))    
    return fieldList


def main():
    lic = License()
    lic.print()
    
    #cmd line argument parsing
    argParser = argparse.ArgumentParser(description='Post processes MOHID Lagrangian outputs. Use -h for help.')
    argParser.add_argument("-i", "--input", dest="caseXML",
                    help=".xml file with the case definition for the MOHID Lagrangian run", metavar=".xml")
    argParser.add_argument("-f", "--force", dest="recipeXML",
                    help=".xml file with the recipe for a post process run - optional", metavar=".xml")
    argParser.add_argument("-o", "--outputDir", dest="outDir",
                    help="output directory", metavar="dir")
    args = argParser.parse_args()
    
    caseXML = getattr(args,'caseXML')
    recipeXML = []
    recipeXML.append(getattr(args,'recipeXML'))
    outDir = getattr(args,'outDir')
    
    print('-> Case definition file is ', caseXML)
    print('-> Main output directory is ', outDir)
    
    #get list of post cicles to run
    if recipeXML == [None]:
        #open caseXML and extract recipeXML names
        recipeXML = getRecipeListFromCase(caseXML)
        
    for recipe in recipeXML:
        outDirLocal = outDir + '/postProcess_' +os.path.basename(os_dir.filename_without_ext(recipe))
        post = GridBasedMeasures(caseXML, recipe, outDir, outDirLocal)
        measures = getFieldsFromRecipe(recipe)
        post.run_postprocessing(outDir, measures)

main()