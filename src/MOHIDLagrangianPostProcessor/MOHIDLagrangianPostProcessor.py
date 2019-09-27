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

import xml.etree.ElementTree as ET
import numpy as np
import xarray as xr
import os
import sys
import argparse
from numba import jit
#from scipy.interpolate import griddata

# This environment variable avoids error on locking file when writing.
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

basePath = os.path.dirname(os.path.realpath(__file__))
commonPath = os.path.abspath(os.path.join(basePath, "../Common"))
sys.path.append(commonPath)

import os_dir
import MDateTime

from about import License
import vtuParser
import MOHIDLagrangianVTUtoHDF5
             
 
class PVDParser:
    
    def __init__(self,pvd_file):
        self.pvd_file = pvd_file
        self.vtu_list = []
        self.files = []
        self.timesteps = []
        self.vtu_data = []
    
    def get_vtu_files(self, outDir):
        self.vtu_list = vtuParser.validVtuFilesList(outDir)
        for vtu_file in self.vtu_list:
            self.files.append(vtu_file)
            self.vtu_data.append(vtuParser.VTUParser(vtu_file))
    

class GridBasedMeasures:
    def __init__(self,xml_file, xml_recipe, outdir, outdirLocal):
        self.xml_file = xml_file
        self.xml_recipe = xml_recipe
        self.pvd_file = outdir +'/'+self.xml_file.replace('.xml','.pvd')
        self.pvd_data = PVDParser(self.pvd_file)
        self.nFiles = len(vtuParser.validVtuFilesList(outdir))
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
        self.beach = []
        self.netcdf_output_file = outdirLocal+'/'+os.path.basename(self.xml_file.replace('.xml','.nc'))
        if os.path.exists(outdirLocal):
            os_dir.deleteDirForce(outdirLocal)
        os.mkdir(outdirLocal)
        
        
    def get_pvd(self, outDir):
        self.pvd_data.get_vtu_files(outDir)
    
    def get_beach(self):
        root = ET.parse(self.xml_recipe).getroot()
        for parameter in root.findall('EulerianMeasures/measures/filters/filter'):
            if parameter.get('key') == 'beaching':
                self.beach = parameter.get('value')
            else:
                self.beach = '0'
                    
    
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
            name_source = source.find('setsource').attrib['name']
            dict_source[id_source] = name_source
        dict_source['global'] = 'global'
        self.sources['id'] = dict_source
            # at the moment we don't require specific information just the ids.
            # In a future releases if we are going to compute measures using 
            # integration time, we will 
    
        
    def get_grid(self):
        root = ET.parse(self.xml_recipe).getroot()
        for parameter in root.findall('EulerianMeasures/gridDefinition/'):
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
            
            self.grid ={'longitude': np.linspace(x_min,x_max,np.int(x_step+1)),
                        'latitude': np.linspace(y_min,y_max,np.int(y_step+1)),
                        'depth': np.linspace(z_min,z_max,np.int(z_step+1))}
    
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

    
    def counts(self,source = 'global'):
        # counts 2d and 3d are splitted in two functions. 
        nz,ny,nx = [np.size(self.centers[key]) for key in ['depth','latitude','longitude']]
        nt = self.time[self.timeMask].size
        
        counts_t = np.zeros((nt,nz,ny,nx))
        i = 0
        t = 0
        bins = [self.grid['depth'],self.grid['latitude'],self.grid['longitude']]
        for vtu_step in self.pvd_data.vtu_data:
            if self.timeMask[t] == True:
                r = vtu_step.points('coords')
                if (source != 'global') and (self.beach == '0'):
                    source_mask = vtu_step.points('source') == int(source)
                    r=r[source_mask]
                elif (source != 'global') and (self.beach == '1'):
                    source_mask = vtu_step.points('source') == int(source)
                    beach_mask = vtu_step.points('state') < 0.5
                    r=r[source_mask*beach_mask]
                elif (source == 'global') and (self.beach == '1'):
                    beach_mask = vtu_step.points('state') < 0.5
                    r = r[beach_mask]
                elif (source != 'global') and (self.beach == '2'):
                    source_mask = vtu_step.points('source') == int(source)
                    beach_mask = vtu_step.points('state') >= 0.5
                    r=r[source_mask*beach_mask]
                elif (source == 'global') and (self.beach == '2'):
                    beach_mask = vtu_step.points('state') >= 0.5
                    r = r[beach_mask]
                counts_t[i], _ = np.histogramdd(r,bins=bins) 
                i=i+1
            t=t+1
        counts_t[:,np.all(counts_t == 0.,axis=0)] = np.nan
        return counts_t
    
    def writeVars(self, measures):
        
        exclusions = ['residence_time', 'concentrations']
        
        for source in self.sources['id'].keys():
            counts_t = self.counts(source=source)
            if 'residence_time' in measures: self.writeResidence_time(counts_t, source)
            if 'concentrations' in measures: self.writeConcentrations(counts_t, source)

            for measure in measures:
                if measure not in exclusions:
                    r_idx=self.var_to_grid(measure, source)
                    self.writeVar(r_idx, measure, source)
                 
                
    def var_to_grid(self,varname,source):
        print('--> Computing '+ varname +' on grid for source:' + self.sources['id'][str(source)])
        nz,ny,nx = [np.size(self.centers[key]) for key in ['depth','latitude','longitude']]
        nt = self.time[self.timeMask].size
        counts_t = np.zeros((nt,nz,ny,nx))
        
        i = 0
        t = 0
        for vtu_step in self.pvd_data.vtu_data:
            if self.timeMask[t] == True:
                r = vtu_step.points('coords')
                var = vtu_step.points(varname)
                if (source != 'global') and (self.beach == '0'):
                    source_mask = vtu_step.points('source') == int(source)
                    r=r[source_mask]
                    var = var[source_mask*beach_mask]
                if (source != 'global') and (self.beach == '1'):
                    source_mask = vtu_step.points('source') == int(source)
                    beach_mask = vtu_step.points('state') < 0.5
                    r=r[source_mask*beach_mask]
                    var = var[source_mask*beach_mask]
                elif (source == 'global') and (self.beach == '1'):
                    beach_mask = vtu_step.points('state') < 0.5
                    r = r[beach_mask]
                    var = var[beach_mask]
                elif (source != 'global') and (self.beach == '2'):
                    source_mask = vtu_step.points('source') == int(source)
                    beach_mask = vtu_step.points('state') >= 0.5
                    r=r[source_mask*beach_mask]
                    var = var[source_mask*beach_mask]
                elif (source == 'global') and (self.beach == '2'):
                    beach_mask = vtu_step.points('state') >= 0.5
                    r = r[beach_mask]
                    var = var[beach_mask]

                
                # Find the indices of the box in each dimension
                z_dig = np.int32((r[:,0] - min((self.grid['depth'])))/abs(self.grid['depth'][1]-self.grid['depth'][0]))
                z_dig[z_dig >= (nz-1)] = nz-1 
                z_dig[0 > z_dig] = 0 
                y_dig = np.int32((r[:,1] - min((self.grid['latitude'])))/abs(self.grid['latitude'][1]-self.grid['latitude'][0]))
                y_dig[y_dig >= (ny-1)] = ny-1 
                y_dig[0 > y_dig] = 0
                x_dig = np.int32((r[:,2] - min((self.grid['longitude'])))/abs(self.grid['longitude'][1]-self.grid['longitude'][0]))
                x_dig[x_dig >= (nx-1)] = nx-1 
                x_dig[0 > x_dig] = 0 
                
                r_idx = np.ravel_multi_index((z_dig,y_dig,x_dig),(nz,ny,nx))
                N = nx*ny*nz    
                counts_t[i]=np.reshape(cellCounting(r_idx,var,N),(nz,ny,nx))
                
                i = i+1
            t = t+1
        counts_t[:,np.all(counts_t == 0.,axis=0)] = np.nan                  
        return counts_t
        

    def writeResidence_time(self,counts_t,source):        
        print('--> Computing residence time on grid for source:' + self.sources['id'][str(source)])        
        ds = xr.open_dataset(self.netcdf_output_file) 
        self.residence_time = np.zeros(counts_t.shape[1:])
        for t_s in range(0,counts_t.shape[0]):
            self.residence_time = (counts_t[t_s] > 0) * self.dt + self.residence_time
        self.residence_time[self.residence_time == 0] = np.nan                   
        var_name ='residence_time_' + self.sources['id'][str(source)]
        ds[var_name] = (self.dims[1:],self.residence_time)
        ds[var_name].attrs = {'long_name':'residence_time', 'units':'s'}
        ds.close()
        ds.to_netcdf(self.netcdf_output_file,'a')
        return
 
    
    def writeConcentrations(self,counts_t,source):        
        print('--> Computing concentrations on grid for source:' + self.sources['id'][str(source)])
        ds = xr.open_dataset(self.netcdf_output_file)
        conc_area = counts_t.sum(axis=1)/self.area
        conc_volume = counts_t/self.volume        
        var_name_as = 'concentration_area_' + self.sources['id'][str(source)]
        ds[var_name_as] = (['time','latitude','longitude'], conc_area)
        ds[var_name_as].attrs = {'long_name':'concentration',
                                  'units':'pp/m^2'}            
        var_name_vs = 'concentration_volume_' + self.sources['id'][str(source)]
        ds[var_name_vs] = (self.dims, conc_volume)
        ds[var_name_vs].attrs = {'long_name':'concentration',
                                  'units':'pp/m^3'}        
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


    def writeVar(self,vardata, varname, source):
        print('--> Writing '+ varname + ' for source:',self.sources['id'][str(source)])
        ds = xr.open_dataset(self.netcdf_output_file)
        ds[varname] = (self.dims, vardata)
        ds.close()
        ds.to_netcdf(self.netcdf_output_file,'a')
    
    
    def checkNc(self):
        ds = xr.open_dataset(self.netcdf_output_file)
        if ds.depth.size == 1:
            print('->The nc has a depth degenerated dimension: Squeezing...' )
            ds = ds.squeeze(dim='depth',drop=True)
            # When you alter the dimensions of the netcdf, you cannot overwrite the
            # netcdf. We ha create a new one without extension, remove the old one,
            # and rename the first one.
            nc_squeezed = self.netcdf_output_file.replace('.nc','_sq.nc')
            ds.to_netcdf(nc_squeezed)
            ds.close()
            try:
                #os.remove(self.netcdf_output_file)
                os.replace(nc_squeezed, nc_squeezed.replace('_sq.nc','.nc'))
            except:
                print('The file cannot be deleted')
                pass
                   

    def run_postprocessing(self, outDir, measures):
        self.get_pvd(outDir)
        self.get_time_axis()
        self.get_grid()
        self.get_sources()
        self.get_beach()
        self.get_areas()
        self.get_netcdf_header()
        #writting fields
        self.writeCount()
        self.writeVolume()
        self.writeVars(measures)
        self.checkNc()

@jit
def cellCounting(r_idx,var,N):
    gridcounts = np.zeros(N)
    for k in range(0,N):
        var_selected = var[k==r_idx]
        if var_selected.size == 0:
            gridcounts[k] = 0
        else:
           gridcounts[k] =np.sum(var_selected)/var_selected.size
    return gridcounts

       
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
    for fieldName in root.findall('EulerianMeasures/measures/field'):
        fieldList.append(fieldName.get('key'))    
    return fieldList


def checkHDF5WriteRecipe(xmlFile):
    convert = False
    root = ET.parse(xmlFile).getroot()    
    for formats in root.findall('convertFiles/format'):
        if formats.get('key') == 'hdf5':
            convert = True
    return convert


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
    
    print('-> Case definition file is', caseXML)
    print('-> Main output directory is', outDir)
    
    #get list of post cicles to run
    if recipeXML == [None]:
        #open caseXML and extract recipeXML names
        recipeXML = getRecipeListFromCase(caseXML)
        
    for recipe in recipeXML:
        print('-> Running recipe', recipe)
        outDirLocal = outDir + '/postProcess_' +os.path.basename(os_dir.filename_without_ext(recipe))
        post = GridBasedMeasures(caseXML, recipe, outDir, outDirLocal)
        measures = getFieldsFromRecipe(recipe)
        post.run_postprocessing(outDir, measures)
        if checkHDF5WriteRecipe(recipe):
            MOHIDLagrangianVTUtoHDF5.run(caseXML, recipe, outDir)

main()