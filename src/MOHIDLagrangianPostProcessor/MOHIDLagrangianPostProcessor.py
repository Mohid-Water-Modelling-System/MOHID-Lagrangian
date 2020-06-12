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

import os
import sys
import argparse

# This environment variable avoids error on locking file when writing.
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

basePath = os.path.dirname(os.path.realpath(__file__))
commonPath = os.path.abspath(os.path.join(basePath, "../Common"))
sys.path.append(commonPath)

import os_dir
import MDateTime
import numpy as np
from about import License
from src.VTUtoHDF5 import vtu2hdf5
from src.PVDParser import PVDParser
from src.XMLReader import *
from src.RectangularGrid import RectangularGridBase
from src.Time import FilesTimesHandler 
from src.NcWriter import NetcdfParser
from src.MeasuresGrid import getResidenceTime, getCountsInCell, getConcentrationsArea,getConcentrationsVolume,getVariableMeanCell
from src.PlotResults import *


class MOHIDLagrangianPostProcessorBase:

    def __init__(self, xml_file, xml_recipe, outdir, outdirLocal):
        self.xml_file = xml_file
        self.xml_recipe = xml_recipe
        self.pvd_file = outdir + '/' + self.xml_file.replace('.xml', '.pvd')
        self.pvd_data = PVDParser(self.pvd_file)
        self.beachCondition = getBeachFromRecipe(xml_recipe)
        self.sources = getSourcesDictFromXML(xml_file)
        self.outdir = outdir
        self.outdirLocal = outdirLocal
        self.time = []
        if os.path.exists(outdirLocal):
            os_dir.deleteDirForce(outdirLocal)
        os.mkdir(outdirLocal)

    def getPVDReader(self):
        self.pvdReader = PVDParser(self.pvd_file)
        self.pvdReader.getVtuFileHandlers(self.outdir)

    def getTimeFileHandler(self):
        self.FileTimeHandler = FilesTimesHandler(self.pvdReader.vtuFilelist)
        self.FileTimeHandler.initializeTimeGrid(self.xml_file, self.xml_recipe)

    def getMOHIDLagrangianBase(self):
        self.getPVDReader()
        self.getTimeFileHandler()
        self.pvdReader.updateVtuFileHandlers(self.FileTimeHandler.cropFileList())


class MOHIDLagrangianGridBasedMeasures:

    def __init__(self):
        self.base = []
        self.grid = []
        self.netcdfWriter = []
        self.outputFile = []
        self.gridBasicMeasures = ['residence_time', 'concentrations']

    def initialize(self,xml_file, xml_recipe, outdir, outdirLocal):
        self.base = MOHIDLagrangianPostProcessorBase(xml_file, xml_recipe, outdir, outdirLocal)
        self.base.getMOHIDLagrangianBase()

        self.outputFile = outdirLocal + self.base.xml_file.replace('.xml', '.nc')

        self.grid = RectangularGridBase(xml_recipe, xml_file)
        self.grid.initializeGrid()

        self.netcdfWriter = NetcdfParser(self.outputFile)
        self.netcdfWriter.initDataset(self.grid, self.base.FileTimeHandler)


    def run(self, measures):
        print('-> Measures to compute: ', measures)

        print('-> Sources: ')
        for i in self.base.sources['id']:
            print('\t', self.base.sources['id'][i])
        t = 0
        for vtu_step in self.base.pvdReader.vtuFileHandlers:
            s = 0
            for source in self.base.sources['id'].keys():
                particlePos = vtu_step.getVtuVariableData('coords', source, beachCondition=self.base.beachCondition) 
                self.grid.getCountsInCell(particlePos)
                nCounts = getCountsInCell(self.grid)
                #
                self.netcdfWriter.appendVariableTimeStepToDataset('n_counts_'+self.base.sources['id'][source],nCounts,t)

                if 'residence_time' in measures:
                    ResidenceTime = getResidenceTime(self.grid, self.base.FileTimeHandler.dt)
                    if t == 0:
                        ResidenceTimeAccum = len(self.base.sources['id'].keys())*[np.zeros_like(ResidenceTime['data'])]
                    ResidenceTimeAccum[s] = ResidenceTime['data'] + ResidenceTimeAccum[s]
                    ResidenceTime['data'] = ResidenceTimeAccum[s]
                    self.netcdfWriter.appendVariableTimeStepToDataset('residence_time_'+self.base.sources['id'][source],ResidenceTime,t)

                if 'concentrations' in measures:
                    ConcentrationsArea = getConcentrationsArea(self.grid)
                    ConcentrationsVolume = getConcentrationsVolume(self.grid)
                    self.netcdfWriter.appendVariableTimeStepToDataset('concentration_volume_'+self.base.sources['id'][source],ConcentrationsArea,t)
                    self.netcdfWriter.appendVariableTimeStepToDataset('concentration_area'+self.base.sources['id'][source],ConcentrationsVolume,t)

                for measure in measures:
                    if measure not in self.gridBasicMeasures:
                        varInParticles = vtu_step.getVtuVariableData(measure, source, beachCondition=self.base.beachCondition)
                        self.grid.getMeanDataInCell(varInParticles)
                        varInCell = getVariableMeanCell(self.grid, measure)
                        self.netcdfWriter.appendVariableTimeStepToDataset(measure+'_'+self.base.sources['id'][source],varInCell,t)
                s = s+1
            t = t + 1
            progress = '-> Progress: %4.2f' %(100*(t/len(self.base.pvdReader.vtuFileHandlers)))
            print(progress, end='\r')


def runPostprocessing(caseXML, recipe, outDir, outDirLocal, measures):
    postProcessor = MOHIDLagrangianGridBasedMeasures()
    postProcessor.initialize(caseXML, recipe, outDir, outDirLocal)
    postProcessor.run(measures)
    return postProcessor


def main():
    lic = License()
    lic.print()

    # cmd line argument parsing
    argParser = argparse.ArgumentParser(description='Post processes MOHID Lagrangian outputs. Use -h for help.')
    argParser.add_argument("-i", "--input", dest="caseXML",
                    help=".xml file with the case definition for the MOHID Lagrangian run", metavar=".xml")
    argParser.add_argument("-f", "--force", dest="recipeXML",
                    help=".xml file with the recipe for a post process run - optional", metavar=".xml")
    argParser.add_argument("-o", "--outputDir", dest="outDir",
                    help="output directory", metavar="dir")
    argParser.add_argument("-po", "--plot-only", dest="plotonly",
                    help="output directory", action='store_true')
    args = argParser.parse_args()


    caseXML = getattr(args, 'caseXML')
    recipeXML = []
    recipeXML.append(getattr(args, 'recipeXML'))
    outDir = getattr(args, 'outDir')

    print('-> Case definition file is', caseXML)
    print('-> Main output directory is', outDir)

    # get list of post cicles to run
    if recipeXML == [None]:
        # open caseXML and extract recipeXML names
        recipeXML = getRecipeListFromCase(caseXML)

    for recipe in recipeXML:
        print('-> Running recipe', recipe)
        outDirLocal = outDir + '/postProcess_' + os.path.basename(os_dir.filename_without_ext(recipe)) + '/' 

        # If plotonly flag > Runs only plotting stage
        if args.plotonly == True:
             plotResultsFromRecipe(outDirLocal, recipe)
             return

        measures = getFieldsFromRecipe(recipe)
        postProcessorBase = runPostprocessing(caseXML, recipe, outDir, outDirLocal, measures)
        if checkHDF5WriteRecipe(recipe):
            vtu2hdf5(postProcessorBase, outDir)
        if checkPlotRecipe(recipe):
            plotResultsFromRecipe(outDirLocal, recipe)

main()
