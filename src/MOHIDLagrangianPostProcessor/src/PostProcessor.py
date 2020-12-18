# -*- coding: utf-8 -*-

from src.XMLReader import *
from src.VTUtoHDF5 import *
from src.VTUParser import VTUParser
from src.Time import FilesTimesHandler
from src.GridBase import GridBase
from src.NcWriter import NetcdfParser
from src.PolygonBase import PolygonBase
import os
import os_dir


class PostProcessor:

    def __init__(self, xml_file, xml_recipe, outdir, outdirLocal):
        self.xml_file = xml_file
        self.xml_recipe = xml_recipe
        self.pvd_file = outdir + '/' + self.xml_file.replace('.xml', '.pvd')
        self.beachCondition = getBeachFromRecipe(xml_recipe)
        self.sources = getSourcesDictFromXML(xml_file)
        self.outdir = outdir
        self.outdirLocal = outdirLocal
        self.time = []
        if os.path.exists(outdirLocal):
            os_dir.deleteDirForce(outdirLocal)
        os.mkdir(outdirLocal)

    def run(self):
        vtuParser = VTUParser(self.outdir)
        fileTimeHandler = FilesTimesHandler(vtuParser.fileList)
        fileTimeHandler.initializeTimeGrid(self.xml_file, self.xml_recipe)
        sliceTimeFileList = fileTimeHandler.cropFileList()
        vtuParser.updateFileList(sliceTimeFileList)

        outputFile = self.outdirLocal + self.xml_file.replace('.xml', '.nc')
        netcdfWriter = NetcdfParser(outputFile)

        measures = getFieldsFromRecipe(self.xml_recipe)
        sources = getSourcesDictFromXML(self.xml_file)

        if isPolygonOrGrid(self.xml_recipe) == 'grid':
            gridBase = GridBase(self.xml_file, self.xml_recipe)
            netcdfWriter.initDataset(gridBase.grid, fileTimeHandler)
            gridBase.run(measures, sources, vtuParser, fileTimeHandler, netcdfWriter)

        elif isPolygonOrGrid(self.xml_recipe) == 'polygon':
            polygonBase = PolygonBase(self.xml_file, self.xml_recipe)
            netcdfWriter.initDataset(polygonBase.polygon, fileTimeHandler)
            polygonBase.run(measures, sources, vtuParser, fileTimeHandler, netcdfWriter)

        if checkHDF5WriteRecipe(self.xml_recipe):
            vtu2hdf5(vtuParser, fileTimeHandler, self.outdirLocal)