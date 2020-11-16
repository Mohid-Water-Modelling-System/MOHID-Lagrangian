# -*- coding: utf-8 -*-
""" Module compute the measures using the grid"""

from src.XMLReader import *
from src.Grid import Grid
import numpy as np
from tqdm import tqdm


def is2D(array):
    """Check if an array is 2D."""
    return array.ndim == 2


def is2Dlayer(array):
    """"Check if an array is 3D array with shape [n_i,n_j,1]"""
    return (array.ndim == 3) and (array.shape[0] == 1)


class MeasuresBase:
    """Parent class to build child measure classes"""

    def __init__(self):
        self.grid = []
        self.base_name = []
        self.long_name = []
        self.units = []
        self.dims = []
        self.coords = []

    def getMeasure(self, nCounts):
        """Measure to """
        pass

    def addSourceName(self, source_name):
        return self.long_name + '_' + source_name

    def toDataArrayDict(self, data):
        d = {}
        d['coords'] = {k: v for k, v in self.coords if k in self.dims}
        d['dims'] = self.dims
        d['data'] = data
        d['attrs'] = {'units': self.units,
                      'long_name': self.long_name}
        return d

    def run(self, nCounts, source_name):
        data = self.getMeasure(nCounts)
        var_name = self.addSourceName(source_name)
        var_dict = self.toDataArrayDict(data)
        return var_dict, var_name


class ConcentrationArea(MeasuresBase):

    def __init__(self, grid):
        MeasuresBase.__init__(self)
        self.grid = grid
        self.base_name = 'concentration'
        self.long_name = 'concentration_area'
        self.units = 'particles/km^2'
        self.dims = ['time', 'depth', 'latitude', 'longitude']
        self.coords = grid.coords.items()

    def getMeasure(self, nCounts):
        data = nCounts/self.grid.cellArea

        if is2Dlayer(data):
            self.dims = ['time', 'latitude', 'longitude']
            data = np.squeeze(data, axis=0)

        return data


class ConcentrationVolume(MeasuresBase):

    def __init__(self, grid):
        MeasuresBase.__init__(self)
        self.grid = grid
        self.base_name = 'concentration'
        self.long_name = 'concentration_volume'
        self.units = 'particles/km^3'
        self.dims = ['time', 'depth', 'latitude', 'longitude']
        self.coords = grid.coords.items()

    def getMeasure(self, nCounts):
        data = nCounts/self.grid.cellArea
        if is2Dlayer(data):
            self.dims = ['time', 'latitude', 'longitude']
            data = np.squeeze(data, axis=0)
        return data


class RawCounts(MeasuresBase):

    def __init__(self, grid):
        MeasuresBase.__init__(self)
        self.grid = grid
        self.base_name = 'n_counts'
        self.long_name = 'n_counts'
        self.units = 'particles-per-box'
        self.dims = ['time', 'depth', 'latitude', 'longitude']
        self.coords = grid.coords.items()

    def getMeasure(self, nCounts):
        data = nCounts
        if is2Dlayer(data):
            self.dims = ['time', 'latitude', 'longitude']
            data = np.squeeze(data, axis=0)
        return data


class ResidenceTime(MeasuresBase):

    def __init__(self, grid, dt, base_name='residence_time', units ='s'):
        MeasuresBase.__init__(self)
        self.grid = grid
        self.base_name = base_name
        self.long_name = base_name
        self.units = units
        self.dims = ['time', 'depth', 'latitude', 'longitude']
        self.coords = grid.coords.items()
        self.accum = np.array([])
        self.dt = np.array(dt)
        self.tidx = 0

    def increase_tidx(self):
        self.tidx += 1

    def getMeasure(self, nCounts):
        data = (nCounts > 0)*self.dt
        if self.tidx == 0:
            self.accum = np.zeros_like(nCounts)
        else:
            self.accum = data + self.accum

        if is2Dlayer(data):
            self.dims = ['time', 'latitude', 'longitude']
            data = np.squeeze(data, axis=0)

        # Increase the counter
        self.increase_tidx()

        return self.accum


class VarInCell(MeasuresBase):

    def __init__(self, grid, base_name='base_name', units='units'):
        MeasuresBase.__init__(self)
        self.grid = grid
        self.base_name = base_name
        self.long_name = base_name
        self.units = units
        self.dims = ['time', 'depth', 'latitude', 'longitude']
        self.coords = grid.coords.items()

    def getMeasure(self, varInCell):
        if is2Dlayer(varInCell):
            self.dims = ['time', 'latitude', 'longitude']
            data = np.squeeze(varInCell, axis=0)
        return data



class GridBase:

    def __init__(self, xml_file, xml_recipe):
        self.grid = Grid(xml_recipe, xml_file)
        self.grid.initializeGrid()
        self.gridBasicMeasures = ['residence_time', 'concentrations']
        self.beachCondition = getBeachFromRecipe(xml_recipe)

    def run(self, measures, sources, vtuParser, fileTimeHandler, netcdfWriter):
        print('-> Measures to compute: ', measures)

        sourcesDict = sources['id']
        for sourceID, sourceName in sourcesDict.items():
            print('\t' + sourceName)

        dt = fileTimeHandler.dt

        # initialize the measures:
        RawCounter = RawCounts(self.grid)
        if 'concentrations' in measures:
            AreaCounter = ConcentrationArea(self.grid)
            VolumeCounter = ConcentrationVolume(self.grid)
        if 'residence_time' in measures:
            ResidenceCounter = ResidenceTime(self.grid, dt)

        VarInCellCounter = {}
        for measure in measures:
            if measure not in self.gridBasicMeasures:
                VarInCellCounter[measure] = VarInCell(self.grid, base_name=measure)

        netcdfWriter.resetTimeIdx()
        vtuFileList = vtuParser.fileList
        for vtuFile in tqdm(vtuFileList, desc='Progress', position=0, leave=True):
            sourceIdx = 0
            vtuParser.updateReaderWithFile(vtuFile)

            for sourceID, sourceName in sourcesDict.items():

                # get particle position
                particlePos = vtuParser.getVariableData('coords', sourceID, beachCondition=self.beachCondition)
                # get counts
                nCountsArray = self.grid.getCountsInCell(particlePos)
                # Update valid cells (to avoid to compute measures on empty cells)
                self.grid.setValidCells(nCountsArray)

                dataArrayDict, dataArrayName = RawCounter.run(nCountsArray, sourceName)
                netcdfWriter.appendVariableTimeStepToDataset(dataArrayName, dataArrayDict)

                if 'residence_time' in measures:
                    dataArrayDict, dataArrayName = ResidenceCounter.run(nCountsArray, sourceName)
                    netcdfWriter.appendVariableTimeStepToDataset(dataArrayName, dataArrayDict)

                if 'concentrations' in measures:
                    dataArrayDict, dataArrayName = AreaCounter.run(nCountsArray, sourceName)
                    netcdfWriter.appendVariableTimeStepToDataset(dataArrayName, dataArrayDict)
                    dataArrayDict, dataArrayName = VolumeCounter.run(nCountsArray, sourceName)

                    netcdfWriter.appendVariableTimeStepToDataset(dataArrayName, dataArrayDict)

                for measure in VarInCellCounter:
                    if measure not in self.gridBasicMeasures:
                        varInParticles = vtuParser.getVariableData(measure, sourceID, beachCondition = self.beachCondition)
                        varInCell = self.grid.getMeanDataInCell(particlePos, varInParticles)
                        dataArrayDict, dataArrayName = VarInCellCounter[measure].run(varInCell, sourceName)
                        netcdfWriter.appendVariableTimeStepToDataset(dataArrayName, dataArrayDict)

                sourceIdx += 1
            netcdfWriter.increaseTimeIdx()