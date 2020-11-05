# -*- coding: utf-8 -*-
import numpy as np

from src.XMLReader import *
from src.Grid import Grid
from tqdm import tqdm


def to_dict(dict_coords, dims, data, units, long_name):
    d = {}
    d['coords'] = {k: v for k, v in dict_coords if k in dims}
    d['dims'] = dims
    d['data'] = data
    d['attrs'] = {'units': units,
                  'long_name': long_name}
    return d


def is2D(array):
    return array.ndim == 2


def is2Dlayer(array):
    return (array.ndim == 3) and (array.shape[0] == 1)


class ConcentrationArea:

    def __init__(self, grid):
        self.grid = grid
        self.base_name = 'concentrations'
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

    def addSourceName(self, source_name):
        return self.base_name + '_' + source_name

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


class ConcentrationVolume:

    def __init__(self, grid):
        self.grid = grid
        self.base_name = 'concentrations'
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

    def addSourceName(self, source_name):
        return self.base_name + '_' + source_name

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


class RawCounts:

    def __init__(self, grid):
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

    def addSourceName(self, source_name):
        return self.base_name + '_' + source_name

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


class ResidenceTime:

    def __init__(self, grid, dt, base_name='', units =''):
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
            data = data + self.accum

        if is2Dlayer(data):
            self.dims = ['time', 'latitude', 'longitude']
            data = np.squeeze(data, axis=0)

        # Increase the counter
        self.increase_tidx()

        return data
    
    def addSourceName(self, source_name):
        return self.base_name + '_' + source_name

    def get_variable_name(self, source_name):
        return self.base_name + source_name

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


class VarInCell:

    def __init__(self, grid, base_name='base_name', units='units'):
        self.grid = grid
        self.base_name = base_name
        self.long_name = base_name
        self.units = units
        self.dims = ['time', 'depth', 'latitude', 'longitude']
        self.coords = grid.coords.items()

    def getMeasure(self, varInCell):
        data = varInCell
        if is2Dlayer(data):
            self.dims = ['time', 'latitude', 'longitude']
            data = np.squeeze(data, axis=0)
        return data

    def addSourceName(self, source_name):
        return self.base_name + '_' + source_name

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

        timeIdx = 0
        vtuFileList = vtuParser.fileList
        for vtuFile in tqdm(vtuFileList, desc='Global'):
            sourceIdx = 0
            vtuParser.updateReaderWithFile(vtuFile)
            for sourceID, sourceName in tqdm(sourcesDict.items(), 'Source'):

                # get particle position
                particlePos = vtuParser.getVariableData('coords', sourceID, beachCondition=self.beachCondition)

                # get raw counts
                nCountsArray = self.grid.getCountsInCell(particlePos)
                dataArrayDict, sourceName = RawCounter.run(nCountsArray, sourceName)
                netcdfWriter.appendVariableTimeStepToDataset(sourceName, dataArrayDict, timeIdx)

                if 'residence_time' in measures:
                    dataArrayDict, sourceName = ResidenceCounter.run(nCountsArray, sourceName)
                    netcdfWriter.appendVariableTimeStepToDataset(sourceName, dataArrayDict, timeIdx)
                if 'concentrations' in measures:
                    dataArrayDict, sourceName = AreaCounter.run(nCountsArray, sourceName)
                    netcdfWriter.appendVariableTimeStepToDataset(sourceName, dataArrayDict, timeIdx)

                    dataArrayDict, sourceName = VolumeCounter.run(nCountsArray, sourceName)
                    netcdfWriter.appendVariableTimeStepToDataset(sourceName, dataArrayDict, timeIdx)

                for measure in VarInCellCounter:
                    if measure not in self.gridBasicMeasures:
                        varInParticles = vtuParser.getVariableData(measure, sourceID, beachCondition = self.beachCondition)
                        self.grid.getMeanDataInCell(varInParticles)
                        varInCell = VarInCellCounter[measure].run(self.grid, measure=measure)
                        dataArrayDict, sourceName = VarInCellCounter[measure].run(nCountsArray)
                        netcdfWriter.appendVariableTimeStepToDataset(sourceName,dataArrayDict,timeIdx)

                sourceIdx += 1
            timeIdx += 1
            progress = '-> Progress: %4.2f' %(100*(timeIdx/len(vtuFileList)))
            print(progress, end='\r')