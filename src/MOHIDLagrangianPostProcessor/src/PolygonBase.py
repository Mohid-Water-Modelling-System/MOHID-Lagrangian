# -*- coding: utf-8 -*-

from src.XMLReader import getBeachFromRecipe
from src.Polygon import Polygon
from tqdm import tqdm


class ConcentrationArea:

    def __init__(self, polygon):
        self.polygon = polygon
        self.base_name = 'concentration_area'
        self.long_name = 'concentration_area'
        self.units = 'particles/km^2'
        self.dims = ['time', 'index']
        self.coords = polygon.coords.items()

    def getMeasure(self, nCounts):
        data = nCounts/self.polygon.geoDataFrame.area
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

    def __init__(self, polygon):
        self.polygon = polygon
        self.base_name = 'n_counts'
        self.long_name = 'n_counts'
        self.units = 'particles-per-polygon'
        self.dims = ['time', 'index']
        self.coords = polygon.coords.items()

    def getMeasure(self, nCounts):
        data = nCounts
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


class PolygonBase:

    def __init__(self, xml_file, xml_recipe):
        self.polygon = Polygon(xml_file, xml_recipe)
        self.gridBasicMeasures = ['residence_time', 'concentrations']
        self.beachCondition = getBeachFromRecipe(xml_recipe)

    def run(self, measures, sources, vtuParser, fileTimeHandler, netcdfWriter):
        print('-> Measures to compute: ', measures)

        sourcesDict = sources['id']
        for sourceID, sourceName in sourcesDict.items():
            print('\t' + sourceName)

        # initialize the measures:
        RawCounter = RawCounts(self.polygon)
        if 'concentrations' in measures:
            AreaCounter = ConcentrationArea(self.polygon)

        # Initialize the measures:
        timeIdx = 0
        for vtuFile in tqdm(vtuParser.fileList, desc='Global'):

            sourceIdx = 0
            vtuParser.updateReaderWithFile(vtuFile)
            for sourceID, sourceName in tqdm(sourcesDict.items(), 'Source'):

                # get particle position
                particlePos = vtuParser.getVariableData('coords', sourceID, beachCondition=self.beachCondition)

                geoseries = Polygon.array_to_geoseries(particlePos)
                points = Polygon.geoseries_to_geodataframe(geoseries)

                # get raw counts
                nCountsArray = self.polygon.getCountsInPolygon(points)

                dataArrayDict, sourceName = RawCounter.run(nCountsArray, sourceName)
                netcdfWriter.appendVariableTimeStepToDataset(sourceName, dataArrayDict, timeIdx)

                if 'concentrations' in measures:
                    dataArrayDict, sourceName = AreaCounter.run(nCountsArray, sourceName)
                    netcdfWriter.appendVariableTimeStepToDataset(sourceName, dataArrayDict, timeIdx)

                sourceIdx += 1
            timeIdx += 1
        progress = '-> Progress: %4.2f' %(100*(timeIdx/len(vtuParser.fileList)))
        print(progress, end='\r')
