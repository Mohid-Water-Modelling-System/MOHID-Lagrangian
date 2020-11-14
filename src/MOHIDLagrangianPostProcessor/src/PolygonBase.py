# -*- coding: utf-8 -*-

from src.XMLReader import getBeachFromRecipe
from src.Polygon import Polygon
from src.GridBase import MeasuresBase
from tqdm import tqdm


class ConcentrationArea(MeasuresBase):

    def __init__(self, polygon):
        MeasuresBase.__init__(self)
        self.polygon = polygon
        self.base_name = 'concentration_area'
        self.long_name = 'concentration_area'
        self.units = 'particles/km^2'
        self.dims = ['time', 'index']
        self.coords = polygon.coords.items()

    def getMeasure(self, nCounts):
        data = nCounts/self.polygon.geoDataFrame.area
        return data


class RawCounts(MeasuresBase):

    def __init__(self, polygon):
        MeasuresBase.__init__(self)
        self.polygon = polygon
        self.base_name = 'n_counts'
        self.long_name = 'n_counts'
        self.units = 'particles-per-polygon'
        self.dims = ['time', 'index']
        self.coords = polygon.coords.items()

    def getMeasure(self, nCounts):
        data = nCounts
        return data



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

        netcdfWriter.resetTimeIdx()
        for vtuFile in tqdm(vtuParser.fileList, desc='Progres'):

            sourceIdx = 0
            vtuParser.updateReaderWithFile(vtuFile)
            for sourceID, sourceName in sourcesDict.items():

                # get particle position
                particlePos = vtuParser.getVariableData('coords', sourceID, beachCondition=self.beachCondition)

                geoseries = Polygon.array_to_geoseries(particlePos)
                points = Polygon.geoseries_to_geodataframe(geoseries)

                # get raw counts
                nCountsArray = self.polygon.getCountsInPolygon(points)

                dataArrayDict, sourceName = RawCounter.run(nCountsArray, sourceName)
                netcdfWriter.appendVariableTimeStepToDataset(sourceName, dataArrayDict)

                if 'concentrations' in measures:
                    dataArrayDict, sourceName = AreaCounter.run(nCountsArray, sourceName)
                    netcdfWriter.appendVariableTimeStepToDataset(sourceName, dataArrayDict)

                sourceIdx += 1
            netcdfWriter.increaseTimeIdx()
