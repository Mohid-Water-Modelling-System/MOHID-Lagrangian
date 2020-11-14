# -*- coding: utf-8 -*-

import numpy as np
import geopandas as gpd
from geopandas import GeoSeries
from shapely.geometry import Point
from src.XMLReader import getPolygonFileFromRecipe


class Polygon:
    """ Class to create a polygon from a shapefile and perform counts.""" 

    def __init__(self, xml_file, xml_recipe):
        self.fileName = getPolygonFileFromRecipe(xml_recipe)[0]
        self.geoDataFrame = gpd.read_file(self.fileName).to_crs({"init": 'EPSG:4326'})
        self.dim = 'index'
        self.ids = np.arange(0, self.geoDataFrame.shape[0])
        self.coords = {self.dim: (self.dim, self.ids)}

    def updateCrsGeoDataFrame(self):
        """ Transform the projection of the shapefile to EPSG:4326"""
        if self.geoDataFrame.crs['init'] != 'EPSG:4326':
            self.geoDataFrame = self.geoDataFrame.to_crs({"init": 'EPSG:4326'})

    @staticmethod
    def array_to_geoseries(array: np.array) -> GeoSeries:
        """ Converts a numpy array x,y into a geoseries.

        Args:
            array (np.array): particle positions [z,y,x]

        Returns:
            GeoSeries: geoseries with position x,y.

        """
        x = array[:, 2]
        y = array[:, 1]

        return GeoSeries(map(Point, zip(x, y)))

    @staticmethod
    def geoseries_to_geodataframe(geoserie):
        """ Converts a geoseries into a geodataFrame.

        Args:
            geoserie (GeoSeries): Input Geoseries

        Returns:
            envgdf (gpd.GeoDataframe): Input Geodataframe

        """

        envgdf = gpd.GeoDataFrame(gpd.GeoSeries(geoserie))
        envgdf = envgdf.rename(columns={0: 'geometry'}).set_geometry('geometry')
        return envgdf

    def getCountsInPolygon(self, points):
        """ Counts the number of points inside each polygon.

            points(Geodataframe): Points position.

        """
        self.geoDataFrame['index'] = self.ids
        if points.shape[0] == 0:
            countsArray = np.zeros_like(self.ids)
        else:
            counts = gpd.sjoin(points, self.geoDataFrame, op='within')
            joins = counts.groupby('index_right').size()
            df = joins.to_frame().reset_index()
            df.columns = ['index', 'counts']
            counts = self.geoDataFrame.merge(df, on='index', how='outer')
            countsArray = counts['counts'].values
        return countsArray