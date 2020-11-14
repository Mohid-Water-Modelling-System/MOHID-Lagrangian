# -*- coding: utf-8 -*-
"""
Module to create the grid and the cell based methods.
"""
import numpy as np
import xml.etree.ElementTree as ET
import src.constants as cte
from numba import jit


class GridJIT:
    """
    Class with grid functions using jit (just in time compiler)
    """

    @jit(nopython=True)
    def cellCounting(rIdCell: np.array, nCells: int) -> np.array:
        """
        Counts the number of particles in each cell identified by an id.

        Particle has an id based on the cell where the particle is.

        Args:
            rIdCell (np.array): Particle id cell position: [[i,j,k]]
            nCells (int): number of cells..

        Returns:
            cellCounts (np.array):counts per cell.

        """
        cellCounts = np.empty(nCells)
        for idCell in range(0, nCells):
            cellCounts[idCell] = np.sum(idCell == rIdCell)
        return cellCounts

    @jit(nopython=True)
    def cellMeanData(rIdCell: np.array, nCells: int, validCells: np.array,
                     varData: np.array) -> np.array:
        """
        Computes the mean value of the particle properties inside a cell.

        Args:
            rIdCell (np.array): Particle id cell position: [[i,j,k]]
            nCells (int): number of cells..
            validCells (np.array): Mask containing particles
            varData (np.array): Particle scalar value [[scalar]]

        Returns:
            cellMean (np.array): mean value of particle property in cell.

        """
        cellMean = np.zeros(nCells)
        for idCell in validCells:
            dataInCell = (idCell == rIdCell)*varData
            if dataInCell.size == 0:
                cellMean[idCell] = 0
            else:
                cellMean[idCell] = np.sum(dataInCell)/dataInCell.size
        return cellMean


class Grid:

    def __init__(self, xml_recipe, xml_file, dims=['depth', 'latitude', 'longitude']):
        self.xml_recipe = xml_recipe
        self.xml_file = xml_file
        self.grid = []
        self.cellCenters = []
        self.cellArea = []
        self.cellVolume = []
        self.dims = dims
        self.coords = {}
        self.countsInCell = []
        self.validCells = []
        self.meanDataInCell = []
        self.rIdCell = []
        self.nCells = []
        self.shape =[]

    def setGrid(self):
        root = ET.parse(self.xml_recipe).getroot()
        self.grid = len(self.dims)*[[]]
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

        if units_value == 'degrees':
            self.grid[2] = np.arange(x_min, x_max, x_step)
            self.grid[1] = np.arange(y_min, y_max, y_step)
            self.grid[0] = np.arange(z_min, z_max, z_step)
        elif units_value == 'relative':
            self.grid[2] = np.linspace(x_min, x_max, np.int(x_step + 1))
            self.grid[1] = np.linspace(y_min, y_max, np.int(y_step + 1))
            self.grid[0] = np.linspace(z_min, z_max, np.int(z_step + 1))
        elif units_value == 'meters':
            y_c = (y_max + y_min)/2.
            dlat = y_step/(cte.degreesToRad*cte.earthRadius)
            dlon = x_step/(cte.degreesToRad*cte.earthRadius * np.cos(cte.degreesToRad*(y_c)))
            self.grid[2] = np.arange(x_min, x_max, dlon)
            self.grid[1] = np.arange(y_min, y_max, dlat)
            self.grid[0] = np.arange(z_min, z_max, z_step)

        print('-> Grid counting domain:',
              ' lon:[', x_min, x_max, ']',
              ' lat:[', y_min, y_max, ']',
              ' depth:[', z_min, z_max, ']')

        print('-> Grid cells:',
              ' lon:[', self.grid[2].size, ']',
              ' lat:[', self.grid[1].size, ']',
              'depth:[', self.grid[0].size, ']'
              )

    def setCellBasics(self):
        "Sets the shape and number of cells"
        nz = self.cellCenters[0].size
        ny = self.cellCenters[1].size
        nx = self.cellCenters[2].size
        self.shape = [nz, ny, nx]
        self.nCells = nz*ny*nx

    def setCellCenters(self):
        "Set the cell centers based on the middle point between extremes."
        for value in self.grid:
            self.cellCenters.append((value[:-1] + value[1:])/2.)

    def setCellAreas(self, units='km'):
        """
        Set the horizontal area of each cell.

        Args:
            units (string, optional): Grid distance units. Defaults to 'km'.

        Returns:
            None.

        """

        dlon = (self.grid[2][1:] - self.grid[2][:-1])
        dlat = (self.grid[1][1:] - self.grid[1][:-1])
        y_c = (self.grid[1][1:] + self.grid[1][:-1])/2.
        # dx and dlon are in meters
        dx = dlon[np.newaxis, :] * (cte.degreesToRad * cte.earthRadius *
             np.cos(cte.degreesToRad * y_c[:, np.newaxis]))
        dy = dlat*cte.degreesToRad*cte.earthRadius
        if units == 'km':
            dx = dx/1000.
            dy = dy/1000.
        self.cellArea = dx*dy[:, np.newaxis]

    def setCellVolumes(self, units='km'):
        """
        Set the volume of each cell.

        Args:
            units (string, optional): Grid distance units. Defaults to 'km'.

        Returns:
            None.

        """
        dz = self.grid[0][1:]-self.grid[0][:-1]
        if units == 'km':
            dz = dz/1000.
        self.cellVolume = dz[:, np.newaxis, np.newaxis] * self.cellArea[np.newaxis, :, :]

    def setCoords(self):
        "Set the coords dictionary of the grid using the cell centers."
        self.coords = {self.dims[0]: ([self.dims[0]], self.cellCenters[0]),
                       self.dims[1]: ([self.dims[1]], self.cellCenters[1]),
                       self.dims[2]: ([self.dims[2]], self.cellCenters[2])}

    def PositionsToIdCell(self, particlePositions: np.array):
        """
        Turns [z,y,x] array position into [id] cell where the particle is.

        Args:
            particlePositions (np.array): [z,y,x] array position.

        Returns:
            None.

        """

        if np.size(particlePositions.shape) == 1:
            particlePositions = particlePositions[np.newaxis, :]

        z_dig = np.digitize(particlePositions[:, 0], self.grid[0], right=True)
        y_dig = np.digitize(particlePositions[:, 1], self.grid[1], right=True)
        x_dig = np.digitize(particlePositions[:, 2], self.grid[2], right=True)

        self.rIdCell = np.ravel_multi_index((z_dig, y_dig, x_dig),
                                            self.shape,
                                            mode='clip')

    def getCountsInCell(self, particlePositions: np.array) -> np.array:
        """
        Get the number of particles in each cell.

        Args:
            particlePositions (np.array): [z,y,x] array position.

        Returns:
            countsInCell (np.array): number of particles in each cell.

        """
        countsInCell, _ = np.histogramdd(particlePositions, self.grid)

        # self.validCells = self.countsInCell > 0
        # self.validCells = np.where(IdCounts > 0)[0]
        # self.countsInCell = np.reshape(IdCounts, (nz, ny, nx))
        return countsInCell

    def getMeanDataInCell(self, varData: np.array):
        """
        Computes the mean value of the particle properties inside a cell.

        Args:
            varData (np.array): scalar value in particle. [[scalar]]

        Returns:
            None.

        """
        cellMean = GridJIT.cellMeanData(self.rIdCell, self.nCells, 
                                        self.validCells, varData)
        self.meanDataInCell = np.reshape(cellMean, self.shape)

    def shape(self) -> list:
        """
        Get the shape of the grid

        Returns:
            list: size of each grid dimension.

        """
        return list(map(len, self.grid))

    def initializeGrid(self):
        """ Initialize the grid and associated properties."""
        self.setGrid()
        self.setCellCenters()
        self.setCellBasics()
        self.setCoords()
        self.setCellAreas()
        self.setCellVolumes()
