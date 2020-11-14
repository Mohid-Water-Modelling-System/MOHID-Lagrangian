# -*- coding: utf-8 -*-

import xml.etree.ElementTree as ET
import numpy as np
import MDateTime
from src.constants import daysToSeconds


class FilesTimesHandler:
    """ Class to manage the interations between the solution filelist with
    particle positions and timeinstants"""

    def __init__(self, fileList):
        self.fileList = fileList
        self.ISO_time_origin = MDateTime.getDateStringFromDateTime(MDateTime.BaseDateTime())
        self.timeAxis = []
        self.timeMask = []
        self.coords = []
        self.startTime = []
        self.endTime = []
        self.dt = []

    def setTimeAxisFromXML(self, xmlFile:str):
        """Sets a time-axis based on the file xml file and the vtu-files written.

        Args:
            xmlFile (str): path to the xml file.

        Returns:
            None.

        """
        root = ET.parse(xmlFile).getroot()
        for parameter in root.findall('execution/parameters/parameter'):
            if parameter.get('key') == 'Start':
                self.startTime = parameter.get('value')
            if parameter.get('key') == 'End':
                self.endTime = parameter.get('value')
            if parameter.get('key') == 'OutputWriteTime':
                self.dt = np.float(parameter.get('value'))
        nFiles = len(self.fileList)
        startTimeStamp = MDateTime.getTimeStampFromISODateString(self.startTime)
        self.timeAxis = np.array([startTimeStamp + i*self.dt/daysToSeconds for i in range(0,nFiles)])

    def setTimeMaskFromRecipe(self, xmlRecipe: str):
        """
        Sets a time-mask from the recipe xml file and vtu 'step' key.

        Args:
            xml xmlRecipe (str): path to the xml file.

        Returns:
            None.

        """
        root = ET.parse(xmlRecipe).getroot()
        self.timeMask = np.ones(self.timeAxis.size, dtype=np.bool)
        for parameter in root.findall('time/'):
            if parameter.tag == 'start':
                timeRecipeStart = MDateTime.getTimeStampFromISODateString(parameter.get('value'))
                self.timeMask = self.timeMask & (self.timeAxis > timeRecipeStart)
            if parameter.tag == 'end':
                timeRecipeEnd = MDateTime.getTimeStampFromISODateString(parameter.get('value'))
                self.timeMask = self.timeMask & (self.timeAxis < timeRecipeEnd)
            if parameter.tag == 'step':
                step = np.int64(parameter.get('value'))
                bufferTimeMask = np.zeros(self.timeAxis.size, dtype=np.bool)
                bufferTimeMask[::step] = True
                self.timeMask = self.timeMask & bufferTimeMask

    def updateTimeAxisWithMask(self):
        """ Updates the time axis with mask."""
        self.timeAxis = self.timeAxis[self.timeMask]

    def setTimeCoords(self):
        """Sets coords variable to pass to xarray methods"""
        self.coords = {'time': ('time', self.timeAxis.round(decimals=10))}

    def initializeTimeGrid(self, xmlFile, xmlRecipe):
        """ Initialized the time grid"""
        self.setTimeAxisFromXML(xmlFile)
        self.setTimeMaskFromRecipe(xmlRecipe)
        self.updateTimeAxisWithMask()
        self.setTimeCoords()

    def cropFileList(self) -> list:
        """ Crops the filelist based on time-mask"""
        tidx = 0
        croppedFileList = []
        for fileName in self.fileList:
            if self.timeMask[tidx] == True:
                croppedFileList.append(fileName)
            tidx += 1

        self.fileList = croppedFileList
        return croppedFileList