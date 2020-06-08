# -*- coding: utf-8 -*-
import glob
from src.VTUParser import VTUParser


class PVDParser:

    def __init__(self, pvdFile):
        self.pvdFile = pvdFile
        self.vtuFilelist = []
        self.vtuFileHandlers = []
        self.vtuParentFile = []

    def getVtuFileList(self, outDir):
        vtu_list = glob.glob(outDir+'/*_?????.vtu')
        vtu_list.sort()
        self.vtuParentFile = vtu_list[0]
        self.vtuFilelist = vtu_list[1:]

    def getVtuFileHandlers(self, outDir):
        self.getVtuFileList(outDir)
        VTUParser.getVtuVarsFromInitialFile(self.vtuParentFile)
        for vtuFile in self.vtuFilelist:
            self.vtuFileHandlers.append(VTUParser(vtuFile))

    def updateVtuFileHandlers(self, fileList):
        self.vtuFileHandlers = []
        self.vtuFilelist = fileList
        for vtuFile in self.vtuFilelist:
            self.vtuFileHandlers.append(VTUParser(vtuFile))
