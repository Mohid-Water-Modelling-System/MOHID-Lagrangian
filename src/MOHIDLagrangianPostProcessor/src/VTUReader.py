#!/usr/bin/env python3

import numpy as np
import glob
from vtk.util.numpy_support import vtk_to_numpy


def getVtuFileList(outDir):
    return glob.glob(outDir+'/*_?????.vtu').sort()[1:]


def getVtuParentFile(outDir):
    return glob.glob(outDir+'/*_?????.vtu').sort()[0]


def getVariableFromVTU(VTKReader, variableName) -> np.array:
    if variableName == 'coords':
        vtu_vars = vtk_to_numpy(VTKReader.GetOutput().GetPoints().GetData())[:, ::-1]
    elif variableName == 'velocity':
        vtu_vars = vtk_to_numpy(VTKReader.GetOutput().GetPointData().GetArray('velocity'))
        vtu_vars = np.sqrt(vtu_vars[:, 0]**2 + vtu_vars[:, 1]**2 + vtu_vars[:, 2]**2)
    else:
        vtk_array = VTKReader.GetOutput().GetPointData().GetArray(variableName)
        if vtk_array is None:
            vtu_vars = np.zeros(1)
        else:
            vtu_vars = vtk_to_numpy(vtk_array)
    return vtu_vars


def getBeachMaskFromVTU(VTKReader, beachCondition) -> np.bool:
    if beachCondition:
        state = getVariableFromVTU(VTKReader, 'state')
        if beachCondition == '0':
            beachMask = True
        elif beachCondition == '1':
            beachMask = state < 0.5
        elif beachCondition == '2':
            beachMask = state >= 0.5
    else:
        beachMask = np.bool(True)
    return beachMask


def getSourceMaskFromVTU(VTKReader, source: int) -> np.bool:
    if source:
        if source == 'global':
            sourceMask = np.bool(True)
        else:
            sourceMask = getVariableFromVTU(VTKReader, 'source') == np.int(source)
    return sourceMask
