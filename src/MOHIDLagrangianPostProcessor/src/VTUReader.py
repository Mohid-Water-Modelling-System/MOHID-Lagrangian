#!/usr/bin/env python3

import numpy as np
import glob
from vtk.util.numpy_support import vtk_to_numpy


def getVtuFileList(outDir):
    return glob.glob(outDir+'/*_?????.vtu').sort()[1:]


def getVtuParentFile(outDir):
    return glob.glob(outDir+'/*_?????.vtu').sort()[0]


def getSourceArrayFromVTU(VTKReader) -> np.array:
    return vtk_to_numpy(VTKReader.GetOutput().GetPointData().GetArray('source'))


def getStateArrayFromVTU(VTKReader) -> np.array:
    state_array = VTKReader.GetOutput().GetPointData().GetArray('state')
    if state_array == None:
        return np.zeros(1)
    else:
        return vtk_to_numpy(state_array)

def getCoordsArrayFromVTU(VTKReader) -> np.array:
    return vtk_to_numpy(VTKReader.GetOutput().GetPoints().GetData())[:, ::-1]


def getVelocityModuleFromVTU(VTKReader) -> np.array:
    vtu_vars = vtk_to_numpy(VTKReader.GetOutput().GetPointData().GetArray('velocity'))
    vtu_vars = np.sqrt(vtu_vars[:, 0]**2 + vtu_vars[:, 1]**2 + vtu_vars[:, 2]**2)
    return vtu_vars


def getVariableArrayFromVTU(VTKReader, variableName) -> np.array:
    if variableName == 'coords':
        vtu_vars = getCoordsArrayFromVTU(VTKReader)
    elif variableName == 'velocity':
        vtu_vars = getVelocityModuleFromVTU(VTKReader)
    elif variableName in VTKReader.availableVtuVars:
        vtu_vars = getVariableArrayFromVTU(VTKReader, variableName)
    return vtu_vars


def getBeachMaskFromVTU(VTKReader, beachCondition = 0) -> np.bool:
    if beachCondition:
        state = getStateArrayFromVTU(VTKReader)
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
            sourceMask = getSourceArrayFromVTU(VTKReader) == np.int(source)
    return sourceMask
