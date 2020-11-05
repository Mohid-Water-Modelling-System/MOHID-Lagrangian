# -*- coding: utf-8 -*-

import h5py
import os
import sys
import os_dir
import MDateTime
from src.VTUReader import *

basePath = os.path.dirname(os.path.realpath(__file__))
commonPath = os.path.abspath(os.path.join(basePath, "../Common"))
sys.path.append(commonPath)


def vtu2hdf5(vtuParser, FileTimeHandler, outDirLocal):
    outdirhdf5 = outDirLocal[:-1] + '_hdf5'
    if os.path.exists(outdirhdf5):
        os_dir.deleteDirForce(outdirhdf5)
    os.mkdir(outdirhdf5)
    f = 0
    print('-> Converting .vtu to .hdf5, MOHID formated')
    vtuFileList = vtuParser.fileList
    for vtuFile in vtuFileList:
        vtuParser.updateReaderWithFile(vtuFile)
        # get particle position
        r = vtuParser.getVariableData('coords')
        hdf5FileName = os_dir.filename_without_ext(os.path.basename(vtuFile))+'.hdf5'
        print('--> ' + os.path.basename(vtuFile) + ' -> ' + hdf5FileName)
        with h5py.File(outdirhdf5 + '/' + hdf5FileName ,'a') as hdf5File:
            # main groups
            grid = hdf5File.create_group("Grid")
            results = hdf5File.create_group("Results")
            time = hdf5File.create_group("Time")
            # subgroups
            group1 = results.create_group("Group_1/Data_1D")
            # writing data
            lon = group1.create_dataset('Longitude/Longitude_00001', data=r[:, 2], dtype='f')
            lon.attrs['Maximum'] = max(r[:, 2])
            lon.attrs['Minimum'] = min(r[:, 2])
            lon.attrs['Units'] = 'º'
            lat = group1.create_dataset('Latitude/Latitude_00001', data=r[:, 1], dtype='f')
            lat.attrs['Maximum'] = max(r[:, 1])
            lat.attrs['Minimum'] = min(r[:, 1])
            lat.attrs['Units'] = 'º'
            zz = group1.create_dataset('Z Pos/Z Position_00001', data=r[:, 0], dtype='f')
            zz.attrs['Maximum'] = max(r[:, 0])
            zz.attrs['Minimum'] = min(r[:, 0])
            zz.attrs['Units'] = 'm'
            r = vtuParser.getVariableData('source')
            source = group1.create_dataset('Origin ID/Origin ID_00001', data=r, dtype='f')
            source.attrs['Maximum'] = max(r)
            source.attrs['Minimum'] = min(r)
            source.attrs['Units'] = '-'
            # writing time
            dateArray = MDateTime.getMOHIDDateFromTimeStamp(FileTimeHandler.timeAxis[f])
            date = time.create_dataset('Time_00001', data=dateArray, dtype='f')
            date.attrs['Maximum'] = max(dateArray)
            date.attrs['Minimum'] = min(dateArray)
            date.attrs['Units'] = 'YYYY/MM/DD HH:MM:SS'
            f = f+1

