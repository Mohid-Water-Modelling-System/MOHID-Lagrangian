# -*- coding: utf-8 -*-
#    
#    MIT License
#    
#    Copyright (c) 2018 RBCanelas
#    
#    Permission is hereby granted, free of charge, to any person obtaining a copy
#    of this software and associated documentation files (the "Software"), to deal
#    in the Software without restriction, including without limitation the rights
#    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#    copies of the Software, and to permit persons to whom the Software is
#    furnished to do so, subject to the following conditions:
#    
#    The above copyright notice and this permission notice shall be included in all
#    copies or substantial portions of the Software.
#    
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#    SOFTWARE.

from datetime import datetime, timedelta
import xarray as xr
import numpy as np
import netCDF4 as nc

class hdf5Metadata:
    def __init__(self, fileName, baseTime):
        self.fileName = []
        self.startTime = []
        self.endTime = []
        self.startDate = []
        self.endDate = []

        self.fileName = fileName
        self.time = []
        
        ncf = nc.Dataset(fileName, diskless=True, persist=False)
        nch = ncf.groups.get('Time')
        xds = xr.open_dataset(xr.backends.NetCDF4DataStore(nch))
        
        time_min = xds['Time_00001']
        number_of_instants = len(xds)
        time_max = xds[self.getTimeString(number_of_instants)]
        
        for i in range(1, number_of_instants+1):
            self.time.append(xds[self.getTimeString(i)].data)
            
        time_start = time_min.data
        time_end = time_max.data
        
        self.startDate = datetime(int(time_start[0]), int(time_start[1]), int(time_start[2]), int(time_start[3]), int(time_start[4]), int(time_start[5]))
        self.endDate = datetime(int(time_end[0]), int(time_end[1]), int(time_end[2]), int(time_end[3]), int(time_end[4]), int(time_end[5]))
        self.startTime = (self.startDate - baseTime).total_seconds()
        self.endTime = (self.endDate - baseTime).total_seconds()
        
    def getTimeString(self, i):
        if i > 9:
            str_time = 'Time_000' + str(i)
        elif i > 99:
            str_time = 'Time_00' + str(i)
        elif i > 999:
            str_time = 'Time_0' + str(i)
        else:
            str_time = 'Time_0000' + str(i)
        return str_time

    def getName(self):
        return self.fileName

    def getstartTime(self):
        return self.startTime

    def getendTime(self):
        return self.endTime

    def getstartDate(self):
        return self.startDate

    def getendDate(self):
        return self.endDate


class hdf5DimParser:
    """
    Class with functions to check that data in the preprocesing stage is ok.

    """

    def checkTime(hdf5MetadataList: list):
        """
        Check that data from all the netcdf files has a good time dimension.
        Prints problematic files as warning. It continues in any case.

        Args:
            hdf5MetadataList (list): hdf5Metadata sorted list (by startTime)

        Returns:
            None.

        """

        time_axis = []
        time_axis_filename = []
        for hdf5_meta in hdf5MetadataList:
            nsteps = len(hdf5_meta.time)
            time_axis.append(hdf5_meta.time)
            time_axis_filename.append([hdf5_meta.fileName for i in range(0, nsteps)])

        # build one dimension time axis
        time_axis = np.hstack(time_axis)
        time_axis_filename = np.hstack(time_axis_filename)

        print('-> Checking time integrity through files... ')

        # test #1: Seek for repeated values.
        # They ill produce gaps -> If exist Return and skip the second test.
        mask_repeated = (np.array([np.sum(time == time_axis) for time in time_axis])) > 1
        if np.any(mask_repeated):
            print(' -> There are repeated values in your time axis')
            problem_files = time_axis_filename[mask_repeated]
            problem_steps = time_axis[mask_repeated]
            problem_type = ['repeated-values' for i in range(0, len(problem_steps))]
            for idx in range(0, len(problem_files)):
                print('->', problem_files[idx],'|', problem_steps[idx], '|', problem_type[idx])
            return

        # test #2: Seek for wholes in data. 
        dt = np.diff(time_axis)
        unique_dt = np.unique(dt)
        most_repeated_dt = unique_dt[-1] # most common dt value is the last index.

        mask_gap = np.zeros_like(dt, dtype=bool)
        #Seek for values where the 'dt' is NOT the most rcommon 
        for non_common_values in unique_dt[:-1]:
            mask_gap[dt == non_common_values] = True

        # Append a value at the end - match dimension with time_axis
        mask_gap = np.append(mask_gap, False)

        if np.any(mask_gap):
            print('-> There are time gaps in your nc-files.')
            problem_files = time_axis_filename[mask_gap]
            problem_steps = time_axis[mask_gap]
            problem_type = ['gaps-in-data' for i in range(0, len(problem_steps))]
            for i in range(0, len(problem_files)):
                print('->', problem_files[i], '|', problem_steps[i], '|', problem_type[i])
            return

