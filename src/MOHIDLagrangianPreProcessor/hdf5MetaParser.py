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

_EPOCH = datetime(1970, 1, 1)


def instantToSeconds(t):
    # Reduce one time entry to a scalar timestamp. MOHID stores each instant as a
    # 6-element [Y,M,D,h,m,s] array; already-scalar formats are passed through.
    a = np.asarray(t).ravel()
    if a.size == 6:
        d = datetime(int(a[0]), int(a[1]), int(a[2]), int(a[3]), int(a[4]), int(a[5]))
        return (d - _EPOCH).total_seconds()
    return float(a)


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

        # use only Time_XXXXX datasets, ordered by their numeric suffix
        time_keys = sorted([k for k in xds.variables if k.startswith('Time_')],
                           key=lambda k: int(k.split('_')[1]))
        number_of_instants = len(time_keys)
        time_min = xds[time_keys[0]]
        time_max = xds[time_keys[-1]]

        for k in time_keys:
            self.time.append(xds[k].data)

        time_start = time_min.data
        time_end = time_max.data

        self.startDate = datetime(int(time_start[0]), int(time_start[1]), int(time_start[2]), int(time_start[3]), int(time_start[4]), int(time_start[5]))
        self.endDate = datetime(int(time_end[0]), int(time_end[1]), int(time_end[2]), int(time_end[3]), int(time_end[4]), int(time_end[5]))
        self.startTime = (self.startDate - baseTime).total_seconds()
        self.endTime = (self.endDate - baseTime).total_seconds()

        xds.close()                      # release the opened datasets
        try:
            ncf.close()
        except RuntimeError:
            pass

    def getTimeString(self, i):

        if i > 999:
            str_time = 'Time_0' + str(i)
        elif i > 99:
            str_time = 'Time_00' + str(i)
        elif i > 9:
            str_time = 'Time_000' + str(i)
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
            # collapse each [Y,M,D,h,m,s] instant to a scalar so the time axis is 1-D
            # and stays the same length as the filename axis
            for t in hdf5_meta.time:
                time_axis.append(instantToSeconds(t))
                time_axis_filename.append(hdf5_meta.fileName)

        # build one dimension time axis
        time_axis = np.asarray(time_axis, dtype=float)
        time_axis_filename = np.asarray(time_axis_filename)

        print('-> Checking time integrity through files... ')

        # test #1: Seek for repeated values.
        # They ill produce gaps -> If exist Return and skip the second test.
        values, counts = np.unique(time_axis, return_counts=True)
        repeated_values = values[counts > 1]
        if repeated_values.size > 0:
            print(' -> There are repeated values in your time axis')
            mask_repeated = np.isin(time_axis, repeated_values)
            problem_files = time_axis_filename[mask_repeated]
            problem_steps = time_axis[mask_repeated]
            problem_type = ['repeated-values' for i in range(0, len(problem_steps))]
            for idx in range(0, len(problem_files)):
                print('->', problem_files[idx],'|', problem_steps[idx], '|', problem_type[idx])
            return

        # test #2: Seek for wholes in data.
        # sort by time so dt is computed between consecutive timestamps
        order = np.argsort(time_axis)
        time_axis = time_axis[order]
        time_axis_filename = time_axis_filename[order]
        dt = np.diff(time_axis)
        unique_dt, counts_dt = np.unique(dt, return_counts=True)
        most_repeated_dt = unique_dt[np.argmax(counts_dt)] # most common dt value is the most frequent one.

        mask_gap = np.zeros_like(dt, dtype=bool)
        #Seek for values where the 'dt' is NOT the most rcommon
        mask_gap[dt != most_repeated_dt] = True

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
