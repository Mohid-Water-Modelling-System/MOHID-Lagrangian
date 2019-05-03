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

class ncMetadata:
    def __init__(self, fileName, baseTime):
        self.fileName = []
        self.startTime = []
        self.endTime = []
        
        self.fileName = fileName        
        ds = xr.open_dataset(self.fileName)
        tMin = ds.time.min()
        tMax = ds.time.max()
        dMin = datetime(tMin.dt.year, tMin.dt.month, tMin.dt.day, tMin.dt.hour, tMin.dt.minute, tMin.dt.second)
        dMax = datetime(tMax.dt.year, tMax.dt.month, tMax.dt.day, tMax.dt.hour, tMax.dt.minute, tMax.dt.second)
        self.startTime = (dMin - baseTime).total_seconds()
        self.endTime = (dMax - baseTime).total_seconds()
        
        ds.close
        
    def getName(self):
        return self.fileName
    
    def startTime(self):
        return self.startTime
    
    def endTime(self):
        return self.endTime
        