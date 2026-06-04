# -*- coding: utf-8 -*-

#    !------------------------------------------------------------------------------
#    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
#    !------------------------------------------------------------------------------
#    !
#    ! TITLE         : MDateTime
#    ! PROJECT       : Mohid python tools
#    ! MODULE        : background
#    ! URL           : http://www.mohid.com
#    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
#    ! DATE          : January 2019
#    ! REVISION      : Canelas 0.1
#    !> @author
#    !> Ricardo Birjukovs Canelas
#    !
#    ! DESCRIPTION:
#    !this encondes/decodes the MOHID date/time array to/from a float of total time since a reference date, in days 
#    !------------------------------------------------------------------------------
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

# Public API

def getTimeStampFromMOHIDDate(MohidDate):
    date = datetime(int(MohidDate[0]), int(MohidDate[1]), int(MohidDate[2]), int(MohidDate[3]), int(MohidDate[4]), int(MohidDate[5]))
    timeStamp = getTimeStampFromDateTime(date)
    return timeStamp

def getMOHIDDateFromTimeStamp(timeStamp):
    MD = getDateTimeFromTimeStamp(timeStamp)
    MohidDate = [MD.year, MD.month, MD.day, MD.hour, MD.minute, MD.second]
    return MohidDate

def getDateStringFromTimeStamp(timeStamp):
    return getDateTimeFromTimeStamp(timeStamp).strftime("%Y-%m-%d %H:%M:%S")

def getDateStringFromMOHIDDate(MohidDate):
    timeStamp = getTimeStampFromMOHIDDate(MohidDate)
    return getDateTimeFromTimeStamp(timeStamp).strftime("%Y-%m-%d %H:%M:%S")

def getTimeStampFromDateString(DateString):
    #string of the type '2000-08-19 01:01:37'
    date = datetime.strptime(DateString, "%Y-%m-%d %H:%M:%S")
    timeStamp = getTimeStampFromDateTime(date)
    return timeStamp


###private functions

def BaseDateTime():
    return datetime(1950, 1, 1, 0, 0, 0)

def getDateTimeFromTimeStamp(timeStamp):
    delta = timedelta(seconds=timeStamp*timedelta (days=1).total_seconds())
    return BaseDateTime() + delta

def getTimeStampFromDateTime(Date):
    delta = Date - BaseDateTime()
    timeStamp = delta.total_seconds()/timedelta(days=1).total_seconds()
    return timeStamp

#API examples
#MOHIDate = [2000, 8, 19, 1, 1, 37]
#print(getTimeStampFromMOHIDDate(MOHIDate))
#print(getMOHIDDateFromTimeStamp(getTimeStampFromMOHIDDate(MOHIDate)))
#print(getDateStringFromTimeStamp(getTimeStampFromMOHIDDate(MOHIDate)))
#print(getDateStringFromMOHIDDate(MOHIDate))
#print(getTimeStampFromDateString(getDateStringFromMOHIDDate(MOHIDate)))