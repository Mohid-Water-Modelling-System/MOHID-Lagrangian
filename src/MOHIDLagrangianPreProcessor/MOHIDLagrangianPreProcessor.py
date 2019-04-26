# -*- coding: utf-8 -*-

#    !------------------------------------------------------------------------------
#    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
#    !------------------------------------------------------------------------------
#    !
#    ! TITLE         : MOHIDLagrangianPreProcessor
#    ! PROJECT       : MOHIDLagrangian
#    ! URL           : http://www.mohid.com
#    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
#    ! DATE          : April 2019
#    ! REVISION      : Canelas 0.1
#    !> @author
#    !> Ricardo Birjukovs Canelas
#    !
#    ! DESCRIPTION:
#    !Preprocessing script for MOHID Lagrangian. Lists input files, composes config 
#    !files, etc 
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

import os
import sys
import argparse

sys.path.append('Common')
import os_dir


def run():
    
    #cmd line argument parsing---------------------------
    argParser = argparse.ArgumentParser(description='Indexes MOHID outputs in xmdf files. Use -h for help.')
    argParser.add_argument("-i", "--input", dest="datadir",
                    help="input directory containing .hdf5 files or subdirectories with them", metavar="dir")
    argParser.add_argument("-g", "--glue", dest="glueFiles", default=False,
                    help="option to atempt to produce a master indexer, that 'glues' all of the files")
    argParser.add_argument("-fd", "--firstdate", dest="firstDate", default='',
                    help="option to control the first date for the master indexer, format as 2000-08-19 01:01:37")
    argParser.add_argument("-ld", "--lastdate", dest="lastDate", default='',
                    help="option to control the last date for the master indexer, format as 2000-08-19 01:01:37")
    args = argParser.parse_args()

    datadir = getattr(args,'datadir')
    if datadir == None: #reverting to the test files
        basepath = os.path.dirname(__file__)
        datadir = os.path.abspath(os.path.join(basepath, "..", "TestFiles"))        
    glueFiles = getattr(args,'glueFiles')
    if glueFiles != False:
        glueFiles = True
    firstDate = str(getattr(args,'firstDate'))
    lastDate = str(getattr(args,'lastDate'))
    
    #datadir="C:\Users\RBC_workhorse\Documents\GitHub\MOHID_python_tools\ConvertCSV2HDF5\testFiles"
    #glueFiles = True
        
    print('-> Main directory is', datadir)
    if glueFiles:
        print('-> Attempting to glue files')
    if firstDate != '':
        print('-> First date to be indexed is', firstDate)
    if lastDate != '':
        print('-> First date to be indexed is', lastDate)
    #------------------------------------------------------    
    
            
run()