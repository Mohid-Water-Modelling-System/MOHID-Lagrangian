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
    argParser = argparse.ArgumentParser(description='Indexes input files for MOHID Lagrangian to parse. Use -h for help.')
    argParser.add_argument("-case", "--caseXML", dest="caseXML",
                    help=".xml file with the case definition for the MOHID Lagrangian run", metavar=".xml")
    argParser.add_argument("-i", "--input", dest="dataDir",
                    help="input directory containing .hdf5 files or subdirectories with them", metavar="dir")
    
    args = argParser.parse_args()

    dataDir = getattr(args,'dataDir')
    if dataDir == None: #reverting to the default
        basepath = os.path.dirname(__file__)
        dataDir = os.path.abspath(os.path.join(basepath, "inputData"))
    caseXML = getattr(args,'caseXML')
    
    #datadir="C:\Users\RBC_workhorse\Documents\GitHub\MOHID_python_tools\ConvertCSV2HDF5\testFiles"
    
    print('-> Case definition file is ', caseXML)
    print('-> Data files directory is', dataDir)
    
    #------------------------------------------------------    
    
            
run()