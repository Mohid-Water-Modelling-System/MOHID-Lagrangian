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
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
import glob

basePath = os.path.dirname(os.path.realpath(__file__))
commonPath = os.path.abspath(os.path.join(basePath, "Common"))
sys.path.append(commonPath)
import os_dir
import about
import ncMetaParser
import xmlWriter


def run():
    
    lic = about.Licence()
    lic.print()
    
    #cmd line argument parsing---------------------------
    argParser = argparse.ArgumentParser(description='Indexes input files for MOHID Lagrangian to parse. Use -h for help.')
    argParser.add_argument("-i", "--input", dest="caseXML",
                    help=".xml file with the case definition for the MOHID Lagrangian run", metavar=".xml")
    argParser.add_argument("-o", "--outputDir", dest="outDir",
                    help="output directory", metavar=".xml")
    args = argParser.parse_args()
    
    caseXML = getattr(args,'caseXML')
    outDir = getattr(args,'outDir')
    print('-> Case definition file is ', caseXML)
    #---------------------------------------------------
    
    #parsing case definition file
    root = ET.parse(caseXML).getroot()
    
    dataDir = []
    for type_tag in root.findall('casedef/inputData/inputDataDir'):
        dataDir.append(type_tag.get('name'))
    
    for type_tag in root.findall('execution/parameters/parameter'):
        if type_tag.get('key') == 'StartTime':
            StartTime = datetime.strptime(type_tag.get('value'), "%Y %m %d %H %M %S")            
        if type_tag.get('key') == 'EndTime':
            EndTime = datetime.strptime(type_tag.get('value'), "%Y %m %d %H %M %S")            
        
    #dataDir="C:\Users\RBC_workhorse\Documents\GitHub\MOHID_python_tools\ConvertCSV2HDF5\testFiles"
    
    if len(dataDir) > 1:
        print('-> Input data directories are', dataDir)
    else:
        print('-> Input data directory is', dataDir)
    #------------------------------------------------------
    fileExtensions = ['.nc', '.nc4']  
    
    #going for each input directory and indexing its files
    inputFiles = []
    for idir in dataDir:
        for ext in fileExtensions:
            inputFiles.append(glob.glob(idir+ '/**/*'+ext, recursive=True))
    #cleaning list of empty values
    inputFiles = list(filter(None, inputFiles))
	
    if not inputFiles:

        print('No input files found. Supported files are ', fileExtensions)
    
    else:
    
        indexerFileName = os_dir.filename_without_ext(caseXML)+'_inputs'
        indexer = xmlWriter.xmlWriter(indexerFileName)
		
		#going trough every file, extracting some metadata and writting in the indexer file
        ncMeta = []
        for idir in inputFiles:
            for ifile in idir:
                print('--> reading file', ifile)
                ncMeta.append(ncMetaParser.ncMetadata(ifile, StartTime))
		
        ncMeta.sort(key=lambda x: x.startTime)
		
        indexer.openCurrentsCollection()
		
        print('--> indexing currents data')
        for ncfile in ncMeta:
            indexer.writeFile(ncfile.getName(), ncfile.getstartTime(), ncfile.getendTime(), ncfile.getstartDate().strftime("%Y/%m/%d, %H:%M:%S"), ncfile.getendDate().strftime("%Y/%m/%d, %H:%M:%S"))
		
        indexer.closeCurrentsCollection()
        indexer.closeFile()
            
run()