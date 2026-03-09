# -*- coding: utf-8 -*-

#    !------------------------------------------------------------------------------
#    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
#    !------------------------------------------------------------------------------
#    !
#    ! TITLE         : MOHIDLagrangianPreProcessor
#    ! PROJECT       : MOHIDLagrangian
#    ! URL           : http://www.mohid.com
#    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
#    ! DATE          : June 2019
#    ! REVISION      : Garaboa 0.1
#    !> @author
#    !> Angel Daniel Garaboa Paz
#    !
#    ! DESCRIPTION:
#    ! Preprocessing script for MOHID Lagrangian. Lists input files, composes config 
#    ! files, etc
#    !------------------------------------------------------------------------------
#
#    MIT License
#
#    Copyright (c) 2018 DGaraboa
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


# This environment variable avoids error on locking file when writing.
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

basePath = os.path.dirname(os.path.realpath(__file__))
commonPath = os.path.abspath(os.path.join(basePath, "../Common"))
sys.path.append(commonPath)

import os_dir
import MDateTime
from about import License
from src.XMLReader import *
from src.PlotResults import *
from src.PostProcessor import PostProcessor
from src.GridBase import GridBase


def main():
    lic = License()
    lic.print()

    # cmd line argument parsing
    argParser = argparse.ArgumentParser(description='Post processes MOHID Lagrangian outputs. Use -h for help.')
    argParser.add_argument("-i", "--input", dest="caseXML",
                    help=".xml file with the case definition for the MOHID Lagrangian run", metavar=".xml")
    argParser.add_argument("-f", "--force", dest="recipeXML",
                    help=".xml file with the recipe for a post process run - optional", metavar=".xml")
    argParser.add_argument("-o", "--outputDir", dest="outDir",
                    help="output directory", metavar="dir")
    argParser.add_argument("-po", "--plot-only", dest="plotonly",
                    help="output directory", action='store_true')
    args = argParser.parse_args()


    caseXML = getattr(args, 'caseXML')
    recipeXML = []
    recipeXML.append(getattr(args, 'recipeXML'))
    outDir = getattr(args, 'outDir')

    print('-> Case definition file is', caseXML)
    print('-> Main output directory is', outDir)

    # get list of post cicles to run
    if recipeXML == [None]:
        # open caseXML and extract recipeXML names
        recipeXML = getRecipeListFromCase(caseXML)

    for recipe in recipeXML:
        print('-> Running recipe', recipe)
        outDirLocal = outDir + '/postProcess_' + os.path.basename(os_dir.filename_without_ext(recipe)) + '/' 

        # If plotonly flag > Runs only plotting stage
        if args.plotonly is True:
            plotResultsFromRecipe(outDirLocal, recipe)
            return

        postProcessor = PostProcessor(caseXML, recipe, outDir, outDirLocal)
        postProcessor.run()

        if checkPlotRecipe(recipe):
            plotResultsFromRecipe(outDirLocal, recipe)

main()
