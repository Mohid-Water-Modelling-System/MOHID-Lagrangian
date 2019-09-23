# -*- coding: utf-8 -*-

import os
import sys
import glob

basePath = os.path.dirname(os.path.realpath(__file__))
commonPath = os.path.abspath(os.path.join(basePath, "../../Common"))
sys.path.append(commonPath)
import os_dir

import MOHIDCurvilinear2Simple


def run():
    
    directories = []    
    directories.append('C:/Users/RBC_workhorse/Documents/GitHub/MOHID-Lagrangian/run_local_cases/PCOMS_fitted_domains/base_data/pcoms')
    directories.append('C:/Users/RBC_workhorse/Documents/GitHub/MOHID-Lagrangian/run_local_cases/PCOMS_fitted_domains/base_data/tagus')
    
    fileExtensions = ['.nc', '.nc4']
    
    for idir in directories:
        outdirLocal = idir+'/simple_grid'
        if os.path.exists(outdirLocal):
            os_dir.deleteDirForce(outdirLocal)
        os.mkdir(outdirLocal)
        for ext in fileExtensions:
            files = glob.glob(idir+'/*'+ext, recursive=True)
            for file in files:
                MOHIDCurvilinear2Simple.MOHIDCurvilinearCorrect(file, outdirLocal)
                








run()