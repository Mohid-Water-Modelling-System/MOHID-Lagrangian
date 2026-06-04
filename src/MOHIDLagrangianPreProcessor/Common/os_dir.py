# -*- coding: utf-8 -*-

import os

def get_immediate_subdirectories(a_dir):
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]
    
def get_contained_files(a_dir,extension):
    return [file for file in os.listdir(a_dir)
            if file.endswith(extension)]

def mkdir_safe(a_dir):
    if not os.path.exists(a_dir):
        os.makedirs(a_dir)
        
def filename_without_ext(path_to_file): 
    return os.path.splitext(path_to_file)[0]