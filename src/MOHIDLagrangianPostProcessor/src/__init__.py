# -*- coding: utf-8 -*-

import os
import sys

basePath = os.path.dirname(os.path.realpath(__file__))
commonPath = os.path.abspath(os.path.join(basePath, "../../Common"))
sys.path.append(commonPath)

import os_dir
import MDateTime
from about import License