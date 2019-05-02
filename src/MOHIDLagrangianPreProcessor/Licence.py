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
#    Copyright (c) 2019 RBCanelas
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


class Licence:
    def __init__(self):
        self.art = '''
  __  __  ____  _    _ _____ _____    _                                       _             
 |  \/  |/ __ \| |  | |_   _|  __ \  | |                                     (_)            
 | \  / | |  | | |__| | | | | |  | | | |     __ _  __ _ _ __ __ _ _ __   __ _ _  __ _ _ __  
 | |\/| | |  | |  __  | | | | |  | | | |    / _` |/ _` | '__/ _` | '_ \ / _` | |/ _` | '_ \ 
 | |  | | |__| | |  | |_| |_| |__| | | |___| (_| | (_| | | | (_| | | | | (_| | | (_| | | | |
 |______|\____/|______|_____|_____/  |______\__,_|\__, |_|  \__,_|_| |_|\__, |_|\__,_|_| |_|
 |  __ \        |  __ \                            __/ |                 __/ |              
 | |__) _ __ ___| |__) _ __ ___   ___ ___ ___ ___ |___/ _ __            |___/               
 |  ___| '__/ _ |  ___| '__/ _ \ / __/ _ / __/ __|/ _ \| '__|                               
 | |   | | |  __| |   | | | (_) | (_|  __\__ \__ | (_) | |                                  
 |_|   |_|  \___|_|   |_|  \___/ \___\___|___|___/\___/|_|                                  
'''
        self.author = 'R. Birjukovs Canelas - MARETEC'
        self.version = ' version 0.1'
        self.date = ' 2-05-2019'
        self.versionInfo = self.author + self.version + self.date
        self.lic = '''
    MIT License
    
    Copyright (c) 2019 RBCanelas
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

'''
        
    def print(self):
        print(self.art)
        print(self.versionInfo)
        print(self.lic)
