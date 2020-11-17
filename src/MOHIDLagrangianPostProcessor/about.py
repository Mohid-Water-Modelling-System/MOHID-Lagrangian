# -*- coding: utf-8 -*-

#    !------------------------------------------------------------------------------
#    !        IST/MARETEC, Water Modelling Group, Mohid modelling system
#    !------------------------------------------------------------------------------
#    !
#    ! TITLE         : MOHIDLagrangianPreProcessor
#    ! PROJECT       : MOHIDLagrangian
#    ! URL           : http://www.mohid.com
#    ! AFFILIATION   : USC/GFNL, Group of Non Linear Physics
#    ! DATE          : August 2020
#    ! REVISION      : Garaboa 0.2
#    !> @author
#    !> Daniel Garaboa Paz, Ricardo Birjukovs Canelas
#    !
#    ! DESCRIPTION:
#    !PostProcessing script for MOHID Lagrangian. Computes concentrations,
#    !residence times, etc
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


class License:
    def __init__(self):
        self.art = '''
 __  __  ____  _    _ _____ _____  _                                       _             
|  \/  |/ __ \| |  | |_   _|  __ \| |                                     (_)            
| \  / | |  | | |__| | | | | |  | | |     __ _  __ _ _ __ __ _ _ __   __ _ _  __ _ _ __  
| |\/| | |  | |  __  | | | | |  | | |    / _` |/ _` | '__/ _` | '_ \ / _` | |/ _` | '_ \ 
| |  | | |__| | |  | |_| |_| |__| | |___| (_| | (_| | | | (_| | | | | (_| | | (_| | | | |
|_|  |_|\____/|_|  |_|_____|_____/|______\__,_|\__, |_|  \__,_|_| |_|\__, |_|\__,_|_| |_|
                                                __/ |                 __/ |              
                                               |___/                 |___/               
 _____          _   _____                                        
|  __ \        | | |  __ \                                       
| |__) |__  ___| |_| |__) | __ ___   ___ ___  ___ ___  ___  _ __ 
|  ___/ _ \/ __| __|  ___/ '__/ _ \ / __/ _ \/ __/ __|/ _ \| '__|
| |  | (_) \__ \ |_| |   | | | (_) | (_|  __/\__ \__ \ (_) | |   
|_|   \___/|___/\__|_|   |_|  \___/ \___\___||___/___/\___/|_|   
                                                                 
                                                                                        
'''
        self.author = 'D. Garaboa Paz - USC-GFNL'
        self.version = ' version 0.1'
        self.date = ' 2-05-2019'
        self.versionInfo = self.author + self.version + self.date
        self.lic = '''
    MIT License

    Copyright (c) 2019 DGaraboaPaz

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
