# -*- coding: utf-8 -*-
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

class xmlWriter:
    def __init__(self, fileName):
        self.filename = fileName
        
        self.openFile()
        
    def openFile(self):
        self.f = open(self.filename + '.xml', 'w')
        self.writeHeader()
    
    def closeFile(self):
        self.f.write('''</file_collection>''')
        self.f.close()
        
    def writeHeader(self):
        self.f.write('''<?xml version="1.0" encoding="UTF-8" ?>
<file_collection>
''')

    def openCurrentsCollection(self):
        self.f.write('''	<currents>
''')
		
    def closeCurrentsCollection(self):
        self.f.write('''	</currents>
''')

    def writeFile(self,fileName,startTime,endTime,startDateStr,endDateStr):                        
        toWrite = '''    	<file>
			<name value="'''+fileName+'''" />
			<startTime value="'''+str(startTime)+'''" />	<!-- '''+startDateStr+'''-->
			<endTime value="'''+str(endTime)+'''" />	<!-- '''+endDateStr+'''-->
		</file>
'''
        self.f.write(toWrite)
        
