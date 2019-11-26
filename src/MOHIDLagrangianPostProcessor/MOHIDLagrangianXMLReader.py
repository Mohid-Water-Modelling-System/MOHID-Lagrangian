# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 21:16:10 2019

@author: Toshiba
"""
import xml.etree.ElementTree as ET

def getSourcesDictFromXML(xmlFile):
    tree = ET.parse(xmlFile)
    sources={'id':{}}
    dict_source = {}
    for source in tree.find('caseDefinitions/sourceDefinitions')[:]:
        id_source = source.find('setsource').attrib['id']
        name_source = source.find('setsource').attrib['name']
        dict_source[id_source] = name_source
    dict_source['global'] = 'global'
    sources['id'] = dict_source
    return sources

def getBeachFromRecipe(xmlFile):
    root = ET.parse(xmlFile).getroot()
    for parameter in root.findall('EulerianMeasures/measures/filters/filter'):
        if parameter.get('key') == 'beaching':
            beach = parameter.get('value')
        else:
            beach = '0'
    return beach


def getRecipeListFromCase(xmlFile):
    recipeList =[]        
    #parsing case definition file
    root = ET.parse(xmlFile).getroot()                
    for fileName in root.findall('execution/postProcessing/file'):
        recipeList.append(fileName.get('name'))    
    return recipeList


def getFieldsFromRecipe(xmlFile):
    fieldList = []
    root = ET.parse(xmlFile).getroot()    
    for fieldName in root.findall('EulerianMeasures/measures/field'):
        fieldList.append(fieldName.get('key'))    
    return fieldList