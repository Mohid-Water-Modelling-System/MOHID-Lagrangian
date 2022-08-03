# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 21:16:10 2019

@author: Toshiba
"""
import xml.etree.ElementTree as ET


def getSourcesDictFromXML(xmlFile):
    tree = ET.parse(xmlFile)
    sources = {'id': {}}
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
    beach = '0'
    for parameter in root.findall('EulerianMeasures/measures/filters/filter'):
        if parameter.get('key') == 'beaching':
            beach = parameter.get('value')
        else:
            beach = '0'
    return beach


def getRecipeListFromCase(xmlFile):
    recipeList = []
    # parsing case definition file
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


def checkHDF5WriteRecipe(xmlFile):
    convert = False
    root = ET.parse(xmlFile).getroot()
    for formats in root.findall('convertFiles/format'):
        if formats.get('key') == 'hdf5':
            convert = True
    return convert


def checkPlotRecipe(xmlFile):
    plot = False
    root = ET.parse(xmlFile).getroot()
    for plots in root.findall('plot'):
        if len(plots) > 1:
            plot = True
    return plot


def getPlotTimeFromRecipe(xmlFile):
    root = ET.parse(xmlFile).getroot()
    fieldName = root.findall('plot/time')[0]
    freq = fieldName.get('value')
    kind = fieldName.get('key')
    return freq, kind


def getPlotMeasuresFromRecipe(xmlFile):
    fieldList = []
    root = ET.parse(xmlFile).getroot()
    for fieldName in root.findall('plot/measure'):
        fieldList.append(fieldName.get('key'))
    return fieldList


def getPolygonFileFromRecipe(xmlFile):
    fieldList = []
    root = ET.parse(xmlFile).getroot()
    for fieldName in root.findall('EulerianMeasures/polygonDefinition/file'):
        fieldList.append(fieldName.get('name'))
    return fieldList


def isPolygonOrGrid(xmlFile):
    if len(getPolygonFileFromRecipe(xmlFile)) >= 1:
        postType = 'polygon'
    else:
        postType = 'grid'
    return postType


def getPlotTypeFromRecipe(xmlFile):
    fieldList = []
    root = ET.parse(xmlFile).getroot()
    for fieldName in root.findall('plot/type'):
        fieldList.append(fieldName.get('value'))
    return fieldList


def getPlotWeightFromRecipe(xmlFile):
    fieldList = []
    root = ET.parse(xmlFile).getroot()
    for fieldName in root.findall('plot/weight'):
        fieldList.append(fieldName.get('file'))
    return fieldList


def getGroupFromRecipe(xmlFile):
    fieldList = []
    root = ET.parse(xmlFile).getroot()
    for fieldName in root.findall('plot/group'):
        fieldList.append(fieldName.get('key'))
    return fieldList


def getNormalizeFromRecipe(xmlFile):
    fieldList = []
    root = ET.parse(xmlFile).getroot()
    for fieldName in root.findall('plot/normalize'):
            fieldList.append(fieldName.get('key'))
    return fieldList

