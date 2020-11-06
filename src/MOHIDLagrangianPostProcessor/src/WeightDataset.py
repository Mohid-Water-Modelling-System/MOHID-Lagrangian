# -*- coding: utf-8 -*-


import pandas as pd
import xarray as xr


def get_source_from_variable(variable: str):
    """
    get the source name from variable

    Args:
        variable (str): variable name from postprocessor

    Returns:
        source_name (str): source name

    """

    post_measures = ['concentration_volume_',
                     'concentration_area_',
                     'n_counts_',
                     'residence_time_']

    for post_measure in post_measures:
        if post_measure in variable:
            source_name = variable.replace(post_measure, '')

    return source_name


def get_material_from_variable(variable: str):
    """
    get the source name from variable

    Args:
        variable (str): variable name from postprocessor

    Returns:
        source_name (str): source name

    """

    post_measures = ['concentration_volume_',
                     'concentration_area_',
                     'n_counts_',
                     'residence_time_']

    material_names = ['light', 'heavy', 'medium']

    for post_measure in post_measures:
        if post_measure in variable:
            source_name = variable.replace(post_measure, '')

    for material_name in material_names:
        if material_name in source_name:
            material_name = source_name.replace(material_name, '')

    return material_name


def get_measure_from_variable(variable: str):
    """
    get the measure from variable

    Args:
        variable (str): DESCRIPTION.

    Returns:
        post_measure (TYPE): DESCRIPTION.

    """

    post_measures = ['concentration_volume_', 'concentration_area_',
                     'n_counts_', 'residence_time_']

    for post_measure in post_measures:
        if post_measure in variable:
            return post_measure


def weight_dataset_with_csv(dataset:xr.Dataset, weight_file:str) -> xr.Dataset:
    """

    Args:
        dataset (xr.Dataset): DESCRIPTION.
        weight_file (str): DESCRIPTION.

    Returns:
        dataset (TYPE): DESCRIPTION.

    """

    df = pd.read_csv(weight_file)
    df = df.set_index('source')
    print("-> %-30s | %9s" % ('Source', 'Weight'))
    print("")
    for datasetVar in list(dataset.keys()):
        for indexVar in df.index:
            if indexVar in datasetVar:
                weight = df.loc[indexVar, 'weight']
                dataset[datasetVar] = dataset[datasetVar]*df.loc[indexVar, 'weight']
                print("-> %-30s | %4.1f" % (datasetVar, weight))

    return dataset


def weight_dataarray_with_csv(dataArray:xr.DataArray, weight_file:str) -> xr.DataArray:
    """


    Args:
        dataset (xr.Dataset): DESCRIPTION.
        weight_file (str): DESCRIPTION.

    Returns:
        dataset (TYPE): DESCRIPTION.

    """

    df = pd.read_csv(weight_file)
    df = df.set_index('source')
    print("-> %-30s | %9s" % ('Source', 'Weight'))
    print("")
    for indexVar in df.index:
        if indexVar in dataArray:
            weight = df.loc[indexVar, 'weight']
            dataArray = dataArray.values*df.loc[indexVar, 'weight']
            print("-> %-30s | %4.1f" % (dataArray.name, weight))

    return dataArray
