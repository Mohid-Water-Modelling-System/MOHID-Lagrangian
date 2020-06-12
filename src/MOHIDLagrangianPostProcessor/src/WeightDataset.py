# -*- coding: utf-8 -*-


import pandas as pd
import xarray as xr

def WeightDataset(dataset_file, weight_file):

    df = pd.read_csv(weight_file)
    ds = xr.open_dataset(dataset_file)

    for datasetVar in list(ds.keys()):
        for indexVar in df.index:
            if indexVar in datasetVar:
                ds[datasetVar] = ds[datasetVar]*df.loc[indexVar, 'weight']

    return ds
