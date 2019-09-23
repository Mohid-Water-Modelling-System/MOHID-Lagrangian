#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 18:27:19 2019

@author: daniel
"""
import os
import xarray as xr


def MOHIDCurvilinearCorrect(MOHID_file, outdir, modified_suffix='_m'):
    # set decode_cf to False
    # Value error related to multiple _fill_value appears.
    ds = xr.open_dataset(MOHID_file,decode_cf=False)
    ds['line_c'] = ds.lat[:,0]
    ds['column_c'] = ds.lon[0,:]
    ds=ds.drop(['lat','lon'])
    #ds=ds.drop(['lat_staggered','lon_staggered'])
    ds=ds.rename({'column_c':'lon','line_c':'lat'})
    ds.close()
    file = os.path.basename(MOHID_file)
    ds.to_netcdf(outdir+'/'+file)
    #ds.to_netcdf(MOHID_file.replace('.nc', modified_suffix + '.nc'))
    # We must write with a different name or the nc file doesn't write.
    return