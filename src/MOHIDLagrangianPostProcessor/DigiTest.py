#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 19:27:27 2019

Flat[x + WIDTH * (y + DEPTH * z)] = Original[x, y, z]

@author: daniel
"""

import numpy as np

def to1D(ridx, ni,nj,nk):
    return (ridx[:,2]* ni * nj) + (ridx[:,1] * ni) + ridx[:,0]


def to3D(idx, ni,nj,nk):
    k = np.int32(idx / (ni * nj))
    idx = idx - (k * ni * nj)
    j = np.int32(idx / ni)
    i = idx % ni
    return  np.c_[i, j, k]


r = np.random.random((10000,3))
t = np.arange(10000)

x_edges = np.linspace(0,1,10)
y_edges = np.linspace(0,1,10)
z_edges = np.linspace(0,1,10)

nx = x_edges.size
ny = y_edges.size
nz = z_edges.size

x_dig = np.digitize(r[:,0],x_edges)
y_dig = np.digitize(r[:,1],y_edges)
z_dig = np.digitize(r[:,2],z_edges)

r_d = np.c_[x_dig,y_dig,z_dig]

r_idx = to1D(r_d,nx,ny,nz)
p_idx = np.arange(0,)
n_counts = np.zeros(nx*ny*nz)

for i in range(0,nx*ny*nz):
    n_counts[i] = np.mean(var[i==r_idx])
    
    
    
    

np.unique(r_d,axis=0,return_counts=True)

