"""
Functions for mapping other functions over (large) data sets, including multiprocessing capabilities.
"""

import numpy as np

from parmap import map as pmap, starmap as smap
from multiprocessing import cpu_count

from astropy.table.column import MaskedColumn, Column
from astropy.table.table import Table

def map_np(function, *args, **kwargs):
    return np.array(map(function, *args, **kwargs))

def pmap_np(function, iterable, *args, **kwargs):
    return np.array(pmap(function, iterable, *args, processes = cpu_count(), **kwargs))

def smap_np(function, iterables, *args, **kwargs):
    return np.array(smap(function, iterables, *args, processes = cpu_count(), **kwargs))

def bin_groups(bo_orig, min_in_group, minrange):
    bo = bo_orig.copy()
    bo.sort()
    current_right_edge = 0
    bin_left_edges = [bo[0]]
    while len(bo[current_right_edge:]) >= min_in_group and np.abs(bo[current_right_edge] - bo[-1]) >= minrange:
        within_range = np.where(bo[current_right_edge+1:] - bo[current_right_edge] <= minrange)[0]
        if len(within_range) <= min_in_group:
            current_right_edge += min_in_group
        else:
            current_right_edge += len(within_range)
        same_values = np.where(bo[current_right_edge+1:] == bo[current_right_edge])[0]
        current_right_edge += len(same_values)+1
        bin_left_edges.append(bo[current_right_edge])
    if len(np.where(bo > bin_left_edges[-1])[0]) < min_in_group:
        _ = bin_left_edges.pop()
    return bin_left_edges

def map_over_bins(function, binover, bins, mapover, linlog="lin", range_mapover=None, range_binover=None, *args, **kwargs):
    """
    Bin data and apply a function to the bins
    """
    if range_mapover is not None:
        indices_range = np.where((mapover >= range_mapover[0]) & (mapover < range_mapover[1]))[0]
        mapover = mapover[indices_range]
        binover = binover[indices_range]
    if range_binover is not None:
        indices_range = np.where((binover >= range_binover[0]) & (binover < range_binover[1]))[0]
        mapover = mapover[indices_range]
        binover = binover[indices_range]
    if isinstance(binover, MaskedColumn):
        bo = binover.copy()
        mapover = np.array(mapover[-bo.mask])
        bo = np.array(bo[-bo.mask])
    elif isinstance(binover, Table):
        bo = binover.copy()
        if len(bo.colnames) == 1:
            bo = bo.columns[0]
            mapover = np.array(mapover[-bo.mask])
            bo = np.array(bo[-bo.mask])
        else:
            mask = np.array([any(row) for row in bo.mask])
            mapover = np.array(mapover[-mask])
            bo = bo[-mask]
            bo = np.vstack(np.array(bo[key]) for key in bo.columns)
    elif isinstance(binover, Column):
        bo = np.array(binover.copy())
        mapover = np.array(mapover)
    elif isinstance(binover, np.ndarray):
        bo = binover.copy()
        mapover = np.array(mapover)
    elif isinstance(binover, list):
        bo = np.array(binover)
        mapover = np.array(mapover)
    else:
        raise NotImplementedError("Parameter binover cannot be of type {0}".format(type(binover)))

    assert len(mapover) == len(bo.T), "mapover has length {0} | bo has length {1}".format(len(mapover), len(bo.T))

    if mapover.ndim > 2:
        raise NotImplementedError("Parameter mapover cannot be of a dimension greater than 2; is {0}".format(mapover.ndim))

    if bo.ndim == 1:
        if isinstance(bins, int):
            if linlog == "lin":
                bin_left_edges = np.linspace(bo.min(), bo.max(), bins)
            elif linlog == "log":
                bin_left_edges = np.logspace(bo.min(), bo.max(), bins)
            else:
                raise ValueError("Parameter linlog can only be `lin` or `log`, not {0}".format(linlog))
        elif isinstance(bins, list):
            bin_left_edges = bins
        elif isinstance(bins, tuple):
            bin_left_edges = bin_groups(bo, *bins)
        else:
            raise NotImplementedError("Parameter bins cannot be of type {0}".format(type(bins)))
        bin_left_edges = np.array(bin_left_edges)
        digitized = np.digitize(bo, bin_left_edges)-1
        grouped_indices = [np.where(digitized == i)[0] for i in range(len(bin_left_edges))]
        x_sigma = [np.std(binover[indices])/np.sqrt(len(binover[indices])) for indices in grouped_indices]
        mapover_rows = [mapover[indices] for indices in grouped_indices]
        try:
            functions = iter(function)
            result = [map_np(f, mapover_rows, *args, **kwargs) for f in functions]
        except TypeError: # only one function given
            result = map_np(function, mapover_rows, *args, **kwargs)
        bin_left_edges_list = bin_left_edges
        highest_value = bo.max()
    elif bo.ndim == 2:
        if isinstance(bins, int):
            if linlog == "lin":
                bin_left_edges_list = [np.linspace(bo_col.min(), bo_col.max(), bins) for bo_col in bo]
            elif linlog == "log":
                bin_left_edges_list = [np.logspace(bo_col.min(), bo_col.max(), bins) for bo_col in bo]
            else:
                raise ValueError("Parameter linlog can only be `lin` or `log`, not {0}".format(linlog))
        elif isinstance(bins, tuple):
            if linlog == "lin":
                bin_left_edges_list = [np.linspace(bo_col.min(), bo_col.max(), bin_nr) for bo_col,bin_nr in zip(bo,bins)]
            elif linlog == "log":
                bin_left_edges_list = [np.logspace(bo_col.min(), bo_col.max(), bin_nr) for bo_col,bin_nr in zip(bo,bins)]
            else:
                raise ValueError("Parameter linlog can only be `lin` or `log`, not {0}".format(linlog))
        else:
            raise NotImplementedError("Parameter bins cannot be of type {0}".format(type(bins)))
        bin_left_edges_list = [np.array(bin_left_edges) for bin_left_edges in bin_left_edges_list]
        digitized_list = [np.digitize(bo_col, bin_left_edges)-1 for bo_col, bin_left_edges in zip(bo, bin_left_edges_list)]
        if len(bo) == 2:
            grouped_indices = [[np.where(np.logical_and(digitized_list[0] == i, digitized_list[1] == j))[0] for i in range(len(bin_left_edges_list[0]))] for j in range(len(bin_left_edges_list[1]))]
            highest_value = (bo[0].max(), bo[1].max())
        else:
            raise NotImplementedError("Cannot map 2-dimensional data with {0} columns".format(len(bo)))
        mapover_rows = [[mapover[indices] for indices in grouped_indices_X] for grouped_indices_X in grouped_indices]
        x_sigma = [[np.std(binover[indices])/np.sqrt(len(binover[indices])) for indices in grouped_indices_X] for grouped_indices_X in grouped_indices]
        try:
            functions = iter(function)
            result = [np.array([[f(mr, *args, **kwargs) for mr in mapover_rows_X] for mapover_rows_X in mapover_rows]) for f in functions]
        except TypeError: # only one function given
            result = np.array([[function(mr, *args, **kwargs) for mr in mapover_rows_X] for mapover_rows_X in mapover_rows])
    else:
        raise NotImplementedError("Cannot map over bins for data of dimension {0}".format(bo.ndim))

    return result, bin_left_edges_list, highest_value, x_sigma
