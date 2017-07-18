"""
Functions for input and output
"""

from astropy import table, units
from astroquery.simbad import Simbad
from astropy.io.ascii import read as read_ascii
from astropy.io import votable
import os
from glob import glob
from sys import stdout
from numpy import save, load

def clear():
    _ = os.system("clear"); del _

def flush():
    stdout.flush()

def fetch_v_r_from_simbad(object_name):
    """
    Fetch an observed radial velocity from SIMBAD
    
    Parameters
    ----------
    object_name: str
        Name of the object to query SIMBAD for

    Returns
    -------
    v_r: float
        Observed radial velocity, in km/s
    v_r_err: float
        Error in observed radial velocity, in km/s
    """
    try:
        Simbad.add_votable_fields("rvz_radvel", "rvz_error")
    except KeyError:
        pass # already present in query
    res = Simbad.query_object(object_name)
    if res is None:
        raise Exception("Cannot fetch radial velocity from SIMBAD for object {0}".format(object_name))
    v_r, v_r_err = res["RVZ_RADVEL"][0], res["RVZ_ERROR"][0]
    return v_r, v_r_err

def read_votable(*args, **kwargs):
    data = votable.parse(*args, **kwargs)
    t = data.get_first_table().to_table()
    return t

def read_csv(*args, **kwargs):
    t = read_ascii(*args, format="csv", **kwargs)
    return t

def find_array_columns(t):
    return [col for col in t.keys() if t[col].shape[1:]]
def remove_array_columns_from_table(t):
    new_t = t.copy()
    remove_columns = find_array_columns(t)
    new_t.remove_columns(remove_columns)
    return new_t

def write_table_without_arrays(t, saveto, format="ascii.fast_csv", exclude_cols=[], *args, **kwargs):
    new_t = t.copy()
    assert all(col in new_t.keys() for col in exclude_cols), "gaia_fc.general.write_table_without_arrays: columns {0} cannot be excluded because they are not in the table".format([col for col in exclude_cols if col not in new_t.keys()])
    new_t.remove_columns(exclude_cols)
    new_t = remove_array_columns_from_table(new_t)
    new_t.write(saveto, format=format, *args, **kwargs)

def write_table_with_separate_arrays(t, saveto_folder, format="ascii.fast_csv", exclude_cols=[], verbose = True, *args, **kwargs):
    t_no_arr = remove_array_columns_from_table(t)
    arrays   = find_array_columns(t)
    assert "table" not in arrays, "gaia_fc.io.write_table_with_separate_arrays: There is a column named `table` in your table, which cannot be written out"
    if not os.path.isdir(saveto_folder):
        os.mkdir(saveto_folder)
    if not os.path.isdir(saveto_folder + "/np/"):
        os.mkdir(saveto_folder + "/np/")
    for arr_key in arrays:
        save("{0}/np/{1}.npy".format(saveto_folder, arr_key), t[arr_key]._data)
        if verbose:
            print arr_key,
            flush()
    t_no_arr.write("{0}/table.csv".format(saveto_folder), format=format, overwrite=True, *args, **kwargs)
    if verbose:
        print ""

def load_table_with_separate_arrays(saveto_folder, format="csv", verbose = True, *args, **kwargs):
    t = read_ascii("{0}/table.dat".format(saveto_folder), format=format, *args, **kwargs)
    array_files = glob("{0}/np/*.npy".format(saveto_folder))
    arrays = [load(f) for f in array_files]
    for f, arr in zip(array_files, arrays):
        name = os.path.splitext(os.path.basename(f))[0]
        as_column = table.Column(data = arr, name = name)
        t.add_column(as_column)
        if verbose:
            print name,
            flush()
    if verbose:
        print ""
    return t

def save_PDFs(amplitudes, means, covariances, saveto_folder, *args, **kwargs):
    if not os.path.isdir(saveto_folder):
        os.mkdir(saveto_folder)
    save("{0}/amplitudes.npy".format(saveto_folder), amplitudes, *args, **kwargs)
    save("{0}/means.npy".format(saveto_folder), means, *args, **kwargs)
    save("{0}/covariances.npy".format(saveto_folder), covariances, *args, **kwargs)

def load_PDFs(folder, *args, **kwargs):
    amplitudes = load("{0}/amplitudes.npy".format(folder), *args, **kwargs)
    means = load("{0}/means.npy".format(folder), *args, **kwargs)
    covariances = load("{0}/covariances.npy".format(folder), *args, **kwargs)
    return amplitudes, means, covariances
