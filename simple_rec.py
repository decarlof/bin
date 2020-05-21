import numpy as np
import matplotlib.pyplot as plt
from tomopy import recon, circ_mask
import h5py
import dxchange
import tomopy 
import time

def get_dx_dims(file_name):
    """
    Read array size of a specific group of Data Exchange file.
    Parameters
    ----------
    fname : str
        String defining the path of file or file name.
    dataset : str
        Path to the dataset inside hdf5 file where data is located.
    Returns
    -------
    ndarray
        Data set size.
    """

    dataset='data'

    grp = '/'.join(['exchange', dataset])

    with h5py.File(file_name, "r") as f:
        try:
            data = f[grp]
        except KeyError:
            return None

        shape = data.shape

    return shape

def binning(proj, flat, dark, binning):

    print('  *** *** binning: %d' % binning)
    proj = _binning(proj, binning)
    flat = _binning(flat, binning)
    dark = _binning(dark, binning)

    return proj, flat, dark

def _binning(data, binning):

    data = tomopy.downsample(data, level=int(binning), axis=2) 
    data = tomopy.downsample(data, level=int(binning), axis=1)

    return data

def main():

    file_name = '/local/data/2020-02/Stock/099_B949_81_84_B2.h5'

    tic_01 =  time.time()

    data_size = get_dx_dims(file_name)

    print(data_size)
    ssino = int(data_size[1] * 0.5)
    detector_center = int(data_size[2] * 0.5)

    # Select sinogram range to reconstruct
    sino_start = ssino
    sino_end = sino_start + 1

    sino = (int(sino_start), int(sino_end))


    # Read APS 2-BM raw data
    projf, flatf, darkf, theta = dxchange.read_aps_32id(file_name, sino=sino)
    proj, flat, dark = binning(projf, flatf, darkf, 2)

    tomo_ind = tomopy.normalize(proj, flat, dark)

    # data = tomopy.normalize_bg(proj, air=10)
    tomo_ind = tomopy.minus_log(tomo_ind)
    tic_02 =  time.time()
    print('prep time', tic_02-tic_01)
    rec = recon(tomo_ind, theta, center=detector_center, sinogram_order=True,
                algorithm='gridrec', filter_name='None')
    rec = circ_mask(rec, axis=0)
    tic_03 =  time.time()
    print('recon time:', tic_03-tic_02)
if __name__ == "__main__":
    main()