import numpy as np
import matplotlib.pyplot as plt
from tomopy import recon, circ_mask
import h5py
import dxchange
import tomopy 
import time
import os

def get_dx_dims(full_file_name):
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

    with h5py.File(full_file_name, "r") as f:
        try:
            data = f[grp]
        except KeyError:
            return None

        shape = data.shape

    return shape

def main():

    data_top = '/local/data/2020-02/Stock/'
    file_name = '099_B949_81_84_B2'
    # data_top = '/local/data/'
    # file_name = 'tomo_00001'
    top = '/local/data/2020-02/Stock/'

    # log.info(os.listdir(top))
    h5_file_list = list(filter(lambda x: x.endswith(('.h5', '.hdf')), os.listdir(top)))
    h5_file_list.sort()

    print("Found: %s" % h5_file_list)

    for fname in h5_file_list:

        full_file_name = data_top + fname
        data_size = get_dx_dims(full_file_name)

        print(data_size)
        ssino = int(data_size[1] * 0.5)
        detector_center = int(data_size[2] * 0.5)

        # Select sinogram range to reconstruct
        sino_start = ssino
        sino_end = sino_start + 10

        sino = (int(sino_start), int(sino_end))


        # Read APS 2-BM raw data
        proj, flat, dark, theta = dxchange.read_aps_32id(full_file_name, sino=sino)

        tomo_ind = tomopy.normalize(proj, flat, dark)
        print(os.path.splitext(full_file_name)[0]+'_sino')
        dxchange.write_tiff_stack(tomo_ind,fname=os.path.splitext(full_file_name)[0]+'_sino', axis=1)
    
        # rec = recon(tomo_ind, theta, center=detector_center, sinogram_order=False,
        #             algorithm='gridrec', filter_name='None')
        # rec = circ_mask(rec, axis=0)
        # dxchange.write_tiff_stack(rec,fname=data_top+file_name+'_rec')

if __name__ == "__main__":
    main()