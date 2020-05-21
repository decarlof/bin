import numpy as np
import matplotlib.pyplot as plt
from tomopy import recon, circ_mask
import h5py
import dxchange
import tomopy 

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

def find_center(
        tomo, theta, ind=None, init=None,
        tol=0.5, mask=True, ratio=1., sinogram_order=False):
    """
    Find rotation axis location.

    The function exploits systematic artifacts in reconstructed images
    due to shifts in the rotation center. It uses image entropy
    as the error metric and ''Nelder-Mead'' routine (of the scipy
    optimization module) as the optimizer :cite:`Donath:06`.

    Parameters
    ----------
    tomo : ndarray
        3D tomographic data.
    theta : array
        Projection angles in radian.
    ind : int, optional
        Index of the slice to be used for reconstruction.
    init : float
        Initial guess for the center.
    tol : scalar
        Desired sub-pixel accuracy.
    mask : bool, optional
        If ``True``, apply a circular mask to the reconstructed image to
        limit the analysis into a circular region.
    ratio : float, optional
        The ratio of the radius of the circular mask to the edge of the
        reconstructed image.
    sinogram_order: bool, optional
        Determins whether data is a stack of sinograms (True, y-axis first axis)
        or a stack of radiographs (False, theta first axis).

    Returns
    -------
    float
        Rotation axis location.
    """
    tomo = dtype.as_float32(tomo)
    theta = dtype.as_float32(theta)

    if sinogram_order:
        dy, dt, dx = tomo.shape
    else:
        dt, dy, dx = tomo.shape

    if ind is None:
        ind = dy // 2
    if init is None:
        init = dx // 2

    # extract slice we are using to find center
    if sinogram_order:
        tomo_ind = tomo[ind:ind + 1]
    else:
        tomo_ind = tomo[:, ind:ind + 1, :]

    hmin, hmax = _adjust_hist_limits(
        tomo_ind, theta, mask, sinogram_order)

    # Magic is ready to happen...
    res = minimize(
        _find_center_cost, init,
        args=(tomo_ind, theta, hmin, hmax, mask, ratio, sinogram_order),
        method='Nelder-Mead',
        tol=tol)
    return res.x

def _find_center_cost(
        center, tomo_ind, theta, hmin, hmax, mask, ratio,
        sinogram_order=False):
    """
    Cost function used for the ``find_center`` routine.
    """
    rec = recon(
        tomo_ind, theta, center,
        sinogram_order=sinogram_order,
        algorithm='gridrec',
        filter_name='shepp')

    if mask is True:
        rec = circ_mask(rec, axis=0)

    hist, e = np.histogram(rec, bins=64, range=[hmin, hmax])
    hist = hist.astype('float32') / rec.size + 1e-12
    val = -np.dot(hist, np.log2(hist))
    # from Matt
    bval = -((rec - rec.mean())**2).sum()

    return val, bval

def _adjust_hist_min(val):
    if val < 0:
        val = 2 * val
    elif val >= 0:
        val = 0.5 * val
    return val

def _adjust_hist_max(val):
    if val < 0:
        val = 0.5 * val
    elif val >= 0:
        val = 2 * val
    return val

def _adjust_hist_limits(tomo_ind, theta, mask, sinogram_order):
    # Make an initial reconstruction to adjust histogram limits.
    rec = recon(tomo_ind,
                theta,
                sinogram_order=sinogram_order,
                algorithm='gridrec')

    # Apply circular mask.
    if mask is True:
        rec = circ_mask(rec, axis=0)

    # Adjust histogram boundaries according to reconstruction.
    return _adjust_hist_min(rec.min()), _adjust_hist_max(rec.max())

def rescale(s):
    s = np.array(s)
    return (s-s.min()) / (s.max() - s.min())

def main():

    # file_name = '/local/data/2020-02/Stock/100_B949_81_84_B2.h5'
    file_name = '/local/data/tomo_00001.h5'

    data_size = get_dx_dims(file_name)

    ssino = int(data_size[1] * 0.5)
    detector_center = int(data_size[2] * 0.5)

    # Select sinogram range to reconstruct
    sino_start = ssino
    sino_end = sino_start + 1

    sino = (int(sino_start), int(sino_end))


    # Read APS 2-BM raw data
    proj, flat, dark, theta = dxchange.read_aps_32id(file_name, sino=sino)

    tomo_ind = tomopy.normalize(proj, flat, dark)

    # data = tomopy.normalize_bg(proj, air=10)
    tomo_ind = tomopy.minus_log(tomo_ind)


    rec = recon(tomo_ind, theta, center=detector_center, sinogram_order=False,
                algorithm='gridrec', filter_name='shepp')
    rec = circ_mask(rec, axis=0)

    # tomopy score, simplified rescaling of histogram bounds 
    hmin, hmax = _adjust_hist_limits(
        tomo_ind, theta, mask=True, sinogram_order=False)

    print(hmin, hmax)

    centers = np.linspace(-25, 25, 101) + detector_center
    print(centers)

    tpscore = []
    blur = []
    for center in centers:
        val, bval = _find_center_cost(
                    center, tomo_ind, theta, hmin, hmax, mask=True, ratio=1.,
                    sinogram_order=False)
        tpscore.append(val)
        blur.append(bval)

    plt.plot(centers, rescale(tpscore), label='tomopy score')
    plt.plot(centers, rescale(blur), label='blurriness')
    plt.legend()
    plt.title('centering scores')
    plt.xlabel('center')
    plt.show()

if __name__ == "__main__":
    main()