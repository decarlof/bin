import numpy as np
import dxchange
import h5py
import argparse
import tomopy
import scipy.signal
import concurrent.futures as cf
import threading
from scipy import ndimage
from functools import partial
dbg = 0

def apply_shift(res, psi, p, id):
    """Apply shift for one projection."""
    # padding to avoid signal wrapping
    [nz,n] = psi[id].shape
    tmp = np.zeros([nz+nz//4, n+n//4], dtype='float32')
    tmp[nz//8:nz+nz//8, n//8:n+n//8] = psi[id]
    res0 = np.fft.ifft2(
        ndimage.fourier_shift(np.fft.fft2(tmp), p))
    res[id] = np.real(res0[nz//8:nz+nz//8, n//8:n+n//8])
    return res[id]

def apply_shift_batch(psi, p):
    """Apply shift for all projections in parallel"""
    res = np.zeros(psi.shape, dtype='float32')
    with cf.ThreadPoolExecutor() as e:
        shift = 0
        for res0 in e.map(partial(apply_shift, res, psi, p), range(0, psi.shape[0])):
            res[shift] = res0
            shift += 1
    return res

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "fname", help="Directory containing an input file name: /data/sample.h5")
    parser.add_argument(
        "foutname", help="Directory containing an output file name: /data/sample.h5")
    parser.add_argument("--axis", nargs='?', type=str, default="100", help="Approximate rotation axis location (pixel): 10.0 (default 10 image horizontal size)")
           
    args = parser.parse_args()
    apr_center = np.int(args.axis)

    # Read data
    proj, flat, dark, theta = dxchange.read_aps_32id(
        args.fname, sino=(0, 512),proj = (0,3000,1500))        
    print(theta[np.arange(0,3000,1500)])
    # filter data        
    data= proj
    data = tomopy.normalize(proj, flat, dark)            
    #data[data<0] = 0
    #data[np.isinf(data)] = 0
    #data[np.isnan(data)] = 0
    # remove stripes
    #data = tomopy.remove_stripe_fw(data,level=7,wname='sym16',sigma=1,pad=True)
    #data = tomopy.remove_stripe_sf(data, size=150)
    #data = tomopy.minus_log(data)        
    
    # stitched data            
    w = apr_center*2
    [ntheta,nz,n] = data.shape    
    datap1 = data[:ntheta//2]    
    datap2 = data[ntheta//2:,:,::-1]    
    if(dbg):
        center = 10
        datap1 = data[:ntheta//2,:,n//2-center:-center]    
        datap2 = data[:ntheta//2,:,center:n//2+center]    
        n=n//2
    shiftza = np.zeros(ntheta//2)
    shiftxa = np.zeros(ntheta//2)
    for k in np.arange(0,ntheta//2):
        # find location datap1 in datap2 by cross-correlation
        corr = scipy.signal.correlate2d(datap1[k,:,0:w]-np.mean(datap1[k,:,0:w]), datap2[k,:,-w:]-np.mean(datap2[k,:,-w:]), boundary='symm', mode='same')
        dxchange.write_tiff(datap1[k],'t1.tiff',overwrite=True)
        dxchange.write_tiff(datap2[k],'t2.tiff',overwrite=True)
        dxchange.write_tiff(datap1[k,:,0:w]-np.mean(datap1[k,:,0:w]),'tc1.tiff',overwrite=True)
        dxchange.write_tiff(datap2[k,:,-w:]-np.mean(datap2[k,:,-w:]),'tc2.tiff',overwrite=True)
        #dxchange.write_tiff(corr.astype('float32'),'t3.tiff',overwrite=True)
        shiftz, shiftx = np.unravel_index(np.argmax(corr),corr.shape)
        shiftza[k] = shiftz-nz/2+1
        shiftxa[k] = shiftx-w/2+1        
        print(k,shiftza[k],shiftxa[k])   
        
    shiftz = np.mean(shiftza)
    shiftx = np.mean(shiftxa)
    ishiftx = np.int(shiftx)
    fshiftx = shiftx - ishiftx
    
    # resulting data
    datanew = np.zeros([ntheta//2,nz,2*n],dtype='float32')
    # make smooth border between data sets
    fw1 = np.ones(n)
    fw2 = np.ones(n)
    fw1[0:ishiftx+w] = np.linspace(0,1,(ishiftx+w))
    fw2[-(ishiftx+w):] = np.linspace(1,0,(ishiftx+w))
    datanew[:,:,n:] = datap1*fw1
    datanew[:,:,ishiftx+w:n+ishiftx+w] = apply_shift_batch(datap2,[0,fshiftx])*fw1    
    dxchange.write_tiff_stack(datanew,'t/t.tiff')    

    fid = h5py.File(args.foutname,'w')
    fid.create_dataset('/exchange/data', (3000//2,1024,n*2), chunks=(3000//2,1,2*n),dtype='float32')
    fid.create_dataset('/exchange/data_white', (flat.shape[0],1024,n*2), chunks=(flat.shape[0],1,2*n),dtype='float32')
    fid.create_dataset('/exchange/data_dark', (dark.shape[0],1024,n*2), chunks=(dark.shape[0],1,2*n),dtype='float32')    
    for k in range(0,8):
        print(k)
        proj, flat, dark, theta = dxchange.read_aps_32id(
            args.fname, sino=(k*128, (k+1)*128))                
        [ntheta,nz,n] = proj.shape        
        nflat = flat.shape[0]                
        ndark = dark.shape[0]            
        projp1 = proj[:ntheta//2]    
        projp2 = proj[ntheta//2:,:,::-1]    
        flatp1 = flat
        flatp2 = flat
        darkp1 = dark
        darkp2 = dark
        if(dbg):
            center = 10
            projp1 = proj[:ntheta//2,:,n//2-center:-center]    
            projp2 = proj[:ntheta//2,:,center:n//2+center]                
            flatp1 = flat[:,:,n//2-center:-center]    
            flatp2 = flat[:,:,center:n//2+center]                
            darkp1 = dark[:,:,n//2-center:-center]    
            darkp2 = dark[:,:,center:n//2+center]                
            n=n//2           
        projnew = np.zeros([ntheta//2,nz,2*n],dtype='float32')        
        flatnew = np.zeros([nflat,nz,2*n],dtype='float32')        
        darknew = np.zeros([ndark,nz,2*n],dtype='float32')                
        # make smooth border between data sets
        fw1 = np.ones(n)
        fw2 = np.ones(n)
        fw1[0:2*ishiftx] = np.linspace(0,1,2*ishiftx)
        fw2[-2*ishiftx:] = np.linspace(1,0,2*ishiftx)
        projnew[:,:,n-ishiftx:-ishiftx] = fw1*projp1
        projnew[:,:,ishiftx:n+ishiftx] += fw2*apply_shift_batch(projp2,[0,fshiftx])    
        flatnew[:,:,n-ishiftx:-ishiftx] = fw1*flatp1
        flatnew[:,:,ishiftx:n+ishiftx] += fw2*apply_shift_batch(flatp2,[0,fshiftx])    
        darknew[:,:,n-ishiftx:-ishiftx] = fw1*darkp1
        darknew[:,:,ishiftx:n+ishiftx] += fw2*apply_shift_batch(darkp2,[0,fshiftx])    
        # save into new h5
        fid['/exchange/data'][:,k*128:(k+1)*128] = projnew
        fid['/exchange/data_white'][:,k*128:(k+1)*128] = flatnew
        fid['/exchange/data_dark'][:,k*128:(k+1)*128] = darknew        
        