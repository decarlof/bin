from epics import PV
import time
import numpy as np
import matplotlib.pyplot as plt


def take_image():

    print('  *** taking a single image')
   
    nRow = PV('2bmbSP1:cam1:SizeY_RBV').get()
    nCol = PV('2bmbSP1:cam1:SizeX_RBV').get()
    pixelFormat = PV('2bmbSP1:cam1:PixelFormat_RBV').get(as_string=True)

    if (pixelFormat == "Mono16"):
        pixel_f = 16
    elif (pixelFormat == "Mono8"):
        pixel_f =8
    else:
        print('error: bit format not supported')
        exit()

    print(pixelFormat, pixel_f)
    image_size = nRow * nCol

    PV('2bmbSP1:cam1:Acquire').put(1, wait=True, timeout=1000.0)
    time.sleep(1)

    # Get the image loaded in memory
    img_vect = PV('2bmbSP1:image1:ArrayData').get(count=image_size)

    img = np.reshape(img_vect,[nRow, nCol])
    img_uint = np.mod(img, 2**pixel_f)

    imgplot = plt.imshow(img_uint)
    plt.show()
    return img_uint

take_image()