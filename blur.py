#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions useful for fly-scan
"""

from __future__ import print_function

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt


def restricted_float(x):

    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x

def calc_blur_error(exposure_time, rotation_speed, camera_size_x):

    blur_delta = exposure_time * rotation_speed
    blur_pixel = (camera_size_x / 2.0) * np.sin(blur_delta * np.pi /180.)

    return blur_pixel



def main(arg):


    parser = argparse.ArgumentParser()
    parser.add_argument("--exposure_time", nargs='?', type=float, default=0.1, help=" ")
    parser.add_argument("--readout_time", nargs='?', type=float, default=0.006, help="flir options: 8bit: 0.006; 16-bit: 0.01")
    parser.add_argument("--camera_size_x", nargs='?', type=int, default=2448, help="flir option: 2448")
    parser.add_argument("--num_projections", nargs='?', type=int, default=1500, help=" ")
    parser.add_argument("--rotation_slow_factor", nargs='?', type=restricted_float, default=1.0, help="rotation speed reduction factor from the max allowed")

    args = parser.parse_args()

    angular_range = 180.0         # deg
    angular_step = angular_range / args.num_projections

    min_scan_time = args.num_projections * (args.exposure_time + args.readout_time)
    max_rot_speed = angular_range / min_scan_time

    rotation_speed = max_rot_speed * args.rotation_slow_factor
    scan_time = angular_range / rotation_speed

    frame_rate = args.num_projections / scan_time

    blur = calc_blur_error(args.exposure_time, rotation_speed, args.camera_size_x)
    print('Scan time: %5.2f s, Frame rate: %5.2f fps,  Rot speed: %5.2f deg/s, Blur: %5.2f pixel' % (scan_time, frame_rate, rotation_speed, blur))



# def set_acquisition(blur_pixel, exposure_time, readout_time, camera_size_x, angular_range, num_projections):

#     """
#     Calculate frame rate and rotation speed for a desired blur error t

#     Parameters
#     ----------
#     blur_pixel : float
#         Desired blur error. For good quality reconstruction this should be < 0.2 pixel.
#     exposure_time: float
#         Detector exposure time
#     readout_time : float
#         Detector read out time
#     camera_size_x : int
#         Detector X size
#     angular_range : float
#         Tomographic scan angular range
#     num_projections : int
#         Numember of projections

#     Returns
#     -------
#     float
#         frame_rate, rot_speed
#     """

#     delta_blur  = np.arccos(((camera_size_x / 2.0) - blur_pixel) / (camera_size_x / 2.0)) * 180.0 / np.pi
#     rot_speed = delta_blur  / exposure_time

#     scan_time = angular_range / rot_speed
#     frame_rate = num_projections / scan_time
#     print("*************************************")
#     print("Total # of proj: ", num_projections)
#     print("Exposure Time: ", exposure_time, "s")
#     print("Readout Time: ", readout_time, "s")
#     print("Angular Range: ", angular_range, "degrees")
#     print("Camera X size: ", camera_size_x)
#     print("Blur Error: ", blur_pixel, "pixels")
#     print("*************************************")
#     print("Rot Speed: : ", rot_speed, "degrees/s")
#     print("Scan Time:: ", scan_time, "s")
#     print("Frame Rate: ", frame_rate, "fps")
#     print("*************************************")
  
#     return frame_rate, rot_speed

   
if __name__ == "__main__":
    main(sys.argv[1:])

 