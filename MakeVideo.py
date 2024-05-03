#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
'''
                            MAKEVIDEO.PY
    Independent script that you should run once main has already been executed.
    It makes a video with the different plots of the simulation.

'''

@author: lucia
"""

import cv2
import os
    
image_folder = "Figures"

video_name = "Video.avi"

def myFunc(e):
  return int(e.split(".")[0])

images = [img for img in os.listdir(image_folder) if img.endswith(".png")]

images.sort(key=myFunc) 

frame = cv2.imread(os.path.join(image_folder, images[0]))
height, width, layers = frame.shape
video = cv2.VideoWriter(video_name, 0, 15, (width,height))

for image in images:
    video.write(cv2.imread(os.path.join(image_folder, image)))

cv2.destroyAllWindows()
video.release()


