# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 13:53:47 2023

@author: ram86
"""
import cv2
import numpy as np
import csv
import matplotlib.pyplot as plt
def imageAnalasis (filepath):
    
    img = cv2.imread(filepath)
    datalocale = filepath + ".csv"
    data = []
    with open(datalocale,'r') as file:
        csvreader = csv.reader(file)
        for row in csvreader:
            data.append(row)
    data = np.array(data)
    print(data)
    utm = data.astype(np.float32)[:11,:3]
    minim = np.min(utm,axis=0)
    utm = utm - np.min(utm, axis=0)
    print(utm)
    utm = [utm]
    rcpixels = data.astype(np.float32)[:,3:]
    pixels = [rcpixels[:11,::-1]]
    
    imgsize = img.shape[1::-1]
    print(utm)
    print(pixels)
    camera_matrix = cv2.initCameraMatrix2D(utm, pixels, imgsize)
    print(camera_matrix)
    ret, mtx, dist, rvecs, tvecs = cv2.calibrateCamera(utm, pixels, imgsize, camera_matrix, None, flags=cv2.CALIB_USE_INTRINSIC_GUESS)
    imgpoints, _ = cv2.projectPoints(utm[0], rvecs[0], tvecs[0], mtx, dist) 
    #creates blue points based on matrixes and not actual image
    imgpoints = imgpoints.reshape(-1,2)
    error = cv2.norm(pixels[0], imgpoints, cv2.NORM_L2)/len(imgpoints)
    print(error)
    plt.figure()
    plt.imshow(img[...,::-1])
    plt.plot(rcpixels[:,1], rcpixels[:,0], linestyle='none', marker='+', color='red')
    plt.plot(imgpoints[:,0], imgpoints[:,1], linestyle='none', marker='+', color='blue')
    plt.show()