# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 11:18:32 2023

@author: ram86
"""
from magnetic_code import coordinateChange, IGRFMaster
from magnetic_code import image_anal as img
from magnetic_code import error_detect as err
import matplotlib.pyplot as plt
import numpy as np

def cellMagneticVectors(lat, lon, alt, date, rolli, pitch, yaw, eyecells = 1000000, isdeg = True):
    '''
    

    Parameters
    ----------
    lat : FLOAT
        Latitude of the location using the WGS-84 elipsoid in degrees.
    lon : FLOAT
        Longitude of the location using the WGS-84 elipsoid in degrees.
    alt : FLOAT
        Altitude of the point off the ground in km.
    date : FLOAT
        Date of intrest between 1900-2030 where the day of the year is the fraction portion.
    rolli : FLOAT
        roll of the bird relative to the horison.
    yaw : FLOAT
        yawing of bird from 'true' geographic north. Anticlockwise
    pitch : FLOAT
        angle the bird makes with the vertical (90+angle off horizon).
    eyecells : INT, optional
        The number of retinal cells you wish to model
    isdeg : BOOL, optional
        Whether the value output are in degrees(True) or radians(False). The default is True.
    
    Returns
    -------
    ix: np.array
        The x coordinate in the image plane
    iy: np.array
        The y coordinate in the image plane
    legendre : np.array
        The normalisesd yield based off of the spherical harmonic function
    heading : FLOAT
        The heading of the bird off of geographic north

    '''
    ## Conversion from pitch being off of horizon to being from z axis
    ## A roll of 0 leads to y pointing down so it needs to be modified by 180
    roll = rolli + 0
    roll = roll % 360
    ## This function gets the magnetic field based on the lat, long, altitude and date
    eff, hoz, dec, inc, mX, mY, mZ = IGRFMaster.getIGRF(lon, lat, alt, date)
    ## Changes coordinates to East->X, North -> Y, Up -> Z
    mX, mY, mZ = coordinateChange.zDownToZUp(mX, mY, mZ)
    ## Goes from World Coordinates to Camera Coordinates
    mX, mY, mZ,heading = coordinateChange.eulerChangeToCamera(mX, mY, mZ, roll, pitch, yaw)
    ## Creates an eye with the given number of cells 
    xc, yc, zc = coordinateChange.createEye(eyecells)
    ## Finds the location on the pixel plane that the eye cell goes to projected from pupil
    iy,ix = coordinateChange.findISpace(xc, yc, zc, 1)
    ## Calculates the angle between the vector to the cell and the magnetic field
    abtwl  = coordinateChange.cartAngleBetween(xc,yc,zc,mX,mY,mZ)
    ## Calculates legendre spherical harmonics and normalises it based on given weightings
    legendre = coordinateChange.legendre(abtwl)
    legendre = np.real(legendre)
    legendre = (legendre - legendre.min())/(legendre.max() - legendre.min())*255
    plt.figure(figsize = (7,7))
    plt.scatter(ix,iy,c=legendre)
    plt.xlabel('ix')
    plt.ylabel('iy')
    plt.title(f'roll = {rolli}, angle off horison = {pitch-90}, heading = {heading}\n{lat},{lon}')
    plt.colorbar()
    coords = []
    iy = -iy
    size = 2700
    ix = (ix - ix.min())/(ix.max() - ix.min())*size
    iy = (iy - iy.min())/(iy.max() - iy.min())*size
    ix = ix.astype('int')
    iy = iy.astype('int')

    print (legendre)
    
    return ix,iy,legendre,heading







#                cellMagneticVectors(latitude,longitude,altitude(km),date,  roll, pitch, yaw)
ix1, iy1, col1,heading = cellMagneticVectors(50.734531, -3.53, .2, 2023+((31+28)/365), 0, 90, 0, 1000000)
