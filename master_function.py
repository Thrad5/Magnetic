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

def cellMagneticVectors(lat, lon, alt, date, rolli, pitch, yaw, eyecells = 600000, w2 = 1, w4 = 0, w6 = 0, w8 = 0, haveimg = False, imgpth = None):
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
    yaw : FLOAT
        yawing of bird from 'true' geographic north. Anticlockwise
    pitch : FLOAT
        angle the bird makes with the vertical (90+angle off horizon).

    Returns
    -------
    None.

    '''
    debug = True
    ## Conversion from pitch being off of horizon to being from z axis
    ## A roll of 0 leads to y pointing down so it needs to be modified by 180
    roll = rolli + 0
    roll = roll % 360
    ## This function gets the magnetic field based on the lat, long, altitude and date
    eff, hoz, dec, inc, mX, mY, mZ = IGRFMaster.getIGRF(lon, lat, alt, date)
    ## Changes coordinates to East->X, North -> Y, Up -> Z
    mX, mY, mZ = coordinateChange.zDownToZUp(mX, mY, mZ)
    # print(mX, mY, mZ)
    if debug == True:
        mX = 0
        mY = 1
        mZ = 0
    ## Goes from World Coordinates to Camera Coordinates
    print(yaw)
    mX, mY, mZ,heading = coordinateChange.eulerChangeToCamera(mX, mY, mZ, roll, pitch, yaw)
    # print(mX, mY, mZ)
    ## Finds the spherical coordinates of the Magnetic field
    if haveimg == False:
        ## Creates an eye with the given number of cells 
        xc, yc, zc, cphi, ctheta,lr, lphi, ltheta = coordinateChange.createEye(eyecells)
        ## Finds the location on the pixel plane that the eye cell goes to projected from pupil
        iy,ix = coordinateChange.findISpace(xc, yc, zc, 1)
        # print('ix min and max')
        # print(ix.min(),ix.max())
        # print('iy min and max')
        # print(yc.min(),yc.max())
    else:
        print('Need to implement')
    ## Calculates the angle between the vector to the cell and the magnetic field
    # abtwl = coordinateChange.angleBetween(ctheta, cphi, thetam, phim)
    abtwl  = coordinateChange.cartAngleBetween(xc,yc,zc,mX,mY,mZ)
    ## Calculates legendre spherical harmonics and normalises it based on given weightings
    legendre = coordinateChange.legendre(abtwl,w2,w4,w6,w8)
    legendre = np.real(legendre)
    legendre = (legendre - legendre.min())/(legendre.max() - legendre.min())*255
    plt.figure(figsize = (7,7))
    plt.scatter(ix,iy,c=legendre)
    plt.xlabel('ix')
    plt.ylabel('iy')
    plt.title(f'L = [{w2,w4,w6,w8}], roll = {rolli}, pitch = {pitch}, yaw = {yaw}')
    # plt.colorbar()
    # print(len(lphi),len(ltheta))
    # print('ltheta max, ltheta min')
    # print(ltheta.max(),ltheta.min())
    # print('lphi min, lphi max')
    # print(lphi.min(),lphi.max())
    coords = []
    iy = -iy
    size = 2700
    ix = (ix - ix.min())/(ix.max() - ix.min())*size
    iy = (iy - iy.min())/(iy.max() - iy.min())*size
    ix = ix.astype('int')
    iy = iy.astype('int')
    legendre = legendre.astype('int')
    print (legendre)
    for i in range(len(xc)):
        coords.append((ix[i],iy[i],xc[i],yc[i],zc[i],legendre[i]))
    x_coords = sorted(list(set(i[0] for i in coords)))
    y_coords = sorted(list(set(i[1] for i in coords)))
    result = np.zeros([len(y_coords),len(x_coords)])
    xsaty0 = []
    ysatx0 = []
    for i in coords:
        x_final = x_coords.index(i[0])
        y_final = y_coords.index(i[1])
        result[y_final][x_final] = i[5]
        if i[1] == size/2:
            xsaty0.append(i[0])
        if i[0] == size/2:
            ysatx0.append(i[1])
    
    if len(set(xsaty0)) == len(xsaty0):
        print("X's Unique")
    else:
        print(len(xsaty0) - len(set(xsaty0)))
        print(len(xsaty0))
        print("X's Not Unique")
    if len(set(ysatx0)) == len(ysatx0):
        print("Y's Unique")
    else:
        print(len(ysatx0) - len(set(ysatx0)))
        print(len(ysatx0))
        print("Y's Not Unique")
    print(result)
    plt.figure()
    plt.imshow(result)
    return ix,iy,legendre,heading







#                cellMagneticVectors(latitude,longitude,altitude(km),date,  yaw, pitch, roll)
#                                                                           roll,

ix1, iy1, col1,heading = cellMagneticVectors(50.734531, -3.53, .2, 2023+((31+28)/365), 0, 0, 0, 1000000)
ix1, iy1, col2,heading = cellMagneticVectors(50.734531, -3.53, .2, 2023+((31+28)/365), 90, 0, 0, 1000000)
coldif = col2-col1
# plt.figure(figsize = (8,7))
# plt.scatter(ix1,iy1,c = coldif)
# plt.title('Difference between the two graphs')
# plt.colorbar()
#numpoints = 10000000
#xc, yc, zc, cphi, ctheta,lr, lphi, ltheta = coordinateChange.createEye(numpoints)
#print(len(xc)/numpoints)
ix1n = (ix1 - ix1.min())/(ix1.max()-ix1.min())
iy1n = (iy1 - iy1.min())/(iy1.max()-iy1.min())
col1n = col1*255
col2n = col2*255
col1n = col1n.astype('int')
col2n = col2n.astype('int')
ix1n,iy1n = ix1n *2000,iy1n*2000
ix1n = ix1n.astype('int')
iy1n = iy1n.astype('int')
print(ix1n.max())
print(iy1n.max())
coldifn = col2n-col1n
uniquexy =[]
pic = np.zeros([2001,2001])
pic[ix1n,iy1n]=col1n
# for i in range(len (ix1n)):
    
#     print(i)
#     if [ix1n[i],iy1n[i]] in uniquexy:
#         print('nonunique')
#         print(ix1n[i],iy1n[i])
#         break
#     else:
#         uniquexy.append([ix1n[i],iy1n[i]])

# plt.figure(figsize = (8,7))
# plt.scatter(ix1n,iy1n,c = coldifn)
# plt.title('Difference between the two graphs')
# plt.colorbar()
print(ix1.min(),ix1.max())
print(iy1.min(),iy1.max())
print('done')
# img.imageAnalasis(r'D:\!University\University_year_3\Master_Project\Fully_Self_Made_Code\MicrosoftTeams-image.png')

coordinateChange.eulerChangeToCamera(0 , 1, 0, 0, 0, 0)
coordinateChange.eulerChangeToCamera(0 , 1, 0, 90, 0, 0) #Rotates around Z axis by 90deg anticlockwize roll
coordinateChange.eulerChangeToCamera(0 , 1, 0, 0, 90, 0) #Rotates around X axis by 90deg anticlockwize pitch bird faces north
coordinateChange.eulerChangeToCamera(0 , 1, 0, 0, 0, 90) #Rotates around Y axis by 90deg anticlockwize yaw
# ix1, iy1, col1,heading = cellMagneticVectors(50.734531, -3.53, .2, 2023+((31+28)/365), 0, 0, 0)