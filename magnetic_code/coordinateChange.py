# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 11:18:28 2023

@author: ram86
"""
import numpy as np
import scipy as sci
def zDownToZUp (X,Y,Z):
    '''
    Converts from X:North, Y:East, Z:Down to X:East, Y:North, Z:Up

    Parameters
    ----------
    X : Float
        North Component.
    Y : Float
        East Component.
    Z : Float
        Down Component.

    Returns
    -------
    X : Float
        East Component.
    Y : Float
        North Component.
    Z : Float
        Up Component.

    '''
    temp = X
    X = Y
    Y  = temp
    Z = -Z
    
    return X, Y, Z

def eulerChangeToCamera(X,Y,Z,roll,pitch,yaw, isdeg = True):
    '''

    Parameters
    ----------
    X : Float
        East Component.
    Y : Float
        North Component.
    Z : Float
        Up Component.
    roll : Float
        The angle off of the ground of the bird
    pitch : Float
        The angle off of the horizon of the bird
    yaw : Float
        The heading of the bird relative to the world's north
    isdeg : BOOL, optional
        Whether the value output are in degrees(True) or radians(False). The default is True.
        
    Returns
    -------
    X : Float
        Left component from the camera.
    Y : Float
        Up component from the camera.
    Z : Float
        Out of lens component of the camera.
    '''
    #Converts from degrees to radians
    if isdeg == True:
        roll = deg2rad(roll)
        pitch = deg2rad(pitch)
        yaw = deg2rad(yaw)

    #Creates the Euler Angle matrix ZYX
    Xr = [[1,0,0],[0,np.cos(pitch),-np.sin(pitch)],[0,np.sin(pitch),np.cos(pitch)]]
    Yr = [[np.cos(yaw),0,np.sin(yaw)],[0,1,0],[-np.sin(yaw),0,np.cos(yaw)]]
    Zr = [[np.cos(roll),-np.sin(roll),0],[np.sin(roll),np.cos(roll),0],[0,0,1]]
    mat = np.matmul(np.matmul(Zr,Yr),Xr)
    vec = np.array([X,Y,Z])
    #Converts the input magnetic vector into the coordinates of the bird's eye
    changed = mat @ vec
    #Calculates the heading of the bird by calculating the inverted matrix multiplied by the vector going out of the pupil [0,0,1]
    inv = np.linalg.inv(mat)
    birdfacing = inv @ [0,0,1]
    phi = np.arctan2(birdfacing[0],birdfacing[1])
    #Converts the heading int degrees if they are being used
    if isdeg == True:
        phi = phi *180/np.pi
        heading = phi%360
    else:
        heading = phi%(2*np.pi)
    return changed[0], changed[1], changed[2],heading


def toSpherical (X,Y,Z,isdeg = True):
    '''
    

    Parameters
    ----------
    X : REAL
        DESCRIPTION.
    Y : REAL
        DESCRIPTION.
    Z : REAL
        DESCRIPTION.
    isdeg : BOOL, optional
        Whether the value output are in degrees(True) or radians(False). The default is True.

    Returns
    -------
    r : TYPE
        DESCRIPTION.
    theta : TYPE
        DESCRIPTION.
    phi : TYPE
        DESCRIPTION.

    '''
    r = np.sqrt(X**2+Y**2+Z**2)
    theta = np.arccos(Z/r)
    phi = np.sign(Y)*np.arccos(X/(np.sqrt(X**2+Y**2)))
    if isdeg == True:
        theta = deg2rad(theta)
        phi = deg2rad(phi)
    return r, theta, phi



def createEye(num_pnts, isdeg = True):
    '''
    

    Parameters
    ----------
    num_pnts : TYPE
        DESCRIPTION.
    isdeg : BOOL, optional
        Whether the value output are in degrees(True) or radians(False). The default is True.

    Returns
    -------
    xc : TYPE
        DESCRIPTION.
    yc : TYPE
        DESCRIPTION.
    zc : TYPE
        DESCRIPTION.

    '''
    #Assume that the eye is covered in retina exept for areas with an angle less than 50º from the vector to the pupil
    #This will use the fibonacci spiral method for generating the points on the retina
    n=round(num_pnts*(1/0.8213946))
    i = np.arange(0,n,1)
    gld = (1 + 5**0.5)/2
    x = (i/gld)%1
    y = i/n
    rad=1#5mm radius normalised for units of bird's retina
    phi = 2*np.pi*y
    theta = np.arccos(1 - 2*x)
    
    
    
    #This calculates the cartesian coordinates of the eye.
    xc1,yc1,zc1 = rad*np.cos(phi)*np.sin(theta),rad*np.sin(phi)*np.sin(theta),rad*np.cos(theta[i])
    selec = cartAngleBetween(xc1,yc1,zc1,1,0,0)


    xc,yc,zc = [],[],[]
    for i in range(0,len(theta)):
        #Using the angle between the point created and the pupil we can see if it is part of the retina 
        if (selec[i]>(50*np.pi/180)):
            xc.append(rad*np.cos(phi[i])*np.sin(theta[i]))
            yc.append(rad*np.sin(phi[i])*np.sin(theta[i]))
            zc.append(rad*np.cos(theta[i]))
    #The function that creates the eye has the pupil as in the positive x, verticaly upwards as positive z and, a right hand system of coordinates.
    #This needs to be changed to be the pupil as positive z, vertically upwards as positive y, and a right hand system of coordinates.
    tem = xc
    temp = yc
    temp2 = zc
    xc = np.array(temp)
    yc = np.array(temp2)
    zc = np.array(tem)
    lz = zc - rad
   
    return xc, yc, zc



def deg2rad(angle):
    '''
    

    Parameters
    ----------
    angle : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    return angle*np.pi/180

def cartAngleBetween(x1,y1,z1,x2,y2,z2,isdeg = False):
    '''
    

    Parameters
    ----------
    x1 : TYPE
        DESCRIPTION.
    y1 : TYPE
        DESCRIPTION.
    z1 : TYPE
        DESCRIPTION.
    x2 : TYPE
        DESCRIPTION.
    y2 : TYPE
        DESCRIPTION.
    z2 : TYPE
        DESCRIPTION.
    isdeg : BOOL, optional
        Whether the value output are in degrees(True) or radians(False). The default is False.

    Returns
    -------
    angle : TYPE
        DESCRIPTION.

    '''
    xx = x1*x2
    yy = y1*y2
    zz = z1*z2
    mag1 = np.sqrt(x1**2 + y1**2 + z1**2)
    mag2 = np.sqrt(x2**2 + y2**2 + z2**2)
    dot = xx + yy + zz
    angle = np.arccos(dot/(mag1*mag2))
    if isdeg == True:
        angle = deg2rad(angle)
    return angle



def legendre(angle, isdeg = True):
    '''
    Calculates the legendre polynomials for an angle

    Parameters
    ----------
    angle : FLOAT or np.array
        DESCRIPTION.
    isdeg : BOOL, optional
        Whether the angle is in degrees or radians. The default is True.

    Returns
    -------
    leg : FLOAT or np.array (same as angle)
        Weighted sum of spherical harmonics from l=2 to l=4.

    '''
    if isdeg == True:
        angle = deg2rad(angle)
        
    with open(r'./magnetic_code/weightings.txt','r') as f:
        weightings = []
        for line in f.readlines():
            weightings.append(float(line))
        weightings = np.asarray(weightings)
    sum = 0
    for i in range(len(weightings)):
        sum += weightings[i]*sci.special.sph_harm(0,i*2,0,angle)
    
    return sum


def findISpace (ex,ey,ez,r):
    '''
    

    Parameters
    ----------
    ex : TYPE
        DESCRIPTION.
    ey : TYPE
        DESCRIPTION.
    ez : TYPE
        DESCRIPTION.
    r : TYPE
        DESCRIPTION.

    Returns
    -------
    iy : TYPE
        DESCRIPTION.
    ix : TYPE
        DESCRIPTION.

    '''
    
    pupil = [0,0,r]
    positions = []
    for i in range(len(ex)):
        positions.append([ex[i],ey[i],ez[i]])
    positions = np.array(positions)
    direction_vector = positions - pupil
    t = (-r - pupil[2])/direction_vector[:,2]
    iy = pupil[1]+ t*direction_vector[:,1]
    ix = pupil[0] + t * direction_vector[:,0]
    return iy, ix