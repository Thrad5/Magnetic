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
    isdeg : Bool, optional
        Changes the angles' measurement units between degrees(True) and radians(False). The default value is True
        
    Returns
    -------
    X : Float
        Left component from the camera.
    Y : Float
        Up component from the camera.
    Z : Float
        Out of lens component of the camera.
    '''
    if isdeg == True:
        roll = deg2rad(roll)
        pitch = deg2rad(pitch)
        yaw = deg2rad(yaw)
    cll00 = np.cos(roll)*np.cos(yaw) - (np.sin(pitch)*np.sin(roll)*np.sin(yaw))
    cll01 = - (np.cos(pitch)*np.sin(roll))
    cll02 = np.cos(pitch)*np.sin(yaw) + np.cos(yaw) * np.sin(pitch)*np.sin(roll)
    cll10 = np.cos(yaw)*np.sin(roll) + np.cos(roll)*np.sin(pitch)*np.sin(yaw)
    cll11 = np.cos(roll)*np.cos(pitch)
    cll12 = (np.sin(roll)*np.sin(yaw)) - (np.cos(roll) * np.cos(yaw) * np.sin(pitch))
    cll20 = - (np.cos(pitch) * np.sin(yaw))
    cll21 = np.sin(pitch)
    cll22 = np.cos(pitch) * np.cos(yaw)
    
    mat = np.array([[cll00,cll01,cll02],[cll10,cll11,cll12],[cll20,cll21,cll22]])
    #print(mat)
    vec = np.array([X,Y,Z])
    #print('Original Vector:',vec)
    changed = mat @ vec
    #print('Changed Vector:',changed)
    inv = np.linalg.inv(mat)
    birdfacing = inv @ [0,0,1]
    #print(birdfacing)
    phi = np.arctan2(birdfacing[0],birdfacing[1])
    heading = 360-((phi+360)%360)
    return changed[0], changed[1], changed[2],heading

def rebound(theta,phi):
    '''
    Rebounds the polar oordinates so they are in the ranges of 0-180 and 0-360

    Parameters
    ----------
    theta : float
        Angle from the positive z axis.
    phi : float
        Angle around the positive z axis.

    Returns
    -------
    theta : float
        Angle from the positive z axis.
    phi : float
        Angle around the positive z axis.

    '''
    if 0<theta or theta>180:
        theta = abs(theta)
        phi = phi +180
        
    phi = phi%360
    return theta,phi

def toSpherical (X,Y,Z):
    r = np.sqrt(X**2+Y**2+Z**2)
    theta = np.arccos(Z/r)
    phi = np.sign(Y)*np.arccos(X/(np.sqrt(X**2+Y**2)))
    theta = theta*180/np.pi
    phi = phi *180/np.pi
    return r, theta, phi

def heading(dec, theta, heading, aoh):
    '''
    Changes the declination / phi term in spherical polars to the phi term of 

    Parameters
    ----------
    dec : FLOAT or np.array
        The Declination off of geographic north that the magnetic field is.
    theta : Float or np.array
        The angle off of up perpendicular to earth of the magnetic field.
    heading : FLOAT
        The heading of the bird relative to geographic north.
    aoh : FLOAT
        The angle off of the horizon of the bird with negative being down

    Returns
    -------
    dec : TYPE
        DESCRIPTION.

    '''
    dec = heading + dec
    theta = theta +aoh
    theta,dec = rebound(theta,dec)
    return dec, theta

def createEye(num_pnts):
    
    #Assume that the eye is covered in the back half by retina
    #This will use the fibonacci spiral method
    n=round(num_pnts*(1/0.8213946))
    i = np.arange(0,n,1)
    gld = (1 + 5**0.5)/2
    x = (i/gld)%1
    y = i/n
    rad=1#5mm radius
    phi = 2*np.pi*y
    theta = np.arccos(1 - 2*x)
    
    # print('Eye debug:')
    # print('Theta min, theta max')
    # print(theta.min(),theta.max())
    # print('Phi min, Phi max')
    # print(phi.min(),phi.max())
    
    
    
    xc1,yc1,zc1 = rad*np.cos(phi)*np.sin(theta),rad*np.sin(phi)*np.sin(theta),rad*np.cos(theta[i])
    selec = cartAngleBetween(xc1,yc1,zc1,1,0,0)


    xc,yc,zc = [],[],[]
    for i in range(0,len(theta)):
        
        # print(selec[i])
        if (selec[i]>(50*np.pi/180)):
            xc.append(rad*np.cos(phi[i])*np.sin(theta[i]))
            yc.append(rad*np.sin(phi[i])*np.sin(theta[i]))
            zc.append(rad*np.cos(theta[i]))
    tem = xc
    temp = yc
    temp2 = zc
    xc = np.array(temp)
    yc = np.array(temp2)
    zc = np.array(tem)
    lz = zc - rad
    lr,ltheta,lphi = toSpherical(xc,yc,lz)
    cr,ctheta,cphi = toSpherical(xc,yc,zc)
    # print('xc min, xc max')
    # print(xc.min(),xc.max())
    # print('yc min, yc max')
    # print(yc.min(),yc.max())
    # print('zc min, zc max')
    # print(zc.min(),zc.max())
    # print('Ltheta min, Ltheta max')
    # print(ltheta.min(),ltheta.max())
    # print('Lphi min, Lphi max')
    # print(lphi.min(),lphi.max())
    # print('Lr min, Lr max')
    # print(lr.min(),lr.max())
    lphi = lphi * 180/np.pi
    ltheta = ltheta * 180/np.pi
    cphi = cphi * 180/np.pi
    ctheta = ctheta * 180/np.pi
    return xc, yc, zc, cphi, ctheta,lr, lphi, ltheta

def recalcXYZ(r, theta, phi, isdeg=True):
    if isdeg == True:
        theta = deg2rad(theta)
        phi = deg2rad(phi)
    return r*np.cos(phi)*np.sin(theta),r*np.sin(phi)*np.sin(theta),r*np.cos(theta)

def deg2rad(angle):
    return angle*np.pi/180

def cartAngleBetween(x1,y1,z1,x2,y2,z2):
    xx = x1*x2
    yy = y1*y2
    zz = z1*z2
    mag1 = np.sqrt(x1**2 + y1**2 + z1**2)
    mag2 = np.sqrt(x2**2 + y2**2 + z2**2)
    dot = xx + yy + zz
    angle = np.arccos(dot/(mag1*mag2))
    return angle

def angleBetween(theta1, phi1, theta2, phi2, isdeg = True):
    '''
    Calculates the angle between two vectors
    One of the two vectors can be an array of vectors as long as they are all in the same units.
    
    Parameters
    ----------
    theta1 : Float
        Angle between z-axis and first vector, degrees default.
    phi1 : Float
        Angle around the z-axis of the first vector, degrees default.
    theta2 : Float
        Angle between z-axis and second vector, degrees default.
    phi2 : Float
        Angle around the z-axis of the second vector, degrees default.
    isdeg : Bool, optional
        Changes between degrees and radians for both input and output, degrees default.
        
    Returns
    -------
    angle : Float
        The angle between the two vectors as either degrees or radians.

    '''
    if isdeg == True:
        theta1 = theta1 *np.pi/180
        theta2 = theta2 *np.pi/180
        phi1 = phi1 * np.pi/180
        phi2 = phi2 * np.pi/180
    
    diff_phi = phi1 - phi2
    lh = np.sin(theta1)*np.sin(theta2)*np.cos(diff_phi)
    rh = np.cos(theta1)
    inarccos = lh + rh
    
    # print('In arccos min and max')
    # print(inarccos.min(),inarccos.max())
    try:
        len(inarccos)
        for i in range(len(inarccos)):
            if inarccos[i]< -1:
                inarccos[i] = -1
            elif inarccos[i] > 1:
                inarccos[i] = 1
    except TypeError:
        if inarccos < -1:
            inarccos = -1
        elif inarccos >1:
            inarccos = 1

    
    # print(inarccos.min(),inarccos.max())
    angle = np.arccos(inarccos)
    if isdeg == True:
        angle = angle *180/np.pi
    return angle

def legendre(angle,w2 = 1, w4 = 0, w6 = 0, w8 = 0, isdeg = True):
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
        angle = angle*np.pi/180
        
    with open(r'./magnetic_code/weightings.txt','r') as f:
        weightings = []
        for line in f.readlines():
            weightings.append(float(line))
        weightings = np.asarray(weightings)
    sum = 0
    for i in range(len(weightings)):
        sum += weightings[i]*sci.special.sph_harm(0,i*2,0,angle)
    
    return sum

def proj (x,theta,phi,isdeg = True):
    if isdeg == True:
        theta = deg2rad(theta)
        phi = deg2rad(phi)
    
    r = 1/(x/(np.cos(phi)*np.sin(theta)))
    y = x *np.sin(phi) / np.cos(phi)
    z = r * np.cos(theta)
    return y, z

def findISpace (ex,ey,ez,r):
    pupil = [0,0,r]
    positions = []
    for i in range(len(ex)):
        positions.append([ex[i],ey[i],ez[i]])
    positions = np.array(positions)
    direction_vector = positions - pupil
    t = (-r - pupil[2])/direction_vector[:,2]
    iy = pupil[1]+ t*direction_vector[:,1]
    ix = pupil[0] + t * direction_vector[:,0]
    # print(iy,ix)
    return iy, ix