# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 17:12:28 2023

@author: ram86
"""
from magnetic_code import coordinateChange, IGRFMaster
from magnetic_code import error_detect as err
import numpy as np
import importlib


def cellVec (lat,long,alt,date,pitchi,rolli,yaw,points,file = None):
    ## Conversion from pitch being off of horizon to being from z axis
    pitch = pitchi + 90
    pitch = pitch % 180
    ## A roll of 0 leads to y pointing down so it needs to be modified by 180
    roll = rolli + 180
    roll = roll % 360
    ## This function gets the magnetic field based on the lat, long, altitude and date
    eff, hoz, dec, inc, mX, mY, mZ = IGRFMaster.getIGRF(long, lat, alt, date)
    ## Changes coordinates to East->X, North -> Y, Up -> Z
    mX, mY, mZ = coordinateChange.zDownToZUp(mX, mY, mZ)
    ## Goes from World Coordinates to Camera Coordinates
    mX, mY, mZ,heading = coordinateChange.eulerChangeToCamera(mX, mY, mZ, roll, pitch, yaw)
    ## Finds the spherical coordinates of the Magnetic field
    r, thetam, phim = coordinateChange.toSpherical(mX,mY,mZ)
    ## Creates an eye with the given number of cells 
    xc, yc, zc, cphi, ctheta,lr, lphi, ltheta = coordinateChange.createEye(points)
    ## Finds the location on the pixel plane that the eye cell goes to projected from pupil
    iy,ix = coordinateChange.findISpace(xc, yc, zc, 1)
    ## Calculates the angle between the vector to the cell and the magnetic field
    abtwl = coordinateChange.cartAngleBetween(xc,yc,zc,mX,mY,mZ)
    ## Calculates legendre spherical harmonics and normalises it based on given weightings
    if file == None:
        legendre = coordinateChange.legendre(abtwl)
    else:
        i = importlib.import_module(file)
        legendre = file.modulation(abtwl)
    legendre = (legendre - legendre.min())/(legendre.max() - legendre.min())
    coords = []
    for i in range(len(xc)):
        coords.append((ix[i],iy[i],xc[i],yc[i],zc[i],legendre[i]))
    x_coords = sorted(list(set(i[0] for i in points)))
    y_coords = sorted(list(set(i[1] for i in points)))
    result = [[0 for i in range(len(x_coords))] for j in range(len(y_coords))]
    for i in points:
        x_final = x_coords.index(i[0])
        y_final = y_coords.index(i[1])
        result[y_final][x_final] = i[5]
    for i in result:
        print(i)

    return xc, yc,zc, iy, ix, legendre,heading

def mstr():
    print('''
    ********************************************************
    *      Magnetic Field Modulation Pattern Generator     *
    *                                                      *
    * This is a program that will create a modulation      *
    * pattern for a simulated eye based off of a radical   *
    * pair mechanism. This will, by default, use spherical *
    * harmonics (Ylm) with m = 0 and l = 2n. But can be    *
    * changed.                                             *
    *                                                      *
    * To get the magnetic field the IGRF version 13 is     *
    * used. It is valid for 1900-2025 with higher accuracy *
    * and until 2030 with reduced accuracy. Only values    *
    * between 1945 and 2015 are seen as definitive.        *
    *                                                      *
    *------------------------------------------------------*
    *                                                      *
    *                 Variable description:                *
    *                                                      *
    * Latitude: The latitude in the WGS 86 elipsoid in     *
    *     degrees with -ve values being S and +ve N.       *
    * Longitude: The longitude in the WGS 86 elipsoid in   *
    *     degrees with -ve values being W and +ve E.       *
    * Altitude: The altitude in km which the simulated eye *
    *     is located.                                      *
    * Pitch: The angle the simulated eye makes with the    *
    *     horizon between -90 pointing up, +90 pointing    *
    *     down.                                            *
    * Roll: The angle around the centre of the eye the eye *
    *     is rotated with an angle of 0 being the eye      *
    *     parralel to the ground. (0-360)                  *
    * Yaw: The angle off of Geographic North the eye is    *
    *     rotated. (0-360)                                 *
    * Eye Points: The number of points you are simulating  *
    *     within the eye. These are evenly spaced out      *
    * Modulation function: This is the function that takes *
    *     the angle between the magnetic field and the     *
    *     vector from the centre of the eye and the cell   *
    *     location. By default this is the spherical       *
    *     harmonic (Yml) with m = 0 and l being even but   *
    *     this can be changed as long as the function is   *
    *     named 'modulation' within a python file you      *
    *     state and only takes in the angles between the   *
    *     vectors as input. This file must be in the same  *
    *     directory as this one.                           *
    *                                                      *
    *------------------------------------------------------*
    *                                                      *
    *                  Output description:                 *
    * The files will be outputted to a directory called    *
    * data that must be created beforehand. The filename   *
    * is date_lat_lon_alt_pitch_roll_yaw.csv. The csv will *
    * be composed of 6 columns as described below.         *
    *                                                      *
    * The coordinate system used for the eye is a right    *
    * handed coordinate system with z pointing out of the  *
    * pupil, y pointing upwards, and x point to the left.  *
    *
    * Cx: this is the X coordinate of cells in the eye     *
    *     with the centre of the eye being the origin.     *
    * Cy: This is the Y coordinate of cells in the eye.    *
    * Cz: This is the Z coordinate of cells in the eye.    *
    * Ix & Iy: are the xy coordinates of the cells         *
    *     projected onto the xy plane at z = -1 by         *
    *     connecting the pupil to the cell of the eye and  *
    *     following that to the plane mentioned previously.*
    * Modulation: This is the normalised value of the      *
    *     modulation pattern for this cell.                *
    *                                                      *
    ********************************************************
    ''')
    crdvr = err.inpTF('Do you wish to loop betweenn multiple coordinates? (y/n): ')
    altvr = err.inpTF('Do you wish to vary in altitude? (y/n): ')
    ptchvr = err.inpTF('Do you wish to vary in pitch? (y/n): ')
    rllvr = err.inpTF('Do you wish to vary in roll? (y/n): ')
    yawvr = err.inpTF('Do you wish to vary in yaw? (y/n): ')
    ptsnm = err.inpTF('Do you wish to change the number of eye cells from the default (1,000,000)? (y/n): ')
    ownfle = err.inpTF('Do you wish to use your own file for the modulation pattern? (y/n): ')
    date = err.inpDate()
    if crdvr:
        minlat = err.inpFlt('What is the minimum latitude? (-90,90): ',-90,90)
        maxlat = err.inpFlt(f'What is the maximum latitude? ({minlat},90): ',minlat,90)
        numlat = err.inpInt('How many latitudes do you want? (min 2): ',2)
        minlong = err.inpFlt('What is the minimum longitude? (-180,180): ',-180,180)
        maxlong = err.inpFlt(f'What is the maximum longitude? ({minlong},180): ',minlong,180)
        numlong = err.inpInt('How many longitudes do you want? (min 2): ',2)
        lats = np.linspace(minlat,maxlat,numlat)
        longs = np.linspace(minlong,maxlong,numlong)
    else:
        lat = err.inpFlt('What is the latitude of the point? (-90,90): ',-90,90)
        long = err.inpFlt('What is the longitude of the point? (-180,180): ',-180,180)
    
    if altvr:
        minalt = err.inpFlt('What is the minimum altitude you want in km? (min 0): ',0)
        maxalt = err.inpFlt(f'What is the maximum altitude you want in km? (min {minalt}): ',minalt)
        numalt = err.inpInt('How many altitudes do you want? (min 2): ',2)
        alts = np.linspace(minalt,maxalt,numalt)
    else:
        alt = err.inpFlt('What is the altitude of the point in km? (min 0): ',0)
    
    if ptchvr:
        minalt = err.inpFlt('What is the minimum pitch of the camera? (-90,90): ',-90,90)
        maxalt = err.inpFlt(f'What is the maximum pitch of the camera? ({minalt},90): ',minalt,90)
        numalt = err.inpInt('How many pitches do you want? (min 2): ',2)
        pitches = np.linspace(minalt,maxalt,numalt)
    else:
        pitch = err.inpFlt('What is the pitch of the camera? (-90,90): ',-90,90)


    if rllvr:
        minalt = err.inpFlt('What is the minimum roll of the camera? (0,360): ',0,360)
        maxalt = err.inpFlt(f'What is the maximum roll of the camera? ({minalt},360): ',minalt,360)
        numalt = err.inpInt('How many rolls do you want? (min 2): ',2)
        rolls = np.linspace(minalt,maxalt,numalt)
    else:
        roll = err.inpFlt('What is the roll of the camera? (0,360): ',0,360)

    if yawvr:
        minalt = err.inpFlt('What is the minimum yaw of the camera? (0,360): ',0,360)
        maxalt = err.inpFlt(f'What is the maximum yaw of the camera? ({minalt},360):',minalt,360)
        numalt = err.inpInt('How many yaw do you want? (min 2): ',2)
        yaws = np.linspace(minalt,maxalt,numalt)
    else:
        yaw = err.inpFlt('What is the yaw of the camera? (0,360): ',0,360)
    
    if ptsnm:
        points = err.inpInt('How many points in the eye do you wish to generate? ',1)
    else:
        points = 1000000
    
    if ownfle:
        file = input('What is the name of the file in which the function is stored (without the .py extension)? ')
    else:
        file = None
    date = round(date , 4)
    headings = []
    if crdvr and altvr and ptchvr and rllvr and yawvr:
        for lat in lats:
            for long in longs:
                for alt in alts:
                    for pitch in pitches:
                        for roll in rolls:
                            for yaw in yaws:
                                file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                                xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                                headings.append(heading)
                                makefile(file,xc,yc,zc,iy,ix,legendre)
    elif crdvr and altvr and ptchvr and rllvr:
        for lat in lats:
            for long in longs:
                for alt in alts:
                    for pitch in pitches:
                        for roll in rolls:
                            file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                            xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                            headings.append(heading)
                            makefile(file,xc,yc,zc,iy,ix,legendre)
    elif crdvr and altvr and ptchvr and yawvr:
        for lat in lats:
            for long in longs:
                for alt in alts:
                    for pitch in pitches:
                        for yaw in yaws:
                            file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                            xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                            headings.append(heading)
                            makefile(file,xc,yc,zc,iy,ix,legendre)
    elif crdvr and altvr and rllvr and yawvr:
        for lat in lats:
            for long in longs:
                for alt in alts:
                    for roll in rolls:
                        for yaw in yaws:
                            file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                            xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                            headings.append(heading)
                            makefile(file,xc,yc,zc,iy,ix,legendre)
    elif crdvr and ptchvr and rllvr and yawvr:
        for lat in lats:
            for long in longs:
                for pitch in pitches:
                    for roll in rolls:
                        for yaw in yaws:
                            file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                            xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                            headings.append(heading)
                            makefile(file,xc,yc,zc,iy,ix,legendre)
    elif altvr and ptchvr and rllvr and yawvr:
        for alt in alts:
            for pitch in pitches:
                for roll in rolls:
                    for yaw in yaws:
                        file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                        xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                        headings.append(heading)
                        makefile(file,xc,yc,zc,iy,ix,legendre)
    elif crdvr and altvr and ptchvr:
        for lat in lats:
            for long in longs:
                for alt in alts:
                    for pitch in pitches:
                        file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                        xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                        headings.append(heading)
                        makefile(file,xc,yc,zc,iy,ix,legendre)
    elif crdvr and altvr and rllvr:
        for lat in lats:
            for long in longs:
                for alt in alts:
                    for roll in rolls:
                        file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                        xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                        headings.append(heading)
                        makefile(file,xc,yc,zc,iy,ix,legendre)
    elif crdvr and ptchvr and rllvr:
        for lat in lats:
            for long in longs:
                for pitch in pitches:
                    for roll in rolls:
                        file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                        xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                        headings.append(heading)
                        makefile(file,xc,yc,zc,iy,ix,legendre)
    elif altvr and ptchvr and rllvr:
        for alt in alts:
            for pitch in pitches:
                for roll in rolls:
                    file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                    xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                    headings.append(heading)
                    makefile(file,xc,yc,zc,iy,ix,legendre)
    elif crdvr and altvr and yawvr:
        for lat in lats:
            for long in longs:
                for alt in alts:
                    for yaw in yaws:
                        file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                        xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                        headings.append(heading)
                        makefile(file,xc,yc,zc,iy,ix,legendre)
    elif crdvr and ptchvr and yawvr:
        for lat in lats:
            for long in longs:
                for pitch in pitches:
                    for yaw in yaws:
                        file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                        xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                        headings.append(heading)
                        makefile(file,xc,yc,zc,iy,ix,legendre)
    elif altvr and ptchvr and yawvr:
        for alt in alts:
            for pitch in pitches:
                for yaw in yaws:
                    file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                    xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                    headings.append(heading)
                    makefile(file,xc,yc,zc,iy,ix,legendre)
    elif crdvr and  rllvr and yawvr:
        for lat in lats:
            for long in longs:
                for roll in rolls:
                    for yaw in yaws:
                        file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                        xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                        headings.append(heading)
                        makefile(file,xc,yc,zc,iy,ix,legendre)
    elif altvr and rllvr and yawvr:
        for alt in alts:
            for roll in rolls:
                for yaw in yaws:
                    file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                    xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                    headings.append(heading)
                    makefile(file,xc,yc,zc,iy,ix,legendre)
    elif ptchvr and rllvr and yawvr:
        for pitch in pitches:
            for roll in rolls:
                for yaw in yaws:
                    file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                    xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                    headings.append(heading)
                    makefile(file,xc,yc,zc,iy,ix,legendre)
    elif crdvr and altvr:
        for lat in lats:
            for long in longs:
                for alt in alts:
                    file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                    xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                    headings.append(heading)
                    makefile(file,xc,yc,zc,iy,ix,legendre)
    elif crdvr and ptchvr:
        for lat in lats:
            for long in longs:
                for pitch in pitches:
                    file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                    xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                    headings.append(heading)
                    makefile(file,xc,yc,zc,iy,ix,legendre)
    elif altvr and ptchvr:
        for alt in alts:
            for pitch in pitches:
                file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                headings.append(heading)
                makefile(file,xc,yc,zc,iy,ix,legendre)
    elif crdvr and rllvr:
        for lat in lats:
            for long in longs:
                for roll in rolls:
                    file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                    xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                    headings.append(heading)
                    makefile(file,xc,yc,zc,iy,ix,legendre)
    elif altvr and rllvr:
        for alt in alts:
            for roll in rolls:
                file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                headings.append(heading)
                makefile(file,xc,yc,zc,iy,ix,legendre)
    elif ptchvr and rllvr:
        for pitch in pitches:
            for roll in rolls:
                file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                headings.append(heading)
                makefile(file,xc,yc,zc,iy,ix,legendre)
    elif crdvr and yawvr:
        for lat in lats:
            for long in longs:
                for yaw in yaws:
                    file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                    xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                    headings.append(heading)
                    makefile(file,xc,yc,zc,iy,ix,legendre)
    elif altvr and yawvr:
        for alt in alts:
                for yaw in yaws:
                    file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                    xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                    headings.append(heading)
                    makefile(file,xc,yc,zc,iy,ix,legendre)
    elif rllvr and yawvr:
        for roll in rolls:
            for yaw in yaws:
                file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                headings.append(heading)
                makefile(file,xc,yc,zc,iy,ix,legendre)
    elif ptchvr and yawvr:
        for pitch in pitches:
            for yaw in yaws:
                file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                headings.append(heading)
                makefile(file,xc,yc,zc,iy,ix,legendre)
    elif yawvr:
        for yaw in yaws:
            file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
            xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
            headings.append(heading)
            makefile(file,xc,yc,zc,iy,ix,legendre)
    elif rllvr:
        for roll in rolls:
            file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
            xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
            headings.append(heading)
            makefile(file,xc,yc,zc,iy,ix,legendre)
    elif ptchvr:
        for pitch in pitches:
            file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
            xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
            headings.append(heading)
            makefile(file,xc,yc,zc,iy,ix,legendre)
    elif altvr:
        for alt in alts:
            file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
            xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
            headings.append(heading)
            makefile(file,xc,yc,zc,iy,ix,legendre)
    elif crdvr:
        for lat in lats:
            for long in longs:
                file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
                xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
                headings.append(heading)
                makefile(file,xc,yc,zc,iy,ix,legendre)
    else:
        file = f'./data/{date}_{lat}_{long}_{alt}_{pitch}_{roll}_{yaw}.csv'
        xc, yc,zc, iy, ix, legendre, heading = cellVec(lat,long,alt,date,pitch,roll,yaw,points,file)
        headings.append(heading)
        makefile(file,xc,yc,zc,iy,ix,legendre)
    with open(r'./data/!headings.txt','w') as f:
        for i in range(len(headings)):
            f.write(f'{headings[i]}\n')

    
    
def makefile (name,cx,cy,cz,iy,iz,legend):
    #This function saves the data to a file
    
    row = []
    for i in range(len(cx)):
        row.append([cx[i],cy[i],cz[i],iz[i],iy[i],legend[i]])
    row = np.asarray(row)
    np.savetxt(name,row,delimiter = ',')

mstr()
