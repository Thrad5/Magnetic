# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 14:02:07 2023

"""
from magnetic_code import igrf_utils as util
from scipy import interpolate
def getIGRF (lon, lat, alt, date, wgs = True):
    '''
    Function that uses the IGRF utils python script by:
        Ciaran Beggan(British Geological Survey)
        With acknowledgements to: Clemens Kloss (DTU Space), David Kerridge (BGS),
        William Brown and Grace Cox.
    This converts longitude latitude altitude and the date into the Magnetic
    Field in Inclanation declination and radius and Cartisien coordinates with 
    X being north
    This is gotten from the 13th edition of the IGRF valid from 1900-2030.
    
    The main-field models for 1900.0, 1905.0,..1940.0 and 2020.0 are 
     non-definitive, those for 1945.0, 1950.0,...2015.0 are definitive and
     the secular-variation model for 2020.0 to 2025.0 is non-definitive.

     Main-field models are to degree and order 10 (i.e. 120 coefficients)
     for 1900.0-1995.0 and to 13 (i.e. 195 coefficients) for 2000.0 onwards. 
     The predictive secular-variation model is to degree and order 8 (i.e. 80
     coefficients).
    

    Parameters
    ----------
    lon : FLOAT
        Longitude of the location in decimal degrees.
    lat : FLOAT
        Latitude of the location in decimal degrees.
    alt : FLOAT
        Altitude of the location in km.
        IF wgs is False alt > 3485.
    date : FLOAT
        The date between 1900.0 and 2030.0 
        with the decimal being the progression through the year.
    wgs : BOOL, optional
        This decides whether the WGS-86 elipsoid is used (what is in google maps eg) for coordinates or spherical earth is used. 
        The default is True (WGS used).

    Returns
    -------
    
    eff : FLOAT
        Total Field Strength (nT).
    hoz : FLOAT
        Horizontal Field Strength (projected into xy plane) (nT).
    dec : FLOAT
        Declination of the magnetic field in degrees (angle from X axis to hoz).
    inc : FLOAT
        Inclination of the magnetic field in degrees (angle from xy plane down +ve).
    X : FLOAT
        North component of the magnetic field (nT).
    Y : FLOAT
        East component of the magnetic field (nT).
    Z : FLOAT
        Downwards component of the magnetic field (nT).

    '''
    
    IGRF_FILE = r'./magnetic_code/IGRF13.shc'
    igrf = util.load_shcfile(IGRF_FILE, None)
    colat = 90-lat
    if wgs == True:
        alt, colat, sd, cd = util.gg_to_geo(alt, colat) #Converts the geoditic altitude and colatitude to geocentric
    else:
        sd = 0
        cd = 0
    # Interpolate the geomagnetic coefficients to the desired date(s)
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    f = interpolate.interp1d(igrf.time, igrf.coeffs, fill_value='extrapolate')
    coeffs = f(date)    
    
    # Compute the main field B_r, B_theta and B_phi value for the location(s) 
    Br, Bt, Bp = util.synth_values(coeffs.T, alt, colat, lon, igrf.parameters['nmax'])
    # For the SV, find the 5 year period in which the date lies and compute
    # the SV within that period. IGRF has constant SV between each 5 year period
    # We don't need to subtract 1900 but it makes it clearer:
    epoch = (date-1900)//5    
    epoch_start = epoch*5
    # Add 1900 back on plus 1 year to account for SV in nT per year (nT/yr):
    coeffs_sv = f(1900+epoch_start+1) - f(1900+epoch_start)   
    Brs, Bts, Bps = util.synth_values(coeffs_sv.T, alt, colat, lon, igrf.parameters['nmax'])
    
    # Use the main field coefficients from the start of each five epoch
    # to compute the SV for Dec, Inc, Hor and Total Field (F) 
    # [Note: these are non-linear components of X, Y and Z so treat separately]
    coeffsm = f(1900+epoch_start);
    Brm, Btm, Bpm = util.synth_values(coeffsm.T, alt, colat, lon, igrf.parameters['nmax'])
    
    
    # Rearrange to X, Y, Z components 
    X = -Bt; Y = Bp; Z = -Br
    # For the SV
    dX = -Bts; dY = Bps; dZ = -Brs 
    Xm = -Btm; Ym = Bpm; Zm = -Brm
    # Rotate back to geodetic coords if needed
    if (wgs == True):
        t = X; X = X*cd + Z*sd;  Z = Z*cd - t*sd
        t = dX; dX = dX*cd + dZ*sd;  dZ = dZ*cd - t*sd
        t = Xm; Xm = Xm*cd + Zm*sd;  Zm = Zm*cd - t*sd
        
    # Compute the four non-linear components 
    dec, hoz, inc, eff = util.xyz2dhif(X,Y,Z)
    # The IGRF SV coefficients are relative to the main field components 
    # at the start of each five year epoch e.g. 2010, 2015, 2020
    decs, hozs, incs, effs = util.xyz2dhif_sv(Xm, Ym, Zm, dX, dY, dZ)
    
    return  eff, hoz, dec, inc, X, Y, Z
    
    