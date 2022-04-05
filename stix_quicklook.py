import hissw

import os
import spiceypy as spice

from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits

from datetime import datetime, timedelta

import sunpy
from sunpy import map
from sunpy.map import make_fitswcs_header
from sunpy.coordinates.frames import HeliocentricEarthEcliptic,HeliographicStonyhurst
from sunpy.map.maputils import solar_angular_radius

import numpy as np
import spiceypy
from spiceypy.utils.exceptions import NotFoundError
import warnings

##### Constants
km2AU = 6.6845871226706e-9*u.AU/u.km



def load_SOLO_SPICE(path_kernel):
    """
    Load the SPICE kernel that will be used to get the
    coordinates of the different spacecrafts.
    """

    #get cwd
    cwd=os.getcwd()
    
    # Check if path_kernel has folder format
    if path_kernel[-1] != '/':
        path_kernel = path_kernel+'/'

    # Change the CWD to the given path. Necessary to load correctly all kernels
    os.chdir(path_kernel)

    # Load one (or more) SPICE kernel into the program
    spice_kernel = 'solo_ANC_soc-flown-mk.tm'
    spice.furnsh(spice_kernel)

    print()
    print('SPICE kernels loaded correctly')
    print()

    #change back to original working directory
    os.chdir(cwd)
    
def coordinates_body(date_body,body_name,light_time=False):
    """
    Load the kernel needed in order to derive the
    coordinates of the given celestial body and then return them in
    Heliocentric Earth Ecliptic (HEE) coordinates.
    """
    
    # Observing time
    obstime = spice.datetime2et(date_body)
    
    
    # Obtain the coordinates of Solar Orbiter
    hee_spice, lighttimes = spice.spkpos(body_name, obstime,
                                     'SOLO_HEE_NASA', #'HEE', # Should be SOLO_HEE_NASA but doesn't work with heliopy. This is okay in this example because no light times are needed.
                                     'NONE', 'SUN')
    hee_spice = hee_spice * u.km

    # Convert the coordinates to HEE
    body_hee = HeliocentricEarthEcliptic(hee_spice,
                                          obstime=Time(date_body).isot,
                                          representation_type='cartesian')
    if not light_time:
        # Return the HEE coordinates of the body
        return body_hee
    else:
        return body_hee,lighttimes

def get_observer(date_in,rsun,obs='Earth',wcs=True,sc=False,out_shape=(4096,4096)):
    '''Get observer information. Get WCS object if requested. Return Observer as SkyCoord at (0",0") if requested.'''
    hee = coordinates_body(date_in, obs)
    observer=hee.transform_to(HeliographicStonyhurst(obstime=date_in))
    if sc:
        observer = SkyCoord(0*u.arcsec, 0*u.arcsec,obstime=date_in,
                            rsun=rsun,
                            observer=observer,frame='helioprojective')

    if wcs == False:
        return observer
    else:
        if isinstance(observer, SkyCoord):
            refcoord=observer
        else:
            refcoord=SkyCoord(0*u.arcsec, 0*u.arcsec,obstime=date_in,
                              rsun=rsun,
                              observer=observer,frame='helioprojective')
                              
        out_header = sunpy.map.make_fitswcs_header(out_shape,refcoord,observatory=obs)

        out_wcs = WCS(out_header)
        return observer,out_wcs
    
def get_rsun_apparent(date_in,radius,observer=False, spacecraft='Solar Orbiter',sc=True):
    if observer == False:
        obs = get_observer(date_in,radius,obs=spacecraft,wcs=False,sc=sc)
    else:
        obs=observer
    return solar_angular_radius(obs)

def parse_sunspice_name_py(spacecraft):
    '''python version of sswidl parse_sunspice_name.pro '''
    if type(spacecraft) == str:
        sc = spacecraft.upper()
        n = len(sc)

        #If the string is recognized as one of the STEREO spacecraft, then return the appropriate ID value.
        if 'AHEAD' in sc and 'STEREO' in sc or sc == 'STA':
            return '-234'
        
        if 'BEHIND' in sc and 'STEREO' in sc or sc == 'STB':
            return '-235'

        #If SOHO, then return -21.
        if sc=='SOHO':
            return '-21'

        #If Solar Orbiter then return -144.
        if 'SOLAR' in sc and 'ORBITER' in sc or sc == 'SOLO' or sc == 'ORBITER':
            return '-144'

        #If Solar Probe Plus then return -96.
        if 'PARKER' in sc and 'SOLAR' in sc and 'PROBE' in sc:
            return '-96'
        if 'PROBE' in sc and 'PLUS' in sc:
            return '-96'
        if sc == 'PSP' or sc == 'SPP':
            return '-96'

        #If BepiColombo MPO then return -121.
        if 'BEPICOLOMBO' in sc and 'MPO' in sc:
            return '-121'
        if 'BC' in sc and 'MPO' in sc:
            return '-121'
        #Otherwise, simply return the (trimmed and uppercase) original name.
        return sc
    else:
        raise TypeError("Input spacecraft name must be string")

def get_sunspice_roll_py(datestr, spacecraft, path_kernel, system='SOLO_SUN_RTN',degrees=True, radians=False,tolerance=100):
    '''Python version of (simpler) sswidl get_sunspice_roll.pro . Assumes spice kernel already furnished'''
    units = float(180)/np.pi
    if radians: units = 1
  
    #Determine which spacecraft was requested, and translate it into the proper input for SPICE.

    inst = 0
    sc_ahead  = '-234'
    sc_behind = '-235'
    sc_psp    = '-96'
    sc = parse_sunspice_name_py(spacecraft)
    if sc == sc_ahead or sc == sc_behind:
        sc_stereo=True
    if sc == '-144':
        sc_frame = -144000

    #Start by deriving the C-matrices.  Make sure that DATE is treated as a vector.

    if system != 'RTN' and type(system) == str:
        system=system.upper()
      
    #don't know why it wants to live in the kernel directory in order to check errons but it does
    pwd=os.getcwd()
    os.chdir(path_kernel)
    et=spiceypy.spiceypy.str2et(datestr)
    sclkdp=spiceypy.spiceypy.sce2c(int(sc), et) #Ephemeris time, seconds past J2000.
    #print(sc,sclkdp, tolerance, system)
    try:
        (cmat,clkout)=spiceypy.spiceypy.ckgp(sc_frame, sclkdp, tolerance, system)
    except NotFoundError:
        warnings.warn("Spice returns not found for function: ckgp, returning roll angle of 0: " + datestr)
        os.chdir(pwd)
        return 0.0

    os.chdir(pwd)
    
    twopi  = 2.*np.pi
    halfpi = np.pi/2.
  
    #sci_frame = (system.upper()== 'SCI') and sc_stereo

    #cspice_m2eul, cmat[*,*,i], 1, 2, 3, rroll, ppitch, yyaw
    rroll,ppitch,yyaw=spiceypy.spiceypy.m2eul(cmat, 1, 2, 3)
    ppitch = -ppitch
    if (sc == sc_ahead) or (sc == sc_behind): rroll = rroll - halfpi
    if sc == sc_behind: rroll = rroll + np.pi
    if abs(rroll) > np.pi: rroll = rroll - sign(twopi, rroll)

    #Correct any cases where the pitch is greater than +/- 90 degrees
    if abs(ppitch) > halfpi:
        ppitch = sign(np.pi,ppitch) - ppitch
        yyaw = yyaw - sign(np.pi, yyaw)
        rroll = rroll - sign(np.pi, rroll)
      
    #Apply the units.
    roll  = units * rroll
    pitch = units * ppitch
    yaw  = units * yyaw

    return roll

def stix_flare_image_reconstruction(
    background_L1A_fits_filename: str,
    flaring_data_L1A_fits_filename: str,
    flare_start_UTC: str,
    flare_end_UTC: str,
    energy_range_science_channel_upper_limit: int,
    energy_range_science_channel_lower_limit: int,
    spice_kernel_path: str,
    map_filename: str,
    ssw_home: str,
    idl_home: str,
    stx_quicklook_fwdfit_path: str):
    """
    Arguments:
        background_L1A_fits_filename: str,
            path of the L1A background file
        flaring_data_L1A_fits_filename: str,
            path of the L1A science file
        flare_start_UTC: str,
            start of the time interval to consider for image reconstruction
        flare_end_UTC: str,
            end of the time interval to consider for image reconstruction
        energy_range_science_channel_upper_limit: int,
            lower bound of the energy range to consider for image reconstruction
        energy_range_science_channel_upper_limit: int,
            lower bound of the energy range to consider for image reconstruction
        spice_kernel_path: str,
            path of the Solar Orbiter SPICE kerner to use for obtaining the ephemeris data
        map_filename: str,
            path and name of the FWDFIT map reconstructed from STIX data
        ssw_home: str,
            path of the SSW folder
        idl_home: str,
            path of the IDL installation folder
        stx_quicklook_fwdfit_path: str,
            path of the folder containing 'stx_quicklook_fwdfit.pro'
    Returns:
        Fits file of the FWDFIT_PSO map of the flaring event
    """

    # Observer name
    observer = 'Solar Orbiter'
    # Load SPICE kernel
    load_SOLO_SPICE(spice_kernel_path)
    
    this_time = Time(flare_start_UTC).to_datetime()
    
    # Then, get the coordinated of SO at a given date
    solo_hee = coordinates_body(this_time, observer)

    # Transform the coordinates in Heliographic Stonyhurst
    solo_hgs = solo_hee.transform_to(HeliographicStonyhurst(obstime=this_time))

    B0 = solo_hgs.lat.deg # Heliographic latitude (B0 angle)
    L0 = solo_hgs.lon.deg # Heliographic longitude (L0 angle)
    
    solar_radius = 695700*u.km

    # RSUN
    apparent_radius_sun = get_rsun_apparent(this_time,solar_radius,observer=False, 
                                            spacecraft='Solar Orbiter',sc=True)
    
    #ROLL ANGLE Solar Orbiter
    roll_angle_solo = get_sunspice_roll_py(this_time.strftime("%Y-%m-%d %H:%M:%S"), 'Solar Orbiter',
                                           spice_kernel_path, system='SOLO_SUN_RTN',degrees=True, 
                                           radians=False,tolerance=100)
    
    # Define dictionary containing the inputs for the IDL code
    inputs = {'background_L1A_fits_filename': background_L1A_fits_filename,
              'flaring_data_L1A_fits_filename': flaring_data_L1A_fits_filename,
              'flare_start_UTC': flare_start_UTC,
              'flare_end_UTC': flare_end_UTC,
              'energy_range_science_channel_upper_limit': energy_range_science_channel_upper_limit,
              'energy_range_science_channel_lower_limit': energy_range_science_channel_lower_limit,
              'map_filename': map_filename,
              'L0': L0,
              'B0': B0,
              'apparent_radius_sun': apparent_radius_sun,
              'roll_angle_solo': roll_angle_solo}
    
    # Define SSWIDL environment
    ssw = hissw.Environment(ssw_home=ssw_home,idl_home=idl_home,ssw_packages=['stix'])
    # Run SSWIDL and the procedure 'stx_quicklook_fwdfit.pro'
    ssw_resp = ssw.run(stx_quicklook_fwdfit_path + 'stx_quicklook_fwdfit.pro', args=inputs)
    

def create_STIX_map(path_fits, datetime_map, flare_center):
    """
    Converts the input FITS to a standard Sunpy Solar map, by
    returning the map containing the STIX source.
    
    At the moment, the function is optimized for the STIX FITS
    files generated by Paolo and Emma.
    """
    
    # Open the FITS file
    this_fits = fits.open(path_fits)
    
    # In order to properly set the WCS in the header, we need to find 
    # the HEE coordinates of Solar Orbiter and set the proper frame
    
    # Convert the string date and time to datetime.datetime
    datetime_map = Time(datetime_map).to_datetime()
    # Obtain the HEE coordinates of Solar Orbiter
    solo_hee = coordinates_body(datetime_map,'Solar Orbiter')
    
    # Set the coordinates of the reference pixel of the STIX map
    stix_ref_coord = SkyCoord(flare_center[0]*u.arcsec, 
                              flare_center[1]*u.arcsec,
                              obstime=datetime_map,
                              observer=solo_hee.transform_to(HeliographicStonyhurst(obstime=datetime_map)),
                              frame='helioprojective')
    
    # Get the distance of Solar Orbiter from the Sun (center)
    dist_AU = np.sqrt(solo_hee.x**2+solo_hee.y**2+solo_hee.z**2)*km2AU

    # Create a FITS header containing the World Coordinate System (WCS) information
    out_header = make_fitswcs_header(this_fits[0].data,
                                     stix_ref_coord,
                                     scale=(1., 1.)*dist_AU/u.AU*u.arcsec/u.pixel,
                                     instrument="STIX VIS_FWDFIT_PSO map",
                                     observatory="Solar Orbiter")

    # Create the STIX map
    stix_map = map.Map((this_fits[0].data, out_header))
    
    return stix_map
    
    
    
    
    
'''
# Below, a really DIRTY code to put the STIX map on a full-disk image.
# We just create a full-disk image full of zeros and then stuck in there
# the STIX map at the right location. This way is a kind of Art Attack (but it works)


    ###### Create an "empty" Sunpy solar map
    
    # Set the coordinates of the reference pixel of the full-disk map
    stix_fd_ref_coord = SkyCoord(0*u.arcsec, 
                                 0*u.arcsec,
                                 obstime=datetime_map,
                                 observer=solo_hee.transform_to(HeliographicStonyhurst(obstime=datetime_map)),
                                 frame='helioprojective')

    # Create a FITS header containing the World Coordinate System (WCS) information
    # for the full-disk map
    out_header = make_fitswcs_header(out_shape,
                                     stix_fd_ref_coord,
                                     scale=(1., 1.)*dist_AU/u.AU*u.arcsec/u.pixel,
                                     instrument="STIX (amplitudes only)",
                                     observatory="Solar Orbiter")

    # Create the STIX full disk map
    stix_fulldisk = Map((np.zeros(out_shape), out_header))
    
    
    ###### Finally, put the STIX map on the full-disk map
    
    # Create the variable containing the axis on the Sun
    axis_solar_x = np.linspace(float(-dist_AU/u.AU*(out_shape[0]/2.)), 
                               float(dist_AU/u.AU*(out_shape[0]/2.)), 
                               num=out_shape[0])
    axis_solar_y = np.linspace(float(-dist_AU/u.AU*(out_shape[1]/2.)), 
                               float(dist_AU/u.AU*(out_shape[1]/2.)), 
                               num=out_shape[1])
    
    # Find the indices of the center of the axis
    ind_center_x = np.argmin(abs(axis_solar_x-flare_center[0]))
    ind_center_y = np.argmin(abs(axis_solar_y-flare_center[1]))
    
    # Put the STIX map on the full-disk map, where its center matches
    # the center of the axis previously defined
    shape_submap = stix_map.data.shape
    x_min = int(ind_center_x-shape_submap[0]/2)
    x_max = int(ind_center_x+shape_submap[0]/2)
    y_min = int(ind_center_y-shape_submap[1]/2)
    y_max = int(ind_center_y+shape_submap[1]/2)
    stix_fulldisk.data[y_min:y_max,x_min:x_max] = stix_map.data

    return stix_fulldisk
    '''