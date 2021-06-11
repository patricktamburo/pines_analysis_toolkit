from astropy.time import Time
from astropy import coordinates as coord
from astropy import units as u

def jd_utc_to_bjd_tdb(jd_utc, ra, dec, location='lowell'):
    ''' AUTHORS
            Patrick Tamburo, BU, June 2021
        PURPOSE
            Converts Julian Date in UTC timescale to Barycentric Julian Date in TDB timescale. 
        INPUTS
            jd_utc (float): the julian date in utc 
            ra (str): the ra of the target (e.g., '23:06:30.0')
            dec (str): the ded of the target (e.g, '-05:01:57')
            location (str, optional): the site of observations, must match a location in the astropy .of_site json file. 
        OUTPUTS
            Time in BJD TDB.
        TODO:
            None.
    '''

    site = coord.EarthLocation.of_site(location)
    input_jd_utc = Time(jd_utc, format='jd', scale='utc', location=site)
    target = coord.SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
    ltt_bary = input_jd_utc.light_travel_time(target)
    return (input_jd_utc.tdb + ltt_bary).value
