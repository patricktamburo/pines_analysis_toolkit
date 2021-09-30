from astropy.time import Time
from astropy import coordinates as coord
from astropy import units as u

def bjd_tdb_to_jd_utc(bjd_tdb, ra, dec, location='lowell'):
    ''' AUTHORS
            Patrick Tamburo, BU, August 2021
        PURPOSE
            Converts a Barycentric Julian Date in TDB timescale to Julian Date in UTC timescale.
        INPUTS
            bjd_tdb (float): the barycentric julian date in tdb
            ra (str): the ra of the target (e.g., '23:06:30.0')
            dec (str): the ded of the target (e.g, '-05:01:57')
            location (str, optional): the site of observations, must match a location in the astropy .of_site json file. 
        OUTPUTS
            Time in JD UTC.
        TODO:
            None.
    '''

    site = coord.EarthLocation.of_site(location)
    input_bjd_tdb = Time(bjd_tdb, format='jd', scale='tdb', location=site)
    target = coord.SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
    ltt_bary = input_bjd_tdb.light_travel_time(target)
    return (input_bjd_tdb.utc - ltt_bary).value