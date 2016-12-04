import numpy as np 

def date2jd(yr, mon, day, hr, minute, sec):
    """Convert date to Julian Date
        
           Purpose:
               - Converts UTC date to julian date valid 1 Mar 1900 to 28 Feb 2100
        
           JD = date2jd(yr, mon, day, hr, min, sec)
        
           Inputs:
               - yr - 4 digit year (between 1900 and 2100)
               - mon - month (between 01 and 12)
               - day - day (between 1 and 31)
               - hr - UT hour (between 0 and 23)
               - min - UT min (between 0 and 59)
               - sec - UT sec (between 0 and 59.999)
        
           Outputs:
               - JD - julian date (days from 1 Jan 4713 BC 12 Noon)
               - MJD - modified julian date
        
           Dependencies:
               - none
        
           Author:
               - Shankar Kulumani 21 Oct 2012
                   - list revisions
               - Shankar Kulumani 3 Dec 2016
                   - convert to python

           References
               - Vallado
               - USAFA Astro 321
    """

    JD = 367.0 * yr - np.floor( (7 * (yr + np.floor( (mon + 9) / 12.0) ) ) * 0.25 ) + np.floor( 275 * mon / 9.0 ) + day + 1721013.5 + ( (sec/60.0 + minute ) / 60.0 + hr ) / 24.0

    MJD = JD - 2400000.5

    return (JD, MJD)

if __name__ == "__main__":
    # test JD for J2000
    yr = 2000
    mon = 1
    day = 1 
    hour = 12
    minute = 0
    sec = 0 

    JD,MJD = date2jd(yr,mon,day,hour,minute,sec)

    print("J2000: %9.2f" % JD)