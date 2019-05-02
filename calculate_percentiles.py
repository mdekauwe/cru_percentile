#!/usr/bin/env python

"""
For each flux site calculate the 95th percentile TMAX using data from CRU

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (05.05.2019)"
__email__ = "mdekauwe@gmail.com"

import numpy as np
import sys
import matplotlib.pyplot as plt
import xarray as xr
import os
import glob
import pandas as pd

def calculate_tmax_percentile(fdir):

    nmonths = 12
    nslices = 5     # 1961-2010, in 10 year slices
    nyears = 10
    nrows = 360
    ncols = 720
    latitudes = np.linspace(-89.75, 89.75, nrows)
    longitudes = np.linspace(-179.75, 179.75, ncols)

    data = np.zeros((nslices, nyears, nmonths, nrows, ncols))
    i = 0
    for fname in glob.glob(os.path.join(fdir, '*.nc')):

        st = int(os.path.basename(fname).split(".")[2])
        en = int(os.path.basename(fname).split(".")[3])
        years = np.arange(st, en+1)

        ds = xr.open_dataset(fname)
        temp = ds["tmx"][:,:,:].values
        temp = temp.reshape(nyears, nmonths, nrows, ncols)
        data[i,:,:,:,:] = temp

        i += 1

    data = data.reshape(nslices*nyears*nmonths, nrows, ncols)

    p = np.percentile(data, 95, axis=0)

    # rejig this so we can mask the percentile map otherwise the see will look
    # bonkers
    x = data.reshape(nslices, nyears, nmonths, nrows, ncols)
    x = x[0,0,0,:,:]
    p = np.where(np.isnan(x), np.nan, p)

    return (p)

def get_flux_cru_percentiles(df, percentile, ofname):

    places = df.SiteCode.values
    lats_needed = df["SiteLatitude"].values
    lons_needed = df["SiteLongitude"].values

    # Fix units to match the 0.5 degree data
    lats_neededx = [float(x_round(float(i))) for i in lats_needed]
    lons_neededx = [float(x_round(float(i))) for i in lons_needed]
    lats_neededr = [float(i) for i in lats_needed]
    lons_neededr = [float(i) for i in lons_needed]

    lats_neededx = dict(zip(places, lats_neededx))
    lons_neededx = dict(zip(places, lons_neededx))
    lats_neededr = dict(zip(places, lats_neededr))
    lons_neededr = dict(zip(places, lons_neededr))

    if os.path.isfile(ofname):
        os.remove(ofname)
    fp  = open(ofname, "w")

    s = "%s,%s,%s,%s" % ("site","latitude","longitude","tmax_95th_percentile")
    print(s, end="\n", file=fp)

    nrows = 360
    ncols = 720
    latitudes = np.linspace(-89.75, 89.75, nrows)
    longitudes = np.linspace(-179.75, 179.75, ncols)

    for p in places:

        r = np.where(latitudes==lats_neededx[p])[0][0]
        c = np.where(longitudes==lons_neededx[p])[0][0]

        s = "%s,%s,%s,%s" % (p.strip(), lats_neededr[p],
                             lons_neededr[p], percentile[r,c])
        print(s, end="\n", file=fp)

    fp.close()

def x_round(x):
    # Need to round to nearest .25 or .75 to match the locations in CRU
    val = round(x * 4.0) / 4.0
    valx = str(val).split(".")
    v1 = valx[0]
    v2 = valx[1]

    if v2 <= "25":
        v2 = "25"
    else:
        v2 = "75"
    valx = float("%s.%s" % (v1, v2))

    return (valx)

if __name__ == "__main__":

    df = pd.read_csv("data/site_data/Site_metadata.csv")
    fdir = "data/TMAX"

    p = calculate_tmax_percentile(fdir)
    #plt.imshow(p, origin='lower')
    #plt.colorbar()
    #plt.show()
    #data.tofile("percentile_TMAX_1960_2010.bin")

    ofname = "flux_site_95th_percentile.csv"
    get_flux_cru_percentiles(df, p, ofname)
