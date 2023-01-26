from skyfield.api import load, wgs84
from skyfield.constants import (AU_M, ERAD, DEG2RAD, IERS_2010_INVERSE_EARTH_FLATTENING, tau)
from skyfield.units import Angle
from skyfield.positionlib import ICRF, Geocentric

# Mapping
import geojson
import geog
import shapely.geometry

# Matpotlib
import matplotlib.pyplot as plt
from descartes import PolygonPatch

import numpy as np

from urllib.parse import quote
import logging

satNames = {
    "NOAA 18",
    "GOES 16",
    "GOES 17"
}

colors = {
    "LAND": "#335533",
    "WATER": "#333355",
    "BG": "#333333",
    "SAT_MKR": "#66AA66",
}

dataUrl = 'http://celestrak.org/NORAD/elements/gp.php?'

# Skyfield Timescale
ts = load.timescale()

# Map elements (globals so we only have to load them once)
geoLand = None
with open("geo/WB_Land_10m_lowres.geojson") as json:
    geoLand = geojson.load(json)

def queryTLEs(list):
    logging.info("Loading satellite TLEs")
    # Empty list of satellites
    sats = []
    # Load each satellite, requesting a reload of TLE data if it's out of date
    for satName in list:
        filename = 'tle-{}.txt'.format(satName)
        # Prepare the URL using urllib to replace invalid characters
        url = '{}NAME={}'.format(dataUrl,quote(satName))
        # Load the TLE, either from disk or download
        sat = load.tle_file(url,filename=filename)[0]
        # Check date
        days = ts.now() - sat.epoch
        if abs(days) > 14:
            sat = load.tle_file(url, filename=filename, reload=True)
            logging.warn("Satellite {} TLE was out of date, reloaded".format(satName))
        sats.append(sat)

    return sats

def getLatLonAlt(satellite):

    # Get geocentric position
    geo = satellite.at(ts.now())
    # Get wgs84 coords
    lat, lon = wgs84.latlon_of(geo)
    height = wgs84.height_of(geo).km

    return lat, lon, height

def plotMap():
    # Create mpl figure and get axis
    fig = plt.figure()
    ax = fig.gca()

    # Set map styling
    ax.set_facecolor(colors["WATER"])
    fig.set_facecolor(colors["BG"])
    
    # Margins
    plt.subplots_adjust(left=0,bottom=0,right=1,top=1,wspace=0,hspace=0)
    ax.margins(0)

    # Plot land
    for feature in geoLand['features']:
        ax.add_patch(PolygonPatch(feature['geometry'], fc=colors['LAND'], ec=colors["LAND"]))

    ax.axis('scaled')

    return fig, ax

def plotSats(fig, ax, satellites):

    # Plot Satellite Points & Coverage Footprint
    for i, sat in enumerate(satellites):
        # Get Sat Coords
        lat, lon, alt = getLatLonAlt(sat)
        # Calculate footprint circle

        # Plot & Annotate
        ax.plot(lon.degrees, lat.degrees, marker='2', markersize=15, color=colors["SAT_MKR"])
        ax.annotate(sat.name, (lon.degrees, lat.degrees), textcoords="offset points", xytext=(0,10), ha='center', color=colors["SAT_MKR"])

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Load satellites
    satellites = queryTLEs(satNames)
    logging.info("Loaded satellites:")
    for sat in satellites:
        logging.info("    {}".format(sat))

    # Get satellite positions
    for sat in satellites:
        lat, lon, alt = getLatLonAlt(sat)
        logging.info("Satellite {}: {}, {}, alt: {:.0f} km".format(sat, lat, lon, alt))

    # Display Map
    fig, ax = plotMap()
    plotSats(fig, ax, satellites)
    plt.show()