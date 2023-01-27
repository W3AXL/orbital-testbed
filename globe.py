from turtle import shape
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

import math

from urllib.parse import quote
import logging

sats = [
    {"name": "NOAA 15", "band": "l"},
    {"name": "NOAA 18", "band": "l"},
    {"name": "NOAA 19", "band": "l"},
    #{"name": "GOES 15", "band": "l"},
    {"name": "GOES 16", "band": "l"},
    #{"name": "GOES 17", "band": "l"},
    {"name": "GOES 18", "band": "l"},
    {"name": "METEOR-M 2", "band": "l"},
    {"name": "METEOR-M2 2", "band": "l"},
    {"name": "METOP-B", "band": "l"},
    {"name": "METOP-C", "band": "l"},
    {"name": "DMSP 5D-3 F18", "band": "s"},
    {"name": "DMSP 5D-3 F17", "band": "s"},
]

mapColors = {
    "LAND": "#335533",
    "LAND_LINE": "#224422",
    "WATER": "#333355",
    "BG": "#333333",
    "SAT_CVG_LINE": "#FFFFFF55",
    "SAT_CVG_FILL": "#FFFFFF11",
}

bandColors = {
    "l": "#77CCCC",
    "s": "#77CC77",
    "c": "#CCCC77",
    "x": "#CC7777",
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
    for sat in list:
        filename = 'tle-{}.txt'.format(sat["name"])
        # Prepare the URL using urllib to replace invalid characters
        url = '{}NAME={}'.format(dataUrl,quote(sat["name"]))
        # Load the TLE, either from disk or download
        logging.debug("Querying {}".format(url))
        sat = load.tle_file(url,filename=filename)[0]
        # Check date
        days = ts.now() - sat.epoch
        if abs(days) > 14:
            sat = load.tle_file(url, filename=filename, reload=True)
            logging.warn("Satellite {} TLE was out of date, reloaded".format(sat["name"]))
        sats.append(sat)

    return sats

def getLatLonAlt(satellite):

    # Get geocentric position
    geo = satellite.at(ts.now())
    # Get wgs84 coords
    lat, lon = wgs84.latlon_of(geo)
    alt = wgs84.height_of(geo)

    return lat, lon, alt

def getCoverageRadius(alt):
    # Get c (satellite altitude in m plus earth's radius in m)
    c = alt.m + ERAD
    # Get length of n (small side of triangle)
    n = c - (2*ERAD*alt.m + alt.m**2)/c
    # Get coverage radius in meters
    return math.sqrt((2*n*alt.m*ERAD + n*alt.m)/c)

def generateCircle(lat, lon, radius, n_points=64):

    point = [lon.degrees, lat.degrees]
    angles = np.linspace(0, 360, n_points)
    polygon = geog.propagate(point, angles, radius)
    shape = shapely.geometry.mapping(shapely.geometry.Polygon(polygon))
    #print(shape)
    return shape

def plotMap():
    # Create mpl figure and get axis
    fig = plt.figure()
    ax = fig.gca()

    # Set map styling
    ax.set_facecolor(mapColors["WATER"])
    fig.set_facecolor(mapColors["BG"])
    
    # Margins
    plt.subplots_adjust(left=0,bottom=0,right=1,top=1,wspace=0,hspace=0)
    ax.margins(0)

    # Plot land
    for feature in geoLand['features']:
        ax.add_patch(PolygonPatch(feature['geometry'], fc=mapColors['LAND'], ec="#00000000"))

    ax.axis('scaled')

    return fig, ax

def plotSats(fig, ax, satellites):

    # Plot Satellite Points & Coverage Footprint
    for i, sat in enumerate(satellites):
        # Get Sat Coords
        lat, lon, alt = getLatLonAlt(sat)

        # Get coverage radius
        #radius = getCoverageRadius(alt)
        # Generate circle geo object
        #circle = generateCircle(lat, lon, radius)
        # Plot circle
        #ax.add_patch(PolygonPatch(circle, fc="#00000000", ec=colors["SAT_CVG_LINE"]))

        # Get band-based color
        satColor = bandColors[sats[i]["band"]]

        # Plot & Annotate Points
        ax.plot(lon.degrees, lat.degrees, marker='2', markersize=15, color=satColor)
        ax.annotate(sat.name, (lon.degrees, lat.degrees), textcoords="offset points", xytext=(6,2), ha='left', color=satColor, fontsize=10, fontweight='medium')

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Load satellites
    satellites = queryTLEs(sats)
    logging.info("Loaded satellites:")
    for sat in satellites:
        logging.info("    {}".format(sat))

    # Get satellite positions
    for sat in satellites:
        lat, lon, alt = getLatLonAlt(sat)
        logging.info("Satellite {}: {}, {}, alt: {:.0f} km, classification: {}".format(sat.name, lat.degrees, lon.degrees, alt.km, sat.model.classification))

    # Display Map
    fig, ax = plotMap()
    plotSats(fig, ax, satellites)
    plt.show()