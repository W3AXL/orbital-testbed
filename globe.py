from turtle import shape
from pandas import interval_range
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
import matplotlib as mpl
from matplotlib import animation
from descartes import PolygonPatch

import numpy as np

import math

import json

import sys

from urllib.parse import quote
import logging

sats = []

mapColors = {
    "LAND": "#223322",
    "LAND_LINE": "#112211",
    "WATER": "#222233",
    "BG": "#111111",
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
with open("geo/WB_Land_10m_lowres.geojson") as jsondata:
    geoLand = geojson.load(jsondata)

# Global figure & axis
fig = None
ax = None
satPoints = []  # XY points of the satellites
satAnns = []    # Annotations for satellite names
satPrints = []  # Footprint circles of the satellites

def importSats(filename):
    global sats
    with open(filename) as file:
        parsed = json.load(file)
        for sat in parsed:
            sats.append(sat)

def queryTLEs():
    global sats

    logging.info("Loading satellite TLEs")
    # Empty list of satellites
    # Load each satellite, requesting a reload of TLE data if it's out of date
    
    for idx, sat in enumerate(sats):
        try:
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
            # Update object in list
            sats[idx]["satobj"] = sat
        except IndexError:
            logging.error("Error loading TLE for satllite {}. Check spelling?".format(sat["name"]))
            exit(1)

def getLatLonAlt(satellite):

    # Get geocentric position
    geo = satellite.at(ts.now())
    # Get wgs84 coords
    #lat, lon = wgs84.latlon_of(geo)
    geoPos = wgs84.subpoint_of(geo)
    lat = geoPos.latitude
    lon = geoPos.longitude
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

    point = [lon, lat]
    angles = np.linspace(0, 360, n_points)
    polygon = geog.propagate(point, angles, radius)
    shape = shapely.geometry.mapping(shapely.geometry.Polygon(polygon))
    #print(shape)
    return shape

def plotMap():

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

def in_frame(annotation):
    ann_bbox = annotation.get_window_extent()
    #print(ann_bbox)
    plt_bbox = plt.gcf().get_window_extent()
    #print(plt_bbox)
    return plt_bbox.contains(ann_bbox.x1, ann_bbox.y0)

def updateSatPos():
    # Update satellite positions
    for i, sat in enumerate(sats):
        # Get Sat Coords
        lat, lon, alt = getLatLonAlt(sat["satobj"])
        # Save to list
        sats[i]["pos"] = [lat, lon, alt]

def init():
    global fig
    global ax
    global sats
    global satPrints
    global satPoints
    global satAnns

    # Draw map
    plotMap()

    # Get initial satellite positions
    updateSatPos()

    # Create initial points for each satellite and add to list
    for sat in sats:
        # Get position & band
        lat = sat["pos"][0].degrees
        lon = sat["pos"][1].degrees
        alt = sat["pos"][2]
        band = sat["band"]

        # Plot point and annotation
        point = ax.plot(lon, lat, marker='2', markersize=15, color=bandColors[band])[0]
        ann = ax.annotate(sat["name"], (lon, lat), textcoords="offset points", xytext=(6,2), ha='left', color=bandColors[band], fontsize=10, fontweight='medium')
        satPoints.append(point)
        satAnns.append(ann)

        # Calculate footprint circle
        #radius = getCoverageRadius(alt)
        # Generate circle geo object
        #circle = generateCircle(lat, lon, radius)
        # Plot circle
        #print = ax.add_patch(PolygonPatch(circle, fc="#00000000", ec=mapColors["SAT_CVG_LINE"]))
        #satPrints.append(print)

    return satPoints, satAnns, satPrints

def update(frame):
    global fig
    global ax
    global sats
    global satPoints
    global satAnns
    global satPrints

    # Get new positions
    updateSatPos()

    for idx, point in enumerate(satPoints):
        # Get new position
        lat = sats[idx]["pos"][0].degrees
        lon = sats[idx]["pos"][1].degrees
        alt = sat["pos"][2]
        # Update marker position
        point.set_data(lon, lat)
        # Update annotation position
        satAnns[idx].xy = (lon, lat)
        satAnns[idx].set(ha='left')
        satAnns[idx].set_position((6,2))
        # Flip annotation if needed
        if not in_frame(satAnns[idx]):
            satAnns[idx].set(ha='right')
            satAnns[idx].set_position((-6,2))  

        # Calculate footprint circle
        #radius = getCoverageRadius(alt)
        # Generate circle geo object
        #circle = generateCircle(lat, lon, radius)
        # update circle
        #satPrints[idx] = PolygonPatch(circle)

    return satPoints, satAnns, satPrints

def plotSats():
    # Update satellite positions
    for i, sat in enumerate(sats):
        # Get Sat Coords
        lat, lon, alt = getLatLonAlt(sat["satobj"])
        # Save to list
        sats[i]["pos"] = [lat, lon, alt]

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
        ax.annotate(sat["name"], (lon.degrees, lat.degrees), textcoords="offset points", xytext=(6,2), ha='left', color=satColor, fontsize=10, fontweight='medium')

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Load sat json
    importSats(sys.argv[1])

    # Load satellite TLEs
    queryTLEs()
    logging.info("Loaded satellites:")
    for sat in sats:
        logging.info("    {}".format(sat["satobj"]))

    # Disable MPL UI elements
    mpl.rcParams['toolbar'] = 'None'

    # Create mpl figure and get axis
    fig = plt.figure()
    ax = fig.gca()
    
    # Animate
    anim = animation.FuncAnimation(fig, update, init_func=init, interval=2000)

    plt.show()