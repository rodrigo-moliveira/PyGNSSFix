""" Module for the planimetric plots """

import simplekml
import folium
import webbrowser
import os


def create_kml_file(points, filename):
    """
    Creates a KML file from a list of points, defined in geodetic coordinates (latitude, longitude and height)

    Note: this file may then be loaded in KML application readers, such as the Google Earth.

    Args:
        points (list[tuple[float, float, float]]): geodetic coordinates of the points to be saved in the KML file,
            This is a list of tuples as [latitude (degree), longitude (degree), height (meters)].
        filename (str): The filename of the KML file to be created.
    """
    kml = simplekml.Kml()
    for lat, lon, alt in points:
        pnt = kml.newpoint(coords=[(lon, lat, alt)])
        pnt.altitudemode = simplekml.AltitudeMode.absolute
    kml.save(filename)


def create_map_html_file(points, filename):
    """
    Plots points on a world map and saves as an HTML file.

    Args:
        points (list[tuple[float, float]]): geodetic coordinates of the points to be saved in the KML file,
            This is a list of tuples as [latitude (degree), longitude (degree)].
        filename (str): The filename of the HTML file to be created.
    """
    # Create a map centered around the average latitude and longitude
    avg_lat = sum(point[0] for point in points) / len(points)
    avg_lon = sum(point[1] for point in points) / len(points)
    my_map = folium.Map(location=[avg_lat, avg_lon], zoom_start=2)

    for lla in points:
        folium.Marker([lla[0], lla[1]]).add_to(my_map)

    # Save the map to an HTML file
    my_map.save(filename)

    # Open the HTML file in the default web browser
    webbrowser.open('file://' + os.path.realpath(filename))
