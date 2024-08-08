import simplekml
import folium
import webbrowser
import os


def create_kml(points, filename):
    """
    Creates a KML file from a list of points with absolute altitude.

    :param points: List of tuples (latitude, longitude, altitude).
    :param filename: The filename of the KML file to be created.
    """
    kml = simplekml.Kml()
    for lat, lon, alt in points:
        pnt = kml.newpoint(coords=[(lon, lat, alt)])
        pnt.altitudemode = simplekml.AltitudeMode.absolute
    kml.save(filename)

# Example usage with latitude, longitude, and altitude points
points = [
    (37.7749, -122.4194, 500),  # Example point with altitude
    (34.0522, -118.2437, 1000),
    (40.7128, -74.0060, 1500)
]

# Save to KML
create_kml(points, "../../../algs/planimetric_plot_tests/points.kml")

def plot_points_on_map(points, filename):
    """
    Plots points on a world map and saves as an HTML file.

    :param points: List of tuples (latitude, longitude).
    :param filename: The filename of the HTML file to be created.
    """
    # Create a map centered around the average latitude and longitude
    avg_lat = sum(point[0] for point in points) / len(points)
    avg_lon = sum(point[1] for point in points) / len(points)
    my_map = folium.Map(location=[avg_lat, avg_lon], zoom_start=2)

    for lat, lon in points:
        folium.Marker([lat, lon]).add_to(my_map)

    # Save the map to an HTML file
    my_map.save(filename)

    # Open the HTML file in the default web browser
    webbrowser.open('file://' + os.path.realpath(filename))

# Example usage with latitude and longitude points
latlon_points = [
    (37.7749, -122.4194),  # San Francisco
    (34.0522, -118.2437),  # Los Angeles
    (40.7128, -74.0060)    # New York
]

# Plot and save to HTML
plot_points_on_map(latlon_points, "../../../algs/planimetric_plot_tests/map.html")
