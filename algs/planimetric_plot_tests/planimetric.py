import gmplot

# Replace these coordinates with your GPS data
latitude_list = [37.7749, 37.7749, 37.7749]
longitude_list = [-122.4194, -122.4194, -122.4194]

# Initialize the map
gmap = gmplot.GoogleMapPlotter(latitude_list[0], longitude_list[0], 13)

# Plot the points on the map
gmap.scatter(latitude_list, longitude_list, '#FF0000', size=40, marker=False)

# Draw the map to an HTML file
gmap.draw("map.html")

# Open the map in a web browser
import webbrowser
webbrowser.open_new("map.html")

from simplekml import Kml

# Replace these coordinates with your GPS data
coordinates = [
    (37.7749, -122.4194),
    (37.7751, -122.4196),
    (37.7753, -122.4198)
]

# Create a KML object
kml = Kml()

# Create a KML linestring and add coordinates
linestring = kml.newlinestring(name="GPS Coordinates")
linestring.coords = coordinates

# Save the KML file
kml.save("gps_coordinates.kml")
