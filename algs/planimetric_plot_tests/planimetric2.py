import plotly.express as px
import pandas as pd

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
exit()

# Import data from USGS
data = pd.read_csv('https://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/all_month.csv')


# Drop rows with missing or invalid values in the 'mag' column
data = data.dropna(subset=['mag'])
data = data[data.mag >= 0]


# Create scatter map
fig = px.scatter_geo(data, lat='latitude', lon='longitude', color='mag',
                     hover_name='place', #size='mag',
                     title='Earthquakes Around the World')
fig.show()