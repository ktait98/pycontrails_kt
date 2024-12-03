import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Define the locations
locations = [
    {"name": "NA", "lat": 47.5, "lon": -32.5, "alt": 12500, "ha": "center", "va": "bottom", "x_offset": 0, "y_offset": 5},
    {"name": "US", "lat": 37.5, "lon": -97.5, "alt": 11500, "ha": "left", "va": "top", "x_offset": 3, "y_offset": -3},
    {"name": "EU", "lat": 42.5, "lon": 7.5, "alt": 9500, "ha": "left", "va": "top", "x_offset": 3, "y_offset": -3},
    {"name": "SA", "lat": -27.5, "lon": -67.5, "alt": 13500, "ha": "left", "va": "top", "x_offset": 3, "y_offset": -3},
    {"name": "SEA", "lat": 22.5, "lon": 102.5, "alt": 10500, "ha": "left", "va": "top", "x_offset": 3, "y_offset": -3},
]

# Define professional colors for each location
colors = ['#1f77b4', '#ff7f0e', '#66c2a5', '#d62728', '#9467bd']  # Blue, Orange, Green, Red, Purple

# Create a figure with a cartopy map projection
fig = plt.figure(figsize=(10, 5))
ax = plt.axes(projection=ccrs.PlateCarree())

# Add features to the map
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAKES, alpha=0.5)
ax.add_feature(cfeature.RIVERS)

for loc, color in zip(locations, colors):
    ax.plot(loc["lon"], loc["lat"], marker='o', color=color, markersize=5, transform=ccrs.PlateCarree())
    ax.text(loc["lon"] + loc["x_offset"], loc["lat"] + loc["y_offset"], 
            f"{loc['name']}\n({loc['lat']}, {loc['lon']})\nAlt: {loc['alt']} m",
            horizontalalignment=loc["ha"], verticalalignment=loc["va"], fontsize=8, transform=ccrs.PlateCarree(),
            bbox=dict(facecolor=color, alpha=0.9, edgecolor='none'))

# Set the extent of the map
ax.set_extent([-180, 180, -90, 90])

# Add gridlines
ax.gridlines(draw_labels=True)

# Add a title
plt.title('Hotspot Locations')

# Show the plot
plt.show()

# Save the plot
plt.savefig('/user/home/kt16229/work/pycontrails_kt/pycontrails/models/gpat/outputs/plots/hotspot_locations.png', format='png')
plt.savefig('/user/home/kt16229/work/pycontrails_kt/pycontrails/models/gpat/outputs/plots/hotspot_locations.pdf', format='pdf')