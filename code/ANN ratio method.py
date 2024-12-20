#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, MultiPoint, Polygon
import os 


# In[ ]:


# read the shapefile of the study area
studyArea = gpd.read_file(r'...\shapefile\Cincy_studyArea.shp')

# The CRS information will help verify that the shapefile's spatial coordinates are in the correct projection.
print(studyArea.crs) 

# Display the GeoDataFrame containing the study area data (optional)
# This will show the content of the study area, like its geometries and attributes
studyArea


# In[ ]:


# export each row of polygon to seprate polygon
for i in studyArea.index:
    subdf = studyArea.iloc[[i]] 
    subdf.to_file('...\Output_divide\CincyStudyUnit_{0}.shp'.format(i))


# In[ ]:


# Define the folder containing the shapefiles
shapefile_folder = "...\\Cincy_studyUnit\\Output_POI_Points"

# List all shapefiles in the folder
shapefiles = [os.path.join(shapefile_folder, f) for f in os.listdir(shapefile_folder) if f.endswith('.shp')]

# Initialize a fresh results list for this run
results = []

# Function to calculate the expected average distance
def calculate_expected_distance(area_sq_km, point_count):
    if point_count == 0 or np.isnan(point_count):
        return np.nan
    return 2 * np.sqrt(area_sq_km / (point_count * np.pi))

# Loop through each shapefile
for shapefile in shapefiles:
    # Read the shapefile
    gdf = gpd.read_file(shapefile)
    
    # Ensure the GeoDataFrame has a projected CRS for accurate distance calculations
    if not gdf.crs.is_projected:
        gdf = gdf.to_crs("EPSG:3735")  # Replace with a suitable projected CRS for your region
    
    # Initialize variables for each shapefile to avoid overlap
    coords = []
    average_distance_km = 0
    
    # Check for empty geometries and non-point geometries
    if gdf.geometry.is_empty.all():
        average_distance_km = 0  # Assign zero for shapefiles with empty geometries
    else:
        # Extract coordinates of valid point geometries
        for geom in gdf.geometry:
            if geom.geom_type == 'Point':  # Ensure we're dealing with points
                coords.append([geom.x, geom.y])
            else:
                print(f"Warning: Skipping non-point geometry in {shapefile}.")
        
        if not coords:  # If no valid points exist
            average_distance_km = 0  # Assign zero if no valid points
        else:
            # Convert coordinates list to numpy array
            coords = np.array(coords)
            
            # Use KDTree for nearest neighbor calculation
            tree = cKDTree(coords)
            distances, _ = tree.query(coords, k=2)  # k=2 to include the point itself
            nearest_distances = distances[:, 1]  # Exclude the self-distance (index 0)
            
            # Compute the average nearest neighbor distance in feet
            average_distance_feet = nearest_distances.mean()
   
            # Convert to kilometers
            average_distance_km = average_distance_feet / 3280.84  # Convert feet to kilometers

    # Add area in square kilometers (convert from square feet to square kilometers)
    gdf['area_sq_km'] = gdf.geometry.area * 9.290304e-8  # Square feet to square kilometers

    # Add expected average distance (based on the area and point count)
    point_count = len(coords)
    expected_distance_km = calculate_expected_distance(gdf['area_sq_km'].sum(), point_count)

    # Add the average distance in kilometers to the GeoDataFrame as a new column
    gdf['avg_nearest_distance_km'] = average_distance_km

    # Calculate the ANN-ratio
    if expected_distance_km > 0:
        ANN_ratio = average_distance_km / expected_distance_km
    else:
        ANN_ratio = np.nan

    # Append the results to the list
    shapefile_name = os.path.basename(shapefile)
    results.append({
        'Shapefile': shapefile_name,
        'Average Nearest Neighbor Distance (km)': average_distance_km,
        'Expected Average Distance (km)': expected_distance_km,
        'ANN Ratio': ANN_ratio
    })

    # Print status for each shapefile
    print(f"Processed: {shapefile_name}")
    
# Create a DataFrame to hold the results
results_df = pd.DataFrame(results)

# Export results to an Excel file
output_excel = "...\\Excel\\results_ANN_ratio.xlsx"
if os.path.exists(output_excel):  # Remove existing Excel file to avoid conflicts
    os.remove(output_excel)
results_df.to_excel(output_excel, index=False)

# Final cleanup of variables
del shapefiles, results, gdf

print("Results for ANN_ratio saved to the exported Excel")

