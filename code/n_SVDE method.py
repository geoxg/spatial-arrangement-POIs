#!/usr/bin/env python
# coding: utf-8

# In[ ]:


## operating system 
import os 

## geospatial libraries 
import geopandas as gpd
from shapely.geometry import MultiPoint
from shapely.geometry import LineString, MultiPoint, Polygon
from shapely.ops import cascaded_union

# Mathematical and Statistical Libraries
import math 
import numpy as np
import collections
from array import array 

# Voronoi and Spatial Entropy Calculation
from geovoronoi import coords_to_points, points_to_coords, voronoi_regions_from_coords, calculate_polygon_areas

## data handling
import pandas as pd

## Logging and debugging
import logging
logging.basicConfig(level=logging.INFO)
geovoronoi_log = logging.getLogger('geovoronoi')
geovoronoi_log.setLevel(logging.INFO)
geovoronoi_log.propagate = True

import warnings
warnings.filterwarnings("ignore")

gpd.__version__


# In[ ]:


# Read the shapefile for the study area (e.g., Cincinnati study area)
studyArea = gpd.read_file(r'...\shapefile\Cincy_studyArea.shp')

# Display the GeoDataFrame containing the study area data (optional)
# This will show the content of the study area, like its geometries and attributes
studyArea  


# In[ ]:


'''
# export each row of polygon to seprate polygon (use if necessary)
for i in studyArea.index:
    subdf = studyArea.iloc[[i]] 
    subdf.to_file('...\Output_divide\CincyStudyUnit_{0}.shp'.format(i)) 
'''


# In[ ]:


# Function to calculate the Shannon Voronoi Entropy for a set of polygon areas

def calculate_shannon_Voronoi(poly_areas):
    b2 = []  # List to store individual entropy contributions for each polygon
    # Iterate through each polygon area in the poly_areas dictionary
    for index, (key, value) in enumerate(poly_areas.items()):
        # Compute the Shannon entropy for each polygon area
        # The formula used is: - (p_i * log(p_i)) where p_i is the proportion of the area
        b_value = (poly_areas[index] / sum(poly_areas.values())) * np.log(poly_areas[index] / sum(poly_areas.values()))
        b2.append(b_value)  # Append the calculated entropy value for each polygon

    # Sum the individual entropy values to get the total entropy (w1)
    w1 = np.sum(np.array(b2))
    
    # Multiply by -1 to adjust for the entropy calculation (since entropy is negative)
    w1 = w1 * (-1)
    
    # Return the final Shannon Voronoi entropy value
    return w1


# In[ ]:


### Setup and initializations ###

# Read the shapefile of Points of Interest (POI) related to Public Transportation
POI_merge = gpd.read_file(r'...\POI_reclassify\PublicTransportation.shp')

# Create empty lists to store results
w_list = []  # Placeholder for storing entropy values (not used in the code, possibly removed later)
studyUnit_list = []  # List to store results for each study unit

# Loop through each study unit (e.g.,Cincinnati block groups, indexed from 0 to 286)
for i in range(0, 287):
    # Define the path for the shapefile corresponding to the study unit (Cincinnati Block Group)
    path = r'...\Output_divide\CincyStudyUnit_'
    
    # Construct the file path for the study unit (e.g., CincyStudyUnit_i.shp)
    file_path = path + str(i) + '.shp'
    
    # Read the shapefile for the current study unit (e.g., Census Block Group, Census Tract )
    studyUnit = gpd.read_file(file_path)

    # Clip the POI dataset with the current study unit to extract POIs within the block group boundary
    Point_clipped = gpd.clip(POI_merge, studyUnit)
    
    # Extract the coordinates (x, y) of the clipped points for further analysis
    locations = [(point.x, point.y) for point in Point_clipped.geometry]
    
    # Optionally, save the clipped points to a new shapefile for future use in the ANN_ratio method
    output_clipped_points = f"...\\Cincy_studyUnit\\Output_POI_Points\\output_clipped_points_{i}.shp"
    Point_clipped.to_file(output_clipped_points)  # Save the clipped POIs
    
    ### Handling cases where there are no or lees than three POIs in the study unit ###
    if len(Point_clipped) == 0:
        w1 = -1  # Assign a placeholder value for entropy when no POIs are present
        studyUnit['nSVDE_PubTr'] = w1  # Add the placeholder value to the study unit
        studyUnit_list.append(studyUnit)  # Append the study unit to the result list
        
    elif len(Point_clipped) == 1:
        w1 = -2  # Assign a different placeholder value when there is exactly one POI
        studyUnit['nSVDE_PubTr'] = w1  # Update the study unit with the new value
        studyUnit_list.append(studyUnit)  # Append the study unit to the result list
        
    elif len(Point_clipped) == 2:
        w1 = -3  # Assign another placeholder value when there are exactly two POIs
        studyUnit['nSVDE_PubTr'] = w1  # Update the study unit with the new value
        studyUnit_list.append(studyUnit)  # Append the study unit to the result list
        
    else:
        # If there are more than two POIs, proceed with the spatial analysis
        boundary_shape = cascaded_union(studyUnit.geometry)  # Combine the geometries of the study unit (polygon)
        
        # Convert the clipped POI geometries to a list of coordinates for Voronoi analysis
        coords = points_to_coords(Point_clipped.geometry)
        
        # Calculate the area of the study unit in square kilometers
        studyUnit["area"] = studyUnit['geometry'].area / 10**6  # Convert from square meters to square kilometers
        
        try:
            # Generate Voronoi regions from the POI coordinates within the study unit boundary
            region_polys, region_pts = voronoi_regions_from_coords(coords, boundary_shape, per_geom=False)
            
            # Calculate the areas of the Voronoi polygons (regions) in square kilometers
            poly_areas = calculate_polygon_areas(region_polys, m2_to_km2=True)  # Convert from square meters to square kilometers
            
            # Calculate the Shannon entropy (SVDE) for the Voronoi polygons
            w1 = calculate_shannon_Voronoi(poly_areas)
            
            # Calculate the maximum possible SVDE value assuming the POIs are evenly distributed
            p = 1 / len(Point_clipped)  # Probability of each POI if evenly distributed
            w2 = (-(p * math.log(p))) * len(Point_clipped)  # Maximum entropy formula
            
            # Normalize the SVDE value to remove the imapcts of numbers of POIs to the entropy value
            studyUnit['nSVDE_PubTr'] = w1 / w2  # Normalized entropy (n_SVDE)
            studyUnit_list.append(studyUnit)  # Append the updated study unit to the result list
            
        except:
            # In case of any errors (e.g., if Voronoi analysis fails), print the file path
            print(file_path)


# In[ ]:


# Concatenate all GeoDataFrames in the 'studyUnit_list' into a single DataFrame
# The 'ignore_index=True' ensures that the index is reset in the concatenated DataFrame
studyArea_entropy = pd.concat(studyUnit_list, ignore_index=True)

# Generate summary statistics 
studyArea_entropy.describe()


# In[ ]:


# Save the concatenated GeoDataFrame 'studyArea_entropy' as a shapefile
studyArea_entropy.to_file(
    '...\Output',  # Path to save the output shapefile (replace with actual path)
    driver='ESRI Shapefile',  # Specifies the format of the file to be saved (Shapefile in this case)
    layer='nSVDE for Public Transportation',  # The name of the layer to be saved within the shapefile
    encoding='utf-8'  # Specify the encoding to use for the attribute data (UTF-8)
)

