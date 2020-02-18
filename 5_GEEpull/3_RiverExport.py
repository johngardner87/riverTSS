# -*- coding: utf-8 -*-
"""
Created on Mon May 13 14:16:42 2019

@author: john

"""

import time
import ee
import os
import numpy
import pandas
#import feather

ee.Initialize()

#Note: This script uses python 3 which has syntax differences from python 2.7

#Source necessary functions.
#execfile('D:/Dropbox/projects/TSS/code/GEE_pull_functions_rivers_ST_0513.py')
#execfile('GEE_pull_functions_rivers_ST_0513.py')
exec(open("D:/Dropbox/projects/TSS/code/GEE_pull_functions_rivers.py").read())    

#Load in Pekel water occurance Layer and Landsat Collections.
# choose Pekel = False to have dynamic DWSE water mask
PekelMask = False
pekel = ee.Image('JRC/GSW1_0/GlobalSurfaceWater')

l8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
l7 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
l5 = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')

#Identify collection for use in sourced functions.
collection = 'SR'

#Standardize band names between the various collections and aggregate 
#them into one image collection

bn8 = ['B2','B3', 'B4', 'B5', 'B6','B7', 'pixel_qa']
bn57 = ['B1', 'B2', 'B3', 'B4', 'B5','B7', 'pixel_qa']
bns = ['Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2', 'qa']
  
ls5 = l5.select(bn57, bns)
ls7 = l7.select(bn57, bns)
ls8 = l8.select(bn8, bns)

# set cloud threshold. Should turn down to 50 since I filter to that later
ls = ee.ImageCollection(ls5.merge(ls7).merge(ls8))\
.filter(ee.Filter.lt('CLOUD_COVER', 60))

#Select the occurence layer in the pekel mask, which is just the 
#percentage of water occurence over a given pixel from 1985-2015.
#Set the percent occurance threshold and create a watermask from the result.
threshold = 80
water = pekel.select('occurrence').gt(threshold)
water = water.updateMask(water)

#lakes = ee.FeatureCollection('ft:1LkWMEc13PpYjUKZiSbUaboLGugOLAmPj3KyjCadk') #

# load reach polygons with same ID and footprint as NHD reach, but wider.
lakes = ee.FeatureCollection("users/johngardner87/CONUS_nhd_reach_poly_collapse");
#load nhd centerlines
centerline = ee.FeatureCollection("users/johngardner87/nhd_grwl_collapse");

# Filter lakes that have already been downloaded
lakesort = lakes.sort('ID')

lakeID = lakesort.aggregate_array('ID').getInfo() 

# make a folder in your google drive manually to output data
dlDir = 'D:/GoogleDrive/WQP_RiverExport_v2' 
filesDown = os.listdir(dlDir)  # -->
filesDown = [int(i.replace(".csv", "")) for i in filesDown]

lakeID  = [i for i in lakeID if i not in filesDown]

print(len(lakeID))

                    
for x in range(0,len(lakeID)):

# when splitting over multiple accounts I change the loop range
#for x in range(0,15000):
  lake = ee.Feature(lakes.filter(ee.Filter.eq('ID', lakeID[x])).first())
  reach = ee.Feature(centerline.filter(ee.Filter.eq('ID', lakeID[x])).first())

  shell = ee.Feature(None)
  
  #FilterBounds for lake, update masks for water occurence, clouds, roads, etc.
  #Remove any images with clouds directly over the waterbody
  lsover = ls.filterBounds(reach.geometry())\
  .map(clipImage)
  
  # Load dem to attach elevation to each reach. IF using GRWL lines, do NOT do this.
  # GRWL has good elevation and this is slow to run. But would have to edit all code below if not extracting elev.
  # And this only needs to be run once if you want elevation. I need to make a if/else option here.
  DEMout = ee.Image("users/eeProject/MERIT").reduceRegion(ee.Reducer.min(), lake.geometry(), 30)
  #DEMdata = DEMshell.set({'elevation': DEMout.get('elevation')})\
  #.set({'ID':lake.get('ID')})
  
  DEMdata = ee.Feature(None, {'elevation': DEMout.get('elevation'),'ID':lake.get('ID')})

  ## Map over sites within specific path/row and pull reflectance values, change cSCore for rivers threshold to 1000
  data = lsover.map(lakePull).filter(ee.Filter.lt('cScore', 1000))
  
  left = ee.Join.inner().apply(data, DEMdata, ee.Filter.equals("ID", None, "ID")).map(join) 
  
  dataOut = ee.batch.Export.table.toDrive(collection = left, \
                                            description = str(lakeID[x]),\
                                            folder = 'WQP_RiverExport_v2',\
                                            selectors = ['system:index','Cloud_Cover', 'ID','azimuth', 'blue', 'cScore',
                                                         'date', 'dswe', 'dswe_sd', 'elevation', 'green', 'hillshade',
                                                         'hillshadow', 'hillshadow_sd', 'nir', 'path', 'pixelCount',
                                                         'qa', 'qa_sd', 'red', 'row', 'sat', 'swir1', 'swir1_sd', 'swir2',
                                                         'zenith'],\
                                            fileFormat = 'csv')

  print(lakeID[x])

  #Check how many existing tasks are running and take a break if it's >15  
  maximum_no_of_tasks(15, 60)
  #Send next task.
  dataOut.start()
  print('done')

#Make sure all Earth engine tasks are completed prior to moving on.  
maximum_no_of_tasks(1,300)
print('done')
