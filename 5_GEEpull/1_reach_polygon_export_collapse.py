# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 14:06:13 2019
Exports reach polygon shapefiles given input of centerlines
add watershed boundaries: NHDplus trimmed and HUC04
@author: john
"""

# Import the Earth Engine Python Package
import time
import ee
#import os

# Initialize the Earth Engine object, using the authentication credentials.
ee.Initialize()

#load centerlines
nhd_gee = ee.FeatureCollection("users/johngardner87/nhd_grwl_collapse");

#load some boundary to iterate over so because it cannot do many reaches at once
# maybe can do 5-15 reaches each iteration. Could just use subsets of connected reaches instead
# of watershed boundaries and would be faster

HUC04 = ee.FeatureCollection("USGS/WBD/2017/HUC04")

#HUC08 = ee.FeatureCollection("USGS/WBD/2017/HUC08")

### define functions for finding pixels  closest to each reach in centerline file
# change the "ID" to whatever the unique ID is of each reach

def fdtFun(f):
    dt = ee.Image.constant(0).byte().paint(ee.FeatureCollection(f), 1).fastDistanceTransform()
    return dt.copyProperties(f, ['ID'])

def mfun(c, p):
    return ee.Image(c).min(p)

def minFun(f):
  return ee.Image(f.iterate(mfun, ee.Image(f.first())))

def idFun(c, p):
    c = ee.Image(c)
    p = ee.Image(p)
    return p.where(c.eq(minMap), ee.Number(c.get('ID')))

# functions for bacth processing
def maximum_no_of_tasks(MaxNActive, waitingPeriod):
  ##maintain a maximum number of active tasks
  time.sleep(10)
  ## initialize submitting jobs
  ts = list(ee.batch.Task.list())

  NActive = 0
  for task in ts:
       if ('RUNNING' in str(task) or 'READY' in str(task)):
           NActive += 1
  ## wait if the number of current active tasks reach the maximum number
  ## defined in MaxNActive
  while (NActive >= MaxNActive):
      time.sleep(waitingPeriod) # if reach or over maximum no. of active tasks, wait for 2min and check again
      ts = list(ee.batch.Task.list())
      NActive = 0
      for task in ts:
        if ('RUNNING' in str(task) or 'READY' in str(task)):
          NActive += 1
  return()

## loop through all HUC04 and output reach 
hucID = HUC04.aggregate_array('huc4').getInfo()


for x in range(0,len(hucID)):
#for x in range(0,2):
    
    #huc = HUC04.filter(ee.Filter.eq('huc4' , hucID[x]))
    huc = HUC04.filter(ee.Filter.stringStartsWith('huc4' , hucID[x]))
    
    nhdClip = nhd_gee.filterBounds(huc.geometry())

# set maximum buffer from center line to assign pixels a reach ID (2000 m)
    # when I do global pull I am going to jack this up to 5000-10000
    nhdClip_buffer = nhdClip.geometry().buffer(2000)

##
    fdt = nhdClip.map(fdtFun)
    
    minMap = minFun(fdt)
    
    idImage = fdt.iterate(idFun, ee.Image(fdt.first()))
    
# make image and clip to max distance buffer from centerline for vector export
    COMID = ee.Image(idImage).rename('ID').clip(nhdClip_buffer)
##
    vectors = COMID.addBands(COMID).reduceToVectors(
        crs = COMID.projection(), 
        scale =  60, 
        geometry = huc, 
        geometryType = 'polygon',
        eightConnected= False,
        tileScale = 16,
        bestEffort = True,
        labelProperty = 'ID',
        maxPixels = 3000000000000000, 
        reducer = ee.Reducer.first()
        )

    
    dataOut = ee.batch.Export.table.toDrive(collection = vectors, \
                                              description = str(hucID[x]),\
                                              folder = 'reachPolygons_collapse',\
                                              fileFormat = 'shp')
    
    print(hucID[x])
  #Check how many existing tasks are running and take a break if it's >15  
    maximum_no_of_tasks(15, 60)
  #Send next task.
    dataOut.start()

#Make sure all Earth engine tasks are completed prior to moving on.  
maximum_no_of_tasks(1,300)
print('done')

### NOTES
## watersheds missed due to exceeeded memory = 1701, 1704, 1706, 1018, 1101, 0708, 0601

### Scratch text
#    fdt = nhdClip.map(function(f) {
#    dt = ee.Image.constant(0).byte().paint(ee.FeatureCollection(f), 1).fastDistanceTransform()
#    return(dt.copyProperties(f, ['COMID']))})

#    minMap = ee.Image(fdt.iterate(function(c, p) {return(ee.Image(c).min(p))}, ee.Image(fdt.first())))

#    idImage = fdt.iterate(function(c, p) {
#    c = ee.Image(c)
#    p = ee.Image(p)
#    return(p.where(c.eq(minMap), ee.Number(c.get('COMID'))))},
#    ee.Image(fdt.first()))
