#!/usr/bin/env python3

#from scipy.io import netcdf  #For NetCDF file reading
import h5py                  #For HDF5 file reading
import sys                   #For command line arguments and exiting
import numpy as np
#import io
import glob                   #Used for listing files in a directory
import os                     #For gettting the server root directory
import bisect                 #For finding where lat/lon are in a sorted list
import itertools              #For mashing years with months
import datetime               #For managing time manipulation
import time                   #For managing time manipulation
import fnmatch                #For filtering file names
import argparse               #For handling command line magic
import code

#import code #For debugging with: code.interact(local=locals())

def binaryfind(arr, val):
  if val<arr[0] or val>arr[-1]:
    raise Exception('Value not in range of array')

  i = bisect.bisect_left(arr, val)
  if arr[i]==val:                         #Obviously should return this index
    return i
  elif i==0:                              #We're at the left edge of list, only one index to return
    return i
  elif abs(arr[i]-val)<abs(arr[i-1]-val): #Decide which is best based on difference
    return i
  else:
    return i-1



class ClimateGrid():
  def __init__(self):
    self.start_time = datetime.datetime(year=1950, month=1, day=1, hour=0, minute=0, second=0)

  def getGridByTime(self, year, month):
    time = self.yearMonthToTime(year, month)
    return self.data[time]

  def yxClosestToLatLon(self, lat, lon):
    y = binaryfind(self.lat[:],float(lat))
    x = binaryfind(self.lon[:],float(lon))
    return y, x

  def timesToUnix(self):
    return map(lambda x: time.mktime((self.start_time+datetime.timedelta(days=x)).timetuple()),self.time)

  def timesToDateTime(self):
    return map(lambda x: self.start_time+datetime.timedelta(days=x),self.time)

  def yearMonthsToTimes(self, startyear, endyear, months):
    years     = range(int(startyear),int(endyear)+1)
    yms       = itertools.product(years,months)
    yms       = map(lambda ym: self.yearMonthToTime(ym[0],ym[1]), yms)
    return yms #TODO: Verify that this is correct

  def yearRangeToTimeRange(self, startyear, endyear):
    start_time = self.yearMonthToTime(startyear,1)
    end_time   = self.yearMonthToTime(endyear,12)
    return (start_time,end_time)

  def pointMean(self, lat, lon, startyear, endyear, months):
    y,x      = self.yxClosestToLatLon(lat,lon)
    times    = self.yearMonthsToTimes(startyear, endyear, months)
    return np.mean(self.data[times,y,x])

  def yearMonthToTime(self, year, month):
    year      = int(year)
    month     = int(month)
    this_time = datetime.datetime(year=year, month=month, day=16, hour=12, minute=0) #Data is labeled in the middle of the month on the 15th/16th @ midnight
    days      = (this_time-self.start_time).days
    return binaryfind(self.time, days)

  def maxVals(self, startyear, endyear, months=None):
    if months:
      times = self.yearMonthsToTimes(startyear, endyear, months)
      tgrid = self.data[times]
    else:
      start,end = self.yearRangeToTimeRange(startyear,endyear)
      tgrid     = self.data[start:end]

    tgrid[tgrid>1e15] = np.nan
    max_grid          = np.nanmax(tgrid, axis=0)
    return max_grid

  def minVals(self, startyear, endyear, months=None):
    if months:
      times = self.yearMonthsToTimes(startyear, endyear, months)
      tgrid = self.data[times]
    else:
      start,end = self.yearRangeToTimeRange(startyear,endyear)
      tgrid     = self.data[start:end]

    tgrid[tgrid>1e15] = np.nan
    min_grid          = np.nanmin(tgrid, axis=0)
    return min_grid

  def _indexOf3(self, startyear, func):
    start_time = self.yearMonthToTime(startyear,1)
    best       = start_time*np.ones(shape=self.data.shape[1:3]) #Index of best place to start
    previous   = None #Previous grouping
    for t in range(start_time,start_time+12+1): #To make range inclusive
      current               = np.sum(self.data[start_time:start_time+3],axis=0)
      current[current>1e15] = np.nan
      if previous is None:
        previous = current
      else:
        better       = func(current,previous)
        best[better] = start_time
    return start_time

  def indexMinOf3(self, startyear):
    return self._indexOf3(startyear, lambda current,previous: current<previous)

  def indexMaxOf3(self, startyear):
    return self._indexOf3(startyear, lambda current,previous: current>previous)

  def meanVals(self, startyear, endyear, months=None):
    if months:
      times = self.yearMonthsToTimes(startyear, endyear, months)
      tgrid = self.data[times]
    else:
      start,end = self.yearRangeToTimeRange(startyear,endyear)
      tgrid     = self.data[start:end]

    tgrid[tgrid>1e15] = np.nan
    avg_grid          = np.mean(tgrid, axis=0)
    return avg_grid

  def stdVals(self, startyear, endyear, months=None):
    if months:
      times = self.yearMonthsToTimes(startyear, endyear, months)
      tgrid = self.data[times]
    else:
      start,end = self.yearRangeToTimeRange(startyear,endyear)
      tgrid     = self.data[start:end]

    tgrid[tgrid>1e15] = np.nan
    std_grid          = np.std(tgrid, axis=0)
    return std_grid

  def timeSeries(self, lat, lon):
    y = binaryfind(self.lat[:],float(lat))
    x = binaryfind(self.lon[:],float(lon))
    return self.data[:,y,x]


class NetCDFClimateGrid(ClimateGrid):
  def __init__(self, filename, varname):
    ClimateGrid.__init__(self)
    self.fin  = netcdf.netcdf_file(filename,'r')
    self.data = self.fin.variables[varname]
    self.lat  = self.fin.variables['latitude'][:]
    self.lon  = self.fin.variables['longitude'][:]
    self.time = self.fin.variables['time'][:]

  def varNames(self):
    return self.fin.variables

  def varShape(self, var):
    return self.fin.variables[var].shape

  def varUnits(self, var):
    return self.fin.variables[var].units

  def varDims(self, var):
    return self.fin.variables[var].dimensions

class HDFClimateGrid(ClimateGrid):
  def __init__(self, filename, varname):
    ClimateGrid.__init__(self)
    self.fin  = h5py.File(filename,'r')
    self.data = self.fin[varname]
    self.lat  = self.fin['latitude'][:]
    self.lon  = -(360-self.fin['longitude'][:]) #Convert from degrees East [0, 360)
    self.time = self.fin['time'][:]

  def varNames(self):
    return self.fin.keys()








#Mean Temperature of Wettest Quarter
def ModelAccum(models, func):
  accum   = None
  for m in models:
    sys.stderr.write('.')
    val = func(models[m])
    if accum is None:
      accum = val
    else:
      accum += val
  sys.stderr.write("\n")
  return accum

#Annual Mean Temperature
def AnnualMeanTemperature(models, startyear, endyear):
  return ModelAccum(models,lambda x: x['tas'].meanVals(startyear,endyear)) /len(models)

#Temperature Seasonality (variation across 12 months)
def TemperatureSeasonality(models, startyear, endyear):
  return ModelAccum(models,lambda x: x['tas'].stdVals(startyear,endyear)) /len(models)

def MaxTemp(models, startyear, endyear):
  return ModelAccum(models,lambda x: x['tasmax'].maxVals(startyear,endyear)) /len(models)

def MinTemp(models, startyear, endyear):
  return ModelAccum(models,lambda x: x['tasmin'].minVals(startyear,endyear)) /len(models)

def Maxpr(models, startyear, endyear):
  return ModelAccum(models,lambda x: x['pr'].maxVals(startyear,endyear)) /len(models)

def Minpr(models, startyear, endyear):
  return ModelAccum(models,lambda x: x['pr'].minVals(startyear,endyear)) /len(models)

#Temperature Seasonality (variation across 12 months)
def PrecipitationSeasonality(models, startyear, endyear):
  return ModelAccum(models,lambda x: x['pr'].stdVals(startyear,endyear)) /len(models)



#Mean Diurnal Range (Mean of monthly (max temp - min temp))
def MeanDiurnalRange(models, startyear, endyear):
  accum                = None
  firstmodel           = list(models.values())[0]
  start_time, end_time = firstmodel['tasmax'].yearRangeToTimeRange(startyear,endyear)
  for m in models:
    for t in range(start_time,end_time+1):
      maxdat               = models[m]['tasmax'].data[t]
      mindat               = models[m]['tasmin'].data[t]
      maxdat[maxdat>1e15]  = np.nan
      mindat[mindat>1e15] = np.nan
      val                  = maxdat-mindat
      if accum is None:
        accum = val
      else:
        accum += val
  accum /= len(models)*(start_time-end_time+1)
  return accum

def _indAccum(models, startyear, endyear, indvar, sumvar, maxmin, mean):
  accum                = None
  code.interact(local=locals())
  firstmodel           = list(models.values())[0]
  start_time, end_time = firstmodel[indvar].yearRangeToTimeRange(startyear,endyear)
  for m in models:
    sys.stderr.write('.')
    for t in range(start_time,end_time+1,12):
      if maxmin=='max':
        ind = models[m][indvar].indexMaxOf3(startyear)
      elif maxmin=='min':
        ind = models[m][indvar].indexMinOf3(startyear)
      val = np.sum(models[m][sumvar].data[ind:ind+3],axis=0)
      if mean:
        val/=3
      if accum is None:
        accum = val
      else:
        accum += val
  sys.stderr.write("\n")
  return accum

def MeanTempWettest(models, startyear, endyear):
  return _indAccum(models,startyear,endyear,'pr','tas','max',mean=True)/len(models)

def MeanTempDriest(models, startyear, endyear):
  return _indAccum(models,startyear,endyear,'pr','tas','min',mean=True)/len(models)

def MeanTempWarmest(models, startyear, endyear):
  return _indAccum(models,startyear,endyear,'tas','tas','max',mean=True)/len(models)

def MeanTempCoolest(models, startyear, endyear):
  return _indAccum(models,startyear,endyear,'tas','tas','min',mean=True)/len(models)

def AnnualPrecip(models, startyear, endyear):
  return ModelAccum(models,lambda x: x['pr'].meanVals(startyear,endyear))

def prWesttest(models, startyear, endyear):
  return _indAccum(models,startyear,endyear,'pr','pr','max',mean=False)/len(models)

def prDriest(models, startyear, endyear):
  return _indAccum(models,startyear,endyear,'pr','pr','min',mean=False)/len(models)

def prWarmest(models, startyear, endyear):
  return _indAccum(models,startyear,endyear,'tas','pr','max',mean=False)/len(models)

def prCoolest(models, startyear, endyear):
  return _indAccum(models,startyear,endyear,'tas','pr','min',mean=False)/len(models)

def writeArrayToArcGrid(filename,arr,exmodel):
  arr                = np.copy(arr)
  arr[np.isnan(arr)] = -9999
  arr                = np.flipud(arr)
  fout               = open(filename,'wb')
  headerstring       = bytes('NCOLS %d\nNROWS %d\nXLLCENTER %f\nYLLCENTER %f\nCELLSIZE %f\nNODATA_value %f\n' % (arr.shape[1], arr.shape[0], exmodel.lon.min(), exmodel.lat.min(), abs(exmodel.lat[1]-exmodel.lat[0]),-9999), 'UTF-8')
  fout.write(headerstring)
  np.savetxt(fout,arr,'%5.2f')





parser = argparse.ArgumentParser(description='Aggregates gridded climate info for use in MaxEnt')
parser.add_argument('basedir',       type=str, help='Base directory of climate data')
parser.add_argument('rcp',           type=str, help='historical/rcp26/rcp45/rcp60/rcp85')
parser.add_argument('startyear',     type=int, help='Start year of averaging interval (inclusive)')
parser.add_argument('endyear',       type=int, help='End year of averaging interval (inclusive)')
parser.add_argument('output_prefix', type=str, help='Prefix to output files')
args = parser.parse_args()

files = []
for root, dirnames, filenames in os.walk(args.basedir):
  for filename in fnmatch.filter(filenames, 'BCSD*.nc'):
    files.append(os.path.join(root, filename))

data = {}
for fname in files:
  fnameparts = fname.split('_')
  variable   = fnameparts[2]
  model      = fnameparts[4]
  rcp        = fnameparts[5]
  if not rcp in data:
    data[rcp] = {}
  if not model in data[rcp]:
    data[rcp][model] = {}
  print(rcp,model,variable)
  data[rcp][model][variable] = HDFClimateGrid(fname, variable)
  dts_of_model               = list(data[rcp][model][variable].timesToDateTime())
  print("\tStart: %s, End %s" % (dts_of_model[0 ].strftime("%Y-%m-%d"), dts_of_model[-1].strftime("%Y-%m-%d")))
  examplemodel               = data[rcp][model][variable]

varstocalculate = [AnnualMeanTemperature,TemperatureSeasonality,MaxTemp,MinTemp,Maxpr,Minpr,PrecipitationSeasonality,MeanDiurnalRange,MeanTempWettest,MeanTempDriest,MeanTempWarmest,MeanTempCoolest,AnnualPrecip,prWesttest,prDriest,prWarmest,prCoolest]
#varstocalculate = [AnnualMeanTemperature]
for v in varstocalculate:
  print("Running %s..." % (v.__name__))
  res = v(data[args.rcp],args.startyear,args.endyear)
  writeArrayToArcGrid(args.output_prefix + v.__name__+'.asc',res,examplemodel)
  np.save(v.__name__,res)
