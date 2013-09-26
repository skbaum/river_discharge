#!/usr/bin/python2.7

#  Program:  discharge_year.py
#
#  This program creates a NetCDF file from a combination of USGS and USACE river discharge
#  values for a given year.  The USGS values are automatically downloaded from a USGS
#  server using the station numbers in the 'stationid' array (except for the first two).
#  The USACE river discharge values are read from files of the form usace_xxxx_yyyy.txt
#  where xxxx is either 'miss' or 'atch' and yyyy is the year.  These files are obtained
#  via the interactive forms at:
#
#   miss- http://www2.mvn.usace.army.mil/cgi-bin/wcmanual.pl?01100
#   atch - http://www2.mvn.usace.army.mil/cgi-bin/watercontrol.pl?03045
#
#  A year is specified on the form, the submit button is pressed, and the resulting
#  text file is downloaded.  The file is then hand edited to remove the extraneous
#  lines from the top and bottom.
#
#  More details about the procedure can be found at:
#
#   http://stommel.tamu.edu/~baum/river.html
#
#   

import urllib,re
import numpy as np
import netCDF4
from netCDF4 import num2date,date2num,datetime
import time,datetime
from time import strftime
from datetime import timedelta,date
import itertools
import sys

year = "2011"
dstart = year + "-01-01"
dend = year + "-12-31"

no_stations = 55
cfs2cms = 0.028316847
fillvalue = -999.

miss_file = "usace_miss_" + year + ".txt"
atch_file = "usace_atch_" + year + ".txt"

try:
	missfile = open(miss_file,"r")
except:
	print " Mississippi discharge file " + miss_file + " not available.  Exiting program."
	sys.exit()

try:
	atchfile = open(atch_file,"r")
except:
	print " Atchafalaya discharge file " + atch_file + " not available.  Exiting program."
	sys.exit()

miss_lines = missfile.readlines()
atch_lines = atchfile.readlines()
len_miss = len(miss_lines)
len_atch = len(atch_lines)

miss = np.zeros(len_miss)
atch = np.zeros(len_atch)

n = 0
for line in miss_lines:
	stuff = line.split(";")
	dis = float((stuff[1]).strip(' \n'))
	miss[n] = dis*1000.*cfs2cms
#	print " n = ",n," miss = ",miss[n]
	n = n + 1

n = 0
for line in atch_lines:
        stuff = line.split(";")
        dis = float((stuff[1]).strip(' \n'))
        atch[n] = dis*1000.*cfs2cms
#	print " n = ",n," miss = ",atch[n]
	n = n + 1

out = "gom_discharge_" + year + ".nc"

stationid = ['01100mis','03045atc','02292010','02296750','02298880','02310525','02310747','02313100',
             '02320500','02324000','02324400','02325000','02326000','02330000','02330150','02359000',
             '02359170','02366500','02368000','02369000','02369600','02370000','02375500','02376033',
             '02376500','02428400','02469761','02479000','02479160','02479300','02479310','02481510',
             '02489500','02492000','07375500','07378500','07385500','07386980','08012000','08015500',
             '08030500','08041000','08041500','08066500','08068090','08075000','08116650','08117500',
             '08162500','08164000','08176500','08188500','08189500','08189700','08211000']

#stationid = ['01100mis','03045atc','02292010','02296750','02298880','02310525','02310747','02313100']

units = "days since 1900-01-01 00:00:00"

#  Establish time information.

yyyy = strftime("%Y")
mm = strftime("%m")
dd = strftime("%d")

dtime = datetime.datetime(int(yyyy),int(mm),int(dd),12,0,0)
dm = datetime.timedelta(days=1)
today_time = date2num(dtime,units)
yesterday_time = date2num(dtime - dm,units)

yyyymmdd = (dtime - dm).strftime("%Y") + (dtime - dm).strftime("%m") + (dtime - dm).strftime("%d")

#  Open input file and read data.

infile = "/home/baum/RIVER/river_programs/gom_river_discharge.nc"
inf = netCDF4.Dataset(infile,'r')

End_Date = inf.variables['End_Date'][:]
Start_Date = inf.variables['Start_Date'][:]
Gage_Latitude = inf.variables['Gage_Latitude'][:]
Gage_Longitude = inf.variables['Gage_Longitude'][:]
Mouth_Latitude = inf.variables['Mouth_Latitude'][:]
Mouth_Longitude = inf.variables['Mouth_Longitude'][:]
Station_ID = inf.variables['Station_ID'][:]
River_Name = inf.variables['River_Name'][:]

inf.close()

#  Open output file.

outfile = "/home/baum/RIVER/river_programs/" + out
ts = netCDF4.Dataset(outfile,'w',format="NETCDF3_CLASSIC")

ts.createDimension('time',None)
ts.createDimension('station',no_stations)
ts.createDimension('id_strlen',8)
ts.createDimension('river_name_strlen',24)

discharge = ts.createVariable('discharge','f4',('time','station',),fill_value="-999.")
time = ts.createVariable('time','i4',('time',))

time[:] = yesterday_time
#time[:] = date2num(dtime,units)

station_id = ts.createVariable('station_id','S1',('station','id_strlen',))
river_name = ts.createVariable('river_name','S1',('station','river_name_strlen',))
end_date = ts.createVariable('end_date','f4',('station',))
start_date = ts.createVariable('start_date','f4',('station',))
gauge_latitude = ts.createVariable('gauge_latitude','f4',('station',))
gauge_longitude = ts.createVariable('gauge_longitude','f4',('station',))
mouth_latitude = ts.createVariable('mouth_latitude','f4',('station',))
mouth_longitude = ts.createVariable('mouth_longitude','f4',('station',))

ts.title = "Time-series discharge data for rivers emptying in the Gulf of Mexico."
ts.source = "USGS and U.S. Army Corps of Engineers"
ts.source_url1 = "http://waterdata.usgs.gov/nwis"
ts.source_url2 = "http://www2.mvn.usace.army.mil/eng/edhd/dailystagedisplay.asp"
ts.featureType = "timeSeries"
ts.rivers = "U.S. Gulf of Mexico Rivers"
ts.contact = "mkhoward@tamu.edu"
ts.conventions = "CF-1.6"
ts.processing_information_url = "http://stommel.tamu.edu/~baum/river.html"

time.standard_name = "time"
time.long_name = "time of measurement"
time.units = "days since 1900-01-01 00:00:00"

discharge.units = "cubic meters per second"
discharge.long_name = "discharge"
discharge.standard_name = "discharge"
discharge.valid_min = "0.0"
discharge.valid_max = "43041.5"
discharge.missing_value = -999.
discharge.coordinates = "gauge_latitude gauge_longitude"

end_date.units = "days since 1900-01-01 00:00:00"
start_date.units = "days since 1900-01-01 00:00:00"
gauge_latitude.units = "degrees_north"
gauge_latitude.long_name = "latitude of gauge station"
gauge_longitude.units = "degrees_east"
gauge_longitude.long_name = "longitude of gauge station"
mouth_latitude.units = "degrees_north"
mouth_latitude.long_name = "latitute of river mouth"
mouth_longitude.units = "degrees_east"
mouth_longitude.long_name = "longitude of river mouth"
station_id.long_name = "gauging station ID number"
station_id.cf_role = "timeseries_id"
river_name.long_name = "name of river"
start_date.long_name = "date of first gauge station data"
end_date.long_name = "date of most recent gauge station data"

end_date[:] = End_Date
start_date[:] = Start_Date
station_id[:] = Station_ID
river_name[:] = River_Name
gauge_latitude[:] = Gage_Latitude
gauge_longitude[:] = Gage_Longitude
mouth_latitude[:] = Mouth_Latitude
mouth_longitude[:] = Mouth_Longitude

timeslice = "&begin_date=" + dstart + "&end_date=" + dend

print " timeslice = ",timeslice

n = 0
for sno in stationid:

	print " Processing station " + sno

       	urlstring = "http://waterdata.usgs.gov/nwis/dv?cb_00060=on&format=rdb" + timeslice + "&site_no=" + sno

       	usgs = urllib.urlretrieve(urlstring,"tmp.txt")

#  The USGS files have a variable number of boilerplate lines at the beginning, so
#  the number must be found each time.

	nskip = 0
	with open('tmp.txt') as f:
		for line in itertools.islice(f,0,None):
			if not line.startswith("USGS"):
				nskip = nskip + 1
			

	with open('tmp.txt') as f:
		nt = 0
		for line in itertools.islice(f,nskip,None):
#			print " n = ",n," nt = ",nt," line = ",line
			try:
				stuff = line.split("\t")
				cfs = float(stuff[3])*cfs2cms
				ds = stuff[2].split("-")
				dtime = datetime.datetime(int(ds[0]),int(ds[1]),int(ds[2]))
				today_time = date2num(dtime,units)
#				print " discharge = ",float(stuff[3])," today = ",today_time

				discharge[nt,n] = cfs
				time[nt] = today_time
				nt = nt + 1

#  If discharge value unavailable, put in a fill value.
			except:

				print " Couldn't obtain desired discharge information for station ",stationid[n]," at time ",dtime
				cfs = fillvalue
				discharge[nt,n] = cfs
				nt = nt + 1


	n = n + 1

#  Add Mississippi and Atchafalaya discharge values.

for nt in range(0,len_miss):
	discharge[nt,0] = miss[nt]
	discharge[nt,1] = atch[nt]


ts.close()

