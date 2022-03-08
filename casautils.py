# casa utilities, some adapted from Todd Hunter's analysisUtils

import numpy as np
from astroquery.simbad import Simbad
import shutil
import os

SECONDS_PER_YEAR = 31556925.
ARCSEC_PER_RAD=206264.80624709636

"""
correctForPM

Adjust the coordinates of the ALMA data to replicate what the image would look like if obsreved at a 
different time based on the proper motion of the source.

*required inputs*
inms - The original .ms directory.
outms - The name of the new .ms directory for the updated ms.
field - The field or fields to update. This can be given as a source name or field IDs.

*optional inputs*
years - Change the coordinates to what they would be if taken this many years in the future
		(or past if negative).
newtime - Change the coordinates to what they would be if taken at this time.
			Must be given as a Modified Julian Date in seconds.
updatepm - If True, get the latest proper motion from Simbad and update the .ms.
			NB: Only works if the source name is given for field and if this name is recognised by Simbad.
			
*example usage*
import casautils as cu

# Cycle 3 dataset
firstms = 'cycle3_calibrated.ms'
# Cycle 6 dataset
secondms = 'cycle6_calibrated.ms'
# Cycle 3 dataset shifted to source position at time of cycle 6 observations
shiftedms = 'cycle3to6_calibrated.ms'

# Extract observation date from cycle 6 data 
tb.open(secondms, nomodify=True)
time6 = tb.getcol('TIME').mean()	# MJD of this ms in seconds
tb.close()

cu.correctForPM(firstms,shiftedms,'hr8799',newtime=time6,updatepm=True)

*improvements wishlist*
Allow newtime to be given in a variety of formats.
Check if user actually wants to overwrite a pre-existing output file.
Instead of inputting a new time, maybe allow user to input an ms and select which field they want to use 
the time from.
"""

def correctForPM(inms,outms,field,years=0,newtime=None,updatepm=False):
	# If output exists, overwrite with data from input
	if os.path.exists(outms): shutil.rmtree(outms)
	shutil.copytree(inms, outms)
	if updatepm == True: updatePM(outms,field)
	# Check if field is a name or ID
	tb.open('%s/FIELD'%(outms),nomodify=True)
	if (isinstance(field,str)):
		names = tb.getcol('NAME')
		fieldid = np.where(names == field)[0]
	elif (isinstance(field,int)):
		fieldid = np.array([field])
	else:
		fieldid = np.array(field)
	sources = tb.getcol('SOURCE_ID')[fieldid]
	times = tb.getcol('TIME')[fieldid]
	tb.close()
	# Get the proper motion from the source table using the first instance of the source
	tb.open('%s/SOURCE'%(outms),nomodify=True)
	sourceids = tb.getcol('SOURCE_ID')
	firstinst = np.where(sourceids==sources[0])[0][0]
	properMotion = tb.getcol('PROPER_MOTION')[:,firstinst]
	tb.close()
	# Determine the difference in seconds between original observation and new time
	seconds=np.zeros_like(fieldid)
	if isinstance(newtime, float):
		for i in range(len(seconds)):
			timestamp = times[i]
			yrs = (newtime - timestamp)/SECONDS_PER_YEAR
			seconds[i] = (newtime - timestamp)
			#print("%f years (%.1f seconds) of proper motion since J2000." % (yrs,seconds[i]))
	elif (years == 0):
		for i in range(len(seconds)):
			timestamp = times[i]
			j2000 = 4453401600
			yrs = (timestamp - j2000)/SECONDS_PER_YEAR
			seconds[i] = (timestamp - j2000)
			#print("%f years (%.1f seconds) of proper motion since J2000." % (yrs,seconds[i]))
	else:
		seconds[:] = years*SECONDS_PER_YEAR
	foundSource = 0
	# Update the coordinates for each field.
	tb.open('%s/FIELD'%(outms),nomodify=False)
	for i,j in enumerate(fieldid):
		print("Calculating new direction for field %i" % j)
		for direct in ['REFERENCE_DIR','PHASE_DIR','DELAY_DIR']:
			direction = tb.getcell(direct,j)
			new_dec = direction[1] + properMotion[1]*seconds[i]
			new_ra = direction[0] + properMotion[0]*seconds[i]/np.cos(new_dec)
			new_direction = np.array([new_ra, new_dec])
			tb.putcell(direct,j,new_direction)
		foundSource = 1
	tb.close()
	if (foundSource == 0):
		print("Did not find source = ", field)
	return

def updatePM(inms,source):
	"""
	Search Simbad for the given name, download latest pm measurement and update in ms
	-Mark Booth (May 2021)
	"""
	Simbad.reset_votable_fields()
	Simbad.add_votable_fields('pmra','pmdec')
	result_table = Simbad.query_object(source)
	if result_table is None:
		print("%s not found" %(source))
	else:
		print("Proper motion found in Simbad: pmra = %.3fmas/yr, pmdec = %.3fmas/yr"
			%(result_table['PMRA'][0],result_table['PMDEC'][0]))
		# Convert to rad/s
		newpm=np.array([result_table['PMRA'][0],result_table['PMDEC'][0]])/(1000*ARCSEC_PER_RAD*SECONDS_PER_YEAR)
		tb.open('%s/SOURCE'%(inms),nomodify=False)
		names=tb.getcol('NAME')
		sid=np.where(names==source)
		pm=tb.getcol('PROPER_MOTION')[:,sid[0][0]]*(1000*ARCSEC_PER_RAD*SECONDS_PER_YEAR)
		print("Current proper motion in ms: pmra = %.3fmas/yr, pmdec = %.3fmas/yr"
			%(pm[0],pm[1]))
		for s in sid:
			tb.putcell('PROPER_MOTION',s,newpm)
		print("Proper motion updated")
		tb.close()
