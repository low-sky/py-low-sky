import sys
import numpy as np
from astropy.table import Table, Column
import csv
import os.path as path
import re
filelist = sys.argv[1:]
scriptloc = sys.argv[0]
tagloc = scriptloc.replace('neatparser.py','neat_tags.txt')
names = ('MaskName','ObsDate','Slit','Object','DBKey')
dtypes = ('S10','S10','S10','S10','S40')
units = (None,None,None,None,None)
with open(tagloc, 'rb') as csvfile:
        reader = csv.DictReader(csvfile)
        dix = []
        for row in reader:
            if row['column_name']:
                names = names+(row['column_name'],
                               row['column_name']+'+err',
                               row['column_name']+'-err',
                               row['column_name']+'_flag')
                dtypes = dtypes + ('f4','f4','f4','S8')
                units = units+(row['units'],row['units'],row['units'],None)
                dix.append(row)
t=Table()
for idx,thisname in enumerate(names):
    c = Column(name=thisname,unit=units[idx],dtype=dtypes[idx])
    t.add_column(c)
    
for thisfile in filelist:
    f = open(thisfile)
    lines = f.readlines()
    f.close()
    t.add_row()
    method='NONE'
    splitpath = path.abspath(thisfile).split('/')
    splitname = re.split(r'(_|\.)',thisfile)[::2]
    t['MaskName'][-1] = splitpath[-2]
    t['ObsDate'][-1] = splitpath[-3]
    t['Slit'][-1] = splitname[1]
    t['Object'][-1] = splitname[2]
    t['DBKey'][-1] = t['MaskName'][-1]+'-'+\
	t['ObsDate'][-1]+'-'+\
	t['Slit'][-1]+'-'+\
	t['Object'][-1]

    for idx,thisline in enumerate(lines):
        if 'ORL' in thisline:
            method = 'ORL'
        if 'Strong line' in thisline:
            method = 'NONE'
        if 'CEL' in thisline:
            method = 'CEL'
        substrs = thisline.split('  ')
        substrs = filter(lambda a: a != '',substrs)
        substrs = [s.strip() for s in substrs]
        for row in dix:
            if (row['neat_tag'] == substrs[0]) and ((method in row['column_name']) or (method == 'NONE')):
#                try:
                    if len(substrs)>3:
                        start = len(substrs)-3
			try:
				values = np.array(substrs[start:],dtype='float')
				t[row['column_name']][-1] = values[0]
				t[row['column_name']+'+err'][-1] = values[1]
				t[row['column_name']+'-err'][-1] = values[2]
			except ValueError:
				t[row['column_name']][-1] = 0.0
				t[row['column_name']+'+err'][-1] = 0.0
				t[row['column_name']+'-err'][-1] = 0.0
                        if 'Warning' in lines[idx-1]:
                            t[row['column_name']+'_flag'][-1] = 'W'
t.write('neat_table.fits',overwrite=True)
