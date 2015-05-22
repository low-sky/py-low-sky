import sys
import numpy as np
from astropy.table import Table, Column
import csv
import os.path as path
import re
filelist = sys.argv[1:]
scriptloc = sys.argv[0]
tagloc = scriptloc.replace('fluxparser.py','line_list.txt')
names = ('MaskName','ObsDate','Slit','Object','DBKey')
dtypes = ('S10','S10','S10','S10','S40')
units = (None,None,None,None,None)
with open(tagloc, 'rb') as csvfile:
        reader = csv.DictReader(csvfile)
        dix = []
        for row in reader:
            if row['line_name']:
                names = names+(row['line_name'],
                               row['line_name']+'_err')
                dtypes = dtypes + ('f4','f4')
                units = units+('W/(m**2)','W/(m**2)')
                dix.append(row)
t=Table()
for idx,thisname in enumerate(names):
    c = Column(name=thisname,unit=units[idx],dtype=dtypes[idx])
    t.add_column(c)

for thisfile in filelist:
    f = open(thisfile)
    textlines = f.readlines()
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
    linelist = Table.read(tagloc,format='csv')
    for idx,thisline in enumerate(textlines):
	    values = np.array(thisline.split(' '),dtype='float')
	    index = np.where(np.abs(linelist['wavelength']-values[0])<0.5)
	    print(thisline,index)
	    if index:
		    index = index[0]
		    t[linelist['line_name'][index][0]][-1] = values[1]
		    t[(linelist['line_name'][index][0])+'_err'][-1] = values[2]
t.write('flux_table.fits',overwrite=True)
