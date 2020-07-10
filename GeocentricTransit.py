from astroquery.jplhorizons import Horizons
from astropy.table import Table, join
import numpy as np

## https://docs.astropy.org/en/stable/api/astropy.table.Table.html

def dl_flats(xtag, idx, locationx, epochsx, quantitiesx):
	obj = Horizons(id = idx, location = locationx, epochs = epochsx, id_type = 'majorbody')
	eph = obj.ephemerides(quantities = quantitiesx, extra_precision = False) # Precision increases only for the RA DEC only, not the angles
	del eph['targetname']
	oldcol = list(eph.columns)
	oldcol = [s for s in oldcol if not s.startswith('datetime')]
	newcol = [xtag + '_' + s for s in oldcol]
	eph.rename_columns(oldcol, newcol)
	return eph

# 1(astrometric RA is RA, DEC), 2(apparent RA is RA_app, DEC_app), 4(AziEl is AZ, EL) ,
# 12(Angular Separation is sat_sep, sat_vis), 13(Target Angular Diameter is ang_width), 23 (SOT is elong, elongFlag)

# epochs_use = {'start':'2084-11-10 02:00', 'stop':'2084-11-10 02:01', 'step':'60'}
epochs_use = {'start':'2084-11-10 02:00', 'stop':'2084-11-10 13:05', 'step':'39900'}
ephT1 = dl_flats('T', '399','500@499', epochs_use, '1,2,12,13,23')
ephS1 = dl_flats('S', '10','500@499', epochs_use, '1,2,13')
ephL1 = dl_flats('L', '301','500@499', epochs_use, '1,2,13,23')

eph = join(ephT1, ephS1, keys = ['datetime_str','datetime_jd'])
eph = join(eph, ephL1, keys = ['datetime_str','datetime_jd'])

# type(eph)
print(eph)
print(eph.colnames)

eph.write('ephTSL1.ecsv', format = 'ascii.ecsv')
eph.write('ephTSL1.csv', format = 'csv')

# eph2 = Table.read('ephTSL1.ecsv', format = 'ascii.ecsv')

