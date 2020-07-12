from astroquery.jplhorizons import Horizons
from astropy.table import Table, join
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
import numpy as np

## https://docs.astropy.org/en/stable/api/astropy.table.Table.html

def dl_flats(xtag, idx, locationx, epochsx, quantitiesx):
    '''
    Standard set-up for downloading flatfiles. Minimal preprocessing.
    '''
    obj = Horizons(id = idx, location = locationx, epochs = epochsx, id_type = 'majorbody')
    eph = obj.ephemerides(quantities = quantitiesx, extra_precision = True) # Precision increases only for the RA DEC only, not the angles
    dropthese = ['targetname']
    del eph[dropthese]
    oldcol = list(eph.columns)
    oldcol = [s for s in oldcol if not s.startswith('datetime')]
    newcol = [xtag + '_' + s for s in oldcol]
    eph.rename_columns(oldcol, newcol)
    return eph

# 1(astrometric RA is RA, DEC), 2(apparent RA is RA_app, DEC_app), 4(AziEl is AZ, EL) ,
# 12(Angular Separation is sat_sep, sat_vis), 13(Target Angular Diameter is ang_width), 23 (SOT is elong, elongFlag)

def downloadSolarDiskFlat():
    '''
    Download and preprocess files for Overview Disc: note that variables are hardcoded.
    '''
    print('downloading Overview Disk files needed from Horizons')
    # epochs_use = {'start':'2084-11-10 02:00', 'stop':'2084-11-10 13:05', 'step':'39900'}
    epochs_use = {'start':'2084-11-10 02:03', 'stop':'2084-11-10 13:02', 'step':'39540'}
    ephT1 = dl_flats('T', '399','500@499', epochs_use, '2,12,13,23') # Earth from Mars geocentric
    ephS1 = dl_flats('S', '10','500@499', epochs_use, '2,13') # Sun from Mars geocentric
    ephL1 = dl_flats('L', '301','500@499', epochs_use, '2,13,23') # Moon from Mars geocentric
    eph = join(ephT1, ephS1, keys = ['datetime_str','datetime_jd'])
    eph = join(eph, ephL1, keys = ['datetime_str','datetime_jd'])
    # type(eph); print(eph.colnames)
    dropthese = ['_flags', '_solar_presence', '_RA', '_DEC']
    dropcol = [s for s in list(eph.columns) if any(map(s.endswith, dropthese))]
    del eph[dropcol]
    print('writing files')
    eph.write('ephTSL1.ecsv', format = 'ascii.ecsv', overwrite=True)
    eph.write('ephTSL1.csv', format = 'csv', overwrite=True)
    print('download complete')
    return eph

def readSolarDiskFlat():
    print('read downloaded files for Overview Disc')
    return Table.read('ephTSL1.ecsv', format = 'ascii.ecsv')

def angular_dist_arcsec( ra1, dec1, ra2, dec2):
    '''
    Takes RA and DEC in radian form and return the angular distance in ~radian~ arcsec form
    d-ang = 2 arcsin * sqrt(A + B), where
    A = sin2(abs(DEC1 - DEC2)/2) and B = cos(DEC1)cos(DEC2)sin2(abs(RA1 - RA2)/2)
    Also, 1 radian = (180/np.pi) degree = 3600 * (180/np.pi) arcsec
    '''
    A = np.sin(np.abs(dec1 - dec2)/2)**2
    B = np.cos(dec1)*np.cos(dec2)*np.sin(np.abs(ra1 - ra2)/2)**2
    AngDist = 2*np.arcsin(np.sqrt(A + B))*(180*3600)/np.pi
    return AngDist

def cal_angdist(eph, transiting_objects):
    '''
    Calculate Angular Distances between Sun and Objects, in arcsec
    '''
    ephnew = eph.copy()
    Solar_RA = np.radians(ephnew['S_RA_app'])
    Solar_DEC = np.radians(ephnew['S_DEC_app'])
    for tobject in transiting_objects:
        this_RA = tobject + '_RA_app'
        this_RA = np.radians(ephnew[this_RA])
        this_DEC = tobject + '_DEC_app'
        this_DEC = np.radians(ephnew[this_DEC])
        this_AngDist = tobject + 'S_angsep'
        ephnew[this_AngDist] = angular_dist_arcsec(Solar_RA, Solar_DEC, this_RA, this_DEC) * u.arcsec
    return ephnew

def solar_InEgress(eph, transiting_objects):
    '''
    Identifying when an object transits the solar disk
    Method
    # Outer Ingress/Egress: xS_angsep < object ang radius + solar ang radius
    # Inner Ingress/Egress: xS_angsep < object ang radius - solar and radius
    # Point of Greatest Transit:
    ## For Terra, Angular Separation at minimum as this info is available
    ## In general, where SOT is minimum
    # Photo Opportunity: where both Luna and Terra are in the frame: obviously subjective
    # Where multiple rows are returned, take the median row (or the earlier median row if number of rows are even)
    '''
    print('Identifying interesting events')
    ### Terra Direct Masking by way of T_sat_vis = t
    ephnew = eph.copy()
    mask = ephnew['T_sat_vis'] == 't'
    tA = []; tB = []; tC = []; tD = []
    tA.append('T'); tB.append('Tflag_start'); tC.append(ephnew['datetime_str'][mask][0]); tD.append(ephnew['datetime_jd'][mask][0])
    tA.append('T'); tB.append('Tflag_end'); tC.append(ephnew['datetime_str'][mask][len(ephnew[mask])-1]); tD.append(ephnew['datetime_jd'][mask][len(ephnew[mask])-1])
    ### My Calculations
    ephnew = eph.copy()
    for tobject in transiting_objects:
        this_width = tobject + '_ang_width'
        this_angsep = tobject + 'S_angsep'
        #this_SOT = tobject + '_elong'
        ### Outer Transit boundaries
        Touter = tobject + 'S_Touter'
        # ephnew[Touter] = 3600*ephnew[this_SOT] < (ephnew['S_ang_width'] + ephnew[this_width])/2
        ephnew[Touter] = ephnew[this_angsep] < (ephnew['S_ang_width'] + ephnew[this_width])/2
        nptTouter = ephnew['datetime_str','datetime_jd'][ephnew[Touter]].copy()
        a = tobject; b = 'Outer_Ingress'; c = nptTouter['datetime_str'][0]; d = nptTouter['datetime_jd'][0]
        tA.append(a); tB.append(b); tC.append(c); tD.append(d)
        a = tobject; b = 'Outer_Egress'; c = nptTouter['datetime_str'][len(nptTouter)-1]; d = nptTouter['datetime_jd'][len(nptTouter)-1]
        tA.append(a); tB.append(b); tC.append(c); tD.append(d)
        ### Inner Transit boundaries
        Tinner = tobject + 'S_Tinner'
        # ephnew[Tinner] = 3600*ephnew[this_SOT] < (ephnew['S_ang_width'] - ephnew[this_width])/2
        ephnew[Tinner] = ephnew[this_angsep] < (ephnew['S_ang_width'] - ephnew[this_width])/2
        nptTinner = ephnew['datetime_str','datetime_jd'][ephnew[Tinner]].copy()
        a = tobject; b = 'Inner_Ingress'; c = nptTinner['datetime_str'][0]; d = nptTinner['datetime_jd'][0]
        tA.append(a); tB.append(b); tC.append(c); tD.append(d)
        a = tobject; b = 'Inner_Egress'; c = nptTinner['datetime_str'][len(nptTinner)-1]; d = nptTinner['datetime_jd'][len(nptTinner)-1]
        tA.append(a); tB.append(b); tC.append(c); tD.append(d)
        ### Greatest Transit (i.e. min SOT)
        #nptmSOT = ephnew['datetime_str','datetime_jd'][ephnew[this_SOT] == np.min(ephnew[this_SOT])].copy()
        nptmSOT = ephnew['datetime_str','datetime_jd'][ephnew[this_angsep] == np.min(ephnew[this_angsep])].copy()
        i = max(int(np.floor(len(nptmSOT)/2))-1,0)
        a = tobject; b = 'Greatest_Transit'; c = nptmSOT['datetime_str'][i]; d = nptmSOT['datetime_jd'][i]
        tA.append(a); tB.append(b); tC.append(c); tD.append(d)
    ### Photo Opportunity
    ephnew = eph.copy()
    ephnew['TL_adist'] = ephnew['TS_angsep'] + ephnew['LS_angsep']
    nptPhoto = ephnew['datetime_str','datetime_jd'][ephnew['TL_adist'] == np.min(ephnew['TL_adist'])].copy()
    i = max(int(np.floor(len(nptPhoto)/2))-1,0)
    a = 'TL'; b = 'Photo_Op'; c = nptPhoto['datetime_str'][i]; d = nptPhoto['datetime_jd'][i]
    tA.append(a); tB.append(b); tC.append(c); tD.append(d)
    ### Sort output table
    Tout = Table([tA, tB, tC, tD], names=('Object','Event','datetime_str','datetime_jd'))
    Tout.sort('datetime_jd')
    returnthese = ['datetime_jd','datetime_str','_ang_width','_angsep','_RA_app','_DEC_app']
    returncol = [s for s in list(ephnew.columns) if any(map(s.endswith, returnthese))]
    plottable = join(Tout, ephnew[returncol])
    print('Writing Event Table')
    plottable.write('GeocentricPlotTable.ecsv', format = 'ascii.ecsv', overwrite=True)
    # return Tout, plottable, ephnew[returncol]
    return Tout['Object','Event','datetime_str'], plottable

#############

### Download or Read from earlier downloads
# eph = downloadSolarDiskFlat()
eph = readSolarDiskFlat()

### Calculate angular distances between sun and transiting objects
# eph2 = cal_angdist(eph, ['T','L']) # your own equations
# use build in from SkyCoord instead
eph2 = eph.copy()
Solar = SkyCoord(eph['S_RA_app'], eph['S_DEC_app'], frame = 'icrs', unit = 'deg')
Terra = SkyCoord(eph['T_RA_app'], eph['T_DEC_app'], frame = 'icrs', unit = 'deg')
Luna = SkyCoord(eph['L_RA_app'], eph['L_DEC_app'], frame = 'icrs', unit = 'deg')
eph2['TS_angsep'] = Solar.separation(Terra).arcsec * u.arcsec
eph2['LS_angsep'] = Solar.separation(Luna).arcsec * u.arcsec

### Generate Summary Events
SummaryDiskEvents, plottable = solar_InEgress(eph2, ['T', 'L'])
SummaryDiskEvents.pprint_all()
plottable.pprint_all()

# eyeball values: print(eph[eph['datetime_jd'] == 2482540.039189815])
