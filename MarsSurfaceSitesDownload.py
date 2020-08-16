# This module downloads data to analyse the transit visibility from various
# interesting sites on the Martian Surface

from astroquery.jplhorizons import Horizons
from astropy.table import Column, Table, join, vstack
from astropy import units as u
from astropy.coordinates import SkyCoord
# from astropy.coordinates import Angle
import numpy as np


def dl_flats(xtag, idx, locationx, epochsx, quantitiesx):
    '''
    Standard set-up for downloading flatfiles. Minimal preprocessing.
    '''
    obj = Horizons(id=idx, location=locationx, epochs=epochsx,
                   id_type='majorbody')
    # Precision increases only for the RA DEC only, not the angles
    eph = obj.ephemerides(quantities=quantitiesx, extra_precision=True)
    dropthese = ['targetname']
    del eph[dropthese]
    oldcol = list(eph.columns)
    oldcol = [s for s in oldcol if not s.startswith('datetime')]
    newcol = [xtag + '_' + s for s in oldcol]
    eph.rename_columns(oldcol, newcol)
    return eph


# 4(apparent Aximuth and Elevation is AZ, EL),
# 12(Angular Separation is sat_sep, sat_vis) of target and target parent
# 13(Target Angular Diameter is ang_width),
# 23(Sun-Observer-Target SOT elongation angle is elong, elongFlag)
# 25(Target-Observer+Moon TOM angle)
# body code: Sun = 10, Terra = 399, Luna = 301, Phobos = 401


def downloadGroundView(sitename, locationx):
    '''
    Download and save raw files for a given location on the Martial surface.
    Angular separations are calculated here.
    Files are saved in /Data/G/
    We don't quite need items in seconds: minutes are good for now to get a
    gauge of where Phobos and Deimos are likely to be.
    (We can find tune later)
    '''
    print('downloading Groundview for '+sitename)
    epochs_use = {'start': '2084-11-10 02:00', 'stop': '2084-11-10 13:10',
                  'step': '1m'}
    # epochs_use = {'start': '2084-11-10 02:03', 'stop': '2084-11-10 13:02',
    #               'step': '39540'}
    # epochs_use = {'start':'2084-11-10 02:00', 'stop':'2084-11-10 13:05',
    #                'step':'39900'}
    # epochs_use = {'start': '2084-11-10 02:00', 'stop': '2084-11-10 13:10',
    #              'step': '40200'}
    # Earth from Mars surface
    ephT1 = dl_flats('T', '399', locationx, epochs_use, '4,12,13,23,25')
    # Sun from Mars surface
    ephS1 = dl_flats('S', '10', locationx, epochs_use, '4,13,25')
    # Moon from Mars surface
    ephL1 = dl_flats('L', '301', locationx, epochs_use, '4,13,23,25')
    # Phobos from Mars surface
    ephP1 = dl_flats('P', '401', locationx, epochs_use, '4,13,23')
    # Deimos from Mars surface
    ephD1 = dl_flats('D', '402', locationx, epochs_use, '4,13,23')
    # join tables
    eph = join(ephS1, ephT1, keys=['datetime_str', 'datetime_jd'])
    eph = join(eph, ephL1, keys=['datetime_str', 'datetime_jd'])
    eph = join(eph, ephP1, keys=['datetime_str', 'datetime_jd'])
    eph = join(eph, ephD1, keys=['datetime_str', 'datetime_jd'])
    # type(eph); print(eph.colnames)
    dropthese = ['_illum']
    dropcol = [s for s in list(eph.columns) if any(map(s.endswith, dropthese))]
    del eph[dropcol]
    eph.add_column(Column(sitename, name='site'), index=0)

    '''
    Calculate angular separation
    '''
    Sol = SkyCoord(eph['S_AZ'], eph['S_EL'], frame='altaz', unit='deg')
    Terra = SkyCoord(eph['T_AZ'], eph['T_EL'], frame='altaz', unit='deg')
    Luna = SkyCoord(eph['L_AZ'], eph['L_EL'], frame='altaz', unit='deg')
    Phobos = SkyCoord(eph['P_AZ'], eph['P_EL'], frame='altaz', unit='deg')
    Deimos = SkyCoord(eph['D_AZ'], eph['D_EL'], frame='altaz', unit='deg')
    eph['TS_angsep'] = Sol.separation(Terra).arcsec * u.arcsec
    eph['LS_angsep'] = Sol.separation(Luna).arcsec * u.arcsec
    eph['PS_angsep'] = Sol.separation(Phobos).arcsec * u.arcsec
    eph['DS_angsep'] = Sol.separation(Deimos).arcsec * u.arcsec

    print('writing raw files')
    eph.write('Data/G/'+sitename+'.ecsv', format='ascii.ecsv', overwrite=True)
    # eph.write('Data/G/'+sitename+'.csv', format='csv', overwrite=True)
    print('download for '+sitename+' complete')
    return eph


def makeEventTable(sitename, eph, transiting_objects):
    '''
    Identifying when an object transits the solar disk
    Method
    # Outer Ingress/Egress: xS_angsep < object ang radius + solar ang radius
    # Inner Ingress/Egress: xS_angsep < object ang radius - solar and radius
    # Point of Greatest Transit:
    ## For Terra, Angular Separation at minimum as this info is available
    ## In general, where SOT is minimum
    # Photo Opportunity: where both Luna and Terra are in frame - subjective
    # Where multiple rows are returned, take the median row
    (or the earlier median row if number of rows are even)
    '''
    print('Identifying interesting events')
    ephN = eph.copy()
    tA = []
    tB = []
    tC = []
    tD = []
    # sunrise by way of S_solar_presence = *
    mask = ephN['S_solar_presence'] == '*'
    if any(mask):
        tA.append('S')
        tB.append('sunrise')
        tC.append(ephN['datetime_str'][mask][0])
        tD.append(ephN['datetime_jd'][mask][0])
        tA.append('S')
        tB.append('sunset')
        tC.append(ephN['datetime_str'][mask][len(ephN[mask])-1])
        tD.append(ephN['datetime_jd'][mask][len(ephN[mask])-1])
        ephN = eph.copy()
    mask = ephN['PS_angsep'] == np.min(ephN['PS_angsep'])
    tA.append('P')
    tB.append('Min_Angsep')
    tC.append(ephN['datetime_str'][mask][0])
    tD.append(ephN['datetime_jd'][mask][0])
    mask = ephN['DS_angsep'] == np.min(ephN['DS_angsep'])
    tA.append('D')
    tB.append('Min_Angsep')
    tC.append(ephN['datetime_str'][mask][0])
    tD.append(ephN['datetime_jd'][mask][0])
    # My Calculations
    ephN = eph.copy()
    for tobject in transiting_objects:
        this_width = tobject + '_ang_width'
        this_angsep = tobject + 'S_angsep'
        # Outer Transit boundaries
        TOut = tobject + 'S_TOut'
        ephN[TOut] = ephN[this_angsep] < ((ephN['S_ang_width'] +
                                           ephN[this_width])/2)
        nptTOut = ephN['datetime_str', 'datetime_jd'][ephN[TOut]].copy()
        tA.append(tobject)
        tB.append('Outer_Ingress')
        tC.append(nptTOut['datetime_str'][0])
        tD.append(nptTOut['datetime_jd'][0])
        tA.append(tobject)
        tB.append('Outer_Egress')
        tC.append(nptTOut['datetime_str'][len(nptTOut)-1])
        tD.append(nptTOut['datetime_jd'][len(nptTOut)-1])
        # Inner Transit boundaries
        TIn = tobject + 'S_TIn'
        ephN[TIn] = ephN[this_angsep] < ((ephN['S_ang_width'] -
                                          ephN[this_width])/2)
        nptTIn = ephN['datetime_str', 'datetime_jd'][ephN[TIn]].copy()
        tA.append(tobject)
        tB.append('Inner_Ingress')
        tC.append(nptTIn['datetime_str'][0])
        tD.append(nptTIn['datetime_jd'][0])
        tA.append(tobject)
        tB.append('Inner_Egress')
        tC.append(nptTIn['datetime_str'][len(nptTIn)-1])
        tD.append(nptTIn['datetime_jd'][len(nptTIn)-1])
        # Greatest Transit (i.e. min SOT)
        mask = ephN[this_angsep] == np.min(ephN[this_angsep])
        nptmSOT = ephN['datetime_str', 'datetime_jd'][mask].copy()
        i = max(int(np.floor(len(nptmSOT)/2))-1, 0)
        tA.append(tobject)
        tB.append('Greatest_Transit')
        tC.append(nptmSOT['datetime_str'][i])
        tD.append(nptmSOT['datetime_jd'][i])
    # Photo Opportunity
    ephN = eph.copy()
    ephN['TL_adist'] = ephN['TS_angsep'] + ephN['LS_angsep']
    mask = ephN['TL_adist'] == np.min(ephN['TL_adist'])
    nptPhoto = ephN['datetime_str', 'datetime_jd'][mask].copy()
    i = max(int(np.floor(len(nptPhoto)/2))-1, 0)
    tA.append('TL')
    tB.append('Photo_Op')
    tC.append(nptPhoto['datetime_str'][i])
    tD.append(nptPhoto['datetime_jd'][i])
    # Sort output table
    Tout = Table([tA, tB, tC, tD],
                 names=('Object', 'Event', 'datetime_str', 'datetime_jd'))
    Tout.sort('datetime_jd')
    xcol = ['datetime_jd', 'datetime_str', 'S_solar_presence', '_ang_width',
            '_angsep', '_AZ', '_EL']
    returncol = [s for s in list(ephN.columns) if any(map(s.endswith, xcol))]
    plottable = join(Tout, ephN[returncol])
    plottable.add_column(Column(sitename, name='site'), index=0)
    print('Writing Event Table\n')
    plottable.write('Data/G/'+sitename+'ETable.ecsv',
                    format='ascii.ecsv', overwrite=True)
    # plottable.write('Data/G/'+sitename+'ETable.csv',
    #                 format='csv', overwrite=True)
    return plottable


def dlAll(sites):
    '''
    For each specified location, download flat files.
    Makes interesting event table for each location.
    '''
    BigTable = Table(names=('site', 'Object', 'Event',
                            'datetime_str', 'datetime_jd', 'S_solar_presence',
                            'S_AZ', 'S_EL', 'S_ang_width',
                            'T_AZ', 'T_EL', 'T_ang_width',
                            'L_AZ', 'L_EL', 'L_ang_width',
                            'P_AZ', 'P_EL', 'P_ang_width',
                            'D_AZ', 'D_EL', 'D_ang_width',
                            'TS_angsep', 'LS_angsep', 'PS_angsep',
                            'DS_angsep'),
                     dtype=('S4', 'S4', 'S4', 'S4', 'f4', 'S4',
                            'f4', 'f4', 'f4', 'f4', 'f4', 'f4',
                            'f4', 'f4', 'f4', 'f4', 'f4', 'f4',
                            'f4', 'f4', 'f4',
                            'f4', 'f4', 'f4', 'f4'))

    for Site, HL, lon, lat, ele in sites.iterrows('Site', 'HorizonsLocation',
                                                  'Elondeg', 'Latdeg',
                                                  'Altkm'):
        print(Site, HL, lon, lat, ele)
        if HL.endswith('499'):
            loc = HL
        else:
            print('making dictionary')
            loc = {'body': '499', 'lon': lon, 'lat': lat, 'elevation': ele}
            print(loc)
        eph = downloadGroundView(Site, loc)
        PT = makeEventTable(Site, eph, ['T', 'L'])
        # print(PT)
        BigTable = vstack([BigTable, PT])
    BigTable.write('Data/AresSurfaceEventTable.ecsv',
                   format='ascii.ecsv', overwrite=True)
    BigTable.write('Data/AresSurfaceEventTable.csv',
                   format='csv', overwrite=True)
    print('Surface Event Table compiled!')
    return BigTable


'''
parse locations to download
'''
sites = Table.read('MarsSurfaceSites.csv', format='csv')
# sites = sites[[0, 46]]  # used if testing code
# fill the masked i.e. missing value with any value not ending with 499
sites['HorizonsLocation'].fill_value = 'not_provided'
sites = sites.filled()

checkthis = dlAll(sites)
print(checkthis)
