# make pretty Solar Disk overview plots

from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
import numpy as np
import matplotlib.pyplot as plt

import T2084Helpers.misc as misc

'''
Purpose: This file plots the relative positions of the Earth and the Moon as
they transit across the Solar Disk.
'''


def readEPH():
    '''Load and transform Event Table'''
    print('read Event Table')
    eph = Table.read('Data/ArescentricPlotTable.ecsv', format='ascii.ecsv')

    '''
    remove rows starting with Tflag - keep calculated rows only
    '''
    ijk = []
    for ii in range(len(eph)):
        if eph['Event'][ii].startswith('Tflag'):
            ijk.append(ii)
    eph.remove_rows(ijk)

    '''
    get event timings
    '''
    zTime = eph['datetime_str']
    zTime = [s.replace('2084-Nov-10 ', '') for s in zTime]
    zTime = [s.replace('.000', '') for s in zTime]

    '''Object radii in degrees'''
    SunR = Angle((eph['S_ang_width']/2))
    TerraR = Angle((eph['T_ang_width']/2))
    LunaR = Angle((eph['L_ang_width']/2))

    '''
    Objects in SkyCoords objects:
    Note to self: RA & DEC are already reported in planetary ICRF.
    Further coordinate transformations not necessary.
    '''
    Sun = SkyCoord(eph['S_RA_app'], eph['S_DEC_app'],
                   frame='icrs', unit='deg')
    Terra = SkyCoord(eph['T_RA_app'], eph['T_DEC_app'],
                     frame='icrs', unit='deg')
    Luna = SkyCoord(eph['L_RA_app'], eph['L_DEC_app'],
                    frame='icrs', unit='deg')

    Obj1frame = Sun.skyoffset_frame()
    TerraOff = Terra.transform_to(Obj1frame)
    LunaOff = Luna.transform_to(Obj1frame)

    return (eph, zTime, SunR, TerraR, LunaR, Sun, Terra, Luna,
            TerraOff, LunaOff)


def plotDisks(r, ra, dec):
    theta = np.linspace(0, 2*np.pi, 100)
    x = r*np.cos(theta) + ra
    y = r*np.sin(theta) + dec
    return x, y


def makeSolarDisk(xpdfname, xdrift=True):
    '''
    Transit across the solar disk
    Options for plots
    '''
    # update labelling options

    xtag = ('Viewed from the Martian-Geocenter.' +
            'RA and DEC in planetary ICRF frame')
    xlabtag = 'RA (degrees)'
    ylabtag = 'DEC (degrees)'

    if xdrift:
        xtitleappend = '\nSolar disk moves across the sky\n' + xtag
        xlabstr = xlabtag
        ylabstr = ylabtag
        legendloc = 'upper right'
        Centre = Sun
    else:
        xtitleappend = '\nSolar disk centered\n' + xtag
        xlabstr = r'$\Delta$'+' '+xlabtag
        ylabstr = r'$\Delta$'+' '+ylabtag
        legendloc = 'lower left'
        Centre = SkyCoord([0]*len(Sun), [0]*len(Sun),
                          frame='icrs', unit='deg')

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection=None)
    # each event gets its own colour
    colormap = iter(plt.cm.rainbow(np.linspace(0, 1, len(eph))))
    for ii in range(len(eph)):
        color = next(colormap)
        ax.scatter(Centre[ii].ra.deg, Centre[ii].dec.deg,
                   color=color, alpha=0.9, marker='+',
                   label=eph['Object'][ii] + " " +
                   eph['Event'][ii] + " " + zTime[ii])
        x, y = plotDisks(r=SunR.deg[ii],
                         ra=Centre[ii].ra.deg,
                         dec=Centre[ii].dec.deg)
        ax.plot(x, y, color=color, alpha=0.7)
        x, y = plotDisks(r=TerraR.deg[ii],
                         ra=TerraOff[ii].lon.deg + Centre[ii].ra.deg,
                         dec=TerraOff[ii].lat.deg + Centre[ii].dec.deg)
        ax.plot(x, y, color=color, alpha=0.9)
        x, y = plotDisks(r=LunaR.deg[ii],
                         ra=LunaOff[ii].lon.deg + Centre[ii].ra.deg,
                         dec=LunaOff[ii].lat.deg + Centre[ii].dec.deg)
        ax.plot(x, y, color=color, alpha=0.9)

    ax.set_aspect(1)
    ax.legend(fontsize='x-small', loc=legendloc)
    ax.invert_xaxis()
    ax.set(xlabel=xlabstr, ylabel=ylabstr)
    # plt.figtext(0.5, 0.025, xtag, wrap=True, horizontalalignment='center')
    plt.figtext(0.5, 0.01, misc.xfooter, size=6,
                horizontalalignment='center')
    plt.title('Terra & Luna Transit from Mars 2084-11-10 UTC' + xtitleappend)

    plt.savefig(xpdfname + '.pdf')
    print(xpdfname + '.pdf saved')
    plt.close()


(eph, zTime, SunR, TerraR, LunaR, Sun, Terra,
 Luna, TerraOff, LunaOff) = readEPH()
makeSolarDisk(xpdfname='Output/AresGeocentricDrift', xdrift=True)
makeSolarDisk(xpdfname='Output/AresGeocentricOffset', xdrift=False)
