# make pretty Solar Disk overview plots

from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u

import numpy as np
import matplotlib.pyplot as plt


'''To Do: Coordinate transformation to Mars rotational frame.'''


def readEPH():
    '''Load and transform Event Table'''
    print('read Event Table')
    eph = Table.read('Data/GeocentricPlotTable.ecsv', format='ascii.ecsv')

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

    '''Mars rotational axis & J2000 rotational axis'''
    Tref = 2451545  # reference JD for '2000-01-01T12:00:00'
    # print(Tref.jd)
    MarsRRA = (317.681 - 0.106*(eph['datetime_jd']-Tref)/36525)
    MarsRDEC = (52.887 - 0.061*(eph['datetime_jd']-Tref)/36525)
    MarsRRA.unit = u.deg
    MarsRDEC.unit = u.deg
    MarsR = SkyCoord(MarsRRA, MarsRDEC, frame='icrs', unit='deg')
    # J2000R = SkyCoord([0]*len(eph), [90]*len(eph), frame='icrs', unit='deg')
    # RotSept = J2000R.separation(MarsR).deg * u.deg
    # RotPosAng = J2000R.position_angle(MarsR).deg * u.deg
    # eph['RotSept'] = RotSept
    # eph['RotPosAng'] = RotPosAng

    '''
    Objects in SkyCoords objects:
    Object1: J2000 frame
    Object2: Mars rotational frame (WIP)
    '''
    Sun1 = SkyCoord(eph['S_RA_app'], eph['S_DEC_app'],
                    frame='icrs', unit='deg')
    Terra1 = SkyCoord(eph['T_RA_app'], eph['T_DEC_app'],
                      frame='icrs', unit='deg')
    Luna1 = SkyCoord(eph['L_RA_app'], eph['L_DEC_app'],
                     frame='icrs', unit='deg')
    # SunJ200RPosAng = Sun1.position_angle(J2000R).deg * u.deg
    # eph['SunJ200RPosAng'] = SunJ200RPosAng  # as expected! no angle diff
    SunMarsRPosAng = Sun1.position_angle(MarsR).deg * u.deg
    eph['SunMarsRPosAng'] = SunMarsRPosAng  # ~36.7

    '''
    Make small offset by 0.1 degree to point in the direction of Mars North
    '''
    # SMN1 = Sun1.directional_offset_by(0*u.deg, 0.1*u.deg)  # North!
    SMN1 = Sun1.directional_offset_by(SunMarsRPosAng, 0.1*u.deg)
    eph['SMN_RA1'] = SMN1.ra.deg
    eph['SMN_DEC1'] = SMN1.dec.deg

    '''Offsets from center of solar disk'''
    Obj1frame = Sun1.skyoffset_frame()
    Terra1Off = Terra1.transform_to(Obj1frame)
    Luna1Off = Luna1.transform_to(Obj1frame)
    SMN1Off = SMN1.transform_to(Obj1frame)

    return (eph, zTime, SunR, TerraR, LunaR, Sun1, Terra1, Luna1, SMN1,
            Terra1Off, Luna1Off, SMN1Off)


def plotDisks(r, ra, dec):
    theta = np.linspace(0, 2*np.pi, 100)
    x = r*np.cos(theta) + ra
    y = r*np.sin(theta) + dec
    return x, y


def makeSolarDisk(xpdfname, xdrift=True):
    '''Transit across the solar disk'''

    '''
    Options for coordinate frames: only J2000 (1) for now
    '''
    Sun = Sun1
    TerraOff = Terra1Off
    LunaOff = Luna1Off
    SMNOff = SMN1Off

    '''
    Options for plots
    '''
    if xdrift:
        xtitleappend = '\nSolar disk moves across the sky'
        xlabstr = 'RA (degrees)'
        ylabstr = 'DEC (degrees)'
        legendloc = 'upper right'
        Centre = Sun
    else:
        xtitleappend = '\nSolar disk centered'
        xlabstr = r'$\Delta$'+' RA (degrees)'
        ylabstr = r'$\Delta$'+' DEC (degrees)'
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
                   color=color, marker='+',
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
        x = SMNOff[ii].lon.deg + Centre[ii].ra.deg
        y = SMNOff[ii].lat.deg + Centre[ii].dec.deg
        ax.scatter(x, y, color=color, alpha=0.9, marker=10)

    ax.set_aspect(1)
    ax.legend(fontsize='x-small', loc=legendloc)
    ax.invert_xaxis()
    ax.set(xlabel=xlabstr, ylabel=ylabstr)
    plt.figtext(0.5, 0.005,
                ('Note: In J2000 frame. +--' + r'$\Delta$'
                 " points to Mars' north (rotational axis)"),
                wrap=True, horizontalalignment='center')
    plt.title('Terra & Luna Transit from Mars 2084-11-10 UTC' + xtitleappend)

    plt.savefig(xpdfname + '.pdf')
    print(xpdfname + '.pdf saved')


def makeWholeSky(xpdfname):
    ''' whole sky view for the transit events '''

    '''
    Options for sets: only 1 for now
    '''
    Sun = Sun1
    Terra = Terra1
    Luna = Luna1
    SMN = SMN1

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='mollweide')
    # each event gets its own colour
    colormap = iter(plt.cm.rainbow(np.linspace(0, 1, len(eph))))
    for ii in range(len(eph)):
        color = next(colormap)
        x = Sun[ii].ra.wrap_at(180*u.deg).rad
        y = Sun[ii].dec.rad
        ax.scatter(x, y, color=color, marker='+', label=eph['Object'][ii] +
                   " " + eph['Event'][ii] + " " + zTime[ii])
        x, y = plotDisks(r=SunR.rad[ii],
                         ra=Sun[ii].ra.wrap_at(180*u.deg).rad,
                         dec=Sun[ii].dec.rad)
        ax.plot(x, y, color=color, alpha=0.5)
        x, y = plotDisks(r=TerraR.rad[ii],
                         ra=Terra[ii].ra.wrap_at(180*u.deg).rad,
                         dec=Terra[ii].dec.rad)
        ax.plot(x, y, color=color, alpha=0.9)
        x, y = plotDisks(r=LunaR.rad[ii],
                         ra=Luna[ii].ra.wrap_at(180*u.deg).rad,
                         dec=Luna[ii].dec.rad)
        ax.plot(x, y, color=color, alpha=0.9)
        x = SMN[ii].ra.wrap_at(180*u.deg).rad
        y = SMN[ii].dec.rad
        ax.scatter(x, y, color=color, alpha=0.3, marker=10)

    ax.legend(fontsize='x-small', loc='upper right')
    ax.grid(True)
    plt.figtext(0.5, 0.15, 'Note: whole sky view. In J2000 frame',
                wrap=True, horizontalalignment='center')
    plt.title('Terra & Luna Transit from Mars 2084-11-10 UTC',
              fontweight='bold')

    plt.savefig(xpdfname+'.pdf')
    print(xpdfname + '.pdf saved')


(eph, zTime, SunR, TerraR, LunaR, Sun1, Terra1,
 Luna1, SMN1, Terra1Off, Luna1Off, SMN1Off) = readEPH()

makeSolarDisk(xpdfname='Output/SolarDiskTravel', xdrift=True)
makeSolarDisk(xpdfname='Output/SolarDiskCentered', xdrift=False)

makeWholeSky(xpdfname='Output/WholeSky')
