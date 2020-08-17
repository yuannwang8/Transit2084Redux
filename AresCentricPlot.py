# make pretty Solar Disk overview plots

from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

import T2084Helpers.MarsNAlign as MNA
import T2084Helpers.misc as misc

'''
Purpose: This file plots the relative positions of the Earth and the Moon as
they transit across the Solar Disk.
Commentary: The direction of Mars's north (rotational axis) looks odd
when considering Earth's orbital plane.
To do: To consider plotting the invariant plane for further clarity.
'''


def readEPH(useMNA, rused):
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

    '''Mars rotational axis & J2000 rotational axis'''
    # use MarsRRA and MarsRDEC defined in MNA module
    MarsRRA = MNA.MarsRRA
    MarsRDEC = MNA.MarsRDEC
    # MarsRRA.unit = u.deg
    # MarsRDEC.unit = u.deg
    MarsR = SkyCoord(MarsRRA, MarsRDEC, frame='icrs', unit='deg')
    # J2000R = SkyCoord([0]*len(eph), [90]*len(eph), frame='icrs', unit='deg')
    # RotSept = J2000R.separation(MarsR).deg * u.deg
    # RotPosAng = J2000R.position_angle(MarsR).deg * u.deg
    # eph['RotSept'] = RotSept
    # eph['RotPosAng'] = RotPosAng

    '''
    Objects in SkyCoords objects:
    Option to switch between J2000 frame or Mars Rotational Axis aligned frame
    '''
    Sun = SkyCoord(eph['S_RA_app'], eph['S_DEC_app'],
                   frame='icrs', unit='deg')
    Terra = SkyCoord(eph['T_RA_app'], eph['T_DEC_app'],
                     frame='icrs', unit='deg')
    Luna = SkyCoord(eph['L_RA_app'], eph['L_DEC_app'],
                    frame='icrs', unit='deg')

    if useMNA:
        Sun = Sun.transform_to(MNA.MarsNAlign)
        Terra = Terra.transform_to(MNA.MarsNAlign)
        Luna = Luna.transform_to(MNA.MarsNAlign)
        MarsR = MarsR.transform_to(MNA.MarsNAlign)

    '''
    Make small offset by 0.1 degree to point in the direction of Mars North
    '''
    # SunJ200RPosAng = Sun.position_angle(J2000R).deg * u.deg
    # eph['SunJ200RPosAng'] = SunJ200RPosAng  # as expected! no angle diff
    SunMarsRPosAng = Sun.position_angle(MarsR).deg * u.deg
    eph['SunMarsRPosAng'] = SunMarsRPosAng  # ~36.7
    # SMN = Sun.directional_offset_by(0*u.deg, 0.1*u.deg)  # North!
    SMN = Sun.directional_offset_by(SunMarsRPosAng, 0.1*u.deg)  # good&weird

    '''
    Position angle offsets from center of solar disk
    To add constraints to rused = [0,1,-1 only]
    '''
    if useMNA or rused == 0:
        '''
        This is for the original J2000 original frame or
        of the MNA.MarsNAlign transformation is used
        '''
        rused = 0
    elif rused == 1:
        '''
        This rotates the offsets.
        When using the J2000 frame, this is functionlly equivalent to
        MNA.MarsAlign without a position angle offset.
        '''
        rused = SunMarsRPosAng
    elif rused == -1:
        '''
        This rotates the offsets the otherway... This is not very logical but
        exists as an option for future considerations.
        '''
        rused = -SunMarsRPosAng
    else:
        print('Your selected option for angle offset is not [0, 1, -1].' +
              'Defaulting to 0')
        rused = 0

    Obj1frame = Sun.skyoffset_frame(rused)
    TerraOff = Terra.transform_to(Obj1frame)
    LunaOff = Luna.transform_to(Obj1frame)
    SMNOff = SMN.transform_to(Obj1frame)

    return (eph, zTime, SunR, TerraR, LunaR, Sun, Terra, Luna, SMN,
            TerraOff, LunaOff, SMNOff)


def plotDisks(r, ra, dec):
    theta = np.linspace(0, 2*np.pi, 100)
    x = r*np.cos(theta) + ra
    y = r*np.sin(theta) + dec
    return x, y


def makeSolarDisk(xpdfname, xdrift=True, useMNAl=True):
    '''
    Transit across the solar disk
    Options for plots
    '''
    # update labelling options
    if useMNAl:
        xtag = ("Note: Frame aligned to Mars' rotational axis. +--" +
                r'$\Delta$'+" points to Mars' north")
        xlabtag = 'Lon (degrees)'
        ylabtag = 'Lat (degrees)'
    else:
        xtag = ('Note: In J2000 frame. +--' + r'$\Delta$'
                " points to Mars' north")
        xlabtag = 'RA (degrees)'
        ylabtag = 'DEC (degrees)'

    if xdrift:
        xtitleappend = '\nSolar disk moves across the sky'
        xlabstr = xlabtag
        ylabstr = ylabtag
        legendloc = 'upper right'
        Centre = Sun
    else:
        xtitleappend = '\nSolar disk centered'
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
    # if xdrift:  # this doesn't work
    #     overlay = ax.get_coords_overlay('icrs')
    #     overlay.grid(color='blue', ls='dotted')
    #     overlay[0].set_axislabel('RA-like')
    #     overlay[1].set_axislabel('DEC-like')
    plt.figtext(0.5, 0.025, xtag, wrap=True, horizontalalignment='center')
    plt.figtext(0.5, 0.01, misc.xfooter, size=6,
                horizontalalignment='center')
    plt.title('Terra & Luna Transit from Mars 2084-11-10 UTC' + xtitleappend)

    plt.savefig(xpdfname + '.pdf')
    print(xpdfname + '.pdf saved')
    plt.close()


def makeWholeSky(xpdfname, useMNAl):
    ''' whole sky view for the transit events '''

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='mollweide')
    # update labelling options
    if useMNAl:
        xtag = 'Note: whole sky view. In rotated J2000 frame'
    else:
        xtag = 'Note: whole sky view. In J2000 frame'

    # each event gets its own colour
    colormap = iter(plt.cm.plasma(np.linspace(0, 1, len(eph))))
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
    plt.figtext(0.5, 0.15, xtag, wrap=True, horizontalalignment='center')
    plt.figtext(0.5, 0.1, misc.xfooter, size=6,
                horizontalalignment='center')
    plt.title('Terra & Luna Transit from Mars 2084-11-10 UTC',
              fontweight='bold')

    plt.savefig(xpdfname+'.pdf')
    print(xpdfname + '.pdf saved')
    plt.close()


(eph, zTime, SunR, TerraR, LunaR, Sun, Terra,
 Luna, SMN, TerraOff, LunaOff, SMNOff) = readEPH(useMNA=False, rused=0)
makeSolarDisk(xpdfname='Output/SolarTravelJ2000', xdrift=True, useMNAl=False)
makeSolarDisk(xpdfname='Output/SolarOffsetJ2000', xdrift=False, useMNAl=False)
makeWholeSky(xpdfname='Output/WholeSkyJ2000', useMNAl=False)

# use with MNA.MarsNAlign
(eph, zTime, SunR, TerraR, LunaR, Sun, Terra,
 Luna, SMN, TerraOff, LunaOff, SMNOff) = readEPH(useMNA=True, rused=0)
makeSolarDisk(xpdfname='Output/SolarTravelMarsN', xdrift=True, useMNAl=True)
makeSolarDisk(xpdfname='Output/SolarOffsetMarsN', xdrift=False, useMNAl=True)
makeWholeSky(xpdfname='Output/WholeSkyMars', useMNAl=True)

# use skyframe offset 1: only the offset chart makes sense!
(eph, zTime, SunR, TerraR, LunaR, Sun, Terra,
 Luna, SMN, TerraOff, LunaOff, SMNOff) = readEPH(useMNA=False, rused=1)
# makeSolarDisk(xpdfname='Output/SolarTravel2', xdrift=True, useMNAl=True)
makeSolarDisk(xpdfname='Output/SolarOffset2', xdrift=False, useMNAl=True)

# use skyframe offset -1: this option is here for curiosity purposes only.
(eph, zTime, SunR, TerraR, LunaR, Sun, Terra,
 Luna, SMN, TerraOff, LunaOff, SMNOff) = readEPH(useMNA=False, rused=-1)
# makeSolarDisk(xpdfname='Output/SolarTravel3wrong', xdrift=True, useMNAl=True)
makeSolarDisk(xpdfname='Output/SolarOffset3wrong', xdrift=False, useMNAl=True)
