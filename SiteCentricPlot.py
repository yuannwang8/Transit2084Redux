# make pretty Solar Disk overview plots from the surface

from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

import T2084Helpers.misc as misc

'''
Purpose: This file plots the relative positions of the Earth and the Moon as
they transit across the Solar Disk.
Yep, the AltAz frame rotation is weird and may need reconsidering.
'''


def readEPH(sitename):
    '''Load and transform Event Table'''
    eph = Table.read('Data/G/'+sitename+'ETable.ecsv', format='ascii.ecsv')

    '''
    remove rows for Phobos and Deimos
    '''
    ijk = []
    for ii in range(len(eph)):
        if eph['Object'][ii].startswith(('D', 'P')):
            ijk.append(ii)
    eph.remove_rows(ijk)

    '''
    Remove sunrise and sunset rows if the occur before T-Ingress and
    after L-Egress
    '''
    if eph['Object'][0].startswith('S'):
        eph.remove_rows(0)

    if eph['Object'][len(eph)-1].startswith('S'):
        eph.remove_rows(len(eph)-1)

    '''
    get event timings & labels
    '''
    zTime = eph['datetime_str']
    zTime = [s.replace('2084-Nov-10 ', '') for s in zTime]

    '''Object radii in degrees'''
    SolR = Angle((eph['S_ang_width']/2))
    TerraR = Angle((eph['T_ang_width']/2))
    LunaR = Angle((eph['L_ang_width']/2))

    '''
    Objects in SkyCoords objects
    (Can't do offset in altaz, it seems...) It really should be altaz instead
    '''
    Sol = SkyCoord(eph['S_AZ'], eph['S_EL'], frame='icrs', unit='deg')
    Terra = SkyCoord(eph['T_AZ'], eph['T_EL'], frame='icrs', unit='deg')
    Luna = SkyCoord(eph['L_AZ'], eph['L_EL'], frame='icrs', unit='deg')
    # Phobos = SkyCoord(eph['P_AZ'], eph['P_EL'], frame='altaz', unit='deg')
    # Deimos = SkyCoord(eph['D_AZ'], eph['D_EL'], frame='altaz', unit='deg')

    '''
    Make small offset by 0.1 degree to point in the direction of Mars North
    '''
    MarsN = SkyCoord([0]*len(eph), [0]*len(eph), frame='icrs', unit='deg')
    MarsNPoint = Sol.position_angle(MarsN).deg * u.deg
    SMN = Sol.directional_offset_by(MarsNPoint, 0.1*u.deg)  # literally 0.1N

    Obj1frame = Sol.skyoffset_frame(MarsNPoint)
    TerraOff = Terra.transform_to(Obj1frame)
    LunaOff = Luna.transform_to(Obj1frame)
    SMNOff = SMN.transform_to(Obj1frame)

    return (eph, zTime, SolR, TerraR, LunaR, Sol, Terra, Luna, SMN,
            TerraOff, LunaOff, SMNOff)


def plotDisks(r, ra, dec):
    theta = np.linspace(0, 2*np.pi, 100)
    x = r*np.cos(theta) + ra
    y = r*np.sin(theta) + dec
    return x, y


def makeSolarDisk(ax, sitename, locale, eph, zTime, SolR, TerraR, LunaR,
                  Sol, TerraOff, LunaOff, SMNOff, xdrift=True):
    '''
    subset observations to those between sunrise and sunset only
    if no observations, return an empty box
    '''
    mask = eph['S_solar_presence'] == '*'

    if any(mask):
        eph = eph[mask]
        zTime = [i for i, m in zip(zTime, mask.data) if m]
        SolR = SolR[mask]
        TerraR = TerraR[mask]
        LunaR = LunaR[mask]
        Sol = Sol[mask]
        TerraOff = TerraOff[mask]
        LunaOff = LunaOff[mask]
        SMNOff = SMNOff[mask]

        '''
        Transit across the solar disk
        '''
        xlabtag = 'Azimuth (degrees)'
        ylabtag = 'Altitude (degrees)'

        if xdrift:
            xtitleappend = 'Solar disk moves across the sky'
            xlabstr = xlabtag
            ylabstr = ylabtag
            legendloc = 'best'  # 'upper right'
            Centre = Sol
            xtag = ' looking South-ish'
        else:
            xtitleappend = 'Solar disk centered'
            xlabstr = r'$\Delta$'+' '+xlabtag
            ylabstr = r'$\Delta$'+' '+ylabtag
            legendloc = 'best'  # 'lower left'
            Centre = SkyCoord([0]*len(Sol), [0]*len(Sol),
                              frame='icrs', unit='deg')
            xtag = r'$\Delta$'+"points to Mars's North"

        # each event gets its own colour
        colormap = iter(plt.cm.rainbow(np.linspace(0, 1, len(eph))))
        for ii in range(len(eph)):
            color = next(colormap)
            ax.scatter(Centre[ii].ra.deg, Centre[ii].dec.deg,
                       color=color, alpha=0.9, marker='+',
                       label=eph['Object'][ii] + " " +
                       eph['Event'][ii] + " " + zTime[ii])
            x, y = plotDisks(r=SolR.deg[ii],
                             ra=Centre[ii].ra.deg,
                             dec=Centre[ii].dec.deg)
            ax.plot(x, y, color=color, alpha=0.3)
            x, y = plotDisks(r=TerraR.deg[ii],
                             ra=TerraOff[ii].lon.deg + Centre[ii].ra.deg,
                             dec=TerraOff[ii].lat.deg + Centre[ii].dec.deg)
            ax.plot(x, y, color=color, alpha=0.9)
            x, y = plotDisks(r=LunaR.deg[ii],
                             ra=LunaOff[ii].lon.deg + Centre[ii].ra.deg,
                             dec=LunaOff[ii].lat.deg + Centre[ii].dec.deg)
            ax.plot(x, y, color=color, alpha=0.9)
            if not xdrift:
                x = SMNOff[ii].lon.deg + Centre[ii].ra.deg
                y = SMNOff[ii].lat.deg + Centre[ii].dec.deg
                ax.scatter(x, y, color=color, alpha=0.3, marker=10)

        ax.set_aspect(1)
        ax.legend(fontsize='x-small', loc=legendloc)
        # ax.invert_xaxis()
        ax.set(xlabel=xlabstr, ylabel=ylabstr)
        ax.set_title(xtitleappend+' '+xtag)
        return ax

    else:
        x = [0, 1]
        y = [0, 1]
        txt = 'Sun not visible during transit:\nsolar disk plot not plotted'
        ax.scatter(x, y, alpha=0.1)
        ax.annotate(txt, (0.27, 0.5))
        ax.set_aspect(1)
        # ax.set_axis_off()
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        return ax


def makeWholeSky(ax, sitename, locale, eph, zTime, Sol):
    '''
    Whole sky view for the transit events
    Looks south - this requires the (angle.ra.deg - 180) as a wrap-around
    Plot only center of solar disk:
    the disks themselves are too small to be seen!
    '''

    # each event gets its own colour
    colormap = iter(plt.cm.plasma(np.linspace(0, 1, len(eph))))
    for ii in range(len(eph)):
        color = next(colormap)
        # x = Sol[ii].ra.wrap_at(180*u.deg).rad
        x = Angle((Sol[ii].ra.deg - 180)*u.deg).wrap_at(180*u.deg).rad
        y = Sol[ii].dec.rad
        ax.scatter(x, y, color=color, marker='+', label=eph['Object'][ii] +
                   " " + eph['Event'][ii] + " " + zTime[ii])

    ax.legend(fontsize='x-small', loc='best')  # 'lower center'
    ax.grid(True)
    ax.set_xticklabels(['30°', '60°', '90°', '120°', '150°', '180°',
                        '210°', '240°', '270°', '300°', '330°'])
    ax.set_title('Whole sky view looking south')
    return ax


def siteCentricPlot(surfViz):
    sitelist = surfViz['Site'].data
    location = []
    for ijk in range(len(surfViz)):
        location.append(str(np.round(-surfViz['Elondeg'][ijk], 2)) + 'E ' +
                        str(np.round(surfViz['Latdeg'][ijk], 2))+'N')

    for site, locale in zip(sitelist, location):
        (eph, zTime, SolR, TerraR, LunaR, Sol, Terra,
         Luna, SMN, TerraOff, LunaOff, SMNOff) = readEPH(site)

        plt.figure(figsize=(16, 6))
        ax1 = plt.subplot(1, 2, 1)
        ax2 = plt.subplot(1, 2, 2, projection='hammer')

        # set xdrift=True to check that the whole sky plot is showing the
        # correct azimuth
        ax1 = makeSolarDisk(ax1, site, locale, eph, zTime, SolR, TerraR,
                            LunaR, Sol, TerraOff, LunaOff, SMNOff,
                            xdrift=False)

        ax2 = makeWholeSky(ax2, site, locale, eph, zTime, Sol)

        # plt.tight_layout()
        plt.suptitle('Terra & Luna Transit from Mars 2084-11-10 UTC at ' +
                     site+' '+locale, fontweight='bold')
        plt.figtext(0.5, 0.05, misc.xfooter,
                    size=6, horizontalalignment='center')
        plt.savefig('Output/G/'+site+'.pdf')
        print(site+' saved')
        plt.close()


surfViz = Table.read('Data/AresSurfaceViz.ecsv', format='ascii.ecsv')
mask = surfViz['Site'] == 'EZ22EZ40GaleCrater'
surfViz = surfViz[np.logical_not(mask)]

siteCentricPlot(surfViz)

# siteCentricPlot(sitename='HCuriosity')
# siteCentricPlot(sitename='EZ11aErebusMontes')
