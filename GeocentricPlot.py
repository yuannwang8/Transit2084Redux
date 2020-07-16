# make pretty Solar Disk overview plots

from astropy.table import Table
from astropy.coordinates import SkyCoord
# from astropy.time import Time
import astropy.units as u

import numpy as np
import matplotlib.pyplot as plt


def readPlotTable():
    print('read Event Table')
    return Table.read('GeocentricPlotTable.ecsv', format='ascii.ecsv')


GeoPlotTable = readPlotTable()


# remove rows starting with Tflag - keep calculated rows only
ijk = []
for ii in range(len(GeoPlotTable)):
    if GeoPlotTable['Event'][ii].startswith('Tflag'):
        ijk.append(ii)

GeoPlotTable.remove_rows(ijk)


'''get time of events'''
zTime = GeoPlotTable['datetime_str']
zTime = [s.replace('2084-Nov-10 ', '') for s in zTime]
zTime = [s.replace('.000', '') for s in zTime]


'''Mars rotational axis & J2000 rotational axis'''
Tref = 2451545  # reference JD for '2000-01-01T12:00:00'
# print(Tref.jd)
MarsRRA = (317.681 - 0.106*(GeoPlotTable['datetime_jd']-Tref)/36525)
MarsRDEC = (52.887 - 0.061*(GeoPlotTable['datetime_jd']-Tref)/36525)
MarsRRA.unit = u.deg
MarsRDEC.unit = u.deg
MarsR = SkyCoord(MarsRRA, MarsRDEC, frame='icrs', unit='deg')

J2000R = SkyCoord([0]*len(GeoPlotTable), [90]*len(GeoPlotTable),
                  frame='icrs', unit='deg')

# RotSept = J2000R.separation(MarsR).deg * u.deg
# RotPosAng = J2000R.position_angle(MarsR).deg * u.deg
# GeoPlotTable['RotSept'] = RotSept
# GeoPlotTable['RotPosAng'] = RotPosAng


'''
Object in SkyCoords objects:
Object1: J2000 frame. Object2: Mars rotational frame.
'''
Sun1 = SkyCoord(GeoPlotTable['S_RA_app'], GeoPlotTable['S_DEC_app'],
                frame='icrs', unit='deg')
SunJ200RPosAng = Sun1.position_angle(J2000R).deg * u.deg
SunMarsRPosAng = Sun1.position_angle(MarsR).deg * u.deg

# GeoPlotTable['SunJ200RPosAng'] = SunJ200RPosAng  # as expected! no angle diff
GeoPlotTable['SunMarsRPosAng'] = SunMarsRPosAng  # ~36.7

Terra1 = SkyCoord(GeoPlotTable['T_RA_app'], GeoPlotTable['T_DEC_app'],
                  frame='icrs', unit='deg')
Luna1 = SkyCoord(GeoPlotTable['L_RA_app'], GeoPlotTable['L_DEC_app'],
                 frame='icrs', unit='deg')


'''Object Radii in degrees'''
SunR = (GeoPlotTable['S_ang_width']/2)/3600
TerraR = (GeoPlotTable['T_ang_width']/2)/3600
LunaR = (GeoPlotTable['L_ang_width']/2)/3600

'''Make small offset by 0.1 degree to point in the direction of Mars North'''
# SMN = Sun1.directional_offset_by(0*u.deg, 0.1*u.deg)  # North!
# SMN = Sun1.directional_offset_by(45*u.deg, 0.1*u.deg)  # test
# SMN = Sun1.directional_offset_by(36.7*u.deg, 0.1*u.deg)
SMN = Sun1.directional_offset_by(SunMarsRPosAng, 0.1*u.deg)
GeoPlotTable['SMN_RA'] = SMN.ra.degree
GeoPlotTable['SMN_DEC'] = SMN.dec.degree

# GeoPlotTable.pprint_all()


'''
Things to update
2a. install ligo.skymap package for small plots
2b. with offsets
'''


def plotDisks(r, ra, dec):
    x = r*np.cos(theta) + ra
    y = r*np.sin(theta) + dec
    return x, y


# Plot 1: solar disk centered
fig = plt.figure(figsize=(8, 6))  # (14,6) for side by side panels
ax = fig.add_subplot(111)  # (121) for side by side panels
# center of solar disk
ax.scatter(0, 0, color='gray', marker='+')
# each event gets its own colour
colormap = iter(plt.cm.rainbow(np.linspace(0, 1, len(GeoPlotTable))))
theta = np.linspace(0, 2*np.pi, 100)
for ii in range(len(GeoPlotTable)):
    color = next(colormap)
    x, y = plotDisks(r=SunR[ii], ra=0, dec=0)
    ax.plot(x, y, color=color, alpha=0.5,
            label=GeoPlotTable['Object'][ii] + " " +
            GeoPlotTable['Event'][ii] + " " + zTime[ii])
    x, y = plotDisks(r=TerraR[ii],
                     ra=Terra1[ii].ra.degree - Sun1[ii].ra.degree,
                     dec=Terra1[ii].dec.degree - Sun1[ii].dec.degree)
    ax.plot(x, y, color=color, alpha=0.9)
    x, y = plotDisks(r=LunaR[ii],
                     ra=Luna1[ii].ra.degree - Sun1[ii].ra.degree,
                     dec=Luna1[ii].dec.degree - Sun1[ii].dec.degree)
    ax.plot(x, y, color=color, alpha=0.9)
    x = SMN[ii].ra.degree - Sun1[ii].ra.degree
    y = SMN[ii].dec.degree - Sun1[ii].dec.degree
    ax.scatter(x, y, color=color, alpha=0.9, marker=10)

ax.set_aspect(1)
ax.legend(fontsize='x-small', loc='lower left')
ax.invert_xaxis()  # invert x-axis
ax.set(xlabel=r'$\Delta$'+' RA (degrees)', ylabel=r'$\Delta$'+' DEC (degrees)')
plt.figtext(0.5, 0.005,
            'Note: plotting error with angles as Cartesian',
            wrap=True, horizontalalignment='center')
plt.title('Terra & Luna Transit from Mars 2084-11-10 UTC\nSolar Disk centered',
          fontweight='bold')
plt.savefig('1.png')


# Plot 2: disks moving across the sky
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection=None)  # (122) for side by side panels
colormap = iter(plt.cm.rainbow(np.linspace(0, 1, len(GeoPlotTable))))
theta = np.linspace(0, 2*np.pi, 100)
for ii in range(len(GeoPlotTable)):
    color = next(colormap)
    ax.scatter(Sun1[ii].ra.degree, Sun1[ii].dec.degree,
               color=color, marker='+')
    x, y = plotDisks(r=SunR[ii], ra=Sun1[ii].ra.degree,
                     dec=Sun1[ii].dec.degree)
    ax.plot(x, y, color=color, alpha=0.5,
            label=GeoPlotTable['Object'][ii] + " " +
            GeoPlotTable['Event'][ii] + " " + zTime[ii])
    x, y = plotDisks(r=TerraR[ii], ra=Terra1[ii].ra.degree,
                     dec=Terra1[ii].dec.degree)
    ax.plot(x, y, color=color, alpha=0.9)
    x, y = plotDisks(r=LunaR[ii], ra=Luna1[ii].ra.degree,
                     dec=Luna1[ii].dec.degree)
    ax.plot(x, y, color=color, alpha=0.9)
    x = SMN.ra.degree[ii]
    y = SMN.dec.degree[ii]
    ax.scatter(x, y, color=color, alpha=0.9, marker=10)

ax.set_aspect(1)
ax.legend(fontsize='x-small', loc='upper right')
ax.invert_xaxis()  # invert x-axis
ax.set(xlabel='RA (degrees)', ylabel='DEC (degrees)')
plt.figtext(0.5, 0.005, 'Note: plotting error with angles as Cartesian',
            wrap=True, horizontalalignment='center')
plt.title('Terra & Luna Transit from Mars 2084-11-10 UTC', fontweight='bold')
plt.savefig('2.png')
