# make pretty Solar Disk overview plots

from astropy.table import Table
import astropy.coordinates as SkyCoord
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

def readPlotTable():
    print('read Event Table')
    return Table.read('GeocentricPlotTable.ecsv', format = 'ascii.ecsv')

GeoPlotTable = readPlotTable()
GeoPlotTable.pprint_all()

# ra = SkyCoord.Angle(GeoPlotTable['S_RA_app'])
# ra = ra.wrap_at(180*u.degree)
# dec = SkyCoord.Angle(GeoPlotTable['S_DEC_app'])
# fig = plt.figure(figsize=(8,6))
# # ax = fig.add_subplot(111,projection='mollweide')
# ax = fig.add_subplot(111)
# ax.scatter(ra.radian, dec.radian)
# plt.xlabel('RA')
# plt.ylabel('DEC')
# plt.show()
# plt.close()

ijk = []
for ii in range(len(GeoPlotTable)):
    if GeoPlotTable['Event'][ii].startswith('Tflag'):
        ijk.append(ii)

GeoPlotTable.remove_rows(ijk)

zTime = GeoPlotTable['datetime_str']
zTime = [s.replace('2084-Nov-10 ','') for s in zTime]
zTime = [s.replace('.000','') for s in zTime]

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)
colormap = iter(plt.cm.rainbow(np.linspace(0,1,len(GeoPlotTable))))
ax.scatter(0,0, color='gray', marker='+')
for ii in range(len(GeoPlotTable)):
    color=next(colormap)
    theta = np.linspace(0, 2*np.pi, 100)
    r = (GeoPlotTable['S_ang_width'][ii]/2)/3600
    x = r*np.cos(theta) + GeoPlotTable['S_RA_app'][ii]-GeoPlotTable['S_RA_app'][ii]
    y = r*np.sin(theta) + GeoPlotTable['S_DEC_app'][ii]-GeoPlotTable['S_DEC_app'][ii]
    ax.plot(-x,y,color=color,alpha=0.5, label = GeoPlotTable['Object'][ii] + " "+ GeoPlotTable['Event'][ii] + " "+ zTime[ii])
    r = (GeoPlotTable['T_ang_width'][ii]/2)/3600
    x = r*np.cos(theta) + GeoPlotTable['T_RA_app'][ii]-GeoPlotTable['S_RA_app'][ii]
    y = r*np.sin(theta) + GeoPlotTable['T_DEC_app'][ii]-GeoPlotTable['S_DEC_app'][ii]
    ax.plot(-x,y,color=color,alpha=0.9)
    r = (GeoPlotTable['L_ang_width'][ii]/2)/3600
    x = r*np.cos(theta) + GeoPlotTable['L_RA_app'][ii]-GeoPlotTable['S_RA_app'][ii]
    y = r*np.sin(theta) + GeoPlotTable['L_DEC_app'][ii]-GeoPlotTable['S_DEC_app'][ii]
    ax.plot(-x,y,color=color,alpha=0.9)

ax.set_aspect(1)
ax.legend(fontsize='x-small')
ax.set(xlabel='-RA (degrees)', ylabel='DEC (degrees)')
plt.figtext(0.5,0.005,'RA flipped to match what we see from ground. Angles as Cartesian distances give plotting error', wrap=True, horizontalalignment='center')
plt.title('Terra & Luna Transit from Mars 2084-11-10 UTC\nSolar Disk centered', fontweight = 'bold')
plt.savefig('1.png')
# plt.show()
# plt.close()

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)
colormap = iter(plt.cm.rainbow(np.linspace(0,1,len(GeoPlotTable))))
for ii in range(len(GeoPlotTable)):
    color=next(colormap)
    ax.scatter(GeoPlotTable['S_RA_app'][ii],GeoPlotTable['S_DEC_app'][ii], color=color, marker='+')
    theta = np.linspace(0, 2*np.pi, 100)
    r = (GeoPlotTable['S_ang_width'][ii]/2)/3600
    x = r*np.cos(theta) + GeoPlotTable['S_RA_app'][ii]
    y = r*np.sin(theta) + GeoPlotTable['S_DEC_app'][ii]
    ax.plot(x,y,color=color,alpha=0.5, label = GeoPlotTable['Object'][ii] + " "+ GeoPlotTable['Event'][ii] + " "+ zTime[ii])
    r = (GeoPlotTable['T_ang_width'][ii]/2)/3600
    x = r*np.cos(theta) + GeoPlotTable['T_RA_app'][ii]
    y = r*np.sin(theta) + GeoPlotTable['T_DEC_app'][ii]
    ax.plot(x,y,color=color,alpha=0.9)
    r = (GeoPlotTable['L_ang_width'][ii]/2)/3600
    x = r*np.cos(theta) + GeoPlotTable['L_RA_app'][ii]
    y = r*np.sin(theta) + GeoPlotTable['L_DEC_app'][ii]
    ax.plot(x,y,color=color,alpha=0.9)

ax.set_aspect(1)
ax.legend(fontsize='x-small')
ax.set(xlabel='RA (degrees)', ylabel='DEC (degrees)')
plt.figtext(0.5,0.005, 'Angles as Cartesian distances give plotting error', wrap=True, horizontalalignment='center')
plt.title('Terra & Luna Transit from Mars 2084-11-10 UTC', fontweight = 'bold')
plt.savefig('2.png')
# plt.show()
# plt.close()
