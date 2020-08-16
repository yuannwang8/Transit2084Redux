# find duration and types of surface viewing

from astropy.table import Table, join
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

import T2084Helpers.misc as misc


def surfaceVisibility():
    '''
    Identifty key events per surface site surveyed
    For the sun: the following events are recorded: sunrise, sunset
    For Earth and the Moon, the following states of transits are recorded:
    Full Viz, Rises during Transit, Sets during Transit, Zero Viz
    '''
    TA = []
    TB = []
    TC = []
    TD = []
    TE = []
    for key, grp in zip(surf.groups.keys, surf.groups):
        # extract times that are not Sol dependent
        # TL Outer Ingress & Egress
        blT = grp['Object'] == 'T'
        blL = grp['Object'] == 'L'
        blOutIng = grp['Event'] == 'Outer_Ingress'
        blOutEg = grp['Event'] == 'Outer_Egress'
        TerraIn = grp[blT & blOutIng]['datetime_jd'][0]
        TerraOut = grp[blT & blOutEg]['datetime_jd'][0]
        LunaIn = grp[blL & blOutIng]['datetime_jd'][0]
        LunaOut = grp[blL & blOutEg]['datetime_jd'][0]

        # Full Viz, Rises during Transit, Sets during Transit, Zero Viz
        TA.append(key[0])
        # Is the Sun visible during the transit time frame?
        if 'S' in grp['Object']:
            TB.append('Visible')
            SunRise = grp[grp['Event'] == 'sunrise']['datetime_jd'][0]
            SunSet = grp[grp['Event'] == 'sunset']['datetime_jd'][0]

            # Terra event options
            if SunRise >= TerraOut:  # sun rises after transit ends
                TC.append('Zero Visibility')
                Tstart = SunRise
                Tend = SunRise
            elif SunRise <= TerraIn:  # sun rises before transit starts
                if SunSet >= TerraOut:
                    TC.append('Full Visibility')
                    Tend = TerraOut
                else:
                    TC.append('Earth sets during transit')
                    Tend = SunSet
                Tstart = TerraIn
            else:
                if SunSet >= TerraOut:
                    TC.append('Earth rises during transit')
                    Tend = TerraOut
                else:
                    TC.append('Earth rises and sets during transit')
                    Tend = SunSet
                Tstart = SunRise

            # Luna event options
            if SunSet <= LunaIn:  # sun sets before transit ends
                TD.append('Zero Visibility')
                Lstart = SunSet
                Lend = SunSet
            elif SunRise <= LunaIn:  # sun rises before transit starts
                if SunSet >= LunaOut:
                    TD.append('Full Visibility')
                    Lend = LunaOut
                else:
                    TD.append('Moon sets during transit')
                    Lend = SunSet
                Lstart = LunaIn
            else:
                if SunSet >= LunaOut:
                    TD.append('Moon rises during transit')
                    Lend = LunaOut
                else:
                    TD.append('Moon rises and sets during transit')
                    Lend = SunSet
                Lstart = SunRise

            # if Earth Zero vis
            if Tend - Tstart == 0:
                TE.append(Lend - Lstart)
            elif Lend - Lstart == 0:
                TE.append(Tend - Tstart)
            else:
                TE.append(Lend - Tstart)

        else:
            TB.append('Not Visible')
            TC.append('Zero Visibility')
            TD.append('Zero Visibility')
            TE.append(0)

    Tout = Table([TA, TB, TC, TD, TE],
                 names=('Site', 'Sun', 'TerraObs', 'LunaObs', 'EventDuration'))

    '''
    Read in longitude and latitude of viewing locations and
    add that to table and returning it as separate object
    '''
    SrfSites = Table.read('MarsSurfaceSites.csv', format='csv')
    SrfSites.remove_columns(['SiteName', 'HorizonsLocation',
                             'Altkm', 'Sources'])
    Tout = join(Tout, SrfSites)
    Tout['Elondeg'] = Tout['Elondeg']*u.deg
    Tout['Latdeg'] = Tout['Latdeg']*u.deg

    SrfLongLat = SkyCoord(-Tout['Elondeg'], Tout['Latdeg'])

    '''Write out table'''
    Tout.write('Data/AresSurfaceViz.ecsv', format='ascii.ecsv', overwrite=True)
    Tout.write('Data/AresSurfaceViz.csv', format='csv', overwrite=True)

    return Tout, SrfLongLat


def makeSurfaceViz(xpdfname, SrfEvts, SrfLongLat):
    ''' plot surface visibility for the transit events '''
    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(111, projection='mollweide')  # mollweide

    '''
    shift center line of geographic projection so that Gale crater and other
    areas with good visibility are at the centre of the plot
    '''
    # x = SrfLongLat.ra.wrap_at(180*u.deg).rad
    x = Angle((SrfLongLat.ra.deg - 180)*u.deg).wrap_at(180*u.deg).rad
    y = SrfLongLat.dec.rad
    SrfSite = SrfEvts['Site'].data

    # colour map
    cm = plt.cm.Reds
    size = SrfEvts['EventDuration'].data * 24  # Julian days to hours

    # markers
    TA = SrfEvts['TerraObs'] == 'Full Visibility'
    LA = SrfEvts['LunaObs'] == 'Full Visibility'
    mask1 = TA & LA
    mask2 = np.logical_not(TA & LA)

    im = ax.scatter(x=x, y=y,
                    c=size, cmap=cm, alpha=0.9)

    # smaller font for partial Viz items
    for i, txt in enumerate(SrfSite[mask2]):
        ax.annotate(txt, (x[mask2][i], y[mask2][i]), alpha=0.7, size=4)

    # larger font for Full Viz items
    for i, txt in enumerate(SrfSite[mask1]):
        ax.annotate(txt, (x[mask1][i]-0.2, y[mask1][i]-0.1),
                    color='red', alpha=0.8, size=8)

    ax.grid(True)
    ax.set_xticklabels(['30°', '60°', '90°', '120°', '150°', '180°',
                        '210°', '240°', '270°', '300°', '330°'])
    plt.figtext(0.5, 0.25,
                'Red labels indicate locations where twin transits are ' +
                'visible in full from object rise till set\n' +
                'Marker colour indicates duration of twin transit in hours',
                wrap=True, horizontalalignment='center')
    plt.figtext(0.5, 0.1, misc.xfooter, size=6,
                horizontalalignment='center')
    plt.title('Preferred Viewing Locations on Mars ' +
              'for Earth & Moon Transit 2018-11-10')
    plt.colorbar(im, orientation='horizontal', shrink=0.5)

    plt.savefig(xpdfname+'.pdf')
    print(xpdfname + '.pdf saved')
    plt.close()


# print('read Event Table')
surf = Table.read('Data/AresSurfaceEventTable.ecsv', format='ascii.ecsv')
# surf.pprint_all()
surf = surf.group_by('site')

SrfEvts, SrfLongLat = surfaceVisibility()
makeSurfaceViz('Output/SurfaceTransitVisibility', SrfEvts, SrfLongLat)
