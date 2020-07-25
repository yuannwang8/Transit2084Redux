# run this as a module from the main working directory
# python3 -m Check.checkMarsNAlign

import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates.representation import (UnitSphericalRepresentation,
                                                CartesianRepresentation)

import T2084Helpers.MarsNAlign as MNA


'''
Check #1: RotY45 frame.
This frame rotates the vector 45 degrees around the y axis.
We expect (0,90) == (0,0,1) rotated 45 degrees around the y axis to yield
(-0.707, 0, 0.707) or (180, 45)
'''
# sph_u = UnitSphericalRepresentation(0*u.deg, 90*u.deg)
# car = sph_u.represent_as(CartesianRepresentation)
# car_rot = car.transform(MNA.rotation_Ry45)
# sph_u_r = car_rot.represent_as(UnitSphericalRepresentation)

C1 = coord.ICRS([0, 90, -180]*u.deg, [90, 0, 45]*u.deg)
C2 = C1.transform_to(MNA.RotY45)

if round(C2[0].lon.deg, 3) == 180.0 and round(C2[0].lat.deg, 3) == 45.0:
    print('\nRotY45 frame: good rotation of (0,90) to (180,45) \n')
else:
    print('\nRotY45 frame: pls check\n')


'''
Check #2: MarsNAlign frame
This frame rotates the J2000 frame to align with Mars's rotational axis, with
North pointing upwards.
The first DEC,RA pair should be a good enough vector for Mars's rotational
axis in the J2000 frame.
'''

vectorR = UnitSphericalRepresentation([317.591, 0, 135, 45]*u.deg,
                                      [52.835, 0, 90, -45]*u.deg)

# carR = vectorR.represent_as(CartesianRepresentation)
# carR_rotated = carR.transform(MNA.rotation_RMarsN)
# carR_rotated_ll = carR_rotated.represent_as(UnitSphericalRepresentation)

C7 = coord.ICRS(vectorR)
C8 = C7.transform_to(MNA.MarsNAlign)
print('Aligning vectors in J2000 frame to Mars North Pole frame\n')
print(C7)
print(C8)

if round(C8[0].dec.deg, 3) == 90:
    print('\nMars rotation axis in J2000 points northwards to 90\n')
else:
    print('\nCheck MarsNAlign frame\n')
