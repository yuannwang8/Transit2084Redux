# creating new coordinate class for Mars North frame

import numpy as np
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates.representation import (
    UnitSphericalRepresentation, CartesianRepresentation)
from astropy.coordinates.matrix_utilities import (rotation_matrix,
                                                  matrix_transpose)
from astropy.coordinates import frame_transform_graph


'''rotates vector 45 degrees around y-axis'''

rotation_Ry45 = rotation_matrix(45*u.deg, axis='y')


class RotY45(coord.BaseCoordinateFrame):
    '''
    This class rotates a vector 45 degrees around the y-axis.
    (0,0,1) in (x,y,z) or (0,90) in long, lat should yield
    (-0.707, 0, 0.707) or (180, 45)
    Parameters:
    longitude-like angle
    latitude-like angle
    '''
    default_representation = UnitSphericalRepresentation

    frame_specific_representation_info = {
        UnitSphericalRepresentation: [
            coord.RepresentationMapping('lon', 'lon'),
            coord.RepresentationMapping('lat', 'lat')]
    }


@frame_transform_graph.transform(coord.StaticMatrixTransform,
                                 coord.ICRS, RotY45)
def icrs_to_rottest():
    '''Compute transformation matrix from ICRS frame to RotTest frame'''
    return rotation_Ry45


@frame_transform_graph.transform(coord.StaticMatrixTransform,
                                 RotY45, coord.ICRS)
def rottest_to_icrs():
    '''Compute trasformation matrix from RotTest frame to ICRS frame'''
    return matrix_transpose(rotation_Ry45)


'''realign Mars rotational axis in J2000 frame to point northwards'''


def makeRotationMatrix(carA, carB):
    '''
    This function creates the 3x3 rotational matrix in Cartesian coordinates to
    rotate vectorA around rotational axis unit vector k to new position
    vector B.
    This is based on Rodrigues' rotation forumla.
    np methods are included as reminders.
    '''
    # convert astropy coordinates to numpy array
    # vectorA = [carA.x, carA.y, carA.z]
    # vectorB = [carB.x, carB.y, carB.z]

    # create rotational vector k
    # A = np.array(vectorA)
    # B = np.array(vectorB)
    # k = np.cross(A, B)
    cark = carA.cross(carB)

    # find modulus of k, A and B
    # kmod = np.sqrt(sum(np.square(k)))
    # Amod = np.sqrt(sum(np.square(A)))
    # Bmod = np.sqrt(sum(np.square(B)))

    # find angle theta
    # sintheta = kmod/(Amod*Bmod)
    sintheta = cark.norm()/(carA.norm()*carB.norm())
    theta = np.arcsin(sintheta)

    # find unit vector k
    # khat = k/(Amod*Bmod*sintheta)
    khat = cark/cark.norm()

    # create the cross product matrix for unit kector k
    # K = np.array([[0, -khat[2], khat[1]],
    #              [khat[2], 0, -khat[0]],
    #              [-khat[1], khat[0], 0]])
    K = np.array([[0, -khat.z, khat.y],
                  [khat.z, 0, -khat.x],
                  [-khat.y, khat.x, 0]])

    # create 3x3 rotational matrix R using Rodrigues' rotation formula
    Iden = np.identity(3)
    K2 = np.matmul(K, K)
    R = Iden + sintheta*K + (1-np.cos(theta))*K2
    # print('The rotational matrix to rotate cartesian vector A to B is')
    # print(R)
    return R


J2000R1 = UnitSphericalRepresentation(0*u.deg, 90*u.deg)
carJ1 = J2000R1.represent_as(CartesianRepresentation)

J2000R2 = UnitSphericalRepresentation(10*u.deg, 90*u.deg)
carJ2 = J2000R2.represent_as(CartesianRepresentation)

# Note that J2000R1 and J2000R2 give the same (x,y,z) : obviously!

# photo op event time used: 2084-11-10T07:53:35
Tref = 2451545  # reference JD for '2000-01-01T12:00:00'
Tevent = 2482539.828877315  # photo op event time used: 2084-11-10T07:53:35
MarsRRA = (317.681 - 0.106*(Tevent-Tref)/36525)
MarsRDEC = (52.887 - 0.061*(Tevent-Tref)/36525)

MarsR = UnitSphericalRepresentation(MarsRRA*u.deg, MarsRDEC*u.deg)
carM = MarsR.represent_as(CartesianRepresentation)

rotation_RMarsN = makeRotationMatrix(carM, carJ1)


class MarsNAlign(coord.BaseCoordinateFrame):
    '''
    A unit spherical representation frame that aligns J2000 (RA,DEC) frame to
    Mars's rotation frame: Mars's rotational axis in J2000 should point
    due north (0,0,1) in Cartesian or (0,90) in long, lat.
    Note:
    This new frame is a simple realignment or rotation of the J2000 frame
    to match that of Mars's rotations axis.
    This doesn't take into consideration the location of Mars's prime meridian.
    To fix this, an additional perimeter for the longitude needs
    to be addressed.
    Parameters:
    longitude-like angle
    latitude-like angle
    '''
    default_representation = UnitSphericalRepresentation

    frame_specific_representation_info = {
        UnitSphericalRepresentation: [
            coord.RepresentationMapping('lon', 'lon'),
            coord.RepresentationMapping('lat', 'lat')]
    }


@frame_transform_graph.transform(coord.StaticMatrixTransform,
                                 coord.ICRS, MarsNAlign)
def icrs_to_MarsRotation():
    '''Compute transformation matrix from ICRS frame to Mars Rotation frame'''
    return rotation_RMarsN


@frame_transform_graph.transform(coord.StaticMatrixTransform,
                                 MarsNAlign, coord.ICRS)
def MarsRotation_to_icrs():
    '''Compute trasformation matrix from Mars Rotation frame to ICRS frame'''
    return matrix_transpose(rotation_RMarsN)
