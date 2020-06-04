"""
    this list of methods contains all relevant convection and diffusion schemes
"""

import numpy as np
import Face

def firstOrderUpwind(i,j,grid,velocity=None,facename="Z"):

    """
        this scheme performs necessary computations for a first-order upwind
        for now it is only really programmed for usage with the convection term
        though usage outside this is unlikely

        inputs:
            i         (int)       -current gridpoint x index
            j         (int)       -current gridpoint y index
            grid      (Grid)      -grid object from Grid.py
            velocity  (ndarray)   -2 element np array containing velocity values
            facename  (string)    -name of face for differencing over (TBLR)
        returns:
            Apn       (ndarray)   -5x1 array of the Ap coefficients around i,j
    """

    # Ap coefficeints we want to finally return
    #   these are ordered LL L P R RR or BB B P T TT
    Apn = np.zeros(5)

    # create a face object
    face = Face.face(i,j,grid,facename,velocity=velocity)

    # set the two Ap coefficients we wish to upstream
    Apn[face.neighbor1Apn] = 1

    #print(str(Apn) + " " + facename)

    return Apn

def secondOrderUpwind(i,j,grid,velocity=None,facename="Z"):

    """
        this scheme performs necessary computations for a second-order upwind
        for now it is only really programmed for usage with the convection term
        though usage outside this is unlikely

        inputs:
            i         (int)       -current gridpoint x index
            j         (int)       -current gridpoint y index
            grid      (Grid)      -grid object from Grid.py
            velocity  (ndarray)   -2 element np array containing velocity values
            facename  (string)    -name of face for differencing over (TBLR)
        returns:
            Apn       (ndarray)   -5x1 array of the Ap coefficients around i,j
    """

    # Ap coefficeints we want to finally return
    #   these are ordered LL L P R RR or BB B P T TT
    Apn = np.zeros(5)

    # create a face object
    face = Face.face(i,j,grid,facename,velocity=velocity)

    neighbor1Size = grid.width\
        (
            face.neighbor1[0],
            face.neighbor1[1],
            face.sizeDirection
        )
    neighbor2Size = grid.width\
        (
            face.neighbor2[0],
            face.neighbor2[1],
            face.sizeDirection
        )

    # create a negative sign if the face normal is against 1,1
    # ensuring usage of the conservation principle
    negative = np.dot(np.array([1,1]),face.normal)
    
    # set the two Ap coefficients we wish to upstream
    Apn[face.neighbor1Apn] = (1 + neighbor1Size / (neighbor1Size + neighbor2Size)) * negative
    Apn[face.neighbor2Apn] = (-neighbor1Size / (neighbor1Size + neighbor2Size)) * negative
     
    # flip all non-P coefficients, as for convection they are put on the RHS
    Apn = Apn * np.array([-1,-1,1,-1,-1])
    print(str(Apn) + " " + facename)

    return Apn

def linear(i,j,grid,facename="Z"):

    """
        this scheme performs necessary computations for a linear scheme
        for now it is only really programmed for usage with the convection term
        it is essentially just calculating linear interpolation weights

        inputs:
            i         (int)       -current gridpoint x index
            j         (int)       -current gridpoint y index
            grid      (Grid)      -grid object from Grid.py
            facename  (string)    -name of face for differencing over (TBLR)
        returns:
            Apn       (ndarray)   -5x1 array of the Ap coefficients around i,j
    """

    # Ap coefficeints we want to finally return
    #   these are ordered LL L P R RR or BB B P T TT
    Apn = np.zeros(5) 

    # create a face object
    face = Face.face(i,j,grid,facename)

    cellSize = grid.width\
        (
            face.cell[0],
            face.cell[1],
            face.sizeDirection
        )
    neighborSize = grid.width\
        (
            face.neighbor1[0],
            face.neighbor1[1],
            face.sizeDirection
        )
    
    # TODO: debug this, it totally doesn't work.

    #negative = np.dot(np.array([1,1]),face.normal)
    negative = 1

    Apn[face.neighbor1Apn] = cellSize     / (cellSize + neighborSize) * negative
    Apn[face.cellApn]      = neighborSize / (cellSize + neighborSize) * negative

    #Apn = Apn * np.array([-1,-1,1,-1,-1])
    print(str(Apn) + " " + facename)

    return Apn

def centralDifference(i,j,grid,facename="Z"):

    """
        this scheme performs necessary computations for a central difference
        for now it is only really programmed for usage with the diffusion term

        inputs:
            i         (int)       -current gridpoint x index
            j         (int)       -current gridpoint y index
            grid      (Grid)      -grid object from Grid.py
            facename  (string)    -name of face for differencing over (TBLR)
        returns:
            Apn       (ndarray)   -5x1 array of the Ap coefficients around i,j
    """

    # Ap coefficeints we want to finally return
    #   these are ordered LL L P R RR or BB B P T TT
    Apn = np.zeros(5)

    # create a face object
    face = Face.face(i,j,grid,facename)

    cellSize = grid.width\
        (
            face.cell[0],
            face.cell[1],
            face.sizeDirection
        )
    neighborSize = grid.width\
        (
            face.neighbor1[0],
            face.neighbor1[1],
            face.sizeDirection
        )

    # TODO: should I allow the area here or put it in outside?
    Apn[face.neighbor1Apn] = 2 * face.area / (cellSize + neighborSize)
    Apn[face.cellApn]      = 2 * face.area / (cellSize + neighborSize)

    print(str(Apn) + " " + facename)

    return Apn


