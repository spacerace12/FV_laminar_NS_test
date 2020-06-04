"""
    this class controls the data needed to produce Ap coefficients from a given
    face on a given cell
"""

import numpy as np

class face:

    def __init__(self,i,j,grid,facename,velocity=None):

        """
        the only goal of this function is to produce information
        on where the discretization scheme should apply the modifications to the
        Ap coefficients, NOT to apply the modifications itself.

        It primarily sets the neighbor1, neighbor 2 properties which are the
        indices in the actual grid of the cells the discretization scheme uses

        In addition it sets neighbor1Apn and neighbor2Apn which are the indices
        in the 5x1 grid used for the computational molecule

        neighbor1, in upwinding schemes, is the cell directly upstream of the face
        neighbor2, then, is the cell upwind of neighbor 1

        for central (linear) schemes, neighbor1 is just the cell across the face
        and neighbor2 has no meaning.

        inputs:
            i         (int)       -current gridpoint x index
            j         (int)       -current gridpoint y index
            grid      (Grid)      -grid object from Grid.py
            facename  (string)    -name of face for differencing over (TBLR)
            velocity  (ndarray)   -2 element np array containing velocity values
        """

        # this face's normal vector (always outward)
        self.normal = np.zeros(2)
        # the direction which is upstream of the current cell
        self.upstream = np.zeros(2)
        # the area of this face
        self.area = 0.0
        # the indicies of the cell that this face belongs to
        self.cell = np.zeros(2,dtype=int)
        # the index of the cell in the Ap vector (0-4)
        self.cellApn = 2
        # the indices of the face position (half indexes)
        self.facePosition = np.zeros(2)
        # the indicies of the cell connected to this face (may be the first upwind cell)
        self.neighbor1 = np.zeros(2,dtype=int)
        # the index of the cell in the Apn vector (0-4)
        self.neighbor1Apn = 0
        # the indicies of the cell connected to neighbor 1 (the second upwind cell)
        self.neighbor2 = np.zeros(2,dtype=int)
        # the index of the cell in the Apn vector (0-4)
        self.neighbor2Apn = 0
        # which direction to pull for grid.width to get the cell's width
        self.sizeDirection = "N"
        # the indicies of the cell downwind of the current face
        self.neighborDownwind = np.zeros(2,dtype=int)
        # the index of the downwind cell in the Apn vector (0-4)
        self.neighborDownwindApn = 0

        self.cell = np.array([i,j])

        # determine other cell value indices for each face
        if facename == "T":
            self.area = grid.width(i,j,"x")
            self.normal = np.array([0,1])
            self.facePosition = self.cell + np.array([0,0.5])
            self.neighbor1Apn = 2.5
            self.sizeDirection = "y"
        elif facename == "B":
            self.area = grid.width(i,j,"x")
            self.normal = np.array([0,-1])
            self.facePosition = self.cell + np.array([0,-0.5])
            self.neighbor1Apn = 1.5
            self.sizeDirection = "y"
        elif facename == "L":
            self.area = grid.width(i,j,"y")
            self.normal = np.array([-1,0])
            self.facePosition = self.cell + np.array([-0.5,0])
            self.neighbor1Apn = 1.5
            self.sizeDirection = "x"
        elif facename == "R":
            self.area = grid.width(i,j,"y")
            self.normal = np.array([1,0])
            self.facePosition = self.cell + np.array([0.5,0])
            self.neighbor1Apn = 2.5
            self.sizeDirection = "x"
        else:
            print("ERROR: Not a face name: {}".format(face))

        # if the velocity is not none, overwrite the face normal
        if velocity is None:
            # if the velocity is none, this is not an upwind scheme
            self.upstream = self.normal
            # assign the neighbor indexes
            self.neighbor1 = self.cell + self.normal
            self.neighbor2 = self.cell + self.normal * 2
            # assign the neighbor Ap indexes
            self.neighbor1Apn = self.cellApn + np.sum(self.normal)
            self.neighbor2Apn = self.cellApn + np.sum(self.normal) * 2
        else:
            # get upstream direction by reversing the velocity's sign
            # and constraining to the absolute face normal direction
            upstreamDirection = np.sign(velocity) * -1 * np.abs(self.normal)
            # the first neighbor position is the 1/2 offset in the upstream
            # direction
            self.neighbor1 = (self.facePosition + 0.5 * upstreamDirection).astype(int)
            # the second neighbor position is just one more 
            # upstream from the previous
            self.neighbor2 = (self.neighbor1 + upstreamDirection).astype(int)

            #the upstream direction (1 or -1) for the Apn index
            ApnUpstream = np.dot(upstreamDirection, np.abs(self.normal))
            # the Apn indices are computed in a similar manner
            self.neighbor1Apn = int(self.neighbor1Apn + ApnUpstream * 0.5)
            self.neighbor2Apn = int(self.neighbor1Apn + ApnUpstream)

            # the downwind position is the 1/2 offset in the downstream dir
            self.neighborDownwind = (self.facePosition + -0.5 * upstreamDirection).astype(int)
            self.neighborDownwindApn = int(self.neighbor1Apn + ApnUpstream * -0.5)





        






