import numpy as np

class grid:

    def __init__(self, nx, ny, verbose=False):

        #number of points in each direction
        self.nx=nx
        self.ny=ny
        self.ncells=nx*ny

        #grid is stored in this numpy array (val)
        ##order is:
        ### 0 centroid x
        ### 1 centroid y
        ### 2 cell width
        ### 3 cell height
        self.val=np.zeros((nx+4,ny+4,4))
        self.points=np.zeros((nx+1,ny+1,2))

        self.initialize(verbose)

    def initialize(self, verbose):
        #get points on which cells are to be made
        self.points=np.zeros((self.nx+1,self.ny+1,2))
        self.points[:,:,0]=np.linspace(0,1,num=self.nx+1)[:,np.newaxis]
        self.points[:,:,1]=np.linspace(0,1,num=self.nx+1)[np.newaxis,:]

        #internal point setup
        #x centroids are the average of the points grid offset in x
        self.val[2:-2,2:-2,0]=(self.points[1:,:-1,0]+self.points[0:-1,:-1,0])/2
        #y centroids are the average of the points grid offset in y
        self.val[2:-2,2:-2,1]=(self.points[:-1,1:,1]+self.points[:-1,0:-1,1])/2
        #x widths are the difference of the points grid offset in x
        self.val[2:-2,2:-2,2]=self.points[1:,:-1,0]-self.points[0:-1,:-1,0]
        #y widths are the difference of the points grid offset in y
        self.val[2:-2,2:-2,3]=self.points[:-1,1:,1]-self.points[:-1,0:-1,1]
        
        #boundary point setup
        # top and bottom
        #the face x quantities are just copies of their associated center values
        self.val[2:-2,1,0] = self.val[2:-2,2,0] # bottom x centroid
        self.val[2:-2,-2,0] = self.val[2:-2,-3,0] # top x centroid
        self.val[2:-2,1,2] = self.val[2:-2,2,2] # bottom x width #
        self.val[2:-2,-2,2] = self.val[2:-2,-3,2] # top x width #
        self.val[2:-2,1,3] = self.val[2:-2,2,3] # bottom x height #
        self.val[2:-2,-2,3] = self.val[2:-2,-3,3] # top x height #
        #self.val[2:-2,1,3] = 0 # bottom x height
        #self.val[2:-2,-2,3] = 0 # top x height
        # the face y quantities are as their adjacent row +- half the grid width
        self.val[2:-2,1,1] = 0 # bottom y centroid
        self.val[2:-2,-2,1] = 1 # top y centroid
        #note no action taken for the y width, as for the face this is zero

        # left and right
        #the face y quantities are just copies of their associated center values
        self.val[1,2:-2,1] = self.val[2,2:-2,1] # left y centroid
        self.val[-2,2:-2,1] = self.val[-3,2:-2,1] # right y centroid
        self.val[1,2:-2,3] = self.val[2,2:-2,3] # left y height #
        self.val[-2,2:-2,3] = self.val[-3,2:-2,3] # right y height #
        self.val[1,2:-2,2] = self.val[2,2:-2,2] # left y width #
        self.val[-2,2:-2,2] = self.val[-3,2:-2,2] # right y width #
        #self.val[1,2:-2,2] = 0 # left y width
        #self.val[-2,2:-2,2] = 0 # right y width
        # the face x quantities are as their adjacent row +- half the grid width
        self.val[1,2:-2,0] = 0 # left x centroid
        self.val[-2,2:-2,0] = 1 # right x centroid
        #note no action taken for the x width, as for the face this is zero   

        #ghost point setup
        self.val[2:-2,0,2] = self.val[2:-2,2,2] # bottom x width
        self.val[2:-2,-1,2] = self.val[2:-2,-3,2] # top x width
        self.val[2:-2,0,3] = self.val[2:-2,2,3] # bottom x height
        self.val[2:-2,-1,3] = self.val[2:-2,-3,3] # top x height
        self.val[0,2:-2,3] = self.val[2,2:-2,3] # left y height
        self.val[-1,2:-2,3] = self.val[-3,2:-2,3] # right y height
        self.val[0,2:-2,2] = self.val[2,2:-2,2] # left y width
        self.val[-1,2:-2,2] = self.val[-3,2:-2,2] # right y width

        # put a value in the corners
        self.val[1,1,0]=0
        self.val[-2,1,0]=1
        self.val[1,-2,0]=0
        self.val[-2,-2,0]=1
        self.val[1,1,1]=0
        self.val[-2,1,1]=0
        self.val[1,-2,1]=1
        self.val[-2,-2,1]=1

        if verbose:
            #for debugging:
            # these will print the x cetroids followed by y centroids
            # note that transpose is necessary to make sure x values appear as 
            # columns and the flip ensures row 0 (x0 appears at the bottom)
            print("centroid position x")
            print(np.flip(np.transpose(self.val[:,:,0]),0))
            print("\n")
            print("centroid position y")
            print(np.flip(np.transpose(self.val[:,:,1]),0))
            print("\n")
            print("cell width x")
            print(np.flip(np.transpose(self.val[:,:,2]),0))
            print("\n")
            print("cell width y")
            print(np.flip(np.transpose(self.val[:,:,3]),0))

    def center(self,i,j,xy):
        #xy=0 for x, 1 for y
        if xy=="x":
            xy=0
        elif xy=="y":
            xy=1
        else:
            print("ERROR: Invalid centroid select: {}".format(xy))
            return
        return self.val[i,j,xy]

    def width(self,i,j,xy):
        #xy=0 for x, 1 for y
        if xy=="x":
            xy=2
        elif xy=="y":
            xy=3
        else:
            print("ERROR: Invalid width select: {}".format(xy))        
            return
        return self.val[i,j,xy]

    def plot(self, saveImage=False):
        import matplotlib.pyplot as mpl
        mpl.plot(
            self.points[:,:,0],
            self.points[:,:,1],
            'k',
            np.transpose(self.points[:,:,0]),
            np.transpose(self.points[:,:,1]),
            'k'
        )
        mpl.scatter(
            self.val[:,:,0],
            self.val[:,:,1],
            c='r',
            s=2
        )

        if saveImage:
            mpl.savefig("grid.png")
        else:
            mpl.show()
        mpl.close()