    """
        this is a main testing program to do diffusion-only and convection-only
        (and perhaps both combined) tests before moving on to a more full CFD
        method
    """

import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as mpl

def main():
    # x gridpoints
    nx=11
    # y gridpoints
    ny=11
    #number of iterations
    iterations = 30

    kinVisc=1/400
    rho=1.18

    #set relaxation factors
    relaxation = 0.99

    #build grid
    import Grid
    grid=Grid.grid(nx,ny,verbose=False)
    #grid.plot()

    #solution ordered u,v,phi
    solution=np.zeros((nx+4,ny+4,3))

    Ap = np.zeros(np.shape(solution))
    #+1, +1
    solution[2:-2,0:2,2]=1.0
    solution[:,:,0]=1
    solution[:,:,1]=1
    #-1, -1
    # solution[1:-1,-1,2]=1.0
    # solution[:,:,0]=-1
    # solution[:,:,1]=-1

    maxResiduals = np.zeros(3)

    # begin GS relaxation loop
    for it in range(1,iterations+1):
        print("Iteration: {}".format(it))

        olditer=np.copy(solution)

        # loop through non-BC cells
        for i in range(2,grid.nx+2):
            for j in range(2,grid.ny+2):
                #print("\t {0} {1}".format(i,j))
                phi, Ap[i,j,2] = solve(i,j,grid,solution,kinVisc)
                phi = solution[i,j,2] * (1 - relaxation) + relaxation * phi

                # +1, +1
                solution[-2,2:-2,2]=solution[-3,2:-2,2]
                solution[-1,2:-2,2]=solution[-3,2:-2,2]
                solution[2:-2,-2,2]=solution[2:-2,-3,2]
                solution[2:-2,-1,2]=solution[2:-2,-3,2]
                # -1, -1
                # solution[0,1:-1,2]=solution[1,1:-1,2]
                # solution[1:-1,0,2]=solution[1:-1,1,2]
                
                solution[i,j,2] = phi

        maxResiduals = residual(olditer, solution, maxResiduals)

        print("\tResiduals: {0:.6f} {1:.6f} {2:.6f}".format(*maxResiduals))
        print("\tMax: phi:{0:.5f} ".format(solution[5,-5,1]))

        # plot the solution each N iterations
        if it % 1 == 0:
            fig,axs = mpl.subplots(1)
            rangev = np.linspace(0, 1.0, 9, endpoint=True)
            cs = axs.contourf\
                (
                    grid.val[2:-2,2:-2,0], 
                    grid.val[2:-2,2:-2,1], 
                    solution[2:-2,2:-2,2], 
                    rangev, 
                    extend='both', 
                    cmap=mpl.get_cmap('bwr')
                )
            fig.colorbar(cs, ax=axs, ticks=rangev, cmap=mpl.get_cmap('bwr'))
            axs.set_title("phi")
            mpl.scatter(
                grid.val[:,:,0],
                grid.val[:,:,1],
                c='k',
                s=2
            )
            #mpl.plot(grid.val[:,15,0],solution[:,15,2],'k.')
            mpl.savefig("sol_{0}".format(it))
            mpl.close()

    fig,axs = mpl.subplots(1)
    mpl.plot(grid.val[5,2:-2,1],np.flip(np.fliplr(solution[2:-2,2:-2,2]).diagonal()),'k.')
    mpl.savefig("test.png")
    mpl.close()

def solve(i,j,grid,solution,kinVisc):

    from interpolate import linear as linInterp
    from fvSchemes import firstOrderUpwind as FOU
    from fvSchemes import secondOrderUpwind as SOU
    from fvSchemes import centralDifference as CD
    from fvSchemes import linear as linear

    #coefficients
    #   these are ordered LL L P R RR (Apx)
    #   or                BB B P T TT (Apy)
    Apx = np.zeros(5)
    Apy = np.zeros(5)

    #distances to adjacent cell centers times 2
    #   this means dx+dx_L for the left cell
    d_X = grid.width(i,j,'x')
    d_Y = grid.width(i,j,'y')
    d_L = grid.width(i,j,'x') + grid.width(i-1,j,'x')
    d_R = grid.width(i,j,'x') + grid.width(i+1,j,'x')
    d_T = grid.width(i,j,'y') + grid.width(i,j+1,'y')
    d_B = grid.width(i,j,'y') + grid.width(i,j-1,'y')

    # get interpolated velocity * area for convection schemes
    Ft = linInterp(i,j,grid,solution,neighbor="T",field=1) * d_X
    Fb = linInterp(i,j,grid,solution,neighbor="B",field=1) * d_X
    Fl = linInterp(i,j,grid,solution,neighbor="L",field=0) * d_Y
    Fr = linInterp(i,j,grid,solution,neighbor="R",field=0) * d_Y

    Apy += Ft * FOU(i,j,grid,velocity=solution[i,j,0:2],facename="T")
    Apy += Fb * FOU(i,j,grid,velocity=solution[i,j,0:2],facename="B")
    Apx += Fl * FOU(i,j,grid,velocity=solution[i,j,0:2],facename="L")
    Apx += Fr * FOU(i,j,grid,velocity=solution[i,j,0:2],facename="R")

    # Apy += Ft * SOU(i,j,grid,velocity=solution[i,j,0:2],facename="T")
    # Apy += Fb * SOU(i,j,grid,velocity=solution[i,j,0:2],facename="B")
    # Apx += Fl * SOU(i,j,grid,velocity=solution[i,j,0:2],facename="L")
    # Apx += Fr * SOU(i,j,grid,velocity=solution[i,j,0:2],facename="R")

    # Apy += Ft * linear(i,j,grid,facename="T")
    # Apy += Fb * linear(i,j,grid,facename="B")
    # Apx += Fl * linear(i,j,grid,facename="L")
    # Apx += Fr * linear(i,j,grid,facename="R")

    ########## diffusion tests below ##########

    # diffusion: CD
    # note: A_P coefficient is backwards, as it shows up in the RHS originally
    # Apx[2] += 2 * kinVisc * d_Y / d_L + \
    #           2 * kinVisc * d_Y / d_R + \
    #           2 * kinVisc * d_X / d_T + \
    #           2 * kinVisc * d_X / d_B
              
    # Apx[1] += 2 * kinVisc * d_Y / d_L
    # Apx[3] += 2 * kinVisc * d_Y / d_R
    # Apy[1] += 2 * kinVisc * d_X / d_B
    # Apy[3] += 2 * kinVisc * d_X / d_T

    # diffusion: CD with new Face class
    # Apy += CD(i,j,grid,facename="T") * kinVisc
    # Apy += CD(i,j,grid,facename="B") * kinVisc
    # Apx += CD(i,j,grid,facename="L") * kinVisc
    # Apx += CD(i,j,grid,facename="R") * kinVisc

    ########## diffusion tests above ##########

    # compute desired value of Up
    print(str(Apx) + " " + str(i) + " " + str(j))
    print(str(Apy) + " " + str(i) + " " + str(j))
    A_P = Apy[2] + Apx[2]

    # TODO: vectorize after setting Apy[2],Apx[2]=0
    value = \
        (
            Apx[0] * solution[i-2,j,2] + \
            Apx[1] * solution[i-1,j,2] + \
            Apx[3] * solution[i+1,j,2] + \
            Apx[4] * solution[i+2,j,2] + \
            Apy[4] * solution[i,j+2,2] + \
            Apy[3] * solution[i,j+1,2] + \
            Apy[1] * solution[i,j-1,2] + \
            Apy[0] * solution[i,j-2,2]
        ) / A_P

    return value, A_P

def residual(old, new, maxResiduals):
    residual = abs(new - old)
    maxResiduals[0] = np.max(residual[:,:,0])
    maxResiduals[1] = np.max(residual[:,:,1])
    maxResiduals[2] = np.max(residual[:,:,2])
    return maxResiduals

if __name__ == "__main__":
    main()