"""
    linearly interpolates a variable between two adjacent cells onto a face
    functional on uniform and non-uniform grids
"""

def linear(i, j, grid, sol, neighbor="Z", field=99):

    """
    inputs:
        i         (int)       -current gridpoint x index
        j         (int)       -current gridpoint y index
        grid      (Grid)      -grid object from Grid.py
        sol       (ndarray)   -numpy array with solution data [x,y,field]
        neighbor  (string)    -name of face to interpolate over (TBLR)
        field     (int)       -field from solution to interpolate
    
    returns:
        value     (float)     -interpolated variable on a face
    """

    #find neighbor to interpolate with
    #find neighbor and this cell's appropriate distance
    if neighbor == "T":
        neighborw = grid.width(i,j+1,'y')
        width = grid.width(i,j,'y')
        neighbori = i
        neighborj = j+1
    elif neighbor == "B":
        neighborw = grid.width(i,j-1,'y')
        width = grid.width(i,j,'y')
        neighbori = i
        neighborj = j-1        
    elif neighbor == "L":
        neighborw = grid.width(i-1,j,'x')
        width = grid.width(i,j,'x')
        neighbori = i-1  
        neighborj = j        
    elif neighbor == "R":
        neighborw = grid.width(i+1,j,'x')
        width = grid.width(i,j,'x')
        neighbori = i+1  
        neighborj = j   
    else:
        print("ERROR: interpolation neighbor not understood: {}".format(neighbor))
        return 0

    if field > 2 or field < 0:
        print("ERROR: field not known: {}".format(field))
        return 0

    value = (
                sol[i,j,field] * neighborw + \
                sol[neighbori,neighborj,field] * width \
            ) / (neighborw + width)
    return value