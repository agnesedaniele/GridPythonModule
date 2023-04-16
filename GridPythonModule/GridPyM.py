#########################################################################
#          V1.1  20-09-22


#          For any question, contact us at agnese.barbensi or 
#          daniele.celoria  followed by @unimelb.edu.au


#########################################################################
#    Created by Agnese Barbensi and Daniele Celoria
#    https://github.com/agnesedaniele
#
#    This module helps with the handling of grid digrams. It is focused on
#    the generation and simplification of grids, rather than on the 
#    computations of invariants. We wish to thank M.Golla and V.Foldvari
#    for helpful comments and suggestions on the code.

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


# Required imports:
from sympy.combinatorics import Permutation
from random import randrange
from matplotlib import pyplot as plt

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def ascending_cusps(input_grid):
    r"""
    Counts the number of ascending cusps in the front. We are following the 
    orientation convention prescribing that O -----> X horizontally.
    
    OUTPUT:
    
    The number of upwards-oriented cusps in the grid.
    
    EXAMPLES::
    
    >> G = load_legendrian_knot('6_2')
    >> draw_grid(G, markings='XO')
    >> descending_cusps(G)
    7
    
    """
    return(descending_cusps(invert_orientation(input_grid)))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def available_knots(cr_number = True):
    r"""
    Writes a list of the names of (knotted) knots whose isotopy classes are pre-loaded. 
    If the variable "cr_number" is specified, it only prints the knots with the given 
    crossing number. In the special case of the unknot, try the "generate_unknot" 
    function instead!

    OUTPUT:

    Prints a list of available knots.

    EXAMPLES::
    
    >> available_knots(5)
    5 crossings:

    5_1  5_2  
    
    """
    if  cr_number < 0 or cr_number > 8:
        raise Exception("We only have links up to 8 crossings!")
    if cr_number not in [True, 3,4,5,6,7,8]:
        raise Exception("Invalid crossing number input")
    list_of_knots = ['3_1', '4_1', '5_1', '5_2', '6_1', '6_2', '6_3', '7_1', '7_2', '7_3', '7_4', '7_5', '7_6', '7_7', '8_1', '8_2', '8_3', '8_4', '8_5', '8_6', '8_7', '8_8', '8_9', '8_10', '8_11', '8_12', '8_13', '8_14', '8_15', '8_16', '8_17', '8_18', '8_19', '8_20', '8_21']
    if cr_number == True or cr_number == 3:
        print('3 crossings:\n',list_of_knots[0],'\n')
    if cr_number == True or cr_number == 4:
        print('4 crossings:\n',list_of_knots[1],'\n')
    if cr_number == True or cr_number == 5:
        print('5 crossings:\n')
        for el in list_of_knots[2:4]:
            print(el,' ', end="")
        print('\n')
    if cr_number == True or cr_number == 6:
        print('6 crossings:\n')
        for el in list_of_knots[4:7]:
            print(el,' ', end="")
        print('\n')
    if cr_number == True or cr_number == 7:
        print('7 crossings:\n')
        for el in list_of_knots[7:14]:
            print(el,' ', end="")
        print('\n')
    if cr_number == True or cr_number == 8:
        print('8 crossings:\n')
        for el in list_of_knots[14:]:
            print(el,' ', end="")
        print('\n')
        
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def available_legendrian_knots(cr_number = True):
    r"""
    Writes a list of the names of (knotted) knots whose Legendrian isotopy classes are 
    pre-loaded. If the variable "cr_number" is specified, it only prints the knots with 
    the given crossing number.

    OUTPUT:

    Prints a list of available Legendrian knots.

    EXAMPLES::
    
    >> available_legendrian_knots(5)
    5 crossings:

    5_1  m(5_1)  5_2  m(5_2)
    
    """
    if  cr_number <= 0 or cr_number > 7:
        raise Exception("We only have nontrivial knots up to 8 crossings!")
    if cr_number not in [True, 3,4,5,6,7]:
        raise Exception("Invalid input")
    list_of_knots = ['3_1', 'm(3_1)', '4_1', '5_1', 'm(5_1)', '5_2', 'm(5_2)', '6_1', 'm(6_1)', '6_2', 'm(6_2)', '6_3', '7_1', 'm(7_1)', '7_2', 'm(7_2)', '7_3', 'm(7_3)', '7_4', 'm(7_4)', '7_5', 'm(7_5)', '7_6', 'm(7_6)', '7_7', 'm(7_7)']
    if cr_number == True or cr_number == 3:
        print('3 crossings:\n')
        for el in list_of_knots[:2]:
            print(el,' ', end="")
        print('\n')
    if cr_number == True or cr_number == 4:
        print('4 crossings:\n')
        for el in list_of_knots[2:3]:
            print(el,' ', end="")
        print('\n')
    if cr_number == True or cr_number == 5:
        print('5 crossings:\n')
        for el in list_of_knots[3:7]:
            print(el,' ', end="")
        print('\n')
    if cr_number == True or cr_number == 6:
        print('6 crossings:\n')
        for el in list_of_knots[7:12]:
            print(el,' ', end="")
        print('\n')
    if cr_number == True or cr_number == 7:
        print('7 crossings:\n')
        for el in list_of_knots[12:]:
            print(el,' ', end="")
        print('\n')

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def check_grid(input_grid):
    r"""
    Checks if the given input represents a valid grid. 

    OUTPUT:

    Return 0 if the grid is valid, and 1 otherwise.

    EXAMPLES::
    
    >> G = [[0,1,2,3,4],[2,3,4,0,1]]
    >> check_grid(G)
    0
    >> G = [[0,1,2,3,4],[0,3,4,0,1]]
    >> check_grid(G)
    1
    
    """
    if len(input_grid) != 2:
        return 1
    A,B = input_grid
    if len(A) != len(B) or set(A) != set(B) or set(A) != set([i for i in range(len(A))]) or set(B) != set([i for i in range(len(A))]):
        return 1
    for j in range(len(A)):
        if A[j] == B[j]:
            return 1
    return 0

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def coherent_bs(input_grid, where, which = 'rows'):
    r"""
    Performs a coherent band attachment on the grid. Can choose between a row/column
    band attachment (using the input 'which'), and specifying the number of the 
    row/column using the 'where' input. The syntax for the 'which' input is either
    'rows'/'r' or 'columns'/'c'.
    
    OUTPUT:
    
    The grid obtained after a coherent band attachment.
    
    EXAMPLES::
    
    >> G = load_knot('3_1')
    >> B_band = coherent_bs(G, 2, 'c')
    >> C
    [[4, 0, 1, 3, 2], [1, 2, 3, 4, 0]]
    
    """
    if which not in ['columns','rows','r','c']:
        raise Exception("Invalid input, the only options are 'rows', 'r', 'columns', 'c'.")
    if check_grid(input_grid) == 1:
        raise Exception("Invalid Input grid")
    if (which == 'r' or which == 'rows') and _check_rows_for_bands(input_grid,where)[0] != 0: 
        raise Exception("The selected rows/columns are not a band attachment site")   
    if (which == 'c' or which == 'columns') and _check_rows_for_bands(rotate(input_grid,1),where)[0] != 0:
        raise Exception("The selected rows/columns are not a band attachment site")
    A = input_grid[0]
    B = input_grid[1]
    if which == 'rows' or which == 'r':    
        remember = [A[where], A[where+1]]
        A[where] = remember[1]
        A[where+1] = remember[0]
        return(A,B)
    if which == 'columns' or which == 'c':
        AA,BB = rotate([A,B],1)
        remember = [AA[where], AA[where+1]]
        AA[where] = remember[1]
        AA[where+1] = remember[0]
        return(rotate([AA,BB],3))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def commute_rows(input_grid, where, interleaving = 'N', verbose = False):
    r"""
    Performs a commutation between two adjacent rows. The "where" parameter (between 
    0 and the dimension of the grid-1) specifies which rows we want to commute.
    If the "interleaving" parameter can be 'Y','N','A'. In the first case we perform the
    commutation only if the markers are interleaving (this produces a crossing change on
    the link type); if it is 'N' only a non-interleaving commutation is performed (this
    corresponds to a planar isotopy or a Reidemeister 2 move). If it is 'A', the commutation
    is always performed (unless the adjacent markings at the same height).

    OUTPUT:

    A grid where the two selected rows are commuted, according to the interleaving
    parameter.

    EXAMPLES::
    
    >> A = [[2, 3, 0, 1, 4], [3, 4, 1, 0, 2]]
    >> B = commute_columns(A,1)
    >> B
    [[1, 3, 0, 2, 4], [3, 4, 2, 0, 1]]
    
    """
    if check_grid(input_grid) == 1:
        raise Exception("Invalid Input grid")
    if where > len(input_grid[0]) or where < 0 or interleaving not in ['A','Y','N']:
        raise Exception("Invalid parameters")
    A = input_grid[0]
    B = input_grid[1]
    if _check_rows(input_grid, where) == 5:
        return 0
    if interleaving == 'A':
        return([A[:where] + [A[where+1]] + [A[where]] + A[where+2:],B[:where] + [B[where+1]] + [B[where]] + B[where+2:]])
    if interleaving == 'Y':
        if _check_rows(input_grid, where) == False:
            return([A[:where] + [A[where+1]] + [A[where]] + A[where+2:],B[:where] + [B[where+1]] + [B[where]] + B[where+2:]])
        else:
            if verbose == True:
                print("These are not interleaved!")
            return 0
    if interleaving == 'N':
        if _check_rows(input_grid, where) == True:
            return([A[:where] + [A[where+1]] + [A[where]] + A[where+2:],B[:where] + [B[where+1]] + [B[where]] + B[where+2:]])
        else:
            if verbose == True:
                print("These rows/columns are interleaved!")
            return 0

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def commute_columns(input_grid,where, interleaving = 'A', verbose = False):
    r"""
    Same as "commute_rows", but on columns.

    OUTPUT:

    A grid where the two selected columns are commuted, according to the interleaving
    parameter.
    
    """
    if check_grid(input_grid) == 1:
        raise Exception("Invalid Input")
    A = input_grid[0]
    B = input_grid[1]
    if commute_rows(rotate([A,B],1),where, interleaving = interleaving, verbose = verbose) == 0:
        return 0
    AA,BB = commute_rows(rotate([A,B],1),where, interleaving = interleaving, verbose = verbose)
    return(rotate([AA,BB],3))    

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def connected_sum(input_grid1,input_grid2):
    r"""
    Creates a grid representing the connected sum of the two inputs. In the case of links,
    it performs the connected sum between the "upper" component of the first grid and the 
    "lower" one in the second (i.e. those components that have a strand in the top/bottom
    rows of the grids).

    OUTPUT:

    A grid for the connected sum.
    
    EXAMPLES::
    
    >> G = generate_torus_link(3,2)
    connected_sum(G, mirror_grid(G))
    [[0, 1, 2, 3, 4, 5, 9, 8, 7, 6], [3, 5, 0, 1, 2, 8, 7, 6, 4, 9]]
    
    """
    if check_grid(input_grid1) == 1 or check_grid(input_grid2) == 1:
        raise Exception("Invalid Input")
    n1 = len(input_grid1[0])
    n2 = len(input_grid2[0])
    A1,B1 = cyclic_shift(input_grid1, horizontal = input_grid1[0][n1-1]+1, vertical = 0)
    A2,B2 = cyclic_shift(input_grid2, horizontal = input_grid2[0][0], vertical= 0)
    Bsum1,Bsum2 = [],[]    
    for i in range(n1):
        if B1[i] != n1 -1:
            Bsum1.append(B1[i])
        else:
            Bsum1.append(B1[i]+1)
    for i in range(n2):
        if B2[i] != 0:
            Bsum2.append(B2[i]+n1)
        else:
            Bsum2.append(B2[i]+n1 -1)
    return([A1+[a+n1 for a in A2], Bsum1 + Bsum2])

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def convert_to_braid(input_grid, optimized = 'Y'):
    r"""
    Return a braid word whose closure represents the same link as the input grid.  
    The "optimized" option can be either "Y", "N", "S". The first gives as output the 
    minimal braid among two possible choices (obtained by suitably rotating the grid),
    the second just output the first braid obtained, while the last option first attempts
    to simplify the grid before producing the smallest braid.
    
    OUTPUT:

    A string of integers representing the product of generators in a braid group.

    EXAMPLES::
    
    >> G = generate_random_grid(8)
    >> print(G)
    [[6, 0, 3, 5, 7, 2, 4, 1], [5, 1, 7, 3, 6, 4, 2, 0]]

    >> convert_to_braid(G)
    array([ 3,  2, -2])
    
    >> F = generate_random_grid(8)
    >> convert_to_braid(F)
    array([-1, -2, -3, -4, -5, -6, -7, -8])
    
    """
    if check_grid(input_grid) == 1:
        raise Exception("Invalid Input")
    if optimized not in ['N','Y','S']:
        raise Exception("Invalid optimization input.")
    if optimized == 'N':
        return _aux_braid(input_grid)
    elif optimized == 'Y':
        firstB = _aux_braid(input_grid)
        secondB = [-a for a in _aux_braid(rotate(input_grid,1))]
        if firstB == [] or secondB == []:
            print('Trivial braid representing the link with %s components'%number_of_components(input_grid))
            return [1]
        if len(firstB) > len(secondB):
            return secondB
        else:
            return firstB
    elif optimized == 'S':
        aux_G = simplify_grid(input_grid)
        firstB = _aux_braid(aux_G)
        secondB = [-a for a in _aux_braid(rotate(input_grid,1))]
        if firstB == [] or secondB == []:
            print('Trivial braid representing the link with %s components'%number_of_components(input_grid))
            return [1]
        if len(firstB) > len(secondB):
            return secondB
        else:
            return firstB

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def convert_to_Sage(input_grid):
    r"""
    It converts a grid diagram into a "Link" type in Sage. For multi-component links
    there might be issues with this function, related to how Sage deals with closures
    of braids.
    
    OUTPUT:

    A Sage Link type object 

    EXAMPLES::
    >> G = load_knot('3_1')
    >> K = convert_to_Sage(G)
    >> K.khovanov_homology()
    {-9: {-3: Z},
     -7: {-3: 0, -2: C2},
     -5: {-3: 0, -2: Z, -1: 0, 0: 0},
     -3: {-3: 0, -2: 0, -1: 0, 0: Z},
     -1: {0: Z}}
     
    """    
    braid_word = convert_to_braid(input_grid)
    br_ind = max([abs(c) for c in braid_word]) + 1
    BG = BraidGroup(br_ind)
    return Link(BG(braid_word))
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def crossing_number(input_grid):
    r"""
    Computes the crossing number of the input grid (note that this is NOT the minimal
    crossing number of the isotopy class!).
    
    OUTPUT:

    A non negative integer.

    EXAMPLES::
    
    >> G = generate_random_grid(8)
    >> crossing_number(G)
    5
    >> G = parallel_copies(generate_torus_link(3,2),11)
    >> crossing_number(G)
    363
    
    """
    A = input_grid[0]
    B = input_grid[1]
    crossings_aux = 0
    for i in range(1,len(A)-1):
        for j in range(min(A[i],B[i]),max(A[i],B[i])):
            if min(A.index(j), B.index(j)) < i < max(A.index(j), B.index(j)):
                crossings_aux += 1
    return crossings_aux

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def cyclic_shift(input_grid, horizontal = 1, vertical = 1):
    r"""
    Performs a cyclic shift (also known as cyclic or torus permutations).
    Takes in input a grid diagram, possibly with the amount of horizontal (to the left) 
    and vertical (downwards) shifts. The defaults values of the shifts are set to 1. 
    
    OUTPUT:

    The shift of the given grid.

    EXAMPLES::
    
    >> G = generate_random_grid(8)
    >> cyclic_shift(G, 3,6)
    [[5, 4, 1, 3, 7, 0, 2, 6], [4, 7, 6, 1, 2, 3, 0, 5]]
    >> cyclic_shift(G)
    [[7, 4, 0, 6, 1, 5, 2, 3], [3, 1, 5, 7, 2, 6, 4, 0]]
    
    """
    if check_grid(input_grid) == 1 or horizontal <0 or vertical <0:
        raise Exception("Invalid Input")
    A = input_grid[0]
    B = input_grid[1]
    if vertical != 0:
        A = A[vertical:]+A[:vertical]
        B = B[vertical:]+B[:vertical]
    if horizontal != 0:
        A = rotate_once(A)[horizontal:] + rotate_once(A)[:horizontal]
        B = rotate_once(B)[horizontal:] + rotate_once(B)[:horizontal]
        A,B = rotate([A,B],3)
    return([A,B])

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def descending_cusps(input_grid):
    r"""
    Counts the number of descending cusps in the front. We are following the 
    orientation convention prescribing that O -----> X horizontally.
    
    OUTPUT:
    
    The number of downwards-oriented cusps in the grid.
    
    EXAMPLES::
    
    >> G = load_legendrian_knot('6_2')
    >> draw_grid(G, markings='XO')
    >> descending_cusps(G)
    7
    
    """
    if check_grid(input_grid) == 1:
        raise Exception("Invalid Input")
    Dcusps = 0
    A = input_grid[0]
    B = input_grid[1]
    for i in range(len(A)):
        if A[i] < B[i]:
            if B.index(A[i]) < i :
                Dcusps += 1
            if A.index(B[i]) > i:
                Dcusps += 1
    return(Dcusps)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def destabilize_all(input_grid):
    r"""
    Performs all possible (generalized)destabilizations on the grid at once.
   
    OUTPUT:

    A grid which is not the stabilisation of another grid.

    EXAMPLES::
   
    >> G = generate_unknot(20)
    >> destabilize_all(G)
   
    [[0, 1], [1, 0]]
   
    """
    A = input_grid[0]
    B = input_grid[1]
    forbidden = 0 
    nn = grid_number(input_grid)
    if nn == 2:
        return input_grid
    while _can_simplify([A,B]) and forbidden ==0:
        if _check_distance_one([A, B]) != []:
            horiz = _check_distance_one([A, B])
            A, B,forbidden = _destabilize_aux([A, B], horiz[0])
        elif _check_distance_one([B, A]) != []:
            horiz = _check_distance_one([B, A])
            B, A,forbidden = _destabilize_aux([B, A], horiz[0])
        elif _check_distance_one(rotate([A,B],1)) != []:
            vert = _check_distance_one(rotate([A,B],1))
            A, B,forbidden = _destabilize_aux(rotate([A,B],1), vert[0])
            A, B = rotate([A,B],3)
        elif _check_distance_one(_invert(rotate([A,B],1))) != []:
            vert = _check_distance_one(_invert(rotate([A,B],1)))
            B, A,forbidden = _destabilize_aux(_invert(rotate([A,B],1)), vert[0])
            B, A = rotate([B,A],3)            
    return([A, B])                

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def destabilize(input_grid, where, selection = 'row', verbose = False):
    r'''
    Performs a single destabilization move on the grid, if possible at the location 
    specified by the 'where' and 'selection' parameters. Takes as input a grid,
    an integer between 0 and the grid's dimension - 2, and the choice of 'row' or
    'column'. 
    Note that in general the number of destabilizations by row and column do not
    necessarily coincide!
    
    OUTPUT: A grid of dimension one less than the input, or 0.
   
    EXAMPLES::
   
    >> G = [[6, 5, 0, 1, 4, 2, 3], [5, 1, 2, 4, 3, 6, 0]]
    >> destabilize(G,0,selection = 'row')
    [[5, 0, 1, 4, 2, 3], [1, 2, 4, 3, 5, 0]]
    
    '''
    A,B = input_grid
    nn = grid_number(input_grid)
    if where > nn-1 or where < 0:
        raise Exception("The index must be between 0 and grid number")
    if nn == 2:
        if verbose == False:
            return 0
        else:
            raise Exception("Grids of size 2 cannot be destabilized")
    if where == nn-1:
        new_out = destabilize(rotate(input_grid,2),0,selection = selection)
        if new_out != 0:
            return rotate(new_out,2)
        else:
            return 0
    if where > nn -1:
        if verbose == False:
             return 0
        else:
            raise Exception("Input index must be less or equal then grid size")
    if selection == 'row' or selection == 'rows':
        if A[where] ==  B[where+1] and A[where] == B[where]+1:
            A,B,check = _destabilize_aux(input_grid,where)
            return [A,B]
        elif A[where] ==  B[where+1] and A[where] == B[where]-1:
            A,B,check = _destabilize_aux(input_grid,where)
            return [A,B]
        elif B[where] ==  A[where+1] and  B[where] == A[where]+1:
            B,A,check = _destabilize_aux([B, A],where)
            return [A,B]
        elif B[where] ==  A[where+1] and  B[where] == A[where]-1:
            B,A,check = _destabilize_aux([B, A],where)
            return [A,B]
        elif A[where] == B[where-1] and A[where] == B[where] +1 and where != 0:
            new_grid = destabilize(rotate([A,B],3),  nn-A[where]-1, selection = 'row')
            A,B = rotate(new_grid,1)
            return [A,B]
        elif A[where] == B[where-1] and A[where] == B[where] -1 and where != 0:
            new_grid = destabilize(rotate([A,B],3),  nn-A[where]-1, selection = 'row')
            A,B = rotate(new_grid,1)
            return [A,B]
        elif B[where] == A[where-1] and B[where] == A[where] +1 and where != 0:
            B,A = destabilize([B,A], where, selection = 'row')
            return [A,B]
        elif B[where] == A[where-1] and B[where] == A[where] -1 and where != 0:
            B,A = destabilize([B,A], where, selection = 'row')
            return [A,B]
        else:  
            if verbose == False:
                return 0
            else:
                raise Exception("There is no destabilization with the given inputs")
    elif selection == 'col' or selection == 'column' or selection == 'columns':
        new_grid = destabilize(rotate([A,B],3),  nn-where-1 , selection = 'row')
        if new_grid == 0:
            if verbose == False:
                return 0
            else:
                raise Exception("There is no destabilization with the given inputs")
        A,B = rotate(new_grid,1)
        return [A,B]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def disjoint_union(input_grid1,input_grid2):
    r"""
    Creates a grid representing the disjoint union of the links represented by the 
    two inputs.
    
    OUTPUT:

    A grid representing G1 U G2.

    EXAMPLES::
    
    >> G = generate_torus_link(3,2)
    >> disjoint_union(G,G)
    [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], [3, 4, 0, 1, 2, 8, 9, 5, 6, 7]]
    
    """
    if check_grid(input_grid1) == 1 or check_grid(input_grid2) == 1:
        raise Exception("Invalid Input")
    A1,B1 = input_grid1
    A2,B2 = input_grid2
    n = len(A1)
    return([A1 + [t + n for t in A2], B1 + [t + n for t in B2]])

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def draw_grid(input_grid, color = 0, markings = False, save_fig = False, fig_name = 'grid_diagram'):
    r"""
    Draws the link associated with the input grid diagram. Crossings are to be resolved 
    with the vertical line as an overcrossing. 
    Options include the color of the grid; the default value is 0, which produces a 
    rainbow grid. Other possible options include standard colors like "red", "green", etc.
    It is possible to draw the position of the markings that determine the grid; the two
    available options are 'XO' (which draws them as small Xs and Os) and 'dots'. If this 
    option is left unselected, the markings are not drawn.
    The 'save_fig' option allows to save the drawing in the current folder, with a name
    specified by the 'fig_name' option.

    OUTPUT:

    A png image of the link represented by the grid.
    
    EXAMPLES:
    
    >> G = generate_random_grid(10)
    >> draw_grid(G, color = 'red')
    >> draw_grid(G, color = 'red', markings = 'dots', save_fig = True)
    Picture saved in the current folder.

    """
    if check_grid(input_grid) == 1:
        raise Exception("Invalid grid Input!")
    if markings not in ['XO',False, 'dots']:
        raise Exception("Invalid markings option.")
    A = input_grid[0]
    B = input_grid[1]
    n = len(A)
    for j in range(n+1):
        plt.plot([j,j ], [0, n], 'k-', lw=0.1)
        plt.plot([0,n], [j, j], 'k-', lw=0.1)
        if j !=n and color != 0:
            plt.plot([B[j] +1/2, A[j] +1/2],[ j +1/2, j +1/2], color = color)
            plt.plot([A[j] +1/2,A[j] +1/2], [ j+1/2,j + 1/2 + _distance_markings(A,B,j)],color = color)
        elif j !=n and color == 0:
            plt.plot([B[j] +1/2, A[j] +1/2],[ j +1/2, j +1/2])
            plt.plot( [A[j] +1/2,A[j] +1/2], [ j+1/2,j + 1/2 + _distance_markings(A,B,j)])
            if markings == 'XO':
                plt.plot([A[i] + 0.5 for i in range(n)], [i + 0.5 for i in range(n)], 'x')
                plt.plot([B[i] + 0.5 for i in range(n)],[i + 0.5 for i in range(n)],  'o')
            if markings == 'dots':
                plt.plot([A[i] + 0.5 for i in range(n)], [i + 0.5 for i in range(n)], '.')
                plt.plot([B[i] + 0.5 for i in range(n)],[i + 0.5 for i in range(n)],  '.')
    plt.axis('off')
    plt.axis('equal')
    if save_fig == True:
        plt.savefig("%s.png" %fig_name)
        print('Picture saved in the current folder.')
    plt.show()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def Gauss_code(input_grid, verbose = False):
    r"""
    Determines the Gauss code of the given grid diagram. At the moment it only works
    for knots.

    OUTPUT: 
    
    A list containing two sublists describing the Gauss code for the input grid.

    EXAMPLES::
    
    >> G = generate_torus_link(3,4)
    >> print(Gauss_code(G))
    [[[-1, -2, 3, 4, -5, -6, 2, 7, -4, 5, -8, 1, -7, -3, 6, 8]], [1, 1, 1, 1, 1, 1, 1, 1]]
    
    """
    if check_grid(input_grid) == 1:
        raise Exception("Invalid grid Input!")
    if number_of_components(input_grid) != 1:
        raise Exception("Right now the Gauss code is only implemented for knots!")
    empty_cpts=0
    A = input_grid[0]
    B = input_grid[1]
    GC = []
    signs = [] 
    c = 1
    coord = [] 
    rows = [i for i in range(len(A))]
    while len(rows) > 1:
        gc = []
        where = min(rows)
        flag = 0
        start = where
        while flag ==0:
            valueB = B[where]
            valueA = A[where]
            if valueB < valueA:
                for kk in range(valueB+1,valueA):
                    v1 = min(A.index(kk),B.index(kk))
                    v2 = max(A.index(kk),B.index(kk))
                    if v1 < where < v2:
                        pos = [where,kk]
                        if pos not in coord:
                            coord.append(pos)
                            gc.append(-c)
                            c = c+1
                            if B.index(kk)> A.index(kk):
                                signs.append(-1)
                            if B.index(kk)< A.index(kk):
                                signs.append(1)
                        else:
                            caux = coord.index(pos)+1
                            gc.append(-caux)
            if valueB > valueA:
                for kk in range(valueB-1,valueA,-1):
                    v1 = min(A.index(kk),B.index(kk))
                    v2 = max(A.index(kk),B.index(kk))
                    if v1 < where < v2:
                        pos = [where,kk]
                        if pos not in coord:
                            coord.append(pos)
                            gc.append(-c)
                            c = c+1
                            if B.index(kk)> A.index(kk):
                                signs.append(1)
                            if B.index(kk)< A.index(kk):
                                signs.append(-1)                        
                        else:
                            caux = coord.index(pos)+1
                            gc.append(-caux)
            if B.index(valueA) > where:
                for kk in range(where+1,B.index(valueA)):
                    if min(A[kk],B[kk]) < valueA < max(A[kk],B[kk]):
                        pos = [kk,valueA]
                        if pos not in coord:
                            coord.append(pos)
                            gc.append(c)
                            c = c+1
                            if B[kk]> A[kk]:
                                signs.append(1)
                            if B[kk]< A[kk]:
                                signs.append(-1)
                        else:   
                            caux = coord.index(pos)+1
                            gc.append(caux)
            if B.index(valueA) < where:
                for kk in range(where-1,B.index(valueA),-1):
                    if min(A[kk],B[kk]) < valueA < max(A[kk],B[kk]):
                        pos = [kk,valueA]
                        if pos not in coord:
                            coord.append(pos)
                            gc.append(c)
                            c = c+1
                            if B[kk]> A[kk]:
                                signs.append(-1)
                            if B[kk]< A[kk]:
                                signs.append(1)
                        else:   
                            caux = coord.index(pos)+1
                            gc.append(caux)
            rows.remove(where)                
            where = B.index(valueA)
            if where ==start:
                flag =1
        if len(gc)>0:        
            GC.append(gc)
        else:
            empty_cpts += 1 
    if empty_cpts==1 and verbose == True:
        print('There is 1 unlinked trivial component')
    if empty_cpts>1 and verbose == True:
        print('There are %s unlinked trivial components' %empty_cpts)
    return([GC,signs])    

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def generate_random_grid(grid_size, components = 0):
    r"""
    Generates a random grid diagram of the specified dimension. If the "components"
    input is set to 0 (the default value), then no check on the number of components
    is performed. If instead the "components" variable is set to a positive integer,
    then the output grid will represent a link with the given number of components 
    (whenever this is possible). The generation is done by randomly sampling two
    collision-free permutations on n elements; so the grid is truly random, but we 
    make no claim of ergodicity of this generation method.
    Furthermore, for large grid sizes, specifying a low number of components might
    result in long computational times.

    OUTPUT:

    A random grid diagram representing a link with a specified number of components.

    EXAMPLES::
    
    >> G = generate_random_grid(10)
    >> print(G)
    [[1, 7, 2, 6, 0, 5, 9, 3, 4, 8], [4, 5, 0, 7, 6, 3, 1, 8, 9, 2]]
    >> F = generate_random_grid(6, components = 2)
    >> print(F)
    [[2, 1, 5, 4, 0, 3], [3, 0, 4, 5, 2, 1]]
    >> number_of_components(F)
    2
    
    """
    if components < 0 or components > grid_size or grid_size <= 2:
        raise Exception("Wrong number of components or incorrect grid size!")
    iters = 0
    while iters < 10000:
        A = Permutation.random(grid_size).list()
        B = Permutation.random(grid_size).list()
        if check_grid([A,B]) == 0:
            if components == 0:
                 return [A,B]
            else:
                if number_of_components([A,B]) == components:
                      return [A,B]
        iters +=1
    raise Exception("something went wrong in the function! This might be due to the absence of grids with the given number of components in the input dimension.")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def generate_torus_link(P,Q):
    r"""
    Creates a grid representing the torus link T(P,Q). The two parameters need to be
    non-zero integers, and the number of components of the link is gcd(P,Q).

    OUTPUT:

    A grid diagram for T(P,Q).

    EXAMPLES::
    
    >> G = generate_torus_link(5,15)
    >> number_of_components(G)
    5
    
    """
    if P == 0 or Q == 0:
        raise Exception("Only non-zero coefficients")
    AA = [r for r in range(abs(P)+abs(Q))]
    BB = [r for r in range(abs(P), abs(P)+abs(Q))] + [r for r in range(abs(P))]
    if P*Q < 0:
        return [AA,BB]
    else:
        return mirror_grid([AA,BB])

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def generate_unknot(grid_size):
    r"""
    Generates a specific (staircase-shaped) unknot diagram of the specified dimension.  

    OUTPUT:

    A grid diagram representing the unknot.

    EXAMPLES::
    
    >> G = generate_unknot(10)
    >> G
    [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], [1, 2, 3, 4, 5, 6, 7, 8, 9, 0]]
    
    """
    if grid_size <=2:
        raise Exception("The grid size needs to be an integer greater than 2!")
    return([[r for r in range(grid_size)],[r + 1 for r in range(grid_size-1)]+[0]])

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def generate_unlink(link_components):
    r"""
    Generates a specific unlink diagram with the specified number of components.  

    OUTPUT:

    A grid diagram representing the unlink.

    EXAMPLES::
    
    >> G = generate_unlink(5)
    >> number_of_components(G)
    5
    
    """
    if link_components < 1:
        raise Exception("The number of components needs to be positive!")
    second_aux = []
    for ii in range(2*link_components):
        if ii%2 == 0:
            second_aux.append((ii+1)%(2*link_components))
        else:
            second_aux.append(ii-1)
    return [[ii for ii in range(2*link_components)], second_aux]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def generate_twist_knot(number_of_twists, clasp = '+'):
    r"""
    Generates a grid diagram representing a twist knot with the given number or twists.
    the sign of the clasp is determined by a '+' or '-' input. If no clasp type is
    selected in the input, a predermined one is assigned.

    OUTPUT:

    A grid diagram representing the given twist knot.

    EXAMPLES::
    
    >> F = generate_twist_knot(5,'-')
    >> F
    [[0, 5, 3, 2, 4, 7, 6, 9, 8, 1], [4, 2, 1, 0, 6, 5, 8, 7, 3, 9]]
    
    """
    if clasp not in ['+','-']:
        raise Exception("Invalid clasp type")
    if number_of_twists == 0:
        raise Exception("There are easier ways of creating unknots!")
    twists = abs(number_of_twists)
    n = 4 + twists
    auxA, auxB = [],[]
    flag1,flag2,flag3,flag4 = 1,1,1,1
    if clasp == '+' and number_of_twists > 0:
        flag1 = 0
    if clasp == '-' and number_of_twists < 0:
        flag2 = 0
    if clasp == '+' and number_of_twists < 0:
        flag3 = 0
    if clasp == '-' and number_of_twists > 0:
        flag4 = 0
    if flag1 == 0 or flag2 == 0:
        auxB.append(0)
        if twists%2 == 0:
            for i in range(int(twists/2) + 1):
                auxA.append(2*i+3)
                auxA.append(2*i+2)
            for i in range(int(twists/2)):
                auxB.append(2*i+4)
                auxB.append(2*i+3)
            auxA.append(0)
            auxA.append(1)
            auxB.append(1)
            auxB.append(2)
            auxB.append(n-1)
        else:
            for i in range(int((twists+1)/2)):
                auxA.append(2*i+3)
                auxA.append(2*i+2)
                auxB.append(2*i+4)
                auxB.append(2*i+3)
            auxA.append(1)
            auxA.append(0)
            auxA.append(n-1)
            auxB.append(2)
            auxB.append(1)
        if flag1 == 0:    
            return([auxA,auxB])
        else:
            return(mirror_grid([auxA,auxB]))
    else:
        if twists%2 == 0:
            auxA.append(4+twists)
            for i in range(int((n-4)/2)-1):
                auxA.append(n - 3 - 2*i )
                auxA.append(n - 2 - 2*i )
                auxB.append(n - 2 - 2*i)
                auxB.append(n - 1 - 2*i)
            for ii in [1,4,n-1, 3,2,0]:
                auxA.append(ii)
            for ii in [4,5,0,2,1,n,3]:
                auxB.append(ii)
        else:
            auxB.append(n)
            for i in range(int((n-3)/2)-1):            
                auxB.append(n-3 -2*i)
                auxB.append(n-2 -2*i)
                auxA.append(n-2*i - 2)
                auxA.append(n-2*i - 1)
            for ii in [0,n-1,1,2,3]:
                auxB.append(ii)
            for ii in [1,4,2,3,n, 0]:
                auxA.append(ii)
        if flag3 == 0:
            return([auxA,auxB])    
        else:
            return(mirror_grid([auxA,auxB]))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def grid_length(input_grid):
    r"""
    Computes the length of the grid. This is computed by summing the length of all the
    segments that compose the grid.
    
    OUTPUT:

    Returns the length of the input grid.

    EXAMPLES::
    
    >> G = generate_twist_knot(5)
    >> grid_length(G)
    56
    
    """
    A = input_grid[0]
    B = input_grid[1]
    Ar, Br = rotate([A,B],1)
    return sum([abs(A[i] - B[i]) for i in range(len(A))]) + sum([abs(Ar[i] - Br[i]) for i in range(len(A))])

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def grid_number(input_grid):
    r"""
    The dimension of the input grid. 
    
    OUTPUT:

    The grid number of the given grid.

    EXAMPLES::
    
    >> G = generate_twist_knot(5)
    >> grid_number(G)
    9
    
    """
    if check_grid(input_grid) == 1:
        raise Exception("Invalid Input")
    return(len(input_grid[0]))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def invert_orientation(input_grid):
    r"""
    Provides a grid representing the same link with the orientation of all of its
    components inverted. 
    
    OUTPUT:

    Returns the reverse of the grid.

    EXAMPLES::
    
    >> G = generate_random_grid(7, components=2)
    >> print(G)
    [[2, 5, 4, 3, 1, 6, 0], [5, 1, 6, 2, 0, 4, 3]]

    >> invert_orientation(G)
    [[5, 1, 6, 2, 0, 4, 3], [2, 5, 4, 3, 1, 6, 0]]
    
    """
    return([input_grid[1],input_grid[0]])

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def load_knot(knot_name, verbose = False):
    r"""
    Produces a minimal grid representing of a given knot type. 
    The pre-loaded knots are for (minimal) crossing numbers less than 9;
    To see which knots are available, use the command "available_knots()". 
    Note that the knot name needs to be between quotes, e.g. '7_2'.
    If the "verbose" option is True, a link to the associated Knotatlas and Knotinfo 
    pages are printed.
    
    OUTPUT:

    Returns the minimal grid representatives of a given knot type.

    EXAMPLES::
    
    >> load_knot('5_2')
    [[1, 5, 6, 2, 3, 4, 0], [6, 0, 4, 5, 1, 2, 3]]
    
    >> load_knot('7_3', verbose = True)
    Check out some of this knot's invariants at the Knot Atlas:
    http://katlas.org/wiki/7_3 
    and Knotinfo:
    https://knotinfo.math.indiana.edu/results.php?searchmode=singleknot&desktopmode=0&mobilemode=0&singleknotprev=&submittype=singleknot&singleknot=7_3

    [[3, 0, 8, 7, 6, 5, 2, 1, 4], [8, 7, 6, 5, 4, 1, 0, 3, 2]]
    
    """
    list_of_knots = {'3_1': [[4, 0, 1, 2, 3], [1, 2, 3, 4, 0]], '4_1': [[5, 0, 3, 4, 2, 1], [2, 4, 5, 1, 0, 3]], '5_1':[[4,5, 6, 0, 1, 2, 3], [6, 2, 3, 4, 5, 0, 1]], '5_2': [[1, 5, 6, 2, 3, 4, 0], [6, 0, 4, 5, 1, 2, 3]],'6_1': [[7, 0, 5, 6, 3, 4, 2, 1], [2, 6, 7, 4, 5, 1, 0, 3]], '6_2': [[7, 2, 4, 5, 6, 1, 0, 3], [5,6, 7, 0, 3, 4, 2, 1]], '6_3': [[7, 0, 3, 5, 6, 4, 2, 1], [3, 4, 6, 7, 2, 1, 0, 5]], '7_1': [[1, 6,7, 8, 0, 2, 3, 4, 5], [8, 0, 1, 4, 5, 6, 7, 2, 3]], '7_2': [[1, 7, 8, 5, 6, 2, 3, 4, 0], [8, 0, 6,7,4, 5, 1, 2, 3]], '7_3': [[3, 0, 8, 7, 6, 5, 2, 1, 4], [8, 7, 6, 5, 4, 1, 0, 3, 2]], '7_4': [[5, 0,8, 2, 1, 7, 4, 3, 6], [8, 7, 1, 0, 6, 3, 2, 5, 4]], '7_5': [[2, 7, 8, 0, 3, 4, 5, 6, 1], [8, 0, 1,6,7, 2, 3, 4, 5]], '7_6': [[1, 7, 8, 3, 2, 4, 5, 6, 0], [8, 0, 2, 1, 6, 7, 3, 4, 5]], '7_7': [[1, 4,7, 8, 6, 3, 2, 5, 0], [8, 0, 2, 5, 1, 7, 4, 3, 6]], '8_1': [[9, 5, 6, 3, 4, 1, 2, 8, 7, 0], [6, 7,4,5, 2, 3, 0, 1, 9, 8]], '8_2': [[9, 0, 8, 1, 7, 2, 3, 4, 5, 6], [1, 7, 2, 3, 9, 4, 5, 6, 8, 0]],'8_3': [[4, 7, 6, 9, 8, 2, 3, 0, 1, 5], [6, 5, 8, 7, 3, 4, 1, 2, 9, 0]], '8_4': [[1, 0, 3, 2, 8,9,4, 5, 6, 7], [9, 2, 1, 4, 3, 5, 6, 7, 8, 0]], '8_5': [[5, 0, 10, 9, 7, 8, 4, 3, 2, 1, 6],[10, 9,8,3, 2, 6, 7, 1, 5, 4, 0]], '8_6': [[9, 3, 2, 0, 7, 8, 1, 4, 5, 6], [2, 1, 8, 3, 9, 4, 5, 6, 7, 0]],'8_7': [[9, 0, 6, 7, 8, 5, 4, 3, 2, 1], [5, 7, 8, 9, 4, 3, 2, 1, 0, 6]], '8_8': [[9, 0, 3, 7, 2,1,5, 4, 8, 6], [1, 7, 8, 9, 4, 3, 2, 6, 5, 0]], '8_9': [[1, 0, 8, 9, 4, 3, 2, 5, 6, 7], [9, 4, 3, 5,2, 1, 6, 7, 8, 0]], '8_10': [[8, 2, 5, 6, 7, 4, 3, 1, 0, 10, 9], [3, 6, 7, 8, 10, 9, 5, 4, 2, 1, 0]],'8_11': [[9, 0, 8, 4, 7, 3, 5, 1, 2, 6], [4, 7, 5, 6, 9, 8, 2, 3, 0, 1]], '8_12': [[2, 3, 0, 1, 5,6, 9, 8, 4, 7], [8, 1, 2, 6, 7, 4, 5, 3, 9, 0]], '8_13': [[9, 0, 5, 7, 8, 4, 2, 3, 1, 6], [4, 7, 8,9, 3, 1, 0, 6, 5, 2]], '8_14': [[9, 2, 8, 3, 7, 4, 0, 1, 5,6], [3, 7, 4, 5, 9, 1, 2, 6, 8, 0]], '8_15': [[10, 1, 6, 8, 7, 2, 3, 4, 5, 9, 0],  [2, 8, 9, 10, 3,6,0, 7, 1, 4, 5]], '8_16': [[2, 1, 7, 8, 3, 5, 4, 9, 0, 10, 6],  [9, 5, 10, 6, 7, 8, 2, 3, 4, 1, 0]],'8_17': [[6, 8, 9, 2, 3, 5, 4, 1, 0, 7], [9, 1, 3, 4, 7, 8, 0, 6, 5, 2]], '8_18': [[2, 1, 3, 4, 0,6, 7, 9, 8, 5], [9, 7, 8, 2, 3, 1, 5, 6, 4, 0]], '8_19': [[2, 1, 0, 6, 5, 4, 3], [6, 5, 4, 3, 2, 1,0]], '8_20': [[7, 1, 0, 3, 4, 6, 5, 2], [3, 6, 4, 5, 7, 2, 1, 0]], '8_21': [[7, 0, 3, 5, 4, 6, 1, 2],[4, 6, 7, 1, 0, 2, 3, 5]] }
    if knot_name not in list_of_knots.keys():
        raise Exception("Invalid input name! Try the command 'available_knots' to see which knots are pre-loaded.")
    if verbose == True:
        print("\nCheck out some of this knot's invariants at the Knot Atlas:\nhttp://katlas.org/wiki/%s \nand Knotinfo:\n https://knotinfo.math.indiana.edu/results.php?searchmode=singleknot&desktopmode=0&mobilemode=0&singleknotprev=&submittype=singleknot&singleknot=%s" %(knot_name,knot_name))
    return list_of_knots[knot_name]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def load_legendrian_knot(knot_name, random_repr = True):
    r"""
    Produces a grid representing a non-destabilizable Legendrian knot type 
    (the XO coordinates are taken from https://services.math.duke.edu/~ng/atlas/). 
    Note that the knot name needs to be between quotes, e.g. 'm(7_7)'.
    The pre-loaded knots are for crossing number less than 8 (and thus arc index 
    less than 10); in case of non-simple knots, the output consists of all the 
    maximal representatives. If instead the "random_repr" option is True, it gives 
    a random representative among the list.
    
    To see which knots are available, use the command "available_legendrian_knots()".
    
    OUTPUT:

    Returns a Legendrian representatives of a given knot type.

    EXAMPLES::
    
    >> load_legendrian_knot('7_2')
    [[1, 7, 8, 5, 6, 2, 3, 4, 0], [8, 0, 6, 7, 4, 5, 1, 2, 3]]
    
    >> load_legendrian_knot('7_2', random_repr = False)
    [[[3, 0, 8, 7, 6, 5, 2, 1, 4], [8, 7, 6, 5, 4, 1, 0, 3, 2]],
    [[3, 1, 8, 7, 6, 5, 0, 4, 2], [8, 7, 6, 5, 4, 2, 3, 1, 0]]]
    
    """
    list_of_legendrian = {'3_1': [[[4, 0, 1, 2, 3], [1, 2, 3, 4, 0]]], 'm(3_1)': [[[4, 3, 2, 1, 0], [1, 0, 4, 3, 2]]], '4_1': [[[5, 0, 3, 4, 2, 1], [2, 4, 5, 1, 0, 3]]], '5_1': [[[4, 5, 6, 0, 1, 2, 3], [6, 2, 3, 4, 5, 0, 1]], [[1, 2, 3, 4, 5, 6, 0], [6, 0, 1, 2, 3, 4, 5]]], 'm(5_1)': [[[1, 0, 6, 5, 4, 3, 2], [6, 5, 4, 3, 2, 1, 0]]], '5_2': [[[1, 5, 6, 2, 3, 4, 0], [6, 0, 4, 5, 1, 2, 3]]], 'm(5_2)': [[[3, 0, 6, 5, 2, 1, 4], [6, 5, 4, 1, 0, 3, 2]], [[3, 1, 6, 5, 0, 4, 2], [6, 5, 4, 2, 3, 1, 0]]], '6_1': [[[7, 0, 5, 6, 3, 4, 2, 1], [2, 6, 7, 4, 5, 1, 0, 3]]], 'm(6_1)': [[[7, 0, 5, 6, 2, 1, 4, 3], [4, 6, 7, 1, 0, 3, 2, 5]], [[7, 2, 5, 6, 4, 1, 0, 3], [4, 6, 7, 3, 0, 5, 2, 1]]], '6_2': [[[7, 2, 4, 5, 6, 1, 0, 3], [5, 6, 7, 0, 3, 4, 2, 1]], [[7, 0, 6, 4, 1, 2, 3, 5], [1, 4, 2, 7, 3, 5, 6, 0]], [[7, 6, 4, 5, 0, 1, 2, 3], [5, 0, 7, 1, 2, 3, 4, 6]]], 'm(6_2)': [[[7, 6, 5, 4, 1, 0, 2, 3], [2, 0, 7, 6, 5, 3, 4, 1]]], '6_3': [[[7, 0, 3, 5, 6, 4, 2, 1], [3, 4, 6, 7, 2, 1, 0, 5]], [[7, 6, 5, 3, 4, 0, 1, 2], [5, 4, 0, 7, 1, 2, 3, 6]]], '7_1': [[[1, 6, 7, 8, 0, 2, 3, 4, 5], [8, 0, 1, 4, 5, 6, 7, 2, 3]], [[6, 7, 8, 0, 1, 2, 3, 4, 5], [8, 4, 5, 6, 7, 0, 1, 2, 3]], [[1, 2, 3, 4, 5, 6, 7, 8, 0], [8, 0, 1, 2, 3, 4, 5, 6, 7]]], 'm(7_1)': [[[1, 0, 8, 7, 6, 5, 4, 3, 2], [8, 7, 6, 5, 4, 3, 2, 1, 0]]], '7_2': [[[1, 7, 8, 5, 6, 2, 3, 4, 0], [8, 0, 6, 7, 4, 5, 1, 2, 3]]], 'm(7_2)': [[[5, 0, 8, 7, 2, 1, 4, 3, 6], [8, 7, 6, 1, 0, 3, 2, 5, 4]], [[5, 2, 8, 7, 0, 1, 6, 3, 4], [8, 7, 6, 4, 5, 3, 2, 0, 1]], [[5, 3, 8, 7, 0, 6, 2, 1, 4], [8, 7, 6, 4, 5, 1, 0, 3, 2]], [[5, 3, 8, 7, 1, 6, 0, 4, 2], [8, 7, 6, 4, 5, 2, 3, 1, 0]]], '7_3': [[[3, 0, 8, 7, 6, 5, 2, 1, 4], [8, 7, 6, 5, 4, 1, 0, 3, 2]], [[3, 1, 8, 7, 6, 5, 0, 4, 2], [8, 7, 6, 5, 4, 2, 3, 1, 0]]], 'm(7_3)': [[[4, 7, 8, 5, 6, 0, 1, 2, 3], [8, 2, 6, 7, 3, 4, 5, 0, 1]], [[1, 7, 8, 2, 3, 4, 5, 6, 0], [8, 0, 6, 7, 1, 2, 3, 4, 5]]], '7_4': [[[5, 0, 8, 2, 1, 7, 4, 3, 6], [8, 7, 1, 0, 6, 3, 2, 5, 4]], [[5, 2, 8, 0, 7, 4, 3, 1, 6], [8, 7, 3, 6, 1, 2, 0, 5, 4]], [[6, 3, 8, 0, 7, 5, 2, 1, 4], [8, 7, 5, 6, 4, 1, 0, 3, 2]], [[6, 3, 8, 1, 7, 5, 0, 4, 2], [8, 7, 5, 6, 4, 2, 3, 1, 0]]], 'm(7_4)': [[[1, 7, 8, 2, 5, 6, 3, 4, 0], [8, 0, 6, 7, 1, 4, 5, 2, 3]]], '7_5': [[[2, 7, 8, 0, 3, 4, 5, 6, 1], [8, 0, 1, 6, 7, 2, 3, 4, 5]], [[1, 7, 8, 2, 4, 3, 5, 6, 0], [8, 0, 4, 6, 7, 1, 2, 3, 5]], [[4, 7, 8, 3, 5, 6, 0, 1, 2], [8, 1, 6, 7, 2, 3, 4, 5, 0]], [[1, 6, 7, 8, 2, 3, 4, 5, 0], [8, 0, 5, 6, 7, 1, 2, 3, 4]]], 'm(7_5)': [[[3, 0, 8, 7, 6, 2, 1, 5, 4], [8, 7, 6, 5, 1, 0, 4, 3, 2]], [[4, 1, 8, 7, 6, 0, 5, 3, 2], [8, 7, 6, 5, 3, 4, 2, 1, 0]]], '7_6': [[[1, 7, 8, 3, 2, 4, 5, 6, 0], [8, 0, 2, 1, 6, 7, 3, 4, 5]], [[5, 7, 8, 0, 1, 3, 2, 6, 4], [8, 2, 6, 4, 7, 0, 5, 3, 1]], [[5, 7, 8, 4, 0, 6, 1, 2, 3], [8, 2, 6, 7, 5, 3, 4, 0, 1]]], 'm(7_6)': [[[2, 0, 8, 7, 5, 6, 4, 3, 1], [8, 7, 6, 1, 0, 3, 2, 5, 4]], [[3, 7, 8, 2, 1, 6, 5, 0, 4], [8, 0, 6, 7, 5, 4, 2, 3, 1]], [[3, 0, 8, 7, 5, 6, 2, 1, 4], [8, 7, 6, 1, 0, 4, 5, 3, 2]]], '7_7': [[[1, 4, 7, 8, 6, 3, 2, 5, 0], [8, 0, 2, 5, 1, 7, 4, 3, 6]], [[1, 5, 7, 8, 4, 3, 2, 6, 0], [8, 0, 2, 6, 7, 1, 5, 4, 3]], [[2, 0, 8, 4, 3, 6, 5, 7, 1], [8, 7, 1, 0, 5, 4, 2, 3, 6]]], 'm(7_7)': [[[3, 7, 8, 1, 0, 6, 4, 5, 2], [8, 0, 6, 7, 5, 2, 1, 3, 4]], [[3, 7, 8, 2, 5, 6, 1, 0, 4], [8, 1, 6, 7, 0, 4, 5, 3, 2]]]}
    if knot_name not in list_of_legendrian.keys():
        raise Exception("Invalid input name! Try the command 'available_knots' to see which knots are pre-loaded.")
    if random_repr == True:
        n = len(list_of_legendrian[knot_name])
        return list_of_legendrian[knot_name][randrange(0,n)]
    if len(list_of_legendrian[knot_name]) == 1:
        return list_of_legendrian[knot_name][0]
    return list_of_legendrian[knot_name]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def mirror_grid(input_grid):
    r"""
    Provides a grid whose associated knot type is the mirror of the given one. 
    
    OUTPUT:

    Returns the mirror of the grid.

    EXAMPLES::
    
    >> G = [[0,1,2,3,4],[2,3,4,0,1]]
    >> mirror_grid(G)
    [[4, 3, 2, 1, 0], [1, 0, 4, 3, 2]]
    
    """
    return(rotate(input_grid, 1))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def number_of_components(input_grid): 
    r"""
    Gives the number of components of the link represented by the input grid diagram.  

    OUTPUT:

    An integer.

    EXAMPLES::
    
    >> G = generate_unknot(15)
    >> print(number_of_components(G))
    1
    >> F = generate_random_grid(6)
    >> F
    [[4, 5, 1, 2, 0, 3], [1, 2, 0, 5, 3, 4]]
    
    print(number_of_components(F))
    2
    
    """
    A = input_grid[0]
    B = input_grid[1]
    VAR = ((Permutation(A)**(-1))*Permutation(B)).cycle_structure
    return sum([VAR[r] for r in VAR if r != 1])

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def parallel_copies(input_grid, strands):
    r"""
    Creates a grid representing the link obtained by taking parallel copies of the link. 
    
    OUTPUT:

    Returns a grid.

    EXAMPLES::
    
    >> G = generate_twist_knot(5)
    number_of_components(parallel_copies(G,7))
    7
        
    """
    if check_grid(input_grid) == 1:
        raise Exception("Invalid Input")
    if strands <0:
        raise Exception("Invalid number of strands")
    A,B = input_grid
    Aout,Bout = [],[]
    for i in range(len(A)):
        if A[i] < B[i] and B.index(A[i]) < i:
            Aout += [strands*(A[i]+1)-t-1 for t in range(strands)]
        elif A[i]> B[i] and B.index(A[i]) < i:
            Aout += [strands*A[i]+t for t in range(strands)]
        elif A[i] < B[i] and B.index(A[i]) > i:
            Aout += [strands*A[i]+t for t in range(strands)]
        elif A[i] > B[i] and B.index(A[i]) > i:
            Aout += [strands*(A[i]+1)-t-1 for t in range(strands)]
        if B[i] < A[i] and A.index(B[i]) < i:
            Bout += [strands*(B[i]+1)-t-1 for t in range(strands)]
        elif B[i]> A[i] and A.index(B[i]) < i:
            Bout += [strands*B[i]+t for t in range(strands)]
        elif B[i] < A[i] and A.index(B[i]) > i:
            Bout += [strands*B[i]+t for t in range(strands)]
        elif B[i] > A[i] and A.index(B[i]) > i:
            Bout += [strands*(B[i]+1)-t-1 for t in range(strands)]
    return([Aout,Bout])

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def perform_all_moves(input_grid):
    r"""
    Performs all possible Cromwell moves on the given grid (not including cyclic shifts).
    In other words, produces all grid diagrams in the radius 1 ball around the input grid
    inside the Reidemeister graph.
        
    OUTPUT:

    Returns a list of all the grid obtained from the input via a single move.

    EXAMPLES::
    
    >> G = [[0,1],[1,0]]
    >> perform_all_moves(G)
    [[[1, 0, 2], [0, 2, 1]],
     [[0, 1, 2], [1, 2, 0]],
     [[0, 1, 2], [2, 0, 1]],
     [[1, 0, 2], [2, 1, 0]],
     [[1, 0, 2], [2, 1, 0]],
     [[2, 0, 1], [1, 2, 0]],
     [[0, 1, 2], [1, 2, 0]],
     [[0, 2, 1], [2, 1, 0]],
     [[0, 2, 1], [2, 1, 0]],
     [[0, 1, 2], [1, 2, 0]],
     [[0, 1, 2], [2, 0, 1]],
     [[0, 2, 1], [1, 0, 2]],
     [[1, 0, 2], [2, 1, 0]],
     [[0, 1, 2], [2, 0, 1]],
     [[1, 2, 0], [2, 0, 1]],
     [[0, 2, 1], [2, 1, 0]]]
     
    """
    if check_grid(input_grid) == 1:
        raise Exception("Invalid Input")
    nn = grid_number(input_grid)
    out_grids = []
    for i in range(nn):
        for stab_type in ['XNE','XNW','XSE','XSW','ONE','ONW','OSE','OSW']:
            out_aux = stabilisation(input_grid, row = i, kind = stab_type)
            if out_aux not in out_grids:
                out_grids.append(out_aux)
    for i in range(nn-1):
        out_comm_cols = commute_columns(input_grid, where = i, interleaving = 'N')
        if out_comm_cols != 0 and out_comm_cols not in out_grids:
            out_grids.append(out_comm_cols)
        out_comm_rows = commute_rows(input_grid, where = i, interleaving = 'N')
        if out_comm_rows != 0 and out_comm_rows not in out_grids:
            out_grids.append(out_comm_rows)
    for i in range(nn):
        out_destab_row = destabilize(input_grid, i, selection = 'row')
        if out_destab_row != 0 and out_destab_row not in out_grids:
            out_grids.append(out_destab_row)
        out_destab_col = destabilize(input_grid, i, selection = 'col')
        if out_destab_col != 0 and out_destab_col not in out_grids:
            out_grids.append(out_destab_col)
    return out_grids

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def rotate(input_grid, number_rotations):
    r"""
    Rotates the grid number_rotations times counter-clockwise. Note that rotating
    by pi/2 returns a grid representing the mirror of the input.
    
    OUTPUT:

    Returns the rotated grid.

    EXAMPLES::
    
    >> G = [[0,1,2,3,4],[2,3,4,0,1]]
    >> rotate(G,1)
    [[4, 3, 2, 1, 0], [1, 0, 4, 3, 2]]
    
    
    >> G = [[0,1,2,3,4],[2,3,4,0,1]]
    >> rotate(G,4)
    [[0,1,2,3,4],[2,3,4,0,1]]
    
    """
    A = input_grid[0]
    B = input_grid[1]
    for i in range(number_rotations):
        A = rotate_once(A)
        B = rotate_once(B)
    return([A, B])
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def rotate_once(input_list):
    r"""
    Auxiliary function that rotates one set of markings once counter-clockwise. 
    
    OUTPUT:

    Returns a list of the rotated markings.

    EXAMPLES::
    
    >> G = [0,1,2,3,4]
    >> rotate_once(G)
    [4, 3, 2, 1, 0]
    
    """
    return [len(input_list)-input_list.index(a)-1 for a in range(len(input_list))] 

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def rotation_number(input_grid):
    r"""
    Computes the rotation number of the Legendrian knot represented by the grid.
    
    OUTPUT:

    The rotation number of the grid as an integer.

    EXAMPLES::
    
    >> G = stabilisation(load_legendrian_knot('4_1'),3,'XNW')
    >> print(rotation_number(G))
    1
    
    """
    if check_grid(input_grid) == 1:
        raise Exception("Invalid Input")
    return int((descending_cusps(input_grid) - ascending_cusps(input_grid))/2)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def scramble_grid(input_grid, effort = 'medium', legendrian = False, transverse = False):
    r"""
    Creates a grid representing the same knot type as the input, but complicated using 
    random Cromwell moves. 
    The "effort" variable can be equal to 'very_low','low', 'medium', 'high' and 'very_high'.
    These choices correspond to the number of random moves applied to the grid (not including 
    destabilisations), which are 20, 50, 100, 1000 or 10000 respectively. 
    It is also possible to specify a custom integer value for the effort.
    If the legendrian or transverse options are set to True, then only random legendrian
    or transverse moves are performed.
    
    OUTPUT:

    A new (more complex) grid in the same isotopy/legendrian/transverse isotopy class as the 
    input.

    EXAMPLES::
    >> G = generate_torus_link(2,3)
    >> scramble_grid(G, effort = 'very_low')
    [[2, 4, 5, 6, 7, 1, 8, 3, 0], [7, 5, 6, 8, 1, 0, 3, 2, 4]]    
    
    >> G = generate_torus_link(2,3)
    >> F = scramble_grid(G, legendrian = True)
    >> [thurston_bennequin(G),rotation_number(G)] == [thurston_bennequin(F) ,rotation_number(F)]
    True
    
    >> G = generate_twist_knot(3, clasp='+')
    >> S = NEW_scramble(G, effort= 42)
    >> grid_number(S)
    14
    
    """
    if check_grid(input_grid) == 1:
        raise Exception("Invalid Input")
    if effort not in ['low', 'high', 'medium', 'very_low','very_high']:
        if type(effort) != type(1):
            raise Exception("Invalid effort!")
        else:
            if effort < 1:
                raise Exception("Invalid effort!")
    if transverse == True and legendrian == True:
        raise Exception("Only one of the two options (legendrian/transverse) at a time")
    if effort == 'very_low':
        tries = 20
    if effort == 'low':
        tries = 50
    if effort == 'medium':
        tries = 100
    if effort == 'high':
        tries = 1000    
    if effort == 'very_high':
        tries = 10000
    elif type(effort) == type(1):
        tries = effort
    A,B = input_grid
    n = len(A)
    count = 0
    while count < tries:
        dice = randrange(0,5)
        n = grid_number([A,B])
        if dice == 0:
            A,B = cyclic_shift([A,B], randrange(0,n+1), randrange(0,n+1))
            count += 1
        elif dice  == 1:
            possible_rows = []
            for iterator in range(n-1):   #why is it < n-1 here? 
                if _check_rows([A,B],iterator) == True:
                    possible_rows.append(iterator)
            if len(possible_rows) > 0:
                A,B = commute_rows([A,B],possible_rows[randrange(0,len(possible_rows))])
                count += 1
        elif dice == 2:
            possible_cols = []
            for iterator in range(n-1):   
                if _check_columns([A,B],iterator) == True:
                    possible_cols.append(iterator)
            if len(possible_cols) > 0:
                A,B = commute_columns([A,B],possible_cols[randrange(0,len(possible_cols))])
                count += 1
        elif dice == 3:
            if legendrian == False and transverse == False:
                possible_stab = ['XNE','XNW','XSE','XSW','ONE','ONW','OSE','OSW']
                A,B = stabilisation([A,B], row = randrange(0,n), kind = possible_stab[randrange(0,8)])
            elif legendrian == True:
                possible_stab = ['XNE','XSW','OSW','ONE']
                A,B = stabilisation([A,B], row = randrange(0,n), kind = possible_stab[randrange(0,4)])
            elif transverse == True:
                possible_stab = ['XNE','XSW','XSE','OSW','ONE','ONW']
                A,B = stabilisation([A,B], row = randrange(0,n), kind = possible_stab[randrange(0,6)])
            count += 1
        elif dice == 4:
            if legendrian == False and transverse == False:
                possible_stabs_row = []
                possible_stabs_col = []
                for iterator in range(n-1):
                    if destabilize([A,B], iterator, selection = 'row', verbose = False) != 0:
                        possible_stabs_row.append(iterator)
                    if destabilize([A,B], iterator, selection = 'col', verbose = False) != 0:
                        possible_stabs_col.append(iterator)
                if len(possible_stabs_row)*len(possible_stabs_col) != 0:
                    coin = randrange(0,2)
                    if coin == 0:
                        A,B = destabilize([A,B], possible_stabs_row[randrange(len(possible_stabs_row))],selection = 'row', verbose = False)
                        count += 1
                    else:
                        A,B = destabilize([A,B], possible_stabs_col[randrange(len(possible_stabs_col))], selection = 'col', verbose = False)
                        count += 1
    return [A,B]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def self_linking(input_grid):
    r"""
    Computes the self-linking number of the transverse knot represented by the grid.
    In the case of multi-component links, it computes the sum over all components.
    
    OUTPUT:

    The self-linking number of the grid.

    EXAMPLES::
    
    >> G = load_legendrian_knot('m(5_2)')
    >> print(self_linking(G))
    1
    
    """
    return(writhe(input_grid) - descending_cusps(input_grid))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def simplify_grid(input_grid, effort = 'medium', verbose = False):
    r"""
    Simplifies the given grid diagram, using an algorithm similar to Gridlink's built-in
    function (see http://homepages.math.uic.edu/~culler/gridlink/). 
    There are two options: the "effort" variable can be equal to 'low', 'medium' or 'high'.
    These choices correspond to the number of random moves applied to the grid (not including 
    destabilizations), which are 200, 1000 or 10000 respectively. 
    If the "verbose" option is True (which is the default value) the amount of simplification
    is printed.
    
    OUTPUT:

    A simplified grid, representing the same knot type.

    EXAMPLES::
    
    >> G = generate_unknot(100)
    >> simplify(G)
    
    Grid simplification from 100 to 2
    [[0, 1], [1, 0]]
    
    >> G = generate_random_grid(20)
    >> print(G)
    [[15, 16, 4, 13, 7, 18, 0, 3, 6, 12, 19, 2, 11, 14, 1, 8, 17, 10, 5, 9], 
    [1, 7, 9, 15, 12, 2, 19, 18, 17, 3, 4, 14, 13, 0, 5, 11, 6, 16, 10, 8]]
    >> simplify(G)
    
    Grid simplification from 20 to 4
    [[3, 0, 1, 2], [1, 2, 3, 0]]
    
    """
    if check_grid(input_grid) == 1:
        raise Exception("Invalid Input")
    count = 0
    A,B = input_grid
    if effort not in ['low', 'high', 'medium']:
        if type(effort) != type(1):
            raise Exception("Invalid effort!")
        else:
            if effort < 1:
                raise Exception("Invalid effort!")
    elif effort == 'low':
        tries = 200
    elif effort == 'medium':
        tries = 1000
    elif effort == 'high':
        tries = 10000    
    else:
        tries = effort
    while count < tries and len(A)>3:
        A,B = destabilize_all([A,B])
        n = len(A)
        dice = randrange(0,3)
        if dice == 0:
            D = cyclic_shift([A,B], randrange(0,n+1), randrange(0,n+1))
            A = D[0]
            B = D[1]
            count += 1
        elif dice ==1:
            done = 0
            iterator = 0
            while done ==0 and iterator<n-1:
                if _check_rows([A,B],iterator) == True:
                    A,B = commute_rows([A,B],iterator)
                    done =1
                    count += 1
                else:
                    iterator += 1
        elif dice ==2:
            done = 0
            iterator = 0
            while done ==0 and iterator<n-1:
                if _check_columns([A,B],iterator) == True:
                    A,B = commute_columns([A,B],iterator)
                    done =1    
                    count += 1
                else:
                    iterator += 1
    if verbose == True:
        print('Grid simplification from %s to %s' %(len(input_grid[0]), len(A)))
    return([A,B])                

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def stabilisation(input_grid, row = 3, kind = 'XNE'):
    r"""
    Performs a single stabilisation on the input grid. The input parameters include the 
    row where the stabilisation is performed, as well as the type of stabilisation.
    These are 'XNE','XNW','XSE','XSW','ONE','ONW','OSE','OSW', where the first letter
    specifies which markings on the row to select, and the other two the position of the
    empty square left after the stabilisation.
    
    OUTPUT:

    A once-stabilised grid digram.

    EXAMPLES::
    
    >> G = generate_torus_link(4,3)
    >> print(G)
    [[0, 1, 2, 3, 4, 5, 6], [4, 5, 6, 0, 1, 2, 3]]
    >> stabilisation(G, 2, kind = 'XNW')
    [[0, 1, 2, 3, 4, 5, 6, 7], [5, 6, 3, 7, 0, 1, 2, 4]]
    
    """
    if check_grid(input_grid) == 1:
        raise Exception("Invalid Input")
    if kind not in ['XNE','XNW','XSE','XSW','ONE','ONW','OSE','OSW'] or 0 > row  or  row > len(input_grid[0]) :
        raise Exception("Invalid kind of stabilisation!")
    n = len(input_grid[0])
    if kind[0] == 'X':
        Oflag = False
        A,B = input_grid
    else:
        Oflag = True
        A,B = invert_orientation(input_grid)
        kind = 'X'+ kind[1:]
    Astab = []
    Bstab = []
    if kind == 'XNW':
        for i in range(n):
            if i == row:
                Astab.append(A[row])
                Astab.append(A[row] +1)
            elif A[i] > A[row]:
                Astab.append(A[i]+1)
            else:
                Astab.append(A[i])
        for i in range(n):
            if i == row:
                Bstab.append(A[row] +1)
                if B[row] < A[row]:
                    Bstab.append(B[row])
                else:
                    Bstab.append(B[row] +1)
            elif B[i] > A[row]:
                Bstab.append(B[i]+1)
            else:
                Bstab.append(B[i])
    elif kind == 'XNE':
        for i in range(n):
            if i == row:
                Astab.append(A[row]+1)
                Astab.append(A[row])
            elif A[i] > A[row]:
                Astab.append(A[i]+1)
            else:
                Astab.append(A[i])
        for i in range(n):
            if i == row:
                Bstab.append(A[row])
                if B[row] < A[row]:
                    Bstab.append(B[row])
                else:
                    Bstab.append(B[row] +1)
            elif B[i] >= A[row]:
                Bstab.append(B[i]+1)
            else:
                Bstab.append(B[i])
    elif kind == 'XSW':
        for i in range(n):
            if i == row:
                Astab.append(A[row]+1)
                Astab.append(A[row])
            elif A[i] > A[row]:
                Astab.append(A[i]+1)
            else:
                Astab.append(A[i])
        for i in range(n):
            if i == row:
                if B[row] < A[row]:
                    Bstab.append(B[row])
                else:
                    Bstab.append(B[row] +1)
                Bstab.append(A[row]+1)
            elif B[i] > A[row]:
                Bstab.append(B[i]+1)
            else:
                Bstab.append(B[i])
    elif kind == 'XSE':
        for i in range(n):
            if i == row:
                Astab.append(A[row])
                Astab.append(A[row]+1)
            elif A[i] > A[row]:
                Astab.append(A[i]+1)
            else:
                Astab.append(A[i])
        for i in range(n):
            if i == row:
                if B[row] < A[row]:
                    Bstab.append(B[row])
                else:
                    Bstab.append(B[row] +1)
                Bstab.append(A[row])
            elif B[i] >= A[row]:
                Bstab.append(B[i]+1)
            else:
                Bstab.append(B[i])
    if Oflag == True:
        return(invert_orientation([Astab, Bstab]))
    else:
        return([Astab, Bstab])
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def thurston_bennequin(input_grid):
    r"""
    Computes the Thurston-Bennequin number of the Legendrian link represented by the grid.
    In the case of multi-component links, it computes the sum over all components.
    
    OUTPUT:

    The Thurston-Bennequin number of the grid.

    EXAMPLES::
    
    >> G = generate_torus_link(4,3)
    >> print(thurston_bennequin(G))
    -12
    >> F = stabilisation(G, 2, kind = 'XNW')
    >> print(thurston_bennequin(F))
    -13    
    
    """
    return int(writhe(input_grid) - (descending_cusps(input_grid) + ascending_cusps(input_grid))/2)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def uncoherent_bs(input_grid, where, which = 'rows'):
    r"""
    Performs a uncoherent band attachment on the grid. Can choose between a row/column
    band attachment (using the input 'which'), and specifying the number of the 
    row/column using the 'where' input. The syntax for the 'which' input is either
    'rows'/'r' or 'columns'/'c'.
    
    OUTPUT:
    
    The grid obtained after a uncoherent band attachment, or an error message if the 
    chosen band attachment site is not available.
    
    EXAMPLES::
    
    >> G = load_knot('4_1')
    >> print(G)
    [[5, 0, 3, 4, 2, 1], [2, 4, 5, 1, 0, 3]]
    >> F = uncoherent_bs(G, 1, 'c')
    >> F
    [[5, 4, 3, 1, 0, 2], [1, 0, 5, 4, 2, 3]]
    
    
    >> G = stabilisation(load_knot('4_1'),3,'XNE')
    >> F = uncoherent_bs(G,4,'c')

    Exception: The selected rows/columns are not a band attachment site
    
    """
    if which not in ['columns','rows','r','c']:
        raise Exception("Invalid input, the only options are 'rows', 'r', 'columns', 'c'.")
    if check_grid(input_grid) == 1:
        raise Exception("Invalid Input grid")
    if (which == 'r' or which == 'rows') and _check_rows_for_bands(input_grid,where)[0] != 0: 
        raise Exception("The selected rows/columns are not a band attachment site")   
    if (which == 'c' or which == 'columns') and _check_rows_for_bands(rotate(input_grid,1),where)[0] != 0:
        raise Exception("The selected rows/columns are not a band attachment site")
    if number_of_components(input_grid) !=1:
        raise Exception("The Input grid does not represent a knot")
    if which == 'columns' or which == 'c':
        input_grid = rotate(input_grid,1)
    case = _check_rows_for_bands(input_grid,where)    
    table = _from_grid_to_table(input_grid)
    aux = [table[where][0][0], table[where][0][1], table[where+1][0][0], table[where+1][0][1]] 
    aux.sort()
    if case[1]== 1 or case[1] == 3:
        table[where][0][0] = aux[2]
        table[where][0][1] = aux[3]
        table[where+1][0][0] = aux[0]
        table[where+1][0][1] = aux[1]
    if case[1] == 2 or case[1] == 6:
        table[where][0][0] = aux[0]
        table[where][0][1] = aux[3]
        table[where+1][0][0] = aux[1]
        table[where+1][0][1] = aux[2]        
    if case[1] == 4 or case[1] == 5:
        table[where][0][0] = aux[1]
        table[where][0][1] = aux[3]
        table[where+1][0][0] = aux[0]
        table[where+1][0][1] = aux[2]        
    if which == 'columns' or which == 'c':
        return(rotate(_from_table_to_grid(table),3))
    return(_from_table_to_grid(table))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def writhe(input_grid):
    r"""
    Computes the writhe of the grid; this is the number of positive crossing minus the
    number of negative crossings in the link diagram represented by the grid.
    In the case of multi-component links, it computes the sum over all components.
    
    OUTPUT:

    An integer.

    EXAMPLES::
    
    >> G = generate_torus_link(2,15)
    >> writhe(G)
    -15
    
    """
    if check_grid(input_grid) == 1:
        raise Exception("Invalid Input")    
    A = input_grid[0]
    B = input_grid[1]
    writhe_aux = 0
    n = len(A)
    for i in range(1,n-1):
        for j in range(min(A[i],B[i]),max(A[i],B[i])):
            if A[i] < B[i] and B.index(j) < i < A.index(j):
                writhe_aux -= 1
            if A[i] < B[i] and A.index(j) < i < B.index(j):
                writhe_aux += 1
            if A[i] > B[i] and B.index(j) < i < A.index(j):
                writhe_aux += 1
            if A[i] > B[i] and A.index(j) < i < B.index(j):
                writhe_aux -= 1
    return writhe_aux


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# What follows is a bunch of auxiliary functions
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _aux_braid(input_grid):    
    A = input_grid[0]
    B = input_grid[1]
    braid_gens = []
    for l in range(len(A)):
        if A[l] < B[l]:
            aux = _count_crossings_braid([A,B],l)
            if aux != 0:
                gens = []
                temp = []
                for s in range(aux):
                    temp.append(s+1+_strandsontheleft([A,B],l))
                gens = gens + temp[::-1]
                braid_gens = braid_gens + gens
        if A[l] > B[l]:
            aux = _count_crossings_braid([A,B],l)
            if aux > 0:
                gens = []
                for s in range(aux):
                    gens.append(-(s+1+_strandsontheleft([A,B],l)))
                braid_gens = braid_gens + gens
    return braid_gens

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _check_rows(input_grid, where):
    #gives True if the commutation is NOT interleaving
    A = input_grid[0]
    B = input_grid[1]
    if where == len(A)-1:
        return 6
    #adding extra check to avoid markers on the same column
    if A[where] == B[where+1] or B[where] == A[where+1]:
        return 5
    if (min(A[where],B[where]) < min(A[where+1],B[where+1]) < max(A[where],B[where]) < max(A[where+1],B[where+1]) ) or  (min(A[where+1],B[where+1]) < min(A[where],B[where]) <  max(A[where+1],B[where+1]) < max(A[where],B[where])):
        return(False)
    else:
        return(True)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _check_distance_one(input_grid):
    A = input_grid[0]
    B = input_grid[1]
    counter = 0
    indici = []
    for i in range(len(A)-1):
        if A[i] == B[i+1] :
            counter += 1
            indici.append(i)
    return(indici)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _check_columns(input_grid,where):
    return(_check_rows(rotate(input_grid,1),where))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _count_crossings_braid(input_grid,segment, writhe = 'True'):  
    ##if writhe == True also keeps track
    #of the signs of the crossings
    A = input_grid[0]
    B = input_grid[1]
    crossings_braid = 0
    if abs(A[segment] - B[segment]) == 1:
        return 0
    elif A[segment]<B[segment]:
        for t in range(len(A)):
            if t != segment:
                if A[segment] < A[t] < B[segment]:
                    if t<segment:
                        if B.index(A[t]) > segment or B.index(A[t]) < t:
                            crossings_braid += 1
                    if t>segment:
                        if segment < B.index(A[t]) < t:
                            crossings_braid += 1
    elif B[segment]<A[segment]:
        for t in range(len(A)):
            if t != segment:
                if B[segment] < A[t] < A[segment]:
                    if t<segment:
                        if B.index(A[t]) > segment or B.index(A[t]) < t:
                            crossings_braid += 1
                    if t>segment:
                        if segment < B.index(A[t]) < t:
                            crossings_braid += 1
    return crossings_braid
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _distance_markings(A,B,J):
    for I in range(len(B)):
        if A[J] == B[I]:
            return I-J

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _howmanystrands(input_grid):
    A = input_grid[0]
    B = input_grid[1]
    counter = 0
    for q in range(len(A)):
        if A.index(q) > B.index(q):
            counter = counter+1
    return counter
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _strandsontheleft(input_grid,segment):
    A = input_grid[0]
    B = input_grid[1]
    strands = 0
    if A[segment]<B[segment]:
        for t in range(len(A)):
            if t != segment:
                if A[t] < A[segment]:
                    if t<segment:
                        if B.index(A[t]) > segment or B.index(A[t]) < t:
                            strands += 1
                    if t>segment:
                        if segment < B.index(A[t]) < t:
                            strands += 1
    if A[segment]>B[segment]:
        for t in range(len(A)):
            if t != segment:
                if A[t] < B[segment]:
                    if t<segment:
                        if B.index(A[t]) > segment or B.index(A[t]) < t:
                            strands += 1
                    if t>segment:
                        if segment < B.index(A[t]) < t:
                            strands += 1
    return strands

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _from_grid_to_table(input_grid):
    #transforms a grid into a table (to perform uncoherent band attachments)
    A = input_grid[0]
    B = input_grid[1]
    out_table = []
    for i in range(len(A)):
        column = []
        column.append([min([A[i],B[i]]), max([A[i],B[i]])])
        out_table.append(column)    
    return(out_table)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _find_value_in_table(val,input_tab):
    for i in range(len(input_tab)):
        if input_tab[i][0][0] == val:
            return(i,0)
        if input_tab[i][0][1] == val:
            return(i,1)
    return(None,None)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _from_table_to_grid(input_table):
    #transforms a table into a grid
    A = len(input_table)*['not_done']
    B = len(input_table)*['not_done']
    finished = 0
    ii = A.index('not_done')
    A[ii] = input_table[ii][0][0]
    B[ii] = input_table[ii][0][1]
    input_table[ii][0][0] = 'done'
    input_table[ii][0][1] = 'done'
    value = A[ii]    
    while 'not_done' in set(A).union(set(B)):
        it,jj = _find_value_in_table(value,input_table)
        if it != None:
            B[it] = input_table[it][0][jj]
            A[it] = input_table[it][0][abs(jj-1)]
            input_table[it][0][jj] = 'done'
            input_table[it][0][abs(jj-1)] = 'done'
            value = A[it]
        else:
            return(A,B)
    return(A,B)    

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _check_rows_for_bands(input_grid, where):
    A = input_grid[0]
    B = input_grid[1]
    if where == len(A)-1:
        raise Exception('Invalid row input!')        
    if A[where] < A[where+1] < B[where] < B[where+1] or B[where] < B[where+1] < A[where] < A[where+1] or B[where+1] < B[where] < A[where+1] < A[where] or A[where+1] < A[where] < B[where+1] < B[where]:
        return([0,1])
    if A[where] < B[where+1] < B[where] < A[where+1] or B[where] < A[where+1] < A[where] < B[where+1] or B[where+1] < A[where] < A[where+1] < B[where] or A[where+1] < B[where] < B[where+1] < A[where]:
        return([0,2])
    if A[where] < A[where+1] < B[where+1] < B[where] or B[where] < B[where+1] < A[where+1] < A[where] or A[where+1] < A[where] < B[where] < B[where+1] or B[where+1] < B[where] < A[where] < A[where+1]: 
        return([0,3])
    if A[where] < B[where+1] < A[where+1] < B[where] or B[where] < A[where+1] < B[where+1] < A[where] or A[where+1] < B[where] < A[where] < B[where+1] or B[where+1] < A[where] < B[where] < A[where+1]:
        return([0,4])
    if A[where+1] < B[where+1] < A[where] < B[where] or B[where+1] < A[where+1] < B[where] < A[where] or A[where] < B[where] < A[where+1] < B[where+1] or B[where] < A[where] < B[where+1] < A[where+1]:
        return([0,5])
    if B[where+1] < A[where+1] < A[where] < B[where] or A[where+1] < B[where+1] < B[where] < A[where] or B[where] < A[where] < A[where+1] < B[where+1] or A[where] < B[where] < B[where+1] < A[where+1]:
        return([0,6])
    else:
        return([1,None])
            
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _invert(M):
    return M[1],M[0]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _destabilize_aux(input_grid,where):   
    if check_grid(input_grid) == 1:
        raise Exception("Invalid Input")
    A = input_grid[0]
    B = input_grid[1]
    n = len(A)
    Anew = (n-1)*[0]
    Bnew = (n-1)*[0]
    j = A[where]   
    if B[where] == A[where+1]: ### check that A,B is not a trivial rectangle
        return(A,B,1)
    for kk in range(where):
        if A[kk] <j:
            Anew[kk] = A[kk]
        else:
            Anew[kk] = A[kk]-1
        if B[kk] <j:    
            Bnew[kk] = B[kk]
        else:
            Bnew[kk] = B[kk]-1
    for kk in range(where+1,n-1):
        if A[kk+1] <j:
            Anew[kk] = A[kk+1]
        else:
            Anew[kk] = A[kk+1]-1
        if B[kk+1] <j:
            Bnew[kk] = B[kk+1]
        else:
            Bnew[kk] = B[kk+1]-1
    if B[where] < A[where] < A[where+1]:
        Bnew[where] = B[where]
        Anew[where] = A[where+1] -1
    if B[where] < A[where] and A[where+1] < A[where]:
        Bnew[where] = B[where]
        Anew[where] = A[where+1]
    if B[where] > A[where] > A[where+1]:
        Bnew[where] = B[where] -1
        Anew[where] = A[where+1]
    if B[where] > A[where] and  A[where+1] > A[where]:
        Bnew[where] = B[where] -1
        Anew[where] = A[where+1] -1
    return(Anew,Bnew,0)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _can_simplify(input_grid):
    A = input_grid[0]
    B = input_grid[1]
    cc = _check_distance_one([A, B])
    ccr = _check_distance_one(rotate([A,B],1))
    cci = _check_distance_one(_invert([A,B]))
    ccir = _check_distance_one(_invert(rotate([A,B],1)))
   
    if len(cc) > 0 or len(ccr) >0 or len(cci) > 0 or len(ccir) >0:
        return True
    else:
        return False

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




