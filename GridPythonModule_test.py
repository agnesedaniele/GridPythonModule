import GridPythonModule 
from GridPythonModule import *
import pytest

test_grid = [[0,1,2,3,4],[2,3,4,0,1]]
no_grid   = [0,1,2,3,4],[2,3,4,0,0]
test_link = [[2, 4, 0, 3, 5, 1], [5, 2, 3, 0, 1, 4]]
moves = ['XNE','XNW','XSE','XSW','ONE','ONW','OSE','OSW']

######################################################################################################################


def test_check_grid():

    assert check_grid(test_grid) == 0

    assert check_grid(no_grid) == 1

    assert check_grid([test_grid[0]]) == 1

######################################################################################################################


def test_available_knots():

    with pytest.raises(Exception) as exc_info:
         available_knots(10)
    assert str(exc_info.value) == 'We only have links up to 8 crossings!'


    with pytest.raises(Exception) as exc_info:
         available_knots(-1)
    assert str(exc_info.value) == 'We only have links up to 8 crossings!'

    with pytest.raises(Exception) as exc_info:
         available_knots(False)
    assert str(exc_info.value) == "Invalid crossing number input"


######################################################################################################################


def test_available_legendrian_knots():

    with pytest.raises(Exception) as exc_info:
         available_legendrian_knots(10)
    assert str(exc_info.value) == 'We only have nontrivial knots up to 7 crossings!'


######################################################################################################################


def test_load_knot():

    assert load_knot('3_1') == [[4, 0, 1, 2, 3], [1, 2, 3, 4, 0]]

    with pytest.raises(Exception) as exc_info:
         load_knot('trefoil')
    assert str(exc_info.value) == "Invalid input name! Try the command 'available_knots' to see which knots are pre-loaded."


######################################################################################################################


def test_load_legendrian_knot():

    assert load_knot('3_1') == [[4, 0, 1, 2, 3], [1, 2, 3, 4, 0]]

    with pytest.raises(Exception) as exc_info:
         load_knot('trefoil')
    assert str(exc_info.value) == "Invalid input name! Try the command 'available_knots' to see which knots are pre-loaded."

######################################################################################################################


def test_generate_torus_link():

    assert generate_torus_link(3,2) == [[4, 3, 2, 1, 0], [2, 1, 0, 4, 3]]

    with pytest.raises(Exception) as exc_info:
         generate_torus_link(0,2)
    assert str(exc_info.value) == "Only non-zero coefficients"

    with pytest.raises(Exception) as exc_info:
         generate_torus_link(2,0)
    assert str(exc_info.value) == "Only non-zero coefficients"

    assert generate_torus_link(3,-2) == [[0, 1, 2, 3, 4], [3, 4, 0, 1, 2]]

######################################################################################################################


def test_generate_twist_knot():

    assert generate_twist_knot(3, '+') == [[3, 2, 5, 4, 1, 0, 6], [0, 4, 3, 6, 5, 2, 1]]

    assert generate_twist_knot(3, '-') == [[0, 5, 3, 2, 4, 7, 6, 1], [4, 2, 1, 0, 6, 5, 3, 7]]

    assert generate_twist_knot(-3, '+') == [[5, 6, 1, 4, 2, 3, 7, 0], [7, 4, 5, 0, 6, 1, 2, 3]]

    assert generate_twist_knot(-3, '-') == [[1, 2, 5, 6, 3, 4, 0], [6, 0, 1, 4, 5, 2, 3]]


    with pytest.raises(Exception) as exc_info:
         generate_twist_knot(0)
    assert str(exc_info.value) == "There are easier ways of creating unknots!"

    with pytest.raises(Exception) as exc_info:
         generate_twist_knot(3, '?')
    assert str(exc_info.value) == "Invalid clasp type"

######################################################################################################################


def test_generate_unknot():

    assert generate_unknot(3) == [[0, 1, 2], [1, 2, 0]]

    with pytest.raises(Exception) as exc_info:
         generate_unknot(2)
    assert str(exc_info.value) == "The grid size needs to be an integer greater than 2!"

######################################################################################################################


def test_generate_unlink():

    assert generate_unlink(2) == [[0, 1, 2, 3], [1, 0, 3, 2]]

    assert generate_unlink(3) == [[0, 1, 2, 3, 4, 5], [1, 0, 3, 2, 5, 4]]

    with pytest.raises(Exception) as exc_info:
         generate_unlink(0)
    assert str(exc_info.value) == "The number of components needs to be positive!"

######################################################################################################################


def test_number_of_components():

    assert number_of_components(test_grid) ==1

    assert number_of_components(test_link) == 2

######################################################################################################################


def test_invert_orientation():

    assert invert_orientation(test_grid) == [test_grid[1],test_grid[0]]

######################################################################################################################

def test_ascending_cusps():

    assert ascending_cusps(test_grid) == 2


######################################################################################################################


def test_descending_cusps():

    assert descending_cusps(test_grid) == 4

    with pytest.raises(Exception) as exc_info:
        descending_cusps(no_grid)
    assert str(exc_info.value) == 'Invalid Input'

######################################################################################################################


def test_cyclic_shift():

    assert cyclic_shift(test_grid, 1,0) == [[4, 0, 1, 2, 3], [1, 2, 3, 4, 0]]

    with pytest.raises(Exception) as exc_info:
        cyclic_shift(no_grid,1,0)
    assert str(exc_info.value) == 'Invalid Input'

    with pytest.raises(Exception) as exc_info:
        cyclic_shift(test_grid,-1,0)
    assert str(exc_info.value) == 'Invalid Input'

    with pytest.raises(Exception) as exc_info:
        cyclic_shift(test_grid,1,-2)
    assert str(exc_info.value) == 'Invalid Input'

######################################################################################################################


def test_simplify_grid():

    assert sorted(simplify_grid(generate_unknot(10), effort = 'high')) == [[0, 1], [1, 0]]

    with pytest.raises(Exception) as exc_info:
        simplify_grid(no_grid, 'high')
    assert str(exc_info.value) == 'Invalid Input'

    with pytest.raises(Exception) as exc_info:
        simplify_grid(generate_unknot(10), effort = 'else')
    assert str(exc_info.value) == 'Invalid effort!'

    with pytest.raises(Exception) as exc_info:
        simplify_grid(generate_unknot(10), effort = 0)
    assert str(exc_info.value) == 'Invalid effort!'

######################################################################################################################


def test_scramble_grid():

    assert check_grid(scramble_grid(test_grid)) == 0

    with pytest.raises(Exception) as exc_info:
        scramble_grid(no_grid, 'high')
    assert str(exc_info.value) == 'Invalid Input'

    with pytest.raises(Exception) as exc_info:
        scramble_grid(generate_unknot(10), effort = 'else')
    assert str(exc_info.value) == 'Invalid effort!'

    with pytest.raises(Exception) as exc_info:
        scramble_grid(generate_unknot(10), effort = 0)
    assert str(exc_info.value) == 'Invalid effort!'

######################################################################################################################


def test_stabilisation():

    for move in moves:
        for pos in range(0,len(test_grid[0])):
            assert len(stabilisation(test_grid, pos, move)[0]) == 6
    assert stabilisation(test_grid, 2, 'XSE') == [[0, 1, 2, 3, 4, 5], [3, 4, 5, 2, 0, 1]]

    with pytest.raises(Exception) as exc_info:
        stabilisation(test_grid, 2, 'w')
    assert str(exc_info.value) == 'Invalid kind of stabilisation!'

    with pytest.raises(Exception) as exc_info:
        stabilisation(no_grid,0, 'XSE')
    assert str(exc_info.value) == 'Invalid Input'


######################################################################################################################


def test_destabilize():

    with pytest.raises(Exception) as exc_info:
        destabilize(no_grid,0)
    assert str(exc_info.value) == 'Invalid Input'

    with pytest.raises(Exception) as exc_info:
        destabilize([[1,0], [0,1]],0, verbose = True)
    assert str(exc_info.value) == 'Grids of size 2 cannot be destabilized'

    assert destabilize([[1,0], [0,1]],0, verbose = False) == 0 

    with pytest.raises(Exception) as exc_info:
        destabilize(test_grid,-1)
    assert str(exc_info.value) == 'The index must be between 0 and grid number'

    for move in moves:
        for pos in range(0,len(test_grid[0])):

            assert destabilize(stabilisation(test_grid, pos, move), pos) == 0 or destabilize(stabilisation(test_grid, pos, move), pos) == test_grid

            assert destabilize(stabilisation(test_grid, pos, move), pos, 'col') == 0 or destabilize(stabilisation(test_grid, pos, move), pos, 'col') == test_grid

######################################################################################################################


def test_destabilize_all():

    for move in moves:
        for pos in range(0,len(test_grid[0])):
            assert destabilize_all(stabilisation(test_grid, pos, move)) == test_grid
            
######################################################################################################################


def test_grid_number():

    assert grid_number(test_grid) == 5

    with pytest.raises(Exception) as exc_info:
        assert grid_number(no_grid)
    assert str(exc_info.value) == 'Invalid Input'

######################################################################################################################


def test_grid_length():

    assert grid_length(test_grid) == 24 

######################################################################################################################


def test_convert_to_braid():

    assert convert_to_braid(test_grid) == [-1,-1,-1]

    assert convert_to_braid(test_grid, optimized = 'N') == [-1,-2,-1, -2]

    assert convert_to_braid(test_grid, optimized = 'S') == [-1,-1,-1]

    with pytest.raises(Exception) as exc_info:
        assert convert_to_braid(no_grid)
    assert str(exc_info.value) == 'Invalid Input'

    with pytest.raises(Exception) as exc_info:
        assert convert_to_braid(test_grid, 'other')
    assert str(exc_info.value) == 'Invalid optimization input.'

######################################################################################################################


def test_crossing_number():

    assert crossing_number(test_grid) == 3

######################################################################################################################


def test_self_linking():

    assert self_linking(invert_orientation(test_grid)) == -5

    assert self_linking(test_grid) == -7

    assert self_linking(mirror_grid(test_grid)) == 1

######################################################################################################################


def test_mirror_grid():

    assert mirror_grid(test_grid) == [[4, 3, 2, 1, 0], [1, 0, 4, 3, 2]]

    assert mirror_grid(mirror_grid(mirror_grid(mirror_grid(test_grid)))) == [test_grid[0], test_grid[1]]

######################################################################################################################

def test_writhe():

    assert writhe(test_grid) == -3

    assert writhe(mirror_grid(test_grid)) == 3

    with pytest.raises(Exception) as exc_info:
        assert writhe(no_grid)
    assert str(exc_info.value) == 'Invalid Input'

######################################################################################################################

def test_thurston_bennequin():

    assert thurston_bennequin(test_grid) == -6

######################################################################################################################

def test_rotation_number():

    assert rotation_number(test_grid) == 1

    assert rotation_number(mirror_grid(test_grid)) ==0

######################################################################################################################

def test_connected_sum():

    assert connected_sum(test_grid, mirror_grid(test_grid)) == [[0, 1, 2, 3, 4, 5, 9, 8, 7, 6], [2, 3, 5, 0, 1, 7, 6, 4, 9, 8]]

    assert connected_sum(test_grid, test_grid) == [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], [2, 3, 5, 0, 1, 7, 8, 9, 4, 6]]

    with pytest.raises(Exception) as exc_info:
        connected_sum(test_grid, no_grid)
    assert str(exc_info.value) == 'Invalid Input'

######################################################################################################################


def test_disjoint_union():

    assert disjoint_union(test_grid, mirror_grid(test_grid)) == [[0, 1, 2, 3, 4, 9, 8, 7, 6, 5], [2, 3, 4, 0, 1, 6, 5, 9, 8, 7]]

    with pytest.raises(Exception) as exc_info:
        disjoint_union(test_grid, no_grid)
    assert str(exc_info.value) == 'Invalid Input'

######################################################################################################################


def test_parallel_copies():
    
    with pytest.raises(Exception) as exc_info:
        parallel_copies(no_grid,3)
    assert str(exc_info.value) == 'Invalid Input'

    with pytest.raises(Exception) as exc_info:
        parallel_copies(test_grid, -1)
    assert str(exc_info.value) == 'Invalid number of strands'

    assert parallel_copies(test_grid,3) == [[0, 1, 2, 3, 4, 5, 8, 7, 6, 9, 10, 11, 12, 13, 14],[8, 7, 6, 11, 10, 9, 14, 13, 12, 2, 1, 0, 5, 4, 3]]

######################################################################################################################


def test_Gauss_code():

    assert Gauss_code(test_grid) == [[[-1, 2, -3, 1, -2, 3]], [-1, -1, -1]]

    with pytest.raises(Exception) as exc_info:
        Gauss_code(no_grid)
    assert str(exc_info.value) == 'Invalid grid Input!'

    with pytest.raises(Exception) as exc_info:
        Gauss_code(test_link)
    assert str(exc_info.value) == 'Right now the Gauss code is only implemented for knots!'

######################################################################################################################

def test_generate_random_grid():

    for i in range(5,20):
        for j in range(1,3):

            assert grid_number(generate_random_grid(i,j)) == i

            assert number_of_components(generate_random_grid(i,j)) == j

    with pytest.raises(Exception) as exc_info:
        generate_random_grid(4,5)
    assert str(exc_info.value) == 'Wrong number of components or incorrect grid size!'

######################################################################################################################


def test_perform_all_moves():

    assert perform_all_moves(test_grid) == [[[1, 0, 2, 3, 4, 5], [0, 3, 4, 5, 1, 2]],
 [[0, 1, 2, 3, 4, 5], [1, 3, 4, 5, 0, 2]],
 [[0, 1, 2, 3, 4, 5], [3, 0, 4, 5, 1, 2]],
 [[1, 0, 2, 3, 4, 5], [3, 1, 4, 5, 0, 2]],
 [[2, 0, 1, 3, 4, 5], [3, 2, 4, 5, 0, 1]],
 [[3, 0, 1, 2, 4, 5], [2, 3, 4, 5, 0, 1]],
 [[0, 2, 1, 3, 4, 5], [2, 3, 4, 5, 0, 1]],
 [[0, 3, 1, 2, 4, 5], [3, 2, 4, 5, 0, 1]],
 [[0, 2, 1, 3, 4, 5], [3, 1, 4, 5, 0, 2]],
 [[0, 1, 2, 3, 4, 5], [3, 2, 4, 5, 0, 1]],
 [[0, 1, 2, 3, 4, 5], [3, 4, 1, 5, 0, 2]],
 [[0, 2, 1, 3, 4, 5], [3, 4, 2, 5, 0, 1]],
 [[0, 3, 1, 2, 4, 5], [2, 4, 3, 5, 0, 1]],
 [[0, 4, 1, 2, 3, 5], [2, 3, 4, 5, 0, 1]],
 [[0, 1, 3, 2, 4, 5], [2, 3, 4, 5, 0, 1]],
 [[0, 1, 4, 2, 3, 5], [2, 4, 3, 5, 0, 1]],
 [[0, 1, 3, 2, 4, 5], [3, 4, 2, 5, 0, 1]],
 [[0, 1, 2, 3, 4, 5], [2, 4, 3, 5, 0, 1]],
 [[0, 1, 2, 3, 4, 5], [3, 4, 5, 2, 0, 1]],
 [[0, 1, 3, 2, 4, 5], [2, 4, 5, 3, 0, 1]],
 [[0, 1, 4, 2, 3, 5], [2, 3, 5, 4, 0, 1]],
 [[0, 1, 5, 2, 3, 4], [2, 3, 4, 5, 0, 1]],
 [[0, 1, 2, 4, 3, 5], [2, 3, 4, 5, 0, 1]],
 [[0, 1, 2, 5, 3, 4], [2, 3, 5, 4, 0, 1]],
 [[0, 1, 2, 4, 3, 5], [2, 4, 5, 3, 0, 1]],
 [[0, 1, 2, 3, 4, 5], [2, 3, 5, 4, 0, 1]],
 [[0, 1, 2, 3, 4, 5], [2, 4, 5, 0, 3, 1]],
 [[0, 1, 2, 4, 3, 5], [2, 3, 5, 0, 4, 1]],
 [[1, 2, 3, 0, 4, 5], [3, 4, 5, 1, 0, 2]],
 [[0, 2, 3, 1, 4, 5], [3, 4, 5, 0, 1, 2]],
 [[1, 2, 3, 4, 0, 5], [3, 4, 5, 0, 1, 2]],
 [[0, 2, 3, 4, 1, 5], [3, 4, 5, 1, 0, 2]],
 [[0, 1, 2, 3, 5, 4], [2, 3, 5, 0, 4, 1]],
 [[0, 1, 2, 3, 4, 5], [2, 3, 4, 0, 5, 1]],
 [[0, 1, 2, 3, 4, 5], [2, 3, 5, 0, 1, 4]],
 [[0, 1, 2, 3, 5, 4], [2, 3, 4, 0, 1, 5]],
 [[0, 2, 3, 4, 1, 5], [3, 4, 5, 0, 2, 1]],
 [[0, 1, 3, 4, 2, 5], [3, 4, 5, 0, 1, 2]],
 [[0, 2, 3, 4, 5, 1], [3, 4, 5, 0, 1, 2]],
 [[0, 1, 3, 4, 5, 2], [3, 4, 5, 0, 2, 1]]]

    with pytest.raises(Exception) as exc_info:
        perform_all_moves(no_grid)
    assert str(exc_info.value) == 'Invalid Input'

######################################################################################################################


def test_rotate():

    assert rotate(test_grid, 1) == [[4, 3, 2, 1, 0], [1, 0, 4, 3, 2]]

    assert rotate(test_grid, 4) == test_grid

    random_grid = generate_random_grid(100)

    assert rotate(random_grid, 1) == mirror_grid(random_grid)

    assert rotate(random_grid, 4) == random_grid

######################################################################################################################


def test_rotate_once():
 
    assert rotate_once(test_grid[0]) == [4, 3, 2, 1, 0]

######################################################################################################################


def test_commute_columns():

   assert commute_columns(test_link, 0, 'A') == [[2, 4, 1, 3, 5, 0], [5, 2, 3, 1, 0, 4]]

   assert commute_columns(test_link, 0, 'A') == commute_columns(test_link, 0, 'N')

   commute_columns(test_link, 0, 'Y') == 0

   assert commute_columns(test_link, 4, 'A') == [[2, 5, 0, 3, 4, 1], [4, 2, 3, 0, 1, 5]]

   assert commute_columns(test_link, 4, 'A') == commute_columns(test_link, 4, 'Y')

   commute_columns(test_link, 4, 'N') == 0


######################################################################################################################


def test_commute_rows():

    assert commute_rows(test_link, 0, 'A') == commute_rows(test_link, 0, 'N') == commute_rows(test_link, 0, 'Y') == 0
    
    assert commute_rows(test_link, 1, 'A') == [[2, 0, 4, 3, 5, 1], [5, 3, 2, 0, 1, 4]]

    assert commute_rows(test_link, 1, 'A') ==  commute_rows(test_link, 1, 'Y') 

    assert commute_rows(test_link, 1, 'N') == 0

    with pytest.raises(Exception) as exc_info:
        commute_rows(no_grid, 1, 'A')
    assert str(exc_info.value) == 'Invalid Input grid'


######################################################################################################################


def test_uncoherent_bs():

    assert uncoherent_bs(test_grid,0, 'c') == [[2, 0, 4, 3, 1], [0, 3, 2, 1, 4]]

    with pytest.raises(Exception) as exc_info:
        uncoherent_bs(no_grid,0)
    assert str(exc_info.value) == 'Invalid Input grid'

    with pytest.raises(Exception) as exc_info:
        uncoherent_bs(test_link,1)
    assert str(exc_info.value) == 'The Input grid does not represent a knot'

    with pytest.raises(Exception) as exc_info:
        uncoherent_bs(generate_unknot(3),0)
    assert str(exc_info.value) == 'The selected rows/columns are not a band attachment site'
######################################################################################################################
  

def test_coherent_bs():
 
    with pytest.raises(Exception) as exc_info:
        coherent_bs(no_grid,0)
    assert str(exc_info.value) == 'Invalid Input grid'

    with pytest.raises(Exception) as exc_info:
        coherent_bs(test_grid,4, 'c') 
    assert str(exc_info.value) == 'Invalid row input!'

    with pytest.raises(Exception) as exc_info:
        coherent_bs(stabilisation(test_grid, 2, 'XNE'), 2) 
    assert str(exc_info.value)== 'The selected rows/columns are not a band attachment site'

    with pytest.raises(Exception) as exc_info:
        coherent_bs(test_grid,1, 'w') 
    assert str(exc_info.value)== "Invalid input, the only options are 'rows', 'r', 'columns', 'c'."

    assert coherent_bs(test_grid,0, 'r') == ([1, 0, 2, 3, 4], [2, 3, 4, 0, 1])  


######################################################################################################################


def test_draw_grid():

    with pytest.raises(Exception) as exc_info:
        draw_grid(no_grid)
    assert str(exc_info.value) == 'Invalid grid Input!'

    with pytest.raises(Exception) as exc_info:
        draw_grid(test_grid, markings = 'other')
    assert str(exc_info.value) == 'Invalid markings option.'




    

