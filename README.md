# GridPythonModule
A Sage compatible Python module to manipulate and simplify grid diagrams.

Examples and instructions can be found in the two Jupyter notebooks "Grid_introduction" and "Crossing_changes_and_bands".

For more details we refer to the associated paper "[GridPyM: a Python module to handle grid diagrams](https://arxiv.org/abs/2210.07399)" by Agnese Barbensi and Daniele Celoria.

**To install**, write `pip install GridPythonModule` in a terminal.

**To import**, write in a Python terminal `import GridPythonModule`, and `from GridPythonModule import *` to import its functions. 

**Requirements**: Python3. 

The following packages are needed: matplotlib, sympy, random2

**List of available functions**: 'ascending_cusps', 'available_knots', 'available_legendrian_knots', 'check_grid', 'coherent_bs', 'commute_columns', 'commute_rows', 'connected_sum', 'convert_to_Sage', 'convert_to_braid', 'crossing_number', 'cyclic_shift', 'descending_cusps', 'destabilize', 'destabilize_all', 'disjoint_union', 'draw_grid', 'Gauss_code', 'generate_random_grid', 'generate_torus_link', 'generate_twist_knot', 'generate_unknot', 'generate_unlink', 'grid_length', 'grid_number', 'invert_orientation', 'load_knot', 'load_legendrian_knot', 'mirror_grid', 'number_of_components', 'parallel_copies', 'perform_all_moves', 'rotate', 'rotate_once', 'rotation_number', 'scramble_grid', 'self_linking', 'simplify_grid', 'stabilisation', 'thurston_bennequin', 'uncoherent_bs', 'writhe'.
