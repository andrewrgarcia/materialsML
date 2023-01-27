from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.local_env import *
import numpy as np
import json
import matplotlib.pyplot as plt


def view__legacy(topol_info, figsize=(8, 6), dpi=80, node_size = 1500, tube_radius = 5, edge = 5 ):
    '''Visualization of Material with designated topology
    Legacy code (older)
    Parameters
    ----------
    topol_info: dict() --> atomic_number, x,y,z
        dictionary containing the above information
    node: int
        atom / cluster size
    edge: int
        bond width
    '''

    plt.figure(figsize=figsize, dpi=dpi)

    ax = plt.axes(projection='3d')
    # colors = np.linspace(2**20,2**24,118,dtype='int') #divide color range into 118 colors (for the 118 chemical elements)
    jmol_colors = {
                'H': '#FFFFFF', 'He': '#D9FFFF', 'Li': '#CC80FF', 'Be': '#C2FF00', 'B': '#FFB5B5', 
                'C': '#909090', 'N': '#3050F8', 'O': '#FF0D0D', 'F': '#90E050', 'Ne': '#B3E3F5',
                'Na': '#AB5CF2', 'Mg': '#8AFF00', 'Al': '#BFA6A6', 'Si': '#F0C8A0', 'P': '#FF8000',
                'S': '#FFFF30', 'Cl': '#1FF01F', 'Ar': '#80D1E3', 'K': '#8F40D4', 'Ca': '#3DFF00', 
                'Sc': '#E6E6E6', 'Ti': '#BFC2C7', 'V': '#A6A6AB', 'Cr': '#8A99C7', 'Mn': '#9C7AC7',
                'Fe': '#E06633', 'Co': '#F090A0', 'Ni': '#50D050', 'Cu': '#C88033', 'Zn': '#7D80B0',
                'Ga': '#C28F8F', 'Ge': '#668F8F', 'As': '#BD80E3', 'Se': '#FFA100', 'Br': '#A62929',
                'Kr': '#5CB8D1', 'Rb': '#702EB0', 'Sr': '#00FF00', 'Y': '#94FFFF', 'Zr': '#94E0E0', 
                'Nb': '#73C2C9', 'Mo': '#54B5B5', 'Tc': '#3B9E9E', 'Ru': '#248F8F', 'Rh': '#0A7D8C', 
                'Pd': '#006985', 'Ag': '#C0C0C0', 'Cd': '#FFD98F', 'In': '#A67573', 'Sn': '#668080', 
                'Sb': '#9E63B5', 'Te': '#D47A00', 'I': '#940094', 'Xe': '#429EB0', 'Cs': '#57178F', 
                'Ba': '#00C900', 'La': '#70D4FF', 'Ce': '#FFFFC7', 'Pr': '#D9FFC7', 'Nd': '#C7FFC7', 
                'Pm': '#A3FFC7', 'Sm': '#8FFFC7', 'Eu': '#61FFC7', 'Gd': '#45FFC7', 'Tb': '#30FFC7', 
                'Dy': '#1FFFC7', 'Ho': '#00FF9C', 'Er': '#00E675', 'Tm': '#00D452', 'Yb': '#00BF38', 
                'Lu': '#00AB24', 'Hf': '#4DC2FF', 'Ta': '#4DA6FF', 'W': '#2194D6', 'Re': '#267DAB', 
                'Os': '#266696', 'Ir': '#175487', 'Pt': '#D0D0E0', 'Au': '#FFD123', 'Hg': '#B8B8D0', 
                'Tl': '#A6544D', 'Pb': '#575961', 'Bi': '#9E4FB5', 'Po': '#AB5C00', 'At': '#754F45', 
                'Rn': '#428296', 'Fr': '#420066', 'Ra': '#007D00', 'Ac': '#70ABFA', 'Th': '#00BAFF', 
                'Pa': '#00A1FF', 'U': '#008FFF', 'Np': '#0080FF', 'Pu': '#006BFF', 'Am': '#545CF2',
                'Cm': '#785CE3', 'Bk': '#8A4FE3', 'Cf': '#A136D4', 'Es': '#B31FD4', 'Fm': '#B31FBA',
                'Md': '#B30DA6', 'No': '#BD0D87', 'Lr': '#C70066', 'Rf': '#CC0059', 'Db': '#D1004F',
                'Sg': '#D90045', 'Bh': '#E00038', 'Hs': '#E6002E', 'Mt': '#EB0026'
                }


    atom = topol_info['atom']

    xyz_arr = np.array([topol_info[j] for j in list('xyz')])

    # print( xyz_arr.T )

    if 'bond_edges' in topol_info.keys():
        for i in range(len(atom)):
            NNidcs = topol_info['bond_edges'][i]
            
            x0,y0,z0 = xyz_arr.T[i]

            x =  [topol_info['x'][k] for k in NNidcs ] 
            y =  [topol_info['y'][k] for k in NNidcs ] 
            z =  [topol_info['z'][k] for k in NNidcs ] 

            [ax.plot((x0,x[k]),(y0,y[k]),(z0,z[k]), linewidth = tube_radius,\
                color=jmol_colors[atom[i]], alpha=0.5) for k in range(len(NNidcs))]
             

    for i in range(len(atom)):
        
        # ax.scatter3D(*xyz_arr.T[i], s=1500, linewidths=3, edgecolors='#dddddd', c="#"+hex(colors[atom[i]])[2:])
        ax.scatter3D(*xyz_arr.T[i], s=node_size, linewidths=3, edgecolors='#dddddd', c=jmol_colors[atom[i]] )

    ax.set_facecolor('#0e0e12')
    # ax.set_facecolor('#626279')

    
    set_axes_equal(ax)           

    plt.axis('off')
    plt.show()


def set_axes_radius(ax, origin, radius):
    '''set_axes_radius and set_axes_equal * * * Credit:
    Mateen Ulhaq (answered Jun 3 '18 at 7:55)
    https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to'''
    ax.set_xlim3d([origin[0] - radius, origin[0] + radius])
    ax.set_ylim3d([origin[1] - radius, origin[1] + radius])
    ax.set_zlim3d([origin[2] - radius, origin[2] + radius])
    
def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.
    Input
      ax: a matplotlib axis, e.g., as output from plt.add_subplot().
    '''
    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])

    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    set_axes_radius(ax, origin, radius)


def save_to_json(filename="graphs.json",class_object={}):
    '''Save object to JSON file
    adapted from: https://pythonspot.com/save-a-dictionary-to-a-file/
    '''
    # load json module
    # python dictionary with key value pairs
    dict = class_object
    # create json object from dictionary
    jsonobj = json.dumps(dict)
    # open file for writing, "w" 
    f = open(filename,"w")
    # write json object to file
    f.write(jsonobj)
    # close file
    f.close()

def load_from_json(filename="graphs.json"):
    '''Load JSON file to object'''
    with open(filename) as f:
        data = f.read()
    return json.loads(data)


mendeleev_elements = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118}

atomic_radii = {'H':  0.53, 'He':  0.31, 'Li':  1.67, 'Be':  1.12, 'B':  0.87, 'C':  0.67, 'N':  0.56, 'O':  0.48, 'F':  0.42, 'Ne':  0.38, 'Na':  1.9, 'Mg':  1.45, 'Al':  1.18, 'Si':  1.11, 'P':  0.98, 'S':  0.88, 'Cl':  0.79, 'Ar':  0.71, 'K':  2.43, 'Ca':  1.94, 'Sc':  1.84, 'Ti':  1.76, 'V':  1.71, 'Cr':  1.66, 'Mn':  1.61, 'Fe':  1.56, 'Co':  1.52, 'Ni':  1.49, 'Cu':  1.45, 'Zn':  1.42, 'Ga':  1.36, 'Ge':  1.25, 'As':  1.14, 'Se':  1.03, 'Br':  0.94, 'Kr':  0.88, 'Rb':  2.65, 'Sr':  2.19, 'Y':  2.12, 'Zr':  2.06, 'Nb':  1.98, 'Mo':  1.9, 'Tc':  1.83, 'Ru':  1.78, 'Rh':  1.73, 'Pd':  1.69, 'Ag':  1.65, 'Cd':  1.61, 'In':  1.56, 'Sn':  1.45, 'Sb':  1.33, 'Te':  1.23, 'I':  1.15, 'Xe':  1.08, 'Cs':  2.98, 'Ba':  2.53, 'La':  1.95, 'Ce':  1.85, 'Pr':  2.47, 'Nd':  2.06, 'Pm':  2.05, 'Sm':  2.38, 'Eu':  2.31, 'Gd':  2.33, 'Tb':  2.25, 'Dy':  2.28, 'Ho':  2.26, 'Er':  2.26, 'Tm':  2.22, 'Yb':  2.22, 'Lu':  2.17, 'Hf':  2.08, 'Ta':  2.0, 'W':  1.93, 'Re':  1.88, 'Os':  1.85, 'Ir':  1.8, 'Pt':  1.77, 'Au':  1.74, 'Hg':  1.71, 'Tl':  1.56, 'Pb':  1.54, 'Bi':  1.43, 'Po':  1.35, 'At':  1.27, 'Rn':  1.2, 'Fr':  0.0001, 'Ra':  0.0001, 'Ac':  1.95, 'Th':  1.8, 'Pa':  1.8, 'U':  1.75, 'Np':  1.75, 'Pu':  1.75, 'Am':  1.75, 'Cm':  0.0001}


def topol(structure, edges="bond_edges", fractional=False, conventional=True):
    '''processes the material topology as a hashmap with an a-xyz coordinate format 
    i.e. [atom_numbers, x_coord, y_coord, z_coord]
    from a MPR[material-id].summary.structure object

    Parameters
    ----------
    fractional: bool
        if False, returns xyz coordinates else returns fractional (a,b,c) coordinates
    conventional: bool
        if True, returns conventional standard structure else returns primitive
    '''
    # structure0 = structure.copy()
    # important to use the conventional structure to ensure
    # that peaks are labelled with the conventional Miller indices
    sga = SpacegroupAnalyzer(structure)
    
    
    if conventional:
        # structure = sga.get_conventional_standard_structure()
        structure = sga.get_refined_structure()

    nodes = {}
    ckeys = list('yxz')
    topol_list = ['atomic_number','atom',*ckeys]

    if edges: 
        topol_list.append('bond_edges')

    for i in topol_list:
        nodes[i] = []     

    
    for i in range(len(structure.sites)):
        # atomic_num = element(str(structure.sites[i].specie)).atomic_number # specie is not a typo

        nodes['atomic_number'].append( mendeleev_elements[str(structure.sites[i].specie)] ) # specie is not a typo
        nodes['atom'].append(str(structure.sites[i].specie))

        if not fractional:
            for k in range(3): nodes[ckeys[k]].append(structure.sites[i].coords[k]) 
        else:
            for k in range(3): nodes[ckeys[k]].append(structure.sites[i].frac_coords[k])


    if edges:
        all_coords = np.array([k.coords for k in structure])

        NEIGHBORS=[]
        for i in range(len(structure)):
            atom = structure.sites[i].specie
            offset = 2.1
            atom_mult = atomic_radii[str(atom)] * offset
            ns = structure.get_neighbors(structure[i],atom_mult)
            neighbors_i = []
            for j in ns:
                locs = np.where((j.coords == all_coords).all(1))[0]
                if len(locs) > 0:
                    neighbors_i.append(int(locs[0]))
            NEIGHBORS.append(neighbors_i)

        nodes['bond_edges'] = NEIGHBORS

    return nodes
