from mp_api.client import MPRester
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# from torch_geometric.data import HeteroData

def load(API_KEY, MATERIAL_ID, conventional=True):
    '''helper function to load structural information of Material

    Parameters
    ----------
    conventional: bool
        if True, returns conventional standard structure else returns primitive
    '''
    with MPRester(API_KEY) as mpr:

        # first retrieve the relevant structure
        structure = mpr.get_structure_by_material_id(MATERIAL_ID)
    
    # important to use the conventional structure to ensure
    # that peaks are labelled with the conventional Miller indices
    sga = SpacegroupAnalyzer(structure)
    conventional_structure = sga.get_conventional_standard_structure()
    
    if conventional:
        structure = conventional_structure

    return structure


class Solid:
    def __init__( self, API_KEY, MATERIAL_ID = 'mp-1103503'):
        '''Material class. To store and process material data

        Parameters
        ----------
        API_KEY : str
            secret API KEY used to run MPRester data requests
        MATERIAL_ID : int
            ID of selected material (without "mp-" prefix)
        graph : object
            Heterogeneous graph to store all information from a `Solid` material
        '''
        self.API_KEY = API_KEY
        self.MATERIAL_ID = MATERIAL_ID
        self.graph = {}

    def load(self, conventional=True):
        '''helper function to load structural information of Material

        Parameters
        ----------
        conventional: bool
            if True, returns conventional standard structure else returns primitive
        '''
        with MPRester(self.API_KEY) as mpr:

            # first retrieve the relevant structure
            structure = mpr.get_structure_by_material_id(self.MATERIAL_ID)
        
        # important to use the conventional structure to ensure
        # that peaks are labelled with the conventional Miller indices
        sga = SpacegroupAnalyzer(structure)
        conventional_structure = sga.get_conventional_standard_structure()
        
        if conventional:
            structure = conventional_structure

        return structure


    def topology(self, fractional=False):
        '''processes the material topology as a hashmap with an a-xyz coordinate format 
        i.e. [atom_numbers, x_coord, y_coord, z_coord]

        Parameters
        ----------
        fractional: bool
            if False, returns xyz coordinates else returns fractional (a,b,c) coordinates
        '''
        structure = Solid(self.API_KEY, self.MATERIAL_ID).load(conventional=True)

        elements = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118}

        Nsites = len(structure.sites)

        # nodes = HeteroData()

        nodes = {}
        ckeys = ['x','y','z']
        for i in ['atomic_number',*ckeys]:
            nodes[i] = [] 

        if not fractional:
            for i in range(Nsites):
                # atomic_num = element(str(structure.sites[i].specie)).atomic_number # specie is not a typo
                atomic_num = elements[str(structure.sites[i].specie)] # specie is not a typo

                # nodes['atomic_number'].x.append(atomic_num)
                # [ nodes[ckeys[k]].x.append(structure.sites[i].coords[k]) for k in range(3) ]

                nodes['atomic_number'].append(atomic_num)
                [ nodes[ckeys[k]].append(structure.sites[i].coords[k]) for k in range(3) ]

        else:
            for i in range(Nsites):
                # atomic_num = element(str(structure.sites[i].specie)).atomic_number 
                atomic_num = elements[str(structure.sites[i].specie)] # specie is not a typo

                # nodes['atomic_number'].append(atomic_num)
                # [ nodes[ckeys[k]].x.append(structure.sites[i].frac_coords[k]) for k in range(3) ]

                nodes['atomic_number'].append(atomic_num)
                [ nodes[ckeys[k]].append(structure.sites[i].frac_coords[k]) for k in range(3) ]

        return nodes




