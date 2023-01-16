from mp_api.client import MPRester
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

mendeleev_elements = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118}

class Solid:
    def __init__( self, API_KEY, MATERIAL_ID = 'mp-1103503'):
        '''Material class. To store and process material data

        Parameters
        ----------
        API_KEY : str
            secret API KEY used to run MPRester data requests
        MATERIAL_ID : str
            ID of selected material 
        graph : dict()
            Heterogeneous graph to store all information from a `Solid` material
        '''
        self.API_KEY = API_KEY
        self.MATERIAL_ID = MATERIAL_ID
        self.graph = {}

    def structure_loader(self, conventional=True):
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
        structure = Solid(self.API_KEY, self.MATERIAL_ID).structure_loader(conventional=True)

        Nsites = len(structure.sites)

        nodes = {}
        ckeys = ['x','y','z']
        for i in ['atomic_number',*ckeys]:
            nodes[i] = [] 


        for i in range(Nsites):
            # atomic_num = element(str(structure.sites[i].specie)).atomic_number # specie is not a typo
            atomic_num = mendeleev_elements[str(structure.sites[i].specie)] # specie is not a typo

            nodes['atomic_number'].append(atomic_num)
            [ nodes[ckeys[k]].append(structure.sites[i].coords[k]) if not fractional else \
              nodes[ckeys[k]].append(structure.sites[i].frac_coords[k]) for k in range(3) ]

        self.graph = nodes

    def addProperty(self, key, value):
        '''Add a specific property label and its corresponding value or set of values. 

        Parameters
        ----------
        key: str
            Label / descriptor for the property to be added
        value: object
            Value(s) for the property to be added. May be any data type and/or be presented as a list()  
        '''
        self.graph[key] = value


class Crate:
    def __init__( self, API_KEY):
        '''Material class. To store and process material data

        Parameters
        ----------
        API_KEY : str
            secret API KEY used to run MPRester data requests
        MATERIALS : list(str)
            List potential materials (IDs) to place in the `graphs` dict / dataset
        graphs : dict(*MATERIALS: *dict())
            Dictionary of Heterogeneous graphs of extracted MATERIALS
        blacklist : list(str)
            List of materials (IDs) rejected from graphs dataset due to lack of information for an added Property. 
            For instance, if a material property to be analyzed is total_magnetization, the materials who lack 
            it will be dropped from the `graphs` variable and placed in this `blacklist`
        num_cores : int
            Number of [ multiprocessing ] cores to be used in parallel for extracting data from the 
            Materials Project server (in progress). 
        '''
        self.API_KEY = API_KEY
        self.MATERIALS = []
        self.graphs = {}
        self.blacklist = []
        self.num_cores = 1

    def topology(self, fractional=False):
        '''Same function as the Materials().topology() method but iterative, with the list of 
        specified material IDs in Crate().MATERIALS
        '''
        for i in self.MATERIALS:
            if i not in self.blacklist:
                if self.graphs.get(i) == None:
                    self.graphs.update({ i: {} })

                material = Solid(self.API_KEY,i)
                try:
                    material.topology()
                except:
                    pass

                if material.graph:
                    [self.graphs[i].update({ j: material.graph[j] }) for j in material.graph.keys()]
                else:
                    self.graphs.pop(i)
                    self.blacklist.append(i)


    def addProperty(self,func):
        '''Decorator used to add several properties iteratively from a REST function

        Usage Example:
        -----------------------
        materials_box = Crate(MY_SECRET_KEY)

        @materials_box.addProperty
        def REST_function_mag(CRATE):      
            with MPRester(api_key=MY_SECRET_KEY) as mpr:
                magnetism_doc = mpr.magnetism.get_data_by_id(CRATE)

            return magnetism_doc.num_unique_magnetic_sites, 'num unique mag. sites'
        
        REST_function_mag(materials_box)      
        '''
        def wrapper(self):
            for i in self.MATERIALS:
                if i not in self.blacklist:
                    if self.graphs.get(i) == None:
                        self.graphs.update({ i: {} })
                    'try if entries for the i material exist else eliminate material from dataset and put on blacklist'
                    x,y = 0,0      # patch for "UnboundLocalError: local variable 'x' referenced before assignment"
                    try:
                        x,y = func(i)
                    except:
                        pass

                    if x:
                        self.graphs[i].update({ y: x })
                    else:
                        self.graphs.pop(i)
                        self.blacklist.append(i)

        return wrapper



