from mp_api.client import MPRester
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.local_env import *

import matplotlib.pyplot as plt

# import stellargraph as sg
# try:
#     sg.utils.validate_notebook_version("1.2.1")
# except AttributeError:
#     raise ValueError(
#         f"This notebook requires StellarGraph version 1.2.1, but a different version {sg.__version__} is installed.  Please see <https://github.com/stellargraph/stellargraph/issues/1172>."
#     ) from None
    
from stellargraph import StellarGraph
import pandas
import json


mendeleev_elements = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118}

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
    structure0 = structure.copy()

    # important to use the conventional structure to ensure
    # that peaks are labelled with the conventional Miller indices
    sga = SpacegroupAnalyzer(structure)
    
    '''attempt to find bonds indirectly (must always specify atomic shell radius)'''
    # persites = structure.sites
    # NN = []
    # for i in range(len(persites)):
    #     NPS = structure.get_neighbors(structure[i],2.6)
    #     Nidcs = [ persites.index(n) for n in NPS ]
    #     NN.append(Nidcs)

    if conventional:
        # structure = sga.get_conventional_standard_structure()
        structure = sga.get_refined_structure()

    nodes = {}
    ckeys = list('xyz')
    topol_list = ['atomic_number',*ckeys]

    if edges: 
        topol_list.append('bond_edges')

    for i in topol_list:
        nodes[i] = []     

    
    for i in range(len(structure.sites)):
        # atomic_num = element(str(structure.sites[i].specie)).atomic_number # specie is not a typo
        atomic_num = mendeleev_elements[str(structure.sites[i].specie)] # specie is not a typo

        nodes['atomic_number'].append(atomic_num)
        [ nodes[ckeys[k]].append(structure.sites[i].coords[k]) if not fractional else \
            nodes[ckeys[k]].append(structure.sites[i].frac_coords[k]) for k in range(3) ]

        if edges:
            'connected sites (indices)'
            sites_for_atom_i = CrystalNN().get_nn_info(structure,i)        # This is a neighbor estimator; it does not display the actual bonds    

            nodes['bond_edges'].append( np.unique( [ l['site_index'] for l in sites_for_atom_i ] ).tolist() ) 

        # nodes['bond_edges'] = nodes['bond_edges'].astype('int')

    return nodes

import numpy as np
def view(topol_info, figsize=(8, 6), dpi=80, node = 400, edge = 5 ):
    '''Visualization of Material with designated topology
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
    colors = np.linspace(2**20,2**24,118,dtype='int') #divide color range into 118 colors (for the 118 chemical elements)


    atom = np.array(topol_info['atomic_number']).astype('int')

    xyz_arr = np.array([topol_info[j] for j in list('xyz')])


    for i in range(len(atom)):
        NNidcs = topol_info['bond_edges'][i]
        
        x0,y0,z0 = xyz_arr.T[i]

        x =  [topol_info['x'][k] for k in NNidcs ] 
        y =  [topol_info['y'][k] for k in NNidcs ] 
        z =  [topol_info['z'][k] for k in NNidcs ] 

        [ax.plot((x0,x[k]),(y0,y[k]),(z0,z[k]), linewidth =5, color="#"+hex(colors[atom[i]])[2:],alpha=0.5) for k in range(len(NNidcs))]

    for i in range(len(atom)):
        
        ax.scatter3D(*xyz_arr.T[i], s=400, c="#"+hex(colors[atom[i]])[2:])
    
    ax.set_facecolor('#5a5d70')
    # set_axes_equal(ax)           

    plt.axis('off')
    plt.show()




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

    def save(self,filename="graph.json"):
        '''Save Crate().graphs as a JSON file
        '''
        save_to_json(filename,self.graph)

    def load(self,filename="graph.json"):
        '''Load JSON file to Crate().graphs
        '''
        self.graph = load_from_json(filename)

    def topology(self, edges="bond_edges", fractional=False, conventional=True):
        '''helper function to load structural information of Material

        Parameters
        ----------
        fractional: bool
            if False, returns xyz coordinates else returns fractional (a,b,c) coordinates
        conventional: bool
            if True, returns conventional standard structure else returns primitive
        '''
        with MPRester(self.API_KEY) as mpr:
            # first retrieve the relevant structure
            structure = mpr.get_structure_by_material_id(self.MATERIAL_ID)
        
        self.graph = topol(structure, edges, fractional, conventional=True)


    def add(self, key, value):
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
        self.num_cores = 1

    def save(self,filename="graphs.json"):
        '''Save Crate().graphs as a JSON file
        '''
        save_to_json(filename,self.graphs)

    def load(self,filename="graphs.json"):
        '''Load JSON file to Crate().graphs
        '''
        self.graphs = load_from_json(filename)


    def queryAdd(self,edges="bond_edges",*args):
        '''
        Example Usage:
        queryAdd("structure","total_magnetization")
        '''
        print(list(args)[0])
        args_pass = list(args)[0]
        with MPRester(api_key=self.API_KEY) as mpr:
            mpresults = mpr.summary.search(material_ids=self.MATERIALS, fields=['material_id',*args_pass])
    
        
        for i in mpresults:
            mid = i.material_id
            self.graphs.update({ mid : {} })
            for j in args_pass:
                topol_list = ["atomic_number",*list('xyz')]

                if j == 'structure':
                    topol_nodes = topol(getattr(i,j), edges, fractional=False, conventional=True)
                    for k in topol_list:
                        self.graphs[mid].update({ k : topol_nodes[k]})

                    if edges:
                        self.graphs[mid].update({ edges : topol_nodes[edges] })
                        # self.graphs[mid].update({ edges : {} })
                        # L = len(topol_nodes['x'])
                        # [self.graphs[mid][edges].update({ l : topol_nodes[edges][l] }) for l in range(L)]

                else:
                    self.graphs[mid].update({j : getattr(i,j)})

    def remove_field(self,field='total_magnetization'):
        '''
        Parameters
        ----------
        field : str
            field / property to remove from graphs data
        '''
        try:
            for i in self.graphs.keys(): 
                    self.graphs[i].pop(field)
        except:
            print(field+" not in graphs")


from tensorflow.keras.optimizers import Adam
from materialsml.neuralnetwork import GraphNet

import stellargraph as sg
from stellargraph.mapper import PaddedGraphGenerator

from sklearn import model_selection
from tensorflow import keras

class Network:
    
    def __init__(self):
        self.nodes = []
        self.edges = []
        self.graphs = []
        self.graph_labels = []

        self.batch_size = 1
        self.epochs = 500
        self.fitted_model = []
        self.test_gen = []

    def importgraphs(self,graphs,label='band_gaps',edges='bond_edges'):
        '''Process graphs for machine learning. 

        Parameters
        ----------
        graphs : dict{*dict: *{}}
            a map of graphs with material properties of N-materials. 
        label : str
            label name of the property to use as the dependent Y variable in the machine
            learning method. Gets removed from graph keys and placed in a new `graphs_label` list.  
        '''
        self.nodes = graphs
        graph_labels = pandas.DataFrame([self.nodes[i].pop(label) for i in self.nodes.keys() ])

        if edges:
            EDGES_ALL_GRAPHS = [self.nodes[i].pop(edges) for i in self.nodes.keys() ]
            EDGES_ALL_MAPS = []
            
            for j in EDGES_ALL_GRAPHS:
                EDGES_I_GRAPH = j
                edges_map = {'source': [], 'target': []}
                source_idx = 0
                for k in EDGES_I_GRAPH:
                    for l in k:
                        edges_map['source'].append(source_idx)
                        edges_map['target'].append(l)
                    source_idx += 1
                EDGES_ALL_MAPS.append(edges_map)

            self.edges = EDGES_ALL_MAPS

        
        self.graph_labels  = graph_labels 

        stellar_graphs = []
        k_count = 0
        for i in self.nodes:
            n_atoms = len(self.nodes[i]['atomic_number'])
            df_nodes = pandas.DataFrame(self.nodes[i], index = [idx for idx in range(n_atoms)])

            if edges:
                df_edges = pandas.DataFrame(EDGES_ALL_MAPS[k_count])
                stellar_graphs.append(StellarGraph(df_nodes,df_edges))
            else:
                stellar_graphs.append(StellarGraph(df_nodes))

            k_count+=1

        self.graphs = stellar_graphs  

        return graph_labels, stellar_graphs

    def save(self,path='path/to/location'):
        
        self.fitted_model.save()

    def load(self,path='path/to/location'):

        model = keras.models.load_model('path/to/location')
        self.fitted_model = model
    
    def description(self):
        summary = pandas.DataFrame(
        [(g.number_of_nodes(), g.number_of_edges()) for g in self.graphs],
        columns=["nodes", "edges"],
        )
        
        print(summary.describe().round(1))

    
    def train(self, plot=True):
        model = GraphNet(self.graphs).build()
        model.compile( optimizer=Adam(lr=0.0001), loss="mean_squared_error")
        
        train_graphs, test_graphs = model_selection.train_test_split(self.graph_labels,\
             train_size=0.7, test_size=None)

        gen = PaddedGraphGenerator(graphs=self.graphs)

        train_gen = gen.flow(
            list(train_graphs.index - 1),
            targets=train_graphs.values,
            batch_size=self.batch_size,
            symmetric_normalization=False,
        )

        test_gen = gen.flow(
            list(test_graphs.index - 1),
            targets=test_graphs.values,
            batch_size=self.batch_size,
            symmetric_normalization=False,
        )


        history = model.fit(
            train_gen, epochs=self.epochs, verbose=1, validation_data=test_gen, shuffle=True,
        )

        self.fitted_model = model

        if plot:
            sg.utils.plot_history(history) 
            plt.show()

    def test(self): return self.fitted_model.evaluate(self.test_gen)

    def predict(self, test_graphs):

        predictions = self.fitted_model.predict(self.test_gen)
        ground_truth = test_graphs.values

        def reshapemodel(x): return x.reshape(x.shape[0])
        
        predictions = reshapemodel(predictions)
        ground_truth = reshapemodel(ground_truth)

        summary= pandas.DataFrame({'ground_truth':ground_truth, 'predictions':predictions})
        print(summary)
