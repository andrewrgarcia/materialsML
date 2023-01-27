from mp_api.client import MPRester

import matplotlib.pyplot as plt
from materialsml.periph import *

# import stellargraph as sg
# try:
#     sg.utils.validate_notebook_version("1.2.1")
# except AttributeError:
#     raise ValueError(
#         f"This notebook requires StellarGraph version 1.2.1, but a different version {sg.__version__} is installed.  Please see <https://github.com/stellargraph/stellargraph/issues/1172>."
#     ) from None
    
from stellargraph import StellarGraph
import pandas
import pyvista as pv
import pickle

import numpy as np

def view(topol_info, window_size=[1024, 768], node_size = 85, tube_radius = 0.1, tube_alpha=0.25, background_color='#696f85' ):
    '''Visualization of Material with designated topology
    Parameters
    ----------
    topol_info: dict() --> atomic_number, x,y,z
        dictionary containing the above information
    window_size: list([int,int])
        Size of PyVista Window default: [1024,768]
    node_size: int
        Size of node. Node represents an atom / cluster of material
    tube_radius: float
        radius of tube. Tube represents bond / edge of material
    tube_alpha: float
        Opacity of tube. 
    background_color: string / hexanumerical 
        Background color of PyVista window
    '''

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


    def Tube(tube_0,tube_f,rad):

        def make_points(crds1=tube_0,crds2=tube_f):
            """Helper to make XYZ points"""

            x0,y0,z0 = crds1
            xf,yf,zf = crds2

            theta = np.linspace(0, 2 * np.pi, 100)
            z = np.linspace(z0,zf, 100)
            x = np.linspace(x0, xf, 100)
            y = np.linspace(y0, yf, 100)
            return np.column_stack((x, y, z))


        def polyline_from_points(points):
            poly = pv.PolyData()
            poly.points = points
            the_cell = np.arange(0, len(points), dtype=np.int_)
            the_cell = np.insert(the_cell, 0, len(points))
            poly.lines = the_cell
            return poly

        points = make_points()
        polyline = polyline_from_points(points)
        polyline["scalars"] = np.arange(polyline.n_points)
        
        return polyline.tube(radius=rad)


    atom = topol_info['atom']

    xyz_arr = np.array([topol_info[j] for j in list('xyz')])


    pl = pv.Plotter(window_size = window_size)
    if background_color != "": 
        pl.background_color = background_color 


    for i in range(len(atom)):

        cloud = pv.wrap(xyz_arr.T[i])
        pl.add_mesh(cloud, color=jmol_colors[atom[i]], render_points_as_spheres=True,point_size=node_size)
    
    if 'bond_edges' in topol_info.keys():
        for i in range(len(atom)):
            NNidcs = topol_info['bond_edges'][i]
            
            x0,y0,z0 = xyz_arr.T[i]

            x =  [topol_info['x'][k] for k in NNidcs ] 
            y =  [topol_info['y'][k] for k in NNidcs ] 
            z =  [topol_info['z'][k] for k in NNidcs ] 

            # [ax.plot((x0,x[k]),(y0,y[k]),(z0,z[k]), linewidth =5,\
            #      color=jmol_colors[atom[i]], alpha=0.5) for k in range(len(NNidcs))]

            for k in range(len(NNidcs)):
                toob = Tube((x0,y0,z0),(x[k],y[k],z[k]),rad=tube_radius)   
                pl.add_mesh(toob,opacity=tube_alpha, smooth_shading=True)
    
        pl.remove_scalar_bar()
        
    pl.show()


class Solid:
    def __init__( self, API_KEY, MATERIAL_ID = 'mp-1103503'):
        '''Solid class. To store and process individual material data

        Parameters
        ----------
        API_KEY : str
            secret API KEY used to run MPRester data requests
        MATERIAL_ID : str
            ID of selected material 
        graph : dict()
            Homogeneous multi-feature graph to store all information from a `Solid` material
        structure: object
            structure object from pymatgen containing all information from the specified Material with MATERIAL_ID
        '''
        self.API_KEY = API_KEY
        self.MATERIAL_ID = MATERIAL_ID
        self.graph = {}
        self.structure = []

    def SolidSave(self,filename="solid.pickle"):
        '''Save Material information loaded from MP server'''
        with open(filename, 'wb') as f: 
            pickle.dump(self.structure, f)

    def SolidLoad(self,filename="solid.pickle"):
        '''Load Material information loaded from MP server'''
        with open(filename, 'rb') as f: 
            self.structure = pickle.load(f,encoding="latin1")

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
        if self.MATERIAL_ID != "":
            with MPRester(self.API_KEY) as mpr:
                # first retrieve the relevant structure
                self.structure = mpr.get_structure_by_material_id(self.MATERIAL_ID)
        
        self.graph = topol(self.structure, edges, fractional, conventional=True)


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
        '''Crate class. To handle multiple materials. 

        Parameters
        ----------
        API_KEY : str
            secret API KEY used to run MPRester data requests
        MATERIALS : list(str)
            List potential materials (IDs) to place in the `graphs` dict / dataset
        graphs : dict(*MATERIALS: *dict())
            Dictionary of Heterogeneous graphs of extracted MATERIALS
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
        # print(list(args)[0])
        # args_pass = list(args)[0]
        args_pass = args
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
        '''Network class. For Machine Learning. Currently uses a DeepGraphCNN graph neural network model to correlate
        material properties. 

        Parameters
        ----------
        nodes : pandas.Dataframe()
            secret API KEY used to run MPRester data requests
        edges : pandas.Dataframe()
            List potential materials (IDs) to place in the `graphs` dict / dataset
        graphs : dict(*MATERIALS: *dict())
            Dictionary of Homogeneous graphs of extracted MATERIALS
        graph_labels : pandas.Dataframe()
            list of quantitative values of property to be used as predicting "Y" variable or "label"
        
        
        batch_size : int
            size of batch for batch training process
        epochs : int
            number of training / learning epochs
        fitted_model : Keras model 
            fitted model after training. This variable is used to move model across functions in the same Network class. 
        test_gen : StellarGraph object 
            Fragment of graph data to be used for testing.      
        '''
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
