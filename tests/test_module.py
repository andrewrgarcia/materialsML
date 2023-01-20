import materialsml as mml
from mp_api.client import MPRester
import random

'''replace line 6 with API key [ in string form ] from Materials Project site (https://materialsproject.org/api#api-key) 
and all mml.SECRET_KEY variable calls with SECRET_KEY  (see lines 9 - 10)'''
# SECRET_KEY = ''    

def tedst_materialprops():
    # material = mml.Solid(SECRET_KEY, 'mp-1103503')
    material = mml.Solid(mml.SECRET_KEY, 'mp-1103503')

    material.topology()
    print(material.graph)

    material.add('yoo',10203)

    print(material.graph)


def tedst_materialviewer():
    # material = mml.Solid(SECRET_KEY, 'mp-1103503')
    material = mml.Solid(mml.SECRET_KEY, 'mp-568671')
    material.topology()
    mml.view(material.graph,figsize=(11, 11))


def test_multiREST():
    crate = mml.Crate(mml.SECRET_KEY)
    crate.MATERIALS =  ['mp-'+str(i) for i in random.sample(range(1,154718), 200)]
    crate.queryAdd("bond_edges", "structure","total_magnetization","band_gap")

    print(crate.graphs)

    crate.save('graphs_test.json')

def tedst_graphs_loadsaveload():
    crate = mml.Crate(mml.SECRET_KEY)
    crate.load('graphs_test.json')
    crate.remove_field('band_gap')
    print(crate.graphs)
    crate.save('graphs_test2.json')
    crate.load('graphs_test2.json')
    crate.remove_field('band_gap')
    print(crate.graphs)



def learn(nodes_only, N_SET=1000):

    bond_edges = False if nodes_only else "bond_edges"

    crate = mml.Crate(mml.SECRET_KEY)
    crate.MATERIALS =  ['mp-'+str(i) for i in random.sample(range(1,154718), N_SET)]
    crate.queryAdd(bond_edges,["structure","total_magnetization"]) 
    # crate.save('graphs.json')

    net = mml.Network()
    net.batch_size = 50
    net.epochs = 100
    net.importgraphs(crate.graphs,'total_magnetization',bond_edges)
    net.description()
    print(net.nodes)
    print(net.edges)
    net.train()

def tedst_learn(): learn(True, 1000)

def tedst_learn_edges(): learn(False, 1000)

def tedst_learn_edges_smaller(): learn(False, 50)

