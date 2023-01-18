import materialsml as mml
from mp_api.client import MPRester
import random

'''replace line 6 with API key [ in string form ] from Materials Project site (https://materialsproject.org/api#api-key) 
and all mml.SECRET_KEY variable calls with SECRET_KEY  (see lines 9 - 10)'''
# SECRET_KEY = ''    

def test_materialprops():
    # material = mml.Solid(SECRET_KEY, 'mp-1103503')
    material = mml.Solid(mml.SECRET_KEY, 'mp-1103503')

    material.topology()
    print(material.graph)

    material.add('yoo',10203)

    print(material.graph)


def test_multiREST():
    crate = mml.Crate(mml.SECRET_KEY)
    crate.MATERIALS =  ['mp-'+str(i) for i in random.sample(range(1,154718), 2000)]
    crate.queryAdd(["structure","total_magnetization","band_gap"])

    print(crate.graphs)

    crate.save('graphs_test.json')

def test_graphs_loadsaveload():
    crate = mml.Crate(mml.SECRET_KEY)
    crate.load('graphs_test.json')
    crate.remove_field('band_gap')
    print(crate.graphs)
    crate.save('graphs_test2.json')
    crate.load('graphs_test2.json')
    crate.remove_field('band_gap')
    print(crate.graphs)



def test_learn():
    crate = mml.Crate(mml.SECRET_KEY)
    crate.MATERIALS =  ['mp-'+str(i) for i in random.sample(range(1,154718), 2000)]
    crate.queryAdd(["structure","total_magnetization"])
    crate.save('graphs.json')

    net = mml.Network()
    net.batch_size = 50
    net.epochs = 100
    net.importgraphs(crate.graphs,'total_magnetization')
    net.description()
    net.train()