import materialsml as mml
from mp_api.client import MPRester

# replace below line with API key [ in string form ] from Materials Project site (https://materialsproject.org/api#api-key)
# SECRET_KEY = ''    

def te3st_materialprops():
    # material = mml.Solid(SECRET_KEY, 'mp-1103503')
    material = mml.Solid(mml.SECRET_KEY, 'mp-1103503')

    material.topology()
    print(material.graph)

    material.add('yoo',10203)

    print(material.graph)


def tes3t_multiREST():
    crate = mml.Crate(mml.SECRET_KEY)
    crate.MATERIALS =  ['mp-'+str(i) for i in range(100)]
    crate.queryAdd(["structure","total_magnetization","band_gap"])

    print(crate.graphs)

    crate.save('graphs_test.json')

def tes3t_graphs_loadsaveload():
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
    crate.load('graphs_test2.json')

    net = mml.Network()

    Y,X = net.importgraphs(crate.graphs,'total_magnetization')
    net.description()

    net.train()