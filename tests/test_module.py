import materialsml as mml
from mp_api.client import MPRester

# replace below line with API key [ in string form ] from Materials Project site (https://materialsproject.org/api#api-key)
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
    crate.MATERIALS =  ['mp-'+str(i) for i in range(10)]
    crate.queryAdd(["structure","total_magnetization","band_gap"])

    print(crate.graphs)

    crate.save()