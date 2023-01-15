import materialsml as mml
from mp_api.client import MPRester
import pandas as pd

# replace below line with API key [ in string form ] from Materials Project site (https://materialsproject.org/api#api-key)
# SECRET_KEY = ''    

def test_materialprops():
    # material = mml.Solid(SECRET_KEY, 'mp-1103503')
    material = mml.Solid(mml.SECRET_KEY, 'mp-1103503')

    # material.topology()
    # print(material.graph)

    material.add_factor('poo',10203)

    print(material.graph)

def test_multirest():
    crate = mml.Crate(mml.SECRET_KEY)

    crate.MATERIALS = ['mp-110350'+str(i) for i in range(9)]
    print(crate.MATERIALS)

    @crate.addFeature
    def rester(CRATE):
        with MPRester(api_key=mml.SECRET_KEY) as mpr:

            bandstructure = mpr.get_bandstructure_by_material_id(CRATE,line_mode=False)
            
            band_gap = bandstructure.get_band_gap()['energy']

        return band_gap, 'band_gap_energy'


    rester(crate)
    print(crate.graphs)

    
    