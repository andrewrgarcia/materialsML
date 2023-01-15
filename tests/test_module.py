import materialsml as mml
from mp_api.client import MPRester as mpr
import pandas as pd

# replace below line with API key [ in string form ] from Materials Project site (https://materialsproject.org/api#api-key)
SECRET_KEY = ''    

def test_materialprops():
    # material = mml.Solid(SECRET_KEY, 'mp-1103503')
    material = mml.Solid(mml.SECRET_KEY, 'mp-1103503')

    # material.topology()
    # print(material.graph)

    material.add_factor('poo',10203)

    print(material.graph)

# def test_wrapper():
    
#     material = mml.Solid(mml.SECRET_KEY, 'mp-1103503')

#     @material.call      
#     def get(material_id):  
#         return mpr.get_charge_density_from_material_id(material_id=material_id)

#     print(get('mp-1103503' ))

    
    