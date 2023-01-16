import materialsml as mml
from mp_api.client import MPRester

# replace below line with API key [ in string form ] from Materials Project site (https://materialsproject.org/api#api-key)
# SECRET_KEY = ''    

def test_materialprops():
    # material = mml.Solid(SECRET_KEY, 'mp-1103503')
    material = mml.Solid(mml.SECRET_KEY, 'mp-1103503')

    material.topology()
    print(material.graph)

    material.addProperty('poo',10203)

    print(material.graph)


def test_multiREST():
    crate = mml.Crate(mml.SECRET_KEY)

    @crate.addProperty
    def bandgap(CRATE):
        with MPRester(api_key=mml.SECRET_KEY) as mpr:
            bandstructure = mpr.get_bandstructure_by_material_id(CRATE,line_mode=False)
            band_gap = bandstructure.get_band_gap()['energy']

        return band_gap, 'band_gap_energy'

    @crate.addProperty
    def totalmag(CRATE):      
        '''magnetic properties data for Material'''
        with MPRester(api_key=mml.SECRET_KEY) as mpr:
            magnetism_doc = mpr.magnetism.get_data_by_id(CRATE)

        return magnetism_doc.total_magnetization, 'total magnetization'

    @crate.addProperty
    def magmoms(CRATE):      
        '''magnetic properties data for Material'''
        with MPRester(api_key=mml.SECRET_KEY) as mpr:
            magnetism_doc = mpr.magnetism.get_data_by_id(CRATE)

        return magnetism_doc.magmoms, 'mag. moments'


    crate.MATERIALS = ['mp-110350'+str(i) for i in range(9)]

    print(crate.MATERIALS)

    crate.topology()

    bandgap(crate)
    totalmag(crate)
    magmoms(crate) 


    print(crate.graphs)
    print(crate.blacklist)
