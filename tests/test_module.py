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
            try:
                bandstructure = mpr.get_bandstructure_by_material_id(CRATE,line_mode=False)
            except:
                pass
            feature = bandstructure.get_band_gap()['energy'] if bandstructure else None
        return feature, 'band_gap_energy'

    @crate.addProperty
    def totalmag(CRATE):      
        '''magnetic properties data for Material'''
        with MPRester(api_key=mml.SECRET_KEY) as mpr:
            try:
                magnetism_doc = mpr.magnetism.get_data_by_id(CRATE)
            except:
                pass
            feature = magnetism_doc.total_magnetization if magnetism_doc else None
        return feature, 'total magnetization'

    @crate.addProperty
    def magmoms(CRATE):      
        '''magnetic properties data for Material'''
        with MPRester(api_key=mml.SECRET_KEY) as mpr:
            try:
                magnetism_doc = mpr.magnetism.get_data_by_id(CRATE)
            except:
                pass
            feature = magnetism_doc.magmoms if magnetism_doc else None 
        return feature, 'mag. moments'

    @crate.addProperty
    def num_magsites(CRATE):      
        '''magnetic properties data for Material'''
        with MPRester(api_key=mml.SECRET_KEY) as mpr:
            try:
                magnetism_doc = mpr.magnetism.get_data_by_id(CRATE)
            except:
                pass
            feature = magnetism_doc.num_unique_magnetic_sites if magnetism_doc else None
        return feature, 'num unique mag. sites'

    @crate.addProperty
    def thermo(CRATE):

        with MPRester(api_key=mml.SECRET_KEY) as mpr:
            # # for a single material
            # thermo_doc = mpr.thermo.get_data_by_id('mp-1103503')

            # # for many materials, it's much faster to use
            # # the `search` method, where additional material_ids can 
            # # be added to this list
            try:
                thermo_docs = mpr.thermo.search(material_ids=[CRATE])
            except:
                pass
            feature = thermo_docs[0].formation_energy_per_atom if thermo_docs else None
        return feature, 'form. energy'

    crate.MATERIALS = ['mp-'+str(i) for i in range(1103500,1103530)]

    print(crate.MATERIALS)

    crate.topology()

    bandgap(crate)
    totalmag(crate)
    # magmoms(crate) 
    # num_magsites(crate)      
    # thermo(crate)


    print(crate.graphs)
    print(crate.blacklist)
