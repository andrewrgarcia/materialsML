import materialsml as mml

# replace below line with API key [ in string form ] from Materials Project site (https://materialsproject.org/api#api-key)
# SECRET_KEY = []

def test_materialprops():
    material = mml.Solid(mml.SECRET_KEY, 'mp-1103503')

    b = material.topology()
    print(b)



