
from os.path import join
import tempfile
import zipfile
from pymol import cmd, finish_launching
from pymol.cgo import *

finish_launching()

dirpath = None
    
def cgo_arrow(atom1='pk1', atom2='pk2', radius=0.07, gap=0.0, hlength=-1, hradius=-1, color='blue red', name=''):
    from chempy import cpv
    radius, gap = float(radius), float(gap)
    hlength, hradius = float(hlength), float(hradius)

    try:
        color1, color2 = color.split()
    except:
        color1 = color2 = color
    color1 = list(cmd.get_color_tuple(color1))
    color2 = list(cmd.get_color_tuple(color2))

    def get_coord(v):
        if not isinstance(v, str):
            return v
        if v.startswith('['):
            return cmd.safe_list_eval(v)
        return cmd.get_atom_coords(v)

    xyz1 = get_coord(atom1)
    xyz2 = get_coord(atom2)
    normal = cpv.normalize(cpv.sub(xyz1, xyz2))

    if hlength < 0:
        hlength = radius * 3.0
    if hradius < 0:
        hradius = hlength * 0.6

    if gap:
        diff = cpv.scale(normal, gap)
        xyz1 = cpv.sub(xyz1, diff)
        xyz2 = cpv.add(xyz2, diff)

    xyz3 = cpv.add(cpv.scale(normal, hlength), xyz2)

    obj = [cgo.CYLINDER] + xyz1 + xyz3 + [radius] + color1 + color2 +           [cgo.CONE] + xyz3 + xyz2 + [hradius, 0.0] + color2 + color2 +           [1.0, 0.0]
    return obj

    
dirpath = tempfile.mkdtemp()
zip_dir = 'hotspot_boundaries.zip'
with zipfile.ZipFile(zip_dir) as hs_zip:
    hs_zip.extractall(dirpath)

cmd.load(join(dirpath,"protein.pdb"), "protein")
cmd.show("cartoon", "protein")

if dirpath:
    f = join(dirpath, "0/label_threshold_30.0.mol2")
else:
    f = "0/label_threshold_30.0.mol2"

cmd.load(f, 'label_threshold_30.0')
cmd.hide('everything', 'label_threshold_30.0')
cmd.label("label_threshold_30.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [30.0]
gfiles = ['0/apolar.grd']
grids = ['apolar']
num = 0
surf_transparency = 0.2

if dirpath:
    gfiles = [join(dirpath, g) for g in gfiles]

for t in threshold_list:
    for i in range(len(grids)):
        try:
            cmd.load(r'%s'%(gfiles[i]), '%s_%s'%(grids[i], str(num)))
            cmd.isosurface('surface_%s_%s_%s'%(grids[i], t, num), '%s_%s'%(grids[i], num), t)
            cmd.set('transparency', surf_transparency, 'surface_%s_%s_%s'%(grids[i], t, num))
            cmd.color(colour_dict['%s'%(grids[i])], 'surface_%s_%s_%s'%(grids[i], t, num))
            cmd.group('threshold_%s'%(t), members = 'surface_%s_%s_%s'%(grids[i],t, num))
            cmd.group('threshold_%s' % (t), members='label_threshold_%s' % (t))
        except:
            continue



    try:
        cmd.group('hotspot_%s' % (num), members='threshold_%s' % (t))
    except:
        continue
    
    for g in grids:
        
        cmd.group('hotspot_%s' % (num), members='%s_%s' % (g,num))


cluster_dict = {"0":[], "0_arrows":[]}

cmd.load_cgo(cluster_dict["0"], "Features_0", 1)
cmd.load_cgo(cluster_dict["0_arrows"], "Arrows_0")
cmd.set("transparency", 0.2,"Features_0")
cmd.group("Pharmacophore_0", members="Features_0")
cmd.group("Pharmacophore_0", members="Arrows_0")

if dirpath:
    f = join(dirpath, "0/label_threshold_0.mol2")
else:
    f = "0/label_threshold_0.mol2"

cmd.load(f, 'label_threshold_0')
cmd.hide('everything', 'label_threshold_0')
cmd.label("label_threshold_0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_0', members= 'label_threshold_0')


if dirpath:
    f = join(dirpath, '0/mesh.grd')
else:
    f = '0/mesh.grd'
cmd.load(f, 'mesh_0')
cmd.isomesh("isomesh_0", "mesh_0", 0.9)
cmd.color("grey80", "isomesh_0")
cmd.set('transparency', 0.4, "isomesh_0")

cmd.group('hotspot_0', "isomesh_0")
cmd.group('hotspot_0', "mesh_0")

if dirpath:
    f = join(dirpath, "1/label_threshold_13.7.mol2")
else:
    f = "1/label_threshold_13.7.mol2"

cmd.load(f, 'label_threshold_13.7')
cmd.hide('everything', 'label_threshold_13.7')
cmd.label("label_threshold_13.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.7]
gfiles = ['1/apolar.grd']
grids = ['apolar']
num = 1
surf_transparency = 0.2

if dirpath:
    gfiles = [join(dirpath, g) for g in gfiles]

for t in threshold_list:
    for i in range(len(grids)):
        try:
            cmd.load(r'%s'%(gfiles[i]), '%s_%s'%(grids[i], str(num)))
            cmd.isosurface('surface_%s_%s_%s'%(grids[i], t, num), '%s_%s'%(grids[i], num), t)
            cmd.set('transparency', surf_transparency, 'surface_%s_%s_%s'%(grids[i], t, num))
            cmd.color(colour_dict['%s'%(grids[i])], 'surface_%s_%s_%s'%(grids[i], t, num))
            cmd.group('threshold_%s'%(t), members = 'surface_%s_%s_%s'%(grids[i],t, num))
            cmd.group('threshold_%s' % (t), members='label_threshold_%s' % (t))
        except:
            continue



    try:
        cmd.group('hotspot_%s' % (num), members='threshold_%s' % (t))
    except:
        continue
    
    for g in grids:
        
        cmd.group('hotspot_%s' % (num), members='%s_%s' % (g,num))


cluster_dict = {"16.375":[], "16.375_arrows":[]}

cluster_dict["16.375"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(57.0), float(24.0), float(1.0)]

cluster_dict["16.375_arrows"] += cgo_arrow([22.0,57.0,24.0], [22.472,54.511,23.015], color="blue red", name="Arrows_16.375_1")

cluster_dict["16.375"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(59.0), float(24.5), float(1.0)]

cluster_dict["16.375_arrows"] += cgo_arrow([27.0,59.0,24.5], [27.568,55.302,26.456], color="blue red", name="Arrows_16.375_2")

cluster_dict["16.375"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(38.0), float(61.5), float(22.0), float(1.0)]

cluster_dict["16.375_arrows"] += cgo_arrow([38.0,61.5,22.0], [38.852,61.496,19.583], color="blue red", name="Arrows_16.375_3")

cluster_dict["16.375"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(39.0), float(60.0), float(25.0), float(1.0)]

cluster_dict["16.375_arrows"] += cgo_arrow([39.0,60.0,25.0], [41.272,61.002,23.719], color="blue red", name="Arrows_16.375_4")

cluster_dict["16.375"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(31.9064534939), float(60.4108489838), float(23.840425297), float(1.0)]


cluster_dict["16.375"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(57.0), float(22.5), float(1.0)]

cluster_dict["16.375_arrows"] += cgo_arrow([18.5,57.0,22.5], [18.518,57.133,19.453], color="red blue", name="Arrows_16.375_5")

cluster_dict["16.375"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(57.0), float(26.5), float(1.0)]

cluster_dict["16.375_arrows"] += cgo_arrow([26.5,57.0,26.5], [27.568,55.302,26.456], color="red blue", name="Arrows_16.375_6")

cluster_dict["16.375"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(60.5), float(26.0), float(1.0)]

cluster_dict["16.375_arrows"] += cgo_arrow([26.5,60.5,26.0], [27.083,58.305,28.603], color="red blue", name="Arrows_16.375_7")

cluster_dict["16.375"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(33.5), float(61.0), float(23.5), float(1.0)]

cluster_dict["16.375_arrows"] += cgo_arrow([33.5,61.0,23.5], [30.398,63.43,24.707], color="red blue", name="Arrows_16.375_8")

cluster_dict["16.375"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(37.0), float(63.5), float(21.5), float(1.0)]

cluster_dict["16.375_arrows"] += cgo_arrow([37.0,63.5,21.5], [40.365,62.817,18.513], color="red blue", name="Arrows_16.375_9")

cmd.load_cgo(cluster_dict["16.375"], "Features_16.375", 1)
cmd.load_cgo(cluster_dict["16.375_arrows"], "Arrows_16.375")
cmd.set("transparency", 0.2,"Features_16.375")
cmd.group("Pharmacophore_16.375", members="Features_16.375")
cmd.group("Pharmacophore_16.375", members="Arrows_16.375")

if dirpath:
    f = join(dirpath, "1/label_threshold_16.375.mol2")
else:
    f = "1/label_threshold_16.375.mol2"

cmd.load(f, 'label_threshold_16.375')
cmd.hide('everything', 'label_threshold_16.375')
cmd.label("label_threshold_16.375", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.375', members= 'label_threshold_16.375')


if dirpath:
    f = join(dirpath, '1/mesh.grd')
else:
    f = '1/mesh.grd'
cmd.load(f, 'mesh_1')
cmd.isomesh("isomesh_1", "mesh_1", 0.9)
cmd.color("grey80", "isomesh_1")
cmd.set('transparency', 0.4, "isomesh_1")

cmd.group('hotspot_1', "isomesh_1")
cmd.group('hotspot_1', "mesh_1")

if dirpath:
    f = join(dirpath, "2/label_threshold_13.0.mol2")
else:
    f = "2/label_threshold_13.0.mol2"

cmd.load(f, 'label_threshold_13.0')
cmd.hide('everything', 'label_threshold_13.0')
cmd.label("label_threshold_13.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.0]
gfiles = ['2/apolar.grd']
grids = ['apolar']
num = 2
surf_transparency = 0.2

if dirpath:
    gfiles = [join(dirpath, g) for g in gfiles]

for t in threshold_list:
    for i in range(len(grids)):
        try:
            cmd.load(r'%s'%(gfiles[i]), '%s_%s'%(grids[i], str(num)))
            cmd.isosurface('surface_%s_%s_%s'%(grids[i], t, num), '%s_%s'%(grids[i], num), t)
            cmd.set('transparency', surf_transparency, 'surface_%s_%s_%s'%(grids[i], t, num))
            cmd.color(colour_dict['%s'%(grids[i])], 'surface_%s_%s_%s'%(grids[i], t, num))
            cmd.group('threshold_%s'%(t), members = 'surface_%s_%s_%s'%(grids[i],t, num))
            cmd.group('threshold_%s' % (t), members='label_threshold_%s' % (t))
        except:
            continue



    try:
        cmd.group('hotspot_%s' % (num), members='threshold_%s' % (t))
    except:
        continue
    
    for g in grids:
        
        cmd.group('hotspot_%s' % (num), members='%s_%s' % (g,num))


cluster_dict = {"16.0054998398":[], "16.0054998398_arrows":[]}

cluster_dict["16.0054998398"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(56.5), float(27.5), float(1.0)]

cluster_dict["16.0054998398_arrows"] += cgo_arrow([24.0,56.5,27.5], [22.078,55.253,29.415], color="blue red", name="Arrows_16.0054998398_1")

cluster_dict["16.0054998398"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(59.0), float(24.5), float(1.0)]

cluster_dict["16.0054998398_arrows"] += cgo_arrow([27.0,59.0,24.5], [27.568,55.302,26.456], color="blue red", name="Arrows_16.0054998398_2")

cluster_dict["16.0054998398"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(37.5), float(61.5), float(26.0), float(1.0)]

cluster_dict["16.0054998398_arrows"] += cgo_arrow([37.5,61.5,26.0], [38.29,60.093,28.64], color="blue red", name="Arrows_16.0054998398_3")

cluster_dict["16.0054998398"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(38.0), float(61.5), float(22.0), float(1.0)]

cluster_dict["16.0054998398_arrows"] += cgo_arrow([38.0,61.5,22.0], [38.852,61.496,19.583], color="blue red", name="Arrows_16.0054998398_4")

cluster_dict["16.0054998398"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(39.0), float(60.0), float(25.0), float(1.0)]

cluster_dict["16.0054998398_arrows"] += cgo_arrow([39.0,60.0,25.0], [41.272,61.002,23.719], color="blue red", name="Arrows_16.0054998398_5")

cluster_dict["16.0054998398"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(33.184386423), float(60.7622870601), float(23.8376706288), float(1.0)]


cluster_dict["16.0054998398"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(57.0), float(26.5), float(1.0)]

cluster_dict["16.0054998398_arrows"] += cgo_arrow([26.5,57.0,26.5], [27.568,55.302,26.456], color="red blue", name="Arrows_16.0054998398_6")

cluster_dict["16.0054998398"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(60.5), float(26.0), float(1.0)]

cluster_dict["16.0054998398_arrows"] += cgo_arrow([26.5,60.5,26.0], [27.083,58.305,28.603], color="red blue", name="Arrows_16.0054998398_7")

cluster_dict["16.0054998398"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(33.5), float(61.0), float(23.5), float(1.0)]

cluster_dict["16.0054998398_arrows"] += cgo_arrow([33.5,61.0,23.5], [30.398,63.43,24.707], color="red blue", name="Arrows_16.0054998398_8")

cluster_dict["16.0054998398"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(37.0), float(63.5), float(21.5), float(1.0)]

cluster_dict["16.0054998398_arrows"] += cgo_arrow([37.0,63.5,21.5], [40.365,62.817,18.513], color="red blue", name="Arrows_16.0054998398_9")

cmd.load_cgo(cluster_dict["16.0054998398"], "Features_16.0054998398", 1)
cmd.load_cgo(cluster_dict["16.0054998398_arrows"], "Arrows_16.0054998398")
cmd.set("transparency", 0.2,"Features_16.0054998398")
cmd.group("Pharmacophore_16.0054998398", members="Features_16.0054998398")
cmd.group("Pharmacophore_16.0054998398", members="Arrows_16.0054998398")

if dirpath:
    f = join(dirpath, "2/label_threshold_16.0054998398.mol2")
else:
    f = "2/label_threshold_16.0054998398.mol2"

cmd.load(f, 'label_threshold_16.0054998398')
cmd.hide('everything', 'label_threshold_16.0054998398')
cmd.label("label_threshold_16.0054998398", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.0054998398', members= 'label_threshold_16.0054998398')


if dirpath:
    f = join(dirpath, '2/mesh.grd')
else:
    f = '2/mesh.grd'
cmd.load(f, 'mesh_2')
cmd.isomesh("isomesh_2", "mesh_2", 0.9)
cmd.color("grey80", "isomesh_2")
cmd.set('transparency', 0.4, "isomesh_2")

cmd.group('hotspot_2', "isomesh_2")
cmd.group('hotspot_2', "mesh_2")

if dirpath:
    f = join(dirpath, "3/label_threshold_0.6.mol2")
else:
    f = "3/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
gfiles = ['3/apolar.grd']
grids = ['apolar']
num = 3
surf_transparency = 0.2

if dirpath:
    gfiles = [join(dirpath, g) for g in gfiles]

for t in threshold_list:
    for i in range(len(grids)):
        try:
            cmd.load(r'%s'%(gfiles[i]), '%s_%s'%(grids[i], str(num)))
            cmd.isosurface('surface_%s_%s_%s'%(grids[i], t, num), '%s_%s'%(grids[i], num), t)
            cmd.set('transparency', surf_transparency, 'surface_%s_%s_%s'%(grids[i], t, num))
            cmd.color(colour_dict['%s'%(grids[i])], 'surface_%s_%s_%s'%(grids[i], t, num))
            cmd.group('threshold_%s'%(t), members = 'surface_%s_%s_%s'%(grids[i],t, num))
            cmd.group('threshold_%s' % (t), members='label_threshold_%s' % (t))
        except:
            continue



    try:
        cmd.group('hotspot_%s' % (num), members='threshold_%s' % (t))
    except:
        continue
    
    for g in grids:
        
        cmd.group('hotspot_%s' % (num), members='%s_%s' % (g,num))


cluster_dict = {"13.4175000191":[], "13.4175000191_arrows":[]}

cluster_dict["13.4175000191"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(41.0), float(33.5), float(26.5), float(1.0)]

cluster_dict["13.4175000191_arrows"] += cgo_arrow([41.0,33.5,26.5], [43.273,34.781,27.189], color="blue red", name="Arrows_13.4175000191_1")

cluster_dict["13.4175000191"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(38.4855778153), float(33.8912117725), float(25.3050585768), float(1.0)]


cluster_dict["13.4175000191"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(33.0), float(24.0), float(1.0)]

cluster_dict["13.4175000191_arrows"] += cgo_arrow([36.5,33.0,24.0], [38.255,32.013,21.714], color="red blue", name="Arrows_13.4175000191_2")

cmd.load_cgo(cluster_dict["13.4175000191"], "Features_13.4175000191", 1)
cmd.load_cgo(cluster_dict["13.4175000191_arrows"], "Arrows_13.4175000191")
cmd.set("transparency", 0.2,"Features_13.4175000191")
cmd.group("Pharmacophore_13.4175000191", members="Features_13.4175000191")
cmd.group("Pharmacophore_13.4175000191", members="Arrows_13.4175000191")

if dirpath:
    f = join(dirpath, "3/label_threshold_13.4175000191.mol2")
else:
    f = "3/label_threshold_13.4175000191.mol2"

cmd.load(f, 'label_threshold_13.4175000191')
cmd.hide('everything', 'label_threshold_13.4175000191')
cmd.label("label_threshold_13.4175000191", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.4175000191', members= 'label_threshold_13.4175000191')


if dirpath:
    f = join(dirpath, '3/mesh.grd')
else:
    f = '3/mesh.grd'
cmd.load(f, 'mesh_3')
cmd.isomesh("isomesh_3", "mesh_3", 0.9)
cmd.color("grey80", "isomesh_3")
cmd.set('transparency', 0.4, "isomesh_3")

cmd.group('hotspot_3', "isomesh_3")
cmd.group('hotspot_3', "mesh_3")

if dirpath:
    f = join(dirpath, "4/label_threshold_9.7.mol2")
else:
    f = "4/label_threshold_9.7.mol2"

cmd.load(f, 'label_threshold_9.7')
cmd.hide('everything', 'label_threshold_9.7')
cmd.label("label_threshold_9.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [9.7]
gfiles = ['4/apolar.grd']
grids = ['apolar']
num = 4
surf_transparency = 0.2

if dirpath:
    gfiles = [join(dirpath, g) for g in gfiles]

for t in threshold_list:
    for i in range(len(grids)):
        try:
            cmd.load(r'%s'%(gfiles[i]), '%s_%s'%(grids[i], str(num)))
            cmd.isosurface('surface_%s_%s_%s'%(grids[i], t, num), '%s_%s'%(grids[i], num), t)
            cmd.set('transparency', surf_transparency, 'surface_%s_%s_%s'%(grids[i], t, num))
            cmd.color(colour_dict['%s'%(grids[i])], 'surface_%s_%s_%s'%(grids[i], t, num))
            cmd.group('threshold_%s'%(t), members = 'surface_%s_%s_%s'%(grids[i],t, num))
            cmd.group('threshold_%s' % (t), members='label_threshold_%s' % (t))
        except:
            continue



    try:
        cmd.group('hotspot_%s' % (num), members='threshold_%s' % (t))
    except:
        continue
    
    for g in grids:
        
        cmd.group('hotspot_%s' % (num), members='%s_%s' % (g,num))


cluster_dict = {"11.951499939":[], "11.951499939_arrows":[]}

cluster_dict["11.951499939"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(11.0), float(58.5), float(20.5), float(1.0)]

cluster_dict["11.951499939_arrows"] += cgo_arrow([11.0,58.5,20.5], [8.854,59.848,19.231], color="blue red", name="Arrows_11.951499939_1")

cluster_dict["11.951499939"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(11.0), float(58.5), float(20.5), float(1.0)]

cluster_dict["11.951499939_arrows"] += cgo_arrow([11.0,58.5,20.5], [8.854,59.848,19.231], color="blue red", name="Arrows_11.951499939_2")

cluster_dict["11.951499939"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(57.0), float(24.0), float(1.0)]

cluster_dict["11.951499939_arrows"] += cgo_arrow([22.0,57.0,24.0], [22.472,54.511,23.015], color="blue red", name="Arrows_11.951499939_3")

cluster_dict["11.951499939"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(57.0), float(24.0), float(1.0)]

cluster_dict["11.951499939_arrows"] += cgo_arrow([22.0,57.0,24.0], [22.472,54.511,23.015], color="blue red", name="Arrows_11.951499939_4")

cluster_dict["11.951499939"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(56.5), float(27.5), float(1.0)]

cluster_dict["11.951499939_arrows"] += cgo_arrow([23.5,56.5,27.5], [22.078,55.253,29.415], color="blue red", name="Arrows_11.951499939_5")

cluster_dict["11.951499939"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(13.5860904825), float(57.8339352864), float(22.0568033382), float(1.0)]


cluster_dict["11.951499939"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(21.8695180269), float(58.8909038818), float(25.1156210555), float(1.0)]


cluster_dict["11.951499939"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(57.0), float(22.5), float(1.0)]

cluster_dict["11.951499939_arrows"] += cgo_arrow([18.5,57.0,22.5], [18.518,57.133,19.453], color="red blue", name="Arrows_11.951499939_6")

cluster_dict["11.951499939"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(57.0), float(26.0), float(1.0)]

cluster_dict["11.951499939_arrows"] += cgo_arrow([25.0,57.0,26.0], [27.568,55.302,26.456], color="red blue", name="Arrows_11.951499939_7")

cluster_dict["11.951499939"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(57.0), float(26.0), float(1.0)]

cluster_dict["11.951499939_arrows"] += cgo_arrow([25.0,57.0,26.0], [27.568,55.302,26.456], color="red blue", name="Arrows_11.951499939_8")

cluster_dict["11.951499939"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(60.5), float(26.0), float(1.0)]

cluster_dict["11.951499939_arrows"] += cgo_arrow([25.0,60.5,26.0], [25.159,63.146,28.81], color="red blue", name="Arrows_11.951499939_9")

cmd.load_cgo(cluster_dict["11.951499939"], "Features_11.951499939", 1)
cmd.load_cgo(cluster_dict["11.951499939_arrows"], "Arrows_11.951499939")
cmd.set("transparency", 0.2,"Features_11.951499939")
cmd.group("Pharmacophore_11.951499939", members="Features_11.951499939")
cmd.group("Pharmacophore_11.951499939", members="Arrows_11.951499939")

if dirpath:
    f = join(dirpath, "4/label_threshold_11.951499939.mol2")
else:
    f = "4/label_threshold_11.951499939.mol2"

cmd.load(f, 'label_threshold_11.951499939')
cmd.hide('everything', 'label_threshold_11.951499939')
cmd.label("label_threshold_11.951499939", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.951499939', members= 'label_threshold_11.951499939')


if dirpath:
    f = join(dirpath, '4/mesh.grd')
else:
    f = '4/mesh.grd'
cmd.load(f, 'mesh_4')
cmd.isomesh("isomesh_4", "mesh_4", 0.9)
cmd.color("grey80", "isomesh_4")
cmd.set('transparency', 0.4, "isomesh_4")

cmd.group('hotspot_4', "isomesh_4")
cmd.group('hotspot_4', "mesh_4")

if dirpath:
    f = join(dirpath, "5/label_threshold_0.6.mol2")
else:
    f = "5/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
gfiles = ['5/apolar.grd']
grids = ['apolar']
num = 5
surf_transparency = 0.2

if dirpath:
    gfiles = [join(dirpath, g) for g in gfiles]

for t in threshold_list:
    for i in range(len(grids)):
        try:
            cmd.load(r'%s'%(gfiles[i]), '%s_%s'%(grids[i], str(num)))
            cmd.isosurface('surface_%s_%s_%s'%(grids[i], t, num), '%s_%s'%(grids[i], num), t)
            cmd.set('transparency', surf_transparency, 'surface_%s_%s_%s'%(grids[i], t, num))
            cmd.color(colour_dict['%s'%(grids[i])], 'surface_%s_%s_%s'%(grids[i], t, num))
            cmd.group('threshold_%s'%(t), members = 'surface_%s_%s_%s'%(grids[i],t, num))
            cmd.group('threshold_%s' % (t), members='label_threshold_%s' % (t))
        except:
            continue



    try:
        cmd.group('hotspot_%s' % (num), members='threshold_%s' % (t))
    except:
        continue
    
    for g in grids:
        
        cmd.group('hotspot_%s' % (num), members='%s_%s' % (g,num))


cluster_dict = {"11.3389997482":[], "11.3389997482_arrows":[]}

cluster_dict["11.3389997482"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(37.5), float(52.0), float(1.5), float(1.0)]

cluster_dict["11.3389997482_arrows"] += cgo_arrow([37.5,52.0,1.5], [34.864,52.782,1.167], color="blue red", name="Arrows_11.3389997482_1")

cluster_dict["11.3389997482"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(37.8139295838), float(52.9249970472), float(3.7280343179), float(1.0)]


cluster_dict["11.3389997482"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(54.5), float(7.0), float(1.0)]

cluster_dict["11.3389997482_arrows"] += cgo_arrow([36.5,54.5,7.0], [37.944,58.205,4.443], color="red blue", name="Arrows_11.3389997482_2")

cmd.load_cgo(cluster_dict["11.3389997482"], "Features_11.3389997482", 1)
cmd.load_cgo(cluster_dict["11.3389997482_arrows"], "Arrows_11.3389997482")
cmd.set("transparency", 0.2,"Features_11.3389997482")
cmd.group("Pharmacophore_11.3389997482", members="Features_11.3389997482")
cmd.group("Pharmacophore_11.3389997482", members="Arrows_11.3389997482")

if dirpath:
    f = join(dirpath, "5/label_threshold_11.3389997482.mol2")
else:
    f = "5/label_threshold_11.3389997482.mol2"

cmd.load(f, 'label_threshold_11.3389997482')
cmd.hide('everything', 'label_threshold_11.3389997482')
cmd.label("label_threshold_11.3389997482", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.3389997482', members= 'label_threshold_11.3389997482')


if dirpath:
    f = join(dirpath, '5/mesh.grd')
else:
    f = '5/mesh.grd'
cmd.load(f, 'mesh_5')
cmd.isomesh("isomesh_5", "mesh_5", 0.9)
cmd.color("grey80", "isomesh_5")
cmd.set('transparency', 0.4, "isomesh_5")

cmd.group('hotspot_5', "isomesh_5")
cmd.group('hotspot_5', "mesh_5")

if dirpath:
    f = join(dirpath, "6/label_threshold_0.6.mol2")
else:
    f = "6/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
gfiles = ['6/apolar.grd']
grids = ['apolar']
num = 6
surf_transparency = 0.2

if dirpath:
    gfiles = [join(dirpath, g) for g in gfiles]

for t in threshold_list:
    for i in range(len(grids)):
        try:
            cmd.load(r'%s'%(gfiles[i]), '%s_%s'%(grids[i], str(num)))
            cmd.isosurface('surface_%s_%s_%s'%(grids[i], t, num), '%s_%s'%(grids[i], num), t)
            cmd.set('transparency', surf_transparency, 'surface_%s_%s_%s'%(grids[i], t, num))
            cmd.color(colour_dict['%s'%(grids[i])], 'surface_%s_%s_%s'%(grids[i], t, num))
            cmd.group('threshold_%s'%(t), members = 'surface_%s_%s_%s'%(grids[i],t, num))
            cmd.group('threshold_%s' % (t), members='label_threshold_%s' % (t))
        except:
            continue



    try:
        cmd.group('hotspot_%s' % (num), members='threshold_%s' % (t))
    except:
        continue
    
    for g in grids:
        
        cmd.group('hotspot_%s' % (num), members='%s_%s' % (g,num))


cluster_dict = {"7.89200019836":[], "7.89200019836_arrows":[]}

cluster_dict["7.89200019836"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(9.16454024581), float(47.1977389199), float(20.6393422917), float(1.0)]


cmd.load_cgo(cluster_dict["7.89200019836"], "Features_7.89200019836", 1)
cmd.load_cgo(cluster_dict["7.89200019836_arrows"], "Arrows_7.89200019836")
cmd.set("transparency", 0.2,"Features_7.89200019836")
cmd.group("Pharmacophore_7.89200019836", members="Features_7.89200019836")
cmd.group("Pharmacophore_7.89200019836", members="Arrows_7.89200019836")

if dirpath:
    f = join(dirpath, "6/label_threshold_7.89200019836.mol2")
else:
    f = "6/label_threshold_7.89200019836.mol2"

cmd.load(f, 'label_threshold_7.89200019836')
cmd.hide('everything', 'label_threshold_7.89200019836')
cmd.label("label_threshold_7.89200019836", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_7.89200019836', members= 'label_threshold_7.89200019836')


if dirpath:
    f = join(dirpath, '6/mesh.grd')
else:
    f = '6/mesh.grd'
cmd.load(f, 'mesh_6')
cmd.isomesh("isomesh_6", "mesh_6", 0.9)
cmd.color("grey80", "isomesh_6")
cmd.set('transparency', 0.4, "isomesh_6")

cmd.group('hotspot_6', "isomesh_6")
cmd.group('hotspot_6', "mesh_6")

if dirpath:
    f = join(dirpath, "7/label_threshold_30.0.mol2")
else:
    f = "7/label_threshold_30.0.mol2"

cmd.load(f, 'label_threshold_30.0')
cmd.hide('everything', 'label_threshold_30.0')
cmd.label("label_threshold_30.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [30.0]
gfiles = ['7/apolar.grd']
grids = ['apolar']
num = 7
surf_transparency = 0.2

if dirpath:
    gfiles = [join(dirpath, g) for g in gfiles]

for t in threshold_list:
    for i in range(len(grids)):
        try:
            cmd.load(r'%s'%(gfiles[i]), '%s_%s'%(grids[i], str(num)))
            cmd.isosurface('surface_%s_%s_%s'%(grids[i], t, num), '%s_%s'%(grids[i], num), t)
            cmd.set('transparency', surf_transparency, 'surface_%s_%s_%s'%(grids[i], t, num))
            cmd.color(colour_dict['%s'%(grids[i])], 'surface_%s_%s_%s'%(grids[i], t, num))
            cmd.group('threshold_%s'%(t), members = 'surface_%s_%s_%s'%(grids[i],t, num))
            cmd.group('threshold_%s' % (t), members='label_threshold_%s' % (t))
        except:
            continue



    try:
        cmd.group('hotspot_%s' % (num), members='threshold_%s' % (t))
    except:
        continue
    
    for g in grids:
        
        cmd.group('hotspot_%s' % (num), members='%s_%s' % (g,num))


cluster_dict = {"0":[], "0_arrows":[]}

cmd.load_cgo(cluster_dict["0"], "Features_0", 1)
cmd.load_cgo(cluster_dict["0_arrows"], "Arrows_0")
cmd.set("transparency", 0.2,"Features_0")
cmd.group("Pharmacophore_0", members="Features_0")
cmd.group("Pharmacophore_0", members="Arrows_0")

if dirpath:
    f = join(dirpath, "7/label_threshold_0.mol2")
else:
    f = "7/label_threshold_0.mol2"

cmd.load(f, 'label_threshold_0')
cmd.hide('everything', 'label_threshold_0')
cmd.label("label_threshold_0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_0', members= 'label_threshold_0')


if dirpath:
    f = join(dirpath, '7/mesh.grd')
else:
    f = '7/mesh.grd'
cmd.load(f, 'mesh_7')
cmd.isomesh("isomesh_7", "mesh_7", 0.9)
cmd.color("grey80", "isomesh_7")
cmd.set('transparency', 0.4, "isomesh_7")

cmd.group('hotspot_7', "isomesh_7")
cmd.group('hotspot_7', "mesh_7")

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
