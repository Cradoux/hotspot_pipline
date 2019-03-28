
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
    f = join(dirpath, "0/label_threshold_8.9.mol2")
else:
    f = "0/label_threshold_8.9.mol2"

cmd.load(f, 'label_threshold_8.9')
cmd.hide('everything', 'label_threshold_8.9')
cmd.label("label_threshold_8.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [8.9]
gfiles = ['0/donor.grd', '0/apolar.grd', '0/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
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


cluster_dict = {"12.0150003433":[], "12.0150003433_arrows":[]}

cluster_dict["12.0150003433"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(36.0), float(58.5), float(31.0), float(1.0)]

cluster_dict["12.0150003433_arrows"] += cgo_arrow([36.0,58.5,31.0], [35.99,56.42,29.615], color="blue red", name="Arrows_12.0150003433_1")

cluster_dict["12.0150003433"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(37.0), float(59.0), float(28.5), float(1.0)]

cluster_dict["12.0150003433_arrows"] += cgo_arrow([37.0,59.0,28.5], [35.99,56.42,29.615], color="blue red", name="Arrows_12.0150003433_2")

cluster_dict["12.0150003433"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(39.5), float(59.5), float(22.0), float(1.0)]

cluster_dict["12.0150003433_arrows"] += cgo_arrow([39.5,59.5,22.0], [38.072,61.49,21.108], color="blue red", name="Arrows_12.0150003433_3")

cluster_dict["12.0150003433"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(38.4254313924), float(59.4446001373), float(25.0120266979), float(1.0)]


cluster_dict["12.0150003433"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(34.0), float(60.0), float(27.5), float(1.0)]

cluster_dict["12.0150003433_arrows"] += cgo_arrow([34.0,60.0,27.5], [35.99,56.42,29.615], color="red blue", name="Arrows_12.0150003433_4")

cluster_dict["12.0150003433"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(39.0), float(59.0), float(29.0), float(1.0)]

cluster_dict["12.0150003433_arrows"] += cgo_arrow([39.0,59.0,29.0], [42.048,60.184,28.406], color="red blue", name="Arrows_12.0150003433_5")

cluster_dict["12.0150003433"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(39.0), float(59.0), float(29.0), float(1.0)]

cluster_dict["12.0150003433_arrows"] += cgo_arrow([39.0,59.0,29.0], [42.048,60.184,28.406], color="red blue", name="Arrows_12.0150003433_6")

cmd.load_cgo(cluster_dict["12.0150003433"], "Features_12.0150003433", 1)
cmd.load_cgo(cluster_dict["12.0150003433_arrows"], "Arrows_12.0150003433")
cmd.set("transparency", 0.2,"Features_12.0150003433")
cmd.group("Pharmacophore_12.0150003433", members="Features_12.0150003433")
cmd.group("Pharmacophore_12.0150003433", members="Arrows_12.0150003433")

if dirpath:
    f = join(dirpath, "0/label_threshold_12.0150003433.mol2")
else:
    f = "0/label_threshold_12.0150003433.mol2"

cmd.load(f, 'label_threshold_12.0150003433')
cmd.hide('everything', 'label_threshold_12.0150003433')
cmd.label("label_threshold_12.0150003433", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.0150003433', members= 'label_threshold_12.0150003433')


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
    f = join(dirpath, "1/label_threshold_17.2.mol2")
else:
    f = "1/label_threshold_17.2.mol2"

cmd.load(f, 'label_threshold_17.2')
cmd.hide('everything', 'label_threshold_17.2')
cmd.label("label_threshold_17.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [17.2]
gfiles = ['1/donor.grd', '1/apolar.grd', '1/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
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


cluster_dict = {"18.9619998932":[], "18.9619998932_arrows":[]}

cluster_dict["18.9619998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(44.0), float(38.0), float(18.5), float(1.0)]

cluster_dict["18.9619998932_arrows"] += cgo_arrow([44.0,38.0,18.5], [45.249,40.474,18.915], color="blue red", name="Arrows_18.9619998932_1")

cluster_dict["18.9619998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(45.0), float(36.5), float(10.0), float(1.0)]

cluster_dict["18.9619998932_arrows"] += cgo_arrow([45.0,36.5,10.0], [46.98,35.513,7.652], color="blue red", name="Arrows_18.9619998932_2")

cluster_dict["18.9619998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(52.0), float(37.0), float(11.0), float(1.0)]

cluster_dict["18.9619998932_arrows"] += cgo_arrow([52.0,37.0,11.0], [54.453,36.097,11.926], color="blue red", name="Arrows_18.9619998932_3")

cluster_dict["18.9619998932"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(47.1533338249), float(34.8521724426), float(13.450621586), float(1.0)]


cluster_dict["18.9619998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(45.0), float(36.5), float(12.0), float(1.0)]

cluster_dict["18.9619998932_arrows"] += cgo_arrow([45.0,36.5,12.0], [45.729,33.684,7.995], color="red blue", name="Arrows_18.9619998932_4")

cluster_dict["18.9619998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(42.0), float(38.0), float(16.5), float(1.0)]

cluster_dict["18.9619998932_arrows"] += cgo_arrow([42.0,38.0,16.5], [41.879,39.772,14.796], color="red blue", name="Arrows_18.9619998932_5")

cluster_dict["18.9619998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(45.0), float(36.5), float(12.0), float(1.0)]

cluster_dict["18.9619998932_arrows"] += cgo_arrow([45.0,36.5,12.0], [45.729,33.684,7.995], color="red blue", name="Arrows_18.9619998932_6")

cluster_dict["18.9619998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(46.0), float(37.5), float(10.0), float(1.0)]

cluster_dict["18.9619998932_arrows"] += cgo_arrow([46.0,37.5,10.0], [45.729,33.684,7.995], color="red blue", name="Arrows_18.9619998932_7")

cluster_dict["18.9619998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(48.5), float(37.5), float(19.5), float(1.0)]

cluster_dict["18.9619998932_arrows"] += cgo_arrow([48.5,37.5,19.5], [46.853,41.904,18.83], color="red blue", name="Arrows_18.9619998932_8")

cluster_dict["18.9619998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(50.5), float(34.0), float(7.5), float(1.0)]

cluster_dict["18.9619998932_arrows"] += cgo_arrow([50.5,34.0,7.5], [48.598,30.873,6.922], color="red blue", name="Arrows_18.9619998932_9")

cmd.load_cgo(cluster_dict["18.9619998932"], "Features_18.9619998932", 1)
cmd.load_cgo(cluster_dict["18.9619998932_arrows"], "Arrows_18.9619998932")
cmd.set("transparency", 0.2,"Features_18.9619998932")
cmd.group("Pharmacophore_18.9619998932", members="Features_18.9619998932")
cmd.group("Pharmacophore_18.9619998932", members="Arrows_18.9619998932")

if dirpath:
    f = join(dirpath, "1/label_threshold_18.9619998932.mol2")
else:
    f = "1/label_threshold_18.9619998932.mol2"

cmd.load(f, 'label_threshold_18.9619998932')
cmd.hide('everything', 'label_threshold_18.9619998932')
cmd.label("label_threshold_18.9619998932", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_18.9619998932', members= 'label_threshold_18.9619998932')


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
    f = join(dirpath, "2/label_threshold_12.4.mol2")
else:
    f = "2/label_threshold_12.4.mol2"

cmd.load(f, 'label_threshold_12.4')
cmd.hide('everything', 'label_threshold_12.4')
cmd.label("label_threshold_12.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.4]
gfiles = ['2/donor.grd', '2/apolar.grd', '2/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
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


cluster_dict = {"15.9745001793":[], "15.9745001793_arrows":[]}

cluster_dict["15.9745001793"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(52.0), float(37.0), float(11.0), float(1.0)]

cluster_dict["15.9745001793_arrows"] += cgo_arrow([52.0,37.0,11.0], [54.453,36.097,11.926], color="blue red", name="Arrows_15.9745001793_1")

cluster_dict["15.9745001793"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(59.0), float(35.5), float(22.5), float(1.0)]

cluster_dict["15.9745001793_arrows"] += cgo_arrow([59.0,35.5,22.5], [59.992,34.123,24.744], color="blue red", name="Arrows_15.9745001793_2")

cluster_dict["15.9745001793"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(49.1343297262), float(33.7087407036), float(16.6115840646), float(1.0)]


cluster_dict["15.9745001793"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(58.8672044278), float(33.1752253398), float(20.2214873288), float(1.0)]


cluster_dict["15.9745001793"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(46.0), float(34.5), float(22.0), float(1.0)]

cluster_dict["15.9745001793_arrows"] += cgo_arrow([46.0,34.5,22.0], [44.211,32.161,23.217], color="red blue", name="Arrows_15.9745001793_3")

cluster_dict["15.9745001793"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(48.5), float(37.5), float(19.5), float(1.0)]

cluster_dict["15.9745001793_arrows"] += cgo_arrow([48.5,37.5,19.5], [46.853,41.904,18.83], color="red blue", name="Arrows_15.9745001793_4")

cluster_dict["15.9745001793"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(49.5), float(36.0), float(22.5), float(1.0)]

cluster_dict["15.9745001793_arrows"] += cgo_arrow([49.5,36.0,22.5], [51.835,36.961,24.694], color="red blue", name="Arrows_15.9745001793_5")

cluster_dict["15.9745001793"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(49.5), float(33.5), float(12.5), float(1.0)]


cluster_dict["15.9745001793"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(55.0), float(35.0), float(20.5), float(1.0)]

cluster_dict["15.9745001793_arrows"] += cgo_arrow([55.0,35.0,20.5], [54.016,36.622,22.042], color="red blue", name="Arrows_15.9745001793_6")

cmd.load_cgo(cluster_dict["15.9745001793"], "Features_15.9745001793", 1)
cmd.load_cgo(cluster_dict["15.9745001793_arrows"], "Arrows_15.9745001793")
cmd.set("transparency", 0.2,"Features_15.9745001793")
cmd.group("Pharmacophore_15.9745001793", members="Features_15.9745001793")
cmd.group("Pharmacophore_15.9745001793", members="Arrows_15.9745001793")

if dirpath:
    f = join(dirpath, "2/label_threshold_15.9745001793.mol2")
else:
    f = "2/label_threshold_15.9745001793.mol2"

cmd.load(f, 'label_threshold_15.9745001793')
cmd.hide('everything', 'label_threshold_15.9745001793')
cmd.label("label_threshold_15.9745001793", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.9745001793', members= 'label_threshold_15.9745001793')


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
    f = join(dirpath, "3/label_threshold_8.9.mol2")
else:
    f = "3/label_threshold_8.9.mol2"

cmd.load(f, 'label_threshold_8.9')
cmd.hide('everything', 'label_threshold_8.9')
cmd.label("label_threshold_8.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [8.9]
gfiles = ['3/donor.grd', '3/apolar.grd', '3/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
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


cluster_dict = {"12.0150003433":[], "12.0150003433_arrows":[]}

cluster_dict["12.0150003433"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(36.0), float(58.5), float(31.0), float(1.0)]

cluster_dict["12.0150003433_arrows"] += cgo_arrow([36.0,58.5,31.0], [35.99,56.42,29.615], color="blue red", name="Arrows_12.0150003433_1")

cluster_dict["12.0150003433"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(37.0), float(59.0), float(28.5), float(1.0)]

cluster_dict["12.0150003433_arrows"] += cgo_arrow([37.0,59.0,28.5], [35.99,56.42,29.615], color="blue red", name="Arrows_12.0150003433_2")

cluster_dict["12.0150003433"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(39.5), float(59.5), float(22.0), float(1.0)]

cluster_dict["12.0150003433_arrows"] += cgo_arrow([39.5,59.5,22.0], [38.072,61.49,21.108], color="blue red", name="Arrows_12.0150003433_3")

cluster_dict["12.0150003433"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(38.4516466836), float(59.4170319545), float(24.9568133932), float(1.0)]


cluster_dict["12.0150003433"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(34.0), float(60.0), float(27.5), float(1.0)]

cluster_dict["12.0150003433_arrows"] += cgo_arrow([34.0,60.0,27.5], [35.99,56.42,29.615], color="red blue", name="Arrows_12.0150003433_4")

cluster_dict["12.0150003433"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(39.0), float(59.0), float(29.0), float(1.0)]

cluster_dict["12.0150003433_arrows"] += cgo_arrow([39.0,59.0,29.0], [42.048,60.184,28.406], color="red blue", name="Arrows_12.0150003433_5")

cluster_dict["12.0150003433"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(39.0), float(59.0), float(29.0), float(1.0)]

cluster_dict["12.0150003433_arrows"] += cgo_arrow([39.0,59.0,29.0], [42.048,60.184,28.406], color="red blue", name="Arrows_12.0150003433_6")

cmd.load_cgo(cluster_dict["12.0150003433"], "Features_12.0150003433", 1)
cmd.load_cgo(cluster_dict["12.0150003433_arrows"], "Arrows_12.0150003433")
cmd.set("transparency", 0.2,"Features_12.0150003433")
cmd.group("Pharmacophore_12.0150003433", members="Features_12.0150003433")
cmd.group("Pharmacophore_12.0150003433", members="Arrows_12.0150003433")

if dirpath:
    f = join(dirpath, "3/label_threshold_12.0150003433.mol2")
else:
    f = "3/label_threshold_12.0150003433.mol2"

cmd.load(f, 'label_threshold_12.0150003433')
cmd.hide('everything', 'label_threshold_12.0150003433')
cmd.label("label_threshold_12.0150003433", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.0150003433', members= 'label_threshold_12.0150003433')


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
    f = join(dirpath, "4/label_threshold_0.6.mol2")
else:
    f = "4/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
gfiles = ['4/donor.grd', '4/apolar.grd', '4/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
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


cluster_dict = {"10.4899997711":[], "10.4899997711_arrows":[]}

cluster_dict["10.4899997711"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(33.4088337831), float(18.2519084606), float(15.5741442743), float(1.0)]


cmd.load_cgo(cluster_dict["10.4899997711"], "Features_10.4899997711", 1)
cmd.load_cgo(cluster_dict["10.4899997711_arrows"], "Arrows_10.4899997711")
cmd.set("transparency", 0.2,"Features_10.4899997711")
cmd.group("Pharmacophore_10.4899997711", members="Features_10.4899997711")
cmd.group("Pharmacophore_10.4899997711", members="Arrows_10.4899997711")

if dirpath:
    f = join(dirpath, "4/label_threshold_10.4899997711.mol2")
else:
    f = "4/label_threshold_10.4899997711.mol2"

cmd.load(f, 'label_threshold_10.4899997711')
cmd.hide('everything', 'label_threshold_10.4899997711')
cmd.label("label_threshold_10.4899997711", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.4899997711', members= 'label_threshold_10.4899997711')


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
gfiles = ['5/donor.grd', '5/apolar.grd', '5/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
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


cluster_dict = {"10.0909996033":[], "10.0909996033_arrows":[]}

cluster_dict["10.0909996033"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(30.2383650899), float(35.9382072446), float(4.67439393538), float(1.0)]


cmd.load_cgo(cluster_dict["10.0909996033"], "Features_10.0909996033", 1)
cmd.load_cgo(cluster_dict["10.0909996033_arrows"], "Arrows_10.0909996033")
cmd.set("transparency", 0.2,"Features_10.0909996033")
cmd.group("Pharmacophore_10.0909996033", members="Features_10.0909996033")
cmd.group("Pharmacophore_10.0909996033", members="Arrows_10.0909996033")

if dirpath:
    f = join(dirpath, "5/label_threshold_10.0909996033.mol2")
else:
    f = "5/label_threshold_10.0909996033.mol2"

cmd.load(f, 'label_threshold_10.0909996033')
cmd.hide('everything', 'label_threshold_10.0909996033')
cmd.label("label_threshold_10.0909996033", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.0909996033', members= 'label_threshold_10.0909996033')


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
    f = join(dirpath, "6/label_threshold_2.8.mol2")
else:
    f = "6/label_threshold_2.8.mol2"

cmd.load(f, 'label_threshold_2.8')
cmd.hide('everything', 'label_threshold_2.8')
cmd.label("label_threshold_2.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [2.8]
gfiles = ['6/donor.grd', '6/apolar.grd', '6/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
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


cluster_dict = {"8.88700008392":[], "8.88700008392_arrows":[]}

cluster_dict["8.88700008392"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(42.0), float(48.5), float(12.5), float(1.0)]

cluster_dict["8.88700008392_arrows"] += cgo_arrow([42.0,48.5,12.5], [42.318,48.694,15.82], color="blue red", name="Arrows_8.88700008392_1")

cluster_dict["8.88700008392"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(41.2515241523), float(46.1984215896), float(12.4385948032), float(1.0)]


cluster_dict["8.88700008392"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(41.5008460895), float(54.3108458561), float(17.3677414386), float(1.0)]


cluster_dict["8.88700008392"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(39.0), float(48.75), float(11.0), float(1.0)]


cluster_dict["8.88700008392"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(40.0), float(47.5), float(10.5), float(1.0)]

cluster_dict["8.88700008392_arrows"] += cgo_arrow([40.0,47.5,10.5], [39.209,46.76,7.629], color="red blue", name="Arrows_8.88700008392_2")

cmd.load_cgo(cluster_dict["8.88700008392"], "Features_8.88700008392", 1)
cmd.load_cgo(cluster_dict["8.88700008392_arrows"], "Arrows_8.88700008392")
cmd.set("transparency", 0.2,"Features_8.88700008392")
cmd.group("Pharmacophore_8.88700008392", members="Features_8.88700008392")
cmd.group("Pharmacophore_8.88700008392", members="Arrows_8.88700008392")

if dirpath:
    f = join(dirpath, "6/label_threshold_8.88700008392.mol2")
else:
    f = "6/label_threshold_8.88700008392.mol2"

cmd.load(f, 'label_threshold_8.88700008392')
cmd.hide('everything', 'label_threshold_8.88700008392')
cmd.label("label_threshold_8.88700008392", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_8.88700008392', members= 'label_threshold_8.88700008392')


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
    f = join(dirpath, "7/label_threshold_0.6.mol2")
else:
    f = "7/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
gfiles = ['7/donor.grd', '7/apolar.grd', '7/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
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


cluster_dict = {"8.47099971771":[], "8.47099971771_arrows":[]}

cluster_dict["8.47099971771"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.0), float(31.0), float(35.0), float(1.0)]

cluster_dict["8.47099971771_arrows"] += cgo_arrow([30.0,31.0,35.0], [32.2,32.335,34.153], color="blue red", name="Arrows_8.47099971771_1")

cluster_dict["8.47099971771"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.0298751826), float(31.0899416908), float(34.7099954518), float(1.0)]


cmd.load_cgo(cluster_dict["8.47099971771"], "Features_8.47099971771", 1)
cmd.load_cgo(cluster_dict["8.47099971771_arrows"], "Arrows_8.47099971771")
cmd.set("transparency", 0.2,"Features_8.47099971771")
cmd.group("Pharmacophore_8.47099971771", members="Features_8.47099971771")
cmd.group("Pharmacophore_8.47099971771", members="Arrows_8.47099971771")

if dirpath:
    f = join(dirpath, "7/label_threshold_8.47099971771.mol2")
else:
    f = "7/label_threshold_8.47099971771.mol2"

cmd.load(f, 'label_threshold_8.47099971771')
cmd.hide('everything', 'label_threshold_8.47099971771')
cmd.label("label_threshold_8.47099971771", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_8.47099971771', members= 'label_threshold_8.47099971771')


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
