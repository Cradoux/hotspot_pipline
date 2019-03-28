
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
    f = join(dirpath, "0/label_threshold_4.0.mol2")
else:
    f = "0/label_threshold_4.0.mol2"

cmd.load(f, 'label_threshold_4.0')
cmd.hide('everything', 'label_threshold_4.0')
cmd.label("label_threshold_4.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [4.0]
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


cluster_dict = {"8.35000038147":[], "8.35000038147_arrows":[]}

cluster_dict["8.35000038147"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(26.5345931319), float(28.1075210604), float(42.6599560213), float(1.0)]


cluster_dict["8.35000038147"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(32.5), float(43.0), float(1.0)]

cluster_dict["8.35000038147_arrows"] += cgo_arrow([27.0,32.5,43.0], [29.312,33.068,45.074], color="red blue", name="Arrows_8.35000038147_1")

cmd.load_cgo(cluster_dict["8.35000038147"], "Features_8.35000038147", 1)
cmd.load_cgo(cluster_dict["8.35000038147_arrows"], "Arrows_8.35000038147")
cmd.set("transparency", 0.2,"Features_8.35000038147")
cmd.group("Pharmacophore_8.35000038147", members="Features_8.35000038147")
cmd.group("Pharmacophore_8.35000038147", members="Arrows_8.35000038147")

if dirpath:
    f = join(dirpath, "0/label_threshold_8.35000038147.mol2")
else:
    f = "0/label_threshold_8.35000038147.mol2"

cmd.load(f, 'label_threshold_8.35000038147')
cmd.hide('everything', 'label_threshold_8.35000038147')
cmd.label("label_threshold_8.35000038147", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_8.35000038147', members= 'label_threshold_8.35000038147')


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
    f = join(dirpath, "1/label_threshold_13.1.mol2")
else:
    f = "1/label_threshold_13.1.mol2"

cmd.load(f, 'label_threshold_13.1')
cmd.hide('everything', 'label_threshold_13.1')
cmd.label("label_threshold_13.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.1]
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


cluster_dict = {"14.2119998932":[], "14.2119998932_arrows":[]}

cluster_dict["14.2119998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(36.5), float(28.0), float(1.0)]

cluster_dict["14.2119998932_arrows"] += cgo_arrow([22.5,36.5,28.0], [25.999,37.887,27.196], color="blue red", name="Arrows_14.2119998932_1")

cluster_dict["14.2119998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(20.5), float(27.5), float(24.0), float(1.0)]

cluster_dict["14.2119998932_arrows"] += cgo_arrow([20.5,27.5,24.0], [20.72,28.331,26.689], color="blue red", name="Arrows_14.2119998932_2")

cluster_dict["14.2119998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(34.5), float(27.5), float(1.0)]

cluster_dict["14.2119998932_arrows"] += cgo_arrow([24.5,34.5,27.5], [26.805,36.028,27.954], color="blue red", name="Arrows_14.2119998932_3")

cluster_dict["14.2119998932"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.4787451371), float(27.2943416197), float(20.829923707), float(1.0)]


cluster_dict["14.2119998932"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(19.9634884903), float(39.0419728048), float(30.2840081011), float(1.0)]


cluster_dict["14.2119998932"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.0213372812), float(33.6998813532), float(26.0375664653), float(1.0)]


cluster_dict["14.2119998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(34.5), float(23.5), float(1.0)]

cluster_dict["14.2119998932_arrows"] += cgo_arrow([26.0,34.5,23.5], [28.493,36.168,24.374], color="red blue", name="Arrows_14.2119998932_4")

cluster_dict["14.2119998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(34.5), float(23.5), float(1.0)]

cluster_dict["14.2119998932_arrows"] += cgo_arrow([26.0,34.5,23.5], [28.493,36.168,24.374], color="red blue", name="Arrows_14.2119998932_5")

cmd.load_cgo(cluster_dict["14.2119998932"], "Features_14.2119998932", 1)
cmd.load_cgo(cluster_dict["14.2119998932_arrows"], "Arrows_14.2119998932")
cmd.set("transparency", 0.2,"Features_14.2119998932")
cmd.group("Pharmacophore_14.2119998932", members="Features_14.2119998932")
cmd.group("Pharmacophore_14.2119998932", members="Arrows_14.2119998932")

if dirpath:
    f = join(dirpath, "1/label_threshold_14.2119998932.mol2")
else:
    f = "1/label_threshold_14.2119998932.mol2"

cmd.load(f, 'label_threshold_14.2119998932')
cmd.hide('everything', 'label_threshold_14.2119998932')
cmd.label("label_threshold_14.2119998932", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.2119998932', members= 'label_threshold_14.2119998932')


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
    f = join(dirpath, "2/label_threshold_12.6.mol2")
else:
    f = "2/label_threshold_12.6.mol2"

cmd.load(f, 'label_threshold_12.6')
cmd.hide('everything', 'label_threshold_12.6')
cmd.label("label_threshold_12.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.6]
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


cluster_dict = {"14.0439996719":[], "14.0439996719_arrows":[]}

cluster_dict["14.0439996719"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(35.5), float(30.5), float(1.0)]

cluster_dict["14.0439996719_arrows"] += cgo_arrow([21.0,35.5,30.5], [20.448,33.722,32.582], color="blue red", name="Arrows_14.0439996719_1")

cluster_dict["14.0439996719"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(20.5), float(27.5), float(24.0), float(1.0)]

cluster_dict["14.0439996719_arrows"] += cgo_arrow([20.5,27.5,24.0], [20.72,28.331,26.689], color="blue red", name="Arrows_14.0439996719_2")

cluster_dict["14.0439996719"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(34.5), float(27.5), float(1.0)]

cluster_dict["14.0439996719_arrows"] += cgo_arrow([24.5,34.5,27.5], [26.805,36.028,27.954], color="blue red", name="Arrows_14.0439996719_3")

cluster_dict["14.0439996719"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.5273798477), float(27.2212120158), float(20.4372334609), float(1.0)]


cluster_dict["14.0439996719"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(21.9612667538), float(33.5929027961), float(25.8923608804), float(1.0)]


cluster_dict["14.0439996719"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(34.5), float(23.5), float(1.0)]

cluster_dict["14.0439996719_arrows"] += cgo_arrow([26.5,34.5,23.5], [28.493,36.168,24.374], color="red blue", name="Arrows_14.0439996719_4")

cluster_dict["14.0439996719"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(34.5), float(23.5), float(1.0)]

cluster_dict["14.0439996719_arrows"] += cgo_arrow([26.5,34.5,23.5], [28.493,36.168,24.374], color="red blue", name="Arrows_14.0439996719_5")

cmd.load_cgo(cluster_dict["14.0439996719"], "Features_14.0439996719", 1)
cmd.load_cgo(cluster_dict["14.0439996719_arrows"], "Arrows_14.0439996719")
cmd.set("transparency", 0.2,"Features_14.0439996719")
cmd.group("Pharmacophore_14.0439996719", members="Features_14.0439996719")
cmd.group("Pharmacophore_14.0439996719", members="Arrows_14.0439996719")

if dirpath:
    f = join(dirpath, "2/label_threshold_14.0439996719.mol2")
else:
    f = "2/label_threshold_14.0439996719.mol2"

cmd.load(f, 'label_threshold_14.0439996719')
cmd.hide('everything', 'label_threshold_14.0439996719')
cmd.label("label_threshold_14.0439996719", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.0439996719', members= 'label_threshold_14.0439996719')


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
    f = join(dirpath, "3/label_threshold_12.4.mol2")
else:
    f = "3/label_threshold_12.4.mol2"

cmd.load(f, 'label_threshold_12.4')
cmd.hide('everything', 'label_threshold_12.4')
cmd.label("label_threshold_12.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.4]
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


cluster_dict = {"13.7580003738":[], "13.7580003738_arrows":[]}

cluster_dict["13.7580003738"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(20.0), float(34.5), float(29.5), float(1.0)]

cluster_dict["13.7580003738_arrows"] += cgo_arrow([20.0,34.5,29.5], [18.088,33.463,29.934], color="blue red", name="Arrows_13.7580003738_1")

cluster_dict["13.7580003738"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(20.5), float(27.5), float(24.0), float(1.0)]

cluster_dict["13.7580003738_arrows"] += cgo_arrow([20.5,27.5,24.0], [20.72,28.331,26.689], color="blue red", name="Arrows_13.7580003738_2")

cluster_dict["13.7580003738"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(34.0), float(27.0), float(1.0)]

cluster_dict["13.7580003738_arrows"] += cgo_arrow([24.0,34.0,27.0], [26.805,36.028,27.954], color="blue red", name="Arrows_13.7580003738_3")

cluster_dict["13.7580003738"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(16.5625), float(37.125), float(27.1875), float(1.0)]


cluster_dict["13.7580003738"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.3945437116), float(28.4503163434), float(21.5234964288), float(1.0)]


cluster_dict["13.7580003738"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(34.5), float(23.5), float(1.0)]

cluster_dict["13.7580003738_arrows"] += cgo_arrow([26.5,34.5,23.5], [28.493,36.168,24.374], color="red blue", name="Arrows_13.7580003738_4")

cluster_dict["13.7580003738"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(34.5), float(23.5), float(1.0)]

cluster_dict["13.7580003738_arrows"] += cgo_arrow([26.5,34.5,23.5], [28.493,36.168,24.374], color="red blue", name="Arrows_13.7580003738_5")

cmd.load_cgo(cluster_dict["13.7580003738"], "Features_13.7580003738", 1)
cmd.load_cgo(cluster_dict["13.7580003738_arrows"], "Arrows_13.7580003738")
cmd.set("transparency", 0.2,"Features_13.7580003738")
cmd.group("Pharmacophore_13.7580003738", members="Features_13.7580003738")
cmd.group("Pharmacophore_13.7580003738", members="Arrows_13.7580003738")

if dirpath:
    f = join(dirpath, "3/label_threshold_13.7580003738.mol2")
else:
    f = "3/label_threshold_13.7580003738.mol2"

cmd.load(f, 'label_threshold_13.7580003738')
cmd.hide('everything', 'label_threshold_13.7580003738')
cmd.label("label_threshold_13.7580003738", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.7580003738', members= 'label_threshold_13.7580003738')


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
    f = join(dirpath, "4/label_threshold_12.2.mol2")
else:
    f = "4/label_threshold_12.2.mol2"

cmd.load(f, 'label_threshold_12.2')
cmd.hide('everything', 'label_threshold_12.2')
cmd.label("label_threshold_12.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.2]
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


cluster_dict = {"13.3950004578":[], "13.3950004578_arrows":[]}

cluster_dict["13.3950004578"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(36.5), float(28.0), float(1.0)]

cluster_dict["13.3950004578_arrows"] += cgo_arrow([22.5,36.5,28.0], [25.999,37.887,27.196], color="blue red", name="Arrows_13.3950004578_1")

cluster_dict["13.3950004578"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(20.0), float(27.5), float(24.0), float(1.0)]

cluster_dict["13.3950004578_arrows"] += cgo_arrow([20.0,27.5,24.0], [20.72,28.331,26.689], color="blue red", name="Arrows_13.3950004578_2")

cluster_dict["13.3950004578"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(36.5), float(28.0), float(1.0)]

cluster_dict["13.3950004578_arrows"] += cgo_arrow([22.5,36.5,28.0], [25.999,37.887,27.196], color="blue red", name="Arrows_13.3950004578_3")

cluster_dict["13.3950004578"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(18.4665395503), float(39.2005501927), float(28.8162888078), float(1.0)]


cluster_dict["13.3950004578"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.6113939144), float(32.7707633539), float(24.6863599517), float(1.0)]


cluster_dict["13.3950004578"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(34.5), float(23.5), float(1.0)]

cluster_dict["13.3950004578_arrows"] += cgo_arrow([26.5,34.5,23.5], [28.493,36.168,24.374], color="red blue", name="Arrows_13.3950004578_4")

cluster_dict["13.3950004578"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(34.5), float(23.5), float(1.0)]

cluster_dict["13.3950004578_arrows"] += cgo_arrow([26.5,34.5,23.5], [28.493,36.168,24.374], color="red blue", name="Arrows_13.3950004578_5")

cmd.load_cgo(cluster_dict["13.3950004578"], "Features_13.3950004578", 1)
cmd.load_cgo(cluster_dict["13.3950004578_arrows"], "Arrows_13.3950004578")
cmd.set("transparency", 0.2,"Features_13.3950004578")
cmd.group("Pharmacophore_13.3950004578", members="Features_13.3950004578")
cmd.group("Pharmacophore_13.3950004578", members="Arrows_13.3950004578")

if dirpath:
    f = join(dirpath, "4/label_threshold_13.3950004578.mol2")
else:
    f = "4/label_threshold_13.3950004578.mol2"

cmd.load(f, 'label_threshold_13.3950004578')
cmd.hide('everything', 'label_threshold_13.3950004578')
cmd.label("label_threshold_13.3950004578", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.3950004578', members= 'label_threshold_13.3950004578')


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
    f = join(dirpath, "5/label_threshold_11.0.mol2")
else:
    f = "5/label_threshold_11.0.mol2"

cmd.load(f, 'label_threshold_11.0')
cmd.hide('everything', 'label_threshold_11.0')
cmd.label("label_threshold_11.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [11.0]
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


cluster_dict = {"12.3780002594":[], "12.3780002594_arrows":[]}

cluster_dict["12.3780002594"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(23.0594471251), float(30.1536810263), float(19.8764765761), float(1.0)]


cluster_dict["12.3780002594"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(30.5), float(25.0), float(1.0)]

cluster_dict["12.3780002594_arrows"] += cgo_arrow([26.0,30.5,25.0], [27.697,33.03,25.727], color="red blue", name="Arrows_12.3780002594_1")

cluster_dict["12.3780002594"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(34.5), float(23.5), float(1.0)]

cluster_dict["12.3780002594_arrows"] += cgo_arrow([26.5,34.5,23.5], [28.493,36.168,24.374], color="red blue", name="Arrows_12.3780002594_2")

cmd.load_cgo(cluster_dict["12.3780002594"], "Features_12.3780002594", 1)
cmd.load_cgo(cluster_dict["12.3780002594_arrows"], "Arrows_12.3780002594")
cmd.set("transparency", 0.2,"Features_12.3780002594")
cmd.group("Pharmacophore_12.3780002594", members="Features_12.3780002594")
cmd.group("Pharmacophore_12.3780002594", members="Arrows_12.3780002594")

if dirpath:
    f = join(dirpath, "5/label_threshold_12.3780002594.mol2")
else:
    f = "5/label_threshold_12.3780002594.mol2"

cmd.load(f, 'label_threshold_12.3780002594')
cmd.hide('everything', 'label_threshold_12.3780002594')
cmd.label("label_threshold_12.3780002594", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.3780002594', members= 'label_threshold_12.3780002594')


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


cluster_dict = {"11.0679998398":[], "11.0679998398_arrows":[]}

cluster_dict["11.0679998398"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(40.9151396481), float(34.3444615162), float(26.3043688816), float(1.0)]


cmd.load_cgo(cluster_dict["11.0679998398"], "Features_11.0679998398", 1)
cmd.load_cgo(cluster_dict["11.0679998398_arrows"], "Arrows_11.0679998398")
cmd.set("transparency", 0.2,"Features_11.0679998398")
cmd.group("Pharmacophore_11.0679998398", members="Features_11.0679998398")
cmd.group("Pharmacophore_11.0679998398", members="Arrows_11.0679998398")

if dirpath:
    f = join(dirpath, "6/label_threshold_11.0679998398.mol2")
else:
    f = "6/label_threshold_11.0679998398.mol2"

cmd.load(f, 'label_threshold_11.0679998398')
cmd.hide('everything', 'label_threshold_11.0679998398')
cmd.label("label_threshold_11.0679998398", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.0679998398', members= 'label_threshold_11.0679998398')


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
    f = join(dirpath, "7/label_threshold_3.9.mol2")
else:
    f = "7/label_threshold_3.9.mol2"

cmd.load(f, 'label_threshold_3.9')
cmd.hide('everything', 'label_threshold_3.9')
cmd.label("label_threshold_3.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [3.9]
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


cluster_dict = {"8.35000038147":[], "8.35000038147_arrows":[]}

cluster_dict["8.35000038147"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(26.5347379086), float(28.1147726138), float(42.6629290786), float(1.0)]


cluster_dict["8.35000038147"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(32.5), float(43.0), float(1.0)]

cluster_dict["8.35000038147_arrows"] += cgo_arrow([27.0,32.5,43.0], [29.312,33.068,45.074], color="red blue", name="Arrows_8.35000038147_1")

cmd.load_cgo(cluster_dict["8.35000038147"], "Features_8.35000038147", 1)
cmd.load_cgo(cluster_dict["8.35000038147_arrows"], "Arrows_8.35000038147")
cmd.set("transparency", 0.2,"Features_8.35000038147")
cmd.group("Pharmacophore_8.35000038147", members="Features_8.35000038147")
cmd.group("Pharmacophore_8.35000038147", members="Arrows_8.35000038147")

if dirpath:
    f = join(dirpath, "7/label_threshold_8.35000038147.mol2")
else:
    f = "7/label_threshold_8.35000038147.mol2"

cmd.load(f, 'label_threshold_8.35000038147')
cmd.hide('everything', 'label_threshold_8.35000038147')
cmd.label("label_threshold_8.35000038147", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_8.35000038147', members= 'label_threshold_8.35000038147')


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
