
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
    f = join(dirpath, "0/label_threshold_6.4.mol2")
else:
    f = "0/label_threshold_6.4.mol2"

cmd.load(f, 'label_threshold_6.4')
cmd.hide('everything', 'label_threshold_6.4')
cmd.label("label_threshold_6.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [6.4]
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


cluster_dict = {"12.9210000038":[], "12.9210000038_arrows":[]}

cluster_dict["12.9210000038"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(66.5), float(19.5), float(-11.0), float(1.0)]

cluster_dict["12.9210000038_arrows"] += cgo_arrow([66.5,19.5,-11.0], [67.519,20.793,-9.958], color="blue red", name="Arrows_12.9210000038_1")

cluster_dict["12.9210000038"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(66.5), float(19.5), float(-11.0), float(1.0)]

cluster_dict["12.9210000038_arrows"] += cgo_arrow([66.5,19.5,-11.0], [67.519,20.793,-9.958], color="blue red", name="Arrows_12.9210000038_2")

cluster_dict["12.9210000038"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(63.8554754203), float(18.1460997512), float(-11.7759993175), float(1.0)]


cluster_dict["12.9210000038"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(62.0), float(18.0), float(-10.0), float(1.0)]

cluster_dict["12.9210000038_arrows"] += cgo_arrow([62.0,18.0,-10.0], [64.661,14.982,-7.165], color="red blue", name="Arrows_12.9210000038_3")

cluster_dict["12.9210000038"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(62.0), float(18.0), float(-10.0), float(1.0)]

cluster_dict["12.9210000038_arrows"] += cgo_arrow([62.0,18.0,-10.0], [64.661,14.982,-7.165], color="red blue", name="Arrows_12.9210000038_4")

cluster_dict["12.9210000038"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(64.5), float(22.5), float(-10.5), float(1.0)]

cluster_dict["12.9210000038_arrows"] += cgo_arrow([64.5,22.5,-10.5], [66.733,24.228,-11.075], color="red blue", name="Arrows_12.9210000038_5")

cluster_dict["12.9210000038"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(67.5), float(19.0), float(-12.5), float(1.0)]

cluster_dict["12.9210000038_arrows"] += cgo_arrow([67.5,19.0,-12.5], [67.519,20.793,-9.958], color="red blue", name="Arrows_12.9210000038_6")

cluster_dict["12.9210000038"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(66.0), float(21.5), float(-12.5), float(1.0)]

cluster_dict["12.9210000038_arrows"] += cgo_arrow([66.0,21.5,-12.5], [67.519,20.793,-9.958], color="red blue", name="Arrows_12.9210000038_7")

cluster_dict["12.9210000038"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(67.5), float(19.0), float(-12.5), float(1.0)]

cluster_dict["12.9210000038_arrows"] += cgo_arrow([67.5,19.0,-12.5], [67.519,20.793,-9.958], color="red blue", name="Arrows_12.9210000038_8")

cluster_dict["12.9210000038"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(68.0), float(21.5), float(-12.5), float(1.0)]

cluster_dict["12.9210000038_arrows"] += cgo_arrow([68.0,21.5,-12.5], [67.519,20.793,-9.958], color="red blue", name="Arrows_12.9210000038_9")

cmd.load_cgo(cluster_dict["12.9210000038"], "Features_12.9210000038", 1)
cmd.load_cgo(cluster_dict["12.9210000038_arrows"], "Arrows_12.9210000038")
cmd.set("transparency", 0.2,"Features_12.9210000038")
cmd.group("Pharmacophore_12.9210000038", members="Features_12.9210000038")
cmd.group("Pharmacophore_12.9210000038", members="Arrows_12.9210000038")

if dirpath:
    f = join(dirpath, "0/label_threshold_12.9210000038.mol2")
else:
    f = "0/label_threshold_12.9210000038.mol2"

cmd.load(f, 'label_threshold_12.9210000038')
cmd.hide('everything', 'label_threshold_12.9210000038')
cmd.label("label_threshold_12.9210000038", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.9210000038', members= 'label_threshold_12.9210000038')


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
    f = join(dirpath, "1/label_threshold_12.5.mol2")
else:
    f = "1/label_threshold_12.5.mol2"

cmd.load(f, 'label_threshold_12.5')
cmd.hide('everything', 'label_threshold_12.5')
cmd.label("label_threshold_12.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.5]
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


cluster_dict = {"15.4989995956":[], "15.4989995956_arrows":[]}

cluster_dict["15.4989995956"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(98.5), float(9.0), float(22.0), float(1.0)]

cluster_dict["15.4989995956_arrows"] += cgo_arrow([98.5,9.0,22.0], [96.849,6.806,22.678], color="blue red", name="Arrows_15.4989995956_1")

cluster_dict["15.4989995956"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(99.5), float(11.0), float(16.0), float(1.0)]

cluster_dict["15.4989995956_arrows"] += cgo_arrow([99.5,11.0,16.0], [101.793,12.269,14.709], color="blue red", name="Arrows_15.4989995956_2")

cluster_dict["15.4989995956"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(105.0), float(18.5), float(28.0), float(1.0)]

cluster_dict["15.4989995956_arrows"] += cgo_arrow([105.0,18.5,28.0], [104.752,19.893,30.768], color="blue red", name="Arrows_15.4989995956_3")

cluster_dict["15.4989995956"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(103.5), float(9.5), float(19.5), float(1.0)]

cluster_dict["15.4989995956_arrows"] += cgo_arrow([103.5,9.5,19.5], [104.648,6.085,21.622], color="blue red", name="Arrows_15.4989995956_4")

cluster_dict["15.4989995956"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(106.5), float(20.0), float(29.0), float(1.0)]

cluster_dict["15.4989995956_arrows"] += cgo_arrow([106.5,20.0,29.0], [104.752,19.893,30.768], color="blue red", name="Arrows_15.4989995956_5")

cluster_dict["15.4989995956"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(110.0), float(9.5), float(20.0), float(1.0)]

cluster_dict["15.4989995956_arrows"] += cgo_arrow([110.0,9.5,20.0], [112.58,10.924,19.833], color="blue red", name="Arrows_15.4989995956_6")

cluster_dict["15.4989995956"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(100.17591383), float(11.6135728884), float(20.166939278), float(1.0)]


cluster_dict["15.4989995956"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(104.920083154), float(17.9933990355), float(25.3440617135), float(1.0)]


cluster_dict["15.4989995956"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(97.0), float(12.5), float(21.0), float(1.0)]

cluster_dict["15.4989995956_arrows"] += cgo_arrow([97.0,12.5,21.0], [94.861,13.295,22.615], color="red blue", name="Arrows_15.4989995956_7")

cluster_dict["15.4989995956"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(98.5), float(14.0), float(19.0), float(1.0)]

cluster_dict["15.4989995956_arrows"] += cgo_arrow([98.5,14.0,19.0], [102.548,16.345,18.531], color="red blue", name="Arrows_15.4989995956_8")

cluster_dict["15.4989995956"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(101.0), float(10.5), float(18.5), float(1.0)]

cluster_dict["15.4989995956_arrows"] += cgo_arrow([101.0,10.5,18.5], [103.102,11.6,16.39], color="red blue", name="Arrows_15.4989995956_9")

cmd.load_cgo(cluster_dict["15.4989995956"], "Features_15.4989995956", 1)
cmd.load_cgo(cluster_dict["15.4989995956_arrows"], "Arrows_15.4989995956")
cmd.set("transparency", 0.2,"Features_15.4989995956")
cmd.group("Pharmacophore_15.4989995956", members="Features_15.4989995956")
cmd.group("Pharmacophore_15.4989995956", members="Arrows_15.4989995956")

if dirpath:
    f = join(dirpath, "1/label_threshold_15.4989995956.mol2")
else:
    f = "1/label_threshold_15.4989995956.mol2"

cmd.load(f, 'label_threshold_15.4989995956')
cmd.hide('everything', 'label_threshold_15.4989995956')
cmd.label("label_threshold_15.4989995956", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.4989995956', members= 'label_threshold_15.4989995956')


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
    f = join(dirpath, "2/label_threshold_12.5.mol2")
else:
    f = "2/label_threshold_12.5.mol2"

cmd.load(f, 'label_threshold_12.5')
cmd.hide('everything', 'label_threshold_12.5')
cmd.label("label_threshold_12.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.5]
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


cluster_dict = {"15.4630002975":[], "15.4630002975_arrows":[]}

cluster_dict["15.4630002975"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(98.5), float(9.0), float(22.0), float(1.0)]

cluster_dict["15.4630002975_arrows"] += cgo_arrow([98.5,9.0,22.0], [96.849,6.806,22.678], color="blue red", name="Arrows_15.4630002975_1")

cluster_dict["15.4630002975"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(99.5), float(11.0), float(16.0), float(1.0)]

cluster_dict["15.4630002975_arrows"] += cgo_arrow([99.5,11.0,16.0], [101.793,12.269,14.709], color="blue red", name="Arrows_15.4630002975_2")

cluster_dict["15.4630002975"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(105.0), float(18.5), float(28.0), float(1.0)]

cluster_dict["15.4630002975_arrows"] += cgo_arrow([105.0,18.5,28.0], [104.752,19.893,30.768], color="blue red", name="Arrows_15.4630002975_3")

cluster_dict["15.4630002975"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(103.5), float(9.5), float(19.5), float(1.0)]

cluster_dict["15.4630002975_arrows"] += cgo_arrow([103.5,9.5,19.5], [104.648,6.085,21.622], color="blue red", name="Arrows_15.4630002975_4")

cluster_dict["15.4630002975"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(109.5), float(9.5), float(20.0), float(1.0)]

cluster_dict["15.4630002975_arrows"] += cgo_arrow([109.5,9.5,20.0], [109.895,6.591,19.309], color="blue red", name="Arrows_15.4630002975_5")

cluster_dict["15.4630002975"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(100.18136948), float(11.6033456184), float(20.1666930137), float(1.0)]


cluster_dict["15.4630002975"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(104.915702895), float(17.9881201888), float(25.3437333109), float(1.0)]


cluster_dict["15.4630002975"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(97.0), float(12.5), float(21.0), float(1.0)]

cluster_dict["15.4630002975_arrows"] += cgo_arrow([97.0,12.5,21.0], [94.861,13.295,22.615], color="red blue", name="Arrows_15.4630002975_6")

cluster_dict["15.4630002975"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(98.5), float(14.0), float(19.0), float(1.0)]

cluster_dict["15.4630002975_arrows"] += cgo_arrow([98.5,14.0,19.0], [102.548,16.345,18.531], color="red blue", name="Arrows_15.4630002975_7")

cluster_dict["15.4630002975"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(101.0), float(10.5), float(18.5), float(1.0)]

cluster_dict["15.4630002975_arrows"] += cgo_arrow([101.0,10.5,18.5], [103.102,11.6,16.39], color="red blue", name="Arrows_15.4630002975_8")

cluster_dict["15.4630002975"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(105.5), float(11.0), float(18.0), float(1.0)]

cluster_dict["15.4630002975_arrows"] += cgo_arrow([105.5,11.0,18.0], [103.102,11.6,16.39], color="red blue", name="Arrows_15.4630002975_9")

cmd.load_cgo(cluster_dict["15.4630002975"], "Features_15.4630002975", 1)
cmd.load_cgo(cluster_dict["15.4630002975_arrows"], "Arrows_15.4630002975")
cmd.set("transparency", 0.2,"Features_15.4630002975")
cmd.group("Pharmacophore_15.4630002975", members="Features_15.4630002975")
cmd.group("Pharmacophore_15.4630002975", members="Arrows_15.4630002975")

if dirpath:
    f = join(dirpath, "2/label_threshold_15.4630002975.mol2")
else:
    f = "2/label_threshold_15.4630002975.mol2"

cmd.load(f, 'label_threshold_15.4630002975')
cmd.hide('everything', 'label_threshold_15.4630002975')
cmd.label("label_threshold_15.4630002975", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.4630002975', members= 'label_threshold_15.4630002975')


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
    f = join(dirpath, "3/label_threshold_12.0.mol2")
else:
    f = "3/label_threshold_12.0.mol2"

cmd.load(f, 'label_threshold_12.0')
cmd.hide('everything', 'label_threshold_12.0')
cmd.label("label_threshold_12.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.0]
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


cluster_dict = {"15.3529996872":[], "15.3529996872_arrows":[]}

cluster_dict["15.3529996872"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(98.5), float(9.0), float(22.0), float(1.0)]

cluster_dict["15.3529996872_arrows"] += cgo_arrow([98.5,9.0,22.0], [96.849,6.806,22.678], color="blue red", name="Arrows_15.3529996872_1")

cluster_dict["15.3529996872"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(99.5), float(11.0), float(16.0), float(1.0)]

cluster_dict["15.3529996872_arrows"] += cgo_arrow([99.5,11.0,16.0], [101.793,12.269,14.709], color="blue red", name="Arrows_15.3529996872_2")

cluster_dict["15.3529996872"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(105.0), float(18.5), float(28.0), float(1.0)]

cluster_dict["15.3529996872_arrows"] += cgo_arrow([105.0,18.5,28.0], [104.752,19.893,30.768], color="blue red", name="Arrows_15.3529996872_3")

cluster_dict["15.3529996872"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(103.5), float(9.5), float(19.5), float(1.0)]

cluster_dict["15.3529996872_arrows"] += cgo_arrow([103.5,9.5,19.5], [104.648,6.085,21.622], color="blue red", name="Arrows_15.3529996872_4")

cluster_dict["15.3529996872"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(106.5), float(20.0), float(29.0), float(1.0)]

cluster_dict["15.3529996872_arrows"] += cgo_arrow([106.5,20.0,29.0], [104.752,19.893,30.768], color="blue red", name="Arrows_15.3529996872_5")

cluster_dict["15.3529996872"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(110.0), float(9.5), float(20.0), float(1.0)]

cluster_dict["15.3529996872_arrows"] += cgo_arrow([110.0,9.5,20.0], [112.58,10.924,19.833], color="blue red", name="Arrows_15.3529996872_6")

cluster_dict["15.3529996872"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(100.46487736), float(11.8321234356), float(20.3782815744), float(1.0)]


cluster_dict["15.3529996872"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(104.700627602), float(17.9744173841), float(25.426892134), float(1.0)]


cluster_dict["15.3529996872"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(97.0), float(12.5), float(21.0), float(1.0)]

cluster_dict["15.3529996872_arrows"] += cgo_arrow([97.0,12.5,21.0], [94.861,13.295,22.615], color="red blue", name="Arrows_15.3529996872_7")

cluster_dict["15.3529996872"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(98.5), float(14.0), float(19.0), float(1.0)]

cluster_dict["15.3529996872_arrows"] += cgo_arrow([98.5,14.0,19.0], [102.548,16.345,18.531], color="red blue", name="Arrows_15.3529996872_8")

cluster_dict["15.3529996872"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(101.0), float(10.5), float(18.5), float(1.0)]

cluster_dict["15.3529996872_arrows"] += cgo_arrow([101.0,10.5,18.5], [103.102,11.6,16.39], color="red blue", name="Arrows_15.3529996872_9")

cmd.load_cgo(cluster_dict["15.3529996872"], "Features_15.3529996872", 1)
cmd.load_cgo(cluster_dict["15.3529996872_arrows"], "Arrows_15.3529996872")
cmd.set("transparency", 0.2,"Features_15.3529996872")
cmd.group("Pharmacophore_15.3529996872", members="Features_15.3529996872")
cmd.group("Pharmacophore_15.3529996872", members="Arrows_15.3529996872")

if dirpath:
    f = join(dirpath, "3/label_threshold_15.3529996872.mol2")
else:
    f = "3/label_threshold_15.3529996872.mol2"

cmd.load(f, 'label_threshold_15.3529996872')
cmd.hide('everything', 'label_threshold_15.3529996872')
cmd.label("label_threshold_15.3529996872", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.3529996872', members= 'label_threshold_15.3529996872')


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
    f = join(dirpath, "4/label_threshold_10.1.mol2")
else:
    f = "4/label_threshold_10.1.mol2"

cmd.load(f, 'label_threshold_10.1')
cmd.hide('everything', 'label_threshold_10.1')
cmd.label("label_threshold_10.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.1]
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


cluster_dict = {"14.9029998779":[], "14.9029998779_arrows":[]}

cluster_dict["14.9029998779"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(67.0), float(4.0), float(6.0), float(1.0)]

cluster_dict["14.9029998779_arrows"] += cgo_arrow([67.0,4.0,6.0], [66.325,1.506,4.625], color="blue red", name="Arrows_14.9029998779_1")

cluster_dict["14.9029998779"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(72.747090928), float(3.45417296345), float(6.02288293499), float(1.0)]


cluster_dict["14.9029998779"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(67.7603550296), float(6.54733727811), float(-0.289940828402), float(1.0)]


cluster_dict["14.9029998779"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(66.0), float(7.0), float(1.5), float(1.0)]

cluster_dict["14.9029998779_arrows"] += cgo_arrow([66.0,7.0,1.5], [65.169,5.124,3.011], color="red blue", name="Arrows_14.9029998779_2")

cluster_dict["14.9029998779"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(69.5), float(6.0), float(5.0), float(1.0)]

cluster_dict["14.9029998779_arrows"] += cgo_arrow([69.5,6.0,5.0], [69.484,8.203,2.825], color="red blue", name="Arrows_14.9029998779_3")

cluster_dict["14.9029998779"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(71.5), float(7.0), float(7.0), float(1.0)]

cluster_dict["14.9029998779_arrows"] += cgo_arrow([71.5,7.0,7.0], [68.285,8.284,7.197], color="red blue", name="Arrows_14.9029998779_4")

cluster_dict["14.9029998779"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(72.0), float(0.5), float(6.5), float(1.0)]

cluster_dict["14.9029998779_arrows"] += cgo_arrow([72.0,0.5,6.5], [70.104,-0.114,3.575], color="red blue", name="Arrows_14.9029998779_5")

cluster_dict["14.9029998779"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(72.0), float(6.5), float(4.0), float(1.0)]

cluster_dict["14.9029998779_arrows"] += cgo_arrow([72.0,6.5,4.0], [70.289,6.272,1.864], color="red blue", name="Arrows_14.9029998779_6")

cmd.load_cgo(cluster_dict["14.9029998779"], "Features_14.9029998779", 1)
cmd.load_cgo(cluster_dict["14.9029998779_arrows"], "Arrows_14.9029998779")
cmd.set("transparency", 0.2,"Features_14.9029998779")
cmd.group("Pharmacophore_14.9029998779", members="Features_14.9029998779")
cmd.group("Pharmacophore_14.9029998779", members="Arrows_14.9029998779")

if dirpath:
    f = join(dirpath, "4/label_threshold_14.9029998779.mol2")
else:
    f = "4/label_threshold_14.9029998779.mol2"

cmd.load(f, 'label_threshold_14.9029998779')
cmd.hide('everything', 'label_threshold_14.9029998779')
cmd.label("label_threshold_14.9029998779", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.9029998779', members= 'label_threshold_14.9029998779')


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
    f = join(dirpath, "5/label_threshold_9.7.mol2")
else:
    f = "5/label_threshold_9.7.mol2"

cmd.load(f, 'label_threshold_9.7')
cmd.hide('everything', 'label_threshold_9.7')
cmd.label("label_threshold_9.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [9.7]
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


cluster_dict = {"14.7770004272":[], "14.7770004272_arrows":[]}

cluster_dict["14.7770004272"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(67.0), float(4.0), float(6.0), float(1.0)]

cluster_dict["14.7770004272_arrows"] += cgo_arrow([67.0,4.0,6.0], [66.325,1.506,4.625], color="blue red", name="Arrows_14.7770004272_1")

cluster_dict["14.7770004272"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(72.422759653), float(3.70317795079), float(5.79647075325), float(1.0)]


cluster_dict["14.7770004272"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(68.1864242304), float(5.82557792023), float(-1.19941436652), float(1.0)]


cluster_dict["14.7770004272"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(66.0), float(7.0), float(1.5), float(1.0)]

cluster_dict["14.7770004272_arrows"] += cgo_arrow([66.0,7.0,1.5], [65.169,5.124,3.011], color="red blue", name="Arrows_14.7770004272_2")

cluster_dict["14.7770004272"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(69.5), float(6.0), float(5.0), float(1.0)]

cluster_dict["14.7770004272_arrows"] += cgo_arrow([69.5,6.0,5.0], [69.484,8.203,2.825], color="red blue", name="Arrows_14.7770004272_3")

cluster_dict["14.7770004272"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(71.5), float(7.0), float(7.0), float(1.0)]

cluster_dict["14.7770004272_arrows"] += cgo_arrow([71.5,7.0,7.0], [68.285,8.284,7.197], color="red blue", name="Arrows_14.7770004272_4")

cluster_dict["14.7770004272"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(71.0), float(5.0), float(0.5), float(1.0)]

cluster_dict["14.7770004272_arrows"] += cgo_arrow([71.0,5.0,0.5], [70.289,6.272,1.864], color="red blue", name="Arrows_14.7770004272_5")

cluster_dict["14.7770004272"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(72.0), float(0.5), float(6.5), float(1.0)]

cluster_dict["14.7770004272_arrows"] += cgo_arrow([72.0,0.5,6.5], [70.104,-0.114,3.575], color="red blue", name="Arrows_14.7770004272_6")

cluster_dict["14.7770004272"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(72.0), float(6.5), float(4.0), float(1.0)]

cluster_dict["14.7770004272_arrows"] += cgo_arrow([72.0,6.5,4.0], [70.289,6.272,1.864], color="red blue", name="Arrows_14.7770004272_7")

cmd.load_cgo(cluster_dict["14.7770004272"], "Features_14.7770004272", 1)
cmd.load_cgo(cluster_dict["14.7770004272_arrows"], "Arrows_14.7770004272")
cmd.set("transparency", 0.2,"Features_14.7770004272")
cmd.group("Pharmacophore_14.7770004272", members="Features_14.7770004272")
cmd.group("Pharmacophore_14.7770004272", members="Arrows_14.7770004272")

if dirpath:
    f = join(dirpath, "5/label_threshold_14.7770004272.mol2")
else:
    f = "5/label_threshold_14.7770004272.mol2"

cmd.load(f, 'label_threshold_14.7770004272')
cmd.hide('everything', 'label_threshold_14.7770004272')
cmd.label("label_threshold_14.7770004272", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.7770004272', members= 'label_threshold_14.7770004272')


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
    f = join(dirpath, "6/label_threshold_6.4.mol2")
else:
    f = "6/label_threshold_6.4.mol2"

cmd.load(f, 'label_threshold_6.4')
cmd.hide('everything', 'label_threshold_6.4')
cmd.label("label_threshold_6.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [6.4]
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


cluster_dict = {"12.9210000038":[], "12.9210000038_arrows":[]}

cluster_dict["12.9210000038"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(66.5), float(19.5), float(-11.0), float(1.0)]

cluster_dict["12.9210000038_arrows"] += cgo_arrow([66.5,19.5,-11.0], [67.519,20.793,-9.958], color="blue red", name="Arrows_12.9210000038_1")

cluster_dict["12.9210000038"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(66.5), float(19.5), float(-11.0), float(1.0)]

cluster_dict["12.9210000038_arrows"] += cgo_arrow([66.5,19.5,-11.0], [67.519,20.793,-9.958], color="blue red", name="Arrows_12.9210000038_2")

cluster_dict["12.9210000038"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(63.8554754203), float(18.1460997512), float(-11.7759993175), float(1.0)]


cluster_dict["12.9210000038"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(62.0), float(18.0), float(-10.0), float(1.0)]

cluster_dict["12.9210000038_arrows"] += cgo_arrow([62.0,18.0,-10.0], [64.661,14.982,-7.165], color="red blue", name="Arrows_12.9210000038_3")

cluster_dict["12.9210000038"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(62.0), float(18.0), float(-10.0), float(1.0)]

cluster_dict["12.9210000038_arrows"] += cgo_arrow([62.0,18.0,-10.0], [64.661,14.982,-7.165], color="red blue", name="Arrows_12.9210000038_4")

cluster_dict["12.9210000038"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(64.5), float(22.5), float(-10.5), float(1.0)]

cluster_dict["12.9210000038_arrows"] += cgo_arrow([64.5,22.5,-10.5], [66.733,24.228,-11.075], color="red blue", name="Arrows_12.9210000038_5")

cluster_dict["12.9210000038"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(67.5), float(19.0), float(-12.5), float(1.0)]

cluster_dict["12.9210000038_arrows"] += cgo_arrow([67.5,19.0,-12.5], [67.519,20.793,-9.958], color="red blue", name="Arrows_12.9210000038_6")

cluster_dict["12.9210000038"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(66.0), float(21.5), float(-12.5), float(1.0)]

cluster_dict["12.9210000038_arrows"] += cgo_arrow([66.0,21.5,-12.5], [67.519,20.793,-9.958], color="red blue", name="Arrows_12.9210000038_7")

cluster_dict["12.9210000038"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(67.5), float(19.0), float(-12.5), float(1.0)]

cluster_dict["12.9210000038_arrows"] += cgo_arrow([67.5,19.0,-12.5], [67.519,20.793,-9.958], color="red blue", name="Arrows_12.9210000038_8")

cluster_dict["12.9210000038"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(68.0), float(21.5), float(-12.5), float(1.0)]

cluster_dict["12.9210000038_arrows"] += cgo_arrow([68.0,21.5,-12.5], [67.519,20.793,-9.958], color="red blue", name="Arrows_12.9210000038_9")

cmd.load_cgo(cluster_dict["12.9210000038"], "Features_12.9210000038", 1)
cmd.load_cgo(cluster_dict["12.9210000038_arrows"], "Arrows_12.9210000038")
cmd.set("transparency", 0.2,"Features_12.9210000038")
cmd.group("Pharmacophore_12.9210000038", members="Features_12.9210000038")
cmd.group("Pharmacophore_12.9210000038", members="Arrows_12.9210000038")

if dirpath:
    f = join(dirpath, "6/label_threshold_12.9210000038.mol2")
else:
    f = "6/label_threshold_12.9210000038.mol2"

cmd.load(f, 'label_threshold_12.9210000038')
cmd.hide('everything', 'label_threshold_12.9210000038')
cmd.label("label_threshold_12.9210000038", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.9210000038', members= 'label_threshold_12.9210000038')


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
    f = join(dirpath, "7/label_threshold_8.8.mol2")
else:
    f = "7/label_threshold_8.8.mol2"

cmd.load(f, 'label_threshold_8.8')
cmd.hide('everything', 'label_threshold_8.8')
cmd.label("label_threshold_8.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [8.8]
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


cluster_dict = {"12.3555002213":[], "12.3555002213_arrows":[]}

cluster_dict["12.3555002213"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(96.0), float(-5.5), float(27.0), float(1.0)]

cluster_dict["12.3555002213_arrows"] += cgo_arrow([96.0,-5.5,27.0], [98.727,-5.395,26.249], color="blue red", name="Arrows_12.3555002213_1")

cluster_dict["12.3555002213"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(97.5), float(-2.5), float(30.0), float(1.0)]

cluster_dict["12.3555002213_arrows"] += cgo_arrow([97.5,-2.5,30.0], [99.521,-4.414,30.008], color="blue red", name="Arrows_12.3555002213_2")

cluster_dict["12.3555002213"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(93.8361372559), float(-4.74785319081), float(29.1308946625), float(1.0)]


cluster_dict["12.3555002213"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(99.4365263609), float(2.99654588085), float(22.8230194206), float(1.0)]


cluster_dict["12.3555002213"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(91.5), float(-3.0), float(31.5), float(1.0)]

cluster_dict["12.3555002213_arrows"] += cgo_arrow([91.5,-3.0,31.5], [88.939,-4.027,31.726], color="red blue", name="Arrows_12.3555002213_3")

cluster_dict["12.3555002213"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(95.5), float(-0.5), float(30.0), float(1.0)]

cluster_dict["12.3555002213_arrows"] += cgo_arrow([95.5,-0.5,30.0], [95.361,2.265,28.552], color="red blue", name="Arrows_12.3555002213_4")

cluster_dict["12.3555002213"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(98.5), float(5.0), float(20.0), float(1.0)]

cluster_dict["12.3555002213_arrows"] += cgo_arrow([98.5,5.0,20.0], [96.932,3.526,18.588], color="red blue", name="Arrows_12.3555002213_5")

cmd.load_cgo(cluster_dict["12.3555002213"], "Features_12.3555002213", 1)
cmd.load_cgo(cluster_dict["12.3555002213_arrows"], "Arrows_12.3555002213")
cmd.set("transparency", 0.2,"Features_12.3555002213")
cmd.group("Pharmacophore_12.3555002213", members="Features_12.3555002213")
cmd.group("Pharmacophore_12.3555002213", members="Arrows_12.3555002213")

if dirpath:
    f = join(dirpath, "7/label_threshold_12.3555002213.mol2")
else:
    f = "7/label_threshold_12.3555002213.mol2"

cmd.load(f, 'label_threshold_12.3555002213')
cmd.hide('everything', 'label_threshold_12.3555002213')
cmd.label("label_threshold_12.3555002213", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.3555002213', members= 'label_threshold_12.3555002213')


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

if dirpath:
    f = join(dirpath, "8/label_threshold_0.6.mol2")
else:
    f = "8/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
gfiles = ['8/donor.grd', '8/apolar.grd', '8/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 8
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


cluster_dict = {"11.6190004349":[], "11.6190004349_arrows":[]}

cluster_dict["11.6190004349"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(53.5), float(16.0), float(16.5), float(1.0)]

cluster_dict["11.6190004349_arrows"] += cgo_arrow([53.5,16.0,16.5], [54.536,17.934,14.833], color="blue red", name="Arrows_11.6190004349_1")

cluster_dict["11.6190004349"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(51.3000091607), float(15.1412010291), float(17.7909066766), float(1.0)]


cluster_dict["11.6190004349"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(51.5), float(14.0), float(20.0), float(1.0)]

cluster_dict["11.6190004349_arrows"] += cgo_arrow([51.5,14.0,20.0], [49.799,16.798,22.539], color="red blue", name="Arrows_11.6190004349_2")

cmd.load_cgo(cluster_dict["11.6190004349"], "Features_11.6190004349", 1)
cmd.load_cgo(cluster_dict["11.6190004349_arrows"], "Arrows_11.6190004349")
cmd.set("transparency", 0.2,"Features_11.6190004349")
cmd.group("Pharmacophore_11.6190004349", members="Features_11.6190004349")
cmd.group("Pharmacophore_11.6190004349", members="Arrows_11.6190004349")

if dirpath:
    f = join(dirpath, "8/label_threshold_11.6190004349.mol2")
else:
    f = "8/label_threshold_11.6190004349.mol2"

cmd.load(f, 'label_threshold_11.6190004349')
cmd.hide('everything', 'label_threshold_11.6190004349')
cmd.label("label_threshold_11.6190004349", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.6190004349', members= 'label_threshold_11.6190004349')


if dirpath:
    f = join(dirpath, '8/mesh.grd')
else:
    f = '8/mesh.grd'
cmd.load(f, 'mesh_8')
cmd.isomesh("isomesh_8", "mesh_8", 0.9)
cmd.color("grey80", "isomesh_8")
cmd.set('transparency', 0.4, "isomesh_8")

cmd.group('hotspot_8', "isomesh_8")
cmd.group('hotspot_8', "mesh_8")

if dirpath:
    f = join(dirpath, "9/label_threshold_7.5.mol2")
else:
    f = "9/label_threshold_7.5.mol2"

cmd.load(f, 'label_threshold_7.5')
cmd.hide('everything', 'label_threshold_7.5')
cmd.label("label_threshold_7.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [7.5]
gfiles = ['9/donor.grd', '9/apolar.grd', '9/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 9
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


cluster_dict = {"10.0369997025":[], "10.0369997025_arrows":[]}

cluster_dict["10.0369997025"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(93.0), float(19.0), float(6.0), float(1.0)]

cluster_dict["10.0369997025_arrows"] += cgo_arrow([93.0,19.0,6.0], [90.356,19.867,5.397], color="blue red", name="Arrows_10.0369997025_1")

cluster_dict["10.0369997025"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(104.5), float(21.5), float(6.5), float(1.0)]

cluster_dict["10.0369997025_arrows"] += cgo_arrow([104.5,21.5,6.5], [101.899,22.062,7.793], color="blue red", name="Arrows_10.0369997025_2")

cluster_dict["10.0369997025"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(97.1106598035), float(17.533384688), float(2.77630118909), float(1.0)]


cluster_dict["10.0369997025"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(103.570339483), float(21.4529868796), float(3.27305253523), float(1.0)]


cmd.load_cgo(cluster_dict["10.0369997025"], "Features_10.0369997025", 1)
cmd.load_cgo(cluster_dict["10.0369997025_arrows"], "Arrows_10.0369997025")
cmd.set("transparency", 0.2,"Features_10.0369997025")
cmd.group("Pharmacophore_10.0369997025", members="Features_10.0369997025")
cmd.group("Pharmacophore_10.0369997025", members="Arrows_10.0369997025")

if dirpath:
    f = join(dirpath, "9/label_threshold_10.0369997025.mol2")
else:
    f = "9/label_threshold_10.0369997025.mol2"

cmd.load(f, 'label_threshold_10.0369997025')
cmd.hide('everything', 'label_threshold_10.0369997025')
cmd.label("label_threshold_10.0369997025", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.0369997025', members= 'label_threshold_10.0369997025')


if dirpath:
    f = join(dirpath, '9/mesh.grd')
else:
    f = '9/mesh.grd'
cmd.load(f, 'mesh_9')
cmd.isomesh("isomesh_9", "mesh_9", 0.9)
cmd.color("grey80", "isomesh_9")
cmd.set('transparency', 0.4, "isomesh_9")

cmd.group('hotspot_9', "isomesh_9")
cmd.group('hotspot_9', "mesh_9")

if dirpath:
    f = join(dirpath, "10/label_threshold_7.2.mol2")
else:
    f = "10/label_threshold_7.2.mol2"

cmd.load(f, 'label_threshold_7.2')
cmd.hide('everything', 'label_threshold_7.2')
cmd.label("label_threshold_7.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [7.2]
gfiles = ['10/donor.grd', '10/apolar.grd', '10/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 10
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


cluster_dict = {"10.0359997749":[], "10.0359997749_arrows":[]}

cluster_dict["10.0359997749"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(93.0), float(19.0), float(6.0), float(1.0)]

cluster_dict["10.0359997749_arrows"] += cgo_arrow([93.0,19.0,6.0], [90.356,19.867,5.397], color="blue red", name="Arrows_10.0359997749_1")

cluster_dict["10.0359997749"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(104.5), float(21.5), float(6.5), float(1.0)]

cluster_dict["10.0359997749_arrows"] += cgo_arrow([104.5,21.5,6.5], [101.899,22.062,7.793], color="blue red", name="Arrows_10.0359997749_2")

cluster_dict["10.0359997749"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(97.1313767055), float(17.5321400102), float(2.77730308865), float(1.0)]


cluster_dict["10.0359997749"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(103.351202763), float(21.3347492155), float(3.24077265433), float(1.0)]


cmd.load_cgo(cluster_dict["10.0359997749"], "Features_10.0359997749", 1)
cmd.load_cgo(cluster_dict["10.0359997749_arrows"], "Arrows_10.0359997749")
cmd.set("transparency", 0.2,"Features_10.0359997749")
cmd.group("Pharmacophore_10.0359997749", members="Features_10.0359997749")
cmd.group("Pharmacophore_10.0359997749", members="Arrows_10.0359997749")

if dirpath:
    f = join(dirpath, "10/label_threshold_10.0359997749.mol2")
else:
    f = "10/label_threshold_10.0359997749.mol2"

cmd.load(f, 'label_threshold_10.0359997749')
cmd.hide('everything', 'label_threshold_10.0359997749')
cmd.label("label_threshold_10.0359997749", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.0359997749', members= 'label_threshold_10.0359997749')


if dirpath:
    f = join(dirpath, '10/mesh.grd')
else:
    f = '10/mesh.grd'
cmd.load(f, 'mesh_10')
cmd.isomesh("isomesh_10", "mesh_10", 0.9)
cmd.color("grey80", "isomesh_10")
cmd.set('transparency', 0.4, "isomesh_10")

cmd.group('hotspot_10', "isomesh_10")
cmd.group('hotspot_10', "mesh_10")

if dirpath:
    f = join(dirpath, "11/label_threshold_7.0.mol2")
else:
    f = "11/label_threshold_7.0.mol2"

cmd.load(f, 'label_threshold_7.0')
cmd.hide('everything', 'label_threshold_7.0')
cmd.label("label_threshold_7.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [7.0]
gfiles = ['11/donor.grd', '11/apolar.grd', '11/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 11
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


cluster_dict = {"9.99899959564":[], "9.99899959564_arrows":[]}

cluster_dict["9.99899959564"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(93.0), float(19.0), float(6.0), float(1.0)]

cluster_dict["9.99899959564_arrows"] += cgo_arrow([93.0,19.0,6.0], [90.356,19.867,5.397], color="blue red", name="Arrows_9.99899959564_1")

cluster_dict["9.99899959564"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(104.5), float(21.5), float(6.5), float(1.0)]

cluster_dict["9.99899959564_arrows"] += cgo_arrow([104.5,21.5,6.5], [101.899,22.062,7.793], color="blue red", name="Arrows_9.99899959564_2")

cluster_dict["9.99899959564"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(97.5274370042), float(17.622131134), float(2.71739058639), float(1.0)]


cluster_dict["9.99899959564"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(103.566925354), float(21.4614145082), float(3.29646284724), float(1.0)]


cmd.load_cgo(cluster_dict["9.99899959564"], "Features_9.99899959564", 1)
cmd.load_cgo(cluster_dict["9.99899959564_arrows"], "Arrows_9.99899959564")
cmd.set("transparency", 0.2,"Features_9.99899959564")
cmd.group("Pharmacophore_9.99899959564", members="Features_9.99899959564")
cmd.group("Pharmacophore_9.99899959564", members="Arrows_9.99899959564")

if dirpath:
    f = join(dirpath, "11/label_threshold_9.99899959564.mol2")
else:
    f = "11/label_threshold_9.99899959564.mol2"

cmd.load(f, 'label_threshold_9.99899959564')
cmd.hide('everything', 'label_threshold_9.99899959564')
cmd.label("label_threshold_9.99899959564", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.99899959564', members= 'label_threshold_9.99899959564')


if dirpath:
    f = join(dirpath, '11/mesh.grd')
else:
    f = '11/mesh.grd'
cmd.load(f, 'mesh_11')
cmd.isomesh("isomesh_11", "mesh_11", 0.9)
cmd.color("grey80", "isomesh_11")
cmd.set('transparency', 0.4, "isomesh_11")

cmd.group('hotspot_11', "isomesh_11")
cmd.group('hotspot_11', "mesh_11")

if dirpath:
    f = join(dirpath, "12/label_threshold_2.6.mol2")
else:
    f = "12/label_threshold_2.6.mol2"

cmd.load(f, 'label_threshold_2.6')
cmd.hide('everything', 'label_threshold_2.6')
cmd.label("label_threshold_2.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [2.6]
gfiles = ['12/donor.grd', '12/apolar.grd', '12/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 12
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


cluster_dict = {"9.79050016403":[], "9.79050016403_arrows":[]}

cluster_dict["9.79050016403"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(87.304872726), float(1.07407413646), float(49.502619574), float(1.0)]


cmd.load_cgo(cluster_dict["9.79050016403"], "Features_9.79050016403", 1)
cmd.load_cgo(cluster_dict["9.79050016403_arrows"], "Arrows_9.79050016403")
cmd.set("transparency", 0.2,"Features_9.79050016403")
cmd.group("Pharmacophore_9.79050016403", members="Features_9.79050016403")
cmd.group("Pharmacophore_9.79050016403", members="Arrows_9.79050016403")

if dirpath:
    f = join(dirpath, "12/label_threshold_9.79050016403.mol2")
else:
    f = "12/label_threshold_9.79050016403.mol2"

cmd.load(f, 'label_threshold_9.79050016403')
cmd.hide('everything', 'label_threshold_9.79050016403')
cmd.label("label_threshold_9.79050016403", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.79050016403', members= 'label_threshold_9.79050016403')


if dirpath:
    f = join(dirpath, '12/mesh.grd')
else:
    f = '12/mesh.grd'
cmd.load(f, 'mesh_12')
cmd.isomesh("isomesh_12", "mesh_12", 0.9)
cmd.color("grey80", "isomesh_12")
cmd.set('transparency', 0.4, "isomesh_12")

cmd.group('hotspot_12', "isomesh_12")
cmd.group('hotspot_12', "mesh_12")

if dirpath:
    f = join(dirpath, "13/label_threshold_2.1.mol2")
else:
    f = "13/label_threshold_2.1.mol2"

cmd.load(f, 'label_threshold_2.1')
cmd.hide('everything', 'label_threshold_2.1')
cmd.label("label_threshold_2.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [2.1]
gfiles = ['13/donor.grd', '13/apolar.grd', '13/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 13
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


cluster_dict = {"9.26799964905":[], "9.26799964905_arrows":[]}

cluster_dict["9.26799964905"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(64.3118878347), float(15.8847154122), float(37.1257201725), float(1.0)]


cluster_dict["9.26799964905"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(66.5), float(18.0), float(36.0), float(1.0)]

cluster_dict["9.26799964905_arrows"] += cgo_arrow([66.5,18.0,36.0], [68.833,18.488,34.538], color="red blue", name="Arrows_9.26799964905_1")

cmd.load_cgo(cluster_dict["9.26799964905"], "Features_9.26799964905", 1)
cmd.load_cgo(cluster_dict["9.26799964905_arrows"], "Arrows_9.26799964905")
cmd.set("transparency", 0.2,"Features_9.26799964905")
cmd.group("Pharmacophore_9.26799964905", members="Features_9.26799964905")
cmd.group("Pharmacophore_9.26799964905", members="Arrows_9.26799964905")

if dirpath:
    f = join(dirpath, "13/label_threshold_9.26799964905.mol2")
else:
    f = "13/label_threshold_9.26799964905.mol2"

cmd.load(f, 'label_threshold_9.26799964905')
cmd.hide('everything', 'label_threshold_9.26799964905')
cmd.label("label_threshold_9.26799964905", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.26799964905', members= 'label_threshold_9.26799964905')


if dirpath:
    f = join(dirpath, '13/mesh.grd')
else:
    f = '13/mesh.grd'
cmd.load(f, 'mesh_13')
cmd.isomesh("isomesh_13", "mesh_13", 0.9)
cmd.color("grey80", "isomesh_13")
cmd.set('transparency', 0.4, "isomesh_13")

cmd.group('hotspot_13', "isomesh_13")
cmd.group('hotspot_13', "mesh_13")

if dirpath:
    f = join(dirpath, "14/label_threshold_0.6.mol2")
else:
    f = "14/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
gfiles = ['14/donor.grd', '14/apolar.grd', '14/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 14
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


cluster_dict = {"9.19999980927":[], "9.19999980927_arrows":[]}

cluster_dict["9.19999980927"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(68.1118304586), float(28.4602299377), float(-2.22879047404), float(1.0)]


cluster_dict["9.19999980927"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(73.1197596589), float(28.1664014504), float(-1.17924165912), float(1.0)]


cmd.load_cgo(cluster_dict["9.19999980927"], "Features_9.19999980927", 1)
cmd.load_cgo(cluster_dict["9.19999980927_arrows"], "Arrows_9.19999980927")
cmd.set("transparency", 0.2,"Features_9.19999980927")
cmd.group("Pharmacophore_9.19999980927", members="Features_9.19999980927")
cmd.group("Pharmacophore_9.19999980927", members="Arrows_9.19999980927")

if dirpath:
    f = join(dirpath, "14/label_threshold_9.19999980927.mol2")
else:
    f = "14/label_threshold_9.19999980927.mol2"

cmd.load(f, 'label_threshold_9.19999980927')
cmd.hide('everything', 'label_threshold_9.19999980927')
cmd.label("label_threshold_9.19999980927", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.19999980927', members= 'label_threshold_9.19999980927')


if dirpath:
    f = join(dirpath, '14/mesh.grd')
else:
    f = '14/mesh.grd'
cmd.load(f, 'mesh_14')
cmd.isomesh("isomesh_14", "mesh_14", 0.9)
cmd.color("grey80", "isomesh_14")
cmd.set('transparency', 0.4, "isomesh_14")

cmd.group('hotspot_14', "isomesh_14")
cmd.group('hotspot_14', "mesh_14")

if dirpath:
    f = join(dirpath, "15/label_threshold_30.0.mol2")
else:
    f = "15/label_threshold_30.0.mol2"

cmd.load(f, 'label_threshold_30.0')
cmd.hide('everything', 'label_threshold_30.0')
cmd.label("label_threshold_30.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [30.0]
gfiles = ['15/donor.grd', '15/apolar.grd', '15/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 15
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
    f = join(dirpath, "15/label_threshold_0.mol2")
else:
    f = "15/label_threshold_0.mol2"

cmd.load(f, 'label_threshold_0')
cmd.hide('everything', 'label_threshold_0')
cmd.label("label_threshold_0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_0', members= 'label_threshold_0')


if dirpath:
    f = join(dirpath, '15/mesh.grd')
else:
    f = '15/mesh.grd'
cmd.load(f, 'mesh_15')
cmd.isomesh("isomesh_15", "mesh_15", 0.9)
cmd.color("grey80", "isomesh_15")
cmd.set('transparency', 0.4, "isomesh_15")

cmd.group('hotspot_15', "isomesh_15")
cmd.group('hotspot_15', "mesh_15")

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
