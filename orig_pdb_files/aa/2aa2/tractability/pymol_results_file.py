
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
    f = join(dirpath, "0/label_threshold_5.0.mol2")
else:
    f = "0/label_threshold_5.0.mol2"

cmd.load(f, 'label_threshold_5.0')
cmd.hide('everything', 'label_threshold_5.0')
cmd.label("label_threshold_5.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [5.0]
gfiles = ['0/apolar.grd', '0/acceptor.grd']
grids = ['apolar', 'acceptor']
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


cluster_dict = {"10.7729997635":[], "10.7729997635_arrows":[]}

cluster_dict["10.7729997635"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(19.4019033114), float(56.3459384754), float(29.5553483161), float(1.0)]


cluster_dict["10.7729997635"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.0), float(60.5), float(32.0), float(1.0)]

cluster_dict["10.7729997635_arrows"] += cgo_arrow([20.0,60.5,32.0], [22.301,61.764,33.421], color="red blue", name="Arrows_10.7729997635_1")

cluster_dict["10.7729997635"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(57.0), float(30.0), float(1.0)]

cluster_dict["10.7729997635_arrows"] += cgo_arrow([21.0,57.0,30.0], [24.957,55.819,31.022], color="red blue", name="Arrows_10.7729997635_2")

cmd.load_cgo(cluster_dict["10.7729997635"], "Features_10.7729997635", 1)
cmd.load_cgo(cluster_dict["10.7729997635_arrows"], "Arrows_10.7729997635")
cmd.set("transparency", 0.2,"Features_10.7729997635")
cmd.group("Pharmacophore_10.7729997635", members="Features_10.7729997635")
cmd.group("Pharmacophore_10.7729997635", members="Arrows_10.7729997635")

if dirpath:
    f = join(dirpath, "0/label_threshold_10.7729997635.mol2")
else:
    f = "0/label_threshold_10.7729997635.mol2"

cmd.load(f, 'label_threshold_10.7729997635')
cmd.hide('everything', 'label_threshold_10.7729997635')
cmd.label("label_threshold_10.7729997635", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.7729997635', members= 'label_threshold_10.7729997635')


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
    f = join(dirpath, "1/label_threshold_11.3.mol2")
else:
    f = "1/label_threshold_11.3.mol2"

cmd.load(f, 'label_threshold_11.3')
cmd.hide('everything', 'label_threshold_11.3')
cmd.label("label_threshold_11.3", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [11.3]
gfiles = ['1/apolar.grd', '1/acceptor.grd']
grids = ['apolar', 'acceptor']
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


cluster_dict = {"17.0569992065":[], "17.0569992065_arrows":[]}

cluster_dict["17.0569992065"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(75.0), float(24.5), float(1.0)]

cluster_dict["17.0569992065_arrows"] += cgo_arrow([18.0,75.0,24.5], [15.611,73.163,23.589], color="blue red", name="Arrows_17.0569992065_1")

cluster_dict["17.0569992065"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.0), float(68.0), float(17.0), float(1.0)]

cluster_dict["17.0569992065_arrows"] += cgo_arrow([19.0,68.0,17.0], [17.176,65.615,17.937], color="blue red", name="Arrows_17.0569992065_2")

cluster_dict["17.0569992065"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(19.1835102977), float(72.372201813), float(19.9057546997), float(1.0)]


cluster_dict["17.0569992065"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(77.0), float(23.0), float(1.0)]

cluster_dict["17.0569992065_arrows"] += cgo_arrow([17.5,77.0,23.0], [18.025,77.713,26.505], color="red blue", name="Arrows_17.0569992065_3")

cluster_dict["17.0569992065"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(68.0), float(15.5), float(1.0)]

cluster_dict["17.0569992065_arrows"] += cgo_arrow([17.5,68.0,15.5], [16.119,65.678,15.964], color="red blue", name="Arrows_17.0569992065_4")

cluster_dict["17.0569992065"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.0), float(70.0), float(22.5), float(1.0)]

cluster_dict["17.0569992065_arrows"] += cgo_arrow([17.0,70.0,22.5], [14.032,67.68,20.514], color="red blue", name="Arrows_17.0569992065_5")

cluster_dict["17.0569992065"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.0), float(72.5), float(17.5), float(1.0)]


cluster_dict["17.0569992065"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(71.0), float(24.0), float(1.0)]

cluster_dict["17.0569992065_arrows"] += cgo_arrow([17.5,71.0,24.0], [20.683,72.995,26.534], color="red blue", name="Arrows_17.0569992065_6")

cluster_dict["17.0569992065"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.0), float(72.0), float(17.0), float(1.0)]


cluster_dict["17.0569992065"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.5), float(73.5), float(18.5), float(1.0)]


cmd.load_cgo(cluster_dict["17.0569992065"], "Features_17.0569992065", 1)
cmd.load_cgo(cluster_dict["17.0569992065_arrows"], "Arrows_17.0569992065")
cmd.set("transparency", 0.2,"Features_17.0569992065")
cmd.group("Pharmacophore_17.0569992065", members="Features_17.0569992065")
cmd.group("Pharmacophore_17.0569992065", members="Arrows_17.0569992065")

if dirpath:
    f = join(dirpath, "1/label_threshold_17.0569992065.mol2")
else:
    f = "1/label_threshold_17.0569992065.mol2"

cmd.load(f, 'label_threshold_17.0569992065')
cmd.hide('everything', 'label_threshold_17.0569992065')
cmd.label("label_threshold_17.0569992065", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.0569992065', members= 'label_threshold_17.0569992065')


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
    f = join(dirpath, "2/label_threshold_11.8.mol2")
else:
    f = "2/label_threshold_11.8.mol2"

cmd.load(f, 'label_threshold_11.8')
cmd.hide('everything', 'label_threshold_11.8')
cmd.label("label_threshold_11.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [11.8]
gfiles = ['2/apolar.grd', '2/acceptor.grd']
grids = ['apolar', 'acceptor']
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


cluster_dict = {"16.4340000153":[], "16.4340000153_arrows":[]}

cluster_dict["16.4340000153"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.0), float(68.0), float(17.0), float(1.0)]

cluster_dict["16.4340000153_arrows"] += cgo_arrow([19.0,68.0,17.0], [17.176,65.615,17.937], color="blue red", name="Arrows_16.4340000153_1")

cluster_dict["16.4340000153"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(70.5), float(22.0), float(1.0)]

cluster_dict["16.4340000153_arrows"] += cgo_arrow([17.5,70.5,22.0], [15.611,73.163,23.589], color="blue red", name="Arrows_16.4340000153_2")

cluster_dict["16.4340000153"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(60.0), float(13.5), float(1.0)]

cluster_dict["16.4340000153_arrows"] += cgo_arrow([21.0,60.0,13.5], [21.096,57.707,12.4], color="blue red", name="Arrows_16.4340000153_3")

cluster_dict["16.4340000153"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(18.5679992943), float(60.7823053207), float(14.817135333), float(1.0)]


cluster_dict["16.4340000153"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(18.8705344336), float(70.0146616287), float(18.5104191233), float(1.0)]


cluster_dict["16.4340000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(61.5), float(14.0), float(1.0)]

cluster_dict["16.4340000153_arrows"] += cgo_arrow([16.0,61.5,14.0], [12.745,61.224,13.4], color="red blue", name="Arrows_16.4340000153_4")

cluster_dict["16.4340000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.5), float(57.5), float(8.5), float(1.0)]

cluster_dict["16.4340000153_arrows"] += cgo_arrow([15.5,57.5,8.5], [13.898,58.668,11.267], color="red blue", name="Arrows_16.4340000153_5")

cluster_dict["16.4340000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(68.0), float(15.5), float(1.0)]

cluster_dict["16.4340000153_arrows"] += cgo_arrow([17.5,68.0,15.5], [16.119,65.678,15.964], color="red blue", name="Arrows_16.4340000153_6")

cluster_dict["16.4340000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.0), float(70.0), float(22.5), float(1.0)]

cluster_dict["16.4340000153_arrows"] += cgo_arrow([17.0,70.0,22.5], [14.032,67.68,20.514], color="red blue", name="Arrows_16.4340000153_7")

cluster_dict["16.4340000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(19.0), float(63.5), float(15.0), float(1.0)]

cluster_dict["16.4340000153_arrows"] += cgo_arrow([19.0,63.5,15.0], [18.607,64.009,12.924], color="red blue", name="Arrows_16.4340000153_8")

cluster_dict["16.4340000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(60.5), float(15.0), float(1.0)]

cluster_dict["16.4340000153_arrows"] += cgo_arrow([21.0,60.5,15.0], [22.773,57.721,15.021], color="red blue", name="Arrows_16.4340000153_9")

cmd.load_cgo(cluster_dict["16.4340000153"], "Features_16.4340000153", 1)
cmd.load_cgo(cluster_dict["16.4340000153_arrows"], "Arrows_16.4340000153")
cmd.set("transparency", 0.2,"Features_16.4340000153")
cmd.group("Pharmacophore_16.4340000153", members="Features_16.4340000153")
cmd.group("Pharmacophore_16.4340000153", members="Arrows_16.4340000153")

if dirpath:
    f = join(dirpath, "2/label_threshold_16.4340000153.mol2")
else:
    f = "2/label_threshold_16.4340000153.mol2"

cmd.load(f, 'label_threshold_16.4340000153')
cmd.hide('everything', 'label_threshold_16.4340000153')
cmd.label("label_threshold_16.4340000153", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.4340000153', members= 'label_threshold_16.4340000153')


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
    f = join(dirpath, "3/label_threshold_2.1.mol2")
else:
    f = "3/label_threshold_2.1.mol2"

cmd.load(f, 'label_threshold_2.1')
cmd.hide('everything', 'label_threshold_2.1')
cmd.label("label_threshold_2.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [2.1]
gfiles = ['3/apolar.grd', '3/acceptor.grd']
grids = ['apolar', 'acceptor']
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


cluster_dict = {"13.9250001907":[], "13.9250001907_arrows":[]}

cluster_dict["13.9250001907"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(67.5), float(25.5), float(1.0)]

cluster_dict["13.9250001907_arrows"] += cgo_arrow([30.5,67.5,25.5], [32.461,65.224,25.43], color="blue red", name="Arrows_13.9250001907_1")

cluster_dict["13.9250001907"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(33.0), float(68.0), float(24.5), float(1.0)]

cluster_dict["13.9250001907_arrows"] += cgo_arrow([33.0,68.0,24.5], [32.461,65.224,25.43], color="blue red", name="Arrows_13.9250001907_2")

cluster_dict["13.9250001907"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(32.4419175365), float(70.5567056117), float(26.8137605897), float(1.0)]


cluster_dict["13.9250001907"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(29.0), float(67.5), float(28.0), float(1.0)]

cluster_dict["13.9250001907_arrows"] += cgo_arrow([29.0,67.5,28.0], [26.136,68.268,28.292], color="red blue", name="Arrows_13.9250001907_3")

cluster_dict["13.9250001907"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(66.5), float(26.0), float(1.0)]

cluster_dict["13.9250001907_arrows"] += cgo_arrow([30.5,66.5,26.0], [32.461,65.224,25.43], color="red blue", name="Arrows_13.9250001907_4")

cluster_dict["13.9250001907"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(31.0), float(68.0), float(24.5), float(1.0)]

cluster_dict["13.9250001907_arrows"] += cgo_arrow([31.0,68.0,24.5], [32.461,65.224,25.43], color="red blue", name="Arrows_13.9250001907_5")

cluster_dict["13.9250001907"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(31.0), float(68.0), float(24.5), float(1.0)]

cluster_dict["13.9250001907_arrows"] += cgo_arrow([31.0,68.0,24.5], [32.461,65.224,25.43], color="red blue", name="Arrows_13.9250001907_6")

cmd.load_cgo(cluster_dict["13.9250001907"], "Features_13.9250001907", 1)
cmd.load_cgo(cluster_dict["13.9250001907_arrows"], "Arrows_13.9250001907")
cmd.set("transparency", 0.2,"Features_13.9250001907")
cmd.group("Pharmacophore_13.9250001907", members="Features_13.9250001907")
cmd.group("Pharmacophore_13.9250001907", members="Arrows_13.9250001907")

if dirpath:
    f = join(dirpath, "3/label_threshold_13.9250001907.mol2")
else:
    f = "3/label_threshold_13.9250001907.mol2"

cmd.load(f, 'label_threshold_13.9250001907')
cmd.hide('everything', 'label_threshold_13.9250001907')
cmd.label("label_threshold_13.9250001907", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.9250001907', members= 'label_threshold_13.9250001907')


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
    f = join(dirpath, "4/label_threshold_4.7.mol2")
else:
    f = "4/label_threshold_4.7.mol2"

cmd.load(f, 'label_threshold_4.7')
cmd.hide('everything', 'label_threshold_4.7')
cmd.label("label_threshold_4.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [4.7]
gfiles = ['4/apolar.grd', '4/acceptor.grd']
grids = ['apolar', 'acceptor']
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


cluster_dict = {"10.7769999504":[], "10.7769999504_arrows":[]}

cluster_dict["10.7769999504"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(19.3910358758), float(56.21023648), float(29.5411388004), float(1.0)]


cluster_dict["10.7769999504"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.0), float(60.5), float(32.0), float(1.0)]

cluster_dict["10.7769999504_arrows"] += cgo_arrow([20.0,60.5,32.0], [22.301,61.764,33.421], color="red blue", name="Arrows_10.7769999504_1")

cluster_dict["10.7769999504"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(57.0), float(30.0), float(1.0)]

cluster_dict["10.7769999504_arrows"] += cgo_arrow([21.0,57.0,30.0], [24.957,55.819,31.022], color="red blue", name="Arrows_10.7769999504_2")

cmd.load_cgo(cluster_dict["10.7769999504"], "Features_10.7769999504", 1)
cmd.load_cgo(cluster_dict["10.7769999504_arrows"], "Arrows_10.7769999504")
cmd.set("transparency", 0.2,"Features_10.7769999504")
cmd.group("Pharmacophore_10.7769999504", members="Features_10.7769999504")
cmd.group("Pharmacophore_10.7769999504", members="Arrows_10.7769999504")

if dirpath:
    f = join(dirpath, "4/label_threshold_10.7769999504.mol2")
else:
    f = "4/label_threshold_10.7769999504.mol2"

cmd.load(f, 'label_threshold_10.7769999504')
cmd.hide('everything', 'label_threshold_10.7769999504')
cmd.label("label_threshold_10.7769999504", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.7769999504', members= 'label_threshold_10.7769999504')


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
    f = join(dirpath, "5/label_threshold_5.0.mol2")
else:
    f = "5/label_threshold_5.0.mol2"

cmd.load(f, 'label_threshold_5.0')
cmd.hide('everything', 'label_threshold_5.0')
cmd.label("label_threshold_5.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [5.0]
gfiles = ['5/apolar.grd', '5/acceptor.grd']
grids = ['apolar', 'acceptor']
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


cluster_dict = {"10.7729997635":[], "10.7729997635_arrows":[]}

cluster_dict["10.7729997635"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(19.4007452381), float(56.3528322981), float(29.5558446954), float(1.0)]


cluster_dict["10.7729997635"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.0), float(60.5), float(32.0), float(1.0)]

cluster_dict["10.7729997635_arrows"] += cgo_arrow([20.0,60.5,32.0], [22.301,61.764,33.421], color="red blue", name="Arrows_10.7729997635_1")

cluster_dict["10.7729997635"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(57.0), float(30.0), float(1.0)]

cluster_dict["10.7729997635_arrows"] += cgo_arrow([21.0,57.0,30.0], [24.957,55.819,31.022], color="red blue", name="Arrows_10.7729997635_2")

cmd.load_cgo(cluster_dict["10.7729997635"], "Features_10.7729997635", 1)
cmd.load_cgo(cluster_dict["10.7729997635_arrows"], "Arrows_10.7729997635")
cmd.set("transparency", 0.2,"Features_10.7729997635")
cmd.group("Pharmacophore_10.7729997635", members="Features_10.7729997635")
cmd.group("Pharmacophore_10.7729997635", members="Arrows_10.7729997635")

if dirpath:
    f = join(dirpath, "5/label_threshold_10.7729997635.mol2")
else:
    f = "5/label_threshold_10.7729997635.mol2"

cmd.load(f, 'label_threshold_10.7729997635')
cmd.hide('everything', 'label_threshold_10.7729997635')
cmd.label("label_threshold_10.7729997635", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.7729997635', members= 'label_threshold_10.7729997635')


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
gfiles = ['6/apolar.grd', '6/acceptor.grd']
grids = ['apolar', 'acceptor']
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


cluster_dict = {"7.89400005341":[], "7.89400005341_arrows":[]}

cluster_dict["7.89400005341"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(33.6174661052), float(46.2438200616), float(23.5251805639), float(1.0)]


cmd.load_cgo(cluster_dict["7.89400005341"], "Features_7.89400005341", 1)
cmd.load_cgo(cluster_dict["7.89400005341_arrows"], "Arrows_7.89400005341")
cmd.set("transparency", 0.2,"Features_7.89400005341")
cmd.group("Pharmacophore_7.89400005341", members="Features_7.89400005341")
cmd.group("Pharmacophore_7.89400005341", members="Arrows_7.89400005341")

if dirpath:
    f = join(dirpath, "6/label_threshold_7.89400005341.mol2")
else:
    f = "6/label_threshold_7.89400005341.mol2"

cmd.load(f, 'label_threshold_7.89400005341')
cmd.hide('everything', 'label_threshold_7.89400005341')
cmd.label("label_threshold_7.89400005341", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_7.89400005341', members= 'label_threshold_7.89400005341')


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
gfiles = ['7/apolar.grd', '7/acceptor.grd']
grids = ['apolar', 'acceptor']
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
