
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
    f = join(dirpath, "0/label_threshold_0.6.mol2")
else:
    f = "0/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
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


cluster_dict = {"13.298500061":[], "13.298500061_arrows":[]}

cluster_dict["13.298500061"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(39.0), float(54.0), float(0.0), float(1.0)]

cluster_dict["13.298500061_arrows"] += cgo_arrow([39.0,54.0,0.0], [38.382,52.515,1.817], color="blue red", name="Arrows_13.298500061_1")

cluster_dict["13.298500061"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(43.5), float(50.5), float(1.0), float(1.0)]

cluster_dict["13.298500061_arrows"] += cgo_arrow([43.5,50.5,1.0], [40.801,50.178,1.888], color="blue red", name="Arrows_13.298500061_2")

cluster_dict["13.298500061"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(42.0194417079), float(54.4453782578), float(-0.326891747616), float(1.0)]


cluster_dict["13.298500061"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(43.1828469715), float(54.1258821272), float(3.04447199637), float(1.0)]


cluster_dict["13.298500061"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(42.5), float(53.0), float(2.0), float(1.0)]

cluster_dict["13.298500061_arrows"] += cgo_arrow([42.5,53.0,2.0], [44.208,54.939,3.202], color="red blue", name="Arrows_13.298500061_3")

cluster_dict["13.298500061"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(42.5), float(55.0), float(-2.5), float(1.0)]

cluster_dict["13.298500061_arrows"] += cgo_arrow([42.5,55.0,-2.5], [38.006,54.184,-3.021], color="red blue", name="Arrows_13.298500061_4")

cmd.load_cgo(cluster_dict["13.298500061"], "Features_13.298500061", 1)
cmd.load_cgo(cluster_dict["13.298500061_arrows"], "Arrows_13.298500061")
cmd.set("transparency", 0.2,"Features_13.298500061")
cmd.group("Pharmacophore_13.298500061", members="Features_13.298500061")
cmd.group("Pharmacophore_13.298500061", members="Arrows_13.298500061")

if dirpath:
    f = join(dirpath, "0/label_threshold_13.298500061.mol2")
else:
    f = "0/label_threshold_13.298500061.mol2"

cmd.load(f, 'label_threshold_13.298500061')
cmd.hide('everything', 'label_threshold_13.298500061')
cmd.label("label_threshold_13.298500061", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.298500061', members= 'label_threshold_13.298500061')


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
    f = join(dirpath, "1/label_threshold_3.2.mol2")
else:
    f = "1/label_threshold_3.2.mol2"

cmd.load(f, 'label_threshold_3.2')
cmd.hide('everything', 'label_threshold_3.2')
cmd.label("label_threshold_3.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [3.2]
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


cluster_dict = {"11.1130003929":[], "11.1130003929_arrows":[]}

cluster_dict["11.1130003929"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(43.7310117538), float(62.3978209539), float(21.5866284149), float(1.0)]


cmd.load_cgo(cluster_dict["11.1130003929"], "Features_11.1130003929", 1)
cmd.load_cgo(cluster_dict["11.1130003929_arrows"], "Arrows_11.1130003929")
cmd.set("transparency", 0.2,"Features_11.1130003929")
cmd.group("Pharmacophore_11.1130003929", members="Features_11.1130003929")
cmd.group("Pharmacophore_11.1130003929", members="Arrows_11.1130003929")

if dirpath:
    f = join(dirpath, "1/label_threshold_11.1130003929.mol2")
else:
    f = "1/label_threshold_11.1130003929.mol2"

cmd.load(f, 'label_threshold_11.1130003929')
cmd.hide('everything', 'label_threshold_11.1130003929')
cmd.label("label_threshold_11.1130003929", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.1130003929', members= 'label_threshold_11.1130003929')


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
    f = join(dirpath, "2/label_threshold_2.0.mol2")
else:
    f = "2/label_threshold_2.0.mol2"

cmd.load(f, 'label_threshold_2.0')
cmd.hide('everything', 'label_threshold_2.0')
cmd.label("label_threshold_2.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [2.0]
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


cluster_dict = {"9.60999965668":[], "9.60999965668_arrows":[]}

cluster_dict["9.60999965668"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(34.5), float(63.0), float(7.0), float(1.0)]

cluster_dict["9.60999965668_arrows"] += cgo_arrow([34.5,63.0,7.0], [37.541,62.959,7.026], color="blue red", name="Arrows_9.60999965668_1")

cluster_dict["9.60999965668"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(36.0), float(62.5), float(2.0), float(1.0)]

cluster_dict["9.60999965668_arrows"] += cgo_arrow([36.0,62.5,2.0], [36.053,60.45,-0.501], color="blue red", name="Arrows_9.60999965668_2")

cluster_dict["9.60999965668"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(32.6925991402), float(62.4146726398), float(4.7620542119), float(1.0)]


cluster_dict["9.60999965668"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(38.8520681456), float(58.0), float(0.0901466618677), float(1.0)]


cluster_dict["9.60999965668"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(65.0), float(3.5), float(1.0)]

cluster_dict["9.60999965668_arrows"] += cgo_arrow([36.5,65.0,3.5], [38.006,67.0,1.792], color="red blue", name="Arrows_9.60999965668_3")

cluster_dict["9.60999965668"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(35.0), float(61.5), float(2.0), float(1.0)]

cluster_dict["9.60999965668_arrows"] += cgo_arrow([35.0,61.5,2.0], [33.884,59.974,-0.824], color="red blue", name="Arrows_9.60999965668_4")

cmd.load_cgo(cluster_dict["9.60999965668"], "Features_9.60999965668", 1)
cmd.load_cgo(cluster_dict["9.60999965668_arrows"], "Arrows_9.60999965668")
cmd.set("transparency", 0.2,"Features_9.60999965668")
cmd.group("Pharmacophore_9.60999965668", members="Features_9.60999965668")
cmd.group("Pharmacophore_9.60999965668", members="Arrows_9.60999965668")

if dirpath:
    f = join(dirpath, "2/label_threshold_9.60999965668.mol2")
else:
    f = "2/label_threshold_9.60999965668.mol2"

cmd.load(f, 'label_threshold_9.60999965668')
cmd.hide('everything', 'label_threshold_9.60999965668')
cmd.label("label_threshold_9.60999965668", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.60999965668', members= 'label_threshold_9.60999965668')


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
    f = join(dirpath, "3/label_threshold_30.0.mol2")
else:
    f = "3/label_threshold_30.0.mol2"

cmd.load(f, 'label_threshold_30.0')
cmd.hide('everything', 'label_threshold_30.0')
cmd.label("label_threshold_30.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [30.0]
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


cluster_dict = {"0":[], "0_arrows":[]}

cmd.load_cgo(cluster_dict["0"], "Features_0", 1)
cmd.load_cgo(cluster_dict["0_arrows"], "Arrows_0")
cmd.set("transparency", 0.2,"Features_0")
cmd.group("Pharmacophore_0", members="Features_0")
cmd.group("Pharmacophore_0", members="Arrows_0")

if dirpath:
    f = join(dirpath, "3/label_threshold_0.mol2")
else:
    f = "3/label_threshold_0.mol2"

cmd.load(f, 'label_threshold_0')
cmd.hide('everything', 'label_threshold_0')
cmd.label("label_threshold_0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_0', members= 'label_threshold_0')


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

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
