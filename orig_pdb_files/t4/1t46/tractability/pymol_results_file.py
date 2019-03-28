
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
    f = join(dirpath, "0/label_threshold_1.1.mol2")
else:
    f = "0/label_threshold_1.1.mol2"

cmd.load(f, 'label_threshold_1.1')
cmd.hide('everything', 'label_threshold_1.1')
cmd.label("label_threshold_1.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [1.1]
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


cluster_dict = {"12.1120004654":[], "12.1120004654_arrows":[]}

cluster_dict["12.1120004654"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.4663398442), float(19.4567198701), float(20.2580733324), float(1.0)]


cluster_dict["12.1120004654"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(20.0), float(21.5), float(1.0)]

cluster_dict["12.1120004654_arrows"] += cgo_arrow([29.5,20.0,21.5], [32.144,18.793,22.139], color="red blue", name="Arrows_12.1120004654_1")

cmd.load_cgo(cluster_dict["12.1120004654"], "Features_12.1120004654", 1)
cmd.load_cgo(cluster_dict["12.1120004654_arrows"], "Arrows_12.1120004654")
cmd.set("transparency", 0.2,"Features_12.1120004654")
cmd.group("Pharmacophore_12.1120004654", members="Features_12.1120004654")
cmd.group("Pharmacophore_12.1120004654", members="Arrows_12.1120004654")

if dirpath:
    f = join(dirpath, "0/label_threshold_12.1120004654.mol2")
else:
    f = "0/label_threshold_12.1120004654.mol2"

cmd.load(f, 'label_threshold_12.1120004654')
cmd.hide('everything', 'label_threshold_12.1120004654')
cmd.label("label_threshold_12.1120004654", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.1120004654', members= 'label_threshold_12.1120004654')


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
    f = join(dirpath, "1/label_threshold_15.1.mol2")
else:
    f = "1/label_threshold_15.1.mol2"

cmd.load(f, 'label_threshold_15.1')
cmd.hide('everything', 'label_threshold_15.1')
cmd.label("label_threshold_15.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.1]
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


cluster_dict = {"16.8024997711":[], "16.8024997711_arrows":[]}

cluster_dict["16.8024997711"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(25.5), float(46.0), float(1.0)]

cluster_dict["16.8024997711_arrows"] += cgo_arrow([24.5,25.5,46.0], [24.992,22.671,45.425], color="blue red", name="Arrows_16.8024997711_1")

cluster_dict["16.8024997711"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.5), float(23.5), float(45.0), float(1.0)]

cluster_dict["16.8024997711_arrows"] += cgo_arrow([27.5,23.5,45.0], [24.992,22.671,45.425], color="blue red", name="Arrows_16.8024997711_2")

cluster_dict["16.8024997711"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(33.0), float(23.5), float(49.0), float(1.0)]

cluster_dict["16.8024997711_arrows"] += cgo_arrow([33.0,23.5,49.0], [33.792,20.626,47.987], color="blue red", name="Arrows_16.8024997711_3")

cluster_dict["16.8024997711"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(26.5971243077), float(26.0312835076), float(40.6655188582), float(1.0)]


cluster_dict["16.8024997711"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(33.7), float(24.2), float(49.5), float(1.0)]


cluster_dict["16.8024997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(26.0), float(31.5), float(1.0)]

cluster_dict["16.8024997711_arrows"] += cgo_arrow([24.0,26.0,31.5], [26.885,27.06,28.588], color="red blue", name="Arrows_16.8024997711_4")

cluster_dict["16.8024997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(24.5), float(44.0), float(1.0)]

cluster_dict["16.8024997711_arrows"] += cgo_arrow([26.0,24.5,44.0], [24.992,22.671,45.425], color="red blue", name="Arrows_16.8024997711_5")

cluster_dict["16.8024997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(26.5), float(41.0), float(1.0)]

cluster_dict["16.8024997711_arrows"] += cgo_arrow([27.0,26.5,41.0], [28.183,26.981,38.218], color="red blue", name="Arrows_16.8024997711_6")

cluster_dict["16.8024997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(24.5), float(44.0), float(1.0)]

cluster_dict["16.8024997711_arrows"] += cgo_arrow([26.0,24.5,44.0], [24.992,22.671,45.425], color="red blue", name="Arrows_16.8024997711_7")

cluster_dict["16.8024997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(28.0), float(43.0), float(1.0)]

cluster_dict["16.8024997711_arrows"] += cgo_arrow([28.0,28.0,43.0], [24.637,29.939,41.447], color="red blue", name="Arrows_16.8024997711_8")

cluster_dict["16.8024997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(31.0), float(23.0), float(47.5), float(1.0)]

cluster_dict["16.8024997711_arrows"] += cgo_arrow([31.0,23.0,47.5], [30.864,20.045,47.642], color="red blue", name="Arrows_16.8024997711_9")

cmd.load_cgo(cluster_dict["16.8024997711"], "Features_16.8024997711", 1)
cmd.load_cgo(cluster_dict["16.8024997711_arrows"], "Arrows_16.8024997711")
cmd.set("transparency", 0.2,"Features_16.8024997711")
cmd.group("Pharmacophore_16.8024997711", members="Features_16.8024997711")
cmd.group("Pharmacophore_16.8024997711", members="Arrows_16.8024997711")

if dirpath:
    f = join(dirpath, "1/label_threshold_16.8024997711.mol2")
else:
    f = "1/label_threshold_16.8024997711.mol2"

cmd.load(f, 'label_threshold_16.8024997711')
cmd.hide('everything', 'label_threshold_16.8024997711')
cmd.label("label_threshold_16.8024997711", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.8024997711', members= 'label_threshold_16.8024997711')


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
    f = join(dirpath, "2/label_threshold_13.4.mol2")
else:
    f = "2/label_threshold_13.4.mol2"

cmd.load(f, 'label_threshold_13.4')
cmd.hide('everything', 'label_threshold_13.4')
cmd.label("label_threshold_13.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.4]
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


cluster_dict = {"16.7150001526":[], "16.7150001526_arrows":[]}

cluster_dict["16.7150001526"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(26.0), float(40.0), float(1.0)]

cluster_dict["16.7150001526_arrows"] += cgo_arrow([24.5,26.0,40.0], [23.381,28.755,39.456], color="blue red", name="Arrows_16.7150001526_1")

cluster_dict["16.7150001526"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(25.5), float(46.0), float(1.0)]

cluster_dict["16.7150001526_arrows"] += cgo_arrow([24.5,25.5,46.0], [24.992,22.671,45.425], color="blue red", name="Arrows_16.7150001526_2")

cluster_dict["16.7150001526"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.5), float(23.5), float(45.0), float(1.0)]

cluster_dict["16.7150001526_arrows"] += cgo_arrow([27.5,23.5,45.0], [24.992,22.671,45.425], color="blue red", name="Arrows_16.7150001526_3")

cluster_dict["16.7150001526"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(33.0), float(23.5), float(49.0), float(1.0)]

cluster_dict["16.7150001526_arrows"] += cgo_arrow([33.0,23.5,49.0], [33.792,20.626,47.987], color="blue red", name="Arrows_16.7150001526_4")

cluster_dict["16.7150001526"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.7599303944), float(25.7926988318), float(42.8509451738), float(1.0)]


cluster_dict["16.7150001526"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.8757634693), float(30.1882811603), float(36.8736727743), float(1.0)]


cluster_dict["16.7150001526"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.7545298107), float(22.5), float(45.0), float(1.0)]


cluster_dict["16.7150001526"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(24.5), float(44.0), float(1.0)]

cluster_dict["16.7150001526_arrows"] += cgo_arrow([26.0,24.5,44.0], [24.992,22.671,45.425], color="red blue", name="Arrows_16.7150001526_5")

cluster_dict["16.7150001526"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(26.5), float(41.0), float(1.0)]

cluster_dict["16.7150001526_arrows"] += cgo_arrow([27.0,26.5,41.0], [28.183,26.981,38.218], color="red blue", name="Arrows_16.7150001526_6")

cluster_dict["16.7150001526"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(24.5), float(44.0), float(1.0)]

cluster_dict["16.7150001526_arrows"] += cgo_arrow([26.0,24.5,44.0], [24.992,22.671,45.425], color="red blue", name="Arrows_16.7150001526_7")

cluster_dict["16.7150001526"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(28.0), float(43.0), float(1.0)]

cluster_dict["16.7150001526_arrows"] += cgo_arrow([28.0,28.0,43.0], [24.637,29.939,41.447], color="red blue", name="Arrows_16.7150001526_8")

cluster_dict["16.7150001526"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(31.0), float(23.0), float(47.5), float(1.0)]

cluster_dict["16.7150001526_arrows"] += cgo_arrow([31.0,23.0,47.5], [30.864,20.045,47.642], color="red blue", name="Arrows_16.7150001526_9")

cmd.load_cgo(cluster_dict["16.7150001526"], "Features_16.7150001526", 1)
cmd.load_cgo(cluster_dict["16.7150001526_arrows"], "Arrows_16.7150001526")
cmd.set("transparency", 0.2,"Features_16.7150001526")
cmd.group("Pharmacophore_16.7150001526", members="Features_16.7150001526")
cmd.group("Pharmacophore_16.7150001526", members="Arrows_16.7150001526")

if dirpath:
    f = join(dirpath, "2/label_threshold_16.7150001526.mol2")
else:
    f = "2/label_threshold_16.7150001526.mol2"

cmd.load(f, 'label_threshold_16.7150001526')
cmd.hide('everything', 'label_threshold_16.7150001526')
cmd.label("label_threshold_16.7150001526", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.7150001526', members= 'label_threshold_16.7150001526')


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
    f = join(dirpath, "3/label_threshold_14.3.mol2")
else:
    f = "3/label_threshold_14.3.mol2"

cmd.load(f, 'label_threshold_14.3')
cmd.hide('everything', 'label_threshold_14.3')
cmd.label("label_threshold_14.3", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [14.3]
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


cluster_dict = {"16.6259994507":[], "16.6259994507_arrows":[]}

cluster_dict["16.6259994507"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.0), float(28.5), float(29.5), float(1.0)]

cluster_dict["16.6259994507_arrows"] += cgo_arrow([19.0,28.5,29.5], [16.337,28.69,28.158], color="blue red", name="Arrows_16.6259994507_1")

cluster_dict["16.6259994507"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(26.0), float(40.0), float(1.0)]

cluster_dict["16.6259994507_arrows"] += cgo_arrow([24.5,26.0,40.0], [23.381,28.755,39.456], color="blue red", name="Arrows_16.6259994507_2")

cluster_dict["16.6259994507"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(25.5), float(46.0), float(1.0)]

cluster_dict["16.6259994507_arrows"] += cgo_arrow([24.5,25.5,46.0], [24.992,22.671,45.425], color="blue red", name="Arrows_16.6259994507_3")

cluster_dict["16.6259994507"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.5), float(23.5), float(45.0), float(1.0)]

cluster_dict["16.6259994507_arrows"] += cgo_arrow([27.5,23.5,45.0], [24.992,22.671,45.425], color="blue red", name="Arrows_16.6259994507_4")

cluster_dict["16.6259994507"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.9009463069), float(26.2122918379), float(38.502418019), float(1.0)]


cluster_dict["16.6259994507"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(28.1671832195), float(21.8346363017), float(22.5836216946), float(1.0)]


cluster_dict["16.6259994507"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(26.0), float(31.5), float(1.0)]

cluster_dict["16.6259994507_arrows"] += cgo_arrow([24.0,26.0,31.5], [26.885,27.06,28.588], color="red blue", name="Arrows_16.6259994507_5")

cluster_dict["16.6259994507"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(24.5), float(44.0), float(1.0)]

cluster_dict["16.6259994507_arrows"] += cgo_arrow([26.0,24.5,44.0], [24.992,22.671,45.425], color="red blue", name="Arrows_16.6259994507_6")

cluster_dict["16.6259994507"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(26.5), float(41.0), float(1.0)]

cluster_dict["16.6259994507_arrows"] += cgo_arrow([27.0,26.5,41.0], [28.183,26.981,38.218], color="red blue", name="Arrows_16.6259994507_7")

cluster_dict["16.6259994507"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(24.5), float(44.0), float(1.0)]

cluster_dict["16.6259994507_arrows"] += cgo_arrow([26.0,24.5,44.0], [24.992,22.671,45.425], color="red blue", name="Arrows_16.6259994507_8")

cluster_dict["16.6259994507"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(28.0), float(43.0), float(1.0)]

cluster_dict["16.6259994507_arrows"] += cgo_arrow([28.0,28.0,43.0], [24.637,29.939,41.447], color="red blue", name="Arrows_16.6259994507_9")

cmd.load_cgo(cluster_dict["16.6259994507"], "Features_16.6259994507", 1)
cmd.load_cgo(cluster_dict["16.6259994507_arrows"], "Arrows_16.6259994507")
cmd.set("transparency", 0.2,"Features_16.6259994507")
cmd.group("Pharmacophore_16.6259994507", members="Features_16.6259994507")
cmd.group("Pharmacophore_16.6259994507", members="Arrows_16.6259994507")

if dirpath:
    f = join(dirpath, "3/label_threshold_16.6259994507.mol2")
else:
    f = "3/label_threshold_16.6259994507.mol2"

cmd.load(f, 'label_threshold_16.6259994507')
cmd.hide('everything', 'label_threshold_16.6259994507')
cmd.label("label_threshold_16.6259994507", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.6259994507', members= 'label_threshold_16.6259994507')


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
    f = join(dirpath, "4/label_threshold_1.1.mol2")
else:
    f = "4/label_threshold_1.1.mol2"

cmd.load(f, 'label_threshold_1.1')
cmd.hide('everything', 'label_threshold_1.1')
cmd.label("label_threshold_1.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [1.1]
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


cluster_dict = {"12.1120004654":[], "12.1120004654_arrows":[]}

cluster_dict["12.1120004654"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.4663398442), float(19.4567198701), float(20.2580733324), float(1.0)]


cluster_dict["12.1120004654"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(20.0), float(21.5), float(1.0)]

cluster_dict["12.1120004654_arrows"] += cgo_arrow([29.5,20.0,21.5], [32.144,18.793,22.139], color="red blue", name="Arrows_12.1120004654_1")

cmd.load_cgo(cluster_dict["12.1120004654"], "Features_12.1120004654", 1)
cmd.load_cgo(cluster_dict["12.1120004654_arrows"], "Arrows_12.1120004654")
cmd.set("transparency", 0.2,"Features_12.1120004654")
cmd.group("Pharmacophore_12.1120004654", members="Features_12.1120004654")
cmd.group("Pharmacophore_12.1120004654", members="Arrows_12.1120004654")

if dirpath:
    f = join(dirpath, "4/label_threshold_12.1120004654.mol2")
else:
    f = "4/label_threshold_12.1120004654.mol2"

cmd.load(f, 'label_threshold_12.1120004654')
cmd.hide('everything', 'label_threshold_12.1120004654')
cmd.label("label_threshold_12.1120004654", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.1120004654', members= 'label_threshold_12.1120004654')


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
    f = join(dirpath, "5/label_threshold_9.4.mol2")
else:
    f = "5/label_threshold_9.4.mol2"

cmd.load(f, 'label_threshold_9.4')
cmd.hide('everything', 'label_threshold_9.4')
cmd.label("label_threshold_9.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [9.4]
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


cluster_dict = {"11.6730003357":[], "11.6730003357_arrows":[]}

cluster_dict["11.6730003357"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(29.0), float(36.5), float(1.0)]

cluster_dict["11.6730003357_arrows"] += cgo_arrow([23.0,29.0,36.5], [23.381,28.755,39.456], color="blue red", name="Arrows_11.6730003357_1")

cluster_dict["11.6730003357"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(39.0), float(40.0), float(1.0)]

cluster_dict["11.6730003357_arrows"] += cgo_arrow([28.0,39.0,40.0], [30.007,37.831,41.593], color="blue red", name="Arrows_11.6730003357_2")

cluster_dict["11.6730003357"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.3719743516), float(31.6991726106), float(34.9598382072), float(1.0)]


cluster_dict["11.6730003357"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(26.1381133615), float(39.7422862788), float(38.208551505), float(1.0)]


cluster_dict["11.6730003357"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.7985709447), float(28.6888260001), float(42.8091877512), float(1.0)]


cluster_dict["11.6730003357"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.5), float(28.0), float(42.5), float(1.0)]

cluster_dict["11.6730003357_arrows"] += cgo_arrow([27.5,28.0,42.5], [24.637,29.939,41.447], color="red blue", name="Arrows_11.6730003357_3")

cmd.load_cgo(cluster_dict["11.6730003357"], "Features_11.6730003357", 1)
cmd.load_cgo(cluster_dict["11.6730003357_arrows"], "Arrows_11.6730003357")
cmd.set("transparency", 0.2,"Features_11.6730003357")
cmd.group("Pharmacophore_11.6730003357", members="Features_11.6730003357")
cmd.group("Pharmacophore_11.6730003357", members="Arrows_11.6730003357")

if dirpath:
    f = join(dirpath, "5/label_threshold_11.6730003357.mol2")
else:
    f = "5/label_threshold_11.6730003357.mol2"

cmd.load(f, 'label_threshold_11.6730003357')
cmd.hide('everything', 'label_threshold_11.6730003357')
cmd.label("label_threshold_11.6730003357", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.6730003357', members= 'label_threshold_11.6730003357')


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

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
