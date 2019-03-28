
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
    f = join(dirpath, "0/label_threshold_12.4.mol2")
else:
    f = "0/label_threshold_12.4.mol2"

cmd.load(f, 'label_threshold_12.4')
cmd.hide('everything', 'label_threshold_12.4')
cmd.label("label_threshold_12.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.4]
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


cluster_dict = {"14.3319997787":[], "14.3319997787_arrows":[]}

cluster_dict["14.3319997787"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(11.0), float(-14.5), float(1.0)]

cluster_dict["14.3319997787_arrows"] += cgo_arrow([26.0,11.0,-14.5], [22.615,11.131,-15.816], color="blue red", name="Arrows_14.3319997787_1")

cluster_dict["14.3319997787"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(15.5), float(-14.5), float(1.0)]

cluster_dict["14.3319997787_arrows"] += cgo_arrow([24.5,15.5,-14.5], [21.438,14.463,-13.601], color="blue red", name="Arrows_14.3319997787_2")

cluster_dict["14.3319997787"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(10.0), float(-13.0), float(1.0)]

cluster_dict["14.3319997787_arrows"] += cgo_arrow([29.5,10.0,-13.0], [30.129,7.657,-14.523], color="blue red", name="Arrows_14.3319997787_3")

cluster_dict["14.3319997787"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(3.5), float(-13.5), float(1.0)]

cluster_dict["14.3319997787_arrows"] += cgo_arrow([27.0,3.5,-13.5], [25.564,1.866,-16.048], color="blue red", name="Arrows_14.3319997787_4")

cluster_dict["14.3319997787"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(10.0), float(-13.0), float(1.0)]

cluster_dict["14.3319997787_arrows"] += cgo_arrow([29.5,10.0,-13.0], [30.129,7.657,-14.523], color="blue red", name="Arrows_14.3319997787_5")

cluster_dict["14.3319997787"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(3.5), float(-13.0), float(1.0)]

cluster_dict["14.3319997787_arrows"] += cgo_arrow([30.5,3.5,-13.0], [30.724,3.981,-15.898], color="blue red", name="Arrows_14.3319997787_6")

cluster_dict["14.3319997787"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(11.9528144616), float(15.6829743685), float(-11.1581897091), float(1.0)]


cluster_dict["14.3319997787"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(14.5556313242), float(3.64270082404), float(-21.9514643567), float(1.0)]


cluster_dict["14.3319997787"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(17.2845646153), float(11.7139041008), float(-11.3087409079), float(1.0)]


cluster_dict["14.3319997787"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(26.3945451802), float(10.1001345614), float(-14.1595758016), float(1.0)]


cluster_dict["14.3319997787"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.0915969252), float(13.0815150041), float(-18.8400166931), float(1.0)]


cluster_dict["14.3319997787"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(28.181050209), float(12.3642123432), float(-17.6852740594), float(1.0)]


cluster_dict["14.3319997787"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(29.8360577735), float(1.93559283288), float(-12.0136131098), float(1.0)]


cluster_dict["14.3319997787"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(9.5), float(-21.0), float(1.0)]

cluster_dict["14.3319997787_arrows"] += cgo_arrow([26.0,9.5,-21.0], [23.813,11.613,-19.373], color="red blue", name="Arrows_14.3319997787_7")

cluster_dict["14.3319997787"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(29.0), float(10.5), float(-12.0), float(1.0)]

cluster_dict["14.3319997787_arrows"] += cgo_arrow([29.0,10.5,-12.0], [29.581,9.727,-9.209], color="red blue", name="Arrows_14.3319997787_8")

cmd.load_cgo(cluster_dict["14.3319997787"], "Features_14.3319997787", 1)
cmd.load_cgo(cluster_dict["14.3319997787_arrows"], "Arrows_14.3319997787")
cmd.set("transparency", 0.2,"Features_14.3319997787")
cmd.group("Pharmacophore_14.3319997787", members="Features_14.3319997787")
cmd.group("Pharmacophore_14.3319997787", members="Arrows_14.3319997787")

if dirpath:
    f = join(dirpath, "0/label_threshold_14.3319997787.mol2")
else:
    f = "0/label_threshold_14.3319997787.mol2"

cmd.load(f, 'label_threshold_14.3319997787')
cmd.hide('everything', 'label_threshold_14.3319997787')
cmd.label("label_threshold_14.3319997787", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.3319997787', members= 'label_threshold_14.3319997787')


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
    f = join(dirpath, "1/label_threshold_16.2.mol2")
else:
    f = "1/label_threshold_16.2.mol2"

cmd.load(f, 'label_threshold_16.2')
cmd.hide('everything', 'label_threshold_16.2')
cmd.label("label_threshold_16.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [16.2]
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


cluster_dict = {"17.3020000458":[], "17.3020000458_arrows":[]}

cluster_dict["17.3020000458"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(1.5), float(-34.0), float(1.0)]

cluster_dict["17.3020000458_arrows"] += cgo_arrow([24.5,1.5,-34.0], [26.171,2.494,-31.486], color="blue red", name="Arrows_17.3020000458_1")

cluster_dict["17.3020000458"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.5), float(10.0), float(-34.0), float(1.0)]

cluster_dict["17.3020000458_arrows"] += cgo_arrow([27.5,10.0,-34.0], [26.257,12.933,-32.801], color="blue red", name="Arrows_17.3020000458_2")

cluster_dict["17.3020000458"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(10.0), float(-37.5), float(1.0)]

cluster_dict["17.3020000458_arrows"] += cgo_arrow([27.0,10.0,-37.5], [25.878,8.199,-38.518], color="blue red", name="Arrows_17.3020000458_3")

cluster_dict["17.3020000458"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(10.5), float(-32.5), float(1.0)]

cluster_dict["17.3020000458_arrows"] += cgo_arrow([28.0,10.5,-32.5], [26.257,12.933,-32.801], color="blue red", name="Arrows_17.3020000458_4")

cluster_dict["17.3020000458"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(31.0), float(10.0), float(-30.5), float(1.0)]

cluster_dict["17.3020000458_arrows"] += cgo_arrow([31.0,10.0,-30.5], [27.36,9.822,-28.982], color="blue red", name="Arrows_17.3020000458_5")

cluster_dict["17.3020000458"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(32.0), float(3.5), float(-37.0), float(1.0)]

cluster_dict["17.3020000458_arrows"] += cgo_arrow([32.0,3.5,-37.0], [32.104,1.695,-34.881], color="blue red", name="Arrows_17.3020000458_6")

cluster_dict["17.3020000458"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(33.5), float(3.5), float(-42.0), float(1.0)]

cluster_dict["17.3020000458_arrows"] += cgo_arrow([33.5,3.5,-42.0], [30.388,4.44,-42.566], color="blue red", name="Arrows_17.3020000458_7")

cluster_dict["17.3020000458"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(30.2497098117), float(-3.02246221681), float(-37.5495152191), float(1.0)]


cluster_dict["17.3020000458"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(32.810687637), float(9.80946495703), float(-36.1106162261), float(1.0)]


cluster_dict["17.3020000458"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(7.5), float(-34.5), float(1.0)]

cluster_dict["17.3020000458_arrows"] += cgo_arrow([30.5,7.5,-34.5], [28.503,5.78,-35.679], color="red blue", name="Arrows_17.3020000458_8")

cluster_dict["17.3020000458"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(33.5), float(3.0), float(-44.0), float(1.0)]

cluster_dict["17.3020000458_arrows"] += cgo_arrow([33.5,3.0,-44.0], [30.965,4.679,-44.712], color="red blue", name="Arrows_17.3020000458_9")

cmd.load_cgo(cluster_dict["17.3020000458"], "Features_17.3020000458", 1)
cmd.load_cgo(cluster_dict["17.3020000458_arrows"], "Arrows_17.3020000458")
cmd.set("transparency", 0.2,"Features_17.3020000458")
cmd.group("Pharmacophore_17.3020000458", members="Features_17.3020000458")
cmd.group("Pharmacophore_17.3020000458", members="Arrows_17.3020000458")

if dirpath:
    f = join(dirpath, "1/label_threshold_17.3020000458.mol2")
else:
    f = "1/label_threshold_17.3020000458.mol2"

cmd.load(f, 'label_threshold_17.3020000458')
cmd.hide('everything', 'label_threshold_17.3020000458')
cmd.label("label_threshold_17.3020000458", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.3020000458', members= 'label_threshold_17.3020000458')


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
    f = join(dirpath, "2/label_threshold_16.0.mol2")
else:
    f = "2/label_threshold_16.0.mol2"

cmd.load(f, 'label_threshold_16.0')
cmd.hide('everything', 'label_threshold_16.0')
cmd.label("label_threshold_16.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [16.0]
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


cluster_dict = {"17.2910003662":[], "17.2910003662_arrows":[]}

cluster_dict["17.2910003662"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(1.5), float(-34.0), float(1.0)]

cluster_dict["17.2910003662_arrows"] += cgo_arrow([24.5,1.5,-34.0], [26.171,2.494,-31.486], color="blue red", name="Arrows_17.2910003662_1")

cluster_dict["17.2910003662"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.5), float(10.0), float(-34.0), float(1.0)]

cluster_dict["17.2910003662_arrows"] += cgo_arrow([27.5,10.0,-34.0], [26.257,12.933,-32.801], color="blue red", name="Arrows_17.2910003662_2")

cluster_dict["17.2910003662"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(10.5), float(-32.5), float(1.0)]

cluster_dict["17.2910003662_arrows"] += cgo_arrow([28.0,10.5,-32.5], [26.257,12.933,-32.801], color="blue red", name="Arrows_17.2910003662_3")

cluster_dict["17.2910003662"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(31.0), float(10.0), float(-30.5), float(1.0)]

cluster_dict["17.2910003662_arrows"] += cgo_arrow([31.0,10.0,-30.5], [27.36,9.822,-28.982], color="blue red", name="Arrows_17.2910003662_4")

cluster_dict["17.2910003662"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(32.0), float(3.5), float(-37.0), float(1.0)]

cluster_dict["17.2910003662_arrows"] += cgo_arrow([32.0,3.5,-37.0], [32.104,1.695,-34.881], color="blue red", name="Arrows_17.2910003662_5")

cluster_dict["17.2910003662"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(33.5), float(3.5), float(-42.0), float(1.0)]

cluster_dict["17.2910003662_arrows"] += cgo_arrow([33.5,3.5,-42.0], [30.388,4.44,-42.566], color="blue red", name="Arrows_17.2910003662_6")

cluster_dict["17.2910003662"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(32.9143862897), float(10.0372706189), float(-36.1100660116), float(1.0)]


cluster_dict["17.2910003662"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(35.0), float(14.5), float(-25.75), float(1.0)]


cluster_dict["17.2910003662"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.5), float(16.0), float(-31.0), float(1.0)]

cluster_dict["17.2910003662_arrows"] += cgo_arrow([20.5,16.0,-31.0], [20.697,19.878,-30.934], color="red blue", name="Arrows_17.2910003662_7")

cluster_dict["17.2910003662"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(7.5), float(-34.5), float(1.0)]

cluster_dict["17.2910003662_arrows"] += cgo_arrow([30.5,7.5,-34.5], [28.503,5.78,-35.679], color="red blue", name="Arrows_17.2910003662_8")

cluster_dict["17.2910003662"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(33.0), float(13.0), float(-28.5), float(1.0)]

cluster_dict["17.2910003662_arrows"] += cgo_arrow([33.0,13.0,-28.5], [32.581,16.518,-27.932], color="red blue", name="Arrows_17.2910003662_9")

cmd.load_cgo(cluster_dict["17.2910003662"], "Features_17.2910003662", 1)
cmd.load_cgo(cluster_dict["17.2910003662_arrows"], "Arrows_17.2910003662")
cmd.set("transparency", 0.2,"Features_17.2910003662")
cmd.group("Pharmacophore_17.2910003662", members="Features_17.2910003662")
cmd.group("Pharmacophore_17.2910003662", members="Arrows_17.2910003662")

if dirpath:
    f = join(dirpath, "2/label_threshold_17.2910003662.mol2")
else:
    f = "2/label_threshold_17.2910003662.mol2"

cmd.load(f, 'label_threshold_17.2910003662')
cmd.hide('everything', 'label_threshold_17.2910003662')
cmd.label("label_threshold_17.2910003662", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.2910003662', members= 'label_threshold_17.2910003662')


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
    f = join(dirpath, "3/label_threshold_15.7.mol2")
else:
    f = "3/label_threshold_15.7.mol2"

cmd.load(f, 'label_threshold_15.7')
cmd.hide('everything', 'label_threshold_15.7')
cmd.label("label_threshold_15.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.7]
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


cluster_dict = {"17.1779994965":[], "17.1779994965_arrows":[]}

cluster_dict["17.1779994965"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.5), float(10.0), float(-34.0), float(1.0)]

cluster_dict["17.1779994965_arrows"] += cgo_arrow([27.5,10.0,-34.0], [26.257,12.933,-32.801], color="blue red", name="Arrows_17.1779994965_1")

cluster_dict["17.1779994965"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(10.0), float(-37.5), float(1.0)]

cluster_dict["17.1779994965_arrows"] += cgo_arrow([27.0,10.0,-37.5], [25.878,8.199,-38.518], color="blue red", name="Arrows_17.1779994965_2")

cluster_dict["17.1779994965"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(31.0), float(10.0), float(-30.5), float(1.0)]

cluster_dict["17.1779994965_arrows"] += cgo_arrow([31.0,10.0,-30.5], [27.36,9.822,-28.982], color="blue red", name="Arrows_17.1779994965_3")

cluster_dict["17.1779994965"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(32.0), float(3.5), float(-37.0), float(1.0)]

cluster_dict["17.1779994965_arrows"] += cgo_arrow([32.0,3.5,-37.0], [32.104,1.695,-34.881], color="blue red", name="Arrows_17.1779994965_4")

cluster_dict["17.1779994965"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(34.5), float(6.5), float(-41.5), float(1.0)]

cluster_dict["17.1779994965_arrows"] += cgo_arrow([34.5,6.5,-41.5], [33.786,9.103,-42.297], color="blue red", name="Arrows_17.1779994965_5")

cluster_dict["17.1779994965"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(33.5), float(15.0), float(-31.5), float(1.0)]

cluster_dict["17.1779994965_arrows"] += cgo_arrow([33.5,15.0,-31.5], [35.064,16.748,-31.789], color="blue red", name="Arrows_17.1779994965_6")

cluster_dict["17.1779994965"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(34.5), float(11.0), float(-30.5), float(1.0)]

cluster_dict["17.1779994965_arrows"] += cgo_arrow([34.5,11.0,-30.5], [36.832,9.96,-29.646], color="blue red", name="Arrows_17.1779994965_7")

cluster_dict["17.1779994965"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(32.879029608), float(10.5553529184), float(-35.7768506863), float(1.0)]


cluster_dict["17.1779994965"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(35.0), float(14.5), float(-25.75), float(1.0)]


cluster_dict["17.1779994965"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(7.5), float(-34.5), float(1.0)]

cluster_dict["17.1779994965_arrows"] += cgo_arrow([30.5,7.5,-34.5], [28.503,5.78,-35.679], color="red blue", name="Arrows_17.1779994965_8")

cluster_dict["17.1779994965"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(33.0), float(13.0), float(-28.5), float(1.0)]

cluster_dict["17.1779994965_arrows"] += cgo_arrow([33.0,13.0,-28.5], [32.581,16.518,-27.932], color="red blue", name="Arrows_17.1779994965_9")

cmd.load_cgo(cluster_dict["17.1779994965"], "Features_17.1779994965", 1)
cmd.load_cgo(cluster_dict["17.1779994965_arrows"], "Arrows_17.1779994965")
cmd.set("transparency", 0.2,"Features_17.1779994965")
cmd.group("Pharmacophore_17.1779994965", members="Features_17.1779994965")
cmd.group("Pharmacophore_17.1779994965", members="Arrows_17.1779994965")

if dirpath:
    f = join(dirpath, "3/label_threshold_17.1779994965.mol2")
else:
    f = "3/label_threshold_17.1779994965.mol2"

cmd.load(f, 'label_threshold_17.1779994965')
cmd.hide('everything', 'label_threshold_17.1779994965')
cmd.label("label_threshold_17.1779994965", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.1779994965', members= 'label_threshold_17.1779994965')


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
    f = join(dirpath, "4/label_threshold_13.7.mol2")
else:
    f = "4/label_threshold_13.7.mol2"

cmd.load(f, 'label_threshold_13.7')
cmd.hide('everything', 'label_threshold_13.7')
cmd.label("label_threshold_13.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.7]
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


cluster_dict = {"16.1760005951":[], "16.1760005951_arrows":[]}

cluster_dict["16.1760005951"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(1.5), float(-34.0), float(1.0)]

cluster_dict["16.1760005951_arrows"] += cgo_arrow([24.5,1.5,-34.0], [26.171,2.494,-31.486], color="blue red", name="Arrows_16.1760005951_1")

cluster_dict["16.1760005951"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(-5.5), float(-36.5), float(1.0)]

cluster_dict["16.1760005951_arrows"] += cgo_arrow([29.5,-5.5,-36.5], [30.398,-6.958,-34.45], color="blue red", name="Arrows_16.1760005951_2")

cluster_dict["16.1760005951"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(32.0), float(3.5), float(-37.0), float(1.0)]

cluster_dict["16.1760005951_arrows"] += cgo_arrow([32.0,3.5,-37.0], [32.104,1.695,-34.881], color="blue red", name="Arrows_16.1760005951_3")

cluster_dict["16.1760005951"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(34.5), float(6.5), float(-41.5), float(1.0)]

cluster_dict["16.1760005951_arrows"] += cgo_arrow([34.5,6.5,-41.5], [33.786,9.103,-42.297], color="blue red", name="Arrows_16.1760005951_4")

cluster_dict["16.1760005951"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(33.5), float(3.5), float(-42.0), float(1.0)]

cluster_dict["16.1760005951_arrows"] += cgo_arrow([33.5,3.5,-42.0], [30.388,4.44,-42.566], color="blue red", name="Arrows_16.1760005951_5")

cluster_dict["16.1760005951"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(2.5), float(-43.5), float(1.0)]

cluster_dict["16.1760005951_arrows"] += cgo_arrow([36.5,2.5,-43.5], [38.722,3.854,-44.392], color="blue red", name="Arrows_16.1760005951_6")

cluster_dict["16.1760005951"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(37.0), float(5.5), float(-40.0), float(1.0)]

cluster_dict["16.1760005951_arrows"] += cgo_arrow([37.0,5.5,-40.0], [38.1,7.79,-41.551], color="blue red", name="Arrows_16.1760005951_7")

cluster_dict["16.1760005951"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(30.0867257132), float(-2.88805050218), float(-37.3847165677), float(1.0)]


cluster_dict["16.1760005951"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(32.9513575424), float(5.51057666049), float(-39.0850243626), float(1.0)]


cluster_dict["16.1760005951"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(7.5), float(-34.5), float(1.0)]

cluster_dict["16.1760005951_arrows"] += cgo_arrow([30.5,7.5,-34.5], [28.503,5.78,-35.679], color="red blue", name="Arrows_16.1760005951_8")

cluster_dict["16.1760005951"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(33.5), float(3.0), float(-44.0), float(1.0)]

cluster_dict["16.1760005951_arrows"] += cgo_arrow([33.5,3.0,-44.0], [30.965,4.679,-44.712], color="red blue", name="Arrows_16.1760005951_9")

cmd.load_cgo(cluster_dict["16.1760005951"], "Features_16.1760005951", 1)
cmd.load_cgo(cluster_dict["16.1760005951_arrows"], "Arrows_16.1760005951")
cmd.set("transparency", 0.2,"Features_16.1760005951")
cmd.group("Pharmacophore_16.1760005951", members="Features_16.1760005951")
cmd.group("Pharmacophore_16.1760005951", members="Arrows_16.1760005951")

if dirpath:
    f = join(dirpath, "4/label_threshold_16.1760005951.mol2")
else:
    f = "4/label_threshold_16.1760005951.mol2"

cmd.load(f, 'label_threshold_16.1760005951')
cmd.hide('everything', 'label_threshold_16.1760005951')
cmd.label("label_threshold_16.1760005951", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.1760005951', members= 'label_threshold_16.1760005951')


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
    f = join(dirpath, "5/label_threshold_13.4.mol2")
else:
    f = "5/label_threshold_13.4.mol2"

cmd.load(f, 'label_threshold_13.4')
cmd.hide('everything', 'label_threshold_13.4')
cmd.label("label_threshold_13.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.4]
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


cluster_dict = {"14.3319997787":[], "14.3319997787_arrows":[]}

cluster_dict["14.3319997787"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(9.0), float(-14.0), float(1.0)]

cluster_dict["14.3319997787_arrows"] += cgo_arrow([24.5,9.0,-14.0], [22.615,11.131,-15.816], color="blue red", name="Arrows_14.3319997787_1")

cluster_dict["14.3319997787"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(3.5), float(-13.5), float(1.0)]

cluster_dict["14.3319997787_arrows"] += cgo_arrow([27.0,3.5,-13.5], [25.564,1.866,-16.048], color="blue red", name="Arrows_14.3319997787_2")

cluster_dict["14.3319997787"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(3.5), float(-13.0), float(1.0)]

cluster_dict["14.3319997787_arrows"] += cgo_arrow([30.5,3.5,-13.0], [30.724,3.981,-15.898], color="blue red", name="Arrows_14.3319997787_3")

cluster_dict["14.3319997787"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(33.0), float(1.5), float(-11.5), float(1.0)]

cluster_dict["14.3319997787_arrows"] += cgo_arrow([33.0,1.5,-11.5], [32.955,-0.155,-14.973], color="blue red", name="Arrows_14.3319997787_4")

cluster_dict["14.3319997787"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(26.9016157003), float(10.353445548), float(-14.9214165016), float(1.0)]


cluster_dict["14.3319997787"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(33.1333824328), float(-0.567564630309), float(-9.71542154965), float(1.0)]


cluster_dict["14.3319997787"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(39.7566668899), float(6.47943923326), float(-8.82091283732), float(1.0)]


cluster_dict["14.3319997787"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(43.1), float(7.85), float(-15.9), float(1.0)]


cluster_dict["14.3319997787"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(28.5), float(10.5), float(-22.0), float(1.0)]

cluster_dict["14.3319997787_arrows"] += cgo_arrow([28.5,10.5,-22.0], [26.191,13.098,-22.189], color="red blue", name="Arrows_14.3319997787_5")

cluster_dict["14.3319997787"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(10.0), float(-13.5), float(1.0)]

cluster_dict["14.3319997787_arrows"] += cgo_arrow([29.5,10.0,-13.5], [32.505,10.662,-12.636], color="red blue", name="Arrows_14.3319997787_6")

cluster_dict["14.3319997787"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(10.0), float(-13.5), float(1.0)]

cluster_dict["14.3319997787_arrows"] += cgo_arrow([29.5,10.0,-13.5], [32.505,10.662,-12.636], color="red blue", name="Arrows_14.3319997787_7")

cluster_dict["14.3319997787"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(6.0), float(-10.5), float(1.0)]

cluster_dict["14.3319997787_arrows"] += cgo_arrow([36.5,6.0,-10.5], [33.693,5.068,-9.683], color="red blue", name="Arrows_14.3319997787_8")

cluster_dict["14.3319997787"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(43.0), float(6.5), float(-14.5), float(1.0)]

cluster_dict["14.3319997787_arrows"] += cgo_arrow([43.0,6.5,-14.5], [43.481,4.46,-13.067], color="red blue", name="Arrows_14.3319997787_9")

cmd.load_cgo(cluster_dict["14.3319997787"], "Features_14.3319997787", 1)
cmd.load_cgo(cluster_dict["14.3319997787_arrows"], "Arrows_14.3319997787")
cmd.set("transparency", 0.2,"Features_14.3319997787")
cmd.group("Pharmacophore_14.3319997787", members="Features_14.3319997787")
cmd.group("Pharmacophore_14.3319997787", members="Arrows_14.3319997787")

if dirpath:
    f = join(dirpath, "5/label_threshold_14.3319997787.mol2")
else:
    f = "5/label_threshold_14.3319997787.mol2"

cmd.load(f, 'label_threshold_14.3319997787')
cmd.hide('everything', 'label_threshold_14.3319997787')
cmd.label("label_threshold_14.3319997787", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.3319997787', members= 'label_threshold_14.3319997787')


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
    f = join(dirpath, "6/label_threshold_9.9.mol2")
else:
    f = "6/label_threshold_9.9.mol2"

cmd.load(f, 'label_threshold_9.9')
cmd.hide('everything', 'label_threshold_9.9')
cmd.label("label_threshold_9.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [9.9]
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


cluster_dict = {"14.2839999199":[], "14.2839999199_arrows":[]}

cluster_dict["14.2839999199"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.5), float(10.0), float(-34.0), float(1.0)]

cluster_dict["14.2839999199_arrows"] += cgo_arrow([27.5,10.0,-34.0], [26.257,12.933,-32.801], color="blue red", name="Arrows_14.2839999199_1")

cluster_dict["14.2839999199"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(10.0), float(-37.5), float(1.0)]

cluster_dict["14.2839999199_arrows"] += cgo_arrow([27.0,10.0,-37.5], [25.878,8.199,-38.518], color="blue red", name="Arrows_14.2839999199_2")

cluster_dict["14.2839999199"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(31.0), float(10.0), float(-30.5), float(1.0)]

cluster_dict["14.2839999199_arrows"] += cgo_arrow([31.0,10.0,-30.5], [27.36,9.822,-28.982], color="blue red", name="Arrows_14.2839999199_3")

cluster_dict["14.2839999199"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(32.5), float(15.0), float(-32.5), float(1.0)]

cluster_dict["14.2839999199_arrows"] += cgo_arrow([32.5,15.0,-32.5], [35.064,16.748,-31.789], color="blue red", name="Arrows_14.2839999199_4")

cluster_dict["14.2839999199"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(15.5602840031), float(4.01430399498), float(-25.1611819989), float(1.0)]


cluster_dict["14.2839999199"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(18.3739704751), float(2.8583845611), float(-25.8925464539), float(1.0)]


cluster_dict["14.2839999199"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(20.5252740104), float(13.6785852336), float(-30.5730124156), float(1.0)]


cluster_dict["14.2839999199"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.8039993419), float(2.39558575163), float(-34.8006883022), float(1.0)]


cluster_dict["14.2839999199"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(30.1899577109), float(11.1510877904), float(-33.7557407051), float(1.0)]


cluster_dict["14.2839999199"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.5), float(16.0), float(-31.0), float(1.0)]

cluster_dict["14.2839999199_arrows"] += cgo_arrow([20.5,16.0,-31.0], [20.697,19.878,-30.934], color="red blue", name="Arrows_14.2839999199_5")

cluster_dict["14.2839999199"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(15.5), float(-33.5), float(1.0)]

cluster_dict["14.2839999199_arrows"] += cgo_arrow([21.0,15.5,-33.5], [23.42,17.057,-33.409], color="red blue", name="Arrows_14.2839999199_6")

cluster_dict["14.2839999199"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(30.0), float(7.5), float(-34.5), float(1.0)]

cluster_dict["14.2839999199_arrows"] += cgo_arrow([30.0,7.5,-34.5], [28.503,5.78,-35.679], color="red blue", name="Arrows_14.2839999199_7")

cluster_dict["14.2839999199"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(30.0), float(15.0), float(-30.5), float(1.0)]

cluster_dict["14.2839999199_arrows"] += cgo_arrow([30.0,15.0,-30.5], [29.452,18.335,-31.502], color="red blue", name="Arrows_14.2839999199_8")

cluster_dict["14.2839999199"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(30.0), float(15.0), float(-30.5), float(1.0)]

cluster_dict["14.2839999199_arrows"] += cgo_arrow([30.0,15.0,-30.5], [29.452,18.335,-31.502], color="red blue", name="Arrows_14.2839999199_9")

cmd.load_cgo(cluster_dict["14.2839999199"], "Features_14.2839999199", 1)
cmd.load_cgo(cluster_dict["14.2839999199_arrows"], "Arrows_14.2839999199")
cmd.set("transparency", 0.2,"Features_14.2839999199")
cmd.group("Pharmacophore_14.2839999199", members="Features_14.2839999199")
cmd.group("Pharmacophore_14.2839999199", members="Arrows_14.2839999199")

if dirpath:
    f = join(dirpath, "6/label_threshold_14.2839999199.mol2")
else:
    f = "6/label_threshold_14.2839999199.mol2"

cmd.load(f, 'label_threshold_14.2839999199')
cmd.hide('everything', 'label_threshold_14.2839999199')
cmd.label("label_threshold_14.2839999199", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.2839999199', members= 'label_threshold_14.2839999199')


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
    f = join(dirpath, "7/label_threshold_12.7.mol2")
else:
    f = "7/label_threshold_12.7.mol2"

cmd.load(f, 'label_threshold_12.7')
cmd.hide('everything', 'label_threshold_12.7')
cmd.label("label_threshold_12.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.7]
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


cluster_dict = {"14.1510000229":[], "14.1510000229_arrows":[]}

cluster_dict["14.1510000229"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(3.5), float(-13.5), float(1.0)]

cluster_dict["14.1510000229_arrows"] += cgo_arrow([27.0,3.5,-13.5], [25.564,1.866,-16.048], color="blue red", name="Arrows_14.1510000229_1")

cluster_dict["14.1510000229"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(10.0), float(-13.0), float(1.0)]

cluster_dict["14.1510000229_arrows"] += cgo_arrow([29.5,10.0,-13.0], [30.129,7.657,-14.523], color="blue red", name="Arrows_14.1510000229_2")

cluster_dict["14.1510000229"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(3.5), float(-13.0), float(1.0)]

cluster_dict["14.1510000229_arrows"] += cgo_arrow([30.5,3.5,-13.0], [30.724,3.981,-15.898], color="blue red", name="Arrows_14.1510000229_3")

cluster_dict["14.1510000229"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(35.0), float(0.5), float(-9.0), float(1.0)]

cluster_dict["14.1510000229_arrows"] += cgo_arrow([35.0,0.5,-9.0], [37.234,0.497,-7.21], color="blue red", name="Arrows_14.1510000229_4")

cluster_dict["14.1510000229"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(33.0), float(1.5), float(-11.5), float(1.0)]

cluster_dict["14.1510000229_arrows"] += cgo_arrow([33.0,1.5,-11.5], [32.955,-0.155,-14.973], color="blue red", name="Arrows_14.1510000229_5")

cluster_dict["14.1510000229"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.4900218989), float(9.3137254336), float(-12.7802427545), float(1.0)]


cluster_dict["14.1510000229"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(32.9909044295), float(-0.557912329366), float(-9.82485130433), float(1.0)]


cluster_dict["14.1510000229"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(39.8654915103), float(6.5498225685), float(-8.78921525031), float(1.0)]


cluster_dict["14.1510000229"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(44.912008325), float(7.7916346073), float(-16.4200540805), float(1.0)]


cluster_dict["14.1510000229"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(10.0), float(-13.5), float(1.0)]

cluster_dict["14.1510000229_arrows"] += cgo_arrow([29.5,10.0,-13.5], [32.505,10.662,-12.636], color="red blue", name="Arrows_14.1510000229_6")

cluster_dict["14.1510000229"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(10.0), float(-13.5), float(1.0)]

cluster_dict["14.1510000229_arrows"] += cgo_arrow([29.5,10.0,-13.5], [32.505,10.662,-12.636], color="red blue", name="Arrows_14.1510000229_7")

cluster_dict["14.1510000229"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(6.0), float(-10.5), float(1.0)]

cluster_dict["14.1510000229_arrows"] += cgo_arrow([36.5,6.0,-10.5], [33.693,5.068,-9.683], color="red blue", name="Arrows_14.1510000229_8")

cluster_dict["14.1510000229"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(43.5), float(7.0), float(-14.5), float(1.0)]

cluster_dict["14.1510000229_arrows"] += cgo_arrow([43.5,7.0,-14.5], [43.481,4.46,-13.067], color="red blue", name="Arrows_14.1510000229_9")

cmd.load_cgo(cluster_dict["14.1510000229"], "Features_14.1510000229", 1)
cmd.load_cgo(cluster_dict["14.1510000229_arrows"], "Arrows_14.1510000229")
cmd.set("transparency", 0.2,"Features_14.1510000229")
cmd.group("Pharmacophore_14.1510000229", members="Features_14.1510000229")
cmd.group("Pharmacophore_14.1510000229", members="Arrows_14.1510000229")

if dirpath:
    f = join(dirpath, "7/label_threshold_14.1510000229.mol2")
else:
    f = "7/label_threshold_14.1510000229.mol2"

cmd.load(f, 'label_threshold_14.1510000229')
cmd.hide('everything', 'label_threshold_14.1510000229')
cmd.label("label_threshold_14.1510000229", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.1510000229', members= 'label_threshold_14.1510000229')


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
    f = join(dirpath, "8/label_threshold_12.7.mol2")
else:
    f = "8/label_threshold_12.7.mol2"

cmd.load(f, 'label_threshold_12.7')
cmd.hide('everything', 'label_threshold_12.7')
cmd.label("label_threshold_12.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.7]
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


cluster_dict = {"14.0":[], "14.0_arrows":[]}

cluster_dict["14.0"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(9.0), float(-14.0), float(1.0)]

cluster_dict["14.0_arrows"] += cgo_arrow([24.5,9.0,-14.0], [22.615,11.131,-15.816], color="blue red", name="Arrows_14.0_1")

cluster_dict["14.0"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(15.5), float(-14.5), float(1.0)]

cluster_dict["14.0_arrows"] += cgo_arrow([24.5,15.5,-14.5], [21.438,14.463,-13.601], color="blue red", name="Arrows_14.0_2")

cluster_dict["14.0"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(18.0), float(-16.5), float(1.0)]

cluster_dict["14.0_arrows"] += cgo_arrow([24.5,18.0,-16.5], [21.965,15.819,-16.753], color="blue red", name="Arrows_14.0_3")

cluster_dict["14.0"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(3.5), float(-13.5), float(1.0)]

cluster_dict["14.0_arrows"] += cgo_arrow([27.0,3.5,-13.5], [25.564,1.866,-16.048], color="blue red", name="Arrows_14.0_4")

cluster_dict["14.0"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.0), float(12.0), float(-23.0), float(1.0)]

cluster_dict["14.0_arrows"] += cgo_arrow([29.0,12.0,-23.0], [28.974,12.038,-25.36], color="blue red", name="Arrows_14.0_5")

cluster_dict["14.0"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(19.0), float(-15.5), float(1.0)]

cluster_dict["14.0_arrows"] += cgo_arrow([29.5,19.0,-15.5], [28.582,16.633,-14.757], color="blue red", name="Arrows_14.0_6")

cluster_dict["14.0"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(3.5), float(-13.0), float(1.0)]

cluster_dict["14.0_arrows"] += cgo_arrow([30.5,3.5,-13.0], [30.724,3.981,-15.898], color="blue red", name="Arrows_14.0_7")

cluster_dict["14.0"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(33.0), float(1.5), float(-11.5), float(1.0)]

cluster_dict["14.0_arrows"] += cgo_arrow([33.0,1.5,-11.5], [32.955,-0.155,-14.973], color="blue red", name="Arrows_14.0_8")

cluster_dict["14.0"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(17.3192351041), float(11.7254059638), float(-11.3102511024), float(1.0)]


cluster_dict["14.0"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.2798068642), float(10.8051640018), float(-16.0911946304), float(1.0)]


cluster_dict["14.0"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(32.2316277565), float(1.66794143404), float(-11.2761829627), float(1.0)]


cluster_dict["14.0"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(37.521438452), float(7.32046695871), float(-11.023390858), float(1.0)]


cluster_dict["14.0"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(29.0), float(10.5), float(-12.0), float(1.0)]

cluster_dict["14.0_arrows"] += cgo_arrow([29.0,10.5,-12.0], [29.581,9.727,-9.209], color="red blue", name="Arrows_14.0_9")

cluster_dict["14.0"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(10.0), float(-13.5), float(1.0)]

cluster_dict["14.0_arrows"] += cgo_arrow([29.5,10.0,-13.5], [32.505,10.662,-12.636], color="red blue", name="Arrows_14.0_10")

cmd.load_cgo(cluster_dict["14.0"], "Features_14.0", 1)
cmd.load_cgo(cluster_dict["14.0_arrows"], "Arrows_14.0")
cmd.set("transparency", 0.2,"Features_14.0")
cmd.group("Pharmacophore_14.0", members="Features_14.0")
cmd.group("Pharmacophore_14.0", members="Arrows_14.0")

if dirpath:
    f = join(dirpath, "8/label_threshold_14.0.mol2")
else:
    f = "8/label_threshold_14.0.mol2"

cmd.load(f, 'label_threshold_14.0')
cmd.hide('everything', 'label_threshold_14.0')
cmd.label("label_threshold_14.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.0', members= 'label_threshold_14.0')


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
    f = join(dirpath, "9/label_threshold_10.6.mol2")
else:
    f = "9/label_threshold_10.6.mol2"

cmd.load(f, 'label_threshold_10.6')
cmd.hide('everything', 'label_threshold_10.6')
cmd.label("label_threshold_10.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.6]
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


cluster_dict = {"13.579000473":[], "13.579000473_arrows":[]}

cluster_dict["13.579000473"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(3.5), float(-13.5), float(1.0)]

cluster_dict["13.579000473_arrows"] += cgo_arrow([27.0,3.5,-13.5], [25.564,1.866,-16.048], color="blue red", name="Arrows_13.579000473_1")

cluster_dict["13.579000473"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(28.5), float(7.5), float(-12.0), float(1.0)]

cluster_dict["13.579000473_arrows"] += cgo_arrow([28.5,7.5,-12.0], [29.116,7.566,-9.101], color="blue red", name="Arrows_13.579000473_2")

cluster_dict["13.579000473"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(3.5), float(-13.0), float(1.0)]

cluster_dict["13.579000473_arrows"] += cgo_arrow([30.5,3.5,-13.0], [30.724,3.981,-15.898], color="blue red", name="Arrows_13.579000473_3")

cluster_dict["13.579000473"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(33.0), float(1.5), float(-11.5), float(1.0)]

cluster_dict["13.579000473_arrows"] += cgo_arrow([33.0,1.5,-11.5], [32.955,-0.155,-14.973], color="blue red", name="Arrows_13.579000473_4")

cluster_dict["13.579000473"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(33.0), float(1.5), float(-11.5), float(1.0)]

cluster_dict["13.579000473_arrows"] += cgo_arrow([33.0,1.5,-11.5], [32.955,-0.155,-14.973], color="blue red", name="Arrows_13.579000473_5")

cluster_dict["13.579000473"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(37.0), float(8.5), float(-11.5), float(1.0)]

cluster_dict["13.579000473_arrows"] += cgo_arrow([37.0,8.5,-11.5], [34.648,10.111,-12.296], color="blue red", name="Arrows_13.579000473_6")

cluster_dict["13.579000473"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(37.0), float(5.5), float(-9.0), float(1.0)]

cluster_dict["13.579000473_arrows"] += cgo_arrow([37.0,5.5,-9.0], [35.325,7.016,-7.626], color="blue red", name="Arrows_13.579000473_7")

cluster_dict["13.579000473"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(26.4997064005), float(6.15373323095), float(-12.7304530467), float(1.0)]


cluster_dict["13.579000473"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(26.8), float(7.0), float(-11.4), float(1.0)]


cluster_dict["13.579000473"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(33.0939473219), float(-0.884017927955), float(-9.63851353317), float(1.0)]


cluster_dict["13.579000473"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(35.5), float(8.0), float(-12.5), float(1.0)]


cluster_dict["13.579000473"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(39.2236921289), float(5.56844039466), float(-9.12958121532), float(1.0)]


cluster_dict["13.579000473"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(6.0), float(-10.5), float(1.0)]

cluster_dict["13.579000473_arrows"] += cgo_arrow([36.5,6.0,-10.5], [33.693,5.068,-9.683], color="red blue", name="Arrows_13.579000473_8")

cmd.load_cgo(cluster_dict["13.579000473"], "Features_13.579000473", 1)
cmd.load_cgo(cluster_dict["13.579000473_arrows"], "Arrows_13.579000473")
cmd.set("transparency", 0.2,"Features_13.579000473")
cmd.group("Pharmacophore_13.579000473", members="Features_13.579000473")
cmd.group("Pharmacophore_13.579000473", members="Arrows_13.579000473")

if dirpath:
    f = join(dirpath, "9/label_threshold_13.579000473.mol2")
else:
    f = "9/label_threshold_13.579000473.mol2"

cmd.load(f, 'label_threshold_13.579000473')
cmd.hide('everything', 'label_threshold_13.579000473')
cmd.label("label_threshold_13.579000473", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.579000473', members= 'label_threshold_13.579000473')


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
    f = join(dirpath, "10/label_threshold_11.8.mol2")
else:
    f = "10/label_threshold_11.8.mol2"

cmd.load(f, 'label_threshold_11.8')
cmd.hide('everything', 'label_threshold_11.8')
cmd.label("label_threshold_11.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [11.8]
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


cluster_dict = {"13.5520000458":[], "13.5520000458_arrows":[]}

cluster_dict["13.5520000458"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(9.0), float(-14.0), float(1.0)]

cluster_dict["13.5520000458_arrows"] += cgo_arrow([24.5,9.0,-14.0], [22.615,11.131,-15.816], color="blue red", name="Arrows_13.5520000458_1")

cluster_dict["13.5520000458"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(18.0), float(-16.5), float(1.0)]

cluster_dict["13.5520000458_arrows"] += cgo_arrow([24.5,18.0,-16.5], [21.965,15.819,-16.753], color="blue red", name="Arrows_13.5520000458_2")

cluster_dict["13.5520000458"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(15.5), float(-14.5), float(1.0)]

cluster_dict["13.5520000458_arrows"] += cgo_arrow([24.5,15.5,-14.5], [21.438,14.463,-13.601], color="blue red", name="Arrows_13.5520000458_3")

cluster_dict["13.5520000458"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(3.5), float(-13.5), float(1.0)]

cluster_dict["13.5520000458_arrows"] += cgo_arrow([27.0,3.5,-13.5], [25.564,1.866,-16.048], color="blue red", name="Arrows_13.5520000458_4")

cluster_dict["13.5520000458"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(3.5), float(-13.0), float(1.0)]

cluster_dict["13.5520000458_arrows"] += cgo_arrow([30.5,3.5,-13.0], [30.724,3.981,-15.898], color="blue red", name="Arrows_13.5520000458_5")

cluster_dict["13.5520000458"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.0), float(12.0), float(-23.0), float(1.0)]

cluster_dict["13.5520000458_arrows"] += cgo_arrow([29.0,12.0,-23.0], [28.974,12.038,-25.36], color="blue red", name="Arrows_13.5520000458_6")

cluster_dict["13.5520000458"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(19.0), float(-15.5), float(1.0)]

cluster_dict["13.5520000458_arrows"] += cgo_arrow([29.5,19.0,-15.5], [28.582,16.633,-14.757], color="blue red", name="Arrows_13.5520000458_7")

cluster_dict["13.5520000458"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.1075844741), float(10.9534561524), float(-15.9279716566), float(1.0)]


cluster_dict["13.5520000458"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(31.2147926795), float(2.62964220704), float(-12.7101720482), float(1.0)]


cluster_dict["13.5520000458"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(32.729390681), float(19.9014336918), float(-15.3870967742), float(1.0)]


cluster_dict["13.5520000458"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(36.6941443325), float(6.92450307752), float(-12.5308123864), float(1.0)]


cluster_dict["13.5520000458"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(37.5), float(7.81049018446), float(-12.7930065437), float(1.0)]


cluster_dict["13.5520000458"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(28.5), float(11.0), float(-22.0), float(1.0)]

cluster_dict["13.5520000458_arrows"] += cgo_arrow([28.5,11.0,-22.0], [26.191,13.098,-22.189], color="red blue", name="Arrows_13.5520000458_8")

cluster_dict["13.5520000458"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(29.0), float(10.5), float(-12.0), float(1.0)]

cluster_dict["13.5520000458_arrows"] += cgo_arrow([29.0,10.5,-12.0], [29.581,9.727,-9.209], color="red blue", name="Arrows_13.5520000458_9")

cmd.load_cgo(cluster_dict["13.5520000458"], "Features_13.5520000458", 1)
cmd.load_cgo(cluster_dict["13.5520000458_arrows"], "Arrows_13.5520000458")
cmd.set("transparency", 0.2,"Features_13.5520000458")
cmd.group("Pharmacophore_13.5520000458", members="Features_13.5520000458")
cmd.group("Pharmacophore_13.5520000458", members="Arrows_13.5520000458")

if dirpath:
    f = join(dirpath, "10/label_threshold_13.5520000458.mol2")
else:
    f = "10/label_threshold_13.5520000458.mol2"

cmd.load(f, 'label_threshold_13.5520000458')
cmd.hide('everything', 'label_threshold_13.5520000458')
cmd.label("label_threshold_13.5520000458", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.5520000458', members= 'label_threshold_13.5520000458')


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
    f = join(dirpath, "11/label_threshold_11.0.mol2")
else:
    f = "11/label_threshold_11.0.mol2"

cmd.load(f, 'label_threshold_11.0')
cmd.hide('everything', 'label_threshold_11.0')
cmd.label("label_threshold_11.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [11.0]
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


cluster_dict = {"13.513999939":[], "13.513999939_arrows":[]}

cluster_dict["13.513999939"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(10.0), float(-14.5), float(1.0)]

cluster_dict["13.513999939_arrows"] += cgo_arrow([25.0,10.0,-14.5], [22.615,11.131,-15.816], color="blue red", name="Arrows_13.513999939_1")

cluster_dict["13.513999939"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(18.0), float(-16.5), float(1.0)]

cluster_dict["13.513999939_arrows"] += cgo_arrow([24.5,18.0,-16.5], [21.965,15.819,-16.753], color="blue red", name="Arrows_13.513999939_2")

cluster_dict["13.513999939"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(15.0), float(-14.0), float(1.0)]

cluster_dict["13.513999939_arrows"] += cgo_arrow([25.0,15.0,-14.0], [21.438,14.463,-13.601], color="blue red", name="Arrows_13.513999939_3")

cluster_dict["13.513999939"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(10.0), float(-14.5), float(1.0)]

cluster_dict["13.513999939_arrows"] += cgo_arrow([25.0,10.0,-14.5], [22.615,11.131,-15.816], color="blue red", name="Arrows_13.513999939_4")

cluster_dict["13.513999939"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(10.0), float(-13.0), float(1.0)]

cluster_dict["13.513999939_arrows"] += cgo_arrow([29.5,10.0,-13.0], [30.129,7.657,-14.523], color="blue red", name="Arrows_13.513999939_5")

cluster_dict["13.513999939"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(19.0), float(-15.5), float(1.0)]

cluster_dict["13.513999939_arrows"] += cgo_arrow([29.5,19.0,-15.5], [28.582,16.633,-14.757], color="blue red", name="Arrows_13.513999939_6")

cluster_dict["13.513999939"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(16.2574086929), float(21.4994176877), float(-10.5639713025), float(1.0)]


cluster_dict["13.513999939"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(17.778336332), float(11.8346477231), float(-11.1327536583), float(1.0)]


cluster_dict["13.513999939"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(20.125), float(9.625), float(-9.875), float(1.0)]


cluster_dict["13.513999939"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(26.7883294568), float(11.6693070681), float(-15.1034155519), float(1.0)]


cluster_dict["13.513999939"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.0879946273), float(11.3333892833), float(-14.4210342246), float(1.0)]


cluster_dict["13.513999939"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(29.1261266924), float(13.4731272176), float(-19.4139966071), float(1.0)]


cluster_dict["13.513999939"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(28.7951499191), float(10.7951499191), float(-19.0167783731), float(1.0)]


cluster_dict["13.513999939"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(32.8040053333), float(20.061087064), float(-15.3229776404), float(1.0)]


cluster_dict["13.513999939"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(10.5), float(-19.5), float(1.0)]

cluster_dict["13.513999939_arrows"] += cgo_arrow([26.0,10.5,-19.5], [23.813,11.613,-19.373], color="red blue", name="Arrows_13.513999939_7")

cluster_dict["13.513999939"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(12.5), float(-21.5), float(1.0)]

cluster_dict["13.513999939_arrows"] += cgo_arrow([28.0,12.5,-21.5], [26.191,13.098,-22.189], color="red blue", name="Arrows_13.513999939_8")

cluster_dict["13.513999939"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(29.0), float(10.5), float(-12.0), float(1.0)]

cluster_dict["13.513999939_arrows"] += cgo_arrow([29.0,10.5,-12.0], [29.581,9.727,-9.209], color="red blue", name="Arrows_13.513999939_9")

cmd.load_cgo(cluster_dict["13.513999939"], "Features_13.513999939", 1)
cmd.load_cgo(cluster_dict["13.513999939_arrows"], "Arrows_13.513999939")
cmd.set("transparency", 0.2,"Features_13.513999939")
cmd.group("Pharmacophore_13.513999939", members="Features_13.513999939")
cmd.group("Pharmacophore_13.513999939", members="Arrows_13.513999939")

if dirpath:
    f = join(dirpath, "11/label_threshold_13.513999939.mol2")
else:
    f = "11/label_threshold_13.513999939.mol2"

cmd.load(f, 'label_threshold_13.513999939')
cmd.hide('everything', 'label_threshold_13.513999939')
cmd.label("label_threshold_13.513999939", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.513999939', members= 'label_threshold_13.513999939')


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
    f = join(dirpath, "12/label_threshold_11.2.mol2")
else:
    f = "12/label_threshold_11.2.mol2"

cmd.load(f, 'label_threshold_11.2')
cmd.hide('everything', 'label_threshold_11.2')
cmd.label("label_threshold_11.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [11.2]
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


cluster_dict = {"13.4709997177":[], "13.4709997177_arrows":[]}

cluster_dict["13.4709997177"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(9.0), float(-14.0), float(1.0)]

cluster_dict["13.4709997177_arrows"] += cgo_arrow([24.5,9.0,-14.0], [22.615,11.131,-15.816], color="blue red", name="Arrows_13.4709997177_1")

cluster_dict["13.4709997177"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(3.5), float(-13.5), float(1.0)]

cluster_dict["13.4709997177_arrows"] += cgo_arrow([27.0,3.5,-13.5], [25.564,1.866,-16.048], color="blue red", name="Arrows_13.4709997177_2")

cluster_dict["13.4709997177"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.5), float(11.0), float(-18.5), float(1.0)]

cluster_dict["13.4709997177_arrows"] += cgo_arrow([27.5,11.0,-18.5], [25.291,13.09,-18.48], color="blue red", name="Arrows_13.4709997177_3")

cluster_dict["13.4709997177"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.0), float(9.5), float(-13.0), float(1.0)]

cluster_dict["13.4709997177_arrows"] += cgo_arrow([29.0,9.5,-13.0], [30.129,7.657,-14.523], color="blue red", name="Arrows_13.4709997177_4")

cluster_dict["13.4709997177"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(3.5), float(-13.0), float(1.0)]

cluster_dict["13.4709997177_arrows"] += cgo_arrow([30.5,3.5,-13.0], [30.724,3.981,-15.898], color="blue red", name="Arrows_13.4709997177_5")

cluster_dict["13.4709997177"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(14.7919923704), float(2.94441467867), float(-23.2558744251), float(1.0)]


cluster_dict["13.4709997177"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(17.9791649224), float(10.6206132726), float(-11.1845092282), float(1.0)]


cluster_dict["13.4709997177"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(18.2100452111), float(2.30523495221), float(-25.8819094404), float(1.0)]


cluster_dict["13.4709997177"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(26.2412757005), float(9.8311117214), float(-14.5112946889), float(1.0)]


cluster_dict["13.4709997177"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.0879946273), float(11.3333892833), float(-14.4210342246), float(1.0)]


cluster_dict["13.4709997177"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.8221112005), float(10.3059968621), float(-23.058168882), float(1.0)]


cluster_dict["13.4709997177"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(30.0554071364), float(1.58982395583), float(-11.9614392694), float(1.0)]


cluster_dict["13.4709997177"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(9.5), float(-21.0), float(1.0)]

cluster_dict["13.4709997177_arrows"] += cgo_arrow([26.0,9.5,-21.0], [23.813,11.613,-19.373], color="red blue", name="Arrows_13.4709997177_6")

cluster_dict["13.4709997177"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(10.5), float(-22.0), float(1.0)]

cluster_dict["13.4709997177_arrows"] += cgo_arrow([28.0,10.5,-22.0], [26.191,13.098,-22.189], color="red blue", name="Arrows_13.4709997177_7")

cmd.load_cgo(cluster_dict["13.4709997177"], "Features_13.4709997177", 1)
cmd.load_cgo(cluster_dict["13.4709997177_arrows"], "Arrows_13.4709997177")
cmd.set("transparency", 0.2,"Features_13.4709997177")
cmd.group("Pharmacophore_13.4709997177", members="Features_13.4709997177")
cmd.group("Pharmacophore_13.4709997177", members="Arrows_13.4709997177")

if dirpath:
    f = join(dirpath, "12/label_threshold_13.4709997177.mol2")
else:
    f = "12/label_threshold_13.4709997177.mol2"

cmd.load(f, 'label_threshold_13.4709997177')
cmd.hide('everything', 'label_threshold_13.4709997177')
cmd.label("label_threshold_13.4709997177", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.4709997177', members= 'label_threshold_13.4709997177')


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
    f = join(dirpath, "13/label_threshold_10.8.mol2")
else:
    f = "13/label_threshold_10.8.mol2"

cmd.load(f, 'label_threshold_10.8')
cmd.hide('everything', 'label_threshold_10.8')
cmd.label("label_threshold_10.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.8]
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


cluster_dict = {"13.3760004044":[], "13.3760004044_arrows":[]}

cluster_dict["13.3760004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(18.0), float(-16.5), float(1.0)]

cluster_dict["13.3760004044_arrows"] += cgo_arrow([24.5,18.0,-16.5], [21.965,15.819,-16.753], color="blue red", name="Arrows_13.3760004044_1")

cluster_dict["13.3760004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(10.0), float(-14.5), float(1.0)]

cluster_dict["13.3760004044_arrows"] += cgo_arrow([25.0,10.0,-14.5], [22.615,11.131,-15.816], color="blue red", name="Arrows_13.3760004044_2")

cluster_dict["13.3760004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(15.0), float(-14.0), float(1.0)]

cluster_dict["13.3760004044_arrows"] += cgo_arrow([25.0,15.0,-14.0], [21.438,14.463,-13.601], color="blue red", name="Arrows_13.3760004044_3")

cluster_dict["13.3760004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(10.0), float(-14.5), float(1.0)]

cluster_dict["13.3760004044_arrows"] += cgo_arrow([25.0,10.0,-14.5], [22.615,11.131,-15.816], color="blue red", name="Arrows_13.3760004044_4")

cluster_dict["13.3760004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.0), float(12.0), float(-23.0), float(1.0)]

cluster_dict["13.3760004044_arrows"] += cgo_arrow([29.0,12.0,-23.0], [28.974,12.038,-25.36], color="blue red", name="Arrows_13.3760004044_5")

cluster_dict["13.3760004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(19.0), float(-15.5), float(1.0)]

cluster_dict["13.3760004044_arrows"] += cgo_arrow([29.5,19.0,-15.5], [28.582,16.633,-14.757], color="blue red", name="Arrows_13.3760004044_6")

cluster_dict["13.3760004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(10.0), float(-13.0), float(1.0)]

cluster_dict["13.3760004044_arrows"] += cgo_arrow([29.5,10.0,-13.0], [30.129,7.657,-14.523], color="blue red", name="Arrows_13.3760004044_7")

cluster_dict["13.3760004044"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.0879946273), float(11.3333892833), float(-14.4210342246), float(1.0)]


cluster_dict["13.3760004044"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.8036372732), float(13.6643119936), float(-16.3908525317), float(1.0)]


cluster_dict["13.3760004044"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(32.5999214052), float(20.5193850444), float(-15.8275500931), float(1.0)]


cluster_dict["13.3760004044"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(35.6960157643), float(8.3153385881), float(-12.6306771762), float(1.0)]


cluster_dict["13.3760004044"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(36.3503127388), float(8.90667399882), float(-10.5184823706), float(1.0)]


cluster_dict["13.3760004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(10.0), float(-20.5), float(1.0)]

cluster_dict["13.3760004044_arrows"] += cgo_arrow([26.5,10.0,-20.5], [23.813,11.613,-19.373], color="red blue", name="Arrows_13.3760004044_8")

cluster_dict["13.3760004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(28.5), float(11.0), float(-22.0), float(1.0)]

cluster_dict["13.3760004044_arrows"] += cgo_arrow([28.5,11.0,-22.0], [26.191,13.098,-22.189], color="red blue", name="Arrows_13.3760004044_9")

cluster_dict["13.3760004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(29.0), float(10.5), float(-12.0), float(1.0)]

cluster_dict["13.3760004044_arrows"] += cgo_arrow([29.0,10.5,-12.0], [29.581,9.727,-9.209], color="red blue", name="Arrows_13.3760004044_10")

cmd.load_cgo(cluster_dict["13.3760004044"], "Features_13.3760004044", 1)
cmd.load_cgo(cluster_dict["13.3760004044_arrows"], "Arrows_13.3760004044")
cmd.set("transparency", 0.2,"Features_13.3760004044")
cmd.group("Pharmacophore_13.3760004044", members="Features_13.3760004044")
cmd.group("Pharmacophore_13.3760004044", members="Arrows_13.3760004044")

if dirpath:
    f = join(dirpath, "13/label_threshold_13.3760004044.mol2")
else:
    f = "13/label_threshold_13.3760004044.mol2"

cmd.load(f, 'label_threshold_13.3760004044')
cmd.hide('everything', 'label_threshold_13.3760004044')
cmd.label("label_threshold_13.3760004044", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.3760004044', members= 'label_threshold_13.3760004044')


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
    f = join(dirpath, "14/label_threshold_10.2.mol2")
else:
    f = "14/label_threshold_10.2.mol2"

cmd.load(f, 'label_threshold_10.2')
cmd.hide('everything', 'label_threshold_10.2')
cmd.label("label_threshold_10.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.2]
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


cluster_dict = {"13.3030004501":[], "13.3030004501_arrows":[]}

cluster_dict["13.3030004501"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(35.5), float(1.0), float(-10.0), float(1.0)]

cluster_dict["13.3030004501_arrows"] += cgo_arrow([35.5,1.0,-10.0], [36.43,-1.147,-11.123], color="blue red", name="Arrows_13.3030004501_1")

cluster_dict["13.3030004501"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(37.0), float(8.5), float(-12.0), float(1.0)]

cluster_dict["13.3030004501_arrows"] += cgo_arrow([37.0,8.5,-12.0], [34.648,10.111,-12.296], color="blue red", name="Arrows_13.3030004501_2")

cluster_dict["13.3030004501"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(37.0), float(5.5), float(-9.0), float(1.0)]

cluster_dict["13.3030004501_arrows"] += cgo_arrow([37.0,5.5,-9.0], [35.325,7.016,-7.626], color="blue red", name="Arrows_13.3030004501_3")

cluster_dict["13.3030004501"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(40.0), float(7.5), float(-5.0), float(1.0)]

cluster_dict["13.3030004501_arrows"] += cgo_arrow([40.0,7.5,-5.0], [37.094,7.651,-4.347], color="blue red", name="Arrows_13.3030004501_4")

cluster_dict["13.3030004501"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(40.9284896708), float(6.6172672037), float(-10.7047319756), float(1.0)]


cluster_dict["13.3030004501"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(6.0), float(-10.5), float(1.0)]

cluster_dict["13.3030004501_arrows"] += cgo_arrow([36.5,6.0,-10.5], [33.693,5.068,-9.683], color="red blue", name="Arrows_13.3030004501_5")

cluster_dict["13.3030004501"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(39.0), float(9.0), float(-6.5), float(1.0)]

cluster_dict["13.3030004501_arrows"] += cgo_arrow([39.0,9.0,-6.5], [34.921,9.085,-6.851], color="red blue", name="Arrows_13.3030004501_6")

cluster_dict["13.3030004501"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(43.5), float(7.0), float(-14.5), float(1.0)]

cluster_dict["13.3030004501_arrows"] += cgo_arrow([43.5,7.0,-14.5], [43.481,4.46,-13.067], color="red blue", name="Arrows_13.3030004501_7")

cluster_dict["13.3030004501"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(45.0), float(9.0), float(-19.0), float(1.0)]

cluster_dict["13.3030004501_arrows"] += cgo_arrow([45.0,9.0,-19.0], [45.862,12.152,-18.907], color="red blue", name="Arrows_13.3030004501_8")

cmd.load_cgo(cluster_dict["13.3030004501"], "Features_13.3030004501", 1)
cmd.load_cgo(cluster_dict["13.3030004501_arrows"], "Arrows_13.3030004501")
cmd.set("transparency", 0.2,"Features_13.3030004501")
cmd.group("Pharmacophore_13.3030004501", members="Features_13.3030004501")
cmd.group("Pharmacophore_13.3030004501", members="Arrows_13.3030004501")

if dirpath:
    f = join(dirpath, "14/label_threshold_13.3030004501.mol2")
else:
    f = "14/label_threshold_13.3030004501.mol2"

cmd.load(f, 'label_threshold_13.3030004501')
cmd.hide('everything', 'label_threshold_13.3030004501')
cmd.label("label_threshold_13.3030004501", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.3030004501', members= 'label_threshold_13.3030004501')


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
    f = join(dirpath, "15/label_threshold_7.9.mol2")
else:
    f = "15/label_threshold_7.9.mol2"

cmd.load(f, 'label_threshold_7.9')
cmd.hide('everything', 'label_threshold_7.9')
cmd.label("label_threshold_7.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [7.9]
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


cluster_dict = {"12.9289999008":[], "12.9289999008_arrows":[]}

cluster_dict["12.9289999008"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(-5.5), float(-36.5), float(1.0)]

cluster_dict["12.9289999008_arrows"] += cgo_arrow([29.5,-5.5,-36.5], [30.398,-6.958,-34.45], color="blue red", name="Arrows_12.9289999008_1")

cluster_dict["12.9289999008"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(32.5), float(3.0), float(-37.0), float(1.0)]

cluster_dict["12.9289999008_arrows"] += cgo_arrow([32.5,3.0,-37.0], [32.104,1.695,-34.881], color="blue red", name="Arrows_12.9289999008_2")

cluster_dict["12.9289999008"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(33.5), float(3.5), float(-37.0), float(1.0)]

cluster_dict["12.9289999008_arrows"] += cgo_arrow([33.5,3.5,-37.0], [34.099,2.438,-34.389], color="blue red", name="Arrows_12.9289999008_3")

cluster_dict["12.9289999008"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(37.0), float(2.0), float(-42.5), float(1.0)]

cluster_dict["12.9289999008_arrows"] += cgo_arrow([37.0,2.0,-42.5], [38.722,3.854,-44.392], color="blue red", name="Arrows_12.9289999008_4")

cluster_dict["12.9289999008"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(37.0), float(2.0), float(-42.5), float(1.0)]

cluster_dict["12.9289999008_arrows"] += cgo_arrow([37.0,2.0,-42.5], [38.722,3.854,-44.392], color="blue red", name="Arrows_12.9289999008_5")

cluster_dict["12.9289999008"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(31.2173887131), float(-3.41112217664), float(-37.8178804902), float(1.0)]


cluster_dict["12.9289999008"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(35.7982789416), float(1.7357446839), float(-40.3198843885), float(1.0)]


cmd.load_cgo(cluster_dict["12.9289999008"], "Features_12.9289999008", 1)
cmd.load_cgo(cluster_dict["12.9289999008_arrows"], "Arrows_12.9289999008")
cmd.set("transparency", 0.2,"Features_12.9289999008")
cmd.group("Pharmacophore_12.9289999008", members="Features_12.9289999008")
cmd.group("Pharmacophore_12.9289999008", members="Arrows_12.9289999008")

if dirpath:
    f = join(dirpath, "15/label_threshold_12.9289999008.mol2")
else:
    f = "15/label_threshold_12.9289999008.mol2"

cmd.load(f, 'label_threshold_12.9289999008')
cmd.hide('everything', 'label_threshold_12.9289999008')
cmd.label("label_threshold_12.9289999008", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.9289999008', members= 'label_threshold_12.9289999008')


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

if dirpath:
    f = join(dirpath, "16/label_threshold_8.5.mol2")
else:
    f = "16/label_threshold_8.5.mol2"

cmd.load(f, 'label_threshold_8.5')
cmd.hide('everything', 'label_threshold_8.5')
cmd.label("label_threshold_8.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [8.5]
gfiles = ['16/donor.grd', '16/apolar.grd', '16/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 16
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


cluster_dict = {"12.1405000687":[], "12.1405000687_arrows":[]}

cluster_dict["12.1405000687"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(11.6056466061), float(16.7199830985), float(-11.4804785106), float(1.0)]


cluster_dict["12.1405000687"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(17.2068563959), float(11.7608595524), float(-11.3663119886), float(1.0)]


cmd.load_cgo(cluster_dict["12.1405000687"], "Features_12.1405000687", 1)
cmd.load_cgo(cluster_dict["12.1405000687_arrows"], "Arrows_12.1405000687")
cmd.set("transparency", 0.2,"Features_12.1405000687")
cmd.group("Pharmacophore_12.1405000687", members="Features_12.1405000687")
cmd.group("Pharmacophore_12.1405000687", members="Arrows_12.1405000687")

if dirpath:
    f = join(dirpath, "16/label_threshold_12.1405000687.mol2")
else:
    f = "16/label_threshold_12.1405000687.mol2"

cmd.load(f, 'label_threshold_12.1405000687')
cmd.hide('everything', 'label_threshold_12.1405000687')
cmd.label("label_threshold_12.1405000687", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.1405000687', members= 'label_threshold_12.1405000687')


if dirpath:
    f = join(dirpath, '16/mesh.grd')
else:
    f = '16/mesh.grd'
cmd.load(f, 'mesh_16')
cmd.isomesh("isomesh_16", "mesh_16", 0.9)
cmd.color("grey80", "isomesh_16")
cmd.set('transparency', 0.4, "isomesh_16")

cmd.group('hotspot_16', "isomesh_16")
cmd.group('hotspot_16', "mesh_16")

if dirpath:
    f = join(dirpath, "17/label_threshold_8.2.mol2")
else:
    f = "17/label_threshold_8.2.mol2"

cmd.load(f, 'label_threshold_8.2')
cmd.hide('everything', 'label_threshold_8.2')
cmd.label("label_threshold_8.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [8.2]
gfiles = ['17/donor.grd', '17/apolar.grd', '17/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 17
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


cluster_dict = {"12.1300001144":[], "12.1300001144_arrows":[]}

cluster_dict["12.1300001144"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(11.6194764837), float(16.7184451404), float(-11.4721746093), float(1.0)]


cluster_dict["12.1300001144"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(17.2074512118), float(11.7644988726), float(-11.3630896949), float(1.0)]


cmd.load_cgo(cluster_dict["12.1300001144"], "Features_12.1300001144", 1)
cmd.load_cgo(cluster_dict["12.1300001144_arrows"], "Arrows_12.1300001144")
cmd.set("transparency", 0.2,"Features_12.1300001144")
cmd.group("Pharmacophore_12.1300001144", members="Features_12.1300001144")
cmd.group("Pharmacophore_12.1300001144", members="Arrows_12.1300001144")

if dirpath:
    f = join(dirpath, "17/label_threshold_12.1300001144.mol2")
else:
    f = "17/label_threshold_12.1300001144.mol2"

cmd.load(f, 'label_threshold_12.1300001144')
cmd.hide('everything', 'label_threshold_12.1300001144')
cmd.label("label_threshold_12.1300001144", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.1300001144', members= 'label_threshold_12.1300001144')


if dirpath:
    f = join(dirpath, '17/mesh.grd')
else:
    f = '17/mesh.grd'
cmd.load(f, 'mesh_17')
cmd.isomesh("isomesh_17", "mesh_17", 0.9)
cmd.color("grey80", "isomesh_17")
cmd.set('transparency', 0.4, "isomesh_17")

cmd.group('hotspot_17', "isomesh_17")
cmd.group('hotspot_17', "mesh_17")

if dirpath:
    f = join(dirpath, "18/label_threshold_8.7.mol2")
else:
    f = "18/label_threshold_8.7.mol2"

cmd.load(f, 'label_threshold_8.7')
cmd.hide('everything', 'label_threshold_8.7')
cmd.label("label_threshold_8.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [8.7]
gfiles = ['18/donor.grd', '18/apolar.grd', '18/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 18
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


cluster_dict = {"10.8459997177":[], "10.8459997177_arrows":[]}

cluster_dict["10.8459997177"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(18.0), float(-16.5), float(1.0)]

cluster_dict["10.8459997177_arrows"] += cgo_arrow([24.5,18.0,-16.5], [21.965,15.819,-16.753], color="blue red", name="Arrows_10.8459997177_1")

cluster_dict["10.8459997177"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(15.5), float(-14.5), float(1.0)]

cluster_dict["10.8459997177_arrows"] += cgo_arrow([24.5,15.5,-14.5], [21.438,14.463,-13.601], color="blue red", name="Arrows_10.8459997177_2")

cluster_dict["10.8459997177"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(18.0), float(-16.5), float(1.0)]

cluster_dict["10.8459997177_arrows"] += cgo_arrow([24.5,18.0,-16.5], [21.965,15.819,-16.753], color="blue red", name="Arrows_10.8459997177_3")

cluster_dict["10.8459997177"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(18.5), float(-16.0), float(1.0)]

cluster_dict["10.8459997177_arrows"] += cgo_arrow([29.5,18.5,-16.0], [28.582,16.633,-14.757], color="blue red", name="Arrows_10.8459997177_4")

cluster_dict["10.8459997177"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(28.8416667258), float(18.6495045342), float(-15.3102002975), float(1.0)]


cluster_dict["10.8459997177"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(29.1443146622), float(14.5124927716), float(-19.4226377962), float(1.0)]


cmd.load_cgo(cluster_dict["10.8459997177"], "Features_10.8459997177", 1)
cmd.load_cgo(cluster_dict["10.8459997177_arrows"], "Arrows_10.8459997177")
cmd.set("transparency", 0.2,"Features_10.8459997177")
cmd.group("Pharmacophore_10.8459997177", members="Features_10.8459997177")
cmd.group("Pharmacophore_10.8459997177", members="Arrows_10.8459997177")

if dirpath:
    f = join(dirpath, "18/label_threshold_10.8459997177.mol2")
else:
    f = "18/label_threshold_10.8459997177.mol2"

cmd.load(f, 'label_threshold_10.8459997177')
cmd.hide('everything', 'label_threshold_10.8459997177')
cmd.label("label_threshold_10.8459997177", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.8459997177', members= 'label_threshold_10.8459997177')


if dirpath:
    f = join(dirpath, '18/mesh.grd')
else:
    f = '18/mesh.grd'
cmd.load(f, 'mesh_18')
cmd.isomesh("isomesh_18", "mesh_18", 0.9)
cmd.color("grey80", "isomesh_18")
cmd.set('transparency', 0.4, "isomesh_18")

cmd.group('hotspot_18', "isomesh_18")
cmd.group('hotspot_18', "mesh_18")

if dirpath:
    f = join(dirpath, "19/label_threshold_5.5.mol2")
else:
    f = "19/label_threshold_5.5.mol2"

cmd.load(f, 'label_threshold_5.5')
cmd.hide('everything', 'label_threshold_5.5')
cmd.label("label_threshold_5.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [5.5]
gfiles = ['19/donor.grd', '19/apolar.grd', '19/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 19
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


cluster_dict = {"10.1459999084":[], "10.1459999084_arrows":[]}

cluster_dict["10.1459999084"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(20.3884858603), float(21.0294068839), float(15.3286451045), float(1.0)]


cluster_dict["10.1459999084"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.0), float(17.0), float(13.0), float(1.0)]

cluster_dict["10.1459999084_arrows"] += cgo_arrow([17.0,17.0,13.0], [18.463,16.048,11.208], color="red blue", name="Arrows_10.1459999084_1")

cmd.load_cgo(cluster_dict["10.1459999084"], "Features_10.1459999084", 1)
cmd.load_cgo(cluster_dict["10.1459999084_arrows"], "Arrows_10.1459999084")
cmd.set("transparency", 0.2,"Features_10.1459999084")
cmd.group("Pharmacophore_10.1459999084", members="Features_10.1459999084")
cmd.group("Pharmacophore_10.1459999084", members="Arrows_10.1459999084")

if dirpath:
    f = join(dirpath, "19/label_threshold_10.1459999084.mol2")
else:
    f = "19/label_threshold_10.1459999084.mol2"

cmd.load(f, 'label_threshold_10.1459999084')
cmd.hide('everything', 'label_threshold_10.1459999084')
cmd.label("label_threshold_10.1459999084", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.1459999084', members= 'label_threshold_10.1459999084')


if dirpath:
    f = join(dirpath, '19/mesh.grd')
else:
    f = '19/mesh.grd'
cmd.load(f, 'mesh_19')
cmd.isomesh("isomesh_19", "mesh_19", 0.9)
cmd.color("grey80", "isomesh_19")
cmd.set('transparency', 0.4, "isomesh_19")

cmd.group('hotspot_19', "isomesh_19")
cmd.group('hotspot_19', "mesh_19")

if dirpath:
    f = join(dirpath, "20/label_threshold_5.5.mol2")
else:
    f = "20/label_threshold_5.5.mol2"

cmd.load(f, 'label_threshold_5.5')
cmd.hide('everything', 'label_threshold_5.5')
cmd.label("label_threshold_5.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [5.5]
gfiles = ['20/donor.grd', '20/apolar.grd', '20/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 20
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


cluster_dict = {"10.140999794":[], "10.140999794_arrows":[]}

cluster_dict["10.140999794"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(20.3898017959), float(21.0289704977), float(15.3277665236), float(1.0)]


cluster_dict["10.140999794"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.0), float(17.0), float(13.0), float(1.0)]

cluster_dict["10.140999794_arrows"] += cgo_arrow([17.0,17.0,13.0], [18.463,16.048,11.208], color="red blue", name="Arrows_10.140999794_1")

cmd.load_cgo(cluster_dict["10.140999794"], "Features_10.140999794", 1)
cmd.load_cgo(cluster_dict["10.140999794_arrows"], "Arrows_10.140999794")
cmd.set("transparency", 0.2,"Features_10.140999794")
cmd.group("Pharmacophore_10.140999794", members="Features_10.140999794")
cmd.group("Pharmacophore_10.140999794", members="Arrows_10.140999794")

if dirpath:
    f = join(dirpath, "20/label_threshold_10.140999794.mol2")
else:
    f = "20/label_threshold_10.140999794.mol2"

cmd.load(f, 'label_threshold_10.140999794')
cmd.hide('everything', 'label_threshold_10.140999794')
cmd.label("label_threshold_10.140999794", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.140999794', members= 'label_threshold_10.140999794')


if dirpath:
    f = join(dirpath, '20/mesh.grd')
else:
    f = '20/mesh.grd'
cmd.load(f, 'mesh_20')
cmd.isomesh("isomesh_20", "mesh_20", 0.9)
cmd.color("grey80", "isomesh_20")
cmd.set('transparency', 0.4, "isomesh_20")

cmd.group('hotspot_20', "isomesh_20")
cmd.group('hotspot_20', "mesh_20")

if dirpath:
    f = join(dirpath, "21/label_threshold_0.6.mol2")
else:
    f = "21/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
gfiles = ['21/donor.grd', '21/apolar.grd', '21/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 21
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


cluster_dict = {"9.94499969482":[], "9.94499969482_arrows":[]}

cluster_dict["9.94499969482"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(20.0934461006), float(35.5768953771), float(-41.6276216917), float(1.0)]


cmd.load_cgo(cluster_dict["9.94499969482"], "Features_9.94499969482", 1)
cmd.load_cgo(cluster_dict["9.94499969482_arrows"], "Arrows_9.94499969482")
cmd.set("transparency", 0.2,"Features_9.94499969482")
cmd.group("Pharmacophore_9.94499969482", members="Features_9.94499969482")
cmd.group("Pharmacophore_9.94499969482", members="Arrows_9.94499969482")

if dirpath:
    f = join(dirpath, "21/label_threshold_9.94499969482.mol2")
else:
    f = "21/label_threshold_9.94499969482.mol2"

cmd.load(f, 'label_threshold_9.94499969482')
cmd.hide('everything', 'label_threshold_9.94499969482')
cmd.label("label_threshold_9.94499969482", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.94499969482', members= 'label_threshold_9.94499969482')


if dirpath:
    f = join(dirpath, '21/mesh.grd')
else:
    f = '21/mesh.grd'
cmd.load(f, 'mesh_21')
cmd.isomesh("isomesh_21", "mesh_21", 0.9)
cmd.color("grey80", "isomesh_21")
cmd.set('transparency', 0.4, "isomesh_21")

cmd.group('hotspot_21', "isomesh_21")
cmd.group('hotspot_21', "mesh_21")

if dirpath:
    f = join(dirpath, "22/label_threshold_8.5.mol2")
else:
    f = "22/label_threshold_8.5.mol2"

cmd.load(f, 'label_threshold_8.5')
cmd.hide('everything', 'label_threshold_8.5')
cmd.label("label_threshold_8.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [8.5]
gfiles = ['22/donor.grd', '22/apolar.grd', '22/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 22
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


cluster_dict = {"9.93999958038":[], "9.93999958038_arrows":[]}

cluster_dict["9.93999958038"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(7.5), float(-13.5), float(1.0)]

cluster_dict["9.93999958038_arrows"] += cgo_arrow([23.5,7.5,-13.5], [23.938,5.884,-15.379], color="blue red", name="Arrows_9.93999958038_1")

cluster_dict["9.93999958038"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(10.8560884192), float(12.4617696739), float(-11.1975486251), float(1.0)]


cluster_dict["9.93999958038"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(12.8064011701), float(1.16553171322), float(-14.6967707156), float(1.0)]


cluster_dict["9.93999958038"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(14.6814267054), float(3.02422068204), float(-23.1804552679), float(1.0)]


cluster_dict["9.93999958038"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(19.3508348859), float(6.17608069851), float(-13.8744613728), float(1.0)]


cluster_dict["9.93999958038"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(17.9851390887), float(2.22967200241), float(-25.8765217755), float(1.0)]


cluster_dict["9.93999958038"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(8.5), float(-14.0), float(1.0)]


cmd.load_cgo(cluster_dict["9.93999958038"], "Features_9.93999958038", 1)
cmd.load_cgo(cluster_dict["9.93999958038_arrows"], "Arrows_9.93999958038")
cmd.set("transparency", 0.2,"Features_9.93999958038")
cmd.group("Pharmacophore_9.93999958038", members="Features_9.93999958038")
cmd.group("Pharmacophore_9.93999958038", members="Arrows_9.93999958038")

if dirpath:
    f = join(dirpath, "22/label_threshold_9.93999958038.mol2")
else:
    f = "22/label_threshold_9.93999958038.mol2"

cmd.load(f, 'label_threshold_9.93999958038')
cmd.hide('everything', 'label_threshold_9.93999958038')
cmd.label("label_threshold_9.93999958038", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.93999958038', members= 'label_threshold_9.93999958038')


if dirpath:
    f = join(dirpath, '22/mesh.grd')
else:
    f = '22/mesh.grd'
cmd.load(f, 'mesh_22')
cmd.isomesh("isomesh_22", "mesh_22", 0.9)
cmd.color("grey80", "isomesh_22")
cmd.set('transparency', 0.4, "isomesh_22")

cmd.group('hotspot_22', "isomesh_22")
cmd.group('hotspot_22', "mesh_22")

if dirpath:
    f = join(dirpath, "23/label_threshold_0.6.mol2")
else:
    f = "23/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
gfiles = ['23/donor.grd', '23/apolar.grd', '23/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 23
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


cluster_dict = {"9.78800010681":[], "9.78800010681_arrows":[]}

cluster_dict["9.78800010681"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(36.563540767), float(25.9684804526), float(-49.7340162443), float(1.0)]


cluster_dict["9.78800010681"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(35.5), float(23.5), float(-50.0), float(1.0)]

cluster_dict["9.78800010681_arrows"] += cgo_arrow([35.5,23.5,-50.0], [30.875,23.44,-49.917], color="red blue", name="Arrows_9.78800010681_1")

cluster_dict["9.78800010681"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(35.5), float(23.5), float(-50.0), float(1.0)]

cluster_dict["9.78800010681_arrows"] += cgo_arrow([35.5,23.5,-50.0], [30.875,23.44,-49.917], color="red blue", name="Arrows_9.78800010681_2")

cmd.load_cgo(cluster_dict["9.78800010681"], "Features_9.78800010681", 1)
cmd.load_cgo(cluster_dict["9.78800010681_arrows"], "Arrows_9.78800010681")
cmd.set("transparency", 0.2,"Features_9.78800010681")
cmd.group("Pharmacophore_9.78800010681", members="Features_9.78800010681")
cmd.group("Pharmacophore_9.78800010681", members="Arrows_9.78800010681")

if dirpath:
    f = join(dirpath, "23/label_threshold_9.78800010681.mol2")
else:
    f = "23/label_threshold_9.78800010681.mol2"

cmd.load(f, 'label_threshold_9.78800010681')
cmd.hide('everything', 'label_threshold_9.78800010681')
cmd.label("label_threshold_9.78800010681", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.78800010681', members= 'label_threshold_9.78800010681')


if dirpath:
    f = join(dirpath, '23/mesh.grd')
else:
    f = '23/mesh.grd'
cmd.load(f, 'mesh_23')
cmd.isomesh("isomesh_23", "mesh_23", 0.9)
cmd.color("grey80", "isomesh_23")
cmd.set('transparency', 0.4, "isomesh_23")

cmd.group('hotspot_23', "isomesh_23")
cmd.group('hotspot_23', "mesh_23")

if dirpath:
    f = join(dirpath, "24/label_threshold_6.5.mol2")
else:
    f = "24/label_threshold_6.5.mol2"

cmd.load(f, 'label_threshold_6.5')
cmd.hide('everything', 'label_threshold_6.5')
cmd.label("label_threshold_6.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [6.5]
gfiles = ['24/donor.grd', '24/apolar.grd', '24/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 24
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


cluster_dict = {"9.32499980927":[], "9.32499980927_arrows":[]}

cluster_dict["9.32499980927"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(13.2805250902), float(1.33244814669), float(-14.7237022282), float(1.0)]


cluster_dict["9.32499980927"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(14.6669559617), float(3.03350590622), float(-23.1978770309), float(1.0)]


cluster_dict["9.32499980927"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(19.4778423594), float(4.02254403081), float(-15.7873119544), float(1.0)]


cluster_dict["9.32499980927"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(16.9463437584), float(11.0517295464), float(-14.2631284158), float(1.0)]


cluster_dict["9.32499980927"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(18.0374943341), float(2.41669314238), float(-25.8962324579), float(1.0)]


cluster_dict["9.32499980927"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.0566270012), float(8.00904214928), float(-14.7719678788), float(1.0)]


cluster_dict["9.32499980927"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.2648660692), float(8.44126059744), float(-17.7621428367), float(1.0)]


cluster_dict["9.32499980927"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(25.1856680723), float(8.32131392649), float(-19.4983332071), float(1.0)]


cmd.load_cgo(cluster_dict["9.32499980927"], "Features_9.32499980927", 1)
cmd.load_cgo(cluster_dict["9.32499980927_arrows"], "Arrows_9.32499980927")
cmd.set("transparency", 0.2,"Features_9.32499980927")
cmd.group("Pharmacophore_9.32499980927", members="Features_9.32499980927")
cmd.group("Pharmacophore_9.32499980927", members="Arrows_9.32499980927")

if dirpath:
    f = join(dirpath, "24/label_threshold_9.32499980927.mol2")
else:
    f = "24/label_threshold_9.32499980927.mol2"

cmd.load(f, 'label_threshold_9.32499980927')
cmd.hide('everything', 'label_threshold_9.32499980927')
cmd.label("label_threshold_9.32499980927", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.32499980927', members= 'label_threshold_9.32499980927')


if dirpath:
    f = join(dirpath, '24/mesh.grd')
else:
    f = '24/mesh.grd'
cmd.load(f, 'mesh_24')
cmd.isomesh("isomesh_24", "mesh_24", 0.9)
cmd.color("grey80", "isomesh_24")
cmd.set('transparency', 0.4, "isomesh_24")

cmd.group('hotspot_24', "isomesh_24")
cmd.group('hotspot_24', "mesh_24")

if dirpath:
    f = join(dirpath, "25/label_threshold_6.7.mol2")
else:
    f = "25/label_threshold_6.7.mol2"

cmd.load(f, 'label_threshold_6.7')
cmd.hide('everything', 'label_threshold_6.7')
cmd.label("label_threshold_6.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [6.7]
gfiles = ['25/donor.grd', '25/apolar.grd', '25/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 25
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


cluster_dict = {"9.08399963379":[], "9.08399963379_arrows":[]}

cluster_dict["9.08399963379"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(18.0), float(-16.5), float(1.0)]

cluster_dict["9.08399963379_arrows"] += cgo_arrow([24.5,18.0,-16.5], [21.965,15.819,-16.753], color="blue red", name="Arrows_9.08399963379_1")

cluster_dict["9.08399963379"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(15.5), float(-14.5), float(1.0)]

cluster_dict["9.08399963379_arrows"] += cgo_arrow([24.5,15.5,-14.5], [21.438,14.463,-13.601], color="blue red", name="Arrows_9.08399963379_2")

cluster_dict["9.08399963379"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(19.0), float(-15.5), float(1.0)]

cluster_dict["9.08399963379_arrows"] += cgo_arrow([29.5,19.0,-15.5], [28.582,16.633,-14.757], color="blue red", name="Arrows_9.08399963379_3")

cluster_dict["9.08399963379"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.7042707357), float(12.5), float(-7.03010676839), float(1.0)]


cluster_dict["9.08399963379"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(26.8277133001), float(19.3130303622), float(-13.3714018781), float(1.0)]


cluster_dict["9.08399963379"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(25.6108091793), float(25.0885116929), float(-18.8728236417), float(1.0)]


cmd.load_cgo(cluster_dict["9.08399963379"], "Features_9.08399963379", 1)
cmd.load_cgo(cluster_dict["9.08399963379_arrows"], "Arrows_9.08399963379")
cmd.set("transparency", 0.2,"Features_9.08399963379")
cmd.group("Pharmacophore_9.08399963379", members="Features_9.08399963379")
cmd.group("Pharmacophore_9.08399963379", members="Arrows_9.08399963379")

if dirpath:
    f = join(dirpath, "25/label_threshold_9.08399963379.mol2")
else:
    f = "25/label_threshold_9.08399963379.mol2"

cmd.load(f, 'label_threshold_9.08399963379')
cmd.hide('everything', 'label_threshold_9.08399963379')
cmd.label("label_threshold_9.08399963379", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.08399963379', members= 'label_threshold_9.08399963379')


if dirpath:
    f = join(dirpath, '25/mesh.grd')
else:
    f = '25/mesh.grd'
cmd.load(f, 'mesh_25')
cmd.isomesh("isomesh_25", "mesh_25", 0.9)
cmd.color("grey80", "isomesh_25")
cmd.set('transparency', 0.4, "isomesh_25")

cmd.group('hotspot_25', "isomesh_25")
cmd.group('hotspot_25', "mesh_25")

if dirpath:
    f = join(dirpath, "26/label_threshold_0.6.mol2")
else:
    f = "26/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
gfiles = ['26/donor.grd', '26/apolar.grd', '26/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 26
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


cluster_dict = {"8.63099956512":[], "8.63099956512_arrows":[]}

cluster_dict["8.63099956512"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(9.0), float(-0.5), float(29.0), float(1.0)]

cluster_dict["8.63099956512_arrows"] += cgo_arrow([9.0,-0.5,29.0], [10.978,-2.707,28.965], color="blue red", name="Arrows_8.63099956512_1")

cluster_dict["8.63099956512"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(8.97777302477), float(0.366303484835), float(30.1340558783), float(1.0)]


cmd.load_cgo(cluster_dict["8.63099956512"], "Features_8.63099956512", 1)
cmd.load_cgo(cluster_dict["8.63099956512_arrows"], "Arrows_8.63099956512")
cmd.set("transparency", 0.2,"Features_8.63099956512")
cmd.group("Pharmacophore_8.63099956512", members="Features_8.63099956512")
cmd.group("Pharmacophore_8.63099956512", members="Arrows_8.63099956512")

if dirpath:
    f = join(dirpath, "26/label_threshold_8.63099956512.mol2")
else:
    f = "26/label_threshold_8.63099956512.mol2"

cmd.load(f, 'label_threshold_8.63099956512')
cmd.hide('everything', 'label_threshold_8.63099956512')
cmd.label("label_threshold_8.63099956512", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_8.63099956512', members= 'label_threshold_8.63099956512')


if dirpath:
    f = join(dirpath, '26/mesh.grd')
else:
    f = '26/mesh.grd'
cmd.load(f, 'mesh_26')
cmd.isomesh("isomesh_26", "mesh_26", 0.9)
cmd.color("grey80", "isomesh_26")
cmd.set('transparency', 0.4, "isomesh_26")

cmd.group('hotspot_26', "isomesh_26")
cmd.group('hotspot_26', "mesh_26")

if dirpath:
    f = join(dirpath, "27/label_threshold_5.9.mol2")
else:
    f = "27/label_threshold_5.9.mol2"

cmd.load(f, 'label_threshold_5.9')
cmd.hide('everything', 'label_threshold_5.9')
cmd.label("label_threshold_5.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [5.9]
gfiles = ['27/donor.grd', '27/apolar.grd', '27/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 27
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


cluster_dict = {"8.55000019073":[], "8.55000019073_arrows":[]}

cluster_dict["8.55000019073"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(17.1045226384), float(1.7233766205), float(-15.02551805), float(1.0)]


cluster_dict["8.55000019073"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(14.5886836763), float(2.16581717415), float(-22.0596139534), float(1.0)]


cmd.load_cgo(cluster_dict["8.55000019073"], "Features_8.55000019073", 1)
cmd.load_cgo(cluster_dict["8.55000019073_arrows"], "Arrows_8.55000019073")
cmd.set("transparency", 0.2,"Features_8.55000019073")
cmd.group("Pharmacophore_8.55000019073", members="Features_8.55000019073")
cmd.group("Pharmacophore_8.55000019073", members="Arrows_8.55000019073")

if dirpath:
    f = join(dirpath, "27/label_threshold_8.55000019073.mol2")
else:
    f = "27/label_threshold_8.55000019073.mol2"

cmd.load(f, 'label_threshold_8.55000019073')
cmd.hide('everything', 'label_threshold_8.55000019073')
cmd.label("label_threshold_8.55000019073", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_8.55000019073', members= 'label_threshold_8.55000019073')


if dirpath:
    f = join(dirpath, '27/mesh.grd')
else:
    f = '27/mesh.grd'
cmd.load(f, 'mesh_27')
cmd.isomesh("isomesh_27", "mesh_27", 0.9)
cmd.color("grey80", "isomesh_27")
cmd.set('transparency', 0.4, "isomesh_27")

cmd.group('hotspot_27', "isomesh_27")
cmd.group('hotspot_27', "mesh_27")

if dirpath:
    f = join(dirpath, "28/label_threshold_0.6.mol2")
else:
    f = "28/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
gfiles = ['28/donor.grd', '28/apolar.grd', '28/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 28
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


cluster_dict = {"8.49899959564":[], "8.49899959564_arrows":[]}

cluster_dict["8.49899959564"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(10.3462642106), float(20.116764814), float(-27.622894219), float(1.0)]


cluster_dict["8.49899959564"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(11.5011030684), float(23.2409131253), float(-28.3382254037), float(1.0)]


cmd.load_cgo(cluster_dict["8.49899959564"], "Features_8.49899959564", 1)
cmd.load_cgo(cluster_dict["8.49899959564_arrows"], "Arrows_8.49899959564")
cmd.set("transparency", 0.2,"Features_8.49899959564")
cmd.group("Pharmacophore_8.49899959564", members="Features_8.49899959564")
cmd.group("Pharmacophore_8.49899959564", members="Arrows_8.49899959564")

if dirpath:
    f = join(dirpath, "28/label_threshold_8.49899959564.mol2")
else:
    f = "28/label_threshold_8.49899959564.mol2"

cmd.load(f, 'label_threshold_8.49899959564')
cmd.hide('everything', 'label_threshold_8.49899959564')
cmd.label("label_threshold_8.49899959564", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_8.49899959564', members= 'label_threshold_8.49899959564')


if dirpath:
    f = join(dirpath, '28/mesh.grd')
else:
    f = '28/mesh.grd'
cmd.load(f, 'mesh_28')
cmd.isomesh("isomesh_28", "mesh_28", 0.9)
cmd.color("grey80", "isomesh_28")
cmd.set('transparency', 0.4, "isomesh_28")

cmd.group('hotspot_28', "isomesh_28")
cmd.group('hotspot_28', "mesh_28")

if dirpath:
    f = join(dirpath, "29/label_threshold_6.3.mol2")
else:
    f = "29/label_threshold_6.3.mol2"

cmd.load(f, 'label_threshold_6.3')
cmd.hide('everything', 'label_threshold_6.3')
cmd.label("label_threshold_6.3", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [6.3]
gfiles = ['29/donor.grd', '29/apolar.grd', '29/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 29
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


cluster_dict = {"7.89699983597":[], "7.89699983597_arrows":[]}

cluster_dict["7.89699983597"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(7.5), float(0.5), float(-0.5), float(1.0)]

cluster_dict["7.89699983597_arrows"] += cgo_arrow([7.5,0.5,-0.5], [7.832,3.488,-1.705], color="blue red", name="Arrows_7.89699983597_1")

cluster_dict["7.89699983597"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(8.83059662616), float(1.50568743153), float(0.605459453919), float(1.0)]


cmd.load_cgo(cluster_dict["7.89699983597"], "Features_7.89699983597", 1)
cmd.load_cgo(cluster_dict["7.89699983597_arrows"], "Arrows_7.89699983597")
cmd.set("transparency", 0.2,"Features_7.89699983597")
cmd.group("Pharmacophore_7.89699983597", members="Features_7.89699983597")
cmd.group("Pharmacophore_7.89699983597", members="Arrows_7.89699983597")

if dirpath:
    f = join(dirpath, "29/label_threshold_7.89699983597.mol2")
else:
    f = "29/label_threshold_7.89699983597.mol2"

cmd.load(f, 'label_threshold_7.89699983597')
cmd.hide('everything', 'label_threshold_7.89699983597')
cmd.label("label_threshold_7.89699983597", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_7.89699983597', members= 'label_threshold_7.89699983597')


if dirpath:
    f = join(dirpath, '29/mesh.grd')
else:
    f = '29/mesh.grd'
cmd.load(f, 'mesh_29')
cmd.isomesh("isomesh_29", "mesh_29", 0.9)
cmd.color("grey80", "isomesh_29")
cmd.set('transparency', 0.4, "isomesh_29")

cmd.group('hotspot_29', "isomesh_29")
cmd.group('hotspot_29', "mesh_29")

if dirpath:
    f = join(dirpath, "30/label_threshold_6.2.mol2")
else:
    f = "30/label_threshold_6.2.mol2"

cmd.load(f, 'label_threshold_6.2')
cmd.hide('everything', 'label_threshold_6.2')
cmd.label("label_threshold_6.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [6.2]
gfiles = ['30/donor.grd', '30/apolar.grd', '30/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 30
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


cluster_dict = {"7.85300016403":[], "7.85300016403_arrows":[]}

cluster_dict["7.85300016403"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(7.5), float(0.5), float(-0.5), float(1.0)]

cluster_dict["7.85300016403_arrows"] += cgo_arrow([7.5,0.5,-0.5], [7.832,3.488,-1.705], color="blue red", name="Arrows_7.85300016403_1")

cluster_dict["7.85300016403"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(0.0), float(0.25), float(-1.77970371302), float(1.0)]


cluster_dict["7.85300016403"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(9.08079176074), float(1.62514962159), float(0.670560128512), float(1.0)]


cmd.load_cgo(cluster_dict["7.85300016403"], "Features_7.85300016403", 1)
cmd.load_cgo(cluster_dict["7.85300016403_arrows"], "Arrows_7.85300016403")
cmd.set("transparency", 0.2,"Features_7.85300016403")
cmd.group("Pharmacophore_7.85300016403", members="Features_7.85300016403")
cmd.group("Pharmacophore_7.85300016403", members="Arrows_7.85300016403")

if dirpath:
    f = join(dirpath, "30/label_threshold_7.85300016403.mol2")
else:
    f = "30/label_threshold_7.85300016403.mol2"

cmd.load(f, 'label_threshold_7.85300016403')
cmd.hide('everything', 'label_threshold_7.85300016403')
cmd.label("label_threshold_7.85300016403", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_7.85300016403', members= 'label_threshold_7.85300016403')


if dirpath:
    f = join(dirpath, '30/mesh.grd')
else:
    f = '30/mesh.grd'
cmd.load(f, 'mesh_30')
cmd.isomesh("isomesh_30", "mesh_30", 0.9)
cmd.color("grey80", "isomesh_30")
cmd.set('transparency', 0.4, "isomesh_30")

cmd.group('hotspot_30', "isomesh_30")
cmd.group('hotspot_30', "mesh_30")

if dirpath:
    f = join(dirpath, "31/label_threshold_4.9.mol2")
else:
    f = "31/label_threshold_4.9.mol2"

cmd.load(f, 'label_threshold_4.9')
cmd.hide('everything', 'label_threshold_4.9')
cmd.label("label_threshold_4.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [4.9]
gfiles = ['31/donor.grd', '31/apolar.grd', '31/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 31
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


cluster_dict = {"7.7610001564":[], "7.7610001564_arrows":[]}

cluster_dict["7.7610001564"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(7.5), float(0.5), float(-0.5), float(1.0)]

cluster_dict["7.7610001564_arrows"] += cgo_arrow([7.5,0.5,-0.5], [7.832,3.488,-1.705], color="blue red", name="Arrows_7.7610001564_1")

cluster_dict["7.7610001564"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(8.20319312234), float(0.976174201195), float(0.633397914245), float(1.0)]


cmd.load_cgo(cluster_dict["7.7610001564"], "Features_7.7610001564", 1)
cmd.load_cgo(cluster_dict["7.7610001564_arrows"], "Arrows_7.7610001564")
cmd.set("transparency", 0.2,"Features_7.7610001564")
cmd.group("Pharmacophore_7.7610001564", members="Features_7.7610001564")
cmd.group("Pharmacophore_7.7610001564", members="Arrows_7.7610001564")

if dirpath:
    f = join(dirpath, "31/label_threshold_7.7610001564.mol2")
else:
    f = "31/label_threshold_7.7610001564.mol2"

cmd.load(f, 'label_threshold_7.7610001564')
cmd.hide('everything', 'label_threshold_7.7610001564')
cmd.label("label_threshold_7.7610001564", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_7.7610001564', members= 'label_threshold_7.7610001564')


if dirpath:
    f = join(dirpath, '31/mesh.grd')
else:
    f = '31/mesh.grd'
cmd.load(f, 'mesh_31')
cmd.isomesh("isomesh_31", "mesh_31", 0.9)
cmd.color("grey80", "isomesh_31")
cmd.set('transparency', 0.4, "isomesh_31")

cmd.group('hotspot_31', "isomesh_31")
cmd.group('hotspot_31', "mesh_31")

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
