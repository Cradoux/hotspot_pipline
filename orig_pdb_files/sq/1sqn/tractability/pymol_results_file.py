
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
    f = join(dirpath, "0/label_threshold_13.1.mol2")
else:
    f = "0/label_threshold_13.1.mol2"

cmd.load(f, 'label_threshold_13.1')
cmd.hide('everything', 'label_threshold_13.1')
cmd.label("label_threshold_13.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.1]
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


cluster_dict = {"16.5699996948":[], "16.5699996948_arrows":[]}

cluster_dict["16.5699996948"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-4.5), float(5.5), float(21.0), float(1.0)]

cluster_dict["16.5699996948_arrows"] += cgo_arrow([-4.5,5.5,21.0], [-5.708,5.87,23.416], color="blue red", name="Arrows_16.5699996948_1")

cluster_dict["16.5699996948"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-2.5), float(-0.5), float(30.5), float(1.0)]

cluster_dict["16.5699996948_arrows"] += cgo_arrow([-2.5,-0.5,30.5], [-0.137,-1.469,31.664], color="blue red", name="Arrows_16.5699996948_2")

cluster_dict["16.5699996948"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-2.6484389603), float(2.00780458088), float(23.7618773147), float(1.0)]


cluster_dict["16.5699996948"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-5.0), float(-1.5), float(30.0), float(1.0)]

cluster_dict["16.5699996948_arrows"] += cgo_arrow([-5.0,-1.5,30.0], [-4.293,-3.037,32.774], color="red blue", name="Arrows_16.5699996948_3")

cluster_dict["16.5699996948"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-4.5), float(0.5), float(24.5), float(1.0)]


cluster_dict["16.5699996948"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-2.0), float(5.5), float(21.5), float(1.0)]

cluster_dict["16.5699996948_arrows"] += cgo_arrow([-2.0,5.5,21.5], [0.09,7.414,22.773], color="red blue", name="Arrows_16.5699996948_4")

cluster_dict["16.5699996948"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(0.0), float(0.5), float(23.5), float(1.0)]


cmd.load_cgo(cluster_dict["16.5699996948"], "Features_16.5699996948", 1)
cmd.load_cgo(cluster_dict["16.5699996948_arrows"], "Arrows_16.5699996948")
cmd.set("transparency", 0.2,"Features_16.5699996948")
cmd.group("Pharmacophore_16.5699996948", members="Features_16.5699996948")
cmd.group("Pharmacophore_16.5699996948", members="Arrows_16.5699996948")

if dirpath:
    f = join(dirpath, "0/label_threshold_16.5699996948.mol2")
else:
    f = "0/label_threshold_16.5699996948.mol2"

cmd.load(f, 'label_threshold_16.5699996948')
cmd.hide('everything', 'label_threshold_16.5699996948')
cmd.label("label_threshold_16.5699996948", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.5699996948', members= 'label_threshold_16.5699996948')


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


cluster_dict = {"16.5699996948":[], "16.5699996948_arrows":[]}

cluster_dict["16.5699996948"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-4.5), float(5.5), float(21.0), float(1.0)]

cluster_dict["16.5699996948_arrows"] += cgo_arrow([-4.5,5.5,21.0], [-5.708,5.87,23.416], color="blue red", name="Arrows_16.5699996948_1")

cluster_dict["16.5699996948"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-2.5), float(-0.5), float(30.5), float(1.0)]

cluster_dict["16.5699996948_arrows"] += cgo_arrow([-2.5,-0.5,30.5], [-0.137,-1.469,31.664], color="blue red", name="Arrows_16.5699996948_2")

cluster_dict["16.5699996948"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-2.64865964506), float(2.00895051665), float(23.7613990629), float(1.0)]


cluster_dict["16.5699996948"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-5.0), float(-1.5), float(30.0), float(1.0)]

cluster_dict["16.5699996948_arrows"] += cgo_arrow([-5.0,-1.5,30.0], [-4.293,-3.037,32.774], color="red blue", name="Arrows_16.5699996948_3")

cluster_dict["16.5699996948"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-4.5), float(0.5), float(24.5), float(1.0)]


cluster_dict["16.5699996948"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-2.0), float(5.5), float(21.5), float(1.0)]

cluster_dict["16.5699996948_arrows"] += cgo_arrow([-2.0,5.5,21.5], [0.09,7.414,22.773], color="red blue", name="Arrows_16.5699996948_4")

cluster_dict["16.5699996948"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(0.0), float(0.5), float(23.5), float(1.0)]


cmd.load_cgo(cluster_dict["16.5699996948"], "Features_16.5699996948", 1)
cmd.load_cgo(cluster_dict["16.5699996948_arrows"], "Arrows_16.5699996948")
cmd.set("transparency", 0.2,"Features_16.5699996948")
cmd.group("Pharmacophore_16.5699996948", members="Features_16.5699996948")
cmd.group("Pharmacophore_16.5699996948", members="Arrows_16.5699996948")

if dirpath:
    f = join(dirpath, "1/label_threshold_16.5699996948.mol2")
else:
    f = "1/label_threshold_16.5699996948.mol2"

cmd.load(f, 'label_threshold_16.5699996948')
cmd.hide('everything', 'label_threshold_16.5699996948')
cmd.label("label_threshold_16.5699996948", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.5699996948', members= 'label_threshold_16.5699996948')


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
    f = join(dirpath, "2/label_threshold_8.6.mol2")
else:
    f = "2/label_threshold_8.6.mol2"

cmd.load(f, 'label_threshold_8.6')
cmd.hide('everything', 'label_threshold_8.6')
cmd.label("label_threshold_8.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [8.6]
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


cluster_dict = {"13.5745000839":[], "13.5745000839_arrows":[]}

cluster_dict["13.5745000839"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-3.5), float(-3.0), float(41.0), float(1.0)]

cluster_dict["13.5745000839_arrows"] += cgo_arrow([-3.5,-3.0,41.0], [-0.462,-2.926,41.812], color="blue red", name="Arrows_13.5745000839_1")

cluster_dict["13.5745000839"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-3.0329623119), float(-1.2383667132), float(38.1405457727), float(1.0)]


cluster_dict["13.5745000839"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-7.0), float(-3.5), float(37.5), float(1.0)]

cluster_dict["13.5745000839_arrows"] += cgo_arrow([-7.0,-3.5,37.5], [-9.781,-0.858,37.828], color="red blue", name="Arrows_13.5745000839_2")

cluster_dict["13.5745000839"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-3.5), float(-6.0), float(38.0), float(1.0)]

cluster_dict["13.5745000839_arrows"] += cgo_arrow([-3.5,-6.0,38.0], [-2.053,-8.109,35.188], color="red blue", name="Arrows_13.5745000839_3")

cluster_dict["13.5745000839"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-3.5), float(2.0), float(38.5), float(1.0)]

cluster_dict["13.5745000839_arrows"] += cgo_arrow([-3.5,2.0,38.5], [-4.852,3.896,39.813], color="red blue", name="Arrows_13.5745000839_4")

cluster_dict["13.5745000839"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-2.0), float(-1.0), float(34.5), float(1.0)]

cluster_dict["13.5745000839_arrows"] += cgo_arrow([-2.0,-1.0,34.5], [-3.664,-3.099,34.971], color="red blue", name="Arrows_13.5745000839_5")

cluster_dict["13.5745000839"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-1.0), float(2.5), float(34.0), float(1.0)]

cluster_dict["13.5745000839_arrows"] += cgo_arrow([-1.0,2.5,34.0], [1.78,1.33,32.478], color="red blue", name="Arrows_13.5745000839_6")

cmd.load_cgo(cluster_dict["13.5745000839"], "Features_13.5745000839", 1)
cmd.load_cgo(cluster_dict["13.5745000839_arrows"], "Arrows_13.5745000839")
cmd.set("transparency", 0.2,"Features_13.5745000839")
cmd.group("Pharmacophore_13.5745000839", members="Features_13.5745000839")
cmd.group("Pharmacophore_13.5745000839", members="Arrows_13.5745000839")

if dirpath:
    f = join(dirpath, "2/label_threshold_13.5745000839.mol2")
else:
    f = "2/label_threshold_13.5745000839.mol2"

cmd.load(f, 'label_threshold_13.5745000839')
cmd.hide('everything', 'label_threshold_13.5745000839')
cmd.label("label_threshold_13.5745000839", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.5745000839', members= 'label_threshold_13.5745000839')


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
    f = join(dirpath, "3/label_threshold_5.4.mol2")
else:
    f = "3/label_threshold_5.4.mol2"

cmd.load(f, 'label_threshold_5.4')
cmd.hide('everything', 'label_threshold_5.4')
cmd.label("label_threshold_5.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [5.4]
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


cluster_dict = {"13.3269996643":[], "13.3269996643_arrows":[]}

cluster_dict["13.3269996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(11.5), float(4.0), float(28.0), float(1.0)]

cluster_dict["13.3269996643_arrows"] += cgo_arrow([11.5,4.0,28.0], [13.264,5.673,31.024], color="blue red", name="Arrows_13.3269996643_1")

cluster_dict["13.3269996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(7.5), float(24.5), float(1.0)]

cluster_dict["13.3269996643_arrows"] += cgo_arrow([15.0,7.5,24.5], [12.681,8.888,23.001], color="blue red", name="Arrows_13.3269996643_2")

cluster_dict["13.3269996643"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(13.5648642643), float(4.65252748597), float(25.5258575102), float(1.0)]


cluster_dict["13.3269996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(11.0), float(5.5), float(23.5), float(1.0)]

cluster_dict["13.3269996643_arrows"] += cgo_arrow([11.0,5.5,23.5], [8.063,5.715,24.743], color="red blue", name="Arrows_13.3269996643_3")

cluster_dict["13.3269996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(11.0), float(3.0), float(27.0), float(1.0)]

cluster_dict["13.3269996643_arrows"] += cgo_arrow([11.0,3.0,27.0], [13.879,2.269,29.166], color="red blue", name="Arrows_13.3269996643_4")

cluster_dict["13.3269996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(11.0), float(5.5), float(29.0), float(1.0)]

cluster_dict["13.3269996643_arrows"] += cgo_arrow([11.0,5.5,29.0], [13.264,5.673,31.024], color="red blue", name="Arrows_13.3269996643_5")

cluster_dict["13.3269996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(11.5), float(8.0), float(27.5), float(1.0)]

cluster_dict["13.3269996643_arrows"] += cgo_arrow([11.5,8.0,27.5], [13.264,5.673,31.024], color="red blue", name="Arrows_13.3269996643_6")

cluster_dict["13.3269996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(11.5), float(8.0), float(27.5), float(1.0)]

cluster_dict["13.3269996643_arrows"] += cgo_arrow([11.5,8.0,27.5], [13.264,5.673,31.024], color="red blue", name="Arrows_13.3269996643_7")

cluster_dict["13.3269996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(14.5), float(6.0), float(28.0), float(1.0)]

cluster_dict["13.3269996643_arrows"] += cgo_arrow([14.5,6.0,28.0], [13.264,5.673,31.024], color="red blue", name="Arrows_13.3269996643_8")

cmd.load_cgo(cluster_dict["13.3269996643"], "Features_13.3269996643", 1)
cmd.load_cgo(cluster_dict["13.3269996643_arrows"], "Arrows_13.3269996643")
cmd.set("transparency", 0.2,"Features_13.3269996643")
cmd.group("Pharmacophore_13.3269996643", members="Features_13.3269996643")
cmd.group("Pharmacophore_13.3269996643", members="Arrows_13.3269996643")

if dirpath:
    f = join(dirpath, "3/label_threshold_13.3269996643.mol2")
else:
    f = "3/label_threshold_13.3269996643.mol2"

cmd.load(f, 'label_threshold_13.3269996643')
cmd.hide('everything', 'label_threshold_13.3269996643')
cmd.label("label_threshold_13.3269996643", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.3269996643', members= 'label_threshold_13.3269996643')


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


cluster_dict = {"12.6520004272":[], "12.6520004272_arrows":[]}

cluster_dict["12.6520004272"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(10.7123759638), float(-8.97957001586), float(26.7726432336), float(1.0)]


cluster_dict["12.6520004272"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(9.0), float(-10.0), float(29.5), float(1.0)]

cluster_dict["12.6520004272_arrows"] += cgo_arrow([9.0,-10.0,29.5], [5.183,-9.495,27.864], color="red blue", name="Arrows_12.6520004272_1")

cluster_dict["12.6520004272"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(12.5), float(-10.0), float(29.0), float(1.0)]

cluster_dict["12.6520004272_arrows"] += cgo_arrow([12.5,-10.0,29.0], [15.215,-8.717,27.965], color="red blue", name="Arrows_12.6520004272_2")

cmd.load_cgo(cluster_dict["12.6520004272"], "Features_12.6520004272", 1)
cmd.load_cgo(cluster_dict["12.6520004272_arrows"], "Arrows_12.6520004272")
cmd.set("transparency", 0.2,"Features_12.6520004272")
cmd.group("Pharmacophore_12.6520004272", members="Features_12.6520004272")
cmd.group("Pharmacophore_12.6520004272", members="Arrows_12.6520004272")

if dirpath:
    f = join(dirpath, "4/label_threshold_12.6520004272.mol2")
else:
    f = "4/label_threshold_12.6520004272.mol2"

cmd.load(f, 'label_threshold_12.6520004272')
cmd.hide('everything', 'label_threshold_12.6520004272')
cmd.label("label_threshold_12.6520004272", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.6520004272', members= 'label_threshold_12.6520004272')


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
    f = join(dirpath, "5/label_threshold_7.6.mol2")
else:
    f = "5/label_threshold_7.6.mol2"

cmd.load(f, 'label_threshold_7.6')
cmd.hide('everything', 'label_threshold_7.6')
cmd.label("label_threshold_7.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [7.6]
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


cluster_dict = {"12.2030000687":[], "12.2030000687_arrows":[]}

cluster_dict["12.2030000687"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-4.5), float(1.5), float(12.5), float(1.0)]

cluster_dict["12.2030000687_arrows"] += cgo_arrow([-4.5,1.5,12.5], [-2.416,-1.047,12.988], color="blue red", name="Arrows_12.2030000687_1")

cluster_dict["12.2030000687"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-2.5), float(4.5), float(14.5), float(1.0)]

cluster_dict["12.2030000687_arrows"] += cgo_arrow([-2.5,4.5,14.5], [-2.48,5.061,16.995], color="blue red", name="Arrows_12.2030000687_2")

cluster_dict["12.2030000687"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-11.7486677302), float(-4.21397701151), float(16.5153391577), float(1.0)]


cluster_dict["12.2030000687"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-3.23477800195), float(3.01497337893), float(12.7865791587), float(1.0)]


cluster_dict["12.2030000687"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-3.20978557745), float(2.93770750819), float(19.922296564), float(1.0)]


cluster_dict["12.2030000687"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-7.0), float(1.0), float(12.5), float(1.0)]

cluster_dict["12.2030000687_arrows"] += cgo_arrow([-7.0,1.0,12.5], [-8.974,0.584,14.486], color="red blue", name="Arrows_12.2030000687_3")

cluster_dict["12.2030000687"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-3.5), float(4.5), float(14.5), float(1.0)]

cluster_dict["12.2030000687_arrows"] += cgo_arrow([-3.5,4.5,14.5], [-2.48,5.061,16.995], color="red blue", name="Arrows_12.2030000687_4")

cluster_dict["12.2030000687"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-3.0), float(1.5), float(12.5), float(1.0)]

cluster_dict["12.2030000687_arrows"] += cgo_arrow([-3.0,1.5,12.5], [-2.416,-1.047,12.988], color="red blue", name="Arrows_12.2030000687_5")

cmd.load_cgo(cluster_dict["12.2030000687"], "Features_12.2030000687", 1)
cmd.load_cgo(cluster_dict["12.2030000687_arrows"], "Arrows_12.2030000687")
cmd.set("transparency", 0.2,"Features_12.2030000687")
cmd.group("Pharmacophore_12.2030000687", members="Features_12.2030000687")
cmd.group("Pharmacophore_12.2030000687", members="Arrows_12.2030000687")

if dirpath:
    f = join(dirpath, "5/label_threshold_12.2030000687.mol2")
else:
    f = "5/label_threshold_12.2030000687.mol2"

cmd.load(f, 'label_threshold_12.2030000687')
cmd.hide('everything', 'label_threshold_12.2030000687')
cmd.label("label_threshold_12.2030000687", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.2030000687', members= 'label_threshold_12.2030000687')


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
    f = join(dirpath, "6/label_threshold_3.8.mol2")
else:
    f = "6/label_threshold_3.8.mol2"

cmd.load(f, 'label_threshold_3.8')
cmd.hide('everything', 'label_threshold_3.8')
cmd.label("label_threshold_3.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [3.8]
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


cluster_dict = {"10.4420003891":[], "10.4420003891_arrows":[]}

cluster_dict["10.4420003891"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-4.5), float(1.5), float(12.5), float(1.0)]

cluster_dict["10.4420003891_arrows"] += cgo_arrow([-4.5,1.5,12.5], [-2.416,-1.047,12.988], color="blue red", name="Arrows_10.4420003891_1")

cluster_dict["10.4420003891"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-4.0), float(4.0), float(15.0), float(1.0)]

cluster_dict["10.4420003891_arrows"] += cgo_arrow([-4.0,4.0,15.0], [-2.48,5.061,16.995], color="blue red", name="Arrows_10.4420003891_2")

cluster_dict["10.4420003891"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-12.2241482961), float(-4.67262333356), float(16.6491217589), float(1.0)]


cluster_dict["10.4420003891"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-5.39100189804), float(2.02412202582), float(12.9223714024), float(1.0)]


cluster_dict["10.4420003891"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-4.12152125861), float(2.09109794182), float(20.1295348815), float(1.0)]


cluster_dict["10.4420003891"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-7.0), float(1.0), float(12.5), float(1.0)]

cluster_dict["10.4420003891_arrows"] += cgo_arrow([-7.0,1.0,12.5], [-8.974,0.584,14.486], color="red blue", name="Arrows_10.4420003891_3")

cluster_dict["10.4420003891"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-4.0), float(4.0), float(14.5), float(1.0)]

cluster_dict["10.4420003891_arrows"] += cgo_arrow([-4.0,4.0,14.5], [-2.48,5.061,16.995], color="red blue", name="Arrows_10.4420003891_4")

cluster_dict["10.4420003891"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-3.0), float(1.5), float(12.5), float(1.0)]

cluster_dict["10.4420003891_arrows"] += cgo_arrow([-3.0,1.5,12.5], [-2.416,-1.047,12.988], color="red blue", name="Arrows_10.4420003891_5")

cmd.load_cgo(cluster_dict["10.4420003891"], "Features_10.4420003891", 1)
cmd.load_cgo(cluster_dict["10.4420003891_arrows"], "Arrows_10.4420003891")
cmd.set("transparency", 0.2,"Features_10.4420003891")
cmd.group("Pharmacophore_10.4420003891", members="Features_10.4420003891")
cmd.group("Pharmacophore_10.4420003891", members="Arrows_10.4420003891")

if dirpath:
    f = join(dirpath, "6/label_threshold_10.4420003891.mol2")
else:
    f = "6/label_threshold_10.4420003891.mol2"

cmd.load(f, 'label_threshold_10.4420003891')
cmd.hide('everything', 'label_threshold_10.4420003891')
cmd.label("label_threshold_10.4420003891", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.4420003891', members= 'label_threshold_10.4420003891')


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

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
