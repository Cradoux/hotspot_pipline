
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
    f = join(dirpath, "0/label_threshold_15.8.mol2")
else:
    f = "0/label_threshold_15.8.mol2"

cmd.load(f, 'label_threshold_15.8')
cmd.hide('everything', 'label_threshold_15.8')
cmd.label("label_threshold_15.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.8]
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


cluster_dict = {"17.4930000305":[], "17.4930000305_arrows":[]}

cluster_dict["17.4930000305"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-52.5), float(-17.0), float(-9.5), float(1.0)]

cluster_dict["17.4930000305_arrows"] += cgo_arrow([-52.5,-17.0,-9.5], [-54.917,-17.538,-10.051], color="blue red", name="Arrows_17.4930000305_1")

cluster_dict["17.4930000305"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-49.5), float(-13.0), float(-6.5), float(1.0)]

cluster_dict["17.4930000305_arrows"] += cgo_arrow([-49.5,-13.0,-6.5], [-47.072,-14.06,-4.739], color="blue red", name="Arrows_17.4930000305_2")

cluster_dict["17.4930000305"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-49.0), float(-16.5), float(-6.5), float(1.0)]

cluster_dict["17.4930000305_arrows"] += cgo_arrow([-49.0,-16.5,-6.5], [-47.072,-14.06,-4.739], color="blue red", name="Arrows_17.4930000305_3")

cluster_dict["17.4930000305"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-42.0), float(-14.5), float(-8.5), float(1.0)]

cluster_dict["17.4930000305_arrows"] += cgo_arrow([-42.0,-14.5,-8.5], [-43.537,-16.731,-9.097], color="blue red", name="Arrows_17.4930000305_4")

cluster_dict["17.4930000305"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-46.2693407668), float(-12.5804843478), float(-8.89678777135), float(1.0)]


cluster_dict["17.4930000305"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-33.3599320875), float(-10.4008982587), float(-15.23607254), float(1.0)]


cluster_dict["17.4930000305"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-52.5), float(-12.5), float(-8.0), float(1.0)]

cluster_dict["17.4930000305_arrows"] += cgo_arrow([-52.5,-12.5,-8.0], [-52.889,-10.6,-6.798], color="red blue", name="Arrows_17.4930000305_5")

cluster_dict["17.4930000305"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-47.5), float(-17.0), float(-5.0), float(1.0)]

cluster_dict["17.4930000305_arrows"] += cgo_arrow([-47.5,-17.0,-5.0], [-45.547,-18.463,-3.903], color="red blue", name="Arrows_17.4930000305_6")

cluster_dict["17.4930000305"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-49.5), float(-10.0), float(-6.5), float(1.0)]

cluster_dict["17.4930000305_arrows"] += cgo_arrow([-49.5,-10.0,-6.5], [-52.889,-10.6,-6.798], color="red blue", name="Arrows_17.4930000305_7")

cluster_dict["17.4930000305"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-48.5), float(-12.0), float(-6.5), float(1.0)]

cluster_dict["17.4930000305_arrows"] += cgo_arrow([-48.5,-12.0,-6.5], [-47.072,-14.06,-4.739], color="red blue", name="Arrows_17.4930000305_8")

cluster_dict["17.4930000305"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-46.0), float(-15.5), float(-8.5), float(1.0)]

cluster_dict["17.4930000305_arrows"] += cgo_arrow([-46.0,-15.5,-8.5], [-43.537,-16.731,-9.097], color="red blue", name="Arrows_17.4930000305_9")

cmd.load_cgo(cluster_dict["17.4930000305"], "Features_17.4930000305", 1)
cmd.load_cgo(cluster_dict["17.4930000305_arrows"], "Arrows_17.4930000305")
cmd.set("transparency", 0.2,"Features_17.4930000305")
cmd.group("Pharmacophore_17.4930000305", members="Features_17.4930000305")
cmd.group("Pharmacophore_17.4930000305", members="Arrows_17.4930000305")

if dirpath:
    f = join(dirpath, "0/label_threshold_17.4930000305.mol2")
else:
    f = "0/label_threshold_17.4930000305.mol2"

cmd.load(f, 'label_threshold_17.4930000305')
cmd.hide('everything', 'label_threshold_17.4930000305')
cmd.label("label_threshold_17.4930000305", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.4930000305', members= 'label_threshold_17.4930000305')


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
    f = join(dirpath, "1/label_threshold_1.0.mol2")
else:
    f = "1/label_threshold_1.0.mol2"

cmd.load(f, 'label_threshold_1.0')
cmd.hide('everything', 'label_threshold_1.0')
cmd.label("label_threshold_1.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [1.0]
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


cluster_dict = {"17.1779994965":[], "17.1779994965_arrows":[]}

cluster_dict["17.1779994965"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-33.0), float(-10.5), float(-12.5), float(1.0)]

cluster_dict["17.1779994965_arrows"] += cgo_arrow([-33.0,-10.5,-12.5], [-33.791,-12.98,-13.233], color="blue red", name="Arrows_17.1779994965_1")

cluster_dict["17.1779994965"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-33.3934795835), float(-10.5980241291), float(-15.2376239091), float(1.0)]


cluster_dict["17.1779994965"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-36.2416347018), float(-8.82770136057), float(-13.1708107385), float(1.0)]


cmd.load_cgo(cluster_dict["17.1779994965"], "Features_17.1779994965", 1)
cmd.load_cgo(cluster_dict["17.1779994965_arrows"], "Arrows_17.1779994965")
cmd.set("transparency", 0.2,"Features_17.1779994965")
cmd.group("Pharmacophore_17.1779994965", members="Features_17.1779994965")
cmd.group("Pharmacophore_17.1779994965", members="Arrows_17.1779994965")

if dirpath:
    f = join(dirpath, "1/label_threshold_17.1779994965.mol2")
else:
    f = "1/label_threshold_17.1779994965.mol2"

cmd.load(f, 'label_threshold_17.1779994965')
cmd.hide('everything', 'label_threshold_17.1779994965')
cmd.label("label_threshold_17.1779994965", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.1779994965', members= 'label_threshold_17.1779994965')


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


cluster_dict = {"16.9039993286":[], "16.9039993286_arrows":[]}

cluster_dict["16.9039993286"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-33.0), float(-47.5), float(-7.0), float(1.0)]

cluster_dict["16.9039993286_arrows"] += cgo_arrow([-33.0,-47.5,-7.0], [-30.7,-45.644,-5.488], color="blue red", name="Arrows_16.9039993286_1")

cluster_dict["16.9039993286"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-37.4006126458), float(-45.4407994057), float(-6.84243676431), float(1.0)]


cluster_dict["16.9039993286"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-39.5), float(-46.0), float(-7.5), float(1.0)]

cluster_dict["16.9039993286_arrows"] += cgo_arrow([-39.5,-46.0,-7.5], [-41.624,-46.288,-5.26], color="red blue", name="Arrows_16.9039993286_2")

cluster_dict["16.9039993286"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(-42.0), float(-5.0), float(1.0)]

cluster_dict["16.9039993286_arrows"] += cgo_arrow([-37.0,-42.0,-5.0], [-38.51,-41.271,-7.117], color="red blue", name="Arrows_16.9039993286_3")

cluster_dict["16.9039993286"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-35.5), float(-44.0), float(-2.5), float(1.0)]

cluster_dict["16.9039993286_arrows"] += cgo_arrow([-35.5,-44.0,-2.5], [-34.702,-46.507,-2.268], color="red blue", name="Arrows_16.9039993286_4")

cluster_dict["16.9039993286"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-35.0), float(-46.0), float(-8.5), float(1.0)]

cluster_dict["16.9039993286_arrows"] += cgo_arrow([-35.0,-46.0,-8.5], [-34.43,-44.132,-11.02], color="red blue", name="Arrows_16.9039993286_5")

cluster_dict["16.9039993286"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-32.5), float(-50.0), float(-7.0), float(1.0)]

cluster_dict["16.9039993286_arrows"] += cgo_arrow([-32.5,-50.0,-7.0], [-33.632,-52.301,-6.73], color="red blue", name="Arrows_16.9039993286_6")

cmd.load_cgo(cluster_dict["16.9039993286"], "Features_16.9039993286", 1)
cmd.load_cgo(cluster_dict["16.9039993286_arrows"], "Arrows_16.9039993286")
cmd.set("transparency", 0.2,"Features_16.9039993286")
cmd.group("Pharmacophore_16.9039993286", members="Features_16.9039993286")
cmd.group("Pharmacophore_16.9039993286", members="Arrows_16.9039993286")

if dirpath:
    f = join(dirpath, "2/label_threshold_16.9039993286.mol2")
else:
    f = "2/label_threshold_16.9039993286.mol2"

cmd.load(f, 'label_threshold_16.9039993286')
cmd.hide('everything', 'label_threshold_16.9039993286')
cmd.label("label_threshold_16.9039993286", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.9039993286', members= 'label_threshold_16.9039993286')


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
    f = join(dirpath, "3/label_threshold_15.4.mol2")
else:
    f = "3/label_threshold_15.4.mol2"

cmd.load(f, 'label_threshold_15.4')
cmd.hide('everything', 'label_threshold_15.4')
cmd.label("label_threshold_15.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.4]
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


cluster_dict = {"16.8540000916":[], "16.8540000916_arrows":[]}

cluster_dict["16.8540000916"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-50.0), float(-43.0), float(-4.0), float(1.0)]

cluster_dict["16.8540000916_arrows"] += cgo_arrow([-50.0,-43.0,-4.0], [-50.188,-40.546,-5.827], color="blue red", name="Arrows_16.8540000916_1")

cluster_dict["16.8540000916"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-45.5), float(-40.0), float(4.5), float(1.0)]

cluster_dict["16.8540000916_arrows"] += cgo_arrow([-45.5,-40.0,4.5], [-43.549,-37.887,3.576], color="blue red", name="Arrows_16.8540000916_2")

cluster_dict["16.8540000916"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-44.5), float(-38.0), float(8.5), float(1.0)]

cluster_dict["16.8540000916_arrows"] += cgo_arrow([-44.5,-38.0,8.5], [-42.285,-36.612,9.774], color="blue red", name="Arrows_16.8540000916_3")

cluster_dict["16.8540000916"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-41.5), float(-45.0), float(6.5), float(1.0)]

cluster_dict["16.8540000916_arrows"] += cgo_arrow([-41.5,-45.0,6.5], [-40.482,-47.349,4.858], color="blue red", name="Arrows_16.8540000916_4")

cluster_dict["16.8540000916"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-47.1429311567), float(-41.4588434719), float(1.84556256702), float(1.0)]


cluster_dict["16.8540000916"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-43.3792214427), float(-42.1909233738), float(5.66313937497), float(1.0)]


cluster_dict["16.8540000916"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-37.6771839458), float(-43.2811010426), float(-5.27550980628), float(1.0)]


cluster_dict["16.8540000916"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-45.5), float(-42.0), float(-0.5), float(1.0)]

cluster_dict["16.8540000916_arrows"] += cgo_arrow([-45.5,-42.0,-0.5], [-45.326,-45.543,-3.379], color="red blue", name="Arrows_16.8540000916_5")

cluster_dict["16.8540000916"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-44.5), float(-39.5), float(10.0), float(1.0)]

cluster_dict["16.8540000916_arrows"] += cgo_arrow([-44.5,-39.5,10.0], [-42.638,-37.974,11.517], color="red blue", name="Arrows_16.8540000916_6")

cluster_dict["16.8540000916"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-39.5), float(-46.0), float(-7.5), float(1.0)]

cluster_dict["16.8540000916_arrows"] += cgo_arrow([-39.5,-46.0,-7.5], [-41.624,-46.288,-5.26], color="red blue", name="Arrows_16.8540000916_7")

cluster_dict["16.8540000916"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(-42.0), float(-5.0), float(1.0)]

cluster_dict["16.8540000916_arrows"] += cgo_arrow([-37.0,-42.0,-5.0], [-38.51,-41.271,-7.117], color="red blue", name="Arrows_16.8540000916_8")

cmd.load_cgo(cluster_dict["16.8540000916"], "Features_16.8540000916", 1)
cmd.load_cgo(cluster_dict["16.8540000916_arrows"], "Arrows_16.8540000916")
cmd.set("transparency", 0.2,"Features_16.8540000916")
cmd.group("Pharmacophore_16.8540000916", members="Features_16.8540000916")
cmd.group("Pharmacophore_16.8540000916", members="Arrows_16.8540000916")

if dirpath:
    f = join(dirpath, "3/label_threshold_16.8540000916.mol2")
else:
    f = "3/label_threshold_16.8540000916.mol2"

cmd.load(f, 'label_threshold_16.8540000916')
cmd.hide('everything', 'label_threshold_16.8540000916')
cmd.label("label_threshold_16.8540000916", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.8540000916', members= 'label_threshold_16.8540000916')


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
    f = join(dirpath, "4/label_threshold_12.5.mol2")
else:
    f = "4/label_threshold_12.5.mol2"

cmd.load(f, 'label_threshold_12.5')
cmd.hide('everything', 'label_threshold_12.5')
cmd.label("label_threshold_12.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.5]
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


cluster_dict = {"16.7430000305":[], "16.7430000305_arrows":[]}

cluster_dict["16.7430000305"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-33.0), float(-47.5), float(-7.0), float(1.0)]

cluster_dict["16.7430000305_arrows"] += cgo_arrow([-33.0,-47.5,-7.0], [-30.7,-45.644,-5.488], color="blue red", name="Arrows_16.7430000305_1")

cluster_dict["16.7430000305"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-37.3773126847), float(-45.1165197701), float(-6.65529551727), float(1.0)]


cluster_dict["16.7430000305"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-39.5), float(-46.0), float(-7.5), float(1.0)]

cluster_dict["16.7430000305_arrows"] += cgo_arrow([-39.5,-46.0,-7.5], [-41.624,-46.288,-5.26], color="red blue", name="Arrows_16.7430000305_2")

cluster_dict["16.7430000305"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(-42.0), float(-5.0), float(1.0)]

cluster_dict["16.7430000305_arrows"] += cgo_arrow([-37.0,-42.0,-5.0], [-38.51,-41.271,-7.117], color="red blue", name="Arrows_16.7430000305_3")

cluster_dict["16.7430000305"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-35.5), float(-44.0), float(-2.5), float(1.0)]

cluster_dict["16.7430000305_arrows"] += cgo_arrow([-35.5,-44.0,-2.5], [-34.702,-46.507,-2.268], color="red blue", name="Arrows_16.7430000305_4")

cluster_dict["16.7430000305"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-35.0), float(-46.0), float(-8.5), float(1.0)]

cluster_dict["16.7430000305_arrows"] += cgo_arrow([-35.0,-46.0,-8.5], [-34.43,-44.132,-11.02], color="red blue", name="Arrows_16.7430000305_5")

cluster_dict["16.7430000305"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-32.5), float(-50.0), float(-7.0), float(1.0)]

cluster_dict["16.7430000305_arrows"] += cgo_arrow([-32.5,-50.0,-7.0], [-33.632,-52.301,-6.73], color="red blue", name="Arrows_16.7430000305_6")

cluster_dict["16.7430000305"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-31.5), float(-41.0), float(-5.0), float(1.0)]

cluster_dict["16.7430000305_arrows"] += cgo_arrow([-31.5,-41.0,-5.0], [-31.294,-43.503,-5.289], color="red blue", name="Arrows_16.7430000305_7")

cmd.load_cgo(cluster_dict["16.7430000305"], "Features_16.7430000305", 1)
cmd.load_cgo(cluster_dict["16.7430000305_arrows"], "Arrows_16.7430000305")
cmd.set("transparency", 0.2,"Features_16.7430000305")
cmd.group("Pharmacophore_16.7430000305", members="Features_16.7430000305")
cmd.group("Pharmacophore_16.7430000305", members="Arrows_16.7430000305")

if dirpath:
    f = join(dirpath, "4/label_threshold_16.7430000305.mol2")
else:
    f = "4/label_threshold_16.7430000305.mol2"

cmd.load(f, 'label_threshold_16.7430000305')
cmd.hide('everything', 'label_threshold_16.7430000305')
cmd.label("label_threshold_16.7430000305", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.7430000305', members= 'label_threshold_16.7430000305')


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
    f = join(dirpath, "5/label_threshold_15.6.mol2")
else:
    f = "5/label_threshold_15.6.mol2"

cmd.load(f, 'label_threshold_15.6')
cmd.hide('everything', 'label_threshold_15.6')
cmd.label("label_threshold_15.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.6]
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


cluster_dict = {"16.6130008698":[], "16.6130008698_arrows":[]}

cluster_dict["16.6130008698"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-52.5), float(-17.0), float(-9.0), float(1.0)]

cluster_dict["16.6130008698_arrows"] += cgo_arrow([-52.5,-17.0,-9.0], [-54.917,-17.538,-10.051], color="blue red", name="Arrows_16.6130008698_1")

cluster_dict["16.6130008698"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-49.5), float(-13.0), float(-6.5), float(1.0)]

cluster_dict["16.6130008698_arrows"] += cgo_arrow([-49.5,-13.0,-6.5], [-47.072,-14.06,-4.739], color="blue red", name="Arrows_16.6130008698_2")

cluster_dict["16.6130008698"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-49.0), float(-16.5), float(-6.5), float(1.0)]

cluster_dict["16.6130008698_arrows"] += cgo_arrow([-49.0,-16.5,-6.5], [-47.072,-14.06,-4.739], color="blue red", name="Arrows_16.6130008698_3")

cluster_dict["16.6130008698"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-42.0), float(-14.5), float(-8.5), float(1.0)]

cluster_dict["16.6130008698_arrows"] += cgo_arrow([-42.0,-14.5,-8.5], [-43.537,-16.731,-9.097], color="blue red", name="Arrows_16.6130008698_4")

cluster_dict["16.6130008698"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-46.0266906834), float(-12.1576438143), float(-8.77046589698), float(1.0)]


cluster_dict["16.6130008698"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-43.9935064935), float(0.233766233766), float(0.811688311688), float(1.0)]


cluster_dict["16.6130008698"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-34.4594749986), float(-10.234695195), float(-14.2218794251), float(1.0)]


cluster_dict["16.6130008698"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-52.5), float(-12.5), float(-8.0), float(1.0)]

cluster_dict["16.6130008698_arrows"] += cgo_arrow([-52.5,-12.5,-8.0], [-52.889,-10.6,-6.798], color="red blue", name="Arrows_16.6130008698_5")

cluster_dict["16.6130008698"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-47.5), float(-17.0), float(-5.0), float(1.0)]

cluster_dict["16.6130008698_arrows"] += cgo_arrow([-47.5,-17.0,-5.0], [-45.547,-18.463,-3.903], color="red blue", name="Arrows_16.6130008698_6")

cluster_dict["16.6130008698"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-49.5), float(-10.0), float(-6.5), float(1.0)]

cluster_dict["16.6130008698_arrows"] += cgo_arrow([-49.5,-10.0,-6.5], [-52.889,-10.6,-6.798], color="red blue", name="Arrows_16.6130008698_7")

cluster_dict["16.6130008698"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-48.5), float(-12.0), float(-6.5), float(1.0)]

cluster_dict["16.6130008698_arrows"] += cgo_arrow([-48.5,-12.0,-6.5], [-47.072,-14.06,-4.739], color="red blue", name="Arrows_16.6130008698_8")

cluster_dict["16.6130008698"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-46.0), float(-15.5), float(-8.5), float(1.0)]

cluster_dict["16.6130008698_arrows"] += cgo_arrow([-46.0,-15.5,-8.5], [-43.537,-16.731,-9.097], color="red blue", name="Arrows_16.6130008698_9")

cmd.load_cgo(cluster_dict["16.6130008698"], "Features_16.6130008698", 1)
cmd.load_cgo(cluster_dict["16.6130008698_arrows"], "Arrows_16.6130008698")
cmd.set("transparency", 0.2,"Features_16.6130008698")
cmd.group("Pharmacophore_16.6130008698", members="Features_16.6130008698")
cmd.group("Pharmacophore_16.6130008698", members="Arrows_16.6130008698")

if dirpath:
    f = join(dirpath, "5/label_threshold_16.6130008698.mol2")
else:
    f = "5/label_threshold_16.6130008698.mol2"

cmd.load(f, 'label_threshold_16.6130008698')
cmd.hide('everything', 'label_threshold_16.6130008698')
cmd.label("label_threshold_16.6130008698", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.6130008698', members= 'label_threshold_16.6130008698')


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
    f = join(dirpath, "6/label_threshold_15.2.mol2")
else:
    f = "6/label_threshold_15.2.mol2"

cmd.load(f, 'label_threshold_15.2')
cmd.hide('everything', 'label_threshold_15.2')
cmd.label("label_threshold_15.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.2]
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


cluster_dict = {"16.3530006409":[], "16.3530006409_arrows":[]}

cluster_dict["16.3530006409"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-52.5), float(-17.0), float(-9.0), float(1.0)]

cluster_dict["16.3530006409_arrows"] += cgo_arrow([-52.5,-17.0,-9.0], [-54.917,-17.538,-10.051], color="blue red", name="Arrows_16.3530006409_1")

cluster_dict["16.3530006409"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-52.0), float(-15.0), float(-3.5), float(1.0)]

cluster_dict["16.3530006409_arrows"] += cgo_arrow([-52.0,-15.0,-3.5], [-49.222,-15.773,-2.99], color="blue red", name="Arrows_16.3530006409_2")

cluster_dict["16.3530006409"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-49.5), float(-13.0), float(-6.5), float(1.0)]

cluster_dict["16.3530006409_arrows"] += cgo_arrow([-49.5,-13.0,-6.5], [-47.072,-14.06,-4.739], color="blue red", name="Arrows_16.3530006409_3")

cluster_dict["16.3530006409"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-49.0), float(-16.5), float(-6.5), float(1.0)]

cluster_dict["16.3530006409_arrows"] += cgo_arrow([-49.0,-16.5,-6.5], [-47.072,-14.06,-4.739], color="blue red", name="Arrows_16.3530006409_4")

cluster_dict["16.3530006409"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-42.0), float(-14.5), float(-8.5), float(1.0)]

cluster_dict["16.3530006409_arrows"] += cgo_arrow([-42.0,-14.5,-8.5], [-43.537,-16.731,-9.097], color="blue red", name="Arrows_16.3530006409_5")

cluster_dict["16.3530006409"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-45.7153583002), float(-12.0019842995), float(-8.63909376135), float(1.0)]


cluster_dict["16.3530006409"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-47.5), float(-17.0), float(-5.0), float(1.0)]

cluster_dict["16.3530006409_arrows"] += cgo_arrow([-47.5,-17.0,-5.0], [-45.547,-18.463,-3.903], color="red blue", name="Arrows_16.3530006409_6")

cluster_dict["16.3530006409"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-49.5), float(-10.0), float(-6.5), float(1.0)]

cluster_dict["16.3530006409_arrows"] += cgo_arrow([-49.5,-10.0,-6.5], [-52.889,-10.6,-6.798], color="red blue", name="Arrows_16.3530006409_7")

cluster_dict["16.3530006409"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-48.5), float(-12.0), float(-6.5), float(1.0)]

cluster_dict["16.3530006409_arrows"] += cgo_arrow([-48.5,-12.0,-6.5], [-47.072,-14.06,-4.739], color="red blue", name="Arrows_16.3530006409_8")

cluster_dict["16.3530006409"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-46.0), float(-15.5), float(-8.5), float(1.0)]

cluster_dict["16.3530006409_arrows"] += cgo_arrow([-46.0,-15.5,-8.5], [-43.537,-16.731,-9.097], color="red blue", name="Arrows_16.3530006409_9")

cmd.load_cgo(cluster_dict["16.3530006409"], "Features_16.3530006409", 1)
cmd.load_cgo(cluster_dict["16.3530006409_arrows"], "Arrows_16.3530006409")
cmd.set("transparency", 0.2,"Features_16.3530006409")
cmd.group("Pharmacophore_16.3530006409", members="Features_16.3530006409")
cmd.group("Pharmacophore_16.3530006409", members="Arrows_16.3530006409")

if dirpath:
    f = join(dirpath, "6/label_threshold_16.3530006409.mol2")
else:
    f = "6/label_threshold_16.3530006409.mol2"

cmd.load(f, 'label_threshold_16.3530006409')
cmd.hide('everything', 'label_threshold_16.3530006409')
cmd.label("label_threshold_16.3530006409", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.3530006409', members= 'label_threshold_16.3530006409')


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
    f = join(dirpath, "7/label_threshold_14.3.mol2")
else:
    f = "7/label_threshold_14.3.mol2"

cmd.load(f, 'label_threshold_14.3')
cmd.hide('everything', 'label_threshold_14.3')
cmd.label("label_threshold_14.3", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [14.3]
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


cluster_dict = {"16.2220001221":[], "16.2220001221_arrows":[]}

cluster_dict["16.2220001221"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-50.0), float(-43.0), float(-4.0), float(1.0)]

cluster_dict["16.2220001221_arrows"] += cgo_arrow([-50.0,-43.0,-4.0], [-50.188,-40.546,-5.827], color="blue red", name="Arrows_16.2220001221_1")

cluster_dict["16.2220001221"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-45.5), float(-40.0), float(4.5), float(1.0)]

cluster_dict["16.2220001221_arrows"] += cgo_arrow([-45.5,-40.0,4.5], [-43.549,-37.887,3.576], color="blue red", name="Arrows_16.2220001221_2")

cluster_dict["16.2220001221"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-44.5), float(-38.0), float(8.5), float(1.0)]

cluster_dict["16.2220001221_arrows"] += cgo_arrow([-44.5,-38.0,8.5], [-42.285,-36.612,9.774], color="blue red", name="Arrows_16.2220001221_3")

cluster_dict["16.2220001221"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-41.5), float(-45.0), float(6.5), float(1.0)]

cluster_dict["16.2220001221_arrows"] += cgo_arrow([-41.5,-45.0,6.5], [-40.482,-47.349,4.858], color="blue red", name="Arrows_16.2220001221_4")

cluster_dict["16.2220001221"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-46.3872081411), float(-41.6800569683), float(3.00972239297), float(1.0)]


cluster_dict["16.2220001221"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-37.2049516512), float(-41.0949358825), float(-3.5566730164), float(1.0)]


cluster_dict["16.2220001221"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-48.5), float(-44.0), float(-5.5), float(1.0)]

cluster_dict["16.2220001221_arrows"] += cgo_arrow([-48.5,-44.0,-5.5], [-48.805,-41.413,-7.401], color="red blue", name="Arrows_16.2220001221_5")

cluster_dict["16.2220001221"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-45.5), float(-42.0), float(-0.5), float(1.0)]

cluster_dict["16.2220001221_arrows"] += cgo_arrow([-45.5,-42.0,-0.5], [-45.326,-45.543,-3.379], color="red blue", name="Arrows_16.2220001221_6")

cluster_dict["16.2220001221"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-44.5), float(-39.5), float(10.0), float(1.0)]

cluster_dict["16.2220001221_arrows"] += cgo_arrow([-44.5,-39.5,10.0], [-42.638,-37.974,11.517], color="red blue", name="Arrows_16.2220001221_7")

cluster_dict["16.2220001221"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(-42.0), float(-5.0), float(1.0)]

cluster_dict["16.2220001221_arrows"] += cgo_arrow([-37.0,-42.0,-5.0], [-38.51,-41.271,-7.117], color="red blue", name="Arrows_16.2220001221_8")

cluster_dict["16.2220001221"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-36.0), float(-44.0), float(-2.5), float(1.0)]

cluster_dict["16.2220001221_arrows"] += cgo_arrow([-36.0,-44.0,-2.5], [-34.702,-46.507,-2.268], color="red blue", name="Arrows_16.2220001221_9")

cmd.load_cgo(cluster_dict["16.2220001221"], "Features_16.2220001221", 1)
cmd.load_cgo(cluster_dict["16.2220001221_arrows"], "Arrows_16.2220001221")
cmd.set("transparency", 0.2,"Features_16.2220001221")
cmd.group("Pharmacophore_16.2220001221", members="Features_16.2220001221")
cmd.group("Pharmacophore_16.2220001221", members="Arrows_16.2220001221")

if dirpath:
    f = join(dirpath, "7/label_threshold_16.2220001221.mol2")
else:
    f = "7/label_threshold_16.2220001221.mol2"

cmd.load(f, 'label_threshold_16.2220001221')
cmd.hide('everything', 'label_threshold_16.2220001221')
cmd.label("label_threshold_16.2220001221", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.2220001221', members= 'label_threshold_16.2220001221')


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
    f = join(dirpath, "8/label_threshold_10.3.mol2")
else:
    f = "8/label_threshold_10.3.mol2"

cmd.load(f, 'label_threshold_10.3')
cmd.hide('everything', 'label_threshold_10.3')
cmd.label("label_threshold_10.3", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.3]
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


cluster_dict = {"15.7869997025":[], "15.7869997025_arrows":[]}

cluster_dict["15.7869997025"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-33.0), float(-47.5), float(-7.0), float(1.0)]

cluster_dict["15.7869997025_arrows"] += cgo_arrow([-33.0,-47.5,-7.0], [-30.7,-45.644,-5.488], color="blue red", name="Arrows_15.7869997025_1")

cluster_dict["15.7869997025"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-36.6020292656), float(-44.6514947556), float(-6.28493028748), float(1.0)]


cluster_dict["15.7869997025"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-39.5), float(-46.0), float(-7.5), float(1.0)]

cluster_dict["15.7869997025_arrows"] += cgo_arrow([-39.5,-46.0,-7.5], [-41.624,-46.288,-5.26], color="red blue", name="Arrows_15.7869997025_2")

cluster_dict["15.7869997025"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-36.5), float(-43.0), float(-5.0), float(1.0)]

cluster_dict["15.7869997025_arrows"] += cgo_arrow([-36.5,-43.0,-5.0], [-34.956,-41.574,-6.87], color="red blue", name="Arrows_15.7869997025_3")

cluster_dict["15.7869997025"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-36.5), float(-43.0), float(-5.0), float(1.0)]

cluster_dict["15.7869997025_arrows"] += cgo_arrow([-36.5,-43.0,-5.0], [-34.956,-41.574,-6.87], color="red blue", name="Arrows_15.7869997025_4")

cluster_dict["15.7869997025"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-32.5), float(-50.0), float(-7.0), float(1.0)]

cluster_dict["15.7869997025_arrows"] += cgo_arrow([-32.5,-50.0,-7.0], [-33.632,-52.301,-6.73], color="red blue", name="Arrows_15.7869997025_5")

cluster_dict["15.7869997025"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-31.5), float(-41.0), float(-5.0), float(1.0)]

cluster_dict["15.7869997025_arrows"] += cgo_arrow([-31.5,-41.0,-5.0], [-31.294,-43.503,-5.289], color="red blue", name="Arrows_15.7869997025_6")

cmd.load_cgo(cluster_dict["15.7869997025"], "Features_15.7869997025", 1)
cmd.load_cgo(cluster_dict["15.7869997025_arrows"], "Arrows_15.7869997025")
cmd.set("transparency", 0.2,"Features_15.7869997025")
cmd.group("Pharmacophore_15.7869997025", members="Features_15.7869997025")
cmd.group("Pharmacophore_15.7869997025", members="Arrows_15.7869997025")

if dirpath:
    f = join(dirpath, "8/label_threshold_15.7869997025.mol2")
else:
    f = "8/label_threshold_15.7869997025.mol2"

cmd.load(f, 'label_threshold_15.7869997025')
cmd.hide('everything', 'label_threshold_15.7869997025')
cmd.label("label_threshold_15.7869997025", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.7869997025', members= 'label_threshold_15.7869997025')


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
    f = join(dirpath, "9/label_threshold_13.0.mol2")
else:
    f = "9/label_threshold_13.0.mol2"

cmd.load(f, 'label_threshold_13.0')
cmd.hide('everything', 'label_threshold_13.0')
cmd.label("label_threshold_13.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.0]
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


cluster_dict = {"15.0550003052":[], "15.0550003052_arrows":[]}

cluster_dict["15.0550003052"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-50.0), float(-43.0), float(-4.0), float(1.0)]

cluster_dict["15.0550003052_arrows"] += cgo_arrow([-50.0,-43.0,-4.0], [-50.188,-40.546,-5.827], color="blue red", name="Arrows_15.0550003052_1")

cluster_dict["15.0550003052"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-45.5), float(-40.0), float(4.5), float(1.0)]

cluster_dict["15.0550003052_arrows"] += cgo_arrow([-45.5,-40.0,4.5], [-43.549,-37.887,3.576], color="blue red", name="Arrows_15.0550003052_2")

cluster_dict["15.0550003052"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-44.5), float(-38.0), float(8.5), float(1.0)]

cluster_dict["15.0550003052_arrows"] += cgo_arrow([-44.5,-38.0,8.5], [-42.285,-36.612,9.774], color="blue red", name="Arrows_15.0550003052_3")

cluster_dict["15.0550003052"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-41.5), float(-45.0), float(6.5), float(1.0)]

cluster_dict["15.0550003052_arrows"] += cgo_arrow([-41.5,-45.0,6.5], [-40.482,-47.349,4.858], color="blue red", name="Arrows_15.0550003052_4")

cluster_dict["15.0550003052"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-45.7555550896), float(-41.9887428509), float(3.62724719548), float(1.0)]


cluster_dict["15.0550003052"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-39.8715765114), float(-47.007820649), float(-6.17774097672), float(1.0)]


cluster_dict["15.0550003052"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-39.6541155621), float(-40.3588332935), float(-2.56636740758), float(1.0)]


cluster_dict["15.0550003052"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-48.5), float(-44.0), float(-5.5), float(1.0)]

cluster_dict["15.0550003052_arrows"] += cgo_arrow([-48.5,-44.0,-5.5], [-48.805,-41.413,-7.401], color="red blue", name="Arrows_15.0550003052_5")

cluster_dict["15.0550003052"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-45.5), float(-42.0), float(-0.5), float(1.0)]

cluster_dict["15.0550003052_arrows"] += cgo_arrow([-45.5,-42.0,-0.5], [-45.326,-45.543,-3.379], color="red blue", name="Arrows_15.0550003052_6")

cluster_dict["15.0550003052"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-44.5), float(-39.5), float(10.0), float(1.0)]

cluster_dict["15.0550003052_arrows"] += cgo_arrow([-44.5,-39.5,10.0], [-42.638,-37.974,11.517], color="red blue", name="Arrows_15.0550003052_7")

cluster_dict["15.0550003052"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-41.0), float(-45.5), float(9.0), float(1.0)]

cluster_dict["15.0550003052_arrows"] += cgo_arrow([-41.0,-45.5,9.0], [-38.284,-43.927,9.57], color="red blue", name="Arrows_15.0550003052_8")

cluster_dict["15.0550003052"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-40.0), float(-41.5), float(3.0), float(1.0)]

cluster_dict["15.0550003052_arrows"] += cgo_arrow([-40.0,-41.5,3.0], [-40.8,-37.246,4.185], color="red blue", name="Arrows_15.0550003052_9")

cmd.load_cgo(cluster_dict["15.0550003052"], "Features_15.0550003052", 1)
cmd.load_cgo(cluster_dict["15.0550003052_arrows"], "Arrows_15.0550003052")
cmd.set("transparency", 0.2,"Features_15.0550003052")
cmd.group("Pharmacophore_15.0550003052", members="Features_15.0550003052")
cmd.group("Pharmacophore_15.0550003052", members="Arrows_15.0550003052")

if dirpath:
    f = join(dirpath, "9/label_threshold_15.0550003052.mol2")
else:
    f = "9/label_threshold_15.0550003052.mol2"

cmd.load(f, 'label_threshold_15.0550003052')
cmd.hide('everything', 'label_threshold_15.0550003052')
cmd.label("label_threshold_15.0550003052", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.0550003052', members= 'label_threshold_15.0550003052')


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
    f = join(dirpath, "10/label_threshold_0.6.mol2")
else:
    f = "10/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
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


cluster_dict = {"14.5399999619":[], "14.5399999619_arrows":[]}

cluster_dict["14.5399999619"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-55.1306636296), float(3.65354011939), float(-20.670186992), float(1.0)]


cmd.load_cgo(cluster_dict["14.5399999619"], "Features_14.5399999619", 1)
cmd.load_cgo(cluster_dict["14.5399999619_arrows"], "Arrows_14.5399999619")
cmd.set("transparency", 0.2,"Features_14.5399999619")
cmd.group("Pharmacophore_14.5399999619", members="Features_14.5399999619")
cmd.group("Pharmacophore_14.5399999619", members="Arrows_14.5399999619")

if dirpath:
    f = join(dirpath, "10/label_threshold_14.5399999619.mol2")
else:
    f = "10/label_threshold_14.5399999619.mol2"

cmd.load(f, 'label_threshold_14.5399999619')
cmd.hide('everything', 'label_threshold_14.5399999619')
cmd.label("label_threshold_14.5399999619", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.5399999619', members= 'label_threshold_14.5399999619')


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
    f = join(dirpath, "11/label_threshold_12.1.mol2")
else:
    f = "11/label_threshold_12.1.mol2"

cmd.load(f, 'label_threshold_12.1')
cmd.hide('everything', 'label_threshold_12.1')
cmd.label("label_threshold_12.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.1]
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


cluster_dict = {"14.5030002594":[], "14.5030002594_arrows":[]}

cluster_dict["14.5030002594"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-46.5), float(1.0), float(-1.0), float(1.0)]

cluster_dict["14.5030002594_arrows"] += cgo_arrow([-46.5,1.0,-1.0], [-46.107,1.785,-3.582], color="blue red", name="Arrows_14.5030002594_1")

cluster_dict["14.5030002594"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-44.0), float(-6.0), float(-8.0), float(1.0)]

cluster_dict["14.5030002594_arrows"] += cgo_arrow([-44.0,-6.0,-8.0], [-45.13,-2.952,-7.642], color="blue red", name="Arrows_14.5030002594_2")

cluster_dict["14.5030002594"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-42.0), float(-1.5), float(2.0), float(1.0)]

cluster_dict["14.5030002594_arrows"] += cgo_arrow([-42.0,-1.5,2.0], [-43.901,-0.73,4.3], color="blue red", name="Arrows_14.5030002594_3")

cluster_dict["14.5030002594"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-41.5553712114), float(-7.22031688956), float(-4.9149323194), float(1.0)]


cluster_dict["14.5030002594"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-43.2190405341), float(1.4246843833), float(0.711779359601), float(1.0)]


cluster_dict["14.5030002594"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-35.4462097429), float(-6.72792468136), float(-4.99769581104), float(1.0)]


cluster_dict["14.5030002594"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-46.0), float(-6.5), float(-10.0), float(1.0)]

cluster_dict["14.5030002594_arrows"] += cgo_arrow([-46.0,-6.5,-10.0], [-48.691,-5.672,-10.851], color="red blue", name="Arrows_14.5030002594_4")

cluster_dict["14.5030002594"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-45.0), float(4.0), float(-0.5), float(1.0)]

cluster_dict["14.5030002594_arrows"] += cgo_arrow([-45.0,4.0,-0.5], [-44.418,4.316,-3.388], color="red blue", name="Arrows_14.5030002594_5")

cluster_dict["14.5030002594"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-43.0), float(-5.5), float(-4.5), float(1.0)]

cluster_dict["14.5030002594_arrows"] += cgo_arrow([-43.0,-5.5,-4.5], [-46.033,-5.35,-4.985], color="red blue", name="Arrows_14.5030002594_6")

cluster_dict["14.5030002594"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-42.5), float(-2.0), float(0.5), float(1.0)]

cluster_dict["14.5030002594_arrows"] += cgo_arrow([-42.5,-2.0,0.5], [-44.511,-3.942,1.429], color="red blue", name="Arrows_14.5030002594_7")

cluster_dict["14.5030002594"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-36.5), float(-6.0), float(-7.0), float(1.0)]

cluster_dict["14.5030002594_arrows"] += cgo_arrow([-36.5,-6.0,-7.0], [-35.008,-3.384,-5.985], color="red blue", name="Arrows_14.5030002594_8")

cmd.load_cgo(cluster_dict["14.5030002594"], "Features_14.5030002594", 1)
cmd.load_cgo(cluster_dict["14.5030002594_arrows"], "Arrows_14.5030002594")
cmd.set("transparency", 0.2,"Features_14.5030002594")
cmd.group("Pharmacophore_14.5030002594", members="Features_14.5030002594")
cmd.group("Pharmacophore_14.5030002594", members="Arrows_14.5030002594")

if dirpath:
    f = join(dirpath, "11/label_threshold_14.5030002594.mol2")
else:
    f = "11/label_threshold_14.5030002594.mol2"

cmd.load(f, 'label_threshold_14.5030002594')
cmd.hide('everything', 'label_threshold_14.5030002594')
cmd.label("label_threshold_14.5030002594", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.5030002594', members= 'label_threshold_14.5030002594')


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
    f = join(dirpath, "12/label_threshold_9.3.mol2")
else:
    f = "12/label_threshold_9.3.mol2"

cmd.load(f, 'label_threshold_9.3')
cmd.hide('everything', 'label_threshold_9.3')
cmd.label("label_threshold_9.3", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [9.3]
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


cluster_dict = {"14.179500103":[], "14.179500103_arrows":[]}

cluster_dict["14.179500103"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-42.5), float(-18.0), float(13.5), float(1.0)]

cluster_dict["14.179500103_arrows"] += cgo_arrow([-42.5,-18.0,13.5], [-44.517,-17.591,11.795], color="blue red", name="Arrows_14.179500103_1")

cluster_dict["14.179500103"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-42.5), float(-29.5), float(15.0), float(1.0)]

cluster_dict["14.179500103_arrows"] += cgo_arrow([-42.5,-29.5,15.0], [-43.481,-30.657,17.066], color="blue red", name="Arrows_14.179500103_2")

cluster_dict["14.179500103"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-41.0), float(-19.5), float(8.5), float(1.0)]

cluster_dict["14.179500103_arrows"] += cgo_arrow([-41.0,-19.5,8.5], [-39.235,-21.84,8.702], color="blue red", name="Arrows_14.179500103_3")

cluster_dict["14.179500103"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-41.5), float(-26.0), float(12.0), float(1.0)]

cluster_dict["14.179500103_arrows"] += cgo_arrow([-41.5,-26.0,12.0], [-43.052,-25.27,9.512], color="blue red", name="Arrows_14.179500103_4")

cluster_dict["14.179500103"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-40.5), float(-21.0), float(11.0), float(1.0)]

cluster_dict["14.179500103_arrows"] += cgo_arrow([-40.5,-21.0,11.0], [-39.235,-21.84,8.702], color="blue red", name="Arrows_14.179500103_5")

cluster_dict["14.179500103"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-40.0), float(-16.0), float(12.0), float(1.0)]

cluster_dict["14.179500103_arrows"] += cgo_arrow([-40.0,-16.0,12.0], [-38.494,-13.8,13.577], color="blue red", name="Arrows_14.179500103_6")

cluster_dict["14.179500103"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-41.0821713454), float(-23.0202621188), float(12.3023608156), float(1.0)]


cluster_dict["14.179500103"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-43.5), float(-27.5), float(9.5), float(1.0)]

cluster_dict["14.179500103_arrows"] += cgo_arrow([-43.5,-27.5,9.5], [-43.614,-30.474,9.596], color="red blue", name="Arrows_14.179500103_7")

cluster_dict["14.179500103"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-39.5), float(-19.5), float(8.0), float(1.0)]

cluster_dict["14.179500103_arrows"] += cgo_arrow([-39.5,-19.5,8.0], [-39.206,-16.839,6.614], color="red blue", name="Arrows_14.179500103_8")

cmd.load_cgo(cluster_dict["14.179500103"], "Features_14.179500103", 1)
cmd.load_cgo(cluster_dict["14.179500103_arrows"], "Arrows_14.179500103")
cmd.set("transparency", 0.2,"Features_14.179500103")
cmd.group("Pharmacophore_14.179500103", members="Features_14.179500103")
cmd.group("Pharmacophore_14.179500103", members="Arrows_14.179500103")

if dirpath:
    f = join(dirpath, "12/label_threshold_14.179500103.mol2")
else:
    f = "12/label_threshold_14.179500103.mol2"

cmd.load(f, 'label_threshold_14.179500103')
cmd.hide('everything', 'label_threshold_14.179500103')
cmd.label("label_threshold_14.179500103", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.179500103', members= 'label_threshold_14.179500103')


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
    f = join(dirpath, "13/label_threshold_7.7.mol2")
else:
    f = "13/label_threshold_7.7.mol2"

cmd.load(f, 'label_threshold_7.7')
cmd.hide('everything', 'label_threshold_7.7')
cmd.label("label_threshold_7.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [7.7]
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


cluster_dict = {"13.8030004501":[], "13.8030004501_arrows":[]}

cluster_dict["13.8030004501"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-42.5), float(-18.0), float(13.5), float(1.0)]

cluster_dict["13.8030004501_arrows"] += cgo_arrow([-42.5,-18.0,13.5], [-44.517,-17.591,11.795], color="blue red", name="Arrows_13.8030004501_1")

cluster_dict["13.8030004501"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-42.5), float(-29.5), float(15.0), float(1.0)]

cluster_dict["13.8030004501_arrows"] += cgo_arrow([-42.5,-29.5,15.0], [-43.481,-30.657,17.066], color="blue red", name="Arrows_13.8030004501_2")

cluster_dict["13.8030004501"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-41.0), float(-19.5), float(8.5), float(1.0)]

cluster_dict["13.8030004501_arrows"] += cgo_arrow([-41.0,-19.5,8.5], [-39.235,-21.84,8.702], color="blue red", name="Arrows_13.8030004501_3")

cluster_dict["13.8030004501"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-41.5), float(-26.0), float(12.0), float(1.0)]

cluster_dict["13.8030004501_arrows"] += cgo_arrow([-41.5,-26.0,12.0], [-43.052,-25.27,9.512], color="blue red", name="Arrows_13.8030004501_4")

cluster_dict["13.8030004501"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-40.5), float(-21.0), float(11.0), float(1.0)]

cluster_dict["13.8030004501_arrows"] += cgo_arrow([-40.5,-21.0,11.0], [-39.235,-21.84,8.702], color="blue red", name="Arrows_13.8030004501_5")

cluster_dict["13.8030004501"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-40.0), float(-16.5), float(12.0), float(1.0)]

cluster_dict["13.8030004501_arrows"] += cgo_arrow([-40.0,-16.5,12.0], [-38.043,-17.594,13.932], color="blue red", name="Arrows_13.8030004501_6")

cluster_dict["13.8030004501"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-41.0498559415), float(-23.9613504427), float(12.4552912728), float(1.0)]


cluster_dict["13.8030004501"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-43.5), float(-27.5), float(9.5), float(1.0)]

cluster_dict["13.8030004501_arrows"] += cgo_arrow([-43.5,-27.5,9.5], [-43.614,-30.474,9.596], color="red blue", name="Arrows_13.8030004501_7")

cluster_dict["13.8030004501"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-39.5), float(-19.5), float(8.0), float(1.0)]

cluster_dict["13.8030004501_arrows"] += cgo_arrow([-39.5,-19.5,8.0], [-39.206,-16.839,6.614], color="red blue", name="Arrows_13.8030004501_8")

cmd.load_cgo(cluster_dict["13.8030004501"], "Features_13.8030004501", 1)
cmd.load_cgo(cluster_dict["13.8030004501_arrows"], "Arrows_13.8030004501")
cmd.set("transparency", 0.2,"Features_13.8030004501")
cmd.group("Pharmacophore_13.8030004501", members="Features_13.8030004501")
cmd.group("Pharmacophore_13.8030004501", members="Arrows_13.8030004501")

if dirpath:
    f = join(dirpath, "13/label_threshold_13.8030004501.mol2")
else:
    f = "13/label_threshold_13.8030004501.mol2"

cmd.load(f, 'label_threshold_13.8030004501')
cmd.hide('everything', 'label_threshold_13.8030004501')
cmd.label("label_threshold_13.8030004501", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.8030004501', members= 'label_threshold_13.8030004501')


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


cluster_dict = {"13.4099998474":[], "13.4099998474_arrows":[]}

cluster_dict["13.4099998474"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-45.5), float(-40.0), float(4.5), float(1.0)]

cluster_dict["13.4099998474_arrows"] += cgo_arrow([-45.5,-40.0,4.5], [-43.549,-37.887,3.576], color="blue red", name="Arrows_13.4099998474_1")

cluster_dict["13.4099998474"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-44.5), float(-38.0), float(8.5), float(1.0)]

cluster_dict["13.4099998474_arrows"] += cgo_arrow([-44.5,-38.0,8.5], [-42.285,-36.612,9.774], color="blue red", name="Arrows_13.4099998474_2")

cluster_dict["13.4099998474"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-41.5), float(-45.0), float(6.5), float(1.0)]

cluster_dict["13.4099998474_arrows"] += cgo_arrow([-41.5,-45.0,6.5], [-40.482,-47.349,4.858], color="blue red", name="Arrows_13.4099998474_3")

cluster_dict["13.4099998474"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-43.2421286992), float(-43.6246792639), float(7.15199059638), float(1.0)]


cluster_dict["13.4099998474"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-44.5), float(-39.5), float(10.0), float(1.0)]

cluster_dict["13.4099998474_arrows"] += cgo_arrow([-44.5,-39.5,10.0], [-42.638,-37.974,11.517], color="red blue", name="Arrows_13.4099998474_4")

cluster_dict["13.4099998474"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-41.0), float(-45.5), float(9.0), float(1.0)]

cluster_dict["13.4099998474_arrows"] += cgo_arrow([-41.0,-45.5,9.0], [-38.284,-43.927,9.57], color="red blue", name="Arrows_13.4099998474_5")

cluster_dict["13.4099998474"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-40.0), float(-41.5), float(3.0), float(1.0)]

cluster_dict["13.4099998474_arrows"] += cgo_arrow([-40.0,-41.5,3.0], [-40.8,-37.246,4.185], color="red blue", name="Arrows_13.4099998474_6")

cmd.load_cgo(cluster_dict["13.4099998474"], "Features_13.4099998474", 1)
cmd.load_cgo(cluster_dict["13.4099998474_arrows"], "Arrows_13.4099998474")
cmd.set("transparency", 0.2,"Features_13.4099998474")
cmd.group("Pharmacophore_13.4099998474", members="Features_13.4099998474")
cmd.group("Pharmacophore_13.4099998474", members="Arrows_13.4099998474")

if dirpath:
    f = join(dirpath, "14/label_threshold_13.4099998474.mol2")
else:
    f = "14/label_threshold_13.4099998474.mol2"

cmd.load(f, 'label_threshold_13.4099998474')
cmd.hide('everything', 'label_threshold_13.4099998474')
cmd.label("label_threshold_13.4099998474", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.4099998474', members= 'label_threshold_13.4099998474')


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
    f = join(dirpath, "15/label_threshold_9.1.mol2")
else:
    f = "15/label_threshold_9.1.mol2"

cmd.load(f, 'label_threshold_9.1')
cmd.hide('everything', 'label_threshold_9.1')
cmd.label("label_threshold_9.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [9.1]
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


cluster_dict = {"13.1870002747":[], "13.1870002747_arrows":[]}

cluster_dict["13.1870002747"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-67.5), float(-32.5), float(-22.0), float(1.0)]

cluster_dict["13.1870002747_arrows"] += cgo_arrow([-67.5,-32.5,-22.0], [-68.022,-35.07,-21.366], color="blue red", name="Arrows_13.1870002747_1")

cluster_dict["13.1870002747"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-66.0), float(-27.5), float(-19.0), float(1.0)]

cluster_dict["13.1870002747_arrows"] += cgo_arrow([-66.0,-27.5,-19.0], [-65.124,-27.855,-15.984], color="blue red", name="Arrows_13.1870002747_2")

cluster_dict["13.1870002747"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-69.6711697986), float(-30.3276067532), float(-21.3986079184), float(1.0)]


cluster_dict["13.1870002747"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-73.0), float(-29.5), float(-22.5), float(1.0)]


cluster_dict["13.1870002747"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-69.0), float(-27.5), float(-20.5), float(1.0)]

cluster_dict["13.1870002747_arrows"] += cgo_arrow([-69.0,-27.5,-20.5], [-66.473,-25.039,-23.597], color="red blue", name="Arrows_13.1870002747_3")

cluster_dict["13.1870002747"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-69.0), float(-27.5), float(-20.5), float(1.0)]

cluster_dict["13.1870002747_arrows"] += cgo_arrow([-69.0,-27.5,-20.5], [-66.473,-25.039,-23.597], color="red blue", name="Arrows_13.1870002747_4")

cluster_dict["13.1870002747"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-65.0), float(-28.5), float(-21.0), float(1.0)]

cluster_dict["13.1870002747_arrows"] += cgo_arrow([-65.0,-28.5,-21.0], [-63.516,-24.906,-21.648], color="red blue", name="Arrows_13.1870002747_5")

cluster_dict["13.1870002747"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-61.5), float(-37.0), float(-23.0), float(1.0)]

cluster_dict["13.1870002747_arrows"] += cgo_arrow([-61.5,-37.0,-23.0], [-61.368,-40.494,-20.992], color="red blue", name="Arrows_13.1870002747_6")

cmd.load_cgo(cluster_dict["13.1870002747"], "Features_13.1870002747", 1)
cmd.load_cgo(cluster_dict["13.1870002747_arrows"], "Arrows_13.1870002747")
cmd.set("transparency", 0.2,"Features_13.1870002747")
cmd.group("Pharmacophore_13.1870002747", members="Features_13.1870002747")
cmd.group("Pharmacophore_13.1870002747", members="Arrows_13.1870002747")

if dirpath:
    f = join(dirpath, "15/label_threshold_13.1870002747.mol2")
else:
    f = "15/label_threshold_13.1870002747.mol2"

cmd.load(f, 'label_threshold_13.1870002747')
cmd.hide('everything', 'label_threshold_13.1870002747')
cmd.label("label_threshold_13.1870002747", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.1870002747', members= 'label_threshold_13.1870002747')


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
    f = join(dirpath, "16/label_threshold_8.2.mol2")
else:
    f = "16/label_threshold_8.2.mol2"

cmd.load(f, 'label_threshold_8.2')
cmd.hide('everything', 'label_threshold_8.2')
cmd.label("label_threshold_8.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [8.2]
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


cluster_dict = {"12.642999649":[], "12.642999649_arrows":[]}

cluster_dict["12.642999649"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-18.0), float(-30.0), float(-19.5), float(1.0)]

cluster_dict["12.642999649_arrows"] += cgo_arrow([-18.0,-30.0,-19.5], [-18.656,-28.605,-16.965], color="blue red", name="Arrows_12.642999649_1")

cluster_dict["12.642999649"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-16.5), float(-26.5), float(-24.5), float(1.0)]

cluster_dict["12.642999649_arrows"] += cgo_arrow([-16.5,-26.5,-24.5], [-16.058,-23.495,-24.824], color="blue red", name="Arrows_12.642999649_2")

cluster_dict["12.642999649"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-16.9418486874), float(-27.3068652995), float(-24.891623744), float(1.0)]


cluster_dict["12.642999649"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-20.5), float(-20.4375), float(-28.625), float(1.0)]


cluster_dict["12.642999649"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-22.0), float(-26.5), float(-28.5), float(1.0)]

cluster_dict["12.642999649_arrows"] += cgo_arrow([-22.0,-26.5,-28.5], [-23.407,-29.252,-29.261], color="red blue", name="Arrows_12.642999649_3")

cluster_dict["12.642999649"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-18.0), float(-29.5), float(-21.0), float(1.0)]

cluster_dict["12.642999649_arrows"] += cgo_arrow([-18.0,-29.5,-21.0], [-20.194,-33.373,-21.039], color="red blue", name="Arrows_12.642999649_4")

cluster_dict["12.642999649"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-14.0), float(-30.5), float(-22.0), float(1.0)]

cluster_dict["12.642999649_arrows"] += cgo_arrow([-14.0,-30.5,-22.0], [-15.12,-33.631,-24.68], color="red blue", name="Arrows_12.642999649_5")

cluster_dict["12.642999649"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-15.5), float(-29.0), float(-19.5), float(1.0)]

cluster_dict["12.642999649_arrows"] += cgo_arrow([-15.5,-29.0,-19.5], [-16.441,-27.124,-16.889], color="red blue", name="Arrows_12.642999649_6")

cmd.load_cgo(cluster_dict["12.642999649"], "Features_12.642999649", 1)
cmd.load_cgo(cluster_dict["12.642999649_arrows"], "Arrows_12.642999649")
cmd.set("transparency", 0.2,"Features_12.642999649")
cmd.group("Pharmacophore_12.642999649", members="Features_12.642999649")
cmd.group("Pharmacophore_12.642999649", members="Arrows_12.642999649")

if dirpath:
    f = join(dirpath, "16/label_threshold_12.642999649.mol2")
else:
    f = "16/label_threshold_12.642999649.mol2"

cmd.load(f, 'label_threshold_12.642999649')
cmd.hide('everything', 'label_threshold_12.642999649')
cmd.label("label_threshold_12.642999649", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.642999649', members= 'label_threshold_12.642999649')


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
    f = join(dirpath, "17/label_threshold_4.9.mol2")
else:
    f = "17/label_threshold_4.9.mol2"

cmd.load(f, 'label_threshold_4.9')
cmd.hide('everything', 'label_threshold_4.9')
cmd.label("label_threshold_4.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [4.9]
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


cluster_dict = {"11.6500000954":[], "11.6500000954_arrows":[]}

cluster_dict["11.6500000954"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-46.0), float(-12.5), float(-12.5), float(1.0)]

cluster_dict["11.6500000954_arrows"] += cgo_arrow([-46.0,-12.5,-12.5], [-47.399,-9.957,-13.153], color="blue red", name="Arrows_11.6500000954_1")

cluster_dict["11.6500000954"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-42.5), float(-15.0), float(-9.0), float(1.0)]

cluster_dict["11.6500000954_arrows"] += cgo_arrow([-42.5,-15.0,-9.0], [-43.537,-16.731,-9.097], color="blue red", name="Arrows_11.6500000954_2")

cluster_dict["11.6500000954"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-43.4646941459), float(-11.5691868735), float(-13.9507856138), float(1.0)]


cluster_dict["11.6500000954"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-42.9598580415), float(-16.1190263759), float(-19.7553026186), float(1.0)]


cluster_dict["11.6500000954"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-46.0), float(-16.0), float(-9.0), float(1.0)]

cluster_dict["11.6500000954_arrows"] += cgo_arrow([-46.0,-16.0,-9.0], [-43.537,-16.731,-9.097], color="red blue", name="Arrows_11.6500000954_3")

cluster_dict["11.6500000954"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-40.5), float(-11.0), float(-12.5), float(1.0)]

cluster_dict["11.6500000954_arrows"] += cgo_arrow([-40.5,-11.0,-12.5], [-39.821,-14.092,-14.343], color="red blue", name="Arrows_11.6500000954_4")

cmd.load_cgo(cluster_dict["11.6500000954"], "Features_11.6500000954", 1)
cmd.load_cgo(cluster_dict["11.6500000954_arrows"], "Arrows_11.6500000954")
cmd.set("transparency", 0.2,"Features_11.6500000954")
cmd.group("Pharmacophore_11.6500000954", members="Features_11.6500000954")
cmd.group("Pharmacophore_11.6500000954", members="Arrows_11.6500000954")

if dirpath:
    f = join(dirpath, "17/label_threshold_11.6500000954.mol2")
else:
    f = "17/label_threshold_11.6500000954.mol2"

cmd.load(f, 'label_threshold_11.6500000954')
cmd.hide('everything', 'label_threshold_11.6500000954')
cmd.label("label_threshold_11.6500000954", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.6500000954', members= 'label_threshold_11.6500000954')


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
    f = join(dirpath, "18/label_threshold_0.6.mol2")
else:
    f = "18/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
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


cluster_dict = {"10.9989995956":[], "10.9989995956_arrows":[]}

cluster_dict["10.9989995956"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-63.5), float(-25.0), float(8.5), float(1.0)]

cluster_dict["10.9989995956_arrows"] += cgo_arrow([-63.5,-25.0,8.5], [-60.695,-25.422,9.101], color="blue red", name="Arrows_10.9989995956_1")

cluster_dict["10.9989995956"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-64.6215491661), float(-23.1440237005), float(6.21166859727), float(1.0)]


cluster_dict["10.9989995956"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-65.5), float(-22.0), float(3.0), float(1.0)]

cluster_dict["10.9989995956_arrows"] += cgo_arrow([-65.5,-22.0,3.0], [-67.529,-19.586,2.278], color="red blue", name="Arrows_10.9989995956_2")

cluster_dict["10.9989995956"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-62.5), float(-26.5), float(7.5), float(1.0)]

cluster_dict["10.9989995956_arrows"] += cgo_arrow([-62.5,-26.5,7.5], [-63.622,-28.096,5.158], color="red blue", name="Arrows_10.9989995956_3")

cmd.load_cgo(cluster_dict["10.9989995956"], "Features_10.9989995956", 1)
cmd.load_cgo(cluster_dict["10.9989995956_arrows"], "Arrows_10.9989995956")
cmd.set("transparency", 0.2,"Features_10.9989995956")
cmd.group("Pharmacophore_10.9989995956", members="Features_10.9989995956")
cmd.group("Pharmacophore_10.9989995956", members="Arrows_10.9989995956")

if dirpath:
    f = join(dirpath, "18/label_threshold_10.9989995956.mol2")
else:
    f = "18/label_threshold_10.9989995956.mol2"

cmd.load(f, 'label_threshold_10.9989995956')
cmd.hide('everything', 'label_threshold_10.9989995956')
cmd.label("label_threshold_10.9989995956", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.9989995956', members= 'label_threshold_10.9989995956')


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
    f = join(dirpath, "19/label_threshold_6.6.mol2")
else:
    f = "19/label_threshold_6.6.mol2"

cmd.load(f, 'label_threshold_6.6')
cmd.hide('everything', 'label_threshold_6.6')
cmd.label("label_threshold_6.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [6.6]
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


cluster_dict = {"10.1180000305":[], "10.1180000305_arrows":[]}

cluster_dict["10.1180000305"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-50.0), float(-43.0), float(-4.0), float(1.0)]

cluster_dict["10.1180000305_arrows"] += cgo_arrow([-50.0,-43.0,-4.0], [-50.188,-40.546,-5.827], color="blue red", name="Arrows_10.1180000305_1")

cluster_dict["10.1180000305"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-50.4920008942), float(-43.4375509877), float(-3.40624511081), float(1.0)]


cluster_dict["10.1180000305"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-45.5), float(-42.5), float(-0.5), float(1.0)]

cluster_dict["10.1180000305_arrows"] += cgo_arrow([-45.5,-42.5,-0.5], [-45.326,-45.543,-3.379], color="red blue", name="Arrows_10.1180000305_2")

cluster_dict["10.1180000305"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-48.5), float(-44.0), float(-5.5), float(1.0)]

cluster_dict["10.1180000305_arrows"] += cgo_arrow([-48.5,-44.0,-5.5], [-48.805,-41.413,-7.401], color="red blue", name="Arrows_10.1180000305_3")

cmd.load_cgo(cluster_dict["10.1180000305"], "Features_10.1180000305", 1)
cmd.load_cgo(cluster_dict["10.1180000305_arrows"], "Arrows_10.1180000305")
cmd.set("transparency", 0.2,"Features_10.1180000305")
cmd.group("Pharmacophore_10.1180000305", members="Features_10.1180000305")
cmd.group("Pharmacophore_10.1180000305", members="Arrows_10.1180000305")

if dirpath:
    f = join(dirpath, "19/label_threshold_10.1180000305.mol2")
else:
    f = "19/label_threshold_10.1180000305.mol2"

cmd.load(f, 'label_threshold_10.1180000305')
cmd.hide('everything', 'label_threshold_10.1180000305')
cmd.label("label_threshold_10.1180000305", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.1180000305', members= 'label_threshold_10.1180000305')


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
    f = join(dirpath, "20/label_threshold_4.1.mol2")
else:
    f = "20/label_threshold_4.1.mol2"

cmd.load(f, 'label_threshold_4.1')
cmd.hide('everything', 'label_threshold_4.1')
cmd.label("label_threshold_4.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [4.1]
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


cluster_dict = {"9.55900001526":[], "9.55900001526_arrows":[]}

cluster_dict["9.55900001526"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-17.0), float(-26.5), float(-25.0), float(1.0)]

cluster_dict["9.55900001526_arrows"] += cgo_arrow([-17.0,-26.5,-25.0], [-16.058,-23.495,-24.824], color="blue red", name="Arrows_9.55900001526_1")

cluster_dict["9.55900001526"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-21.1398296359), float(-24.1252091154), float(-27.5085898322), float(1.0)]


cluster_dict["9.55900001526"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-22.0), float(-26.5), float(-28.5), float(1.0)]

cluster_dict["9.55900001526_arrows"] += cgo_arrow([-22.0,-26.5,-28.5], [-23.407,-29.252,-29.261], color="red blue", name="Arrows_9.55900001526_2")

cluster_dict["9.55900001526"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-20.0), float(-29.5), float(-23.0), float(1.0)]

cluster_dict["9.55900001526_arrows"] += cgo_arrow([-20.0,-29.5,-23.0], [-22.142,-30.633,-23.814], color="red blue", name="Arrows_9.55900001526_3")

cmd.load_cgo(cluster_dict["9.55900001526"], "Features_9.55900001526", 1)
cmd.load_cgo(cluster_dict["9.55900001526_arrows"], "Arrows_9.55900001526")
cmd.set("transparency", 0.2,"Features_9.55900001526")
cmd.group("Pharmacophore_9.55900001526", members="Features_9.55900001526")
cmd.group("Pharmacophore_9.55900001526", members="Arrows_9.55900001526")

if dirpath:
    f = join(dirpath, "20/label_threshold_9.55900001526.mol2")
else:
    f = "20/label_threshold_9.55900001526.mol2"

cmd.load(f, 'label_threshold_9.55900001526')
cmd.hide('everything', 'label_threshold_9.55900001526')
cmd.label("label_threshold_9.55900001526", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.55900001526', members= 'label_threshold_9.55900001526')


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
    f = join(dirpath, "21/label_threshold_2.4.mol2")
else:
    f = "21/label_threshold_2.4.mol2"

cmd.load(f, 'label_threshold_2.4')
cmd.hide('everything', 'label_threshold_2.4')
cmd.label("label_threshold_2.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [2.4]
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


cluster_dict = {"9.48299980164":[], "9.48299980164_arrows":[]}

cluster_dict["9.48299980164"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-16.5), float(-42.5), float(-3.5), float(1.0)]

cluster_dict["9.48299980164_arrows"] += cgo_arrow([-16.5,-42.5,-3.5], [-19.472,-43.28,-3.522], color="blue red", name="Arrows_9.48299980164_1")

cluster_dict["9.48299980164"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-16.5), float(-42.5), float(-3.5), float(1.0)]

cluster_dict["9.48299980164_arrows"] += cgo_arrow([-16.5,-42.5,-3.5], [-19.472,-43.28,-3.522], color="blue red", name="Arrows_9.48299980164_2")

cluster_dict["9.48299980164"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-15.2802692514), float(-43.5308221753), float(-2.089903114), float(1.0)]


cmd.load_cgo(cluster_dict["9.48299980164"], "Features_9.48299980164", 1)
cmd.load_cgo(cluster_dict["9.48299980164_arrows"], "Arrows_9.48299980164")
cmd.set("transparency", 0.2,"Features_9.48299980164")
cmd.group("Pharmacophore_9.48299980164", members="Features_9.48299980164")
cmd.group("Pharmacophore_9.48299980164", members="Arrows_9.48299980164")

if dirpath:
    f = join(dirpath, "21/label_threshold_9.48299980164.mol2")
else:
    f = "21/label_threshold_9.48299980164.mol2"

cmd.load(f, 'label_threshold_9.48299980164')
cmd.hide('everything', 'label_threshold_9.48299980164')
cmd.label("label_threshold_9.48299980164", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.48299980164', members= 'label_threshold_9.48299980164')


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
    f = join(dirpath, "22/label_threshold_2.8.mol2")
else:
    f = "22/label_threshold_2.8.mol2"

cmd.load(f, 'label_threshold_2.8')
cmd.hide('everything', 'label_threshold_2.8')
cmd.label("label_threshold_2.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [2.8]
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


cluster_dict = {"8.64200019836":[], "8.64200019836_arrows":[]}

cluster_dict["8.64200019836"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-68.0), float(-12.0), float(-9.0), float(1.0)]

cluster_dict["8.64200019836_arrows"] += cgo_arrow([-68.0,-12.0,-9.0], [-68.199,-9.877,-11.593], color="blue red", name="Arrows_8.64200019836_1")

cluster_dict["8.64200019836"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-68.7351044673), float(-8.91932182544), float(-8.06156854891), float(1.0)]


cluster_dict["8.64200019836"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-70.0), float(-10.0), float(-6.0), float(1.0)]

cluster_dict["8.64200019836_arrows"] += cgo_arrow([-70.0,-10.0,-6.0], [-68.145,-9.058,-3.925], color="red blue", name="Arrows_8.64200019836_2")

cluster_dict["8.64200019836"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-68.5), float(-10.0), float(-6.5), float(1.0)]

cluster_dict["8.64200019836_arrows"] += cgo_arrow([-68.5,-10.0,-6.5], [-68.145,-9.058,-3.925], color="red blue", name="Arrows_8.64200019836_3")

cmd.load_cgo(cluster_dict["8.64200019836"], "Features_8.64200019836", 1)
cmd.load_cgo(cluster_dict["8.64200019836_arrows"], "Arrows_8.64200019836")
cmd.set("transparency", 0.2,"Features_8.64200019836")
cmd.group("Pharmacophore_8.64200019836", members="Features_8.64200019836")
cmd.group("Pharmacophore_8.64200019836", members="Arrows_8.64200019836")

if dirpath:
    f = join(dirpath, "22/label_threshold_8.64200019836.mol2")
else:
    f = "22/label_threshold_8.64200019836.mol2"

cmd.load(f, 'label_threshold_8.64200019836')
cmd.hide('everything', 'label_threshold_8.64200019836')
cmd.label("label_threshold_8.64200019836", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_8.64200019836', members= 'label_threshold_8.64200019836')


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
    f = join(dirpath, "23/label_threshold_2.9.mol2")
else:
    f = "23/label_threshold_2.9.mol2"

cmd.load(f, 'label_threshold_2.9')
cmd.hide('everything', 'label_threshold_2.9')
cmd.label("label_threshold_2.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [2.9]
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


cluster_dict = {"8.47500038147":[], "8.47500038147_arrows":[]}

cluster_dict["8.47500038147"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-67.5), float(-32.5), float(-22.0), float(1.0)]

cluster_dict["8.47500038147_arrows"] += cgo_arrow([-67.5,-32.5,-22.0], [-68.022,-35.07,-21.366], color="blue red", name="Arrows_8.47500038147_1")

cluster_dict["8.47500038147"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-62.9489671694), float(-35.7018037634), float(-24.2745101674), float(1.0)]


cluster_dict["8.47500038147"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-64.5), float(-29.5), float(-20.0), float(1.0)]

cluster_dict["8.47500038147_arrows"] += cgo_arrow([-64.5,-29.5,-20.0], [-65.877,-33.139,-17.898], color="red blue", name="Arrows_8.47500038147_2")

cluster_dict["8.47500038147"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-61.5), float(-37.0), float(-23.0), float(1.0)]

cluster_dict["8.47500038147_arrows"] += cgo_arrow([-61.5,-37.0,-23.0], [-61.368,-40.494,-20.992], color="red blue", name="Arrows_8.47500038147_3")

cluster_dict["8.47500038147"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-61.5), float(-37.0), float(-23.0), float(1.0)]

cluster_dict["8.47500038147_arrows"] += cgo_arrow([-61.5,-37.0,-23.0], [-61.368,-40.494,-20.992], color="red blue", name="Arrows_8.47500038147_4")

cmd.load_cgo(cluster_dict["8.47500038147"], "Features_8.47500038147", 1)
cmd.load_cgo(cluster_dict["8.47500038147_arrows"], "Arrows_8.47500038147")
cmd.set("transparency", 0.2,"Features_8.47500038147")
cmd.group("Pharmacophore_8.47500038147", members="Features_8.47500038147")
cmd.group("Pharmacophore_8.47500038147", members="Arrows_8.47500038147")

if dirpath:
    f = join(dirpath, "23/label_threshold_8.47500038147.mol2")
else:
    f = "23/label_threshold_8.47500038147.mol2"

cmd.load(f, 'label_threshold_8.47500038147')
cmd.hide('everything', 'label_threshold_8.47500038147')
cmd.label("label_threshold_8.47500038147", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_8.47500038147', members= 'label_threshold_8.47500038147')


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
    f = join(dirpath, "24/label_threshold_0.6.mol2")
else:
    f = "24/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
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


cluster_dict = {"8.09200000763":[], "8.09200000763_arrows":[]}

cluster_dict["8.09200000763"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-18.7128345079), float(-23.9543313624), float(6.08145959903), float(1.0)]


cmd.load_cgo(cluster_dict["8.09200000763"], "Features_8.09200000763", 1)
cmd.load_cgo(cluster_dict["8.09200000763_arrows"], "Arrows_8.09200000763")
cmd.set("transparency", 0.2,"Features_8.09200000763")
cmd.group("Pharmacophore_8.09200000763", members="Features_8.09200000763")
cmd.group("Pharmacophore_8.09200000763", members="Arrows_8.09200000763")

if dirpath:
    f = join(dirpath, "24/label_threshold_8.09200000763.mol2")
else:
    f = "24/label_threshold_8.09200000763.mol2"

cmd.load(f, 'label_threshold_8.09200000763')
cmd.hide('everything', 'label_threshold_8.09200000763')
cmd.label("label_threshold_8.09200000763", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_8.09200000763', members= 'label_threshold_8.09200000763')


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
    f = join(dirpath, "25/label_threshold_0.6.mol2")
else:
    f = "25/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
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


cluster_dict = {"7.16400003433":[], "7.16400003433_arrows":[]}

cluster_dict["7.16400003433"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-49.2429939467), float(-40.7260773162), float(-16.2667229994), float(1.0)]


cluster_dict["7.16400003433"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-45.5762752007), float(-41.0231126617), float(-17.3648807227), float(1.0)]


cmd.load_cgo(cluster_dict["7.16400003433"], "Features_7.16400003433", 1)
cmd.load_cgo(cluster_dict["7.16400003433_arrows"], "Arrows_7.16400003433")
cmd.set("transparency", 0.2,"Features_7.16400003433")
cmd.group("Pharmacophore_7.16400003433", members="Features_7.16400003433")
cmd.group("Pharmacophore_7.16400003433", members="Arrows_7.16400003433")

if dirpath:
    f = join(dirpath, "25/label_threshold_7.16400003433.mol2")
else:
    f = "25/label_threshold_7.16400003433.mol2"

cmd.load(f, 'label_threshold_7.16400003433')
cmd.hide('everything', 'label_threshold_7.16400003433')
cmd.label("label_threshold_7.16400003433", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_7.16400003433', members= 'label_threshold_7.16400003433')


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


cluster_dict = {"6.17899990082":[], "6.17899990082_arrows":[]}

cluster_dict["6.17899990082"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-43.706668707), float(-31.7570001214), float(-25.5139998887), float(1.0)]


cmd.load_cgo(cluster_dict["6.17899990082"], "Features_6.17899990082", 1)
cmd.load_cgo(cluster_dict["6.17899990082_arrows"], "Arrows_6.17899990082")
cmd.set("transparency", 0.2,"Features_6.17899990082")
cmd.group("Pharmacophore_6.17899990082", members="Features_6.17899990082")
cmd.group("Pharmacophore_6.17899990082", members="Arrows_6.17899990082")

if dirpath:
    f = join(dirpath, "26/label_threshold_6.17899990082.mol2")
else:
    f = "26/label_threshold_6.17899990082.mol2"

cmd.load(f, 'label_threshold_6.17899990082')
cmd.hide('everything', 'label_threshold_6.17899990082')
cmd.label("label_threshold_6.17899990082", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_6.17899990082', members= 'label_threshold_6.17899990082')


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

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
