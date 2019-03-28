
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


cluster_dict = {"8.61400032043":[], "8.61400032043_arrows":[]}

cluster_dict["8.61400032043"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(16.9853224264), float(-30.4380859215), float(14.5), float(1.0)]


cluster_dict["8.61400032043"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(20.8122773335), float(-27.696702135), float(7.99495805385), float(1.0)]


cluster_dict["8.61400032043"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(23.1658746464), float(-27.9996146955), float(7.0), float(1.0)]


cmd.load_cgo(cluster_dict["8.61400032043"], "Features_8.61400032043", 1)
cmd.load_cgo(cluster_dict["8.61400032043_arrows"], "Arrows_8.61400032043")
cmd.set("transparency", 0.2,"Features_8.61400032043")
cmd.group("Pharmacophore_8.61400032043", members="Features_8.61400032043")
cmd.group("Pharmacophore_8.61400032043", members="Arrows_8.61400032043")

if dirpath:
    f = join(dirpath, "0/label_threshold_8.61400032043.mol2")
else:
    f = "0/label_threshold_8.61400032043.mol2"

cmd.load(f, 'label_threshold_8.61400032043')
cmd.hide('everything', 'label_threshold_8.61400032043')
cmd.label("label_threshold_8.61400032043", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_8.61400032043', members= 'label_threshold_8.61400032043')


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
    f = join(dirpath, "1/label_threshold_10.8.mol2")
else:
    f = "1/label_threshold_10.8.mol2"

cmd.load(f, 'label_threshold_10.8')
cmd.hide('everything', 'label_threshold_10.8')
cmd.label("label_threshold_10.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.8]
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


cluster_dict = {"17.9939994812":[], "17.9939994812_arrows":[]}

cluster_dict["17.9939994812"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(0.5), float(-37.0), float(18.0), float(1.0)]

cluster_dict["17.9939994812_arrows"] += cgo_arrow([0.5,-37.0,18.0], [0.174,-36.164,15.33], color="blue red", name="Arrows_17.9939994812_1")

cluster_dict["17.9939994812"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(5.0), float(-37.0), float(18.5), float(1.0)]

cluster_dict["17.9939994812_arrows"] += cgo_arrow([5.0,-37.0,18.5], [6.939,-35.33,17.433], color="blue red", name="Arrows_17.9939994812_2")

cluster_dict["17.9939994812"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(1.81156224728), float(-36.4096162446), float(21.4481129302), float(1.0)]


cluster_dict["17.9939994812"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-1.0), float(-36.5), float(21.0), float(1.0)]

cluster_dict["17.9939994812_arrows"] += cgo_arrow([-1.0,-36.5,21.0], [-2.138,-31.939,19.534], color="red blue", name="Arrows_17.9939994812_3")

cluster_dict["17.9939994812"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-2.5), float(-34.5), float(23.0), float(1.0)]

cluster_dict["17.9939994812_arrows"] += cgo_arrow([-2.5,-34.5,23.0], [-5.169,-31.989,21.428], color="red blue", name="Arrows_17.9939994812_4")

cluster_dict["17.9939994812"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-1.0), float(-39.0), float(18.5), float(1.0)]

cluster_dict["17.9939994812_arrows"] += cgo_arrow([-1.0,-39.0,18.5], [-1.77,-43.143,17.759], color="red blue", name="Arrows_17.9939994812_5")

cluster_dict["17.9939994812"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(3.5), float(-37.5), float(17.0), float(1.0)]

cluster_dict["17.9939994812_arrows"] += cgo_arrow([3.5,-37.5,17.0], [2.681,-35.21,14.895], color="red blue", name="Arrows_17.9939994812_6")

cmd.load_cgo(cluster_dict["17.9939994812"], "Features_17.9939994812", 1)
cmd.load_cgo(cluster_dict["17.9939994812_arrows"], "Arrows_17.9939994812")
cmd.set("transparency", 0.2,"Features_17.9939994812")
cmd.group("Pharmacophore_17.9939994812", members="Features_17.9939994812")
cmd.group("Pharmacophore_17.9939994812", members="Arrows_17.9939994812")

if dirpath:
    f = join(dirpath, "1/label_threshold_17.9939994812.mol2")
else:
    f = "1/label_threshold_17.9939994812.mol2"

cmd.load(f, 'label_threshold_17.9939994812')
cmd.hide('everything', 'label_threshold_17.9939994812')
cmd.label("label_threshold_17.9939994812", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.9939994812', members= 'label_threshold_17.9939994812')


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
    f = join(dirpath, "2/label_threshold_16.6.mol2")
else:
    f = "2/label_threshold_16.6.mol2"

cmd.load(f, 'label_threshold_16.6')
cmd.hide('everything', 'label_threshold_16.6')
cmd.label("label_threshold_16.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [16.6]
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


cluster_dict = {"17.7140007019":[], "17.7140007019_arrows":[]}

cluster_dict["17.7140007019"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(-41.5), float(36.0), float(1.0)]

cluster_dict["17.7140007019_arrows"] += cgo_arrow([30.5,-41.5,36.0], [27.575,-41.176,36.727], color="blue red", name="Arrows_17.7140007019_1")

cluster_dict["17.7140007019"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(32.5), float(-35.5), float(36.0), float(1.0)]

cluster_dict["17.7140007019_arrows"] += cgo_arrow([32.5,-35.5,36.0], [34.071,-33.616,34.804], color="blue red", name="Arrows_17.7140007019_2")

cluster_dict["17.7140007019"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(33.0), float(-41.0), float(39.5), float(1.0)]

cluster_dict["17.7140007019_arrows"] += cgo_arrow([33.0,-41.0,39.5], [31.939,-40.222,42.03], color="blue red", name="Arrows_17.7140007019_3")

cluster_dict["17.7140007019"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(-32.4383212106), float(28.9691407786), float(1.0)]


cluster_dict["17.7140007019"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(32.2458735969), float(-41.650127408), float(34.0718433069), float(1.0)]


cluster_dict["17.7140007019"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(30.0), float(-44.0), float(36.5), float(1.0)]

cluster_dict["17.7140007019_arrows"] += cgo_arrow([30.0,-44.0,36.5], [28.683,-42.257,38.348], color="red blue", name="Arrows_17.7140007019_4")

cluster_dict["17.7140007019"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(-36.5), float(38.0), float(1.0)]

cluster_dict["17.7140007019_arrows"] += cgo_arrow([30.5,-36.5,38.0], [27.542,-36.538,38.091], color="red blue", name="Arrows_17.7140007019_5")

cluster_dict["17.7140007019"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(34.0), float(-39.5), float(34.5), float(1.0)]

cluster_dict["17.7140007019_arrows"] += cgo_arrow([34.0,-39.5,34.5], [33.968,-36.953,33.337], color="red blue", name="Arrows_17.7140007019_6")

cmd.load_cgo(cluster_dict["17.7140007019"], "Features_17.7140007019", 1)
cmd.load_cgo(cluster_dict["17.7140007019_arrows"], "Arrows_17.7140007019")
cmd.set("transparency", 0.2,"Features_17.7140007019")
cmd.group("Pharmacophore_17.7140007019", members="Features_17.7140007019")
cmd.group("Pharmacophore_17.7140007019", members="Arrows_17.7140007019")

if dirpath:
    f = join(dirpath, "2/label_threshold_17.7140007019.mol2")
else:
    f = "2/label_threshold_17.7140007019.mol2"

cmd.load(f, 'label_threshold_17.7140007019')
cmd.hide('everything', 'label_threshold_17.7140007019')
cmd.label("label_threshold_17.7140007019", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.7140007019', members= 'label_threshold_17.7140007019')


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
    f = join(dirpath, "3/label_threshold_16.1.mol2")
else:
    f = "3/label_threshold_16.1.mol2"

cmd.load(f, 'label_threshold_16.1')
cmd.hide('everything', 'label_threshold_16.1')
cmd.label("label_threshold_16.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [16.1]
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


cluster_dict = {"16.7420005798":[], "16.7420005798_arrows":[]}

cluster_dict["16.7420005798"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(-3.0), float(12.0), float(1.0)]

cluster_dict["16.7420005798_arrows"] += cgo_arrow([24.0,-3.0,12.0], [21.766,-1.319,13.37], color="blue red", name="Arrows_16.7420005798_1")

cluster_dict["16.7420005798"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(2.0), float(2.0), float(1.0)]

cluster_dict["16.7420005798_arrows"] += cgo_arrow([26.0,2.0,2.0], [23.655,2.241,1.129], color="blue red", name="Arrows_16.7420005798_2")

cluster_dict["16.7420005798"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(2.5), float(-1.0), float(1.0)]

cluster_dict["16.7420005798_arrows"] += cgo_arrow([26.0,2.5,-1.0], [26.236,5.121,-2.387], color="blue red", name="Arrows_16.7420005798_3")

cluster_dict["16.7420005798"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(2.0), float(-2.5), float(1.0)]

cluster_dict["16.7420005798_arrows"] += cgo_arrow([27.0,2.0,-2.5], [26.236,5.121,-2.387], color="blue red", name="Arrows_16.7420005798_4")

cluster_dict["16.7420005798"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(7.0), float(0.0), float(1.0)]

cluster_dict["16.7420005798_arrows"] += cgo_arrow([30.5,7.0,0.0], [31.908,8.68,1.718], color="blue red", name="Arrows_16.7420005798_5")

cluster_dict["16.7420005798"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(44.5), float(-0.5), float(11.5), float(1.0)]

cluster_dict["16.7420005798_arrows"] += cgo_arrow([44.5,-0.5,11.5], [47.549,0.315,10.678], color="blue red", name="Arrows_16.7420005798_6")

cluster_dict["16.7420005798"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(29.4331221431), float(1.09660523711), float(1.11365114511), float(1.0)]


cluster_dict["16.7420005798"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(34.0345533343), float(1.77668193439), float(11.979851349), float(1.0)]


cluster_dict["16.7420005798"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(34.9028224901), float(2.39631919966), float(8.28702885726), float(1.0)]


cluster_dict["16.7420005798"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(41.4901511383), float(2.63355613373), float(-0.683013480244), float(1.0)]


cluster_dict["16.7420005798"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(43.2834367478), float(-4.72043011736), float(10.5970722458), float(1.0)]


cluster_dict["16.7420005798"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(44.5), float(4.0), float(8.5), float(1.0)]


cluster_dict["16.7420005798"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(-2.5), float(11.5), float(1.0)]

cluster_dict["16.7420005798_arrows"] += cgo_arrow([27.0,-2.5,11.5], [25.888,-5.209,9.814], color="red blue", name="Arrows_16.7420005798_7")

cluster_dict["16.7420005798"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(43.0), float(-7.0), float(8.5), float(1.0)]

cluster_dict["16.7420005798_arrows"] += cgo_arrow([43.0,-7.0,8.5], [40.826,-8.37,6.871], color="red blue", name="Arrows_16.7420005798_8")

cluster_dict["16.7420005798"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(44.5), float(-3.5), float(10.0), float(1.0)]

cluster_dict["16.7420005798_arrows"] += cgo_arrow([44.5,-3.5,10.0], [45.601,-1.014,8.286], color="red blue", name="Arrows_16.7420005798_9")

cmd.load_cgo(cluster_dict["16.7420005798"], "Features_16.7420005798", 1)
cmd.load_cgo(cluster_dict["16.7420005798_arrows"], "Arrows_16.7420005798")
cmd.set("transparency", 0.2,"Features_16.7420005798")
cmd.group("Pharmacophore_16.7420005798", members="Features_16.7420005798")
cmd.group("Pharmacophore_16.7420005798", members="Arrows_16.7420005798")

if dirpath:
    f = join(dirpath, "3/label_threshold_16.7420005798.mol2")
else:
    f = "3/label_threshold_16.7420005798.mol2"

cmd.load(f, 'label_threshold_16.7420005798')
cmd.hide('everything', 'label_threshold_16.7420005798')
cmd.label("label_threshold_16.7420005798", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.7420005798', members= 'label_threshold_16.7420005798')


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
    f = join(dirpath, "4/label_threshold_14.9.mol2")
else:
    f = "4/label_threshold_14.9.mol2"

cmd.load(f, 'label_threshold_14.9')
cmd.hide('everything', 'label_threshold_14.9')
cmd.label("label_threshold_14.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [14.9]
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


cluster_dict = {"16.5830001831":[], "16.5830001831_arrows":[]}

cluster_dict["16.5830001831"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(-31.0), float(22.0), float(1.0)]

cluster_dict["16.5830001831_arrows"] += cgo_arrow([17.5,-31.0,22.0], [15.264,-32.638,22.22], color="blue red", name="Arrows_16.5830001831_1")

cluster_dict["16.5830001831"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(-44.0), float(26.5), float(1.0)]

cluster_dict["16.5830001831_arrows"] += cgo_arrow([18.5,-44.0,26.5], [15.475,-45.652,25.95], color="blue red", name="Arrows_16.5830001831_2")

cluster_dict["16.5830001831"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(20.0), float(-33.5), float(26.0), float(1.0)]

cluster_dict["16.5830001831_arrows"] += cgo_arrow([20.0,-33.5,26.0], [20.162,-31.121,28.672], color="blue red", name="Arrows_16.5830001831_3")

cluster_dict["16.5830001831"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(12.7976276996), float(-35.398043533), float(32.7976276996), float(1.0)]


cluster_dict["16.5830001831"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(13.1689145848), float(-35.5211648263), float(33.0), float(1.0)]


cluster_dict["16.5830001831"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(19.6089477009), float(-31.6599193519), float(21.875703526), float(1.0)]


cluster_dict["16.5830001831"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(19.2724252792), float(-41.6899719042), float(28.1003851422), float(1.0)]


cluster_dict["16.5830001831"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(18.5625), float(-45.0), float(24.9375), float(1.0)]


cluster_dict["16.5830001831"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(-28.0), float(13.5), float(1.0)]


cluster_dict["16.5830001831"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(-32.4383212106), float(28.9691407786), float(1.0)]


cluster_dict["16.5830001831"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(-40.5), float(28.5), float(1.0)]

cluster_dict["16.5830001831_arrows"] += cgo_arrow([16.0,-40.5,28.5], [13.403,-39.381,29.79], color="red blue", name="Arrows_16.5830001831_4")

cluster_dict["16.5830001831"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(-29.5), float(22.0), float(1.0)]

cluster_dict["16.5830001831_arrows"] += cgo_arrow([16.0,-29.5,22.0], [13.97,-30.815,22.039], color="red blue", name="Arrows_16.5830001831_5")

cluster_dict["16.5830001831"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(-31.5), float(24.0), float(1.0)]

cluster_dict["16.5830001831_arrows"] += cgo_arrow([17.5,-31.5,24.0], [15.279,-30.499,26.119], color="red blue", name="Arrows_16.5830001831_6")

cluster_dict["16.5830001831"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(-32.0), float(17.0), float(1.0)]

cluster_dict["16.5830001831_arrows"] += cgo_arrow([18.5,-32.0,17.0], [20.191,-30.573,14.342], color="red blue", name="Arrows_16.5830001831_7")

cluster_dict["16.5830001831"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.5), float(-44.0), float(23.5), float(1.0)]

cluster_dict["16.5830001831_arrows"] += cgo_arrow([20.5,-44.0,23.5], [22.794,-42.19,22.011], color="red blue", name="Arrows_16.5830001831_8")

cluster_dict["16.5830001831"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(-31.0), float(24.0), float(1.0)]

cluster_dict["16.5830001831_arrows"] += cgo_arrow([23.0,-31.0,24.0], [24.932,-30.214,21.553], color="red blue", name="Arrows_16.5830001831_9")

cmd.load_cgo(cluster_dict["16.5830001831"], "Features_16.5830001831", 1)
cmd.load_cgo(cluster_dict["16.5830001831_arrows"], "Arrows_16.5830001831")
cmd.set("transparency", 0.2,"Features_16.5830001831")
cmd.group("Pharmacophore_16.5830001831", members="Features_16.5830001831")
cmd.group("Pharmacophore_16.5830001831", members="Arrows_16.5830001831")

if dirpath:
    f = join(dirpath, "4/label_threshold_16.5830001831.mol2")
else:
    f = "4/label_threshold_16.5830001831.mol2"

cmd.load(f, 'label_threshold_16.5830001831')
cmd.hide('everything', 'label_threshold_16.5830001831')
cmd.label("label_threshold_16.5830001831", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.5830001831', members= 'label_threshold_16.5830001831')


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
    f = join(dirpath, "5/label_threshold_15.1.mol2")
else:
    f = "5/label_threshold_15.1.mol2"

cmd.load(f, 'label_threshold_15.1')
cmd.hide('everything', 'label_threshold_15.1')
cmd.label("label_threshold_15.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.1]
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


cluster_dict = {"16.4230003357":[], "16.4230003357_arrows":[]}

cluster_dict["16.4230003357"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(-31.0), float(22.0), float(1.0)]

cluster_dict["16.4230003357_arrows"] += cgo_arrow([17.5,-31.0,22.0], [15.264,-32.638,22.22], color="blue red", name="Arrows_16.4230003357_1")

cluster_dict["16.4230003357"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(-44.0), float(26.5), float(1.0)]

cluster_dict["16.4230003357_arrows"] += cgo_arrow([18.5,-44.0,26.5], [15.475,-45.652,25.95], color="blue red", name="Arrows_16.4230003357_2")

cluster_dict["16.4230003357"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(20.0), float(-33.5), float(26.0), float(1.0)]

cluster_dict["16.4230003357_arrows"] += cgo_arrow([20.0,-33.5,26.0], [20.162,-31.121,28.672], color="blue red", name="Arrows_16.4230003357_3")

cluster_dict["16.4230003357"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.5), float(-30.5), float(21.0), float(1.0)]

cluster_dict["16.4230003357_arrows"] += cgo_arrow([19.5,-30.5,21.0], [22.725,-29.591,19.691], color="blue red", name="Arrows_16.4230003357_4")

cluster_dict["16.4230003357"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(12.1383293424), float(-37.0312123023), float(33.8141889245), float(1.0)]


cluster_dict["16.4230003357"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(19.3261786062), float(-41.7096184277), float(28.0941020175), float(1.0)]


cluster_dict["16.4230003357"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(20.3694078785), float(-32.7798788655), float(22.9498504669), float(1.0)]


cluster_dict["16.4230003357"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(20.1801188625), float(-47.8549009041), float(25.6989864157), float(1.0)]


cluster_dict["16.4230003357"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(21.2), float(-32.8), float(27.0), float(1.0)]


cluster_dict["16.4230003357"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(-32.9084143942), float(28.9471952019), float(1.0)]


cluster_dict["16.4230003357"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(28.5), float(-44.8426303555), float(30.194764049), float(1.0)]


cluster_dict["16.4230003357"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(28.5), float(-39.4995775304), float(34.7862616234), float(1.0)]


cluster_dict["16.4230003357"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(-40.5), float(28.5), float(1.0)]

cluster_dict["16.4230003357_arrows"] += cgo_arrow([16.0,-40.5,28.5], [13.403,-39.381,29.79], color="red blue", name="Arrows_16.4230003357_5")

cluster_dict["16.4230003357"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(-30.0), float(22.0), float(1.0)]

cluster_dict["16.4230003357_arrows"] += cgo_arrow([16.0,-30.0,22.0], [13.97,-30.815,22.039], color="red blue", name="Arrows_16.4230003357_6")

cluster_dict["16.4230003357"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(-31.5), float(24.0), float(1.0)]

cluster_dict["16.4230003357_arrows"] += cgo_arrow([17.5,-31.5,24.0], [15.279,-30.499,26.119], color="red blue", name="Arrows_16.4230003357_7")

cluster_dict["16.4230003357"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.5), float(-44.0), float(23.5), float(1.0)]

cluster_dict["16.4230003357_arrows"] += cgo_arrow([20.5,-44.0,23.5], [22.794,-42.19,22.011], color="red blue", name="Arrows_16.4230003357_8")

cluster_dict["16.4230003357"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(-31.0), float(24.0), float(1.0)]

cluster_dict["16.4230003357_arrows"] += cgo_arrow([23.0,-31.0,24.0], [24.932,-30.214,21.553], color="red blue", name="Arrows_16.4230003357_9")

cmd.load_cgo(cluster_dict["16.4230003357"], "Features_16.4230003357", 1)
cmd.load_cgo(cluster_dict["16.4230003357_arrows"], "Arrows_16.4230003357")
cmd.set("transparency", 0.2,"Features_16.4230003357")
cmd.group("Pharmacophore_16.4230003357", members="Features_16.4230003357")
cmd.group("Pharmacophore_16.4230003357", members="Arrows_16.4230003357")

if dirpath:
    f = join(dirpath, "5/label_threshold_16.4230003357.mol2")
else:
    f = "5/label_threshold_16.4230003357.mol2"

cmd.load(f, 'label_threshold_16.4230003357')
cmd.hide('everything', 'label_threshold_16.4230003357')
cmd.label("label_threshold_16.4230003357", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.4230003357', members= 'label_threshold_16.4230003357')


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
    f = join(dirpath, "6/label_threshold_13.4.mol2")
else:
    f = "6/label_threshold_13.4.mol2"

cmd.load(f, 'label_threshold_13.4')
cmd.hide('everything', 'label_threshold_13.4')
cmd.label("label_threshold_13.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.4]
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


cluster_dict = {"16.1469993591":[], "16.1469993591_arrows":[]}

cluster_dict["16.1469993591"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(2.0), float(-2.5), float(1.0)]

cluster_dict["16.1469993591_arrows"] += cgo_arrow([27.0,2.0,-2.5], [26.236,5.121,-2.387], color="blue red", name="Arrows_16.1469993591_1")

cluster_dict["16.1469993591"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(7.0), float(0.0), float(1.0)]

cluster_dict["16.1469993591_arrows"] += cgo_arrow([30.5,7.0,0.0], [31.908,8.68,1.718], color="blue red", name="Arrows_16.1469993591_2")

cluster_dict["16.1469993591"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(29.8402203895), float(1.24572916546), float(0.684319074505), float(1.0)]


cluster_dict["16.1469993591"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(34.8081563576), float(2.19319259381), float(8.38723598175), float(1.0)]


cluster_dict["16.1469993591"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(41.6672174937), float(2.67726504769), float(-0.764066775479), float(1.0)]


cluster_dict["16.1469993591"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(-2.5), float(4.0), float(1.0)]

cluster_dict["16.1469993591_arrows"] += cgo_arrow([26.0,-2.5,4.0], [26.993,-4.956,2.181], color="red blue", name="Arrows_16.1469993591_3")

cluster_dict["16.1469993591"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(31.5), float(-0.5), float(2.0), float(1.0)]

cluster_dict["16.1469993591_arrows"] += cgo_arrow([31.5,-0.5,2.0], [33.219,-2.042,2.779], color="red blue", name="Arrows_16.1469993591_4")

cmd.load_cgo(cluster_dict["16.1469993591"], "Features_16.1469993591", 1)
cmd.load_cgo(cluster_dict["16.1469993591_arrows"], "Arrows_16.1469993591")
cmd.set("transparency", 0.2,"Features_16.1469993591")
cmd.group("Pharmacophore_16.1469993591", members="Features_16.1469993591")
cmd.group("Pharmacophore_16.1469993591", members="Arrows_16.1469993591")

if dirpath:
    f = join(dirpath, "6/label_threshold_16.1469993591.mol2")
else:
    f = "6/label_threshold_16.1469993591.mol2"

cmd.load(f, 'label_threshold_16.1469993591')
cmd.hide('everything', 'label_threshold_16.1469993591')
cmd.label("label_threshold_16.1469993591", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.1469993591', members= 'label_threshold_16.1469993591')


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
    f = join(dirpath, "7/label_threshold_10.7.mol2")
else:
    f = "7/label_threshold_10.7.mol2"

cmd.load(f, 'label_threshold_10.7')
cmd.hide('everything', 'label_threshold_10.7')
cmd.label("label_threshold_10.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.7]
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


cluster_dict = {"16.1259994507":[], "16.1259994507_arrows":[]}

cluster_dict["16.1259994507"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(0.5), float(-37.0), float(18.0), float(1.0)]

cluster_dict["16.1259994507_arrows"] += cgo_arrow([0.5,-37.0,18.0], [0.174,-36.164,15.33], color="blue red", name="Arrows_16.1259994507_1")

cluster_dict["16.1259994507"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(5.0), float(-37.0), float(18.5), float(1.0)]

cluster_dict["16.1259994507_arrows"] += cgo_arrow([5.0,-37.0,18.5], [6.939,-35.33,17.433], color="blue red", name="Arrows_16.1259994507_2")

cluster_dict["16.1259994507"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(2.36139853989), float(-37.2744824606), float(20.3749822028), float(1.0)]


cluster_dict["16.1259994507"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(14.9584185084), float(-44.6313431733), float(18.3150330733), float(1.0)]


cluster_dict["16.1259994507"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-1.0), float(-37.0), float(20.5), float(1.0)]


cluster_dict["16.1259994507"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-1.0), float(-39.0), float(18.5), float(1.0)]

cluster_dict["16.1259994507_arrows"] += cgo_arrow([-1.0,-39.0,18.5], [-1.77,-43.143,17.759], color="red blue", name="Arrows_16.1259994507_3")

cluster_dict["16.1259994507"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(3.5), float(-37.5), float(17.0), float(1.0)]

cluster_dict["16.1259994507_arrows"] += cgo_arrow([3.5,-37.5,17.0], [2.681,-35.21,14.895], color="red blue", name="Arrows_16.1259994507_4")

cmd.load_cgo(cluster_dict["16.1259994507"], "Features_16.1259994507", 1)
cmd.load_cgo(cluster_dict["16.1259994507_arrows"], "Arrows_16.1259994507")
cmd.set("transparency", 0.2,"Features_16.1259994507")
cmd.group("Pharmacophore_16.1259994507", members="Features_16.1259994507")
cmd.group("Pharmacophore_16.1259994507", members="Arrows_16.1259994507")

if dirpath:
    f = join(dirpath, "7/label_threshold_16.1259994507.mol2")
else:
    f = "7/label_threshold_16.1259994507.mol2"

cmd.load(f, 'label_threshold_16.1259994507')
cmd.hide('everything', 'label_threshold_16.1259994507')
cmd.label("label_threshold_16.1259994507", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.1259994507', members= 'label_threshold_16.1259994507')


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
    f = join(dirpath, "8/label_threshold_13.4.mol2")
else:
    f = "8/label_threshold_13.4.mol2"

cmd.load(f, 'label_threshold_13.4')
cmd.hide('everything', 'label_threshold_13.4')
cmd.label("label_threshold_13.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.4]
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


cluster_dict = {"15.8090000153":[], "15.8090000153_arrows":[]}

cluster_dict["15.8090000153"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(-3.0), float(12.0), float(1.0)]

cluster_dict["15.8090000153_arrows"] += cgo_arrow([24.0,-3.0,12.0], [21.766,-1.319,13.37], color="blue red", name="Arrows_15.8090000153_1")

cluster_dict["15.8090000153"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(2.0), float(-2.5), float(1.0)]

cluster_dict["15.8090000153_arrows"] += cgo_arrow([27.0,2.0,-2.5], [26.236,5.121,-2.387], color="blue red", name="Arrows_15.8090000153_2")

cluster_dict["15.8090000153"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(7.0), float(0.0), float(1.0)]

cluster_dict["15.8090000153_arrows"] += cgo_arrow([30.5,7.0,0.0], [31.908,8.68,1.718], color="blue red", name="Arrows_15.8090000153_3")

cluster_dict["15.8090000153"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(29.8360627635), float(1.24309894995), float(0.688932158277), float(1.0)]


cluster_dict["15.8090000153"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(33.7801552175), float(1.70494645616), float(10.8799148095), float(1.0)]


cluster_dict["15.8090000153"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(40.531605642), float(1.94987620702), float(-0.753589955375), float(1.0)]


cluster_dict["15.8090000153"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(40.5), float(2.72222222222), float(-4.38888888889), float(1.0)]


cluster_dict["15.8090000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(-7.5), float(6.5), float(1.0)]

cluster_dict["15.8090000153_arrows"] += cgo_arrow([22.5,-7.5,6.5], [22.101,-6.363,4.763], color="red blue", name="Arrows_15.8090000153_4")

cluster_dict["15.8090000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(-2.5), float(11.5), float(1.0)]

cluster_dict["15.8090000153_arrows"] += cgo_arrow([27.0,-2.5,11.5], [25.888,-5.209,9.814], color="red blue", name="Arrows_15.8090000153_5")

cluster_dict["15.8090000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(-2.5), float(4.0), float(1.0)]

cluster_dict["15.8090000153_arrows"] += cgo_arrow([26.0,-2.5,4.0], [26.993,-4.956,2.181], color="red blue", name="Arrows_15.8090000153_6")

cluster_dict["15.8090000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(31.5), float(-0.5), float(2.0), float(1.0)]

cluster_dict["15.8090000153_arrows"] += cgo_arrow([31.5,-0.5,2.0], [33.219,-2.042,2.779], color="red blue", name="Arrows_15.8090000153_7")

cmd.load_cgo(cluster_dict["15.8090000153"], "Features_15.8090000153", 1)
cmd.load_cgo(cluster_dict["15.8090000153_arrows"], "Arrows_15.8090000153")
cmd.set("transparency", 0.2,"Features_15.8090000153")
cmd.group("Pharmacophore_15.8090000153", members="Features_15.8090000153")
cmd.group("Pharmacophore_15.8090000153", members="Arrows_15.8090000153")

if dirpath:
    f = join(dirpath, "8/label_threshold_15.8090000153.mol2")
else:
    f = "8/label_threshold_15.8090000153.mol2"

cmd.load(f, 'label_threshold_15.8090000153')
cmd.hide('everything', 'label_threshold_15.8090000153')
cmd.label("label_threshold_15.8090000153", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.8090000153', members= 'label_threshold_15.8090000153')


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
    f = join(dirpath, "9/label_threshold_13.8.mol2")
else:
    f = "9/label_threshold_13.8.mol2"

cmd.load(f, 'label_threshold_13.8')
cmd.hide('everything', 'label_threshold_13.8')
cmd.label("label_threshold_13.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.8]
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


cluster_dict = {"15.8039999008":[], "15.8039999008_arrows":[]}

cluster_dict["15.8039999008"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(-3.0), float(12.5), float(1.0)]

cluster_dict["15.8039999008_arrows"] += cgo_arrow([24.0,-3.0,12.5], [21.766,-1.319,13.37], color="blue red", name="Arrows_15.8039999008_1")

cluster_dict["15.8039999008"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(-11.0), float(9.0), float(1.0)]

cluster_dict["15.8039999008_arrows"] += cgo_arrow([24.5,-11.0,9.0], [26.667,-9.958,10.507], color="blue red", name="Arrows_15.8039999008_2")

cluster_dict["15.8039999008"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(2.0), float(-2.5), float(1.0)]

cluster_dict["15.8039999008_arrows"] += cgo_arrow([27.0,2.0,-2.5], [26.236,5.121,-2.387], color="blue red", name="Arrows_15.8039999008_3")

cluster_dict["15.8039999008"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(-9.5), float(14.0), float(1.0)]

cluster_dict["15.8039999008_arrows"] += cgo_arrow([26.5,-9.5,14.0], [28.783,-7.532,13.582], color="blue red", name="Arrows_15.8039999008_4")

cluster_dict["15.8039999008"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(7.0), float(0.0), float(1.0)]

cluster_dict["15.8039999008_arrows"] += cgo_arrow([30.5,7.0,0.0], [31.908,8.68,1.718], color="blue red", name="Arrows_15.8039999008_5")

cluster_dict["15.8039999008"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.9206195707), float(-11.0512520991), float(8.84098633983), float(1.0)]


cluster_dict["15.8039999008"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(29.8236361515), float(1.21709503082), float(0.723811118156), float(1.0)]


cluster_dict["15.8039999008"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(33.9630721097), float(1.79129273337), float(11.884841271), float(1.0)]


cluster_dict["15.8039999008"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(39.5), float(-1.7633203123), float(-0.560385400252), float(1.0)]


cluster_dict["15.8039999008"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(-8.5), float(7.0), float(1.0)]

cluster_dict["15.8039999008_arrows"] += cgo_arrow([22.5,-8.5,7.0], [22.101,-6.363,4.763], color="red blue", name="Arrows_15.8039999008_6")

cluster_dict["15.8039999008"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(-2.5), float(11.5), float(1.0)]

cluster_dict["15.8039999008_arrows"] += cgo_arrow([27.0,-2.5,11.5], [25.888,-5.209,9.814], color="red blue", name="Arrows_15.8039999008_7")

cluster_dict["15.8039999008"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(-2.5), float(4.0), float(1.0)]

cluster_dict["15.8039999008_arrows"] += cgo_arrow([26.0,-2.5,4.0], [26.993,-4.956,2.181], color="red blue", name="Arrows_15.8039999008_8")

cluster_dict["15.8039999008"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(31.5), float(-0.5), float(2.0), float(1.0)]

cluster_dict["15.8039999008_arrows"] += cgo_arrow([31.5,-0.5,2.0], [33.219,-2.042,2.779], color="red blue", name="Arrows_15.8039999008_9")

cmd.load_cgo(cluster_dict["15.8039999008"], "Features_15.8039999008", 1)
cmd.load_cgo(cluster_dict["15.8039999008_arrows"], "Arrows_15.8039999008")
cmd.set("transparency", 0.2,"Features_15.8039999008")
cmd.group("Pharmacophore_15.8039999008", members="Features_15.8039999008")
cmd.group("Pharmacophore_15.8039999008", members="Arrows_15.8039999008")

if dirpath:
    f = join(dirpath, "9/label_threshold_15.8039999008.mol2")
else:
    f = "9/label_threshold_15.8039999008.mol2"

cmd.load(f, 'label_threshold_15.8039999008')
cmd.hide('everything', 'label_threshold_15.8039999008')
cmd.label("label_threshold_15.8039999008", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.8039999008', members= 'label_threshold_15.8039999008')


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
    f = join(dirpath, "10/label_threshold_14.8.mol2")
else:
    f = "10/label_threshold_14.8.mol2"

cmd.load(f, 'label_threshold_14.8')
cmd.hide('everything', 'label_threshold_14.8')
cmd.label("label_threshold_14.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [14.8]
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


cluster_dict = {"15.7994999886":[], "15.7994999886_arrows":[]}

cluster_dict["15.7994999886"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(43.5), float(-8.5), float(9.5), float(1.0)]

cluster_dict["15.7994999886_arrows"] += cgo_arrow([43.5,-8.5,9.5], [43.181,-9.676,7.212], color="blue red", name="Arrows_15.7994999886_1")

cluster_dict["15.7994999886"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(44.5), float(-0.5), float(11.5), float(1.0)]

cluster_dict["15.7994999886_arrows"] += cgo_arrow([44.5,-0.5,11.5], [47.549,0.315,10.678], color="blue red", name="Arrows_15.7994999886_2")

cluster_dict["15.7994999886"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(48.0), float(-6.0), float(8.0), float(1.0)]

cluster_dict["15.7994999886_arrows"] += cgo_arrow([48.0,-6.0,8.0], [45.592,-6.463,6.489], color="blue red", name="Arrows_15.7994999886_3")

cluster_dict["15.7994999886"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(52.5), float(-7.0), float(6.5), float(1.0)]

cluster_dict["15.7994999886_arrows"] += cgo_arrow([52.5,-7.0,6.5], [50.139,-8.382,8.294], color="blue red", name="Arrows_15.7994999886_4")

cluster_dict["15.7994999886"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(42.0363852893), float(2.50855217932), float(0.800768802169), float(1.0)]


cluster_dict["15.7994999886"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(43.9230041039), float(-5.99446548399), float(11.0629519479), float(1.0)]


cluster_dict["15.7994999886"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(44.5026696863), float(3.84705930573), float(8.5), float(1.0)]


cluster_dict["15.7994999886"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(51.1316842148), float(-8.25151219982), float(4.11539457304), float(1.0)]


cluster_dict["15.7994999886"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(53.0), float(-7.0), float(0.5), float(1.0)]


cluster_dict["15.7994999886"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(47.0), float(-7.0), float(9.0), float(1.0)]

cluster_dict["15.7994999886_arrows"] += cgo_arrow([47.0,-7.0,9.0], [45.592,-6.463,6.489], color="red blue", name="Arrows_15.7994999886_5")

cluster_dict["15.7994999886"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(45.0), float(-3.5), float(10.0), float(1.0)]

cluster_dict["15.7994999886_arrows"] += cgo_arrow([45.0,-3.5,10.0], [45.601,-1.014,8.286], color="red blue", name="Arrows_15.7994999886_6")

cluster_dict["15.7994999886"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(48.0), float(-7.0), float(5.5), float(1.0)]

cluster_dict["15.7994999886_arrows"] += cgo_arrow([48.0,-7.0,5.5], [45.592,-6.463,6.489], color="red blue", name="Arrows_15.7994999886_7")

cluster_dict["15.7994999886"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(51.0), float(-9.0), float(3.5), float(1.0)]

cluster_dict["15.7994999886_arrows"] += cgo_arrow([51.0,-9.0,3.5], [48.96,-11.301,1.715], color="red blue", name="Arrows_15.7994999886_8")

cluster_dict["15.7994999886"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(51.0), float(-9.0), float(3.5), float(1.0)]

cluster_dict["15.7994999886_arrows"] += cgo_arrow([51.0,-9.0,3.5], [48.96,-11.301,1.715], color="red blue", name="Arrows_15.7994999886_9")

cmd.load_cgo(cluster_dict["15.7994999886"], "Features_15.7994999886", 1)
cmd.load_cgo(cluster_dict["15.7994999886_arrows"], "Arrows_15.7994999886")
cmd.set("transparency", 0.2,"Features_15.7994999886")
cmd.group("Pharmacophore_15.7994999886", members="Features_15.7994999886")
cmd.group("Pharmacophore_15.7994999886", members="Arrows_15.7994999886")

if dirpath:
    f = join(dirpath, "10/label_threshold_15.7994999886.mol2")
else:
    f = "10/label_threshold_15.7994999886.mol2"

cmd.load(f, 'label_threshold_15.7994999886')
cmd.hide('everything', 'label_threshold_15.7994999886')
cmd.label("label_threshold_15.7994999886", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.7994999886', members= 'label_threshold_15.7994999886')


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
    f = join(dirpath, "11/label_threshold_14.0.mol2")
else:
    f = "11/label_threshold_14.0.mol2"

cmd.load(f, 'label_threshold_14.0')
cmd.hide('everything', 'label_threshold_14.0')
cmd.label("label_threshold_14.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [14.0]
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


cluster_dict = {"15.7620000839":[], "15.7620000839_arrows":[]}

cluster_dict["15.7620000839"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(16.5), float(-35.5), float(22.0), float(1.0)]

cluster_dict["15.7620000839_arrows"] += cgo_arrow([16.5,-35.5,22.0], [14.405,-37.587,22.27], color="blue red", name="Arrows_15.7620000839_1")

cluster_dict["15.7620000839"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(-44.0), float(26.0), float(1.0)]

cluster_dict["15.7620000839_arrows"] += cgo_arrow([18.5,-44.0,26.0], [15.475,-45.652,25.95], color="blue red", name="Arrows_15.7620000839_2")

cluster_dict["15.7620000839"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(-44.0), float(26.0), float(1.0)]

cluster_dict["15.7620000839_arrows"] += cgo_arrow([18.5,-44.0,26.0], [15.475,-45.652,25.95], color="blue red", name="Arrows_15.7620000839_3")

cluster_dict["15.7620000839"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.5), float(-34.5), float(26.5), float(1.0)]

cluster_dict["15.7620000839_arrows"] += cgo_arrow([19.5,-34.5,26.5], [19.257,-35.376,29.708], color="blue red", name="Arrows_15.7620000839_4")

cluster_dict["15.7620000839"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(6.09902086484), float(-38.2906550109), float(19.7127593343), float(1.0)]


cluster_dict["15.7620000839"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(16.3114139499), float(-44.4330254646), float(21.3558913511), float(1.0)]


cluster_dict["15.7620000839"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(19.7171878551), float(-44.5390458116), float(26.7535597902), float(1.0)]


cluster_dict["15.7620000839"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(19.9734023839), float(-34.092715892), float(22.3072385808), float(1.0)]


cluster_dict["15.7620000839"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(-40.5), float(28.5), float(1.0)]

cluster_dict["15.7620000839_arrows"] += cgo_arrow([16.0,-40.5,28.5], [13.403,-39.381,29.79], color="red blue", name="Arrows_15.7620000839_5")

cluster_dict["15.7620000839"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.5), float(-35.5), float(24.0), float(1.0)]

cluster_dict["15.7620000839_arrows"] += cgo_arrow([20.5,-35.5,24.0], [21.795,-37.163,26.941], color="red blue", name="Arrows_15.7620000839_6")

cluster_dict["15.7620000839"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(-34.5), float(25.5), float(1.0)]

cluster_dict["15.7620000839_arrows"] += cgo_arrow([18.5,-34.5,25.5], [16.555,-32.827,27.873], color="red blue", name="Arrows_15.7620000839_7")

cluster_dict["15.7620000839"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(-32.5), float(17.0), float(1.0)]

cluster_dict["15.7620000839_arrows"] += cgo_arrow([18.0,-32.5,17.0], [17.274,-35.929,18.15], color="red blue", name="Arrows_15.7620000839_8")

cluster_dict["15.7620000839"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.5), float(-44.0), float(23.5), float(1.0)]

cluster_dict["15.7620000839_arrows"] += cgo_arrow([20.5,-44.0,23.5], [22.794,-42.19,22.011], color="red blue", name="Arrows_15.7620000839_9")

cluster_dict["15.7620000839"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(-34.5), float(22.5), float(1.0)]

cluster_dict["15.7620000839_arrows"] += cgo_arrow([22.5,-34.5,22.5], [24.521,-37.391,22.877], color="red blue", name="Arrows_15.7620000839_10")

cmd.load_cgo(cluster_dict["15.7620000839"], "Features_15.7620000839", 1)
cmd.load_cgo(cluster_dict["15.7620000839_arrows"], "Arrows_15.7620000839")
cmd.set("transparency", 0.2,"Features_15.7620000839")
cmd.group("Pharmacophore_15.7620000839", members="Features_15.7620000839")
cmd.group("Pharmacophore_15.7620000839", members="Arrows_15.7620000839")

if dirpath:
    f = join(dirpath, "11/label_threshold_15.7620000839.mol2")
else:
    f = "11/label_threshold_15.7620000839.mol2"

cmd.load(f, 'label_threshold_15.7620000839')
cmd.hide('everything', 'label_threshold_15.7620000839')
cmd.label("label_threshold_15.7620000839", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.7620000839', members= 'label_threshold_15.7620000839')


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
    f = join(dirpath, "12/label_threshold_13.0.mol2")
else:
    f = "12/label_threshold_13.0.mol2"

cmd.load(f, 'label_threshold_13.0')
cmd.hide('everything', 'label_threshold_13.0')
cmd.label("label_threshold_13.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.0]
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


cluster_dict = {"15.1840000153":[], "15.1840000153_arrows":[]}

cluster_dict["15.1840000153"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(43.5), float(-8.5), float(9.5), float(1.0)]

cluster_dict["15.1840000153_arrows"] += cgo_arrow([43.5,-8.5,9.5], [43.181,-9.676,7.212], color="blue red", name="Arrows_15.1840000153_1")

cluster_dict["15.1840000153"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(44.5), float(-1.0), float(11.0), float(1.0)]

cluster_dict["15.1840000153_arrows"] += cgo_arrow([44.5,-1.0,11.0], [47.549,0.315,10.678], color="blue red", name="Arrows_15.1840000153_2")

cluster_dict["15.1840000153"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(48.0), float(-6.0), float(8.0), float(1.0)]

cluster_dict["15.1840000153_arrows"] += cgo_arrow([48.0,-6.0,8.0], [45.592,-6.463,6.489], color="blue red", name="Arrows_15.1840000153_3")

cluster_dict["15.1840000153"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(52.5), float(-7.0), float(6.5), float(1.0)]

cluster_dict["15.1840000153_arrows"] += cgo_arrow([52.5,-7.0,6.5], [50.139,-8.382,8.294], color="blue red", name="Arrows_15.1840000153_4")

cluster_dict["15.1840000153"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(41.7715163298), float(-0.616554517252), float(-1.12806149459), float(1.0)]


cluster_dict["15.1840000153"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(44.1433299831), float(-6.15080711521), float(10.7317697741), float(1.0)]


cluster_dict["15.1840000153"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(51.1627419202), float(-8.22607166055), float(4.17812572322), float(1.0)]


cluster_dict["15.1840000153"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(53.0), float(-7.0), float(0.5), float(1.0)]


cluster_dict["15.1840000153"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(54.0), float(-6.99755319129), float(1.5), float(1.0)]


cluster_dict["15.1840000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(47.0), float(-7.0), float(9.0), float(1.0)]

cluster_dict["15.1840000153_arrows"] += cgo_arrow([47.0,-7.0,9.0], [45.592,-6.463,6.489], color="red blue", name="Arrows_15.1840000153_5")

cluster_dict["15.1840000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(45.0), float(-3.5), float(10.0), float(1.0)]

cluster_dict["15.1840000153_arrows"] += cgo_arrow([45.0,-3.5,10.0], [45.601,-1.014,8.286], color="red blue", name="Arrows_15.1840000153_6")

cluster_dict["15.1840000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(48.0), float(-7.0), float(5.5), float(1.0)]

cluster_dict["15.1840000153_arrows"] += cgo_arrow([48.0,-7.0,5.5], [45.592,-6.463,6.489], color="red blue", name="Arrows_15.1840000153_7")

cluster_dict["15.1840000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(51.0), float(-9.0), float(3.5), float(1.0)]

cluster_dict["15.1840000153_arrows"] += cgo_arrow([51.0,-9.0,3.5], [48.96,-11.301,1.715], color="red blue", name="Arrows_15.1840000153_8")

cluster_dict["15.1840000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(51.0), float(-9.0), float(3.5), float(1.0)]

cluster_dict["15.1840000153_arrows"] += cgo_arrow([51.0,-9.0,3.5], [48.96,-11.301,1.715], color="red blue", name="Arrows_15.1840000153_9")

cmd.load_cgo(cluster_dict["15.1840000153"], "Features_15.1840000153", 1)
cmd.load_cgo(cluster_dict["15.1840000153_arrows"], "Arrows_15.1840000153")
cmd.set("transparency", 0.2,"Features_15.1840000153")
cmd.group("Pharmacophore_15.1840000153", members="Features_15.1840000153")
cmd.group("Pharmacophore_15.1840000153", members="Arrows_15.1840000153")

if dirpath:
    f = join(dirpath, "12/label_threshold_15.1840000153.mol2")
else:
    f = "12/label_threshold_15.1840000153.mol2"

cmd.load(f, 'label_threshold_15.1840000153')
cmd.hide('everything', 'label_threshold_15.1840000153')
cmd.label("label_threshold_15.1840000153", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.1840000153', members= 'label_threshold_15.1840000153')


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
    f = join(dirpath, "13/label_threshold_12.6.mol2")
else:
    f = "13/label_threshold_12.6.mol2"

cmd.load(f, 'label_threshold_12.6')
cmd.hide('everything', 'label_threshold_12.6')
cmd.label("label_threshold_12.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.6]
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


cluster_dict = {"15.1260004044":[], "15.1260004044_arrows":[]}

cluster_dict["15.1260004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(-43.5), float(3.5), float(1.0)]

cluster_dict["15.1260004044_arrows"] += cgo_arrow([29.5,-43.5,3.5], [29.393,-40.692,4.576], color="blue red", name="Arrows_15.1260004044_1")

cluster_dict["15.1260004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(36.0), float(-35.5), float(-1.0), float(1.0)]

cluster_dict["15.1260004044_arrows"] += cgo_arrow([36.0,-35.5,-1.0], [35.162,-31.753,-0.868], color="blue red", name="Arrows_15.1260004044_2")

cluster_dict["15.1260004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(36.0), float(-35.5), float(-1.0), float(1.0)]

cluster_dict["15.1260004044_arrows"] += cgo_arrow([36.0,-35.5,-1.0], [35.162,-31.753,-0.868], color="blue red", name="Arrows_15.1260004044_3")

cluster_dict["15.1260004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(38.5), float(-25.0), float(0.0), float(1.0)]

cluster_dict["15.1260004044_arrows"] += cgo_arrow([38.5,-25.0,0.0], [37.904,-22.855,1.023], color="blue red", name="Arrows_15.1260004044_4")

cluster_dict["15.1260004044"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(31.4169215532), float(-43.3567809865), float(3.67457401371), float(1.0)]


cluster_dict["15.1260004044"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(38.5727701659), float(-26.82189472), float(-1.50137258236), float(1.0)]


cluster_dict["15.1260004044"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(36.0503854534), float(-37.7014227167), float(-0.633929414707), float(1.0)]


cluster_dict["15.1260004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(-28.5), float(-4.5), float(1.0)]

cluster_dict["15.1260004044_arrows"] += cgo_arrow([36.5,-28.5,-4.5], [37.132,-30.946,-6.482], color="red blue", name="Arrows_15.1260004044_5")

cluster_dict["15.1260004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(39.0), float(-23.5), float(-4.0), float(1.0)]

cluster_dict["15.1260004044_arrows"] += cgo_arrow([39.0,-23.5,-4.0], [36.528,-21.922,-3.46], color="red blue", name="Arrows_15.1260004044_6")

cluster_dict["15.1260004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(40.5), float(-27.5), float(2.0), float(1.0)]

cluster_dict["15.1260004044_arrows"] += cgo_arrow([40.5,-27.5,2.0], [43.271,-28.102,3.157], color="red blue", name="Arrows_15.1260004044_7")

cmd.load_cgo(cluster_dict["15.1260004044"], "Features_15.1260004044", 1)
cmd.load_cgo(cluster_dict["15.1260004044_arrows"], "Arrows_15.1260004044")
cmd.set("transparency", 0.2,"Features_15.1260004044")
cmd.group("Pharmacophore_15.1260004044", members="Features_15.1260004044")
cmd.group("Pharmacophore_15.1260004044", members="Arrows_15.1260004044")

if dirpath:
    f = join(dirpath, "13/label_threshold_15.1260004044.mol2")
else:
    f = "13/label_threshold_15.1260004044.mol2"

cmd.load(f, 'label_threshold_15.1260004044')
cmd.hide('everything', 'label_threshold_15.1260004044')
cmd.label("label_threshold_15.1260004044", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.1260004044', members= 'label_threshold_15.1260004044')


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
    f = join(dirpath, "14/label_threshold_13.3.mol2")
else:
    f = "14/label_threshold_13.3.mol2"

cmd.load(f, 'label_threshold_13.3')
cmd.hide('everything', 'label_threshold_13.3')
cmd.label("label_threshold_13.3", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.3]
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


cluster_dict = {"15.0850000381":[], "15.0850000381_arrows":[]}

cluster_dict["15.0850000381"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(-36.0), float(25.0), float(1.0)]

cluster_dict["15.0850000381_arrows"] += cgo_arrow([18.5,-36.0,25.0], [17.839,-38.095,23.656], color="blue red", name="Arrows_15.0850000381_1")

cluster_dict["15.0850000381"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(-44.0), float(26.5), float(1.0)]

cluster_dict["15.0850000381_arrows"] += cgo_arrow([18.5,-44.0,26.5], [15.475,-45.652,25.95], color="blue red", name="Arrows_15.0850000381_2")

cluster_dict["15.0850000381"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(-44.0), float(26.5), float(1.0)]

cluster_dict["15.0850000381_arrows"] += cgo_arrow([18.5,-44.0,26.5], [15.475,-45.652,25.95], color="blue red", name="Arrows_15.0850000381_3")

cluster_dict["15.0850000381"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(19.0332476795), float(-44.5082179855), float(25.4970313898), float(1.0)]


cluster_dict["15.0850000381"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(18.2500280802), float(-36.5), float(22.1), float(1.0)]


cluster_dict["15.0850000381"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(19.3151669054), float(-36.2636536851), float(22.9662785976), float(1.0)]


cluster_dict["15.0850000381"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(-40.5), float(28.5), float(1.0)]

cluster_dict["15.0850000381_arrows"] += cgo_arrow([16.0,-40.5,28.5], [13.403,-39.381,29.79], color="red blue", name="Arrows_15.0850000381_4")

cluster_dict["15.0850000381"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.0), float(-44.0), float(23.5), float(1.0)]

cluster_dict["15.0850000381_arrows"] += cgo_arrow([20.0,-44.0,23.5], [22.794,-42.19,22.011], color="red blue", name="Arrows_15.0850000381_5")

cluster_dict["15.0850000381"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.0), float(-44.0), float(23.5), float(1.0)]

cluster_dict["15.0850000381_arrows"] += cgo_arrow([20.0,-44.0,23.5], [22.794,-42.19,22.011], color="red blue", name="Arrows_15.0850000381_6")

cluster_dict["15.0850000381"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.5), float(-36.0), float(24.0), float(1.0)]

cluster_dict["15.0850000381_arrows"] += cgo_arrow([20.5,-36.0,24.0], [21.795,-37.163,26.941], color="red blue", name="Arrows_15.0850000381_7")

cluster_dict["15.0850000381"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.0), float(-36.5), float(22.0), float(1.0)]

cluster_dict["15.0850000381_arrows"] += cgo_arrow([20.0,-36.5,22.0], [17.086,-39.498,22.053], color="red blue", name="Arrows_15.0850000381_8")

cluster_dict["15.0850000381"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(-49.0), float(23.5), float(1.0)]

cluster_dict["15.0850000381_arrows"] += cgo_arrow([21.0,-49.0,23.5], [21.007,-48.652,20.493], color="red blue", name="Arrows_15.0850000381_9")

cmd.load_cgo(cluster_dict["15.0850000381"], "Features_15.0850000381", 1)
cmd.load_cgo(cluster_dict["15.0850000381_arrows"], "Arrows_15.0850000381")
cmd.set("transparency", 0.2,"Features_15.0850000381")
cmd.group("Pharmacophore_15.0850000381", members="Features_15.0850000381")
cmd.group("Pharmacophore_15.0850000381", members="Arrows_15.0850000381")

if dirpath:
    f = join(dirpath, "14/label_threshold_15.0850000381.mol2")
else:
    f = "14/label_threshold_15.0850000381.mol2"

cmd.load(f, 'label_threshold_15.0850000381')
cmd.hide('everything', 'label_threshold_15.0850000381')
cmd.label("label_threshold_15.0850000381", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.0850000381', members= 'label_threshold_15.0850000381')


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
    f = join(dirpath, "15/label_threshold_12.3.mol2")
else:
    f = "15/label_threshold_12.3.mol2"

cmd.load(f, 'label_threshold_12.3')
cmd.hide('everything', 'label_threshold_12.3')
cmd.label("label_threshold_12.3", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.3]
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


cluster_dict = {"15.0749998093":[], "15.0749998093_arrows":[]}

cluster_dict["15.0749998093"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(36.0), float(-35.5), float(-1.0), float(1.0)]

cluster_dict["15.0749998093_arrows"] += cgo_arrow([36.0,-35.5,-1.0], [35.162,-31.753,-0.868], color="blue red", name="Arrows_15.0749998093_1")

cluster_dict["15.0749998093"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(36.0), float(-35.5), float(-1.0), float(1.0)]

cluster_dict["15.0749998093_arrows"] += cgo_arrow([36.0,-35.5,-1.0], [35.162,-31.753,-0.868], color="blue red", name="Arrows_15.0749998093_2")

cluster_dict["15.0749998093"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(38.5), float(-25.0), float(0.0), float(1.0)]

cluster_dict["15.0749998093_arrows"] += cgo_arrow([38.5,-25.0,0.0], [37.904,-22.855,1.023], color="blue red", name="Arrows_15.0749998093_3")

cluster_dict["15.0749998093"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(38.7182882862), float(-26.7314564146), float(-1.51568532478), float(1.0)]


cluster_dict["15.0749998093"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(36.0529861577), float(-37.7332648403), float(-0.727823894393), float(1.0)]


cluster_dict["15.0749998093"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(-28.5), float(-4.5), float(1.0)]

cluster_dict["15.0749998093_arrows"] += cgo_arrow([36.5,-28.5,-4.5], [37.132,-30.946,-6.482], color="red blue", name="Arrows_15.0749998093_4")

cluster_dict["15.0749998093"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(39.0), float(-23.5), float(-4.0), float(1.0)]

cluster_dict["15.0749998093_arrows"] += cgo_arrow([39.0,-23.5,-4.0], [36.528,-21.922,-3.46], color="red blue", name="Arrows_15.0749998093_5")

cluster_dict["15.0749998093"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(40.5), float(-27.5), float(2.0), float(1.0)]

cluster_dict["15.0749998093_arrows"] += cgo_arrow([40.5,-27.5,2.0], [43.271,-28.102,3.157], color="red blue", name="Arrows_15.0749998093_6")

cmd.load_cgo(cluster_dict["15.0749998093"], "Features_15.0749998093", 1)
cmd.load_cgo(cluster_dict["15.0749998093_arrows"], "Arrows_15.0749998093")
cmd.set("transparency", 0.2,"Features_15.0749998093")
cmd.group("Pharmacophore_15.0749998093", members="Features_15.0749998093")
cmd.group("Pharmacophore_15.0749998093", members="Arrows_15.0749998093")

if dirpath:
    f = join(dirpath, "15/label_threshold_15.0749998093.mol2")
else:
    f = "15/label_threshold_15.0749998093.mol2"

cmd.load(f, 'label_threshold_15.0749998093')
cmd.hide('everything', 'label_threshold_15.0749998093')
cmd.label("label_threshold_15.0749998093", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.0749998093', members= 'label_threshold_15.0749998093')


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
    f = join(dirpath, "16/label_threshold_12.8.mol2")
else:
    f = "16/label_threshold_12.8.mol2"

cmd.load(f, 'label_threshold_12.8')
cmd.hide('everything', 'label_threshold_12.8')
cmd.label("label_threshold_12.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.8]
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


cluster_dict = {"14.7399997711":[], "14.7399997711_arrows":[]}

cluster_dict["14.7399997711"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.0), float(-31.0), float(21.5), float(1.0)]

cluster_dict["14.7399997711_arrows"] += cgo_arrow([19.0,-31.0,21.5], [15.264,-32.638,22.22], color="blue red", name="Arrows_14.7399997711_1")

cluster_dict["14.7399997711"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.0), float(-31.0), float(21.5), float(1.0)]

cluster_dict["14.7399997711_arrows"] += cgo_arrow([19.0,-31.0,21.5], [15.264,-32.638,22.22], color="blue red", name="Arrows_14.7399997711_2")

cluster_dict["14.7399997711"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(20.5), float(-33.0), float(26.0), float(1.0)]

cluster_dict["14.7399997711_arrows"] += cgo_arrow([20.5,-33.0,26.0], [20.162,-31.121,28.672], color="blue red", name="Arrows_14.7399997711_3")

cluster_dict["14.7399997711"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.0), float(-26.5), float(16.0), float(1.0)]

cluster_dict["14.7399997711_arrows"] += cgo_arrow([29.0,-26.5,16.0], [27.469,-24.336,14.801], color="blue red", name="Arrows_14.7399997711_4")

cluster_dict["14.7399997711"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(20.8727115296), float(-30.9851875476), float(22.7942605172), float(1.0)]


cluster_dict["14.7399997711"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(28.479841117), float(-26.0752281469), float(19.5143209604), float(1.0)]


cluster_dict["14.7399997711"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(-32.4383212106), float(28.9691407786), float(1.0)]


cluster_dict["14.7399997711"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(36.3634887665), float(-25.1801162355), float(20.0853889465), float(1.0)]


cluster_dict["14.7399997711"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(37.895945019), float(-31.7493058494), float(15.2686160489), float(1.0)]


cluster_dict["14.7399997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.0), float(-29.5), float(22.0), float(1.0)]

cluster_dict["14.7399997711_arrows"] += cgo_arrow([17.0,-29.5,22.0], [13.97,-30.815,22.039], color="red blue", name="Arrows_14.7399997711_5")

cluster_dict["14.7399997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(-32.0), float(17.0), float(1.0)]

cluster_dict["14.7399997711_arrows"] += cgo_arrow([18.5,-32.0,17.0], [20.191,-30.573,14.342], color="red blue", name="Arrows_14.7399997711_6")

cluster_dict["14.7399997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(-34.0), float(21.0), float(1.0)]

cluster_dict["14.7399997711_arrows"] += cgo_arrow([18.5,-34.0,21.0], [17.274,-35.929,18.15], color="red blue", name="Arrows_14.7399997711_7")

cluster_dict["14.7399997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.5), float(-35.5), float(23.5), float(1.0)]

cluster_dict["14.7399997711_arrows"] += cgo_arrow([20.5,-35.5,23.5], [21.795,-37.163,26.941], color="red blue", name="Arrows_14.7399997711_8")

cluster_dict["14.7399997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(-31.0), float(24.0), float(1.0)]

cluster_dict["14.7399997711_arrows"] += cgo_arrow([23.0,-31.0,24.0], [24.932,-30.214,21.553], color="red blue", name="Arrows_14.7399997711_9")

cmd.load_cgo(cluster_dict["14.7399997711"], "Features_14.7399997711", 1)
cmd.load_cgo(cluster_dict["14.7399997711_arrows"], "Arrows_14.7399997711")
cmd.set("transparency", 0.2,"Features_14.7399997711")
cmd.group("Pharmacophore_14.7399997711", members="Features_14.7399997711")
cmd.group("Pharmacophore_14.7399997711", members="Arrows_14.7399997711")

if dirpath:
    f = join(dirpath, "16/label_threshold_14.7399997711.mol2")
else:
    f = "16/label_threshold_14.7399997711.mol2"

cmd.load(f, 'label_threshold_14.7399997711')
cmd.hide('everything', 'label_threshold_14.7399997711')
cmd.label("label_threshold_14.7399997711", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.7399997711', members= 'label_threshold_14.7399997711')


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
    f = join(dirpath, "17/label_threshold_11.5.mol2")
else:
    f = "17/label_threshold_11.5.mol2"

cmd.load(f, 'label_threshold_11.5')
cmd.hide('everything', 'label_threshold_11.5')
cmd.label("label_threshold_11.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [11.5]
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


cluster_dict = {"14.6309995651":[], "14.6309995651_arrows":[]}

cluster_dict["14.6309995651"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(-43.5), float(4.0), float(1.0)]

cluster_dict["14.6309995651_arrows"] += cgo_arrow([29.5,-43.5,4.0], [29.393,-40.692,4.576], color="blue red", name="Arrows_14.6309995651_1")

cluster_dict["14.6309995651"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(35.5), float(-35.5), float(-1.0), float(1.0)]

cluster_dict["14.6309995651_arrows"] += cgo_arrow([35.5,-35.5,-1.0], [33.351,-35.839,-4.065], color="blue red", name="Arrows_14.6309995651_2")

cluster_dict["14.6309995651"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(35.5), float(-35.5), float(-1.0), float(1.0)]

cluster_dict["14.6309995651_arrows"] += cgo_arrow([35.5,-35.5,-1.0], [33.351,-35.839,-4.065], color="blue red", name="Arrows_14.6309995651_3")

cluster_dict["14.6309995651"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(31.4075174442), float(-43.3475890629), float(3.68558489096), float(1.0)]


cluster_dict["14.6309995651"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(37.7627149749), float(-27.5349111774), float(-1.1222617153), float(1.0)]


cluster_dict["14.6309995651"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(36.0171850212), float(-37.7495738776), float(-0.959856382421), float(1.0)]


cluster_dict["14.6309995651"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(36.0), float(-37.5), float(-5.5), float(1.0)]


cluster_dict["14.6309995651"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(-28.5), float(-4.5), float(1.0)]

cluster_dict["14.6309995651_arrows"] += cgo_arrow([36.5,-28.5,-4.5], [37.132,-30.946,-6.482], color="red blue", name="Arrows_14.6309995651_4")

cluster_dict["14.6309995651"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(40.5), float(-27.5), float(2.0), float(1.0)]

cluster_dict["14.6309995651_arrows"] += cgo_arrow([40.5,-27.5,2.0], [43.271,-28.102,3.157], color="red blue", name="Arrows_14.6309995651_5")

cmd.load_cgo(cluster_dict["14.6309995651"], "Features_14.6309995651", 1)
cmd.load_cgo(cluster_dict["14.6309995651_arrows"], "Arrows_14.6309995651")
cmd.set("transparency", 0.2,"Features_14.6309995651")
cmd.group("Pharmacophore_14.6309995651", members="Features_14.6309995651")
cmd.group("Pharmacophore_14.6309995651", members="Arrows_14.6309995651")

if dirpath:
    f = join(dirpath, "17/label_threshold_14.6309995651.mol2")
else:
    f = "17/label_threshold_14.6309995651.mol2"

cmd.load(f, 'label_threshold_14.6309995651')
cmd.hide('everything', 'label_threshold_14.6309995651')
cmd.label("label_threshold_14.6309995651", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.6309995651', members= 'label_threshold_14.6309995651')


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
    f = join(dirpath, "18/label_threshold_10.5.mol2")
else:
    f = "18/label_threshold_10.5.mol2"

cmd.load(f, 'label_threshold_10.5')
cmd.hide('everything', 'label_threshold_10.5')
cmd.label("label_threshold_10.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.5]
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


cluster_dict = {"14.5200004578":[], "14.5200004578_arrows":[]}

cluster_dict["14.5200004578"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(36.0), float(-34.5), float(-1.0), float(1.0)]

cluster_dict["14.5200004578_arrows"] += cgo_arrow([36.0,-34.5,-1.0], [35.162,-31.753,-0.868], color="blue red", name="Arrows_14.5200004578_1")

cluster_dict["14.5200004578"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(-36.0), float(-4.0), float(1.0)]

cluster_dict["14.5200004578_arrows"] += cgo_arrow([36.5,-36.0,-4.0], [33.351,-35.839,-4.065], color="blue red", name="Arrows_14.5200004578_2")

cluster_dict["14.5200004578"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(38.5), float(-25.0), float(0.0), float(1.0)]

cluster_dict["14.5200004578_arrows"] += cgo_arrow([38.5,-25.0,0.0], [37.904,-22.855,1.023], color="blue red", name="Arrows_14.5200004578_3")

cluster_dict["14.5200004578"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(38.8468691937), float(-26.7496661103), float(-1.47969838646), float(1.0)]


cluster_dict["14.5200004578"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(36.3383055946), float(-34.4811835163), float(-1.56143115676), float(1.0)]


cluster_dict["14.5200004578"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(-28.5), float(-4.5), float(1.0)]

cluster_dict["14.5200004578_arrows"] += cgo_arrow([36.5,-28.5,-4.5], [37.132,-30.946,-6.482], color="red blue", name="Arrows_14.5200004578_4")

cluster_dict["14.5200004578"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(39.0), float(-23.5), float(-4.0), float(1.0)]

cluster_dict["14.5200004578_arrows"] += cgo_arrow([39.0,-23.5,-4.0], [36.528,-21.922,-3.46], color="red blue", name="Arrows_14.5200004578_5")

cluster_dict["14.5200004578"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(40.5), float(-27.5), float(2.0), float(1.0)]

cluster_dict["14.5200004578_arrows"] += cgo_arrow([40.5,-27.5,2.0], [43.271,-28.102,3.157], color="red blue", name="Arrows_14.5200004578_6")

cmd.load_cgo(cluster_dict["14.5200004578"], "Features_14.5200004578", 1)
cmd.load_cgo(cluster_dict["14.5200004578_arrows"], "Arrows_14.5200004578")
cmd.set("transparency", 0.2,"Features_14.5200004578")
cmd.group("Pharmacophore_14.5200004578", members="Features_14.5200004578")
cmd.group("Pharmacophore_14.5200004578", members="Arrows_14.5200004578")

if dirpath:
    f = join(dirpath, "18/label_threshold_14.5200004578.mol2")
else:
    f = "18/label_threshold_14.5200004578.mol2"

cmd.load(f, 'label_threshold_14.5200004578')
cmd.hide('everything', 'label_threshold_14.5200004578')
cmd.label("label_threshold_14.5200004578", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.5200004578', members= 'label_threshold_14.5200004578')


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
    f = join(dirpath, "19/label_threshold_9.5.mol2")
else:
    f = "19/label_threshold_9.5.mol2"

cmd.load(f, 'label_threshold_9.5')
cmd.hide('everything', 'label_threshold_9.5')
cmd.label("label_threshold_9.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [9.5]
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


cluster_dict = {"14.0740003586":[], "14.0740003586_arrows":[]}

cluster_dict["14.0740003586"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(38.5), float(-25.0), float(0.0), float(1.0)]

cluster_dict["14.0740003586_arrows"] += cgo_arrow([38.5,-25.0,0.0], [37.904,-22.855,1.023], color="blue red", name="Arrows_14.0740003586_1")

cluster_dict["14.0740003586"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(39.038131242), float(-26.7692356885), float(-1.56818975172), float(1.0)]


cluster_dict["14.0740003586"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(37.218029479), float(-32.9359460235), float(-2.03699557676), float(1.0)]


cluster_dict["14.0740003586"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(-28.5), float(-4.5), float(1.0)]

cluster_dict["14.0740003586_arrows"] += cgo_arrow([36.5,-28.5,-4.5], [37.132,-30.946,-6.482], color="red blue", name="Arrows_14.0740003586_2")

cluster_dict["14.0740003586"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(39.0), float(-23.5), float(-4.0), float(1.0)]

cluster_dict["14.0740003586_arrows"] += cgo_arrow([39.0,-23.5,-4.0], [36.528,-21.922,-3.46], color="red blue", name="Arrows_14.0740003586_3")

cluster_dict["14.0740003586"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(40.5), float(-27.5), float(2.0), float(1.0)]

cluster_dict["14.0740003586_arrows"] += cgo_arrow([40.5,-27.5,2.0], [43.271,-28.102,3.157], color="red blue", name="Arrows_14.0740003586_4")

cmd.load_cgo(cluster_dict["14.0740003586"], "Features_14.0740003586", 1)
cmd.load_cgo(cluster_dict["14.0740003586_arrows"], "Arrows_14.0740003586")
cmd.set("transparency", 0.2,"Features_14.0740003586")
cmd.group("Pharmacophore_14.0740003586", members="Features_14.0740003586")
cmd.group("Pharmacophore_14.0740003586", members="Arrows_14.0740003586")

if dirpath:
    f = join(dirpath, "19/label_threshold_14.0740003586.mol2")
else:
    f = "19/label_threshold_14.0740003586.mol2"

cmd.load(f, 'label_threshold_14.0740003586')
cmd.hide('everything', 'label_threshold_14.0740003586')
cmd.label("label_threshold_14.0740003586", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.0740003586', members= 'label_threshold_14.0740003586')


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
    f = join(dirpath, "20/label_threshold_11.6.mol2")
else:
    f = "20/label_threshold_11.6.mol2"

cmd.load(f, 'label_threshold_11.6')
cmd.hide('everything', 'label_threshold_11.6')
cmd.label("label_threshold_11.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [11.6]
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


cluster_dict = {"13.9555001259":[], "13.9555001259_arrows":[]}

cluster_dict["13.9555001259"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(6.5), float(-0.5), float(1.0)]

cluster_dict["13.9555001259_arrows"] += cgo_arrow([30.5,6.5,-0.5], [28.734,6.913,-2.582], color="blue red", name="Arrows_13.9555001259_1")

cluster_dict["13.9555001259"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(31.8827578632), float(1.85750875869), float(-0.727195146032), float(1.0)]


cluster_dict["13.9555001259"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(41.4698348192), float(2.84517415584), float(-2.1133533582), float(1.0)]


cluster_dict["13.9555001259"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(31.5), float(-0.5), float(2.0), float(1.0)]

cluster_dict["13.9555001259_arrows"] += cgo_arrow([31.5,-0.5,2.0], [33.219,-2.042,2.779], color="red blue", name="Arrows_13.9555001259_2")

cmd.load_cgo(cluster_dict["13.9555001259"], "Features_13.9555001259", 1)
cmd.load_cgo(cluster_dict["13.9555001259_arrows"], "Arrows_13.9555001259")
cmd.set("transparency", 0.2,"Features_13.9555001259")
cmd.group("Pharmacophore_13.9555001259", members="Features_13.9555001259")
cmd.group("Pharmacophore_13.9555001259", members="Arrows_13.9555001259")

if dirpath:
    f = join(dirpath, "20/label_threshold_13.9555001259.mol2")
else:
    f = "20/label_threshold_13.9555001259.mol2"

cmd.load(f, 'label_threshold_13.9555001259')
cmd.hide('everything', 'label_threshold_13.9555001259')
cmd.label("label_threshold_13.9555001259", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.9555001259', members= 'label_threshold_13.9555001259')


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
    f = join(dirpath, "21/label_threshold_11.6.mol2")
else:
    f = "21/label_threshold_11.6.mol2"

cmd.load(f, 'label_threshold_11.6')
cmd.hide('everything', 'label_threshold_11.6')
cmd.label("label_threshold_11.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [11.6]
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


cluster_dict = {"13.9350004196":[], "13.9350004196_arrows":[]}

cluster_dict["13.9350004196"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(7.0), float(0.0), float(1.0)]

cluster_dict["13.9350004196_arrows"] += cgo_arrow([30.5,7.0,0.0], [31.908,8.68,1.718], color="blue red", name="Arrows_13.9350004196_1")

cluster_dict["13.9350004196"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(31.80193098), float(1.92879676784), float(-0.931573361237), float(1.0)]


cluster_dict["13.9350004196"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(41.4698348192), float(2.84517415584), float(-2.1133533582), float(1.0)]


cluster_dict["13.9350004196"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(31.5), float(-0.5), float(2.0), float(1.0)]

cluster_dict["13.9350004196_arrows"] += cgo_arrow([31.5,-0.5,2.0], [33.219,-2.042,2.779], color="red blue", name="Arrows_13.9350004196_2")

cmd.load_cgo(cluster_dict["13.9350004196"], "Features_13.9350004196", 1)
cmd.load_cgo(cluster_dict["13.9350004196_arrows"], "Arrows_13.9350004196")
cmd.set("transparency", 0.2,"Features_13.9350004196")
cmd.group("Pharmacophore_13.9350004196", members="Features_13.9350004196")
cmd.group("Pharmacophore_13.9350004196", members="Arrows_13.9350004196")

if dirpath:
    f = join(dirpath, "21/label_threshold_13.9350004196.mol2")
else:
    f = "21/label_threshold_13.9350004196.mol2"

cmd.load(f, 'label_threshold_13.9350004196')
cmd.hide('everything', 'label_threshold_13.9350004196')
cmd.label("label_threshold_13.9350004196", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.9350004196', members= 'label_threshold_13.9350004196')


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
    f = join(dirpath, "22/label_threshold_10.5.mol2")
else:
    f = "22/label_threshold_10.5.mol2"

cmd.load(f, 'label_threshold_10.5')
cmd.hide('everything', 'label_threshold_10.5')
cmd.label("label_threshold_10.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.5]
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


cluster_dict = {"13.7419996262":[], "13.7419996262_arrows":[]}

cluster_dict["13.7419996262"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(-9.5), float(14.0), float(1.0)]

cluster_dict["13.7419996262_arrows"] += cgo_arrow([26.5,-9.5,14.0], [28.783,-7.532,13.582], color="blue red", name="Arrows_13.7419996262_1")

cluster_dict["13.7419996262"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(-3.0), float(12.5), float(1.0)]

cluster_dict["13.7419996262_arrows"] += cgo_arrow([24.0,-3.0,12.5], [21.766,-1.319,13.37], color="blue red", name="Arrows_13.7419996262_2")

cluster_dict["13.7419996262"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(2.0), float(2.0), float(1.0)]

cluster_dict["13.7419996262_arrows"] += cgo_arrow([26.0,2.0,2.0], [23.655,2.241,1.129], color="blue red", name="Arrows_13.7419996262_3")

cluster_dict["13.7419996262"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.5867799887), float(-11.5121570333), float(9.15035770299), float(1.0)]


cluster_dict["13.7419996262"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(23.9427420159), float(5.53011168496), float(18.0928887152), float(1.0)]


cluster_dict["13.7419996262"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(28.3490227797), float(-0.0706966591752), float(3.14272382178), float(1.0)]


cluster_dict["13.7419996262"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(33.3326539876), float(1.59999519252), float(12.3137894672), float(1.0)]


cluster_dict["13.7419996262"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(-13.0), float(14.5), float(1.0)]

cluster_dict["13.7419996262_arrows"] += cgo_arrow([18.0,-13.0,14.5], [16.48,-10.765,14.656], color="red blue", name="Arrows_13.7419996262_4")

cluster_dict["13.7419996262"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(-8.5), float(7.0), float(1.0)]

cluster_dict["13.7419996262_arrows"] += cgo_arrow([22.5,-8.5,7.0], [22.101,-6.363,4.763], color="red blue", name="Arrows_13.7419996262_5")

cluster_dict["13.7419996262"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(25.5), float(0.0), float(12.0), float(1.0)]

cluster_dict["13.7419996262_arrows"] += cgo_arrow([25.5,0.0,12.0], [28.5,3.731,11.033], color="red blue", name="Arrows_13.7419996262_6")

cluster_dict["13.7419996262"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(-2.5), float(11.5), float(1.0)]

cluster_dict["13.7419996262_arrows"] += cgo_arrow([27.0,-2.5,11.5], [25.888,-5.209,9.814], color="red blue", name="Arrows_13.7419996262_7")

cluster_dict["13.7419996262"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(-2.5), float(4.0), float(1.0)]

cluster_dict["13.7419996262_arrows"] += cgo_arrow([26.0,-2.5,4.0], [26.993,-4.956,2.181], color="red blue", name="Arrows_13.7419996262_8")

cluster_dict["13.7419996262"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(31.5), float(-0.5), float(2.0), float(1.0)]

cluster_dict["13.7419996262_arrows"] += cgo_arrow([31.5,-0.5,2.0], [33.219,-2.042,2.779], color="red blue", name="Arrows_13.7419996262_9")

cmd.load_cgo(cluster_dict["13.7419996262"], "Features_13.7419996262", 1)
cmd.load_cgo(cluster_dict["13.7419996262_arrows"], "Arrows_13.7419996262")
cmd.set("transparency", 0.2,"Features_13.7419996262")
cmd.group("Pharmacophore_13.7419996262", members="Features_13.7419996262")
cmd.group("Pharmacophore_13.7419996262", members="Arrows_13.7419996262")

if dirpath:
    f = join(dirpath, "22/label_threshold_13.7419996262.mol2")
else:
    f = "22/label_threshold_13.7419996262.mol2"

cmd.load(f, 'label_threshold_13.7419996262')
cmd.hide('everything', 'label_threshold_13.7419996262')
cmd.label("label_threshold_13.7419996262", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.7419996262', members= 'label_threshold_13.7419996262')


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
    f = join(dirpath, "23/label_threshold_11.5.mol2")
else:
    f = "23/label_threshold_11.5.mol2"

cmd.load(f, 'label_threshold_11.5')
cmd.hide('everything', 'label_threshold_11.5')
cmd.label("label_threshold_11.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [11.5]
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


cluster_dict = {"13.6949996948":[], "13.6949996948_arrows":[]}

cluster_dict["13.6949996948"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(-43.0), float(25.0), float(1.0)]

cluster_dict["13.6949996948_arrows"] += cgo_arrow([21.0,-43.0,25.0], [19.688,-40.611,22.652], color="blue red", name="Arrows_13.6949996948_1")

cluster_dict["13.6949996948"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(18.1455637264), float(-46.4642957146), float(21.7121446876), float(1.0)]


cluster_dict["13.6949996948"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(20.3853712315), float(-50.7936682167), float(15.7434501278), float(1.0)]


cluster_dict["13.6949996948"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.0), float(-44.0), float(24.5), float(1.0)]

cluster_dict["13.6949996948_arrows"] += cgo_arrow([17.0,-44.0,24.5], [15.475,-45.652,25.95], color="red blue", name="Arrows_13.6949996948_2")

cluster_dict["13.6949996948"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(-48.5), float(23.5), float(1.0)]

cluster_dict["13.6949996948_arrows"] += cgo_arrow([16.0,-48.5,23.5], [16.145,-50.262,20.162], color="red blue", name="Arrows_13.6949996948_3")

cluster_dict["13.6949996948"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(-47.0), float(20.0), float(1.0)]

cluster_dict["13.6949996948_arrows"] += cgo_arrow([18.5,-47.0,20.0], [21.007,-48.652,20.493], color="red blue", name="Arrows_13.6949996948_4")

cluster_dict["13.6949996948"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.5), float(-44.0), float(23.5), float(1.0)]

cluster_dict["13.6949996948_arrows"] += cgo_arrow([20.5,-44.0,23.5], [22.794,-42.19,22.011], color="red blue", name="Arrows_13.6949996948_5")

cluster_dict["13.6949996948"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.5), float(-48.5), float(16.0), float(1.0)]

cluster_dict["13.6949996948_arrows"] += cgo_arrow([21.5,-48.5,16.0], [21.784,-47.8,13.196], color="red blue", name="Arrows_13.6949996948_6")

cluster_dict["13.6949996948"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(-49.0), float(23.5), float(1.0)]

cluster_dict["13.6949996948_arrows"] += cgo_arrow([21.0,-49.0,23.5], [21.007,-48.652,20.493], color="red blue", name="Arrows_13.6949996948_7")

cmd.load_cgo(cluster_dict["13.6949996948"], "Features_13.6949996948", 1)
cmd.load_cgo(cluster_dict["13.6949996948_arrows"], "Arrows_13.6949996948")
cmd.set("transparency", 0.2,"Features_13.6949996948")
cmd.group("Pharmacophore_13.6949996948", members="Features_13.6949996948")
cmd.group("Pharmacophore_13.6949996948", members="Arrows_13.6949996948")

if dirpath:
    f = join(dirpath, "23/label_threshold_13.6949996948.mol2")
else:
    f = "23/label_threshold_13.6949996948.mol2"

cmd.load(f, 'label_threshold_13.6949996948')
cmd.hide('everything', 'label_threshold_13.6949996948')
cmd.label("label_threshold_13.6949996948", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.6949996948', members= 'label_threshold_13.6949996948')


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
    f = join(dirpath, "24/label_threshold_8.9.mol2")
else:
    f = "24/label_threshold_8.9.mol2"

cmd.load(f, 'label_threshold_8.9')
cmd.hide('everything', 'label_threshold_8.9')
cmd.label("label_threshold_8.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [8.9]
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


cluster_dict = {"13.2060003281":[], "13.2060003281_arrows":[]}

cluster_dict["13.2060003281"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(43.5), float(-8.5), float(9.5), float(1.0)]

cluster_dict["13.2060003281_arrows"] += cgo_arrow([43.5,-8.5,9.5], [43.181,-9.676,7.212], color="blue red", name="Arrows_13.2060003281_1")

cluster_dict["13.2060003281"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(44.5), float(-0.5), float(11.5), float(1.0)]

cluster_dict["13.2060003281_arrows"] += cgo_arrow([44.5,-0.5,11.5], [47.549,0.315,10.678], color="blue red", name="Arrows_13.2060003281_2")

cluster_dict["13.2060003281"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(48.0), float(-6.0), float(8.0), float(1.0)]

cluster_dict["13.2060003281_arrows"] += cgo_arrow([48.0,-6.0,8.0], [45.592,-6.463,6.489], color="blue red", name="Arrows_13.2060003281_3")

cluster_dict["13.2060003281"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(41.4414422947), float(3.68747959748), float(1.70346262478), float(1.0)]


cluster_dict["13.2060003281"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(45.1875699499), float(-5.45839100727), float(9.9474714477), float(1.0)]


cluster_dict["13.2060003281"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(43.5600614954), float(3.07829969033), float(10.5201798413), float(1.0)]


cluster_dict["13.2060003281"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(48.1811893186), float(-7.81881068136), float(8.90940534068), float(1.0)]


cluster_dict["13.2060003281"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(47.0), float(-7.0), float(9.0), float(1.0)]

cluster_dict["13.2060003281_arrows"] += cgo_arrow([47.0,-7.0,9.0], [45.592,-6.463,6.489], color="red blue", name="Arrows_13.2060003281_4")

cluster_dict["13.2060003281"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(45.0), float(-3.5), float(10.0), float(1.0)]

cluster_dict["13.2060003281_arrows"] += cgo_arrow([45.0,-3.5,10.0], [45.601,-1.014,8.286], color="red blue", name="Arrows_13.2060003281_5")

cluster_dict["13.2060003281"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(48.0), float(-7.0), float(5.5), float(1.0)]

cluster_dict["13.2060003281_arrows"] += cgo_arrow([48.0,-7.0,5.5], [45.592,-6.463,6.489], color="red blue", name="Arrows_13.2060003281_6")

cluster_dict["13.2060003281"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(51.0), float(-5.5), float(3.5), float(1.0)]

cluster_dict["13.2060003281_arrows"] += cgo_arrow([51.0,-5.5,3.5], [54.537,-4.976,2.652], color="red blue", name="Arrows_13.2060003281_7")

cmd.load_cgo(cluster_dict["13.2060003281"], "Features_13.2060003281", 1)
cmd.load_cgo(cluster_dict["13.2060003281_arrows"], "Arrows_13.2060003281")
cmd.set("transparency", 0.2,"Features_13.2060003281")
cmd.group("Pharmacophore_13.2060003281", members="Features_13.2060003281")
cmd.group("Pharmacophore_13.2060003281", members="Arrows_13.2060003281")

if dirpath:
    f = join(dirpath, "24/label_threshold_13.2060003281.mol2")
else:
    f = "24/label_threshold_13.2060003281.mol2"

cmd.load(f, 'label_threshold_13.2060003281')
cmd.hide('everything', 'label_threshold_13.2060003281')
cmd.label("label_threshold_13.2060003281", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.2060003281', members= 'label_threshold_13.2060003281')


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
    f = join(dirpath, "25/label_threshold_11.5.mol2")
else:
    f = "25/label_threshold_11.5.mol2"

cmd.load(f, 'label_threshold_11.5')
cmd.hide('everything', 'label_threshold_11.5')
cmd.label("label_threshold_11.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [11.5]
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


cluster_dict = {"12.8179998398":[], "12.8179998398_arrows":[]}

cluster_dict["12.8179998398"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(36.9505341131), float(-26.1058308565), float(20.0452048799), float(1.0)]


cluster_dict["12.8179998398"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(37.6575973535), float(-32.3172103056), float(14.6756673137), float(1.0)]


cluster_dict["12.8179998398"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(45.357990713), float(-30.5366188878), float(17.5257365499), float(1.0)]


cluster_dict["12.8179998398"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(44.4616010855), float(-22.6639540863), float(16.0355033586), float(1.0)]


cluster_dict["12.8179998398"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(52.8353097109), float(-40.7058140956), float(20.7556214403), float(1.0)]


cluster_dict["12.8179998398"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(55.8007709584), float(-39.8108216251), float(12.6786538277), float(1.0)]


cluster_dict["12.8179998398"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(42.5), float(-29.5), float(14.0), float(1.0)]

cluster_dict["12.8179998398_arrows"] += cgo_arrow([42.5,-29.5,14.0], [45.812,-29.166,13.284], color="red blue", name="Arrows_12.8179998398_1")

cluster_dict["12.8179998398"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(42.5), float(-29.5), float(14.0), float(1.0)]

cluster_dict["12.8179998398_arrows"] += cgo_arrow([42.5,-29.5,14.0], [45.812,-29.166,13.284], color="red blue", name="Arrows_12.8179998398_2")

cluster_dict["12.8179998398"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(50.5), float(-44.0), float(21.5), float(1.0)]

cluster_dict["12.8179998398_arrows"] += cgo_arrow([50.5,-44.0,21.5], [48.706,-42.054,18.711], color="red blue", name="Arrows_12.8179998398_3")

cluster_dict["12.8179998398"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(51.5), float(-41.5), float(22.5), float(1.0)]

cluster_dict["12.8179998398_arrows"] += cgo_arrow([51.5,-41.5,22.5], [48.706,-42.054,18.711], color="red blue", name="Arrows_12.8179998398_4")

cluster_dict["12.8179998398"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(53.5), float(-38.5), float(10.5), float(1.0)]

cluster_dict["12.8179998398_arrows"] += cgo_arrow([53.5,-38.5,10.5], [53.831,-40.753,12.95], color="red blue", name="Arrows_12.8179998398_5")

cmd.load_cgo(cluster_dict["12.8179998398"], "Features_12.8179998398", 1)
cmd.load_cgo(cluster_dict["12.8179998398_arrows"], "Arrows_12.8179998398")
cmd.set("transparency", 0.2,"Features_12.8179998398")
cmd.group("Pharmacophore_12.8179998398", members="Features_12.8179998398")
cmd.group("Pharmacophore_12.8179998398", members="Arrows_12.8179998398")

if dirpath:
    f = join(dirpath, "25/label_threshold_12.8179998398.mol2")
else:
    f = "25/label_threshold_12.8179998398.mol2"

cmd.load(f, 'label_threshold_12.8179998398')
cmd.hide('everything', 'label_threshold_12.8179998398')
cmd.label("label_threshold_12.8179998398", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.8179998398', members= 'label_threshold_12.8179998398')


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
    f = join(dirpath, "26/label_threshold_10.2.mol2")
else:
    f = "26/label_threshold_10.2.mol2"

cmd.load(f, 'label_threshold_10.2')
cmd.hide('everything', 'label_threshold_10.2')
cmd.label("label_threshold_10.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.2]
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


cluster_dict = {"12.7989997864":[], "12.7989997864_arrows":[]}

cluster_dict["12.7989997864"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(48.3058393721), float(-33.0094996568), float(17.9704053446), float(1.0)]


cluster_dict["12.7989997864"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(53.1346975785), float(-41.2465312706), float(20.7765980451), float(1.0)]


cluster_dict["12.7989997864"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(56.9419305973), float(-40.6077772871), float(11.9667772685), float(1.0)]


cluster_dict["12.7989997864"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(51.0), float(-44.5), float(21.5), float(1.0)]

cluster_dict["12.7989997864_arrows"] += cgo_arrow([51.0,-44.5,21.5], [48.706,-42.054,18.711], color="red blue", name="Arrows_12.7989997864_1")

cluster_dict["12.7989997864"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(53.5), float(-38.5), float(10.5), float(1.0)]

cluster_dict["12.7989997864_arrows"] += cgo_arrow([53.5,-38.5,10.5], [53.831,-40.753,12.95], color="red blue", name="Arrows_12.7989997864_2")

cluster_dict["12.7989997864"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(55.0), float(-42.0), float(9.5), float(1.0)]

cluster_dict["12.7989997864_arrows"] += cgo_arrow([55.0,-42.0,9.5], [54.318,-45.115,10.737], color="red blue", name="Arrows_12.7989997864_3")

cmd.load_cgo(cluster_dict["12.7989997864"], "Features_12.7989997864", 1)
cmd.load_cgo(cluster_dict["12.7989997864_arrows"], "Arrows_12.7989997864")
cmd.set("transparency", 0.2,"Features_12.7989997864")
cmd.group("Pharmacophore_12.7989997864", members="Features_12.7989997864")
cmd.group("Pharmacophore_12.7989997864", members="Arrows_12.7989997864")

if dirpath:
    f = join(dirpath, "26/label_threshold_12.7989997864.mol2")
else:
    f = "26/label_threshold_12.7989997864.mol2"

cmd.load(f, 'label_threshold_12.7989997864')
cmd.hide('everything', 'label_threshold_12.7989997864')
cmd.label("label_threshold_12.7989997864", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.7989997864', members= 'label_threshold_12.7989997864')


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
    f = join(dirpath, "27/label_threshold_9.0.mol2")
else:
    f = "27/label_threshold_9.0.mol2"

cmd.load(f, 'label_threshold_9.0')
cmd.hide('everything', 'label_threshold_9.0')
cmd.label("label_threshold_9.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [9.0]
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


cluster_dict = {"12.6540002823":[], "12.6540002823_arrows":[]}

cluster_dict["12.6540002823"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(43.5), float(-8.5), float(9.5), float(1.0)]

cluster_dict["12.6540002823_arrows"] += cgo_arrow([43.5,-8.5,9.5], [43.181,-9.676,7.212], color="blue red", name="Arrows_12.6540002823_1")

cluster_dict["12.6540002823"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(44.5), float(-0.5), float(11.5), float(1.0)]

cluster_dict["12.6540002823_arrows"] += cgo_arrow([44.5,-0.5,11.5], [47.549,0.315,10.678], color="blue red", name="Arrows_12.6540002823_2")

cluster_dict["12.6540002823"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(48.0), float(-6.0), float(8.0), float(1.0)]

cluster_dict["12.6540002823_arrows"] += cgo_arrow([48.0,-6.0,8.0], [45.592,-6.463,6.489], color="blue red", name="Arrows_12.6540002823_3")

cluster_dict["12.6540002823"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(44.4534251902), float(-5.75602283645), float(10.8310556277), float(1.0)]


cluster_dict["12.6540002823"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(43.5600614954), float(3.07829969033), float(10.5201798413), float(1.0)]


cluster_dict["12.6540002823"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(47.0), float(-7.0), float(9.0), float(1.0)]

cluster_dict["12.6540002823_arrows"] += cgo_arrow([47.0,-7.0,9.0], [45.592,-6.463,6.489], color="red blue", name="Arrows_12.6540002823_4")

cluster_dict["12.6540002823"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(45.0), float(-3.5), float(10.0), float(1.0)]

cluster_dict["12.6540002823_arrows"] += cgo_arrow([45.0,-3.5,10.0], [45.601,-1.014,8.286], color="red blue", name="Arrows_12.6540002823_5")

cluster_dict["12.6540002823"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(47.5), float(-7.0), float(5.5), float(1.0)]

cluster_dict["12.6540002823_arrows"] += cgo_arrow([47.5,-7.0,5.5], [45.592,-6.463,6.489], color="red blue", name="Arrows_12.6540002823_6")

cmd.load_cgo(cluster_dict["12.6540002823"], "Features_12.6540002823", 1)
cmd.load_cgo(cluster_dict["12.6540002823_arrows"], "Arrows_12.6540002823")
cmd.set("transparency", 0.2,"Features_12.6540002823")
cmd.group("Pharmacophore_12.6540002823", members="Features_12.6540002823")
cmd.group("Pharmacophore_12.6540002823", members="Arrows_12.6540002823")

if dirpath:
    f = join(dirpath, "27/label_threshold_12.6540002823.mol2")
else:
    f = "27/label_threshold_12.6540002823.mol2"

cmd.load(f, 'label_threshold_12.6540002823')
cmd.hide('everything', 'label_threshold_12.6540002823')
cmd.label("label_threshold_12.6540002823", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.6540002823', members= 'label_threshold_12.6540002823')


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
    f = join(dirpath, "28/label_threshold_7.5.mol2")
else:
    f = "28/label_threshold_7.5.mol2"

cmd.load(f, 'label_threshold_7.5')
cmd.hide('everything', 'label_threshold_7.5')
cmd.label("label_threshold_7.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [7.5]
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


cluster_dict = {"12.5010004044":[], "12.5010004044_arrows":[]}

cluster_dict["12.5010004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(43.5), float(-8.5), float(9.5), float(1.0)]

cluster_dict["12.5010004044_arrows"] += cgo_arrow([43.5,-8.5,9.5], [43.181,-9.676,7.212], color="blue red", name="Arrows_12.5010004044_1")

cluster_dict["12.5010004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(44.5), float(-0.5), float(11.5), float(1.0)]

cluster_dict["12.5010004044_arrows"] += cgo_arrow([44.5,-0.5,11.5], [47.549,0.315,10.678], color="blue red", name="Arrows_12.5010004044_2")

cluster_dict["12.5010004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(48.0), float(-6.0), float(8.0), float(1.0)]

cluster_dict["12.5010004044_arrows"] += cgo_arrow([48.0,-6.0,8.0], [45.592,-6.463,6.489], color="blue red", name="Arrows_12.5010004044_3")

cluster_dict["12.5010004044"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(44.1988681307), float(-5.93181655342), float(11.0588883007), float(1.0)]


cluster_dict["12.5010004044"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(43.5561487951), float(2.2103789445), float(10.9821651025), float(1.0)]


cluster_dict["12.5010004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(47.0), float(-7.0), float(9.0), float(1.0)]

cluster_dict["12.5010004044_arrows"] += cgo_arrow([47.0,-7.0,9.0], [45.592,-6.463,6.489], color="red blue", name="Arrows_12.5010004044_4")

cluster_dict["12.5010004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(45.0), float(-3.5), float(10.0), float(1.0)]

cluster_dict["12.5010004044_arrows"] += cgo_arrow([45.0,-3.5,10.0], [45.601,-1.014,8.286], color="red blue", name="Arrows_12.5010004044_5")

cluster_dict["12.5010004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(48.0), float(-6.5), float(7.0), float(1.0)]

cluster_dict["12.5010004044_arrows"] += cgo_arrow([48.0,-6.5,7.0], [45.592,-6.463,6.489], color="red blue", name="Arrows_12.5010004044_6")

cmd.load_cgo(cluster_dict["12.5010004044"], "Features_12.5010004044", 1)
cmd.load_cgo(cluster_dict["12.5010004044_arrows"], "Arrows_12.5010004044")
cmd.set("transparency", 0.2,"Features_12.5010004044")
cmd.group("Pharmacophore_12.5010004044", members="Features_12.5010004044")
cmd.group("Pharmacophore_12.5010004044", members="Arrows_12.5010004044")

if dirpath:
    f = join(dirpath, "28/label_threshold_12.5010004044.mol2")
else:
    f = "28/label_threshold_12.5010004044.mol2"

cmd.load(f, 'label_threshold_12.5010004044')
cmd.hide('everything', 'label_threshold_12.5010004044')
cmd.label("label_threshold_12.5010004044", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.5010004044', members= 'label_threshold_12.5010004044')


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
    f = join(dirpath, "29/label_threshold_9.8.mol2")
else:
    f = "29/label_threshold_9.8.mol2"

cmd.load(f, 'label_threshold_9.8')
cmd.hide('everything', 'label_threshold_9.8')
cmd.label("label_threshold_9.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [9.8]
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


cluster_dict = {"12.1820001602":[], "12.1820001602_arrows":[]}

cluster_dict["12.1820001602"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(-43.5), float(4.0), float(1.0)]

cluster_dict["12.1820001602_arrows"] += cgo_arrow([29.5,-43.5,4.0], [29.393,-40.692,4.576], color="blue red", name="Arrows_12.1820001602_1")

cluster_dict["12.1820001602"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(31.0), float(-42.5), float(3.5), float(1.0)]

cluster_dict["12.1820001602_arrows"] += cgo_arrow([31.0,-42.5,3.5], [29.393,-40.692,4.576], color="blue red", name="Arrows_12.1820001602_2")

cluster_dict["12.1820001602"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(29.5213145668), float(-45.246176176), float(6.29503206958), float(1.0)]


cluster_dict["12.1820001602"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(35.9991862931), float(-39.0421710118), float(0.669795276505), float(1.0)]


cluster_dict["12.1820001602"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(36.2664319983), float(-40.3510078108), float(7.7637118295), float(1.0)]


cluster_dict["12.1820001602"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(35.75), float(-39.5), float(10.5), float(1.0)]


cluster_dict["12.1820001602"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(33.5), float(-43.5), float(9.5), float(1.0)]

cluster_dict["12.1820001602_arrows"] += cgo_arrow([33.5,-43.5,9.5], [36.115,-42.175,10.275], color="red blue", name="Arrows_12.1820001602_3")

cmd.load_cgo(cluster_dict["12.1820001602"], "Features_12.1820001602", 1)
cmd.load_cgo(cluster_dict["12.1820001602_arrows"], "Arrows_12.1820001602")
cmd.set("transparency", 0.2,"Features_12.1820001602")
cmd.group("Pharmacophore_12.1820001602", members="Features_12.1820001602")
cmd.group("Pharmacophore_12.1820001602", members="Arrows_12.1820001602")

if dirpath:
    f = join(dirpath, "29/label_threshold_12.1820001602.mol2")
else:
    f = "29/label_threshold_12.1820001602.mol2"

cmd.load(f, 'label_threshold_12.1820001602')
cmd.hide('everything', 'label_threshold_12.1820001602')
cmd.label("label_threshold_12.1820001602", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.1820001602', members= 'label_threshold_12.1820001602')


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
    f = join(dirpath, "30/label_threshold_9.7.mol2")
else:
    f = "30/label_threshold_9.7.mol2"

cmd.load(f, 'label_threshold_9.7')
cmd.hide('everything', 'label_threshold_9.7')
cmd.label("label_threshold_9.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [9.7]
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


cluster_dict = {"12.1029996872":[], "12.1029996872_arrows":[]}

cluster_dict["12.1029996872"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(53.3301705154), float(-40.0778521353), float(18.7983058637), float(1.0)]


cluster_dict["12.1029996872"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(57.9041726053), float(-40.5769700325), float(10.8592893609), float(1.0)]


cluster_dict["12.1029996872"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(55.2487758074), float(-38.0048967703), float(16.0), float(1.0)]


cluster_dict["12.1029996872"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(55.5), float(-41.5), float(9.0), float(1.0)]

cluster_dict["12.1029996872_arrows"] += cgo_arrow([55.5,-41.5,9.0], [56.45,-42.109,6.174], color="red blue", name="Arrows_12.1029996872_1")

cluster_dict["12.1029996872"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(57.0), float(-42.5), float(8.5), float(1.0)]

cluster_dict["12.1029996872_arrows"] += cgo_arrow([57.0,-42.5,8.5], [56.45,-42.109,6.174], color="red blue", name="Arrows_12.1029996872_2")

cmd.load_cgo(cluster_dict["12.1029996872"], "Features_12.1029996872", 1)
cmd.load_cgo(cluster_dict["12.1029996872_arrows"], "Arrows_12.1029996872")
cmd.set("transparency", 0.2,"Features_12.1029996872")
cmd.group("Pharmacophore_12.1029996872", members="Features_12.1029996872")
cmd.group("Pharmacophore_12.1029996872", members="Arrows_12.1029996872")

if dirpath:
    f = join(dirpath, "30/label_threshold_12.1029996872.mol2")
else:
    f = "30/label_threshold_12.1029996872.mol2"

cmd.load(f, 'label_threshold_12.1029996872')
cmd.hide('everything', 'label_threshold_12.1029996872')
cmd.label("label_threshold_12.1029996872", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.1029996872', members= 'label_threshold_12.1029996872')


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
    f = join(dirpath, "31/label_threshold_9.5.mol2")
else:
    f = "31/label_threshold_9.5.mol2"

cmd.load(f, 'label_threshold_9.5')
cmd.hide('everything', 'label_threshold_9.5')
cmd.label("label_threshold_9.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [9.5]
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


cluster_dict = {"12.0954999924":[], "12.0954999924_arrows":[]}

cluster_dict["12.0954999924"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(-3.0), float(12.5), float(1.0)]

cluster_dict["12.0954999924_arrows"] += cgo_arrow([24.0,-3.0,12.5], [21.766,-1.319,13.37], color="blue red", name="Arrows_12.0954999924_1")

cluster_dict["12.0954999924"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(-9.0), float(14.0), float(1.0)]

cluster_dict["12.0954999924_arrows"] += cgo_arrow([26.5,-9.0,14.0], [28.783,-7.532,13.582], color="blue red", name="Arrows_12.0954999924_2")

cluster_dict["12.0954999924"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(2.0), float(2.5), float(1.0)]

cluster_dict["12.0954999924_arrows"] += cgo_arrow([26.5,2.0,2.5], [26.497,1.535,4.9], color="blue red", name="Arrows_12.0954999924_3")

cluster_dict["12.0954999924"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.0102040816), float(-2.20408163265), float(10.9285714286), float(1.0)]


cluster_dict["12.0954999924"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.6010706186), float(6.8955896563), float(17.269098198), float(1.0)]


cluster_dict["12.0954999924"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(30.002745416), float(0.859399399329), float(3.15462940441), float(1.0)]


cluster_dict["12.0954999924"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(26.0454545455), float(-3.54545454545), float(13.2727272727), float(1.0)]


cluster_dict["12.0954999924"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(34.0298684188), float(1.7691418965), float(12.186743167), float(1.0)]


cluster_dict["12.0954999924"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(25.5), float(0.0), float(12.0), float(1.0)]

cluster_dict["12.0954999924_arrows"] += cgo_arrow([25.5,0.0,12.0], [28.5,3.731,11.033], color="red blue", name="Arrows_12.0954999924_4")

cluster_dict["12.0954999924"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(-2.5), float(11.5), float(1.0)]

cluster_dict["12.0954999924_arrows"] += cgo_arrow([27.0,-2.5,11.5], [25.888,-5.209,9.814], color="red blue", name="Arrows_12.0954999924_5")

cluster_dict["12.0954999924"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(-2.5), float(4.0), float(1.0)]

cluster_dict["12.0954999924_arrows"] += cgo_arrow([26.0,-2.5,4.0], [26.993,-4.956,2.181], color="red blue", name="Arrows_12.0954999924_6")

cluster_dict["12.0954999924"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(31.5), float(-0.5), float(2.0), float(1.0)]

cluster_dict["12.0954999924_arrows"] += cgo_arrow([31.5,-0.5,2.0], [33.219,-2.042,2.779], color="red blue", name="Arrows_12.0954999924_7")

cmd.load_cgo(cluster_dict["12.0954999924"], "Features_12.0954999924", 1)
cmd.load_cgo(cluster_dict["12.0954999924_arrows"], "Arrows_12.0954999924")
cmd.set("transparency", 0.2,"Features_12.0954999924")
cmd.group("Pharmacophore_12.0954999924", members="Features_12.0954999924")
cmd.group("Pharmacophore_12.0954999924", members="Arrows_12.0954999924")

if dirpath:
    f = join(dirpath, "31/label_threshold_12.0954999924.mol2")
else:
    f = "31/label_threshold_12.0954999924.mol2"

cmd.load(f, 'label_threshold_12.0954999924')
cmd.hide('everything', 'label_threshold_12.0954999924')
cmd.label("label_threshold_12.0954999924", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.0954999924', members= 'label_threshold_12.0954999924')


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

if dirpath:
    f = join(dirpath, "32/label_threshold_8.0.mol2")
else:
    f = "32/label_threshold_8.0.mol2"

cmd.load(f, 'label_threshold_8.0')
cmd.hide('everything', 'label_threshold_8.0')
cmd.label("label_threshold_8.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [8.0]
gfiles = ['32/donor.grd', '32/apolar.grd', '32/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 32
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


cluster_dict = {"12.0860004425":[], "12.0860004425_arrows":[]}

cluster_dict["12.0860004425"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(41.0), float(17.0), float(3.0), float(1.0)]

cluster_dict["12.0860004425_arrows"] += cgo_arrow([41.0,17.0,3.0], [39.58,15.088,1.315], color="blue red", name="Arrows_12.0860004425_1")

cluster_dict["12.0860004425"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(46.5), float(11.5), float(10.5), float(1.0)]

cluster_dict["12.0860004425_arrows"] += cgo_arrow([46.5,11.5,10.5], [44.601,13.071,11.785], color="blue red", name="Arrows_12.0860004425_2")

cluster_dict["12.0860004425"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(50.0), float(11.0), float(9.5), float(1.0)]

cluster_dict["12.0860004425_arrows"] += cgo_arrow([50.0,11.0,9.5], [49.291,9.094,7.797], color="blue red", name="Arrows_12.0860004425_3")

cluster_dict["12.0860004425"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(49.5), float(12.0), float(7.5), float(1.0)]

cluster_dict["12.0860004425_arrows"] += cgo_arrow([49.5,12.0,7.5], [49.291,9.094,7.797], color="blue red", name="Arrows_12.0860004425_4")

cluster_dict["12.0860004425"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(40.0), float(16.6093263944), float(2.71865278872), float(1.0)]


cluster_dict["12.0860004425"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(45.8699230083), float(14.1376107004), float(6.41775495018), float(1.0)]


cmd.load_cgo(cluster_dict["12.0860004425"], "Features_12.0860004425", 1)
cmd.load_cgo(cluster_dict["12.0860004425_arrows"], "Arrows_12.0860004425")
cmd.set("transparency", 0.2,"Features_12.0860004425")
cmd.group("Pharmacophore_12.0860004425", members="Features_12.0860004425")
cmd.group("Pharmacophore_12.0860004425", members="Arrows_12.0860004425")

if dirpath:
    f = join(dirpath, "32/label_threshold_12.0860004425.mol2")
else:
    f = "32/label_threshold_12.0860004425.mol2"

cmd.load(f, 'label_threshold_12.0860004425')
cmd.hide('everything', 'label_threshold_12.0860004425')
cmd.label("label_threshold_12.0860004425", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.0860004425', members= 'label_threshold_12.0860004425')


if dirpath:
    f = join(dirpath, '32/mesh.grd')
else:
    f = '32/mesh.grd'
cmd.load(f, 'mesh_32')
cmd.isomesh("isomesh_32", "mesh_32", 0.9)
cmd.color("grey80", "isomesh_32")
cmd.set('transparency', 0.4, "isomesh_32")

cmd.group('hotspot_32', "isomesh_32")
cmd.group('hotspot_32', "mesh_32")

if dirpath:
    f = join(dirpath, "33/label_threshold_9.7.mol2")
else:
    f = "33/label_threshold_9.7.mol2"

cmd.load(f, 'label_threshold_9.7')
cmd.hide('everything', 'label_threshold_9.7')
cmd.label("label_threshold_9.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [9.7]
gfiles = ['33/donor.grd', '33/apolar.grd', '33/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 33
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


cluster_dict = {"12.0419998169":[], "12.0419998169_arrows":[]}

cluster_dict["12.0419998169"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(-43.5), float(4.0), float(1.0)]

cluster_dict["12.0419998169_arrows"] += cgo_arrow([29.5,-43.5,4.0], [29.393,-40.692,4.576], color="blue red", name="Arrows_12.0419998169_1")

cluster_dict["12.0419998169"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(31.0), float(-42.5), float(3.5), float(1.0)]

cluster_dict["12.0419998169_arrows"] += cgo_arrow([31.0,-42.5,3.5], [29.393,-40.692,4.576], color="blue red", name="Arrows_12.0419998169_2")

cluster_dict["12.0419998169"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(29.5283355315), float(-45.2513034402), float(6.28343107095), float(1.0)]


cluster_dict["12.0419998169"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(35.8041464737), float(-39.2312521331), float(0.75414185492), float(1.0)]


cluster_dict["12.0419998169"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(36.2664319983), float(-40.3510078108), float(7.7637118295), float(1.0)]


cluster_dict["12.0419998169"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(35.75), float(-39.5), float(10.5), float(1.0)]


cluster_dict["12.0419998169"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(33.5), float(-43.5), float(9.5), float(1.0)]

cluster_dict["12.0419998169_arrows"] += cgo_arrow([33.5,-43.5,9.5], [36.115,-42.175,10.275], color="red blue", name="Arrows_12.0419998169_3")

cmd.load_cgo(cluster_dict["12.0419998169"], "Features_12.0419998169", 1)
cmd.load_cgo(cluster_dict["12.0419998169_arrows"], "Arrows_12.0419998169")
cmd.set("transparency", 0.2,"Features_12.0419998169")
cmd.group("Pharmacophore_12.0419998169", members="Features_12.0419998169")
cmd.group("Pharmacophore_12.0419998169", members="Arrows_12.0419998169")

if dirpath:
    f = join(dirpath, "33/label_threshold_12.0419998169.mol2")
else:
    f = "33/label_threshold_12.0419998169.mol2"

cmd.load(f, 'label_threshold_12.0419998169')
cmd.hide('everything', 'label_threshold_12.0419998169')
cmd.label("label_threshold_12.0419998169", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.0419998169', members= 'label_threshold_12.0419998169')


if dirpath:
    f = join(dirpath, '33/mesh.grd')
else:
    f = '33/mesh.grd'
cmd.load(f, 'mesh_33')
cmd.isomesh("isomesh_33", "mesh_33", 0.9)
cmd.color("grey80", "isomesh_33")
cmd.set('transparency', 0.4, "isomesh_33")

cmd.group('hotspot_33', "isomesh_33")
cmd.group('hotspot_33', "mesh_33")

if dirpath:
    f = join(dirpath, "34/label_threshold_8.1.mol2")
else:
    f = "34/label_threshold_8.1.mol2"

cmd.load(f, 'label_threshold_8.1')
cmd.hide('everything', 'label_threshold_8.1')
cmd.label("label_threshold_8.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [8.1]
gfiles = ['34/donor.grd', '34/apolar.grd', '34/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 34
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


cluster_dict = {"12.029999733":[], "12.029999733_arrows":[]}

cluster_dict["12.029999733"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(-43.5), float(4.0), float(1.0)]

cluster_dict["12.029999733_arrows"] += cgo_arrow([29.5,-43.5,4.0], [29.393,-40.692,4.576], color="blue red", name="Arrows_12.029999733_1")

cluster_dict["12.029999733"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(31.0), float(-42.5), float(3.5), float(1.0)]

cluster_dict["12.029999733_arrows"] += cgo_arrow([31.0,-42.5,3.5], [29.393,-40.692,4.576], color="blue red", name="Arrows_12.029999733_2")

cluster_dict["12.029999733"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(30.4724924888), float(-44.9272981086), float(6.71223038159), float(1.0)]


cluster_dict["12.029999733"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(35.9794122665), float(-38.9848604213), float(0.802960305056), float(1.0)]


cluster_dict["12.029999733"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(36.2742313175), float(-40.3106254609), float(7.8193281673), float(1.0)]


cluster_dict["12.029999733"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(35.75), float(-39.5), float(10.5), float(1.0)]


cluster_dict["12.029999733"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(33.5), float(-43.5), float(9.5), float(1.0)]

cluster_dict["12.029999733_arrows"] += cgo_arrow([33.5,-43.5,9.5], [36.115,-42.175,10.275], color="red blue", name="Arrows_12.029999733_3")

cmd.load_cgo(cluster_dict["12.029999733"], "Features_12.029999733", 1)
cmd.load_cgo(cluster_dict["12.029999733_arrows"], "Arrows_12.029999733")
cmd.set("transparency", 0.2,"Features_12.029999733")
cmd.group("Pharmacophore_12.029999733", members="Features_12.029999733")
cmd.group("Pharmacophore_12.029999733", members="Arrows_12.029999733")

if dirpath:
    f = join(dirpath, "34/label_threshold_12.029999733.mol2")
else:
    f = "34/label_threshold_12.029999733.mol2"

cmd.load(f, 'label_threshold_12.029999733')
cmd.hide('everything', 'label_threshold_12.029999733')
cmd.label("label_threshold_12.029999733", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.029999733', members= 'label_threshold_12.029999733')


if dirpath:
    f = join(dirpath, '34/mesh.grd')
else:
    f = '34/mesh.grd'
cmd.load(f, 'mesh_34')
cmd.isomesh("isomesh_34", "mesh_34", 0.9)
cmd.color("grey80", "isomesh_34")
cmd.set('transparency', 0.4, "isomesh_34")

cmd.group('hotspot_34', "isomesh_34")
cmd.group('hotspot_34', "mesh_34")

if dirpath:
    f = join(dirpath, "35/label_threshold_8.3.mol2")
else:
    f = "35/label_threshold_8.3.mol2"

cmd.load(f, 'label_threshold_8.3')
cmd.hide('everything', 'label_threshold_8.3')
cmd.label("label_threshold_8.3", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [8.3]
gfiles = ['35/donor.grd', '35/apolar.grd', '35/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 35
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


cluster_dict = {"11.9879999161":[], "11.9879999161_arrows":[]}

cluster_dict["11.9879999161"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(36.4279743046), float(-25.1182023563), float(20.000089215), float(1.0)]


cluster_dict["11.9879999161"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(42.0232888331), float(-27.5142026781), float(17.9064148486), float(1.0)]


cluster_dict["11.9879999161"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(44.3952620346), float(-26.2457956523), float(21.0513184768), float(1.0)]


cluster_dict["11.9879999161"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(42.5), float(-29.5), float(14.0), float(1.0)]

cluster_dict["11.9879999161_arrows"] += cgo_arrow([42.5,-29.5,14.0], [45.812,-29.166,13.284], color="red blue", name="Arrows_11.9879999161_1")

cluster_dict["11.9879999161"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(42.5), float(-29.5), float(14.0), float(1.0)]

cluster_dict["11.9879999161_arrows"] += cgo_arrow([42.5,-29.5,14.0], [45.812,-29.166,13.284], color="red blue", name="Arrows_11.9879999161_2")

cmd.load_cgo(cluster_dict["11.9879999161"], "Features_11.9879999161", 1)
cmd.load_cgo(cluster_dict["11.9879999161_arrows"], "Arrows_11.9879999161")
cmd.set("transparency", 0.2,"Features_11.9879999161")
cmd.group("Pharmacophore_11.9879999161", members="Features_11.9879999161")
cmd.group("Pharmacophore_11.9879999161", members="Arrows_11.9879999161")

if dirpath:
    f = join(dirpath, "35/label_threshold_11.9879999161.mol2")
else:
    f = "35/label_threshold_11.9879999161.mol2"

cmd.load(f, 'label_threshold_11.9879999161')
cmd.hide('everything', 'label_threshold_11.9879999161')
cmd.label("label_threshold_11.9879999161", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.9879999161', members= 'label_threshold_11.9879999161')


if dirpath:
    f = join(dirpath, '35/mesh.grd')
else:
    f = '35/mesh.grd'
cmd.load(f, 'mesh_35')
cmd.isomesh("isomesh_35", "mesh_35", 0.9)
cmd.color("grey80", "isomesh_35")
cmd.set('transparency', 0.4, "isomesh_35")

cmd.group('hotspot_35', "isomesh_35")
cmd.group('hotspot_35', "mesh_35")

if dirpath:
    f = join(dirpath, "36/label_threshold_10.4.mol2")
else:
    f = "36/label_threshold_10.4.mol2"

cmd.load(f, 'label_threshold_10.4')
cmd.hide('everything', 'label_threshold_10.4')
cmd.label("label_threshold_10.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.4]
gfiles = ['36/donor.grd', '36/apolar.grd', '36/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 36
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


cluster_dict = {"11.7349996567":[], "11.7349996567_arrows":[]}

cluster_dict["11.7349996567"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(36.401580872), float(-25.1434731223), float(20.0782661172), float(1.0)]


cluster_dict["11.7349996567"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(37.6531829504), float(-32.3089441279), float(14.7192307483), float(1.0)]


cluster_dict["11.7349996567"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(46.0976829103), float(-30.1991480536), float(17.3162447442), float(1.0)]


cluster_dict["11.7349996567"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(44.4146203318), float(-22.6676985524), float(16.039699147), float(1.0)]


cluster_dict["11.7349996567"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(42.5), float(-29.5), float(14.0), float(1.0)]

cluster_dict["11.7349996567_arrows"] += cgo_arrow([42.5,-29.5,14.0], [45.812,-29.166,13.284], color="red blue", name="Arrows_11.7349996567_1")

cluster_dict["11.7349996567"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(42.5), float(-29.5), float(14.0), float(1.0)]

cluster_dict["11.7349996567_arrows"] += cgo_arrow([42.5,-29.5,14.0], [45.812,-29.166,13.284], color="red blue", name="Arrows_11.7349996567_2")

cmd.load_cgo(cluster_dict["11.7349996567"], "Features_11.7349996567", 1)
cmd.load_cgo(cluster_dict["11.7349996567_arrows"], "Arrows_11.7349996567")
cmd.set("transparency", 0.2,"Features_11.7349996567")
cmd.group("Pharmacophore_11.7349996567", members="Features_11.7349996567")
cmd.group("Pharmacophore_11.7349996567", members="Arrows_11.7349996567")

if dirpath:
    f = join(dirpath, "36/label_threshold_11.7349996567.mol2")
else:
    f = "36/label_threshold_11.7349996567.mol2"

cmd.load(f, 'label_threshold_11.7349996567')
cmd.hide('everything', 'label_threshold_11.7349996567')
cmd.label("label_threshold_11.7349996567", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.7349996567', members= 'label_threshold_11.7349996567')


if dirpath:
    f = join(dirpath, '36/mesh.grd')
else:
    f = '36/mesh.grd'
cmd.load(f, 'mesh_36')
cmd.isomesh("isomesh_36", "mesh_36", 0.9)
cmd.color("grey80", "isomesh_36")
cmd.set('transparency', 0.4, "isomesh_36")

cmd.group('hotspot_36', "isomesh_36")
cmd.group('hotspot_36', "mesh_36")

if dirpath:
    f = join(dirpath, "37/label_threshold_10.0.mol2")
else:
    f = "37/label_threshold_10.0.mol2"

cmd.load(f, 'label_threshold_10.0')
cmd.hide('everything', 'label_threshold_10.0')
cmd.label("label_threshold_10.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.0]
gfiles = ['37/donor.grd', '37/apolar.grd', '37/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 37
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


cluster_dict = {"11.6470003128":[], "11.6470003128_arrows":[]}

cluster_dict["11.6470003128"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(37.0264287201), float(-25.0877567914), float(20.0539717724), float(1.0)]


cluster_dict["11.6470003128"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(38.1231552129), float(-30.9098165759), float(15.8924255082), float(1.0)]


cluster_dict["11.6470003128"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(46.5227995208), float(-29.9727389681), float(17.2870759309), float(1.0)]


cluster_dict["11.6470003128"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(44.4146203318), float(-22.6676985524), float(16.039699147), float(1.0)]


cluster_dict["11.6470003128"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(42.5), float(-29.5), float(14.0), float(1.0)]

cluster_dict["11.6470003128_arrows"] += cgo_arrow([42.5,-29.5,14.0], [45.812,-29.166,13.284], color="red blue", name="Arrows_11.6470003128_1")

cluster_dict["11.6470003128"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(42.5), float(-29.5), float(14.0), float(1.0)]

cluster_dict["11.6470003128_arrows"] += cgo_arrow([42.5,-29.5,14.0], [45.812,-29.166,13.284], color="red blue", name="Arrows_11.6470003128_2")

cmd.load_cgo(cluster_dict["11.6470003128"], "Features_11.6470003128", 1)
cmd.load_cgo(cluster_dict["11.6470003128_arrows"], "Arrows_11.6470003128")
cmd.set("transparency", 0.2,"Features_11.6470003128")
cmd.group("Pharmacophore_11.6470003128", members="Features_11.6470003128")
cmd.group("Pharmacophore_11.6470003128", members="Arrows_11.6470003128")

if dirpath:
    f = join(dirpath, "37/label_threshold_11.6470003128.mol2")
else:
    f = "37/label_threshold_11.6470003128.mol2"

cmd.load(f, 'label_threshold_11.6470003128')
cmd.hide('everything', 'label_threshold_11.6470003128')
cmd.label("label_threshold_11.6470003128", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.6470003128', members= 'label_threshold_11.6470003128')


if dirpath:
    f = join(dirpath, '37/mesh.grd')
else:
    f = '37/mesh.grd'
cmd.load(f, 'mesh_37')
cmd.isomesh("isomesh_37", "mesh_37", 0.9)
cmd.color("grey80", "isomesh_37")
cmd.set('transparency', 0.4, "isomesh_37")

cmd.group('hotspot_37', "isomesh_37")
cmd.group('hotspot_37', "mesh_37")

if dirpath:
    f = join(dirpath, "38/label_threshold_9.6.mol2")
else:
    f = "38/label_threshold_9.6.mol2"

cmd.load(f, 'label_threshold_9.6')
cmd.hide('everything', 'label_threshold_9.6')
cmd.label("label_threshold_9.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [9.6]
gfiles = ['38/donor.grd', '38/apolar.grd', '38/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 38
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


cluster_dict = {"11.5900001526":[], "11.5900001526_arrows":[]}

cluster_dict["11.5900001526"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(45.7953508174), float(-29.4423513193), float(17.3461529304), float(1.0)]


cluster_dict["11.5900001526"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(37.3649189437), float(-25.4135755847), float(19.9874324103), float(1.0)]


cluster_dict["11.5900001526"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(44.3718462905), float(-22.6492874893), float(16.040255507), float(1.0)]


cluster_dict["11.5900001526"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(42.5), float(-29.5), float(14.0), float(1.0)]

cluster_dict["11.5900001526_arrows"] += cgo_arrow([42.5,-29.5,14.0], [45.812,-29.166,13.284], color="red blue", name="Arrows_11.5900001526_1")

cluster_dict["11.5900001526"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(42.5), float(-29.5), float(14.0), float(1.0)]

cluster_dict["11.5900001526_arrows"] += cgo_arrow([42.5,-29.5,14.0], [45.812,-29.166,13.284], color="red blue", name="Arrows_11.5900001526_2")

cmd.load_cgo(cluster_dict["11.5900001526"], "Features_11.5900001526", 1)
cmd.load_cgo(cluster_dict["11.5900001526_arrows"], "Arrows_11.5900001526")
cmd.set("transparency", 0.2,"Features_11.5900001526")
cmd.group("Pharmacophore_11.5900001526", members="Features_11.5900001526")
cmd.group("Pharmacophore_11.5900001526", members="Arrows_11.5900001526")

if dirpath:
    f = join(dirpath, "38/label_threshold_11.5900001526.mol2")
else:
    f = "38/label_threshold_11.5900001526.mol2"

cmd.load(f, 'label_threshold_11.5900001526')
cmd.hide('everything', 'label_threshold_11.5900001526')
cmd.label("label_threshold_11.5900001526", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.5900001526', members= 'label_threshold_11.5900001526')


if dirpath:
    f = join(dirpath, '38/mesh.grd')
else:
    f = '38/mesh.grd'
cmd.load(f, 'mesh_38')
cmd.isomesh("isomesh_38", "mesh_38", 0.9)
cmd.color("grey80", "isomesh_38")
cmd.set('transparency', 0.4, "isomesh_38")

cmd.group('hotspot_38', "isomesh_38")
cmd.group('hotspot_38', "mesh_38")

if dirpath:
    f = join(dirpath, "39/label_threshold_9.6.mol2")
else:
    f = "39/label_threshold_9.6.mol2"

cmd.load(f, 'label_threshold_9.6')
cmd.hide('everything', 'label_threshold_9.6')
cmd.label("label_threshold_9.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [9.6]
gfiles = ['39/donor.grd', '39/apolar.grd', '39/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 39
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


cluster_dict = {"11.5340003967":[], "11.5340003967_arrows":[]}

cluster_dict["11.5340003967"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(46.6385955723), float(-30.0080790032), float(17.2383721812), float(1.0)]


cluster_dict["11.5340003967"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(44.3718462905), float(-22.6492874893), float(16.040255507), float(1.0)]


cluster_dict["11.5340003967"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(42.5), float(-29.5), float(14.0), float(1.0)]

cluster_dict["11.5340003967_arrows"] += cgo_arrow([42.5,-29.5,14.0], [45.812,-29.166,13.284], color="red blue", name="Arrows_11.5340003967_1")

cluster_dict["11.5340003967"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(42.5), float(-29.5), float(14.0), float(1.0)]

cluster_dict["11.5340003967_arrows"] += cgo_arrow([42.5,-29.5,14.0], [45.812,-29.166,13.284], color="red blue", name="Arrows_11.5340003967_2")

cmd.load_cgo(cluster_dict["11.5340003967"], "Features_11.5340003967", 1)
cmd.load_cgo(cluster_dict["11.5340003967_arrows"], "Arrows_11.5340003967")
cmd.set("transparency", 0.2,"Features_11.5340003967")
cmd.group("Pharmacophore_11.5340003967", members="Features_11.5340003967")
cmd.group("Pharmacophore_11.5340003967", members="Arrows_11.5340003967")

if dirpath:
    f = join(dirpath, "39/label_threshold_11.5340003967.mol2")
else:
    f = "39/label_threshold_11.5340003967.mol2"

cmd.load(f, 'label_threshold_11.5340003967')
cmd.hide('everything', 'label_threshold_11.5340003967')
cmd.label("label_threshold_11.5340003967", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.5340003967', members= 'label_threshold_11.5340003967')


if dirpath:
    f = join(dirpath, '39/mesh.grd')
else:
    f = '39/mesh.grd'
cmd.load(f, 'mesh_39')
cmd.isomesh("isomesh_39", "mesh_39", 0.9)
cmd.color("grey80", "isomesh_39")
cmd.set('transparency', 0.4, "isomesh_39")

cmd.group('hotspot_39', "isomesh_39")
cmd.group('hotspot_39', "mesh_39")

if dirpath:
    f = join(dirpath, "40/label_threshold_8.2.mol2")
else:
    f = "40/label_threshold_8.2.mol2"

cmd.load(f, 'label_threshold_8.2')
cmd.hide('everything', 'label_threshold_8.2')
cmd.label("label_threshold_8.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [8.2]
gfiles = ['40/donor.grd', '40/apolar.grd', '40/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 40
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


cluster_dict = {"11.5059995651":[], "11.5059995651_arrows":[]}

cluster_dict["11.5059995651"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(-3.5), float(12.5), float(1.0)]

cluster_dict["11.5059995651_arrows"] += cgo_arrow([24.0,-3.5,12.5], [21.766,-1.319,13.37], color="blue red", name="Arrows_11.5059995651_1")

cluster_dict["11.5059995651"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(-3.5), float(12.5), float(1.0)]

cluster_dict["11.5059995651_arrows"] += cgo_arrow([24.0,-3.5,12.5], [21.766,-1.319,13.37], color="blue red", name="Arrows_11.5059995651_2")

cluster_dict["11.5059995651"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(-0.5), float(3.0), float(1.0)]

cluster_dict["11.5059995651_arrows"] += cgo_arrow([25.0,-0.5,3.0], [23.784,-1.609,0.528], color="blue red", name="Arrows_11.5059995651_3")

cluster_dict["11.5059995651"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(15.5), float(-15.5), float(18.0), float(1.0)]


cluster_dict["11.5059995651"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(18.184603147), float(-16.2256148549), float(16.1963207519), float(1.0)]


cluster_dict["11.5059995651"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.6164940474), float(-11.6013310024), float(9.14526445125), float(1.0)]


cluster_dict["11.5059995651"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(25.0957410712), float(-2.44218658832), float(11.9308281248), float(1.0)]


cluster_dict["11.5059995651"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.2531745281), float(-1.6884982441), float(4.99753601511), float(1.0)]


cluster_dict["11.5059995651"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(-13.0), float(15.0), float(1.0)]

cluster_dict["11.5059995651_arrows"] += cgo_arrow([18.0,-13.0,15.0], [16.48,-10.765,14.656], color="red blue", name="Arrows_11.5059995651_4")

cluster_dict["11.5059995651"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(-8.5), float(7.0), float(1.0)]

cluster_dict["11.5059995651_arrows"] += cgo_arrow([22.5,-8.5,7.0], [22.101,-6.363,4.763], color="red blue", name="Arrows_11.5059995651_5")

cluster_dict["11.5059995651"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(-2.5), float(11.5), float(1.0)]

cluster_dict["11.5059995651_arrows"] += cgo_arrow([26.5,-2.5,11.5], [25.888,-5.209,9.814], color="red blue", name="Arrows_11.5059995651_6")

cluster_dict["11.5059995651"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(-2.5), float(11.5), float(1.0)]

cluster_dict["11.5059995651_arrows"] += cgo_arrow([26.5,-2.5,11.5], [25.888,-5.209,9.814], color="red blue", name="Arrows_11.5059995651_7")

cluster_dict["11.5059995651"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(-2.5), float(4.0), float(1.0)]

cluster_dict["11.5059995651_arrows"] += cgo_arrow([26.0,-2.5,4.0], [26.993,-4.956,2.181], color="red blue", name="Arrows_11.5059995651_8")

cmd.load_cgo(cluster_dict["11.5059995651"], "Features_11.5059995651", 1)
cmd.load_cgo(cluster_dict["11.5059995651_arrows"], "Arrows_11.5059995651")
cmd.set("transparency", 0.2,"Features_11.5059995651")
cmd.group("Pharmacophore_11.5059995651", members="Features_11.5059995651")
cmd.group("Pharmacophore_11.5059995651", members="Arrows_11.5059995651")

if dirpath:
    f = join(dirpath, "40/label_threshold_11.5059995651.mol2")
else:
    f = "40/label_threshold_11.5059995651.mol2"

cmd.load(f, 'label_threshold_11.5059995651')
cmd.hide('everything', 'label_threshold_11.5059995651')
cmd.label("label_threshold_11.5059995651", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.5059995651', members= 'label_threshold_11.5059995651')


if dirpath:
    f = join(dirpath, '40/mesh.grd')
else:
    f = '40/mesh.grd'
cmd.load(f, 'mesh_40')
cmd.isomesh("isomesh_40", "mesh_40", 0.9)
cmd.color("grey80", "isomesh_40")
cmd.set('transparency', 0.4, "isomesh_40")

cmd.group('hotspot_40', "isomesh_40")
cmd.group('hotspot_40', "mesh_40")

if dirpath:
    f = join(dirpath, "41/label_threshold_7.6.mol2")
else:
    f = "41/label_threshold_7.6.mol2"

cmd.load(f, 'label_threshold_7.6')
cmd.hide('everything', 'label_threshold_7.6')
cmd.label("label_threshold_7.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [7.6]
gfiles = ['41/donor.grd', '41/apolar.grd', '41/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 41
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


cluster_dict = {"11.0439996719":[], "11.0439996719_arrows":[]}

cluster_dict["11.0439996719"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(-3.5), float(12.5), float(1.0)]

cluster_dict["11.0439996719_arrows"] += cgo_arrow([24.0,-3.5,12.5], [21.766,-1.319,13.37], color="blue red", name="Arrows_11.0439996719_1")

cluster_dict["11.0439996719"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(-3.5), float(12.5), float(1.0)]

cluster_dict["11.0439996719_arrows"] += cgo_arrow([24.0,-3.5,12.5], [21.766,-1.319,13.37], color="blue red", name="Arrows_11.0439996719_2")

cluster_dict["11.0439996719"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(18.0927916096), float(-16.1571336724), float(16.2778909622), float(1.0)]


cluster_dict["11.0439996719"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.6701222476), float(-11.5772672608), float(9.13904883518), float(1.0)]


cluster_dict["11.0439996719"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(25.3895035037), float(-2.59393837504), float(12.2428014913), float(1.0)]


cluster_dict["11.0439996719"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(28.1894508781), float(-2.48163051478), float(6.75780555983), float(1.0)]


cluster_dict["11.0439996719"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(-13.0), float(15.0), float(1.0)]

cluster_dict["11.0439996719_arrows"] += cgo_arrow([18.0,-13.0,15.0], [16.48,-10.765,14.656], color="red blue", name="Arrows_11.0439996719_3")

cluster_dict["11.0439996719"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(-8.5), float(7.0), float(1.0)]

cluster_dict["11.0439996719_arrows"] += cgo_arrow([22.5,-8.5,7.0], [22.101,-6.363,4.763], color="red blue", name="Arrows_11.0439996719_4")

cluster_dict["11.0439996719"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(-2.5), float(11.5), float(1.0)]

cluster_dict["11.0439996719_arrows"] += cgo_arrow([27.0,-2.5,11.5], [25.888,-5.209,9.814], color="red blue", name="Arrows_11.0439996719_5")

cluster_dict["11.0439996719"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(-2.5), float(11.5), float(1.0)]

cluster_dict["11.0439996719_arrows"] += cgo_arrow([27.0,-2.5,11.5], [25.888,-5.209,9.814], color="red blue", name="Arrows_11.0439996719_6")

cmd.load_cgo(cluster_dict["11.0439996719"], "Features_11.0439996719", 1)
cmd.load_cgo(cluster_dict["11.0439996719_arrows"], "Arrows_11.0439996719")
cmd.set("transparency", 0.2,"Features_11.0439996719")
cmd.group("Pharmacophore_11.0439996719", members="Features_11.0439996719")
cmd.group("Pharmacophore_11.0439996719", members="Arrows_11.0439996719")

if dirpath:
    f = join(dirpath, "41/label_threshold_11.0439996719.mol2")
else:
    f = "41/label_threshold_11.0439996719.mol2"

cmd.load(f, 'label_threshold_11.0439996719')
cmd.hide('everything', 'label_threshold_11.0439996719')
cmd.label("label_threshold_11.0439996719", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.0439996719', members= 'label_threshold_11.0439996719')


if dirpath:
    f = join(dirpath, '41/mesh.grd')
else:
    f = '41/mesh.grd'
cmd.load(f, 'mesh_41')
cmd.isomesh("isomesh_41", "mesh_41", 0.9)
cmd.color("grey80", "isomesh_41")
cmd.set('transparency', 0.4, "isomesh_41")

cmd.group('hotspot_41', "isomesh_41")
cmd.group('hotspot_41', "mesh_41")

if dirpath:
    f = join(dirpath, "42/label_threshold_6.8.mol2")
else:
    f = "42/label_threshold_6.8.mol2"

cmd.load(f, 'label_threshold_6.8')
cmd.hide('everything', 'label_threshold_6.8')
cmd.label("label_threshold_6.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [6.8]
gfiles = ['42/donor.grd', '42/apolar.grd', '42/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 42
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


cluster_dict = {"10.876999855":[], "10.876999855_arrows":[]}

cluster_dict["10.876999855"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(35.5), float(-35.5), float(-1.0), float(1.0)]

cluster_dict["10.876999855_arrows"] += cgo_arrow([35.5,-35.5,-1.0], [33.351,-35.839,-4.065], color="blue red", name="Arrows_10.876999855_1")

cluster_dict["10.876999855"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(35.5), float(-35.5), float(-1.0), float(1.0)]

cluster_dict["10.876999855_arrows"] += cgo_arrow([35.5,-35.5,-1.0], [33.351,-35.839,-4.065], color="blue red", name="Arrows_10.876999855_2")

cluster_dict["10.876999855"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(36.2146059755), float(-37.8309649601), float(-2.11954347092), float(1.0)]


cluster_dict["10.876999855"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(40.6999787531), float(-28.2181823939), float(-3.47437827384), float(1.0)]


cluster_dict["10.876999855"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(43.7326540123), float(-37.6326600629), float(-9.28176599814), float(1.0)]


cluster_dict["10.876999855"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(-29.0), float(-5.0), float(1.0)]

cluster_dict["10.876999855_arrows"] += cgo_arrow([36.5,-29.0,-5.0], [37.132,-30.946,-6.482], color="red blue", name="Arrows_10.876999855_3")

cmd.load_cgo(cluster_dict["10.876999855"], "Features_10.876999855", 1)
cmd.load_cgo(cluster_dict["10.876999855_arrows"], "Arrows_10.876999855")
cmd.set("transparency", 0.2,"Features_10.876999855")
cmd.group("Pharmacophore_10.876999855", members="Features_10.876999855")
cmd.group("Pharmacophore_10.876999855", members="Arrows_10.876999855")

if dirpath:
    f = join(dirpath, "42/label_threshold_10.876999855.mol2")
else:
    f = "42/label_threshold_10.876999855.mol2"

cmd.load(f, 'label_threshold_10.876999855')
cmd.hide('everything', 'label_threshold_10.876999855')
cmd.label("label_threshold_10.876999855", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.876999855', members= 'label_threshold_10.876999855')


if dirpath:
    f = join(dirpath, '42/mesh.grd')
else:
    f = '42/mesh.grd'
cmd.load(f, 'mesh_42')
cmd.isomesh("isomesh_42", "mesh_42", 0.9)
cmd.color("grey80", "isomesh_42")
cmd.set('transparency', 0.4, "isomesh_42")

cmd.group('hotspot_42', "isomesh_42")
cmd.group('hotspot_42', "mesh_42")

if dirpath:
    f = join(dirpath, "43/label_threshold_4.4.mol2")
else:
    f = "43/label_threshold_4.4.mol2"

cmd.load(f, 'label_threshold_4.4')
cmd.hide('everything', 'label_threshold_4.4')
cmd.label("label_threshold_4.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [4.4]
gfiles = ['43/donor.grd', '43/apolar.grd', '43/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 43
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


cluster_dict = {"9.98200035095":[], "9.98200035095_arrows":[]}

cluster_dict["9.98200035095"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(-9.5), float(14.0), float(1.0)]

cluster_dict["9.98200035095_arrows"] += cgo_arrow([26.5,-9.5,14.0], [28.783,-7.532,13.582], color="blue red", name="Arrows_9.98200035095_1")

cluster_dict["9.98200035095"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(18.193311252), float(-16.1239391994), float(16.2652059336), float(1.0)]


cluster_dict["9.98200035095"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.56770901), float(-12.0466505801), float(10.3365330733), float(1.0)]


cluster_dict["9.98200035095"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(-13.0), float(15.0), float(1.0)]

cluster_dict["9.98200035095_arrows"] += cgo_arrow([18.0,-13.0,15.0], [16.48,-10.765,14.656], color="red blue", name="Arrows_9.98200035095_2")

cluster_dict["9.98200035095"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(-9.0), float(7.0), float(1.0)]

cluster_dict["9.98200035095_arrows"] += cgo_arrow([22.0,-9.0,7.0], [22.101,-6.363,4.763], color="red blue", name="Arrows_9.98200035095_3")

cmd.load_cgo(cluster_dict["9.98200035095"], "Features_9.98200035095", 1)
cmd.load_cgo(cluster_dict["9.98200035095_arrows"], "Arrows_9.98200035095")
cmd.set("transparency", 0.2,"Features_9.98200035095")
cmd.group("Pharmacophore_9.98200035095", members="Features_9.98200035095")
cmd.group("Pharmacophore_9.98200035095", members="Arrows_9.98200035095")

if dirpath:
    f = join(dirpath, "43/label_threshold_9.98200035095.mol2")
else:
    f = "43/label_threshold_9.98200035095.mol2"

cmd.load(f, 'label_threshold_9.98200035095')
cmd.hide('everything', 'label_threshold_9.98200035095')
cmd.label("label_threshold_9.98200035095", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.98200035095', members= 'label_threshold_9.98200035095')


if dirpath:
    f = join(dirpath, '43/mesh.grd')
else:
    f = '43/mesh.grd'
cmd.load(f, 'mesh_43')
cmd.isomesh("isomesh_43", "mesh_43", 0.9)
cmd.color("grey80", "isomesh_43")
cmd.set('transparency', 0.4, "isomesh_43")

cmd.group('hotspot_43', "isomesh_43")
cmd.group('hotspot_43', "mesh_43")

if dirpath:
    f = join(dirpath, "44/label_threshold_0.6.mol2")
else:
    f = "44/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
gfiles = ['44/donor.grd', '44/apolar.grd', '44/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 44
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


cluster_dict = {"9.40950012207":[], "9.40950012207_arrows":[]}

cluster_dict["9.40950012207"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-9.0), float(-46.5), float(22.0), float(1.0)]

cluster_dict["9.40950012207_arrows"] += cgo_arrow([-9.0,-46.5,22.0], [-10.053,-48.203,20.0], color="blue red", name="Arrows_9.40950012207_1")

cluster_dict["9.40950012207"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-10.5104567505), float(-43.6616399276), float(21.8990304042), float(1.0)]


cmd.load_cgo(cluster_dict["9.40950012207"], "Features_9.40950012207", 1)
cmd.load_cgo(cluster_dict["9.40950012207_arrows"], "Arrows_9.40950012207")
cmd.set("transparency", 0.2,"Features_9.40950012207")
cmd.group("Pharmacophore_9.40950012207", members="Features_9.40950012207")
cmd.group("Pharmacophore_9.40950012207", members="Arrows_9.40950012207")

if dirpath:
    f = join(dirpath, "44/label_threshold_9.40950012207.mol2")
else:
    f = "44/label_threshold_9.40950012207.mol2"

cmd.load(f, 'label_threshold_9.40950012207')
cmd.hide('everything', 'label_threshold_9.40950012207')
cmd.label("label_threshold_9.40950012207", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.40950012207', members= 'label_threshold_9.40950012207')


if dirpath:
    f = join(dirpath, '44/mesh.grd')
else:
    f = '44/mesh.grd'
cmd.load(f, 'mesh_44')
cmd.isomesh("isomesh_44", "mesh_44", 0.9)
cmd.color("grey80", "isomesh_44")
cmd.set('transparency', 0.4, "isomesh_44")

cmd.group('hotspot_44', "isomesh_44")
cmd.group('hotspot_44', "mesh_44")

if dirpath:
    f = join(dirpath, "45/label_threshold_0.6.mol2")
else:
    f = "45/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
gfiles = ['45/donor.grd', '45/apolar.grd', '45/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 45
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


cluster_dict = {"9.38199996948":[], "9.38199996948_arrows":[]}

cluster_dict["9.38199996948"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(11.4730523383), float(-36.7460288906), float(35.4814821071), float(1.0)]


cluster_dict["9.38199996948"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(11.0), float(-39.0), float(35.5), float(1.0)]

cluster_dict["9.38199996948_arrows"] += cgo_arrow([11.0,-39.0,35.5], [11.718,-39.152,38.215], color="red blue", name="Arrows_9.38199996948_1")

cmd.load_cgo(cluster_dict["9.38199996948"], "Features_9.38199996948", 1)
cmd.load_cgo(cluster_dict["9.38199996948_arrows"], "Arrows_9.38199996948")
cmd.set("transparency", 0.2,"Features_9.38199996948")
cmd.group("Pharmacophore_9.38199996948", members="Features_9.38199996948")
cmd.group("Pharmacophore_9.38199996948", members="Arrows_9.38199996948")

if dirpath:
    f = join(dirpath, "45/label_threshold_9.38199996948.mol2")
else:
    f = "45/label_threshold_9.38199996948.mol2"

cmd.load(f, 'label_threshold_9.38199996948')
cmd.hide('everything', 'label_threshold_9.38199996948')
cmd.label("label_threshold_9.38199996948", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.38199996948', members= 'label_threshold_9.38199996948')


if dirpath:
    f = join(dirpath, '45/mesh.grd')
else:
    f = '45/mesh.grd'
cmd.load(f, 'mesh_45')
cmd.isomesh("isomesh_45", "mesh_45", 0.9)
cmd.color("grey80", "isomesh_45")
cmd.set('transparency', 0.4, "isomesh_45")

cmd.group('hotspot_45', "isomesh_45")
cmd.group('hotspot_45', "mesh_45")

if dirpath:
    f = join(dirpath, "46/label_threshold_4.3.mol2")
else:
    f = "46/label_threshold_4.3.mol2"

cmd.load(f, 'label_threshold_4.3')
cmd.hide('everything', 'label_threshold_4.3')
cmd.label("label_threshold_4.3", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [4.3]
gfiles = ['46/donor.grd', '46/apolar.grd', '46/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 46
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


cluster_dict = {"8.93799972534":[], "8.93799972534_arrows":[]}

cluster_dict["8.93799972534"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(52.5), float(-7.5), float(6.5), float(1.0)]

cluster_dict["8.93799972534_arrows"] += cgo_arrow([52.5,-7.5,6.5], [50.139,-8.382,8.294], color="blue red", name="Arrows_8.93799972534_1")

cluster_dict["8.93799972534"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(46.5), float(-23.2757988631), float(-1.39680454755), float(1.0)]


cluster_dict["8.93799972534"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(52.5708664333), float(-12.3399113307), float(3.29011214206), float(1.0)]


cluster_dict["8.93799972534"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(49.0), float(-9.0), float(1.5), float(1.0)]

cluster_dict["8.93799972534_arrows"] += cgo_arrow([49.0,-9.0,1.5], [48.96,-11.301,1.715], color="red blue", name="Arrows_8.93799972534_2")

cluster_dict["8.93799972534"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(51.0), float(-9.0), float(3.5), float(1.0)]

cluster_dict["8.93799972534_arrows"] += cgo_arrow([51.0,-9.0,3.5], [48.96,-11.301,1.715], color="red blue", name="Arrows_8.93799972534_3")

cluster_dict["8.93799972534"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(51.0), float(-9.0), float(3.5), float(1.0)]

cluster_dict["8.93799972534_arrows"] += cgo_arrow([51.0,-9.0,3.5], [48.96,-11.301,1.715], color="red blue", name="Arrows_8.93799972534_4")

cmd.load_cgo(cluster_dict["8.93799972534"], "Features_8.93799972534", 1)
cmd.load_cgo(cluster_dict["8.93799972534_arrows"], "Arrows_8.93799972534")
cmd.set("transparency", 0.2,"Features_8.93799972534")
cmd.group("Pharmacophore_8.93799972534", members="Features_8.93799972534")
cmd.group("Pharmacophore_8.93799972534", members="Arrows_8.93799972534")

if dirpath:
    f = join(dirpath, "46/label_threshold_8.93799972534.mol2")
else:
    f = "46/label_threshold_8.93799972534.mol2"

cmd.load(f, 'label_threshold_8.93799972534')
cmd.hide('everything', 'label_threshold_8.93799972534')
cmd.label("label_threshold_8.93799972534", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_8.93799972534', members= 'label_threshold_8.93799972534')


if dirpath:
    f = join(dirpath, '46/mesh.grd')
else:
    f = '46/mesh.grd'
cmd.load(f, 'mesh_46')
cmd.isomesh("isomesh_46", "mesh_46", 0.9)
cmd.color("grey80", "isomesh_46")
cmd.set('transparency', 0.4, "isomesh_46")

cmd.group('hotspot_46', "isomesh_46")
cmd.group('hotspot_46', "mesh_46")

if dirpath:
    f = join(dirpath, "47/label_threshold_0.6.mol2")
else:
    f = "47/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
gfiles = ['47/donor.grd', '47/apolar.grd', '47/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 47
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


cluster_dict = {"8.61400032043":[], "8.61400032043_arrows":[]}

cluster_dict["8.61400032043"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(16.9853224264), float(-30.4380859215), float(14.5), float(1.0)]


cluster_dict["8.61400032043"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(20.8122773335), float(-27.696702135), float(7.99495805385), float(1.0)]


cluster_dict["8.61400032043"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(23.1658746464), float(-27.9996146955), float(7.0), float(1.0)]


cmd.load_cgo(cluster_dict["8.61400032043"], "Features_8.61400032043", 1)
cmd.load_cgo(cluster_dict["8.61400032043_arrows"], "Arrows_8.61400032043")
cmd.set("transparency", 0.2,"Features_8.61400032043")
cmd.group("Pharmacophore_8.61400032043", members="Features_8.61400032043")
cmd.group("Pharmacophore_8.61400032043", members="Arrows_8.61400032043")

if dirpath:
    f = join(dirpath, "47/label_threshold_8.61400032043.mol2")
else:
    f = "47/label_threshold_8.61400032043.mol2"

cmd.load(f, 'label_threshold_8.61400032043')
cmd.hide('everything', 'label_threshold_8.61400032043')
cmd.label("label_threshold_8.61400032043", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_8.61400032043', members= 'label_threshold_8.61400032043')


if dirpath:
    f = join(dirpath, '47/mesh.grd')
else:
    f = '47/mesh.grd'
cmd.load(f, 'mesh_47')
cmd.isomesh("isomesh_47", "mesh_47", 0.9)
cmd.color("grey80", "isomesh_47")
cmd.set('transparency', 0.4, "isomesh_47")

cmd.group('hotspot_47', "isomesh_47")
cmd.group('hotspot_47', "mesh_47")

if dirpath:
    f = join(dirpath, "48/label_threshold_1.9.mol2")
else:
    f = "48/label_threshold_1.9.mol2"

cmd.load(f, 'label_threshold_1.9')
cmd.hide('everything', 'label_threshold_1.9')
cmd.label("label_threshold_1.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [1.9]
gfiles = ['48/donor.grd', '48/apolar.grd', '48/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 48
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


cluster_dict = {"8.4359998703":[], "8.4359998703_arrows":[]}

cluster_dict["8.4359998703"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(6.5), float(-8.0), float(23.5), float(1.0)]

cluster_dict["8.4359998703_arrows"] += cgo_arrow([6.5,-8.0,23.5], [5.933,-11.235,24.279], color="blue red", name="Arrows_8.4359998703_1")

cluster_dict["8.4359998703"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(7.09453250767), float(-10.7128155245), float(21.1006044846), float(1.0)]


cluster_dict["8.4359998703"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(9.16639106182), float(0.243903419536), float(18.9119091002), float(1.0)]


cmd.load_cgo(cluster_dict["8.4359998703"], "Features_8.4359998703", 1)
cmd.load_cgo(cluster_dict["8.4359998703_arrows"], "Arrows_8.4359998703")
cmd.set("transparency", 0.2,"Features_8.4359998703")
cmd.group("Pharmacophore_8.4359998703", members="Features_8.4359998703")
cmd.group("Pharmacophore_8.4359998703", members="Arrows_8.4359998703")

if dirpath:
    f = join(dirpath, "48/label_threshold_8.4359998703.mol2")
else:
    f = "48/label_threshold_8.4359998703.mol2"

cmd.load(f, 'label_threshold_8.4359998703')
cmd.hide('everything', 'label_threshold_8.4359998703')
cmd.label("label_threshold_8.4359998703", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_8.4359998703', members= 'label_threshold_8.4359998703')


if dirpath:
    f = join(dirpath, '48/mesh.grd')
else:
    f = '48/mesh.grd'
cmd.load(f, 'mesh_48')
cmd.isomesh("isomesh_48", "mesh_48", 0.9)
cmd.color("grey80", "isomesh_48")
cmd.set('transparency', 0.4, "isomesh_48")

cmd.group('hotspot_48', "isomesh_48")
cmd.group('hotspot_48', "mesh_48")

if dirpath:
    f = join(dirpath, "49/label_threshold_0.6.mol2")
else:
    f = "49/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
gfiles = ['49/donor.grd', '49/apolar.grd', '49/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 49
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


cluster_dict = {"8.22900009155":[], "8.22900009155_arrows":[]}

cluster_dict["8.22900009155"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-2.36140885888), float(-0.0542945881066), float(19.2690039029), float(1.0)]


cmd.load_cgo(cluster_dict["8.22900009155"], "Features_8.22900009155", 1)
cmd.load_cgo(cluster_dict["8.22900009155_arrows"], "Arrows_8.22900009155")
cmd.set("transparency", 0.2,"Features_8.22900009155")
cmd.group("Pharmacophore_8.22900009155", members="Features_8.22900009155")
cmd.group("Pharmacophore_8.22900009155", members="Arrows_8.22900009155")

if dirpath:
    f = join(dirpath, "49/label_threshold_8.22900009155.mol2")
else:
    f = "49/label_threshold_8.22900009155.mol2"

cmd.load(f, 'label_threshold_8.22900009155')
cmd.hide('everything', 'label_threshold_8.22900009155')
cmd.label("label_threshold_8.22900009155", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_8.22900009155', members= 'label_threshold_8.22900009155')


if dirpath:
    f = join(dirpath, '49/mesh.grd')
else:
    f = '49/mesh.grd'
cmd.load(f, 'mesh_49')
cmd.isomesh("isomesh_49", "mesh_49", 0.9)
cmd.color("grey80", "isomesh_49")
cmd.set('transparency', 0.4, "isomesh_49")

cmd.group('hotspot_49', "isomesh_49")
cmd.group('hotspot_49', "mesh_49")

if dirpath:
    f = join(dirpath, "50/label_threshold_2.4.mol2")
else:
    f = "50/label_threshold_2.4.mol2"

cmd.load(f, 'label_threshold_2.4')
cmd.hide('everything', 'label_threshold_2.4')
cmd.label("label_threshold_2.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [2.4]
gfiles = ['50/donor.grd', '50/apolar.grd', '50/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 50
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


cluster_dict = {"7.08400011063":[], "7.08400011063_arrows":[]}

cluster_dict["7.08400011063"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(6.5), float(-8.0), float(23.5), float(1.0)]

cluster_dict["7.08400011063_arrows"] += cgo_arrow([6.5,-8.0,23.5], [5.933,-11.235,24.279], color="blue red", name="Arrows_7.08400011063_1")

cluster_dict["7.08400011063"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(13.5), float(0.0), float(17.5), float(1.0)]

cluster_dict["7.08400011063_arrows"] += cgo_arrow([13.5,0.0,17.5], [16.796,0.404,18.339], color="blue red", name="Arrows_7.08400011063_2")

cluster_dict["7.08400011063"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(7.42110694909), float(-8.36522876583), float(20.6666238093), float(1.0)]


cluster_dict["7.08400011063"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(9.43102063205), float(2.29461275388), float(18.9120245822), float(1.0)]


cmd.load_cgo(cluster_dict["7.08400011063"], "Features_7.08400011063", 1)
cmd.load_cgo(cluster_dict["7.08400011063_arrows"], "Arrows_7.08400011063")
cmd.set("transparency", 0.2,"Features_7.08400011063")
cmd.group("Pharmacophore_7.08400011063", members="Features_7.08400011063")
cmd.group("Pharmacophore_7.08400011063", members="Arrows_7.08400011063")

if dirpath:
    f = join(dirpath, "50/label_threshold_7.08400011063.mol2")
else:
    f = "50/label_threshold_7.08400011063.mol2"

cmd.load(f, 'label_threshold_7.08400011063')
cmd.hide('everything', 'label_threshold_7.08400011063')
cmd.label("label_threshold_7.08400011063", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_7.08400011063', members= 'label_threshold_7.08400011063')


if dirpath:
    f = join(dirpath, '50/mesh.grd')
else:
    f = '50/mesh.grd'
cmd.load(f, 'mesh_50')
cmd.isomesh("isomesh_50", "mesh_50", 0.9)
cmd.color("grey80", "isomesh_50")
cmd.set('transparency', 0.4, "isomesh_50")

cmd.group('hotspot_50', "isomesh_50")
cmd.group('hotspot_50', "mesh_50")

if dirpath:
    f = join(dirpath, "51/label_threshold_0.6.mol2")
else:
    f = "51/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
gfiles = ['51/donor.grd', '51/apolar.grd', '51/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 51
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


cluster_dict = {"6.70699977875":[], "6.70699977875_arrows":[]}

cluster_dict["6.70699977875"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.7743956738), float(6.83684400291), float(19.0774554209), float(1.0)]


cluster_dict["6.70699977875"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(26.7228391938), float(9.36550133575), float(17.8757056845), float(1.0)]


cmd.load_cgo(cluster_dict["6.70699977875"], "Features_6.70699977875", 1)
cmd.load_cgo(cluster_dict["6.70699977875_arrows"], "Arrows_6.70699977875")
cmd.set("transparency", 0.2,"Features_6.70699977875")
cmd.group("Pharmacophore_6.70699977875", members="Features_6.70699977875")
cmd.group("Pharmacophore_6.70699977875", members="Arrows_6.70699977875")

if dirpath:
    f = join(dirpath, "51/label_threshold_6.70699977875.mol2")
else:
    f = "51/label_threshold_6.70699977875.mol2"

cmd.load(f, 'label_threshold_6.70699977875')
cmd.hide('everything', 'label_threshold_6.70699977875')
cmd.label("label_threshold_6.70699977875", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_6.70699977875', members= 'label_threshold_6.70699977875')


if dirpath:
    f = join(dirpath, '51/mesh.grd')
else:
    f = '51/mesh.grd'
cmd.load(f, 'mesh_51')
cmd.isomesh("isomesh_51", "mesh_51", 0.9)
cmd.color("grey80", "isomesh_51")
cmd.set('transparency', 0.4, "isomesh_51")

cmd.group('hotspot_51', "isomesh_51")
cmd.group('hotspot_51', "mesh_51")

if dirpath:
    f = join(dirpath, "52/label_threshold_30.0.mol2")
else:
    f = "52/label_threshold_30.0.mol2"

cmd.load(f, 'label_threshold_30.0')
cmd.hide('everything', 'label_threshold_30.0')
cmd.label("label_threshold_30.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [30.0]
gfiles = ['52/donor.grd', '52/apolar.grd', '52/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 52
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
    f = join(dirpath, "52/label_threshold_0.mol2")
else:
    f = "52/label_threshold_0.mol2"

cmd.load(f, 'label_threshold_0')
cmd.hide('everything', 'label_threshold_0')
cmd.label("label_threshold_0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_0', members= 'label_threshold_0')


if dirpath:
    f = join(dirpath, '52/mesh.grd')
else:
    f = '52/mesh.grd'
cmd.load(f, 'mesh_52')
cmd.isomesh("isomesh_52", "mesh_52", 0.9)
cmd.color("grey80", "isomesh_52")
cmd.set('transparency', 0.4, "isomesh_52")

cmd.group('hotspot_52', "isomesh_52")
cmd.group('hotspot_52', "mesh_52")

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
