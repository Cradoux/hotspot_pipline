
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
    f = join(dirpath, "0/label_threshold_10.9.mol2")
else:
    f = "0/label_threshold_10.9.mol2"

cmd.load(f, 'label_threshold_10.9')
cmd.hide('everything', 'label_threshold_10.9')
cmd.label("label_threshold_10.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.9]
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


cluster_dict = {"13.3310003281":[], "13.3310003281_arrows":[]}

cluster_dict["13.3310003281"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-11.3831813907), float(31.5318793094), float(51.121979609), float(1.0)]


cluster_dict["13.3310003281"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-14.5), float(33.0), float(54.5), float(1.0)]

cluster_dict["13.3310003281_arrows"] += cgo_arrow([-14.5,33.0,54.5], [-15.963,34.648,56.255], color="red blue", name="Arrows_13.3310003281_1")

cluster_dict["13.3310003281"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-11.5), float(33.5), float(54.5), float(1.0)]

cluster_dict["13.3310003281_arrows"] += cgo_arrow([-11.5,33.5,54.5], [-8.203,33.346,53.721], color="red blue", name="Arrows_13.3310003281_2")

cluster_dict["13.3310003281"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-10.5), float(29.5), float(54.5), float(1.0)]

cluster_dict["13.3310003281_arrows"] += cgo_arrow([-10.5,29.5,54.5], [-7.337,28.877,55.918], color="red blue", name="Arrows_13.3310003281_3")

cluster_dict["13.3310003281"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-8.5), float(29.0), float(53.0), float(1.0)]

cluster_dict["13.3310003281_arrows"] += cgo_arrow([-8.5,29.0,53.0], [-7.337,28.877,55.918], color="red blue", name="Arrows_13.3310003281_4")

cluster_dict["13.3310003281"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-7.5), float(32.0), float(49.5), float(1.0)]

cluster_dict["13.3310003281_arrows"] += cgo_arrow([-7.5,32.0,49.5], [-6.6,34.578,50.608], color="red blue", name="Arrows_13.3310003281_5")

cmd.load_cgo(cluster_dict["13.3310003281"], "Features_13.3310003281", 1)
cmd.load_cgo(cluster_dict["13.3310003281_arrows"], "Arrows_13.3310003281")
cmd.set("transparency", 0.2,"Features_13.3310003281")
cmd.group("Pharmacophore_13.3310003281", members="Features_13.3310003281")
cmd.group("Pharmacophore_13.3310003281", members="Arrows_13.3310003281")

if dirpath:
    f = join(dirpath, "0/label_threshold_13.3310003281.mol2")
else:
    f = "0/label_threshold_13.3310003281.mol2"

cmd.load(f, 'label_threshold_13.3310003281')
cmd.hide('everything', 'label_threshold_13.3310003281')
cmd.label("label_threshold_13.3310003281", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.3310003281', members= 'label_threshold_13.3310003281')


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
    f = join(dirpath, "1/label_threshold_0.6.mol2")
else:
    f = "1/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
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


cluster_dict = {"11.5965003967":[], "11.5965003967_arrows":[]}

cluster_dict["11.5965003967"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-28.5), float(22.5), float(66.5), float(1.0)]

cluster_dict["11.5965003967_arrows"] += cgo_arrow([-28.5,22.5,66.5], [-30.793,24.547,65.715], color="blue red", name="Arrows_11.5965003967_1")

cluster_dict["11.5965003967"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-28.383880056), float(22.0952669), float(69.5792271123), float(1.0)]


cmd.load_cgo(cluster_dict["11.5965003967"], "Features_11.5965003967", 1)
cmd.load_cgo(cluster_dict["11.5965003967_arrows"], "Arrows_11.5965003967")
cmd.set("transparency", 0.2,"Features_11.5965003967")
cmd.group("Pharmacophore_11.5965003967", members="Features_11.5965003967")
cmd.group("Pharmacophore_11.5965003967", members="Arrows_11.5965003967")

if dirpath:
    f = join(dirpath, "1/label_threshold_11.5965003967.mol2")
else:
    f = "1/label_threshold_11.5965003967.mol2"

cmd.load(f, 'label_threshold_11.5965003967')
cmd.hide('everything', 'label_threshold_11.5965003967')
cmd.label("label_threshold_11.5965003967", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.5965003967', members= 'label_threshold_11.5965003967')


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
    f = join(dirpath, "2/label_threshold_4.2.mol2")
else:
    f = "2/label_threshold_4.2.mol2"

cmd.load(f, 'label_threshold_4.2')
cmd.hide('everything', 'label_threshold_4.2')
cmd.label("label_threshold_4.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [4.2]
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


cluster_dict = {"10.9099998474":[], "10.9099998474_arrows":[]}

cluster_dict["10.9099998474"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-19.5), float(48.0), float(70.0), float(1.0)]

cluster_dict["10.9099998474_arrows"] += cgo_arrow([-19.5,48.0,70.0], [-20.078,45.596,70.111], color="blue red", name="Arrows_10.9099998474_1")

cluster_dict["10.9099998474"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-18.5), float(54.0), float(67.5), float(1.0)]

cluster_dict["10.9099998474_arrows"] += cgo_arrow([-18.5,54.0,67.5], [-15.43,53.913,68.224], color="blue red", name="Arrows_10.9099998474_2")

cluster_dict["10.9099998474"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-20.0157079944), float(51.3604862025), float(68.8723051472), float(1.0)]


cluster_dict["10.9099998474"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-22.5), float(48.0), float(67.5), float(1.0)]

cluster_dict["10.9099998474_arrows"] += cgo_arrow([-22.5,48.0,67.5], [-23.246,45.59,68.299], color="red blue", name="Arrows_10.9099998474_3")

cmd.load_cgo(cluster_dict["10.9099998474"], "Features_10.9099998474", 1)
cmd.load_cgo(cluster_dict["10.9099998474_arrows"], "Arrows_10.9099998474")
cmd.set("transparency", 0.2,"Features_10.9099998474")
cmd.group("Pharmacophore_10.9099998474", members="Features_10.9099998474")
cmd.group("Pharmacophore_10.9099998474", members="Arrows_10.9099998474")

if dirpath:
    f = join(dirpath, "2/label_threshold_10.9099998474.mol2")
else:
    f = "2/label_threshold_10.9099998474.mol2"

cmd.load(f, 'label_threshold_10.9099998474')
cmd.hide('everything', 'label_threshold_10.9099998474')
cmd.label("label_threshold_10.9099998474", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.9099998474', members= 'label_threshold_10.9099998474')


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

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
