
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
    f = join(dirpath, "0/label_threshold_13.7.mol2")
else:
    f = "0/label_threshold_13.7.mol2"

cmd.load(f, 'label_threshold_13.7')
cmd.hide('everything', 'label_threshold_13.7')
cmd.label("label_threshold_13.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.7]
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


cluster_dict = {"20.9880008698":[], "20.9880008698_arrows":[]}

cluster_dict["20.9880008698"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-8.0), float(3.0), float(28.0), float(1.0)]

cluster_dict["20.9880008698_arrows"] += cgo_arrow([-8.0,3.0,28.0], [-8.314,2.092,26.001], color="blue red", name="Arrows_20.9880008698_1")

cluster_dict["20.9880008698"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-5.0), float(-1.0), float(27.0), float(1.0)]

cluster_dict["20.9880008698_arrows"] += cgo_arrow([-5.0,-1.0,27.0], [-6.273,-3.125,28.294], color="blue red", name="Arrows_20.9880008698_2")

cluster_dict["20.9880008698"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-3.0), float(1.5), float(27.0), float(1.0)]

cluster_dict["20.9880008698_arrows"] += cgo_arrow([-3.0,1.5,27.0], [-1.429,3.89,26.009], color="blue red", name="Arrows_20.9880008698_3")

cluster_dict["20.9880008698"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-3.0), float(1.5), float(27.0), float(1.0)]

cluster_dict["20.9880008698_arrows"] += cgo_arrow([-3.0,1.5,27.0], [-1.429,3.89,26.009], color="blue red", name="Arrows_20.9880008698_4")

cluster_dict["20.9880008698"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-8.5), float(-2.44829862188), float(20.3118334141), float(1.0)]


cluster_dict["20.9880008698"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-3.78102821215), float(0.781370134612), float(27.0099923098), float(1.0)]


cluster_dict["20.9880008698"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-5.0), float(1.5), float(25.5), float(1.0)]

cluster_dict["20.9880008698_arrows"] += cgo_arrow([-5.0,1.5,25.5], [-8.314,2.092,26.001], color="red blue", name="Arrows_20.9880008698_5")

cluster_dict["20.9880008698"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-5.0), float(2.0), float(23.0), float(1.0)]

cluster_dict["20.9880008698_arrows"] += cgo_arrow([-5.0,2.0,23.0], [-6.691,4.697,22.392], color="red blue", name="Arrows_20.9880008698_6")

cluster_dict["20.9880008698"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-4.0), float(4.0), float(32.5), float(1.0)]

cluster_dict["20.9880008698_arrows"] += cgo_arrow([-4.0,4.0,32.5], [-5.472,4.138,35.024], color="red blue", name="Arrows_20.9880008698_7")

cluster_dict["20.9880008698"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-3.5), float(4.5), float(29.5), float(1.0)]

cluster_dict["20.9880008698_arrows"] += cgo_arrow([-3.5,4.5,29.5], [-1.861,6.026,25.478], color="red blue", name="Arrows_20.9880008698_8")

cluster_dict["20.9880008698"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-1.0), float(0.5), float(31.0), float(1.0)]

cluster_dict["20.9880008698_arrows"] += cgo_arrow([-1.0,0.5,31.0], [0.797,-0.594,32.559], color="red blue", name="Arrows_20.9880008698_9")

cmd.load_cgo(cluster_dict["20.9880008698"], "Features_20.9880008698", 1)
cmd.load_cgo(cluster_dict["20.9880008698_arrows"], "Arrows_20.9880008698")
cmd.set("transparency", 0.2,"Features_20.9880008698")
cmd.group("Pharmacophore_20.9880008698", members="Features_20.9880008698")
cmd.group("Pharmacophore_20.9880008698", members="Arrows_20.9880008698")

if dirpath:
    f = join(dirpath, "0/label_threshold_20.9880008698.mol2")
else:
    f = "0/label_threshold_20.9880008698.mol2"

cmd.load(f, 'label_threshold_20.9880008698')
cmd.hide('everything', 'label_threshold_20.9880008698')
cmd.label("label_threshold_20.9880008698", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_20.9880008698', members= 'label_threshold_20.9880008698')


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
    f = join(dirpath, "1/label_threshold_13.4.mol2")
else:
    f = "1/label_threshold_13.4.mol2"

cmd.load(f, 'label_threshold_13.4')
cmd.hide('everything', 'label_threshold_13.4')
cmd.label("label_threshold_13.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.4]
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


cluster_dict = {"20.9880008698":[], "20.9880008698_arrows":[]}

cluster_dict["20.9880008698"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-8.0), float(3.0), float(28.0), float(1.0)]

cluster_dict["20.9880008698_arrows"] += cgo_arrow([-8.0,3.0,28.0], [-8.314,2.092,26.001], color="blue red", name="Arrows_20.9880008698_1")

cluster_dict["20.9880008698"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-5.0), float(-1.0), float(27.0), float(1.0)]

cluster_dict["20.9880008698_arrows"] += cgo_arrow([-5.0,-1.0,27.0], [-6.273,-3.125,28.294], color="blue red", name="Arrows_20.9880008698_2")

cluster_dict["20.9880008698"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-3.0), float(1.5), float(27.0), float(1.0)]

cluster_dict["20.9880008698_arrows"] += cgo_arrow([-3.0,1.5,27.0], [-1.429,3.89,26.009], color="blue red", name="Arrows_20.9880008698_3")

cluster_dict["20.9880008698"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-3.0), float(1.5), float(27.0), float(1.0)]

cluster_dict["20.9880008698_arrows"] += cgo_arrow([-3.0,1.5,27.0], [-1.429,3.89,26.009], color="blue red", name="Arrows_20.9880008698_4")

cluster_dict["20.9880008698"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-8.5), float(-2.44829862188), float(20.3118334141), float(1.0)]


cluster_dict["20.9880008698"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-3.78106683414), float(0.795302398501), float(27.005725637), float(1.0)]


cluster_dict["20.9880008698"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-5.0), float(1.5), float(25.5), float(1.0)]

cluster_dict["20.9880008698_arrows"] += cgo_arrow([-5.0,1.5,25.5], [-8.314,2.092,26.001], color="red blue", name="Arrows_20.9880008698_5")

cluster_dict["20.9880008698"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-5.0), float(2.0), float(23.0), float(1.0)]

cluster_dict["20.9880008698_arrows"] += cgo_arrow([-5.0,2.0,23.0], [-6.691,4.697,22.392], color="red blue", name="Arrows_20.9880008698_6")

cluster_dict["20.9880008698"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-4.0), float(4.0), float(32.5), float(1.0)]

cluster_dict["20.9880008698_arrows"] += cgo_arrow([-4.0,4.0,32.5], [-5.472,4.138,35.024], color="red blue", name="Arrows_20.9880008698_7")

cluster_dict["20.9880008698"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-3.5), float(4.5), float(29.5), float(1.0)]

cluster_dict["20.9880008698_arrows"] += cgo_arrow([-3.5,4.5,29.5], [-1.861,6.026,25.478], color="red blue", name="Arrows_20.9880008698_8")

cluster_dict["20.9880008698"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-1.0), float(0.5), float(31.0), float(1.0)]

cluster_dict["20.9880008698_arrows"] += cgo_arrow([-1.0,0.5,31.0], [0.797,-0.594,32.559], color="red blue", name="Arrows_20.9880008698_9")

cmd.load_cgo(cluster_dict["20.9880008698"], "Features_20.9880008698", 1)
cmd.load_cgo(cluster_dict["20.9880008698_arrows"], "Arrows_20.9880008698")
cmd.set("transparency", 0.2,"Features_20.9880008698")
cmd.group("Pharmacophore_20.9880008698", members="Features_20.9880008698")
cmd.group("Pharmacophore_20.9880008698", members="Arrows_20.9880008698")

if dirpath:
    f = join(dirpath, "1/label_threshold_20.9880008698.mol2")
else:
    f = "1/label_threshold_20.9880008698.mol2"

cmd.load(f, 'label_threshold_20.9880008698')
cmd.hide('everything', 'label_threshold_20.9880008698')
cmd.label("label_threshold_20.9880008698", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_20.9880008698', members= 'label_threshold_20.9880008698')


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
    f = join(dirpath, "2/label_threshold_14.2.mol2")
else:
    f = "2/label_threshold_14.2.mol2"

cmd.load(f, 'label_threshold_14.2')
cmd.hide('everything', 'label_threshold_14.2')
cmd.label("label_threshold_14.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [14.2]
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


cluster_dict = {"17.033000946":[], "17.033000946_arrows":[]}

cluster_dict["17.033000946"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-14.5), float(-7.5), float(21.5), float(1.0)]

cluster_dict["17.033000946_arrows"] += cgo_arrow([-14.5,-7.5,21.5], [-14.228,-10.728,20.922], color="blue red", name="Arrows_17.033000946_1")

cluster_dict["17.033000946"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-14.0), float(-0.5), float(21.5), float(1.0)]

cluster_dict["17.033000946_arrows"] += cgo_arrow([-14.0,-0.5,21.5], [-12.971,1.338,24.129], color="blue red", name="Arrows_17.033000946_2")

cluster_dict["17.033000946"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-5.0), float(-1.0), float(27.0), float(1.0)]

cluster_dict["17.033000946_arrows"] += cgo_arrow([-5.0,-1.0,27.0], [-6.273,-3.125,28.294], color="blue red", name="Arrows_17.033000946_3")

cluster_dict["17.033000946"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-6.0), float(1.0), float(24.0), float(1.0)]

cluster_dict["17.033000946_arrows"] += cgo_arrow([-6.0,1.0,24.0], [-6.966,-1.252,23.733], color="blue red", name="Arrows_17.033000946_4")

cluster_dict["17.033000946"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-2.5), float(1.5), float(27.0), float(1.0)]

cluster_dict["17.033000946_arrows"] += cgo_arrow([-2.5,1.5,27.0], [-1.429,3.89,26.009], color="blue red", name="Arrows_17.033000946_5")

cluster_dict["17.033000946"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-13.3640230673), float(-3.93234932172), float(22.2990220113), float(1.0)]


cluster_dict["17.033000946"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-8.5), float(-2.69201015266), float(20.4015374726), float(1.0)]


cluster_dict["17.033000946"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-3.94383867025), float(0.19413497002), float(27.1381220728), float(1.0)]


cluster_dict["17.033000946"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-15.0), float(-1.0), float(21.5), float(1.0)]

cluster_dict["17.033000946_arrows"] += cgo_arrow([-15.0,-1.0,21.5], [-18.128,-2.744,23.189], color="red blue", name="Arrows_17.033000946_6")

cluster_dict["17.033000946"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-5.0), float(1.5), float(23.0), float(1.0)]

cluster_dict["17.033000946_arrows"] += cgo_arrow([-5.0,1.5,23.0], [-6.691,4.697,22.392], color="red blue", name="Arrows_17.033000946_7")

cluster_dict["17.033000946"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-5.0), float(1.5), float(25.5), float(1.0)]

cluster_dict["17.033000946_arrows"] += cgo_arrow([-5.0,1.5,25.5], [-8.314,2.092,26.001], color="red blue", name="Arrows_17.033000946_8")

cluster_dict["17.033000946"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-1.0), float(0.5), float(31.0), float(1.0)]

cluster_dict["17.033000946_arrows"] += cgo_arrow([-1.0,0.5,31.0], [0.797,-0.594,32.559], color="red blue", name="Arrows_17.033000946_9")

cmd.load_cgo(cluster_dict["17.033000946"], "Features_17.033000946", 1)
cmd.load_cgo(cluster_dict["17.033000946_arrows"], "Arrows_17.033000946")
cmd.set("transparency", 0.2,"Features_17.033000946")
cmd.group("Pharmacophore_17.033000946", members="Features_17.033000946")
cmd.group("Pharmacophore_17.033000946", members="Arrows_17.033000946")

if dirpath:
    f = join(dirpath, "2/label_threshold_17.033000946.mol2")
else:
    f = "2/label_threshold_17.033000946.mol2"

cmd.load(f, 'label_threshold_17.033000946')
cmd.hide('everything', 'label_threshold_17.033000946')
cmd.label("label_threshold_17.033000946", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.033000946', members= 'label_threshold_17.033000946')


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
    f = join(dirpath, "3/label_threshold_10.0.mol2")
else:
    f = "3/label_threshold_10.0.mol2"

cmd.load(f, 'label_threshold_10.0')
cmd.hide('everything', 'label_threshold_10.0')
cmd.label("label_threshold_10.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.0]
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


cluster_dict = {"14.7270002365":[], "14.7270002365_arrows":[]}

cluster_dict["14.7270002365"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-12.5), float(4.5), float(8.5), float(1.0)]

cluster_dict["14.7270002365_arrows"] += cgo_arrow([-12.5,4.5,8.5], [-10.305,3.336,10.139], color="blue red", name="Arrows_14.7270002365_1")

cluster_dict["14.7270002365"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-12.5), float(5.0), float(6.5), float(1.0)]

cluster_dict["14.7270002365_arrows"] += cgo_arrow([-12.5,5.0,6.5], [-14.913,4.845,5.836], color="blue red", name="Arrows_14.7270002365_2")

cluster_dict["14.7270002365"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-9.0), float(9.0), float(10.0), float(1.0)]

cluster_dict["14.7270002365_arrows"] += cgo_arrow([-9.0,9.0,10.0], [-7.009,6.434,10.305], color="blue red", name="Arrows_14.7270002365_3")

cluster_dict["14.7270002365"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-9.0), float(9.0), float(10.0), float(1.0)]

cluster_dict["14.7270002365_arrows"] += cgo_arrow([-9.0,9.0,10.0], [-7.009,6.434,10.305], color="blue red", name="Arrows_14.7270002365_4")

cluster_dict["14.7270002365"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-6.0), float(10.0), float(2.5), float(1.0)]

cluster_dict["14.7270002365_arrows"] += cgo_arrow([-6.0,10.0,2.5], [-7.103,12.548,3.941], color="blue red", name="Arrows_14.7270002365_5")

cluster_dict["14.7270002365"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-8.96656733607), float(8.81008756201), float(6.8382272265), float(1.0)]


cluster_dict["14.7270002365"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-14.0), float(8.5), float(10.5), float(1.0)]

cluster_dict["14.7270002365_arrows"] += cgo_arrow([-14.0,8.5,10.5], [-16.699,9.26,11.346], color="red blue", name="Arrows_14.7270002365_6")

cluster_dict["14.7270002365"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-11.0), float(7.5), float(6.5), float(1.0)]

cluster_dict["14.7270002365_arrows"] += cgo_arrow([-11.0,7.5,6.5], [-13.664,9.104,5.429], color="red blue", name="Arrows_14.7270002365_7")

cmd.load_cgo(cluster_dict["14.7270002365"], "Features_14.7270002365", 1)
cmd.load_cgo(cluster_dict["14.7270002365_arrows"], "Arrows_14.7270002365")
cmd.set("transparency", 0.2,"Features_14.7270002365")
cmd.group("Pharmacophore_14.7270002365", members="Features_14.7270002365")
cmd.group("Pharmacophore_14.7270002365", members="Arrows_14.7270002365")

if dirpath:
    f = join(dirpath, "3/label_threshold_14.7270002365.mol2")
else:
    f = "3/label_threshold_14.7270002365.mol2"

cmd.load(f, 'label_threshold_14.7270002365')
cmd.hide('everything', 'label_threshold_14.7270002365')
cmd.label("label_threshold_14.7270002365", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.7270002365', members= 'label_threshold_14.7270002365')


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
    f = join(dirpath, "4/label_threshold_10.0.mol2")
else:
    f = "4/label_threshold_10.0.mol2"

cmd.load(f, 'label_threshold_10.0')
cmd.hide('everything', 'label_threshold_10.0')
cmd.label("label_threshold_10.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.0]
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


cluster_dict = {"14.0970001221":[], "14.0970001221_arrows":[]}

cluster_dict["14.0970001221"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-14.5), float(-7.5), float(21.5), float(1.0)]

cluster_dict["14.0970001221_arrows"] += cgo_arrow([-14.5,-7.5,21.5], [-14.228,-10.728,20.922], color="blue red", name="Arrows_14.0970001221_1")

cluster_dict["14.0970001221"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-14.0), float(-0.5), float(21.5), float(1.0)]

cluster_dict["14.0970001221_arrows"] += cgo_arrow([-14.0,-0.5,21.5], [-12.971,1.338,24.129], color="blue red", name="Arrows_14.0970001221_2")

cluster_dict["14.0970001221"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-13.0), float(1.0), float(21.0), float(1.0)]

cluster_dict["14.0970001221_arrows"] += cgo_arrow([-13.0,1.0,21.0], [-12.971,1.338,24.129], color="blue red", name="Arrows_14.0970001221_3")

cluster_dict["14.0970001221"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-11.0), float(-5.0), float(23.5), float(1.0)]

cluster_dict["14.0970001221_arrows"] += cgo_arrow([-11.0,-5.0,23.5], [-10.199,-6.72,22.235], color="blue red", name="Arrows_14.0970001221_4")

cluster_dict["14.0970001221"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-6.0), float(1.5), float(24.0), float(1.0)]

cluster_dict["14.0970001221_arrows"] += cgo_arrow([-6.0,1.5,24.0], [-6.966,-1.252,23.733], color="blue red", name="Arrows_14.0970001221_5")

cluster_dict["14.0970001221"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-13.4334821607), float(-4.35849451844), float(23.4710020673), float(1.0)]


cluster_dict["14.0970001221"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-8.5), float(-2.63222170249), float(20.4306236026), float(1.0)]


cluster_dict["14.0970001221"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-6.6672141948), float(0.450498108582), float(26.7523386979), float(1.0)]


cluster_dict["14.0970001221"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-14.0), float(-2.0), float(24.5), float(1.0)]

cluster_dict["14.0970001221_arrows"] += cgo_arrow([-14.0,-2.0,24.5], [-15.564,-2.598,27.017], color="red blue", name="Arrows_14.0970001221_6")

cluster_dict["14.0970001221"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-15.0), float(-1.0), float(21.5), float(1.0)]

cluster_dict["14.0970001221_arrows"] += cgo_arrow([-15.0,-1.0,21.5], [-18.128,-2.744,23.189], color="red blue", name="Arrows_14.0970001221_7")

cmd.load_cgo(cluster_dict["14.0970001221"], "Features_14.0970001221", 1)
cmd.load_cgo(cluster_dict["14.0970001221_arrows"], "Arrows_14.0970001221")
cmd.set("transparency", 0.2,"Features_14.0970001221")
cmd.group("Pharmacophore_14.0970001221", members="Features_14.0970001221")
cmd.group("Pharmacophore_14.0970001221", members="Arrows_14.0970001221")

if dirpath:
    f = join(dirpath, "4/label_threshold_14.0970001221.mol2")
else:
    f = "4/label_threshold_14.0970001221.mol2"

cmd.load(f, 'label_threshold_14.0970001221')
cmd.hide('everything', 'label_threshold_14.0970001221')
cmd.label("label_threshold_14.0970001221", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.0970001221', members= 'label_threshold_14.0970001221')


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
    f = join(dirpath, "5/label_threshold_8.1.mol2")
else:
    f = "5/label_threshold_8.1.mol2"

cmd.load(f, 'label_threshold_8.1')
cmd.hide('everything', 'label_threshold_8.1')
cmd.label("label_threshold_8.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [8.1]
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


cluster_dict = {"12.0484995842":[], "12.0484995842_arrows":[]}

cluster_dict["12.0484995842"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-11.0), float(-5.0), float(23.5), float(1.0)]

cluster_dict["12.0484995842_arrows"] += cgo_arrow([-11.0,-5.0,23.5], [-10.199,-6.72,22.235], color="blue red", name="Arrows_12.0484995842_1")

cluster_dict["12.0484995842"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-4.5), float(-1.0), float(27.0), float(1.0)]

cluster_dict["12.0484995842_arrows"] += cgo_arrow([-4.5,-1.0,27.0], [-4.584,-1.966,29.695], color="blue red", name="Arrows_12.0484995842_2")

cluster_dict["12.0484995842"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-4.0), float(-1.0), float(22.5), float(1.0)]

cluster_dict["12.0484995842_arrows"] += cgo_arrow([-4.0,-1.0,22.5], [-1.617,-2.921,22.521], color="blue red", name="Arrows_12.0484995842_3")

cluster_dict["12.0484995842"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-0.5), float(-1.0), float(24.0), float(1.0)]

cluster_dict["12.0484995842_arrows"] += cgo_arrow([-0.5,-1.0,24.0], [-1.617,-2.921,22.521], color="blue red", name="Arrows_12.0484995842_4")

cluster_dict["12.0484995842"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-12.8808239872), float(-5.64732724572), float(26.988004426), float(1.0)]


cluster_dict["12.0484995842"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-5.02291724418), float(-8.02114801849), float(33.4552867368), float(1.0)]


cluster_dict["12.0484995842"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-3.7200038411), float(-1.26762223337), float(26.5612831936), float(1.0)]


cluster_dict["12.0484995842"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-14.0), float(-2.5), float(24.5), float(1.0)]

cluster_dict["12.0484995842_arrows"] += cgo_arrow([-14.0,-2.5,24.5], [-15.564,-2.598,27.017], color="red blue", name="Arrows_12.0484995842_5")

cluster_dict["12.0484995842"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-3.5), float(-2.0), float(24.0), float(1.0)]

cluster_dict["12.0484995842_arrows"] += cgo_arrow([-3.5,-2.0,24.0], [-3.866,-4.266,23.606], color="red blue", name="Arrows_12.0484995842_6")

cmd.load_cgo(cluster_dict["12.0484995842"], "Features_12.0484995842", 1)
cmd.load_cgo(cluster_dict["12.0484995842_arrows"], "Arrows_12.0484995842")
cmd.set("transparency", 0.2,"Features_12.0484995842")
cmd.group("Pharmacophore_12.0484995842", members="Features_12.0484995842")
cmd.group("Pharmacophore_12.0484995842", members="Arrows_12.0484995842")

if dirpath:
    f = join(dirpath, "5/label_threshold_12.0484995842.mol2")
else:
    f = "5/label_threshold_12.0484995842.mol2"

cmd.load(f, 'label_threshold_12.0484995842')
cmd.hide('everything', 'label_threshold_12.0484995842')
cmd.label("label_threshold_12.0484995842", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.0484995842', members= 'label_threshold_12.0484995842')


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
    f = join(dirpath, "6/label_threshold_2.4.mol2")
else:
    f = "6/label_threshold_2.4.mol2"

cmd.load(f, 'label_threshold_2.4')
cmd.hide('everything', 'label_threshold_2.4')
cmd.label("label_threshold_2.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [2.4]
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


cluster_dict = {"11.7089996338":[], "11.7089996338_arrows":[]}

cluster_dict["11.7089996338"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-10.5), float(-7.0), float(9.0), float(1.0)]

cluster_dict["11.7089996338_arrows"] += cgo_arrow([-10.5,-7.0,9.0], [-9.655,-10.024,9.811], color="blue red", name="Arrows_11.7089996338_1")

cluster_dict["11.7089996338"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-9.0), float(-3.0), float(7.0), float(1.0)]

cluster_dict["11.7089996338_arrows"] += cgo_arrow([-9.0,-3.0,7.0], [-8.077,-1.117,8.926], color="blue red", name="Arrows_11.7089996338_2")

cluster_dict["11.7089996338"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-9.17209664208), float(-5.33342170345), float(5.81352939116), float(1.0)]


cluster_dict["11.7089996338"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-9.5), float(-8.0), float(5.5), float(1.0)]

cluster_dict["11.7089996338_arrows"] += cgo_arrow([-9.5,-8.0,5.5], [-10.385,-9.862,3.535], color="red blue", name="Arrows_11.7089996338_3")

cmd.load_cgo(cluster_dict["11.7089996338"], "Features_11.7089996338", 1)
cmd.load_cgo(cluster_dict["11.7089996338_arrows"], "Arrows_11.7089996338")
cmd.set("transparency", 0.2,"Features_11.7089996338")
cmd.group("Pharmacophore_11.7089996338", members="Features_11.7089996338")
cmd.group("Pharmacophore_11.7089996338", members="Arrows_11.7089996338")

if dirpath:
    f = join(dirpath, "6/label_threshold_11.7089996338.mol2")
else:
    f = "6/label_threshold_11.7089996338.mol2"

cmd.load(f, 'label_threshold_11.7089996338')
cmd.hide('everything', 'label_threshold_11.7089996338')
cmd.label("label_threshold_11.7089996338", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.7089996338', members= 'label_threshold_11.7089996338')


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
