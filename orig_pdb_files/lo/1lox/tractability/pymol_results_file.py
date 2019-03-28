
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
    f = join(dirpath, "0/label_threshold_17.6.mol2")
else:
    f = "0/label_threshold_17.6.mol2"

cmd.load(f, 'label_threshold_17.6')
cmd.hide('everything', 'label_threshold_17.6')
cmd.label("label_threshold_17.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [17.6]
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


cluster_dict = {"22.420999527":[], "22.420999527_arrows":[]}

cluster_dict["22.420999527"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-26.8166365629), float(150.561683856), float(54.2155399368), float(1.0)]


cluster_dict["22.420999527"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-28.0), float(154.0), float(57.0), float(1.0)]

cluster_dict["22.420999527_arrows"] += cgo_arrow([-28.0,154.0,57.0], [-28.645,151.11,60.419], color="red blue", name="Arrows_22.420999527_1")

cluster_dict["22.420999527"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-25.0), float(154.5), float(58.0), float(1.0)]


cmd.load_cgo(cluster_dict["22.420999527"], "Features_22.420999527", 1)
cmd.load_cgo(cluster_dict["22.420999527_arrows"], "Arrows_22.420999527")
cmd.set("transparency", 0.2,"Features_22.420999527")
cmd.group("Pharmacophore_22.420999527", members="Features_22.420999527")
cmd.group("Pharmacophore_22.420999527", members="Arrows_22.420999527")

if dirpath:
    f = join(dirpath, "0/label_threshold_22.420999527.mol2")
else:
    f = "0/label_threshold_22.420999527.mol2"

cmd.load(f, 'label_threshold_22.420999527')
cmd.hide('everything', 'label_threshold_22.420999527')
cmd.label("label_threshold_22.420999527", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_22.420999527', members= 'label_threshold_22.420999527')


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
    f = join(dirpath, "1/label_threshold_17.6.mol2")
else:
    f = "1/label_threshold_17.6.mol2"

cmd.load(f, 'label_threshold_17.6')
cmd.hide('everything', 'label_threshold_17.6')
cmd.label("label_threshold_17.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [17.6]
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


cluster_dict = {"22.420999527":[], "22.420999527_arrows":[]}

cluster_dict["22.420999527"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-26.8166365629), float(150.561683856), float(54.2155399368), float(1.0)]


cluster_dict["22.420999527"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-28.0), float(154.0), float(57.0), float(1.0)]

cluster_dict["22.420999527_arrows"] += cgo_arrow([-28.0,154.0,57.0], [-28.645,151.11,60.419], color="red blue", name="Arrows_22.420999527_1")

cluster_dict["22.420999527"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-25.0), float(154.5), float(58.0), float(1.0)]


cmd.load_cgo(cluster_dict["22.420999527"], "Features_22.420999527", 1)
cmd.load_cgo(cluster_dict["22.420999527_arrows"], "Arrows_22.420999527")
cmd.set("transparency", 0.2,"Features_22.420999527")
cmd.group("Pharmacophore_22.420999527", members="Features_22.420999527")
cmd.group("Pharmacophore_22.420999527", members="Arrows_22.420999527")

if dirpath:
    f = join(dirpath, "1/label_threshold_22.420999527.mol2")
else:
    f = "1/label_threshold_22.420999527.mol2"

cmd.load(f, 'label_threshold_22.420999527')
cmd.hide('everything', 'label_threshold_22.420999527')
cmd.label("label_threshold_22.420999527", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_22.420999527', members= 'label_threshold_22.420999527')


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
    f = join(dirpath, "2/label_threshold_17.4.mol2")
else:
    f = "2/label_threshold_17.4.mol2"

cmd.load(f, 'label_threshold_17.4')
cmd.hide('everything', 'label_threshold_17.4')
cmd.label("label_threshold_17.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [17.4]
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


cluster_dict = {"22.4039993286":[], "22.4039993286_arrows":[]}

cluster_dict["22.4039993286"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-26.820143328), float(150.568563396), float(54.2452055912), float(1.0)]


cluster_dict["22.4039993286"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-28.0), float(154.0), float(57.0), float(1.0)]

cluster_dict["22.4039993286_arrows"] += cgo_arrow([-28.0,154.0,57.0], [-28.645,151.11,60.419], color="red blue", name="Arrows_22.4039993286_1")

cluster_dict["22.4039993286"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-25.0), float(154.5), float(58.0), float(1.0)]


cmd.load_cgo(cluster_dict["22.4039993286"], "Features_22.4039993286", 1)
cmd.load_cgo(cluster_dict["22.4039993286_arrows"], "Arrows_22.4039993286")
cmd.set("transparency", 0.2,"Features_22.4039993286")
cmd.group("Pharmacophore_22.4039993286", members="Features_22.4039993286")
cmd.group("Pharmacophore_22.4039993286", members="Arrows_22.4039993286")

if dirpath:
    f = join(dirpath, "2/label_threshold_22.4039993286.mol2")
else:
    f = "2/label_threshold_22.4039993286.mol2"

cmd.load(f, 'label_threshold_22.4039993286')
cmd.hide('everything', 'label_threshold_22.4039993286')
cmd.label("label_threshold_22.4039993286", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_22.4039993286', members= 'label_threshold_22.4039993286')


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
    f = join(dirpath, "3/label_threshold_16.9.mol2")
else:
    f = "3/label_threshold_16.9.mol2"

cmd.load(f, 'label_threshold_16.9')
cmd.hide('everything', 'label_threshold_16.9')
cmd.label("label_threshold_16.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [16.9]
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


cluster_dict = {"22.2099990845":[], "22.2099990845_arrows":[]}

cluster_dict["22.2099990845"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-26.0), float(152.0), float(59.5), float(1.0)]

cluster_dict["22.2099990845_arrows"] += cgo_arrow([-26.0,152.0,59.5], [-25.679,152.997,62.198], color="blue red", name="Arrows_22.2099990845_1")

cluster_dict["22.2099990845"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-26.8594349657), float(150.613972569), float(54.4799358648), float(1.0)]


cluster_dict["22.2099990845"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-28.0), float(154.0), float(57.0), float(1.0)]

cluster_dict["22.2099990845_arrows"] += cgo_arrow([-28.0,154.0,57.0], [-28.645,151.11,60.419], color="red blue", name="Arrows_22.2099990845_2")

cluster_dict["22.2099990845"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-25.0), float(154.5), float(58.0), float(1.0)]


cmd.load_cgo(cluster_dict["22.2099990845"], "Features_22.2099990845", 1)
cmd.load_cgo(cluster_dict["22.2099990845_arrows"], "Arrows_22.2099990845")
cmd.set("transparency", 0.2,"Features_22.2099990845")
cmd.group("Pharmacophore_22.2099990845", members="Features_22.2099990845")
cmd.group("Pharmacophore_22.2099990845", members="Arrows_22.2099990845")

if dirpath:
    f = join(dirpath, "3/label_threshold_22.2099990845.mol2")
else:
    f = "3/label_threshold_22.2099990845.mol2"

cmd.load(f, 'label_threshold_22.2099990845')
cmd.hide('everything', 'label_threshold_22.2099990845')
cmd.label("label_threshold_22.2099990845", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_22.2099990845', members= 'label_threshold_22.2099990845')


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
    f = join(dirpath, "4/label_threshold_10.7.mol2")
else:
    f = "4/label_threshold_10.7.mol2"

cmd.load(f, 'label_threshold_10.7')
cmd.hide('everything', 'label_threshold_10.7')
cmd.label("label_threshold_10.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.7]
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


cluster_dict = {"20.6819992065":[], "20.6819992065_arrows":[]}

cluster_dict["20.6819992065"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-41.1099824188), float(147.719453279), float(47.9683607039), float(1.0)]


cluster_dict["20.6819992065"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-43.5), float(150.0), float(50.5), float(1.0)]

cluster_dict["20.6819992065_arrows"] += cgo_arrow([-43.5,150.0,50.5], [-44.889,152.307,50.103], color="red blue", name="Arrows_20.6819992065_1")

cmd.load_cgo(cluster_dict["20.6819992065"], "Features_20.6819992065", 1)
cmd.load_cgo(cluster_dict["20.6819992065_arrows"], "Arrows_20.6819992065")
cmd.set("transparency", 0.2,"Features_20.6819992065")
cmd.group("Pharmacophore_20.6819992065", members="Features_20.6819992065")
cmd.group("Pharmacophore_20.6819992065", members="Arrows_20.6819992065")

if dirpath:
    f = join(dirpath, "4/label_threshold_20.6819992065.mol2")
else:
    f = "4/label_threshold_20.6819992065.mol2"

cmd.load(f, 'label_threshold_20.6819992065')
cmd.hide('everything', 'label_threshold_20.6819992065')
cmd.label("label_threshold_20.6819992065", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_20.6819992065', members= 'label_threshold_20.6819992065')


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
    f = join(dirpath, "5/label_threshold_15.7.mol2")
else:
    f = "5/label_threshold_15.7.mol2"

cmd.load(f, 'label_threshold_15.7')
cmd.hide('everything', 'label_threshold_15.7')
cmd.label("label_threshold_15.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.7]
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


cluster_dict = {"17.0769996643":[], "17.0769996643_arrows":[]}

cluster_dict["17.0769996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-12.5), float(142.5), float(46.0), float(1.0)]

cluster_dict["17.0769996643_arrows"] += cgo_arrow([-12.5,142.5,46.0], [-9.665,144.622,45.987], color="blue red", name="Arrows_17.0769996643_1")

cluster_dict["17.0769996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-10.0), float(138.5), float(49.5), float(1.0)]

cluster_dict["17.0769996643_arrows"] += cgo_arrow([-10.0,138.5,49.5], [-7.298,137.092,49.558], color="blue red", name="Arrows_17.0769996643_2")

cluster_dict["17.0769996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-7.0), float(139.5), float(48.0), float(1.0)]

cluster_dict["17.0769996643_arrows"] += cgo_arrow([-7.0,139.5,48.0], [-7.298,137.092,49.558], color="blue red", name="Arrows_17.0769996643_3")

cluster_dict["17.0769996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(1.0), float(138.5), float(46.5), float(1.0)]

cluster_dict["17.0769996643_arrows"] += cgo_arrow([1.0,138.5,46.5], [-0.043,135.193,45.922], color="blue red", name="Arrows_17.0769996643_4")

cluster_dict["17.0769996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-3.5), float(146.5), float(46.0), float(1.0)]

cluster_dict["17.0769996643_arrows"] += cgo_arrow([-3.5,146.5,46.0], [-4.668,144.707,44.258], color="blue red", name="Arrows_17.0769996643_5")

cluster_dict["17.0769996643"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-4.97006010744), float(142.422856048), float(47.1333514378), float(1.0)]


cluster_dict["17.0769996643"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-11.3135629733), float(145.718996015), float(40.587923436), float(1.0)]


cluster_dict["17.0769996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-9.5), float(137.0), float(51.0), float(1.0)]

cluster_dict["17.0769996643_arrows"] += cgo_arrow([-9.5,137.0,51.0], [-7.298,137.092,49.558], color="red blue", name="Arrows_17.0769996643_6")

cluster_dict["17.0769996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-5.5), float(139.0), float(48.0), float(1.0)]

cluster_dict["17.0769996643_arrows"] += cgo_arrow([-5.5,139.0,48.0], [-7.298,137.092,49.558], color="red blue", name="Arrows_17.0769996643_7")

cluster_dict["17.0769996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-1.0), float(145.0), float(48.0), float(1.0)]

cluster_dict["17.0769996643_arrows"] += cgo_arrow([-1.0,145.0,48.0], [-0.44,144.992,50.646], color="red blue", name="Arrows_17.0769996643_8")

cluster_dict["17.0769996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-1.0), float(145.0), float(48.0), float(1.0)]

cluster_dict["17.0769996643_arrows"] += cgo_arrow([-1.0,145.0,48.0], [-0.44,144.992,50.646], color="red blue", name="Arrows_17.0769996643_9")

cluster_dict["17.0769996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(0.0), float(140.5), float(47.5), float(1.0)]

cluster_dict["17.0769996643_arrows"] += cgo_arrow([0.0,140.5,47.5], [0.695,140.446,50.066], color="red blue", name="Arrows_17.0769996643_10")

cmd.load_cgo(cluster_dict["17.0769996643"], "Features_17.0769996643", 1)
cmd.load_cgo(cluster_dict["17.0769996643_arrows"], "Arrows_17.0769996643")
cmd.set("transparency", 0.2,"Features_17.0769996643")
cmd.group("Pharmacophore_17.0769996643", members="Features_17.0769996643")
cmd.group("Pharmacophore_17.0769996643", members="Arrows_17.0769996643")

if dirpath:
    f = join(dirpath, "5/label_threshold_17.0769996643.mol2")
else:
    f = "5/label_threshold_17.0769996643.mol2"

cmd.load(f, 'label_threshold_17.0769996643')
cmd.hide('everything', 'label_threshold_17.0769996643')
cmd.label("label_threshold_17.0769996643", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.0769996643', members= 'label_threshold_17.0769996643')


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


cluster_dict = {"16.8369998932":[], "16.8369998932_arrows":[]}

cluster_dict["16.8369998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-12.0), float(142.5), float(46.0), float(1.0)]

cluster_dict["16.8369998932_arrows"] += cgo_arrow([-12.0,142.5,46.0], [-9.665,144.622,45.987], color="blue red", name="Arrows_16.8369998932_1")

cluster_dict["16.8369998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-7.0), float(139.5), float(48.0), float(1.0)]

cluster_dict["16.8369998932_arrows"] += cgo_arrow([-7.0,139.5,48.0], [-7.298,137.092,49.558], color="blue red", name="Arrows_16.8369998932_2")

cluster_dict["16.8369998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(1.0), float(138.5), float(46.5), float(1.0)]

cluster_dict["16.8369998932_arrows"] += cgo_arrow([1.0,138.5,46.5], [-0.043,135.193,45.922], color="blue red", name="Arrows_16.8369998932_3")

cluster_dict["16.8369998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-3.5), float(146.5), float(46.0), float(1.0)]

cluster_dict["16.8369998932_arrows"] += cgo_arrow([-3.5,146.5,46.0], [-4.668,144.707,44.258], color="blue red", name="Arrows_16.8369998932_4")

cluster_dict["16.8369998932"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-3.55981374539), float(143.007358506), float(47.1101520032), float(1.0)]


cluster_dict["16.8369998932"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-10.6928335597), float(146.422748694), float(39.9988980016), float(1.0)]


cluster_dict["16.8369998932"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-10.5), float(148.578731042), float(39.7572137962), float(1.0)]


cluster_dict["16.8369998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-5.5), float(139.0), float(48.0), float(1.0)]

cluster_dict["16.8369998932_arrows"] += cgo_arrow([-5.5,139.0,48.0], [-7.298,137.092,49.558], color="red blue", name="Arrows_16.8369998932_5")

cluster_dict["16.8369998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-1.0), float(145.0), float(48.0), float(1.0)]

cluster_dict["16.8369998932_arrows"] += cgo_arrow([-1.0,145.0,48.0], [-0.44,144.992,50.646], color="red blue", name="Arrows_16.8369998932_6")

cluster_dict["16.8369998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-1.0), float(145.0), float(48.0), float(1.0)]

cluster_dict["16.8369998932_arrows"] += cgo_arrow([-1.0,145.0,48.0], [-0.44,144.992,50.646], color="red blue", name="Arrows_16.8369998932_7")

cluster_dict["16.8369998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(0.5), float(140.5), float(47.5), float(1.0)]

cluster_dict["16.8369998932_arrows"] += cgo_arrow([0.5,140.5,47.5], [0.695,140.446,50.066], color="red blue", name="Arrows_16.8369998932_8")

cluster_dict["16.8369998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(9.0), float(141.5), float(41.5), float(1.0)]

cluster_dict["16.8369998932_arrows"] += cgo_arrow([9.0,141.5,41.5], [9.366,139.015,42.432], color="red blue", name="Arrows_16.8369998932_9")

cmd.load_cgo(cluster_dict["16.8369998932"], "Features_16.8369998932", 1)
cmd.load_cgo(cluster_dict["16.8369998932_arrows"], "Arrows_16.8369998932")
cmd.set("transparency", 0.2,"Features_16.8369998932")
cmd.group("Pharmacophore_16.8369998932", members="Features_16.8369998932")
cmd.group("Pharmacophore_16.8369998932", members="Arrows_16.8369998932")

if dirpath:
    f = join(dirpath, "6/label_threshold_16.8369998932.mol2")
else:
    f = "6/label_threshold_16.8369998932.mol2"

cmd.load(f, 'label_threshold_16.8369998932')
cmd.hide('everything', 'label_threshold_16.8369998932')
cmd.label("label_threshold_16.8369998932", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.8369998932', members= 'label_threshold_16.8369998932')


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
    f = join(dirpath, "7/label_threshold_15.5.mol2")
else:
    f = "7/label_threshold_15.5.mol2"

cmd.load(f, 'label_threshold_15.5')
cmd.hide('everything', 'label_threshold_15.5')
cmd.label("label_threshold_15.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.5]
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


cluster_dict = {"16.811000824":[], "16.811000824_arrows":[]}

cluster_dict["16.811000824"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-12.5), float(142.5), float(46.0), float(1.0)]

cluster_dict["16.811000824_arrows"] += cgo_arrow([-12.5,142.5,46.0], [-9.665,144.622,45.987], color="blue red", name="Arrows_16.811000824_1")

cluster_dict["16.811000824"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-10.0), float(138.5), float(49.5), float(1.0)]

cluster_dict["16.811000824_arrows"] += cgo_arrow([-10.0,138.5,49.5], [-7.298,137.092,49.558], color="blue red", name="Arrows_16.811000824_2")

cluster_dict["16.811000824"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-7.0), float(139.5), float(48.0), float(1.0)]

cluster_dict["16.811000824_arrows"] += cgo_arrow([-7.0,139.5,48.0], [-7.298,137.092,49.558], color="blue red", name="Arrows_16.811000824_3")

cluster_dict["16.811000824"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(1.0), float(138.5), float(46.5), float(1.0)]

cluster_dict["16.811000824_arrows"] += cgo_arrow([1.0,138.5,46.5], [-0.043,135.193,45.922], color="blue red", name="Arrows_16.811000824_4")

cluster_dict["16.811000824"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-3.5), float(146.5), float(46.0), float(1.0)]

cluster_dict["16.811000824_arrows"] += cgo_arrow([-3.5,146.5,46.0], [-4.668,144.707,44.258], color="blue red", name="Arrows_16.811000824_5")

cluster_dict["16.811000824"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-12.5610262097), float(146.81505447), float(39.8148740758), float(1.0)]


cluster_dict["16.811000824"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-6.40372189211), float(141.595149575), float(47.5527450446), float(1.0)]


cluster_dict["16.811000824"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-9.5), float(137.0), float(51.0), float(1.0)]

cluster_dict["16.811000824_arrows"] += cgo_arrow([-9.5,137.0,51.0], [-7.298,137.092,49.558], color="red blue", name="Arrows_16.811000824_6")

cluster_dict["16.811000824"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-5.5), float(139.0), float(48.0), float(1.0)]

cluster_dict["16.811000824_arrows"] += cgo_arrow([-5.5,139.0,48.0], [-7.298,137.092,49.558], color="red blue", name="Arrows_16.811000824_7")

cluster_dict["16.811000824"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-1.0), float(145.0), float(48.0), float(1.0)]

cluster_dict["16.811000824_arrows"] += cgo_arrow([-1.0,145.0,48.0], [-0.44,144.992,50.646], color="red blue", name="Arrows_16.811000824_8")

cluster_dict["16.811000824"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-1.0), float(145.0), float(48.0), float(1.0)]

cluster_dict["16.811000824_arrows"] += cgo_arrow([-1.0,145.0,48.0], [-0.44,144.992,50.646], color="red blue", name="Arrows_16.811000824_9")

cluster_dict["16.811000824"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(0.0), float(140.5), float(47.5), float(1.0)]

cluster_dict["16.811000824_arrows"] += cgo_arrow([0.0,140.5,47.5], [0.695,140.446,50.066], color="red blue", name="Arrows_16.811000824_10")

cmd.load_cgo(cluster_dict["16.811000824"], "Features_16.811000824", 1)
cmd.load_cgo(cluster_dict["16.811000824_arrows"], "Arrows_16.811000824")
cmd.set("transparency", 0.2,"Features_16.811000824")
cmd.group("Pharmacophore_16.811000824", members="Features_16.811000824")
cmd.group("Pharmacophore_16.811000824", members="Arrows_16.811000824")

if dirpath:
    f = join(dirpath, "7/label_threshold_16.811000824.mol2")
else:
    f = "7/label_threshold_16.811000824.mol2"

cmd.load(f, 'label_threshold_16.811000824')
cmd.hide('everything', 'label_threshold_16.811000824')
cmd.label("label_threshold_16.811000824", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.811000824', members= 'label_threshold_16.811000824')


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
    f = join(dirpath, "8/label_threshold_12.6.mol2")
else:
    f = "8/label_threshold_12.6.mol2"

cmd.load(f, 'label_threshold_12.6')
cmd.hide('everything', 'label_threshold_12.6')
cmd.label("label_threshold_12.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.6]
gfiles = ['8/apolar.grd', '8/acceptor.grd']
grids = ['apolar', 'acceptor']
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


cluster_dict = {"15.1850004196":[], "15.1850004196_arrows":[]}

cluster_dict["15.1850004196"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-2.5), float(145.5), float(43.0), float(1.0)]

cluster_dict["15.1850004196_arrows"] += cgo_arrow([-2.5,145.5,43.0], [-4.668,144.707,44.258], color="blue red", name="Arrows_15.1850004196_1")

cluster_dict["15.1850004196"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-2.5), float(143.5), float(46.5), float(1.0)]

cluster_dict["15.1850004196_arrows"] += cgo_arrow([-2.5,143.5,46.5], [-4.668,144.707,44.258], color="blue red", name="Arrows_15.1850004196_2")

cluster_dict["15.1850004196"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-3.5), float(146.5), float(46.0), float(1.0)]

cluster_dict["15.1850004196_arrows"] += cgo_arrow([-3.5,146.5,46.0], [-4.668,144.707,44.258], color="blue red", name="Arrows_15.1850004196_3")

cluster_dict["15.1850004196"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-2.0), float(140.0), float(47.0), float(1.0)]

cluster_dict["15.1850004196_arrows"] += cgo_arrow([-2.0,140.0,47.0], [-3.278,138.005,45.4], color="blue red", name="Arrows_15.1850004196_4")

cluster_dict["15.1850004196"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(2.0), float(136.5), float(46.5), float(1.0)]

cluster_dict["15.1850004196_arrows"] += cgo_arrow([2.0,136.5,46.5], [-0.043,135.193,45.922], color="blue red", name="Arrows_15.1850004196_5")

cluster_dict["15.1850004196"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(1.69875958824), float(144.332388005), float(45.2705399776), float(1.0)]


cluster_dict["15.1850004196"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-1.0), float(145.0), float(48.0), float(1.0)]

cluster_dict["15.1850004196_arrows"] += cgo_arrow([-1.0,145.0,48.0], [-0.44,144.992,50.646], color="red blue", name="Arrows_15.1850004196_6")

cluster_dict["15.1850004196"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-1.0), float(145.0), float(48.0), float(1.0)]

cluster_dict["15.1850004196_arrows"] += cgo_arrow([-1.0,145.0,48.0], [-0.44,144.992,50.646], color="red blue", name="Arrows_15.1850004196_7")

cluster_dict["15.1850004196"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(0.5), float(140.5), float(47.5), float(1.0)]

cluster_dict["15.1850004196_arrows"] += cgo_arrow([0.5,140.5,47.5], [0.695,140.446,50.066], color="red blue", name="Arrows_15.1850004196_8")

cluster_dict["15.1850004196"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(9.0), float(141.5), float(41.5), float(1.0)]

cluster_dict["15.1850004196_arrows"] += cgo_arrow([9.0,141.5,41.5], [9.366,139.015,42.432], color="red blue", name="Arrows_15.1850004196_9")

cmd.load_cgo(cluster_dict["15.1850004196"], "Features_15.1850004196", 1)
cmd.load_cgo(cluster_dict["15.1850004196_arrows"], "Arrows_15.1850004196")
cmd.set("transparency", 0.2,"Features_15.1850004196")
cmd.group("Pharmacophore_15.1850004196", members="Features_15.1850004196")
cmd.group("Pharmacophore_15.1850004196", members="Arrows_15.1850004196")

if dirpath:
    f = join(dirpath, "8/label_threshold_15.1850004196.mol2")
else:
    f = "8/label_threshold_15.1850004196.mol2"

cmd.load(f, 'label_threshold_15.1850004196')
cmd.hide('everything', 'label_threshold_15.1850004196')
cmd.label("label_threshold_15.1850004196", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.1850004196', members= 'label_threshold_15.1850004196')


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
    f = join(dirpath, "9/label_threshold_12.8.mol2")
else:
    f = "9/label_threshold_12.8.mol2"

cmd.load(f, 'label_threshold_12.8')
cmd.hide('everything', 'label_threshold_12.8')
cmd.label("label_threshold_12.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.8]
gfiles = ['9/apolar.grd', '9/acceptor.grd']
grids = ['apolar', 'acceptor']
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


cluster_dict = {"15.1739997864":[], "15.1739997864_arrows":[]}

cluster_dict["15.1739997864"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-7.0), float(139.5), float(48.0), float(1.0)]

cluster_dict["15.1739997864_arrows"] += cgo_arrow([-7.0,139.5,48.0], [-7.298,137.092,49.558], color="blue red", name="Arrows_15.1739997864_1")

cluster_dict["15.1739997864"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-11.5), float(142.5), float(45.5), float(1.0)]

cluster_dict["15.1739997864_arrows"] += cgo_arrow([-11.5,142.5,45.5], [-9.665,144.622,45.987], color="blue red", name="Arrows_15.1739997864_2")

cluster_dict["15.1739997864"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-11.5), float(136.5), float(43.0), float(1.0)]


cluster_dict["15.1739997864"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-10.0), float(138.5), float(49.5), float(1.0)]

cluster_dict["15.1739997864_arrows"] += cgo_arrow([-10.0,138.5,49.5], [-7.298,137.092,49.558], color="blue red", name="Arrows_15.1739997864_3")

cluster_dict["15.1739997864"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-8.0), float(143.5), float(51.0), float(1.0)]

cluster_dict["15.1739997864_arrows"] += cgo_arrow([-8.0,143.5,51.0], [-8.269,146.312,49.744], color="blue red", name="Arrows_15.1739997864_4")

cluster_dict["15.1739997864"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-26.9405511871), float(140.39145096), float(46.0658136482), float(1.0)]


cluster_dict["15.1739997864"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-28.5029664984), float(146.0), float(46.2588821555), float(1.0)]


cluster_dict["15.1739997864"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-12.9455994644), float(145.605752198), float(40.7263110561), float(1.0)]


cluster_dict["15.1739997864"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-10.1879739038), float(139.689088791), float(47.3621475679), float(1.0)]


cluster_dict["15.1739997864"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-11.5), float(136.5), float(42.0), float(1.0)]


cluster_dict["15.1739997864"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-26.0), float(142.5), float(47.0), float(1.0)]

cluster_dict["15.1739997864_arrows"] += cgo_arrow([-26.0,142.5,47.0], [-23.081,143.083,47.078], color="red blue", name="Arrows_15.1739997864_5")

cluster_dict["15.1739997864"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-11.5), float(143.5), float(49.5), float(1.0)]

cluster_dict["15.1739997864_arrows"] += cgo_arrow([-11.5,143.5,49.5], [-14.448,143.769,50.64], color="red blue", name="Arrows_15.1739997864_6")

cluster_dict["15.1739997864"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-9.5), float(137.0), float(51.0), float(1.0)]

cluster_dict["15.1739997864_arrows"] += cgo_arrow([-9.5,137.0,51.0], [-7.298,137.092,49.558], color="red blue", name="Arrows_15.1739997864_7")

cluster_dict["15.1739997864"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-9.0), float(137.5), float(45.5), float(1.0)]

cluster_dict["15.1739997864_arrows"] += cgo_arrow([-9.0,137.5,45.5], [-7.298,137.092,49.558], color="red blue", name="Arrows_15.1739997864_8")

cluster_dict["15.1739997864"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-6.5), float(138.5), float(48.0), float(1.0)]

cluster_dict["15.1739997864_arrows"] += cgo_arrow([-6.5,138.5,48.0], [-7.298,137.092,49.558], color="red blue", name="Arrows_15.1739997864_9")

cmd.load_cgo(cluster_dict["15.1739997864"], "Features_15.1739997864", 1)
cmd.load_cgo(cluster_dict["15.1739997864_arrows"], "Arrows_15.1739997864")
cmd.set("transparency", 0.2,"Features_15.1739997864")
cmd.group("Pharmacophore_15.1739997864", members="Features_15.1739997864")
cmd.group("Pharmacophore_15.1739997864", members="Arrows_15.1739997864")

if dirpath:
    f = join(dirpath, "9/label_threshold_15.1739997864.mol2")
else:
    f = "9/label_threshold_15.1739997864.mol2"

cmd.load(f, 'label_threshold_15.1739997864')
cmd.hide('everything', 'label_threshold_15.1739997864')
cmd.label("label_threshold_15.1739997864", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.1739997864', members= 'label_threshold_15.1739997864')


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
    f = join(dirpath, "10/label_threshold_13.0.mol2")
else:
    f = "10/label_threshold_13.0.mol2"

cmd.load(f, 'label_threshold_13.0')
cmd.hide('everything', 'label_threshold_13.0')
cmd.label("label_threshold_13.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.0]
gfiles = ['10/apolar.grd', '10/acceptor.grd']
grids = ['apolar', 'acceptor']
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


cluster_dict = {"14.8260002136":[], "14.8260002136_arrows":[]}

cluster_dict["14.8260002136"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-26.0), float(154.0), float(77.0), float(1.0)]

cluster_dict["14.8260002136_arrows"] += cgo_arrow([-26.0,154.0,77.0], [-28.521,154.807,77.064], color="blue red", name="Arrows_14.8260002136_1")

cluster_dict["14.8260002136"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-26.0), float(154.5), float(75.5), float(1.0)]

cluster_dict["14.8260002136_arrows"] += cgo_arrow([-26.0,154.5,75.5], [-28.521,154.807,77.064], color="blue red", name="Arrows_14.8260002136_2")

cluster_dict["14.8260002136"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-19.5), float(147.5), float(76.5), float(1.0)]

cluster_dict["14.8260002136_arrows"] += cgo_arrow([-19.5,147.5,76.5], [-20.783,145.068,75.824], color="blue red", name="Arrows_14.8260002136_3")

cluster_dict["14.8260002136"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-19.5), float(150.5), float(74.5), float(1.0)]

cluster_dict["14.8260002136_arrows"] += cgo_arrow([-19.5,150.5,74.5], [-17.276,152.355,73.761], color="blue red", name="Arrows_14.8260002136_4")

cluster_dict["14.8260002136"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-17.5), float(147.0), float(78.5), float(1.0)]

cluster_dict["14.8260002136_arrows"] += cgo_arrow([-17.5,147.0,78.5], [-16.688,145.544,80.941], color="blue red", name="Arrows_14.8260002136_5")

cluster_dict["14.8260002136"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-22.4969991752), float(152.657217589), float(76.4505571012), float(1.0)]


cluster_dict["14.8260002136"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-17.0714285714), float(149.142857143), float(79.5), float(1.0)]


cluster_dict["14.8260002136"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-27.0), float(155.5), float(74.5), float(1.0)]

cluster_dict["14.8260002136_arrows"] += cgo_arrow([-27.0,155.5,74.5], [-26.779,159.137,75.242], color="red blue", name="Arrows_14.8260002136_6")

cluster_dict["14.8260002136"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-25.5), float(152.5), float(72.0), float(1.0)]

cluster_dict["14.8260002136_arrows"] += cgo_arrow([-25.5,152.5,72.0], [-22.98,153.045,70.263], color="red blue", name="Arrows_14.8260002136_7")

cluster_dict["14.8260002136"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-20.5), float(152.0), float(72.0), float(1.0)]

cluster_dict["14.8260002136_arrows"] += cgo_arrow([-20.5,152.0,72.0], [-22.98,153.045,70.263], color="red blue", name="Arrows_14.8260002136_8")

cluster_dict["14.8260002136"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-17.5), float(150.5), float(77.0), float(1.0)]

cluster_dict["14.8260002136_arrows"] += cgo_arrow([-17.5,150.5,77.0], [-15.287,149.378,76.0], color="red blue", name="Arrows_14.8260002136_9")

cmd.load_cgo(cluster_dict["14.8260002136"], "Features_14.8260002136", 1)
cmd.load_cgo(cluster_dict["14.8260002136_arrows"], "Arrows_14.8260002136")
cmd.set("transparency", 0.2,"Features_14.8260002136")
cmd.group("Pharmacophore_14.8260002136", members="Features_14.8260002136")
cmd.group("Pharmacophore_14.8260002136", members="Arrows_14.8260002136")

if dirpath:
    f = join(dirpath, "10/label_threshold_14.8260002136.mol2")
else:
    f = "10/label_threshold_14.8260002136.mol2"

cmd.load(f, 'label_threshold_14.8260002136')
cmd.hide('everything', 'label_threshold_14.8260002136')
cmd.label("label_threshold_14.8260002136", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.8260002136', members= 'label_threshold_14.8260002136')


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
    f = join(dirpath, "11/label_threshold_6.5.mol2")
else:
    f = "11/label_threshold_6.5.mol2"

cmd.load(f, 'label_threshold_6.5')
cmd.hide('everything', 'label_threshold_6.5')
cmd.label("label_threshold_6.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [6.5]
gfiles = ['11/apolar.grd', '11/acceptor.grd']
grids = ['apolar', 'acceptor']
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


cluster_dict = {"12.5530004501":[], "12.5530004501_arrows":[]}

cluster_dict["12.5530004501"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-37.5), float(136.5), float(76.5), float(1.0)]

cluster_dict["12.5530004501_arrows"] += cgo_arrow([-37.5,136.5,76.5], [-37.706,136.19,79.478], color="blue red", name="Arrows_12.5530004501_1")

cluster_dict["12.5530004501"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-40.5558641236), float(135.672426359), float(74.8464323826), float(1.0)]


cluster_dict["12.5530004501"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-38.0), float(132.5), float(73.5), float(1.0)]

cluster_dict["12.5530004501_arrows"] += cgo_arrow([-38.0,132.5,73.5], [-40.176,130.752,70.647], color="red blue", name="Arrows_12.5530004501_2")

cluster_dict["12.5530004501"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-39.0), float(139.0), float(77.0), float(1.0)]

cluster_dict["12.5530004501_arrows"] += cgo_arrow([-39.0,139.0,77.0], [-38.713,140.071,79.854], color="red blue", name="Arrows_12.5530004501_3")

cmd.load_cgo(cluster_dict["12.5530004501"], "Features_12.5530004501", 1)
cmd.load_cgo(cluster_dict["12.5530004501_arrows"], "Arrows_12.5530004501")
cmd.set("transparency", 0.2,"Features_12.5530004501")
cmd.group("Pharmacophore_12.5530004501", members="Features_12.5530004501")
cmd.group("Pharmacophore_12.5530004501", members="Arrows_12.5530004501")

if dirpath:
    f = join(dirpath, "11/label_threshold_12.5530004501.mol2")
else:
    f = "11/label_threshold_12.5530004501.mol2"

cmd.load(f, 'label_threshold_12.5530004501')
cmd.hide('everything', 'label_threshold_12.5530004501')
cmd.label("label_threshold_12.5530004501", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.5530004501', members= 'label_threshold_12.5530004501')


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
    f = join(dirpath, "12/label_threshold_9.1.mol2")
else:
    f = "12/label_threshold_9.1.mol2"

cmd.load(f, 'label_threshold_9.1')
cmd.hide('everything', 'label_threshold_9.1')
cmd.label("label_threshold_9.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [9.1]
gfiles = ['12/apolar.grd', '12/acceptor.grd']
grids = ['apolar', 'acceptor']
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


cluster_dict = {"11.7399997711":[], "11.7399997711_arrows":[]}

cluster_dict["11.7399997711"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-36.3590398537), float(134.591224067), float(60.7540361526), float(1.0)]


cluster_dict["11.7399997711"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-40.4135250821), float(144.555990698), float(50.9116078561), float(1.0)]


cluster_dict["11.7399997711"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-27.6240596058), float(140.881142199), float(46.7660256512), float(1.0)]


cluster_dict["11.7399997711"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-28.9778507746), float(146.0), float(46.8208066201), float(1.0)]


cluster_dict["11.7399997711"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-26.7445671909), float(146.0), float(51.2675887404), float(1.0)]


cluster_dict["11.7399997711"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-20.397016053), float(140.28619316), float(50.1016262645), float(1.0)]


cluster_dict["11.7399997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-40.5), float(135.5), float(64.5), float(1.0)]

cluster_dict["11.7399997711_arrows"] += cgo_arrow([-40.5,135.5,64.5], [-40.92,137.025,61.157], color="red blue", name="Arrows_11.7399997711_1")

cluster_dict["11.7399997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-26.0), float(142.5), float(47.0), float(1.0)]

cluster_dict["11.7399997711_arrows"] += cgo_arrow([-26.0,142.5,47.0], [-23.081,143.083,47.078], color="red blue", name="Arrows_11.7399997711_2")

cmd.load_cgo(cluster_dict["11.7399997711"], "Features_11.7399997711", 1)
cmd.load_cgo(cluster_dict["11.7399997711_arrows"], "Arrows_11.7399997711")
cmd.set("transparency", 0.2,"Features_11.7399997711")
cmd.group("Pharmacophore_11.7399997711", members="Features_11.7399997711")
cmd.group("Pharmacophore_11.7399997711", members="Arrows_11.7399997711")

if dirpath:
    f = join(dirpath, "12/label_threshold_11.7399997711.mol2")
else:
    f = "12/label_threshold_11.7399997711.mol2"

cmd.load(f, 'label_threshold_11.7399997711')
cmd.hide('everything', 'label_threshold_11.7399997711')
cmd.label("label_threshold_11.7399997711", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.7399997711', members= 'label_threshold_11.7399997711')


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
    f = join(dirpath, "13/label_threshold_9.3.mol2")
else:
    f = "13/label_threshold_9.3.mol2"

cmd.load(f, 'label_threshold_9.3')
cmd.hide('everything', 'label_threshold_9.3')
cmd.label("label_threshold_9.3", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [9.3]
gfiles = ['13/apolar.grd', '13/acceptor.grd']
grids = ['apolar', 'acceptor']
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


cluster_dict = {"11.3380002975":[], "11.3380002975_arrows":[]}

cluster_dict["11.3380002975"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-30.9553105114), float(137.724998238), float(56.7253355522), float(1.0)]


cluster_dict["11.3380002975"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-24.005307814), float(139.72520643), float(46.4594932182), float(1.0)]


cluster_dict["11.3380002975"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-28.8453138418), float(146.0), float(46.4567950638), float(1.0)]


cluster_dict["11.3380002975"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-26.4583218129), float(146.388952014), float(51.2510394384), float(1.0)]


cluster_dict["11.3380002975"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-26.0), float(142.5), float(47.0), float(1.0)]

cluster_dict["11.3380002975_arrows"] += cgo_arrow([-26.0,142.5,47.0], [-23.081,143.083,47.078], color="red blue", name="Arrows_11.3380002975_1")

cmd.load_cgo(cluster_dict["11.3380002975"], "Features_11.3380002975", 1)
cmd.load_cgo(cluster_dict["11.3380002975_arrows"], "Arrows_11.3380002975")
cmd.set("transparency", 0.2,"Features_11.3380002975")
cmd.group("Pharmacophore_11.3380002975", members="Features_11.3380002975")
cmd.group("Pharmacophore_11.3380002975", members="Arrows_11.3380002975")

if dirpath:
    f = join(dirpath, "13/label_threshold_11.3380002975.mol2")
else:
    f = "13/label_threshold_11.3380002975.mol2"

cmd.load(f, 'label_threshold_11.3380002975')
cmd.hide('everything', 'label_threshold_11.3380002975')
cmd.label("label_threshold_11.3380002975", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.3380002975', members= 'label_threshold_11.3380002975')


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
    f = join(dirpath, "14/label_threshold_6.3.mol2")
else:
    f = "14/label_threshold_6.3.mol2"

cmd.load(f, 'label_threshold_6.3')
cmd.hide('everything', 'label_threshold_6.3')
cmd.label("label_threshold_6.3", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [6.3]
gfiles = ['14/apolar.grd', '14/acceptor.grd']
grids = ['apolar', 'acceptor']
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


cluster_dict = {"11.0869998932":[], "11.0869998932_arrows":[]}

cluster_dict["11.0869998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(1.0), float(123.5), float(49.5), float(1.0)]

cluster_dict["11.0869998932_arrows"] += cgo_arrow([1.0,123.5,49.5], [2.853,123.406,47.025], color="blue red", name="Arrows_11.0869998932_1")

cluster_dict["11.0869998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(5.0), float(125.0), float(52.0), float(1.0)]

cluster_dict["11.0869998932_arrows"] += cgo_arrow([5.0,125.0,52.0], [5.551,127.715,52.611], color="blue red", name="Arrows_11.0869998932_2")

cluster_dict["11.0869998932"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(1.75786417079), float(122.698830399), float(50.205445228), float(1.0)]


cmd.load_cgo(cluster_dict["11.0869998932"], "Features_11.0869998932", 1)
cmd.load_cgo(cluster_dict["11.0869998932_arrows"], "Arrows_11.0869998932")
cmd.set("transparency", 0.2,"Features_11.0869998932")
cmd.group("Pharmacophore_11.0869998932", members="Features_11.0869998932")
cmd.group("Pharmacophore_11.0869998932", members="Arrows_11.0869998932")

if dirpath:
    f = join(dirpath, "14/label_threshold_11.0869998932.mol2")
else:
    f = "14/label_threshold_11.0869998932.mol2"

cmd.load(f, 'label_threshold_11.0869998932')
cmd.hide('everything', 'label_threshold_11.0869998932')
cmd.label("label_threshold_11.0869998932", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.0869998932', members= 'label_threshold_11.0869998932')


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
    f = join(dirpath, "15/label_threshold_0.9.mol2")
else:
    f = "15/label_threshold_0.9.mol2"

cmd.load(f, 'label_threshold_0.9')
cmd.hide('everything', 'label_threshold_0.9')
cmd.label("label_threshold_0.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.9]
gfiles = ['15/apolar.grd', '15/acceptor.grd']
grids = ['apolar', 'acceptor']
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


cluster_dict = {"10.220000267":[], "10.220000267_arrows":[]}

cluster_dict["10.220000267"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-13.4127539059), float(148.255730326), float(39.1706906184), float(1.0)]


cmd.load_cgo(cluster_dict["10.220000267"], "Features_10.220000267", 1)
cmd.load_cgo(cluster_dict["10.220000267_arrows"], "Arrows_10.220000267")
cmd.set("transparency", 0.2,"Features_10.220000267")
cmd.group("Pharmacophore_10.220000267", members="Features_10.220000267")
cmd.group("Pharmacophore_10.220000267", members="Arrows_10.220000267")

if dirpath:
    f = join(dirpath, "15/label_threshold_10.220000267.mol2")
else:
    f = "15/label_threshold_10.220000267.mol2"

cmd.load(f, 'label_threshold_10.220000267')
cmd.hide('everything', 'label_threshold_10.220000267')
cmd.label("label_threshold_10.220000267", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.220000267', members= 'label_threshold_10.220000267')


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
    f = join(dirpath, "16/label_threshold_0.6.mol2")
else:
    f = "16/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
gfiles = ['16/apolar.grd', '16/acceptor.grd']
grids = ['apolar', 'acceptor']
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


cluster_dict = {"10.1850004196":[], "10.1850004196_arrows":[]}

cluster_dict["10.1850004196"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(13.0), float(121.5), float(43.0), float(1.0)]

cluster_dict["10.1850004196_arrows"] += cgo_arrow([13.0,121.5,43.0], [12.604,122.713,40.402], color="blue red", name="Arrows_10.1850004196_1")

cluster_dict["10.1850004196"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(13.6399607035), float(118.338722938), float(42.0016045726), float(1.0)]


cluster_dict["10.1850004196"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(121.5), float(41.5), float(1.0)]

cluster_dict["10.1850004196_arrows"] += cgo_arrow([15.0,121.5,41.5], [12.604,122.713,40.402], color="red blue", name="Arrows_10.1850004196_2")

cmd.load_cgo(cluster_dict["10.1850004196"], "Features_10.1850004196", 1)
cmd.load_cgo(cluster_dict["10.1850004196_arrows"], "Arrows_10.1850004196")
cmd.set("transparency", 0.2,"Features_10.1850004196")
cmd.group("Pharmacophore_10.1850004196", members="Features_10.1850004196")
cmd.group("Pharmacophore_10.1850004196", members="Arrows_10.1850004196")

if dirpath:
    f = join(dirpath, "16/label_threshold_10.1850004196.mol2")
else:
    f = "16/label_threshold_10.1850004196.mol2"

cmd.load(f, 'label_threshold_10.1850004196')
cmd.hide('everything', 'label_threshold_10.1850004196')
cmd.label("label_threshold_10.1850004196", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.1850004196', members= 'label_threshold_10.1850004196')


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
    f = join(dirpath, "17/label_threshold_5.5.mol2")
else:
    f = "17/label_threshold_5.5.mol2"

cmd.load(f, 'label_threshold_5.5')
cmd.hide('everything', 'label_threshold_5.5')
cmd.label("label_threshold_5.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [5.5]
gfiles = ['17/apolar.grd', '17/acceptor.grd']
grids = ['apolar', 'acceptor']
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


cluster_dict = {"10.1160001755":[], "10.1160001755_arrows":[]}

cluster_dict["10.1160001755"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-37.0679347947), float(134.158597856), float(61.0610843517), float(1.0)]


cluster_dict["10.1160001755"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-41.0), float(135.0), float(64.5), float(1.0)]

cluster_dict["10.1160001755_arrows"] += cgo_arrow([-41.0,135.0,64.5], [-40.92,137.025,61.157], color="red blue", name="Arrows_10.1160001755_1")

cmd.load_cgo(cluster_dict["10.1160001755"], "Features_10.1160001755", 1)
cmd.load_cgo(cluster_dict["10.1160001755_arrows"], "Arrows_10.1160001755")
cmd.set("transparency", 0.2,"Features_10.1160001755")
cmd.group("Pharmacophore_10.1160001755", members="Features_10.1160001755")
cmd.group("Pharmacophore_10.1160001755", members="Arrows_10.1160001755")

if dirpath:
    f = join(dirpath, "17/label_threshold_10.1160001755.mol2")
else:
    f = "17/label_threshold_10.1160001755.mol2"

cmd.load(f, 'label_threshold_10.1160001755')
cmd.hide('everything', 'label_threshold_10.1160001755')
cmd.label("label_threshold_10.1160001755", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.1160001755', members= 'label_threshold_10.1160001755')


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
    f = join(dirpath, "18/label_threshold_30.0.mol2")
else:
    f = "18/label_threshold_30.0.mol2"

cmd.load(f, 'label_threshold_30.0')
cmd.hide('everything', 'label_threshold_30.0')
cmd.label("label_threshold_30.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [30.0]
gfiles = ['18/apolar.grd', '18/acceptor.grd']
grids = ['apolar', 'acceptor']
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


cluster_dict = {"0":[], "0_arrows":[]}

cmd.load_cgo(cluster_dict["0"], "Features_0", 1)
cmd.load_cgo(cluster_dict["0_arrows"], "Arrows_0")
cmd.set("transparency", 0.2,"Features_0")
cmd.group("Pharmacophore_0", members="Features_0")
cmd.group("Pharmacophore_0", members="Arrows_0")

if dirpath:
    f = join(dirpath, "18/label_threshold_0.mol2")
else:
    f = "18/label_threshold_0.mol2"

cmd.load(f, 'label_threshold_0')
cmd.hide('everything', 'label_threshold_0')
cmd.label("label_threshold_0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_0', members= 'label_threshold_0')


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

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
