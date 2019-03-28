
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


cluster_dict = {"10.295999527":[], "10.295999527_arrows":[]}

cluster_dict["10.295999527"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(6.40536656594), float(-20.9986897702), float(2.78902959972), float(1.0)]


cluster_dict["10.295999527"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(10.7914084108), float(-26.9377889718), float(-1.44661388676), float(1.0)]


cmd.load_cgo(cluster_dict["10.295999527"], "Features_10.295999527", 1)
cmd.load_cgo(cluster_dict["10.295999527_arrows"], "Arrows_10.295999527")
cmd.set("transparency", 0.2,"Features_10.295999527")
cmd.group("Pharmacophore_10.295999527", members="Features_10.295999527")
cmd.group("Pharmacophore_10.295999527", members="Arrows_10.295999527")

if dirpath:
    f = join(dirpath, "0/label_threshold_10.295999527.mol2")
else:
    f = "0/label_threshold_10.295999527.mol2"

cmd.load(f, 'label_threshold_10.295999527')
cmd.hide('everything', 'label_threshold_10.295999527')
cmd.label("label_threshold_10.295999527", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.295999527', members= 'label_threshold_10.295999527')


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
    f = join(dirpath, "1/label_threshold_19.0.mol2")
else:
    f = "1/label_threshold_19.0.mol2"

cmd.load(f, 'label_threshold_19.0')
cmd.hide('everything', 'label_threshold_19.0')
cmd.label("label_threshold_19.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [19.0]
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


cluster_dict = {"21.2399997711":[], "21.2399997711_arrows":[]}

cluster_dict["21.2399997711"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-4.0), float(-44.5), float(6.5), float(1.0)]

cluster_dict["21.2399997711_arrows"] += cgo_arrow([-4.0,-44.5,6.5], [-7.136,-43.562,7.403], color="blue red", name="Arrows_21.2399997711_1")

cluster_dict["21.2399997711"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-2.0), float(-39.5), float(0.5), float(1.0)]

cluster_dict["21.2399997711_arrows"] += cgo_arrow([-2.0,-39.5,0.5], [-4.407,-37.137,0.543], color="blue red", name="Arrows_21.2399997711_2")

cluster_dict["21.2399997711"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-2.12318469976), float(-41.5472484872), float(3.11201188107), float(1.0)]


cluster_dict["21.2399997711"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(9.79620486447), float(-44.4395949043), float(3.05920881905), float(1.0)]


cluster_dict["21.2399997711"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(4.89715316027), float(-38.9933158471), float(-0.257979922911), float(1.0)]


cluster_dict["21.2399997711"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(8.7), float(-27.0), float(2.9), float(1.0)]


cluster_dict["21.2399997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-5.0), float(-44.0), float(4.0), float(1.0)]

cluster_dict["21.2399997711_arrows"] += cgo_arrow([-5.0,-44.0,4.0], [-3.912,-44.155,1.825], color="red blue", name="Arrows_21.2399997711_3")

cluster_dict["21.2399997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-4.5), float(-38.0), float(2.5), float(1.0)]

cluster_dict["21.2399997711_arrows"] += cgo_arrow([-4.5,-38.0,2.5], [-4.407,-37.137,0.543], color="red blue", name="Arrows_21.2399997711_4")

cluster_dict["21.2399997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-3.0), float(-39.5), float(0.5), float(1.0)]

cluster_dict["21.2399997711_arrows"] += cgo_arrow([-3.0,-39.5,0.5], [-5.681,-39.347,0.026], color="red blue", name="Arrows_21.2399997711_5")

cluster_dict["21.2399997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(2.5), float(-37.0), float(-5.5), float(1.0)]

cluster_dict["21.2399997711_arrows"] += cgo_arrow([2.5,-37.0,-5.5], [1.834,-38.561,-7.841], color="red blue", name="Arrows_21.2399997711_6")

cluster_dict["21.2399997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(4.0), float(-37.5), float(0.0), float(1.0)]

cluster_dict["21.2399997711_arrows"] += cgo_arrow([4.0,-37.5,0.0], [2.326,-39.824,2.141], color="red blue", name="Arrows_21.2399997711_7")

cluster_dict["21.2399997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(10.5), float(-29.0), float(-2.0), float(1.0)]

cluster_dict["21.2399997711_arrows"] += cgo_arrow([10.5,-29.0,-2.0], [10.383,-26.287,-2.934], color="red blue", name="Arrows_21.2399997711_8")

cluster_dict["21.2399997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(10.5), float(-46.5), float(4.5), float(1.0)]

cluster_dict["21.2399997711_arrows"] += cgo_arrow([10.5,-46.5,4.5], [9.393,-45.27,7.582], color="red blue", name="Arrows_21.2399997711_9")

cmd.load_cgo(cluster_dict["21.2399997711"], "Features_21.2399997711", 1)
cmd.load_cgo(cluster_dict["21.2399997711_arrows"], "Arrows_21.2399997711")
cmd.set("transparency", 0.2,"Features_21.2399997711")
cmd.group("Pharmacophore_21.2399997711", members="Features_21.2399997711")
cmd.group("Pharmacophore_21.2399997711", members="Arrows_21.2399997711")

if dirpath:
    f = join(dirpath, "1/label_threshold_21.2399997711.mol2")
else:
    f = "1/label_threshold_21.2399997711.mol2"

cmd.load(f, 'label_threshold_21.2399997711')
cmd.hide('everything', 'label_threshold_21.2399997711')
cmd.label("label_threshold_21.2399997711", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_21.2399997711', members= 'label_threshold_21.2399997711')


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
    f = join(dirpath, "2/label_threshold_14.4.mol2")
else:
    f = "2/label_threshold_14.4.mol2"

cmd.load(f, 'label_threshold_14.4')
cmd.hide('everything', 'label_threshold_14.4')
cmd.label("label_threshold_14.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [14.4]
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


cluster_dict = {"20.3339996338":[], "20.3339996338_arrows":[]}

cluster_dict["20.3339996338"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(10.1137004876), float(-44.5667574296), float(3.11379246808), float(1.0)]


cluster_dict["20.3339996338"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(4.80409933888), float(-39.2850012741), float(-0.142745572068), float(1.0)]


cluster_dict["20.3339996338"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(5.0), float(-44.5), float(4.5), float(1.0)]

cluster_dict["20.3339996338_arrows"] += cgo_arrow([5.0,-44.5,4.5], [5.589,-47.905,4.4], color="red blue", name="Arrows_20.3339996338_1")

cluster_dict["20.3339996338"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(10.5), float(-46.5), float(4.5), float(1.0)]

cluster_dict["20.3339996338_arrows"] += cgo_arrow([10.5,-46.5,4.5], [9.393,-45.27,7.582], color="red blue", name="Arrows_20.3339996338_2")

cmd.load_cgo(cluster_dict["20.3339996338"], "Features_20.3339996338", 1)
cmd.load_cgo(cluster_dict["20.3339996338_arrows"], "Arrows_20.3339996338")
cmd.set("transparency", 0.2,"Features_20.3339996338")
cmd.group("Pharmacophore_20.3339996338", members="Features_20.3339996338")
cmd.group("Pharmacophore_20.3339996338", members="Arrows_20.3339996338")

if dirpath:
    f = join(dirpath, "2/label_threshold_20.3339996338.mol2")
else:
    f = "2/label_threshold_20.3339996338.mol2"

cmd.load(f, 'label_threshold_20.3339996338')
cmd.hide('everything', 'label_threshold_20.3339996338')
cmd.label("label_threshold_20.3339996338", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_20.3339996338', members= 'label_threshold_20.3339996338')


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
    f = join(dirpath, "3/label_threshold_17.6.mol2")
else:
    f = "3/label_threshold_17.6.mol2"

cmd.load(f, 'label_threshold_17.6')
cmd.hide('everything', 'label_threshold_17.6')
cmd.label("label_threshold_17.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [17.6]
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


cluster_dict = {"18.9740009308":[], "18.9740009308_arrows":[]}

cluster_dict["18.9740009308"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-4.0), float(-44.5), float(7.0), float(1.0)]

cluster_dict["18.9740009308_arrows"] += cgo_arrow([-4.0,-44.5,7.0], [-7.136,-43.562,7.403], color="blue red", name="Arrows_18.9740009308_1")

cluster_dict["18.9740009308"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-2.0), float(-39.5), float(0.5), float(1.0)]

cluster_dict["18.9740009308_arrows"] += cgo_arrow([-2.0,-39.5,0.5], [-4.407,-37.137,0.543], color="blue red", name="Arrows_18.9740009308_2")

cluster_dict["18.9740009308"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-2.05431136089), float(-41.6408186531), float(2.96766686819), float(1.0)]


cluster_dict["18.9740009308"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(5.74235335986), float(-35.7793357531), float(-2.29525944957), float(1.0)]


cluster_dict["18.9740009308"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(8.11574611749), float(-44.6980671328), float(2.83163684302), float(1.0)]


cluster_dict["18.9740009308"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-5.0), float(-44.0), float(4.5), float(1.0)]

cluster_dict["18.9740009308_arrows"] += cgo_arrow([-5.0,-44.0,4.5], [-3.912,-44.155,1.825], color="red blue", name="Arrows_18.9740009308_3")

cluster_dict["18.9740009308"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-4.5), float(-43.0), float(8.0), float(1.0)]

cluster_dict["18.9740009308_arrows"] += cgo_arrow([-4.5,-43.0,8.0], [-7.161,-44.188,9.572], color="red blue", name="Arrows_18.9740009308_4")

cluster_dict["18.9740009308"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-4.5), float(-38.0), float(2.5), float(1.0)]

cluster_dict["18.9740009308_arrows"] += cgo_arrow([-4.5,-38.0,2.5], [-4.407,-37.137,0.543], color="red blue", name="Arrows_18.9740009308_5")

cluster_dict["18.9740009308"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-3.0), float(-39.5), float(0.5), float(1.0)]

cluster_dict["18.9740009308_arrows"] += cgo_arrow([-3.0,-39.5,0.5], [-5.681,-39.347,0.026], color="red blue", name="Arrows_18.9740009308_6")

cluster_dict["18.9740009308"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(2.5), float(-37.0), float(-5.5), float(1.0)]

cluster_dict["18.9740009308_arrows"] += cgo_arrow([2.5,-37.0,-5.5], [1.834,-38.561,-7.841], color="red blue", name="Arrows_18.9740009308_7")

cluster_dict["18.9740009308"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(4.0), float(-37.5), float(0.0), float(1.0)]

cluster_dict["18.9740009308_arrows"] += cgo_arrow([4.0,-37.5,0.0], [2.326,-39.824,2.141], color="red blue", name="Arrows_18.9740009308_8")

cmd.load_cgo(cluster_dict["18.9740009308"], "Features_18.9740009308", 1)
cmd.load_cgo(cluster_dict["18.9740009308_arrows"], "Arrows_18.9740009308")
cmd.set("transparency", 0.2,"Features_18.9740009308")
cmd.group("Pharmacophore_18.9740009308", members="Features_18.9740009308")
cmd.group("Pharmacophore_18.9740009308", members="Arrows_18.9740009308")

if dirpath:
    f = join(dirpath, "3/label_threshold_18.9740009308.mol2")
else:
    f = "3/label_threshold_18.9740009308.mol2"

cmd.load(f, 'label_threshold_18.9740009308')
cmd.hide('everything', 'label_threshold_18.9740009308')
cmd.label("label_threshold_18.9740009308", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_18.9740009308', members= 'label_threshold_18.9740009308')


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
    f = join(dirpath, "4/label_threshold_15.8.mol2")
else:
    f = "4/label_threshold_15.8.mol2"

cmd.load(f, 'label_threshold_15.8')
cmd.hide('everything', 'label_threshold_15.8')
cmd.label("label_threshold_15.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.8]
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


cluster_dict = {"18.6779994965":[], "18.6779994965_arrows":[]}

cluster_dict["18.6779994965"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-8.0), float(-37.5), float(6.0), float(1.0)]

cluster_dict["18.6779994965_arrows"] += cgo_arrow([-8.0,-37.5,6.0], [-9.247,-38.712,8.549], color="blue red", name="Arrows_18.6779994965_1")

cluster_dict["18.6779994965"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-4.0), float(-44.5), float(7.0), float(1.0)]

cluster_dict["18.6779994965_arrows"] += cgo_arrow([-4.0,-44.5,7.0], [-7.136,-43.562,7.403], color="blue red", name="Arrows_18.6779994965_2")

cluster_dict["18.6779994965"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-2.0), float(-39.5), float(0.5), float(1.0)]

cluster_dict["18.6779994965_arrows"] += cgo_arrow([-2.0,-39.5,0.5], [-4.407,-37.137,0.543], color="blue red", name="Arrows_18.6779994965_3")

cluster_dict["18.6779994965"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(2.49444228124), float(-40.3160216237), float(1.43586988129), float(1.0)]


cluster_dict["18.6779994965"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(6.53935850538), float(-43.8659859969), float(3.17725234624), float(1.0)]


cluster_dict["18.6779994965"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-4.5), float(-38.0), float(2.5), float(1.0)]

cluster_dict["18.6779994965_arrows"] += cgo_arrow([-4.5,-38.0,2.5], [-4.407,-37.137,0.543], color="red blue", name="Arrows_18.6779994965_4")

cluster_dict["18.6779994965"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-5.0), float(-44.0), float(4.5), float(1.0)]

cluster_dict["18.6779994965_arrows"] += cgo_arrow([-5.0,-44.0,4.5], [-3.912,-44.155,1.825], color="red blue", name="Arrows_18.6779994965_5")

cluster_dict["18.6779994965"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-4.5), float(-45.5), float(9.0), float(1.0)]

cluster_dict["18.6779994965_arrows"] += cgo_arrow([-4.5,-45.5,9.0], [-4.214,-46.08,11.547], color="red blue", name="Arrows_18.6779994965_6")

cluster_dict["18.6779994965"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-3.0), float(-39.5), float(0.5), float(1.0)]

cluster_dict["18.6779994965_arrows"] += cgo_arrow([-3.0,-39.5,0.5], [-5.681,-39.347,0.026], color="red blue", name="Arrows_18.6779994965_7")

cluster_dict["18.6779994965"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(2.5), float(-37.0), float(-5.5), float(1.0)]

cluster_dict["18.6779994965_arrows"] += cgo_arrow([2.5,-37.0,-5.5], [1.834,-38.561,-7.841], color="red blue", name="Arrows_18.6779994965_8")

cluster_dict["18.6779994965"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(4.0), float(-37.5), float(0.0), float(1.0)]

cluster_dict["18.6779994965_arrows"] += cgo_arrow([4.0,-37.5,0.0], [2.326,-39.824,2.141], color="red blue", name="Arrows_18.6779994965_9")

cmd.load_cgo(cluster_dict["18.6779994965"], "Features_18.6779994965", 1)
cmd.load_cgo(cluster_dict["18.6779994965_arrows"], "Arrows_18.6779994965")
cmd.set("transparency", 0.2,"Features_18.6779994965")
cmd.group("Pharmacophore_18.6779994965", members="Features_18.6779994965")
cmd.group("Pharmacophore_18.6779994965", members="Arrows_18.6779994965")

if dirpath:
    f = join(dirpath, "4/label_threshold_18.6779994965.mol2")
else:
    f = "4/label_threshold_18.6779994965.mol2"

cmd.load(f, 'label_threshold_18.6779994965')
cmd.hide('everything', 'label_threshold_18.6779994965')
cmd.label("label_threshold_18.6779994965", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_18.6779994965', members= 'label_threshold_18.6779994965')


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
    f = join(dirpath, "5/label_threshold_15.4.mol2")
else:
    f = "5/label_threshold_15.4.mol2"

cmd.load(f, 'label_threshold_15.4')
cmd.hide('everything', 'label_threshold_15.4')
cmd.label("label_threshold_15.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.4]
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


cluster_dict = {"17.3990001678":[], "17.3990001678_arrows":[]}

cluster_dict["17.3990001678"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(12.0), float(-31.5), float(-1.0), float(1.0)]

cluster_dict["17.3990001678_arrows"] += cgo_arrow([12.0,-31.5,-1.0], [13.926,-29.478,0.721], color="blue red", name="Arrows_17.3990001678_1")

cluster_dict["17.3990001678"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(12.5), float(-29.5), float(-3.0), float(1.0)]

cluster_dict["17.3990001678_arrows"] += cgo_arrow([12.5,-29.5,-3.0], [10.383,-26.287,-2.934], color="blue red", name="Arrows_17.3990001678_2")

cluster_dict["17.3990001678"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(8.31206521703), float(-33.8074602349), float(-2.58367793808), float(1.0)]


cluster_dict["17.3990001678"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(8.03348511527), float(-25.9319927391), float(2.47118941425), float(1.0)]


cluster_dict["17.3990001678"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(11.2236057806), float(-41.7766649153), float(2.16609473326), float(1.0)]


cluster_dict["17.3990001678"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(10.852736094), float(-27.0), float(-1.69081454744), float(1.0)]


cluster_dict["17.3990001678"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(2.5), float(-37.0), float(-5.5), float(1.0)]

cluster_dict["17.3990001678_arrows"] += cgo_arrow([2.5,-37.0,-5.5], [1.834,-38.561,-7.841], color="red blue", name="Arrows_17.3990001678_3")

cluster_dict["17.3990001678"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(4.0), float(-37.5), float(0.0), float(1.0)]

cluster_dict["17.3990001678_arrows"] += cgo_arrow([4.0,-37.5,0.0], [2.326,-39.824,2.141], color="red blue", name="Arrows_17.3990001678_4")

cluster_dict["17.3990001678"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(10.5), float(-29.0), float(-2.0), float(1.0)]

cluster_dict["17.3990001678_arrows"] += cgo_arrow([10.5,-29.0,-2.0], [10.383,-26.287,-2.934], color="red blue", name="Arrows_17.3990001678_5")

cluster_dict["17.3990001678"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(11.5), float(-35.0), float(-2.5), float(1.0)]

cluster_dict["17.3990001678_arrows"] += cgo_arrow([11.5,-35.0,-2.5], [10.478,-37.477,-0.934], color="red blue", name="Arrows_17.3990001678_6")

cmd.load_cgo(cluster_dict["17.3990001678"], "Features_17.3990001678", 1)
cmd.load_cgo(cluster_dict["17.3990001678_arrows"], "Arrows_17.3990001678")
cmd.set("transparency", 0.2,"Features_17.3990001678")
cmd.group("Pharmacophore_17.3990001678", members="Features_17.3990001678")
cmd.group("Pharmacophore_17.3990001678", members="Arrows_17.3990001678")

if dirpath:
    f = join(dirpath, "5/label_threshold_17.3990001678.mol2")
else:
    f = "5/label_threshold_17.3990001678.mol2"

cmd.load(f, 'label_threshold_17.3990001678')
cmd.hide('everything', 'label_threshold_17.3990001678')
cmd.label("label_threshold_17.3990001678", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.3990001678', members= 'label_threshold_17.3990001678')


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
    f = join(dirpath, "6/label_threshold_7.4.mol2")
else:
    f = "6/label_threshold_7.4.mol2"

cmd.load(f, 'label_threshold_7.4')
cmd.hide('everything', 'label_threshold_7.4')
cmd.label("label_threshold_7.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [7.4]
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


cluster_dict = {"11.7550001144":[], "11.7550001144_arrows":[]}

cluster_dict["11.7550001144"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(6.5), float(-33.0), float(26.5), float(1.0)]

cluster_dict["11.7550001144_arrows"] += cgo_arrow([6.5,-33.0,26.5], [4.28,-33.309,24.424], color="blue red", name="Arrows_11.7550001144_1")

cluster_dict["11.7550001144"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(6.5), float(-33.0), float(26.5), float(1.0)]

cluster_dict["11.7550001144_arrows"] += cgo_arrow([6.5,-33.0,26.5], [4.28,-33.309,24.424], color="blue red", name="Arrows_11.7550001144_2")

cluster_dict["11.7550001144"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(9.5), float(-37.5), float(28.0), float(1.0)]

cluster_dict["11.7550001144_arrows"] += cgo_arrow([9.5,-37.5,28.0], [8.025,-39.643,26.61], color="blue red", name="Arrows_11.7550001144_3")

cluster_dict["11.7550001144"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(12.5), float(-33.5), float(29.5), float(1.0)]

cluster_dict["11.7550001144_arrows"] += cgo_arrow([12.5,-33.5,29.5], [15.021,-33.063,27.885], color="blue red", name="Arrows_11.7550001144_4")

cluster_dict["11.7550001144"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(8.71418621207), float(-33.8371919608), float(27.4878051013), float(1.0)]


cluster_dict["11.7550001144"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(7.49644974953), float(-30.7267821791), float(25.7732178209), float(1.0)]


cluster_dict["11.7550001144"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(6.5), float(-34.0), float(23.5), float(1.0)]

cluster_dict["11.7550001144_arrows"] += cgo_arrow([6.5,-34.0,23.5], [4.619,-32.866,22.306], color="red blue", name="Arrows_11.7550001144_5")

cluster_dict["11.7550001144"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(9.5), float(-37.5), float(26.5), float(1.0)]

cluster_dict["11.7550001144_arrows"] += cgo_arrow([9.5,-37.5,26.5], [11.665,-39.071,26.398], color="red blue", name="Arrows_11.7550001144_6")

cluster_dict["11.7550001144"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(12.5), float(-35.5), float(29.0), float(1.0)]

cluster_dict["11.7550001144_arrows"] += cgo_arrow([12.5,-35.5,29.0], [15.074,-35.309,27.791], color="red blue", name="Arrows_11.7550001144_7")

cmd.load_cgo(cluster_dict["11.7550001144"], "Features_11.7550001144", 1)
cmd.load_cgo(cluster_dict["11.7550001144_arrows"], "Arrows_11.7550001144")
cmd.set("transparency", 0.2,"Features_11.7550001144")
cmd.group("Pharmacophore_11.7550001144", members="Features_11.7550001144")
cmd.group("Pharmacophore_11.7550001144", members="Arrows_11.7550001144")

if dirpath:
    f = join(dirpath, "6/label_threshold_11.7550001144.mol2")
else:
    f = "6/label_threshold_11.7550001144.mol2"

cmd.load(f, 'label_threshold_11.7550001144')
cmd.hide('everything', 'label_threshold_11.7550001144')
cmd.label("label_threshold_11.7550001144", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.7550001144', members= 'label_threshold_11.7550001144')


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
    f = join(dirpath, "7/label_threshold_5.0.mol2")
else:
    f = "7/label_threshold_5.0.mol2"

cmd.load(f, 'label_threshold_5.0')
cmd.hide('everything', 'label_threshold_5.0')
cmd.label("label_threshold_5.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [5.0]
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


cluster_dict = {"10.295999527":[], "10.295999527_arrows":[]}

cluster_dict["10.295999527"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(6.40536656594), float(-20.9986897702), float(2.78902959972), float(1.0)]


cluster_dict["10.295999527"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(10.7671588202), float(-26.9486122171), float(-1.00894519217), float(1.0)]


cmd.load_cgo(cluster_dict["10.295999527"], "Features_10.295999527", 1)
cmd.load_cgo(cluster_dict["10.295999527_arrows"], "Arrows_10.295999527")
cmd.set("transparency", 0.2,"Features_10.295999527")
cmd.group("Pharmacophore_10.295999527", members="Features_10.295999527")
cmd.group("Pharmacophore_10.295999527", members="Arrows_10.295999527")

if dirpath:
    f = join(dirpath, "7/label_threshold_10.295999527.mol2")
else:
    f = "7/label_threshold_10.295999527.mol2"

cmd.load(f, 'label_threshold_10.295999527')
cmd.hide('everything', 'label_threshold_10.295999527')
cmd.label("label_threshold_10.295999527", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.295999527', members= 'label_threshold_10.295999527')


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
    f = join(dirpath, "8/label_threshold_0.6.mol2")
else:
    f = "8/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
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


cluster_dict = {"5.40899991989":[], "5.40899991989_arrows":[]}

cluster_dict["5.40899991989"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(2.82790552902), float(-17.9995657092), float(18.9869188077), float(1.0)]


cmd.load_cgo(cluster_dict["5.40899991989"], "Features_5.40899991989", 1)
cmd.load_cgo(cluster_dict["5.40899991989_arrows"], "Arrows_5.40899991989")
cmd.set("transparency", 0.2,"Features_5.40899991989")
cmd.group("Pharmacophore_5.40899991989", members="Features_5.40899991989")
cmd.group("Pharmacophore_5.40899991989", members="Arrows_5.40899991989")

if dirpath:
    f = join(dirpath, "8/label_threshold_5.40899991989.mol2")
else:
    f = "8/label_threshold_5.40899991989.mol2"

cmd.load(f, 'label_threshold_5.40899991989')
cmd.hide('everything', 'label_threshold_5.40899991989')
cmd.label("label_threshold_5.40899991989", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_5.40899991989', members= 'label_threshold_5.40899991989')


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

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
