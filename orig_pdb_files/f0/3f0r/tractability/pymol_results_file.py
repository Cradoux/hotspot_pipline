
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
    f = join(dirpath, "0/label_threshold_7.3.mol2")
else:
    f = "0/label_threshold_7.3.mol2"

cmd.load(f, 'label_threshold_7.3')
cmd.hide('everything', 'label_threshold_7.3')
cmd.label("label_threshold_7.3", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [7.3]
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


cluster_dict = {"14.0630002022":[], "14.0630002022_arrows":[]}

cluster_dict["14.0630002022"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(32.5), float(-6.5), float(6.5), float(1.0)]

cluster_dict["14.0630002022_arrows"] += cgo_arrow([32.5,-6.5,6.5], [31.839,-6.19,4.213], color="blue red", name="Arrows_14.0630002022_1")

cluster_dict["14.0630002022"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(34.0), float(-7.0), float(2.5), float(1.0)]

cluster_dict["14.0630002022_arrows"] += cgo_arrow([34.0,-7.0,2.5], [31.839,-6.19,4.213], color="blue red", name="Arrows_14.0630002022_2")

cluster_dict["14.0630002022"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(33.5), float(-3.5), float(2.5), float(1.0)]

cluster_dict["14.0630002022_arrows"] += cgo_arrow([33.5,-3.5,2.5], [31.597,-4.137,0.858], color="blue red", name="Arrows_14.0630002022_3")

cluster_dict["14.0630002022"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(-4.0), float(2.0), float(1.0)]

cluster_dict["14.0630002022_arrows"] += cgo_arrow([36.5,-4.0,2.0], [37.036,-1.381,0.495], color="blue red", name="Arrows_14.0630002022_4")

cluster_dict["14.0630002022"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(39.5), float(-4.5), float(3.0), float(1.0)]

cluster_dict["14.0630002022_arrows"] += cgo_arrow([39.5,-4.5,3.0], [38.748,-6.913,1.552], color="blue red", name="Arrows_14.0630002022_5")

cluster_dict["14.0630002022"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(36.6766121811), float(-3.51100683087), float(5.26062508525), float(1.0)]


cluster_dict["14.0630002022"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(35.0), float(-4.5), float(0.0), float(1.0)]

cluster_dict["14.0630002022_arrows"] += cgo_arrow([35.0,-4.5,0.0], [36.469,-5.476,-2.507], color="red blue", name="Arrows_14.0630002022_6")

cluster_dict["14.0630002022"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(34.0), float(-6.5), float(2.0), float(1.0)]

cluster_dict["14.0630002022_arrows"] += cgo_arrow([34.0,-6.5,2.0], [31.839,-6.19,4.213], color="red blue", name="Arrows_14.0630002022_7")

cluster_dict["14.0630002022"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(34.0), float(-6.5), float(2.0), float(1.0)]

cluster_dict["14.0630002022_arrows"] += cgo_arrow([34.0,-6.5,2.0], [31.839,-6.19,4.213], color="red blue", name="Arrows_14.0630002022_8")

cluster_dict["14.0630002022"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(37.0), float(-6.0), float(8.0), float(1.0)]

cluster_dict["14.0630002022_arrows"] += cgo_arrow([37.0,-6.0,8.0], [39.742,-6.281,9.454], color="red blue", name="Arrows_14.0630002022_9")

cluster_dict["14.0630002022"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(40.5), float(-2.5), float(4.0), float(1.0)]

cluster_dict["14.0630002022_arrows"] += cgo_arrow([40.5,-2.5,4.0], [43.985,-5.458,4.808], color="red blue", name="Arrows_14.0630002022_10")

cmd.load_cgo(cluster_dict["14.0630002022"], "Features_14.0630002022", 1)
cmd.load_cgo(cluster_dict["14.0630002022_arrows"], "Arrows_14.0630002022")
cmd.set("transparency", 0.2,"Features_14.0630002022")
cmd.group("Pharmacophore_14.0630002022", members="Features_14.0630002022")
cmd.group("Pharmacophore_14.0630002022", members="Arrows_14.0630002022")

if dirpath:
    f = join(dirpath, "0/label_threshold_14.0630002022.mol2")
else:
    f = "0/label_threshold_14.0630002022.mol2"

cmd.load(f, 'label_threshold_14.0630002022')
cmd.hide('everything', 'label_threshold_14.0630002022')
cmd.label("label_threshold_14.0630002022", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.0630002022', members= 'label_threshold_14.0630002022')


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
    f = join(dirpath, "1/label_threshold_10.6.mol2")
else:
    f = "1/label_threshold_10.6.mol2"

cmd.load(f, 'label_threshold_10.6')
cmd.hide('everything', 'label_threshold_10.6')
cmd.label("label_threshold_10.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.6]
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


cluster_dict = {"14.0080003738":[], "14.0080003738_arrows":[]}

cluster_dict["14.0080003738"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(32.5), float(-6.5), float(6.5), float(1.0)]

cluster_dict["14.0080003738_arrows"] += cgo_arrow([32.5,-6.5,6.5], [31.839,-6.19,4.213], color="blue red", name="Arrows_14.0080003738_1")

cluster_dict["14.0080003738"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(34.0), float(-7.0), float(2.5), float(1.0)]

cluster_dict["14.0080003738_arrows"] += cgo_arrow([34.0,-7.0,2.5], [31.839,-6.19,4.213], color="blue red", name="Arrows_14.0080003738_2")

cluster_dict["14.0080003738"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(34.0), float(-3.5), float(2.5), float(1.0)]

cluster_dict["14.0080003738_arrows"] += cgo_arrow([34.0,-3.5,2.5], [31.597,-4.137,0.858], color="blue red", name="Arrows_14.0080003738_3")

cluster_dict["14.0080003738"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(-4.0), float(2.0), float(1.0)]

cluster_dict["14.0080003738_arrows"] += cgo_arrow([36.5,-4.0,2.0], [37.036,-1.381,0.495], color="blue red", name="Arrows_14.0080003738_4")

cluster_dict["14.0080003738"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(39.5), float(-4.5), float(3.0), float(1.0)]

cluster_dict["14.0080003738_arrows"] += cgo_arrow([39.5,-4.5,3.0], [38.748,-6.913,1.552], color="blue red", name="Arrows_14.0080003738_5")

cluster_dict["14.0080003738"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(36.795633817), float(-3.43865303097), float(5.49714043308), float(1.0)]


cluster_dict["14.0080003738"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(39.962087635), float(-4.27949997829), float(13.8465463125), float(1.0)]


cluster_dict["14.0080003738"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(33.5), float(-6.5), float(2.0), float(1.0)]

cluster_dict["14.0080003738_arrows"] += cgo_arrow([33.5,-6.5,2.0], [31.839,-6.19,4.213], color="red blue", name="Arrows_14.0080003738_6")

cluster_dict["14.0080003738"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(37.0), float(-6.0), float(8.0), float(1.0)]

cluster_dict["14.0080003738_arrows"] += cgo_arrow([37.0,-6.0,8.0], [39.742,-6.281,9.454], color="red blue", name="Arrows_14.0080003738_7")

cluster_dict["14.0080003738"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(40.5), float(-2.5), float(4.0), float(1.0)]

cluster_dict["14.0080003738_arrows"] += cgo_arrow([40.5,-2.5,4.0], [43.985,-5.458,4.808], color="red blue", name="Arrows_14.0080003738_8")

cmd.load_cgo(cluster_dict["14.0080003738"], "Features_14.0080003738", 1)
cmd.load_cgo(cluster_dict["14.0080003738_arrows"], "Arrows_14.0080003738")
cmd.set("transparency", 0.2,"Features_14.0080003738")
cmd.group("Pharmacophore_14.0080003738", members="Features_14.0080003738")
cmd.group("Pharmacophore_14.0080003738", members="Arrows_14.0080003738")

if dirpath:
    f = join(dirpath, "1/label_threshold_14.0080003738.mol2")
else:
    f = "1/label_threshold_14.0080003738.mol2"

cmd.load(f, 'label_threshold_14.0080003738')
cmd.hide('everything', 'label_threshold_14.0080003738')
cmd.label("label_threshold_14.0080003738", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.0080003738', members= 'label_threshold_14.0080003738')


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
    f = join(dirpath, "2/label_threshold_0.6.mol2")
else:
    f = "2/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
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


cluster_dict = {"11.5279998779":[], "11.5279998779_arrows":[]}

cluster_dict["11.5279998779"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(-10.0), float(-7.5), float(1.0)]

cluster_dict["11.5279998779_arrows"] += cgo_arrow([23.0,-10.0,-7.5], [23.21,-9.109,-10.141], color="blue red", name="Arrows_11.5279998779_1")

cluster_dict["11.5279998779"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.9532888127), float(-9.11713701356), float(-4.40654102579), float(1.0)]


cluster_dict["11.5279998779"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(-8.5), float(-3.0), float(1.0)]


cluster_dict["11.5279998779"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(-11.5), float(-7.0), float(1.0)]

cluster_dict["11.5279998779_arrows"] += cgo_arrow([23.5,-11.5,-7.0], [26.003,-10.371,-8.463], color="red blue", name="Arrows_11.5279998779_2")

cluster_dict["11.5279998779"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(-8.5), float(-3.0), float(1.0)]

cluster_dict["11.5279998779_arrows"] += cgo_arrow([24.5,-8.5,-3.0], [26.804,-8.147,-6.458], color="red blue", name="Arrows_11.5279998779_3")

cmd.load_cgo(cluster_dict["11.5279998779"], "Features_11.5279998779", 1)
cmd.load_cgo(cluster_dict["11.5279998779_arrows"], "Arrows_11.5279998779")
cmd.set("transparency", 0.2,"Features_11.5279998779")
cmd.group("Pharmacophore_11.5279998779", members="Features_11.5279998779")
cmd.group("Pharmacophore_11.5279998779", members="Arrows_11.5279998779")

if dirpath:
    f = join(dirpath, "2/label_threshold_11.5279998779.mol2")
else:
    f = "2/label_threshold_11.5279998779.mol2"

cmd.load(f, 'label_threshold_11.5279998779')
cmd.hide('everything', 'label_threshold_11.5279998779')
cmd.label("label_threshold_11.5279998779", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.5279998779', members= 'label_threshold_11.5279998779')


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
    f = join(dirpath, "3/label_threshold_0.6.mol2")
else:
    f = "3/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
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


cluster_dict = {"10.9890003204":[], "10.9890003204_arrows":[]}

cluster_dict["10.9890003204"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(40.0), float(-30.5), float(-2.5), float(1.0)]

cluster_dict["10.9890003204_arrows"] += cgo_arrow([40.0,-30.5,-2.5], [42.411,-31.449,-1.202], color="blue red", name="Arrows_10.9890003204_1")

cluster_dict["10.9890003204"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(40.0), float(-30.5), float(-2.5), float(1.0)]

cluster_dict["10.9890003204_arrows"] += cgo_arrow([40.0,-30.5,-2.5], [42.411,-31.449,-1.202], color="blue red", name="Arrows_10.9890003204_2")

cluster_dict["10.9890003204"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(42.5), float(-29.0), float(-6.5), float(1.0)]

cluster_dict["10.9890003204_arrows"] += cgo_arrow([42.5,-29.0,-6.5], [42.81,-25.521,-7.429], color="blue red", name="Arrows_10.9890003204_3")

cluster_dict["10.9890003204"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(39.7375000458), float(-30.5855743051), float(-5.76213425097), float(1.0)]


cluster_dict["10.9890003204"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(40.3548501636), float(-30.2674052996), float(-0.338699864971), float(1.0)]


cmd.load_cgo(cluster_dict["10.9890003204"], "Features_10.9890003204", 1)
cmd.load_cgo(cluster_dict["10.9890003204_arrows"], "Arrows_10.9890003204")
cmd.set("transparency", 0.2,"Features_10.9890003204")
cmd.group("Pharmacophore_10.9890003204", members="Features_10.9890003204")
cmd.group("Pharmacophore_10.9890003204", members="Arrows_10.9890003204")

if dirpath:
    f = join(dirpath, "3/label_threshold_10.9890003204.mol2")
else:
    f = "3/label_threshold_10.9890003204.mol2"

cmd.load(f, 'label_threshold_10.9890003204')
cmd.hide('everything', 'label_threshold_10.9890003204')
cmd.label("label_threshold_10.9890003204", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.9890003204', members= 'label_threshold_10.9890003204')


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


cluster_dict = {"10.4090003967":[], "10.4090003967_arrows":[]}

cluster_dict["10.4090003967"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(19.0153747938), float(-22.9368185628), float(7.38310950564), float(1.0)]


cmd.load_cgo(cluster_dict["10.4090003967"], "Features_10.4090003967", 1)
cmd.load_cgo(cluster_dict["10.4090003967_arrows"], "Arrows_10.4090003967")
cmd.set("transparency", 0.2,"Features_10.4090003967")
cmd.group("Pharmacophore_10.4090003967", members="Features_10.4090003967")
cmd.group("Pharmacophore_10.4090003967", members="Arrows_10.4090003967")

if dirpath:
    f = join(dirpath, "4/label_threshold_10.4090003967.mol2")
else:
    f = "4/label_threshold_10.4090003967.mol2"

cmd.load(f, 'label_threshold_10.4090003967')
cmd.hide('everything', 'label_threshold_10.4090003967')
cmd.label("label_threshold_10.4090003967", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.4090003967', members= 'label_threshold_10.4090003967')


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
    f = join(dirpath, "5/label_threshold_1.8.mol2")
else:
    f = "5/label_threshold_1.8.mol2"

cmd.load(f, 'label_threshold_1.8')
cmd.hide('everything', 'label_threshold_1.8')
cmd.label("label_threshold_1.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [1.8]
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


cluster_dict = {"10.0380001068":[], "10.0380001068_arrows":[]}

cluster_dict["10.0380001068"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(33.0), float(-6.5), float(6.0), float(1.0)]

cluster_dict["10.0380001068_arrows"] += cgo_arrow([33.0,-6.5,6.0], [31.839,-6.19,4.213], color="blue red", name="Arrows_10.0380001068_1")

cluster_dict["10.0380001068"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(34.0), float(-7.0), float(2.5), float(1.0)]

cluster_dict["10.0380001068_arrows"] += cgo_arrow([34.0,-7.0,2.5], [31.839,-6.19,4.213], color="blue red", name="Arrows_10.0380001068_2")

cluster_dict["10.0380001068"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(33.5), float(-3.5), float(2.5), float(1.0)]

cluster_dict["10.0380001068_arrows"] += cgo_arrow([33.5,-3.5,2.5], [31.597,-4.137,0.858], color="blue red", name="Arrows_10.0380001068_3")

cluster_dict["10.0380001068"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(-4.0), float(2.0), float(1.0)]

cluster_dict["10.0380001068_arrows"] += cgo_arrow([36.5,-4.0,2.0], [37.036,-1.381,0.495], color="blue red", name="Arrows_10.0380001068_4")

cluster_dict["10.0380001068"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(39.5), float(-4.5), float(3.0), float(1.0)]

cluster_dict["10.0380001068_arrows"] += cgo_arrow([39.5,-4.5,3.0], [38.748,-6.913,1.552], color="blue red", name="Arrows_10.0380001068_5")

cluster_dict["10.0380001068"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(34.4830329577), float(-11.5681629118), float(-4.93749746018), float(1.0)]


cluster_dict["10.0380001068"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(35.5998584569), float(-5.2805775259), float(2.75216303626), float(1.0)]


cluster_dict["10.0380001068"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(35.3222405683), float(-8.46711037727), float(-1.5), float(1.0)]


cluster_dict["10.0380001068"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(-10.7111234635), float(-1.5), float(1.0)]


cluster_dict["10.0380001068"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(38.0), float(-12.4566357403), float(-4.41327148064), float(1.0)]


cluster_dict["10.0380001068"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(33.0), float(-12.5), float(-3.0), float(1.0)]

cluster_dict["10.0380001068_arrows"] += cgo_arrow([33.0,-12.5,-3.0], [33.491,-9.923,-1.371], color="red blue", name="Arrows_10.0380001068_6")

cluster_dict["10.0380001068"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(32.5), float(-11.0), float(-3.5), float(1.0)]

cluster_dict["10.0380001068_arrows"] += cgo_arrow([32.5,-11.0,-3.5], [33.491,-9.923,-1.371], color="red blue", name="Arrows_10.0380001068_7")

cluster_dict["10.0380001068"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(35.0), float(-4.5), float(0.0), float(1.0)]

cluster_dict["10.0380001068_arrows"] += cgo_arrow([35.0,-4.5,0.0], [36.469,-5.476,-2.507], color="red blue", name="Arrows_10.0380001068_8")

cluster_dict["10.0380001068"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(34.0), float(-6.5), float(2.0), float(1.0)]

cluster_dict["10.0380001068_arrows"] += cgo_arrow([34.0,-6.5,2.0], [31.839,-6.19,4.213], color="red blue", name="Arrows_10.0380001068_9")

cluster_dict["10.0380001068"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(34.0), float(-6.5), float(2.0), float(1.0)]

cluster_dict["10.0380001068_arrows"] += cgo_arrow([34.0,-6.5,2.0], [31.839,-6.19,4.213], color="red blue", name="Arrows_10.0380001068_10")

cluster_dict["10.0380001068"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(35.5), float(-14.0), float(-4.0), float(1.0)]

cluster_dict["10.0380001068_arrows"] += cgo_arrow([35.5,-14.0,-4.0], [35.409,-14.894,-1.405], color="red blue", name="Arrows_10.0380001068_11")

cluster_dict["10.0380001068"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(35.5), float(-12.0), float(-2.5), float(1.0)]

cluster_dict["10.0380001068_arrows"] += cgo_arrow([35.5,-12.0,-2.5], [34.64,-12.416,-0.25], color="red blue", name="Arrows_10.0380001068_12")

cluster_dict["10.0380001068"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(36.0), float(-10.5), float(-3.0), float(1.0)]

cluster_dict["10.0380001068_arrows"] += cgo_arrow([36.0,-10.5,-3.0], [33.491,-9.923,-1.371], color="red blue", name="Arrows_10.0380001068_13")

cmd.load_cgo(cluster_dict["10.0380001068"], "Features_10.0380001068", 1)
cmd.load_cgo(cluster_dict["10.0380001068_arrows"], "Arrows_10.0380001068")
cmd.set("transparency", 0.2,"Features_10.0380001068")
cmd.group("Pharmacophore_10.0380001068", members="Features_10.0380001068")
cmd.group("Pharmacophore_10.0380001068", members="Arrows_10.0380001068")

if dirpath:
    f = join(dirpath, "5/label_threshold_10.0380001068.mol2")
else:
    f = "5/label_threshold_10.0380001068.mol2"

cmd.load(f, 'label_threshold_10.0380001068')
cmd.hide('everything', 'label_threshold_10.0380001068')
cmd.label("label_threshold_10.0380001068", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.0380001068', members= 'label_threshold_10.0380001068')


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
