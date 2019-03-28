
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
    f = join(dirpath, "0/label_threshold_13.6.mol2")
else:
    f = "0/label_threshold_13.6.mol2"

cmd.load(f, 'label_threshold_13.6')
cmd.hide('everything', 'label_threshold_13.6')
cmd.label("label_threshold_13.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.6]
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


cluster_dict = {"16.5279998779":[], "16.5279998779_arrows":[]}

cluster_dict["16.5279998779"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(15.5), float(-4.5), float(19.5), float(1.0)]

cluster_dict["16.5279998779_arrows"] += cgo_arrow([15.5,-4.5,19.5], [12.907,-4.46,18.205], color="blue red", name="Arrows_16.5279998779_1")

cluster_dict["16.5279998779"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(-7.5), float(31.0), float(1.0)]

cluster_dict["16.5279998779_arrows"] += cgo_arrow([24.5,-7.5,31.0], [21.289,-7.22,30.788], color="blue red", name="Arrows_16.5279998779_2")

cluster_dict["16.5279998779"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(17.2820920599), float(-7.41809412719), float(17.8584068059), float(1.0)]


cluster_dict["16.5279998779"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(14.6612025271), float(-2.60338962482), float(22.1612025271), float(1.0)]


cluster_dict["16.5279998779"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.6603494549), float(-8.1243797371), float(29.8063053593), float(1.0)]


cluster_dict["16.5279998779"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(-8.5), float(13.5), float(1.0)]

cluster_dict["16.5279998779_arrows"] += cgo_arrow([18.5,-8.5,13.5], [20.712,-4.579,11.589], color="red blue", name="Arrows_16.5279998779_3")

cluster_dict["16.5279998779"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(-3.5), float(23.0), float(1.0)]

cluster_dict["16.5279998779_arrows"] += cgo_arrow([15.0,-3.5,23.0], [13.814,-5.448,21.846], color="red blue", name="Arrows_16.5279998779_4")

cluster_dict["16.5279998779"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(-2.5), float(18.0), float(1.0)]

cluster_dict["16.5279998779_arrows"] += cgo_arrow([16.0,-2.5,18.0], [15.247,-0.45,15.978], color="red blue", name="Arrows_16.5279998779_5")

cluster_dict["16.5279998779"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(-7.5), float(17.0), float(1.0)]

cluster_dict["16.5279998779_arrows"] += cgo_arrow([18.0,-7.5,17.0], [20.681,-4.461,18.07], color="red blue", name="Arrows_16.5279998779_6")

cluster_dict["16.5279998779"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(-8.5), float(13.5), float(1.0)]

cluster_dict["16.5279998779_arrows"] += cgo_arrow([18.5,-8.5,13.5], [20.712,-4.579,11.589], color="red blue", name="Arrows_16.5279998779_7")

cluster_dict["16.5279998779"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(-10.0), float(28.0), float(1.0)]

cluster_dict["16.5279998779_arrows"] += cgo_arrow([23.0,-10.0,28.0], [22.302,-11.378,26.314], color="red blue", name="Arrows_16.5279998779_8")

cluster_dict["16.5279998779"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(-6.5), float(32.5), float(1.0)]

cluster_dict["16.5279998779_arrows"] += cgo_arrow([23.5,-6.5,32.5], [22.178,-4.461,34.632], color="red blue", name="Arrows_16.5279998779_9")

cmd.load_cgo(cluster_dict["16.5279998779"], "Features_16.5279998779", 1)
cmd.load_cgo(cluster_dict["16.5279998779_arrows"], "Arrows_16.5279998779")
cmd.set("transparency", 0.2,"Features_16.5279998779")
cmd.group("Pharmacophore_16.5279998779", members="Features_16.5279998779")
cmd.group("Pharmacophore_16.5279998779", members="Arrows_16.5279998779")

if dirpath:
    f = join(dirpath, "0/label_threshold_16.5279998779.mol2")
else:
    f = "0/label_threshold_16.5279998779.mol2"

cmd.load(f, 'label_threshold_16.5279998779')
cmd.hide('everything', 'label_threshold_16.5279998779')
cmd.label("label_threshold_16.5279998779", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.5279998779', members= 'label_threshold_16.5279998779')


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
    f = join(dirpath, "1/label_threshold_12.6.mol2")
else:
    f = "1/label_threshold_12.6.mol2"

cmd.load(f, 'label_threshold_12.6')
cmd.hide('everything', 'label_threshold_12.6')
cmd.label("label_threshold_12.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.6]
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


cluster_dict = {"16.1299991608":[], "16.1299991608_arrows":[]}

cluster_dict["16.1299991608"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(-3.0), float(23.5), float(1.0)]

cluster_dict["16.1299991608_arrows"] += cgo_arrow([15.0,-3.0,23.5], [12.932,-3.265,24.968], color="blue red", name="Arrows_16.1299991608_1")

cluster_dict["16.1299991608"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(15.5), float(-4.5), float(19.5), float(1.0)]

cluster_dict["16.1299991608_arrows"] += cgo_arrow([15.5,-4.5,19.5], [12.907,-4.46,18.205], color="blue red", name="Arrows_16.1299991608_2")

cluster_dict["16.1299991608"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(17.3342556183), float(-7.50462159966), float(17.7951420183), float(1.0)]


cluster_dict["16.1299991608"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(-0.5), float(19.5), float(1.0)]


cluster_dict["16.1299991608"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(23.0784932487), float(-7.63340659948), float(27.5), float(1.0)]


cluster_dict["16.1299991608"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(-8.5), float(13.5), float(1.0)]

cluster_dict["16.1299991608_arrows"] += cgo_arrow([18.5,-8.5,13.5], [20.712,-4.579,11.589], color="red blue", name="Arrows_16.1299991608_3")

cluster_dict["16.1299991608"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(-3.5), float(23.0), float(1.0)]

cluster_dict["16.1299991608_arrows"] += cgo_arrow([15.0,-3.5,23.0], [13.814,-5.448,21.846], color="red blue", name="Arrows_16.1299991608_4")

cluster_dict["16.1299991608"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(-2.5), float(18.0), float(1.0)]

cluster_dict["16.1299991608_arrows"] += cgo_arrow([16.0,-2.5,18.0], [15.247,-0.45,15.978], color="red blue", name="Arrows_16.1299991608_5")

cluster_dict["16.1299991608"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(-7.5), float(17.0), float(1.0)]

cluster_dict["16.1299991608_arrows"] += cgo_arrow([18.0,-7.5,17.0], [20.681,-4.461,18.07], color="red blue", name="Arrows_16.1299991608_6")

cluster_dict["16.1299991608"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(-8.5), float(13.5), float(1.0)]

cluster_dict["16.1299991608_arrows"] += cgo_arrow([18.5,-8.5,13.5], [20.712,-4.579,11.589], color="red blue", name="Arrows_16.1299991608_7")

cluster_dict["16.1299991608"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(-8.0), float(12.0), float(1.0)]

cluster_dict["16.1299991608_arrows"] += cgo_arrow([21.0,-8.0,12.0], [21.313,-6.25,9.361], color="red blue", name="Arrows_16.1299991608_8")

cluster_dict["16.1299991608"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(-10.0), float(27.5), float(1.0)]

cluster_dict["16.1299991608_arrows"] += cgo_arrow([22.5,-10.0,27.5], [22.302,-11.378,26.314], color="red blue", name="Arrows_16.1299991608_9")

cmd.load_cgo(cluster_dict["16.1299991608"], "Features_16.1299991608", 1)
cmd.load_cgo(cluster_dict["16.1299991608_arrows"], "Arrows_16.1299991608")
cmd.set("transparency", 0.2,"Features_16.1299991608")
cmd.group("Pharmacophore_16.1299991608", members="Features_16.1299991608")
cmd.group("Pharmacophore_16.1299991608", members="Arrows_16.1299991608")

if dirpath:
    f = join(dirpath, "1/label_threshold_16.1299991608.mol2")
else:
    f = "1/label_threshold_16.1299991608.mol2"

cmd.load(f, 'label_threshold_16.1299991608')
cmd.hide('everything', 'label_threshold_16.1299991608')
cmd.label("label_threshold_16.1299991608", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.1299991608', members= 'label_threshold_16.1299991608')


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
    f = join(dirpath, "2/label_threshold_7.1.mol2")
else:
    f = "2/label_threshold_7.1.mol2"

cmd.load(f, 'label_threshold_7.1')
cmd.hide('everything', 'label_threshold_7.1')
cmd.label("label_threshold_7.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [7.1]
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


cluster_dict = {"14.8190002441":[], "14.8190002441_arrows":[]}

cluster_dict["14.8190002441"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(-7.0), float(28.0), float(1.0)]

cluster_dict["14.8190002441_arrows"] += cgo_arrow([23.0,-7.0,28.0], [21.289,-7.22,30.788], color="blue red", name="Arrows_14.8190002441_1")

cluster_dict["14.8190002441"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(-7.0), float(28.0), float(1.0)]

cluster_dict["14.8190002441_arrows"] += cgo_arrow([23.0,-7.0,28.0], [21.289,-7.22,30.788], color="blue red", name="Arrows_14.8190002441_2")

cluster_dict["14.8190002441"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(-5.5), float(34.5), float(1.0)]

cluster_dict["14.8190002441_arrows"] += cgo_arrow([26.0,-5.5,34.5], [23.841,-3.603,35.944], color="blue red", name="Arrows_14.8190002441_3")

cluster_dict["14.8190002441"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(26.9964771213), float(-7.68033226517), float(32.2168000869), float(1.0)]


cluster_dict["14.8190002441"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(-10.0), float(28.0), float(1.0)]

cluster_dict["14.8190002441_arrows"] += cgo_arrow([23.0,-10.0,28.0], [22.302,-11.378,26.314], color="red blue", name="Arrows_14.8190002441_4")

cluster_dict["14.8190002441"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(-6.5), float(33.0), float(1.0)]

cluster_dict["14.8190002441_arrows"] += cgo_arrow([24.0,-6.5,33.0], [22.178,-4.461,34.632], color="red blue", name="Arrows_14.8190002441_5")

cmd.load_cgo(cluster_dict["14.8190002441"], "Features_14.8190002441", 1)
cmd.load_cgo(cluster_dict["14.8190002441_arrows"], "Arrows_14.8190002441")
cmd.set("transparency", 0.2,"Features_14.8190002441")
cmd.group("Pharmacophore_14.8190002441", members="Features_14.8190002441")
cmd.group("Pharmacophore_14.8190002441", members="Arrows_14.8190002441")

if dirpath:
    f = join(dirpath, "2/label_threshold_14.8190002441.mol2")
else:
    f = "2/label_threshold_14.8190002441.mol2"

cmd.load(f, 'label_threshold_14.8190002441')
cmd.hide('everything', 'label_threshold_14.8190002441')
cmd.label("label_threshold_14.8190002441", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.8190002441', members= 'label_threshold_14.8190002441')


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
    f = join(dirpath, "3/label_threshold_8.3.mol2")
else:
    f = "3/label_threshold_8.3.mol2"

cmd.load(f, 'label_threshold_8.3')
cmd.hide('everything', 'label_threshold_8.3')
cmd.label("label_threshold_8.3", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [8.3]
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


cluster_dict = {"12.9890003204":[], "12.9890003204_arrows":[]}

cluster_dict["12.9890003204"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.5), float(-5.0), float(22.5), float(1.0)]

cluster_dict["12.9890003204_arrows"] += cgo_arrow([19.5,-5.0,22.5], [21.795,-3.156,22.131], color="blue red", name="Arrows_12.9890003204_1")

cluster_dict["12.9890003204"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(19.8739521316), float(-8.21670902118), float(15.0221315962), float(1.0)]


cluster_dict["12.9890003204"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(-8.0), float(13.5), float(1.0)]

cluster_dict["12.9890003204_arrows"] += cgo_arrow([18.5,-8.0,13.5], [20.712,-4.579,11.589], color="red blue", name="Arrows_12.9890003204_2")

cluster_dict["12.9890003204"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.5), float(-3.0), float(17.5), float(1.0)]

cluster_dict["12.9890003204_arrows"] += cgo_arrow([16.5,-3.0,17.5], [15.247,-0.45,15.978], color="red blue", name="Arrows_12.9890003204_3")

cluster_dict["12.9890003204"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(-7.5), float(17.0), float(1.0)]

cluster_dict["12.9890003204_arrows"] += cgo_arrow([18.0,-7.5,17.0], [20.681,-4.461,18.07], color="red blue", name="Arrows_12.9890003204_4")

cluster_dict["12.9890003204"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(-8.0), float(13.5), float(1.0)]

cluster_dict["12.9890003204_arrows"] += cgo_arrow([18.5,-8.0,13.5], [20.712,-4.579,11.589], color="red blue", name="Arrows_12.9890003204_5")

cluster_dict["12.9890003204"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(-6.0), float(13.5), float(1.0)]

cluster_dict["12.9890003204_arrows"] += cgo_arrow([21.0,-6.0,13.5], [20.712,-4.579,11.589], color="red blue", name="Arrows_12.9890003204_6")

cluster_dict["12.9890003204"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(-6.0), float(13.5), float(1.0)]

cluster_dict["12.9890003204_arrows"] += cgo_arrow([21.0,-6.0,13.5], [20.712,-4.579,11.589], color="red blue", name="Arrows_12.9890003204_7")

cmd.load_cgo(cluster_dict["12.9890003204"], "Features_12.9890003204", 1)
cmd.load_cgo(cluster_dict["12.9890003204_arrows"], "Arrows_12.9890003204")
cmd.set("transparency", 0.2,"Features_12.9890003204")
cmd.group("Pharmacophore_12.9890003204", members="Features_12.9890003204")
cmd.group("Pharmacophore_12.9890003204", members="Arrows_12.9890003204")

if dirpath:
    f = join(dirpath, "3/label_threshold_12.9890003204.mol2")
else:
    f = "3/label_threshold_12.9890003204.mol2"

cmd.load(f, 'label_threshold_12.9890003204')
cmd.hide('everything', 'label_threshold_12.9890003204')
cmd.label("label_threshold_12.9890003204", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.9890003204', members= 'label_threshold_12.9890003204')


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


cluster_dict = {"12.1490001678":[], "12.1490001678_arrows":[]}

cluster_dict["12.1490001678"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-2.0), float(9.0), float(24.0), float(1.0)]

cluster_dict["12.1490001678_arrows"] += cgo_arrow([-2.0,9.0,24.0], [-0.803,7.183,25.687], color="blue red", name="Arrows_12.1490001678_1")

cluster_dict["12.1490001678"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(2.5), float(11.0), float(22.0), float(1.0)]

cluster_dict["12.1490001678_arrows"] += cgo_arrow([2.5,11.0,22.0], [3.385,11.103,19.6], color="blue red", name="Arrows_12.1490001678_2")

cluster_dict["12.1490001678"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-0.39239359998), float(11.6198780014), float(22.9163385745), float(1.0)]


cluster_dict["12.1490001678"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(1.5), float(9.0), float(22.0), float(1.0)]

cluster_dict["12.1490001678_arrows"] += cgo_arrow([1.5,9.0,22.0], [2.304,6.061,22.03], color="red blue", name="Arrows_12.1490001678_3")

cmd.load_cgo(cluster_dict["12.1490001678"], "Features_12.1490001678", 1)
cmd.load_cgo(cluster_dict["12.1490001678_arrows"], "Arrows_12.1490001678")
cmd.set("transparency", 0.2,"Features_12.1490001678")
cmd.group("Pharmacophore_12.1490001678", members="Features_12.1490001678")
cmd.group("Pharmacophore_12.1490001678", members="Arrows_12.1490001678")

if dirpath:
    f = join(dirpath, "4/label_threshold_12.1490001678.mol2")
else:
    f = "4/label_threshold_12.1490001678.mol2"

cmd.load(f, 'label_threshold_12.1490001678')
cmd.hide('everything', 'label_threshold_12.1490001678')
cmd.label("label_threshold_12.1490001678", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.1490001678', members= 'label_threshold_12.1490001678')


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
    f = join(dirpath, "5/label_threshold_4.0.mol2")
else:
    f = "5/label_threshold_4.0.mol2"

cmd.load(f, 'label_threshold_4.0')
cmd.hide('everything', 'label_threshold_4.0')
cmd.label("label_threshold_4.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [4.0]
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


cluster_dict = {"10.7869997025":[], "10.7869997025_arrows":[]}

cluster_dict["10.7869997025"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.0), float(10.5), float(15.0), float(1.0)]

cluster_dict["10.7869997025_arrows"] += cgo_arrow([29.0,10.5,15.0], [29.506,7.736,14.929], color="blue red", name="Arrows_10.7869997025_1")

cluster_dict["10.7869997025"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.1758019789), float(9.17696502529), float(13.4116421202), float(1.0)]


cmd.load_cgo(cluster_dict["10.7869997025"], "Features_10.7869997025", 1)
cmd.load_cgo(cluster_dict["10.7869997025_arrows"], "Arrows_10.7869997025")
cmd.set("transparency", 0.2,"Features_10.7869997025")
cmd.group("Pharmacophore_10.7869997025", members="Features_10.7869997025")
cmd.group("Pharmacophore_10.7869997025", members="Arrows_10.7869997025")

if dirpath:
    f = join(dirpath, "5/label_threshold_10.7869997025.mol2")
else:
    f = "5/label_threshold_10.7869997025.mol2"

cmd.load(f, 'label_threshold_10.7869997025')
cmd.hide('everything', 'label_threshold_10.7869997025')
cmd.label("label_threshold_10.7869997025", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.7869997025', members= 'label_threshold_10.7869997025')


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
