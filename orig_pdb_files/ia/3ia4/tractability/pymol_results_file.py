
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
    f = join(dirpath, "0/label_threshold_0.7.mol2")
else:
    f = "0/label_threshold_0.7.mol2"

cmd.load(f, 'label_threshold_0.7')
cmd.hide('everything', 'label_threshold_0.7')
cmd.label("label_threshold_0.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.7]
gfiles = ['0/apolar.grd']
grids = ['apolar']
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


cluster_dict = {"23.8479995728":[], "23.8479995728_arrows":[]}

cluster_dict["23.8479995728"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(46.3934331569), float(-15.2053076247), float(-27.3129684872), float(1.0)]


cmd.load_cgo(cluster_dict["23.8479995728"], "Features_23.8479995728", 1)
cmd.load_cgo(cluster_dict["23.8479995728_arrows"], "Arrows_23.8479995728")
cmd.set("transparency", 0.2,"Features_23.8479995728")
cmd.group("Pharmacophore_23.8479995728", members="Features_23.8479995728")
cmd.group("Pharmacophore_23.8479995728", members="Arrows_23.8479995728")

if dirpath:
    f = join(dirpath, "0/label_threshold_23.8479995728.mol2")
else:
    f = "0/label_threshold_23.8479995728.mol2"

cmd.load(f, 'label_threshold_23.8479995728')
cmd.hide('everything', 'label_threshold_23.8479995728')
cmd.label("label_threshold_23.8479995728", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_23.8479995728', members= 'label_threshold_23.8479995728')


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
    f = join(dirpath, "1/label_threshold_16.0.mol2")
else:
    f = "1/label_threshold_16.0.mol2"

cmd.load(f, 'label_threshold_16.0')
cmd.hide('everything', 'label_threshold_16.0')
cmd.label("label_threshold_16.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [16.0]
gfiles = ['1/apolar.grd']
grids = ['apolar']
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


cluster_dict = {"17.0990009308":[], "17.0990009308_arrows":[]}

cluster_dict["17.0990009308"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(38.5), float(-2.0), float(-19.5), float(1.0)]

cluster_dict["17.0990009308_arrows"] += cgo_arrow([38.5,-2.0,-19.5], [38.057,-3.157,-21.901], color="blue red", name="Arrows_17.0990009308_1")

cluster_dict["17.0990009308"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(39.5), float(-0.5), float(-17.0), float(1.0)]

cluster_dict["17.0990009308_arrows"] += cgo_arrow([39.5,-0.5,-17.0], [39.544,-0.566,-14.152], color="blue red", name="Arrows_17.0990009308_2")

cluster_dict["17.0990009308"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(42.0), float(-0.5), float(-21.0), float(1.0)]

cluster_dict["17.0990009308_arrows"] += cgo_arrow([42.0,-0.5,-21.0], [42.983,-2.845,-21.818], color="blue red", name="Arrows_17.0990009308_3")

cluster_dict["17.0990009308"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(41.8267724498), float(3.5983583096), float(-21.4543893521), float(1.0)]


cluster_dict["17.0990009308"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(36.0), float(-3.0), float(-17.5), float(1.0)]

cluster_dict["17.0990009308_arrows"] += cgo_arrow([36.0,-3.0,-17.5], [36.188,-3.85,-14.011], color="red blue", name="Arrows_17.0990009308_4")

cluster_dict["17.0990009308"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(36.0), float(-3.0), float(-17.5), float(1.0)]

cluster_dict["17.0990009308_arrows"] += cgo_arrow([36.0,-3.0,-17.5], [36.188,-3.85,-14.011], color="red blue", name="Arrows_17.0990009308_5")

cluster_dict["17.0990009308"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(38.0), float(5.0), float(-18.0), float(1.0)]

cluster_dict["17.0990009308_arrows"] += cgo_arrow([38.0,5.0,-18.0], [36.971,3.647,-15.506], color="red blue", name="Arrows_17.0990009308_6")

cluster_dict["17.0990009308"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(42.0), float(-2.5), float(-17.5), float(1.0)]

cluster_dict["17.0990009308_arrows"] += cgo_arrow([42.0,-2.5,-17.5], [43.999,-4.142,-19.151], color="red blue", name="Arrows_17.0990009308_7")

cluster_dict["17.0990009308"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(42.0), float(0.5), float(-22.0), float(1.0)]


cluster_dict["17.0990009308"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(46.5), float(6.0), float(-25.5), float(1.0)]

cluster_dict["17.0990009308_arrows"] += cgo_arrow([46.5,6.0,-25.5], [47.586,4.444,-27.546], color="red blue", name="Arrows_17.0990009308_8")

cmd.load_cgo(cluster_dict["17.0990009308"], "Features_17.0990009308", 1)
cmd.load_cgo(cluster_dict["17.0990009308_arrows"], "Arrows_17.0990009308")
cmd.set("transparency", 0.2,"Features_17.0990009308")
cmd.group("Pharmacophore_17.0990009308", members="Features_17.0990009308")
cmd.group("Pharmacophore_17.0990009308", members="Arrows_17.0990009308")

if dirpath:
    f = join(dirpath, "1/label_threshold_17.0990009308.mol2")
else:
    f = "1/label_threshold_17.0990009308.mol2"

cmd.load(f, 'label_threshold_17.0990009308')
cmd.hide('everything', 'label_threshold_17.0990009308')
cmd.label("label_threshold_17.0990009308", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.0990009308', members= 'label_threshold_17.0990009308')


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
    f = join(dirpath, "2/label_threshold_5.0.mol2")
else:
    f = "2/label_threshold_5.0.mol2"

cmd.load(f, 'label_threshold_5.0')
cmd.hide('everything', 'label_threshold_5.0')
cmd.label("label_threshold_5.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [5.0]
gfiles = ['2/apolar.grd']
grids = ['apolar']
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


cluster_dict = {"10.3554997444":[], "10.3554997444_arrows":[]}

cluster_dict["10.3554997444"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(33.5), float(-4.5), float(-18.5), float(1.0)]

cluster_dict["10.3554997444_arrows"] += cgo_arrow([33.5,-4.5,-18.5], [32.663,-1.547,-19.815], color="blue red", name="Arrows_10.3554997444_1")

cluster_dict["10.3554997444"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(-15.0), float(-26.5), float(1.0)]

cluster_dict["10.3554997444_arrows"] += cgo_arrow([36.5,-15.0,-26.5], [38.138,-15.746,-24.174], color="blue red", name="Arrows_10.3554997444_2")

cluster_dict["10.3554997444"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(-3.5), float(-20.0), float(1.0)]

cluster_dict["10.3554997444_arrows"] += cgo_arrow([36.5,-3.5,-20.0], [38.057,-3.157,-21.901], color="blue red", name="Arrows_10.3554997444_3")

cluster_dict["10.3554997444"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(35.1347682069), float(-10.5274588972), float(-24.4289468181), float(1.0)]


cluster_dict["10.3554997444"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(37.0), float(-14.5), float(-28.0), float(1.0)]

cluster_dict["10.3554997444_arrows"] += cgo_arrow([37.0,-14.5,-28.0], [40.026,-15.294,-30.305], color="red blue", name="Arrows_10.3554997444_4")

cluster_dict["10.3554997444"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(36.0), float(-5.5), float(-21.0), float(1.0)]

cluster_dict["10.3554997444_arrows"] += cgo_arrow([36.0,-5.5,-21.0], [38.057,-3.157,-21.901], color="red blue", name="Arrows_10.3554997444_5")

cluster_dict["10.3554997444"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(37.0), float(-14.5), float(-28.0), float(1.0)]

cluster_dict["10.3554997444_arrows"] += cgo_arrow([37.0,-14.5,-28.0], [40.026,-15.294,-30.305], color="red blue", name="Arrows_10.3554997444_6")

cmd.load_cgo(cluster_dict["10.3554997444"], "Features_10.3554997444", 1)
cmd.load_cgo(cluster_dict["10.3554997444_arrows"], "Arrows_10.3554997444")
cmd.set("transparency", 0.2,"Features_10.3554997444")
cmd.group("Pharmacophore_10.3554997444", members="Features_10.3554997444")
cmd.group("Pharmacophore_10.3554997444", members="Arrows_10.3554997444")

if dirpath:
    f = join(dirpath, "2/label_threshold_10.3554997444.mol2")
else:
    f = "2/label_threshold_10.3554997444.mol2"

cmd.load(f, 'label_threshold_10.3554997444')
cmd.hide('everything', 'label_threshold_10.3554997444')
cmd.label("label_threshold_10.3554997444", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.3554997444', members= 'label_threshold_10.3554997444')


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
