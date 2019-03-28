
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
    f = join(dirpath, "0/label_threshold_8.5.mol2")
else:
    f = "0/label_threshold_8.5.mol2"

cmd.load(f, 'label_threshold_8.5')
cmd.hide('everything', 'label_threshold_8.5')
cmd.label("label_threshold_8.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [8.5]
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


cluster_dict = {"11.7354998589":[], "11.7354998589_arrows":[]}

cluster_dict["11.7354998589"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(13.0), float(32.0), float(1.0)]

cluster_dict["11.7354998589_arrows"] += cgo_arrow([23.0,13.0,32.0], [20.399,12.771,31.489], color="blue red", name="Arrows_11.7354998589_1")

cluster_dict["11.7354998589"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(23.1563290441), float(16.9766663534), float(34.0816572156), float(1.0)]


cluster_dict["11.7354998589"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(20.0), float(21.2417573617), float(34.0), float(1.0)]


cluster_dict["11.7354998589"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(12.5), float(33.5), float(1.0)]

cluster_dict["11.7354998589_arrows"] += cgo_arrow([22.0,12.5,33.5], [20.399,12.771,31.489], color="red blue", name="Arrows_11.7354998589_2")

cluster_dict["11.7354998589"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(6.0), float(30.0), float(1.0)]

cluster_dict["11.7354998589_arrows"] += cgo_arrow([28.0,6.0,30.0], [30.726,8.722,27.538], color="red blue", name="Arrows_11.7354998589_3")

cmd.load_cgo(cluster_dict["11.7354998589"], "Features_11.7354998589", 1)
cmd.load_cgo(cluster_dict["11.7354998589_arrows"], "Arrows_11.7354998589")
cmd.set("transparency", 0.2,"Features_11.7354998589")
cmd.group("Pharmacophore_11.7354998589", members="Features_11.7354998589")
cmd.group("Pharmacophore_11.7354998589", members="Arrows_11.7354998589")

if dirpath:
    f = join(dirpath, "0/label_threshold_11.7354998589.mol2")
else:
    f = "0/label_threshold_11.7354998589.mol2"

cmd.load(f, 'label_threshold_11.7354998589')
cmd.hide('everything', 'label_threshold_11.7354998589')
cmd.label("label_threshold_11.7354998589", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.7354998589', members= 'label_threshold_11.7354998589')


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
    f = join(dirpath, "1/label_threshold_8.3.mol2")
else:
    f = "1/label_threshold_8.3.mol2"

cmd.load(f, 'label_threshold_8.3')
cmd.hide('everything', 'label_threshold_8.3')
cmd.label("label_threshold_8.3", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [8.3]
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


cluster_dict = {"11.6099996567":[], "11.6099996567_arrows":[]}

cluster_dict["11.6099996567"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(13.0), float(32.0), float(1.0)]

cluster_dict["11.6099996567_arrows"] += cgo_arrow([23.0,13.0,32.0], [20.399,12.771,31.489], color="blue red", name="Arrows_11.6099996567_1")

cluster_dict["11.6099996567"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(23.4099255862), float(16.4748294713), float(34.0025782316), float(1.0)]


cluster_dict["11.6099996567"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(12.5), float(33.5), float(1.0)]

cluster_dict["11.6099996567_arrows"] += cgo_arrow([22.0,12.5,33.5], [20.399,12.771,31.489], color="red blue", name="Arrows_11.6099996567_2")

cluster_dict["11.6099996567"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(6.0), float(30.0), float(1.0)]

cluster_dict["11.6099996567_arrows"] += cgo_arrow([28.0,6.0,30.0], [30.726,8.722,27.538], color="red blue", name="Arrows_11.6099996567_3")

cmd.load_cgo(cluster_dict["11.6099996567"], "Features_11.6099996567", 1)
cmd.load_cgo(cluster_dict["11.6099996567_arrows"], "Arrows_11.6099996567")
cmd.set("transparency", 0.2,"Features_11.6099996567")
cmd.group("Pharmacophore_11.6099996567", members="Features_11.6099996567")
cmd.group("Pharmacophore_11.6099996567", members="Arrows_11.6099996567")

if dirpath:
    f = join(dirpath, "1/label_threshold_11.6099996567.mol2")
else:
    f = "1/label_threshold_11.6099996567.mol2"

cmd.load(f, 'label_threshold_11.6099996567')
cmd.hide('everything', 'label_threshold_11.6099996567')
cmd.label("label_threshold_11.6099996567", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.6099996567', members= 'label_threshold_11.6099996567')


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

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
