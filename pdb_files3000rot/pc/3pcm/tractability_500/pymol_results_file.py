
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
zip_dir = 'out.zip'
with zipfile.ZipFile(zip_dir) as hs_zip:
    hs_zip.extractall(dirpath)

cmd.load(join(dirpath,"protein.pdb"), "protein")
cmd.show("cartoon", "protein")

if dirpath:
    f = join(dirpath, "label_threshold_10.mol2")
else:
    f = "label_threshold_10.mol2"

cmd.load(f, 'label_threshold_10')
cmd.hide('everything', 'label_threshold_10')
cmd.label("label_threshold_10", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


if dirpath:
    f = join(dirpath, "label_threshold_14.mol2")
else:
    f = "label_threshold_14.mol2"

cmd.load(f, 'label_threshold_14')
cmd.hide('everything', 'label_threshold_14')
cmd.label("label_threshold_14", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


if dirpath:
    f = join(dirpath, "label_threshold_17.mol2")
else:
    f = "label_threshold_17.mol2"

cmd.load(f, 'label_threshold_17')
cmd.hide('everything', 'label_threshold_17')
cmd.label("label_threshold_17", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10, 14, 17]
gfiles = ['donor.grd', 'apolar.grd', 'acceptor.grd']
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


cluster_dict = {"14.1840000153":[], "14.1840000153_arrows":[]}

cluster_dict["14.1840000153"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-44.0), float(-14.5), float(23.5), float(1.0)]

cluster_dict["14.1840000153_arrows"] += cgo_arrow([-44.0,-14.5,23.5], [-45.333,-13.303,25.639], color="blue red", name="Arrows_14.1840000153_1")

cluster_dict["14.1840000153"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-41.5), float(-23.0), float(25.0), float(1.0)]

cluster_dict["14.1840000153_arrows"] += cgo_arrow([-41.5,-23.0,25.0], [-41.401,-25.107,23.037], color="blue red", name="Arrows_14.1840000153_2")

cluster_dict["14.1840000153"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-38.5), float(-19.0), float(25.0), float(1.0)]

cluster_dict["14.1840000153_arrows"] += cgo_arrow([-38.5,-19.0,25.0], [-36.914,-21.867,27.084], color="blue red", name="Arrows_14.1840000153_3")

cluster_dict["14.1840000153"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-39.0), float(-17.5), float(21.0), float(1.0)]

cluster_dict["14.1840000153_arrows"] += cgo_arrow([-39.0,-17.5,21.0], [-39.22,-17.997,18.326], color="blue red", name="Arrows_14.1840000153_4")

cluster_dict["14.1840000153"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-41.1631287375), float(-16.4309674001), float(22.79440645), float(1.0)]


cluster_dict["14.1840000153"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-38.9375), float(-19.9375), float(27.5), float(1.0)]


cluster_dict["14.1840000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-43.0), float(-15.0), float(25.0), float(1.0)]

cluster_dict["14.1840000153_arrows"] += cgo_arrow([-43.0,-15.0,25.0], [-45.577,-16.324,26.096], color="red blue", name="Arrows_14.1840000153_5")

cluster_dict["14.1840000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-41.0), float(-22.5), float(25.0), float(1.0)]

cluster_dict["14.1840000153_arrows"] += cgo_arrow([-41.0,-22.5,25.0], [-39.053,-23.756,26.626], color="red blue", name="Arrows_14.1840000153_6")

cluster_dict["14.1840000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-41.0), float(-16.5), float(17.5), float(1.0)]

cluster_dict["14.1840000153_arrows"] += cgo_arrow([-41.0,-16.5,17.5], [-39.22,-17.997,18.326], color="red blue", name="Arrows_14.1840000153_7")

cluster_dict["14.1840000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-40.0), float(-18.5), float(21.0), float(1.0)]

cluster_dict["14.1840000153_arrows"] += cgo_arrow([-40.0,-18.5,21.0], [-39.22,-17.997,18.326], color="red blue", name="Arrows_14.1840000153_8")

cluster_dict["14.1840000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-39.0), float(-12.5), float(23.5), float(1.0)]

cluster_dict["14.1840000153_arrows"] += cgo_arrow([-39.0,-12.5,23.5], [-35.867,-11.387,23.771], color="red blue", name="Arrows_14.1840000153_9")

cmd.load_cgo(cluster_dict["14.1840000153"], "Features_14.1840000153", 1)
cmd.load_cgo(cluster_dict["14.1840000153_arrows"], "Arrows_14.1840000153")
cmd.set("transparency", 0.2,"Features_14.1840000153")
cmd.group("Pharmacophore_14.1840000153", members="Features_14.1840000153")
cmd.group("Pharmacophore_14.1840000153", members="Arrows_14.1840000153")

if dirpath:
    f = join(dirpath, "label_threshold_14.1840000153.mol2")
else:
    f = "label_threshold_14.1840000153.mol2"

cmd.load(f, 'label_threshold_14.1840000153')
cmd.hide('everything', 'label_threshold_14.1840000153')
cmd.label("label_threshold_14.1840000153", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.1840000153', members= 'label_threshold_14.1840000153')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
