
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


cluster_dict = {"15.1700000763":[], "15.1700000763_arrows":[]}

cluster_dict["15.1700000763"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(39.0), float(33.0), float(42.0), float(1.0)]

cluster_dict["15.1700000763_arrows"] += cgo_arrow([39.0,33.0,42.0], [36.909,31.21,41.407], color="blue red", name="Arrows_15.1700000763_1")

cluster_dict["15.1700000763"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(43.0), float(37.5), float(47.0), float(1.0)]

cluster_dict["15.1700000763_arrows"] += cgo_arrow([43.0,37.5,47.0], [45.075,37.341,45.381], color="blue red", name="Arrows_15.1700000763_2")

cluster_dict["15.1700000763"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(41.0), float(25.0), float(42.5), float(1.0)]

cluster_dict["15.1700000763_arrows"] += cgo_arrow([41.0,25.0,42.5], [40.137,22.295,42.592], color="blue red", name="Arrows_15.1700000763_3")

cluster_dict["15.1700000763"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(42.0), float(28.0), float(40.5), float(1.0)]

cluster_dict["15.1700000763_arrows"] += cgo_arrow([42.0,28.0,40.5], [44.787,26.78,38.927], color="blue red", name="Arrows_15.1700000763_4")

cluster_dict["15.1700000763"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(42.5), float(25.0), float(40.5), float(1.0)]

cluster_dict["15.1700000763_arrows"] += cgo_arrow([42.5,25.0,40.5], [44.787,26.78,38.927], color="blue red", name="Arrows_15.1700000763_5")

cluster_dict["15.1700000763"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(44.0), float(33.5), float(38.0), float(1.0)]

cluster_dict["15.1700000763_arrows"] += cgo_arrow([44.0,33.5,38.0], [46.638,32.89,38.91], color="blue red", name="Arrows_15.1700000763_6")

cluster_dict["15.1700000763"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(44.0), float(34.0), float(45.5), float(1.0)]

cluster_dict["15.1700000763_arrows"] += cgo_arrow([44.0,34.0,45.5], [45.159,34.188,48.359], color="blue red", name="Arrows_15.1700000763_7")

cluster_dict["15.1700000763"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(45.5), float(27.5), float(42.5), float(1.0)]

cluster_dict["15.1700000763_arrows"] += cgo_arrow([45.5,27.5,42.5], [47.301,26.62,44.963], color="blue red", name="Arrows_15.1700000763_8")

cluster_dict["15.1700000763"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(43.3321819954), float(30.4863403479), float(43.0690563014), float(1.0)]


cluster_dict["15.1700000763"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(41.0), float(37.5), float(46.5), float(1.0)]

cluster_dict["15.1700000763_arrows"] += cgo_arrow([41.0,37.5,46.5], [39.572,37.561,45.648], color="red blue", name="Arrows_15.1700000763_9")

cluster_dict["15.1700000763"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(46.5), float(23.0), float(43.0), float(1.0)]

cluster_dict["15.1700000763_arrows"] += cgo_arrow([46.5,23.0,43.0], [48.915,22.058,41.981], color="red blue", name="Arrows_15.1700000763_10")

cmd.load_cgo(cluster_dict["15.1700000763"], "Features_15.1700000763", 1)
cmd.load_cgo(cluster_dict["15.1700000763_arrows"], "Arrows_15.1700000763")
cmd.set("transparency", 0.2,"Features_15.1700000763")
cmd.group("Pharmacophore_15.1700000763", members="Features_15.1700000763")
cmd.group("Pharmacophore_15.1700000763", members="Arrows_15.1700000763")

if dirpath:
    f = join(dirpath, "label_threshold_15.1700000763.mol2")
else:
    f = "label_threshold_15.1700000763.mol2"

cmd.load(f, 'label_threshold_15.1700000763')
cmd.hide('everything', 'label_threshold_15.1700000763')
cmd.label("label_threshold_15.1700000763", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.1700000763', members= 'label_threshold_15.1700000763')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
