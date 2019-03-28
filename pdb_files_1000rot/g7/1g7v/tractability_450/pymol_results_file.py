
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


cluster_dict = {"15.9200000763":[], "15.9200000763_arrows":[]}

cluster_dict["15.9200000763"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(5.5), float(102.0), float(82.5), float(1.0)]

cluster_dict["15.9200000763_arrows"] += cgo_arrow([5.5,102.0,82.5], [6.438,101.52,85.208], color="blue red", name="Arrows_15.9200000763_1")

cluster_dict["15.9200000763"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(9.0), float(106.0), float(79.5), float(1.0)]

cluster_dict["15.9200000763_arrows"] += cgo_arrow([9.0,106.0,79.5], [11.728,106.388,80.182], color="blue red", name="Arrows_15.9200000763_2")

cluster_dict["15.9200000763"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(9.0), float(106.0), float(79.5), float(1.0)]

cluster_dict["15.9200000763_arrows"] += cgo_arrow([9.0,106.0,79.5], [11.728,106.388,80.182], color="blue red", name="Arrows_15.9200000763_3")

cluster_dict["15.9200000763"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(11.0), float(110.5), float(81.0), float(1.0)]

cluster_dict["15.9200000763_arrows"] += cgo_arrow([11.0,110.5,81.0], [13.419,110.042,82.914], color="blue red", name="Arrows_15.9200000763_4")

cluster_dict["15.9200000763"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(7.46469903546), float(107.966396345), float(78.2984931418), float(1.0)]


cluster_dict["15.9200000763"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(6.0), float(103.5), float(79.0), float(1.0)]

cluster_dict["15.9200000763_arrows"] += cgo_arrow([6.0,103.5,79.0], [5.046,101.486,78.2], color="red blue", name="Arrows_15.9200000763_5")

cluster_dict["15.9200000763"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(8.5), float(115.0), float(78.0), float(1.0)]

cluster_dict["15.9200000763_arrows"] += cgo_arrow([8.5,115.0,78.0], [6.555,117.389,78.438], color="red blue", name="Arrows_15.9200000763_6")

cluster_dict["15.9200000763"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(10.5), float(109.0), float(79.5), float(1.0)]

cluster_dict["15.9200000763_arrows"] += cgo_arrow([10.5,109.0,79.5], [11.883,106.786,77.978], color="red blue", name="Arrows_15.9200000763_7")

cluster_dict["15.9200000763"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(11.5), float(110.0), float(81.0), float(1.0)]

cluster_dict["15.9200000763_arrows"] += cgo_arrow([11.5,110.0,81.0], [13.419,110.042,82.914], color="red blue", name="Arrows_15.9200000763_8")

cluster_dict["15.9200000763"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(11.5), float(113.5), float(78.0), float(1.0)]

cluster_dict["15.9200000763_arrows"] += cgo_arrow([11.5,113.5,78.0], [13.328,115.357,77.05], color="red blue", name="Arrows_15.9200000763_9")

cmd.load_cgo(cluster_dict["15.9200000763"], "Features_15.9200000763", 1)
cmd.load_cgo(cluster_dict["15.9200000763_arrows"], "Arrows_15.9200000763")
cmd.set("transparency", 0.2,"Features_15.9200000763")
cmd.group("Pharmacophore_15.9200000763", members="Features_15.9200000763")
cmd.group("Pharmacophore_15.9200000763", members="Arrows_15.9200000763")

if dirpath:
    f = join(dirpath, "label_threshold_15.9200000763.mol2")
else:
    f = "label_threshold_15.9200000763.mol2"

cmd.load(f, 'label_threshold_15.9200000763')
cmd.hide('everything', 'label_threshold_15.9200000763')
cmd.label("label_threshold_15.9200000763", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.9200000763', members= 'label_threshold_15.9200000763')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
