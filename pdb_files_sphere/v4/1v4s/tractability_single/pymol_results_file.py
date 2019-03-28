
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


cluster_dict = {"15.7950000763":[], "15.7950000763_arrows":[]}

cluster_dict["15.7950000763"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(19.0), float(57.0), float(1.0)]

cluster_dict["15.7950000763_arrows"] += cgo_arrow([36.5,19.0,57.0], [34.748,16.322,56.81], color="blue red", name="Arrows_15.7950000763_1")

cluster_dict["15.7950000763"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(38.5), float(16.0), float(60.5), float(1.0)]

cluster_dict["15.7950000763_arrows"] += cgo_arrow([38.5,16.0,60.5], [40.059,15.519,58.856], color="blue red", name="Arrows_15.7950000763_2")

cluster_dict["15.7950000763"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(40.0), float(23.0), float(60.0), float(1.0)]

cluster_dict["15.7950000763_arrows"] += cgo_arrow([40.0,23.0,60.0], [42.68,21.965,59.695], color="blue red", name="Arrows_15.7950000763_3")

cluster_dict["15.7950000763"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(42.5), float(9.5), float(64.0), float(1.0)]

cluster_dict["15.7950000763_arrows"] += cgo_arrow([42.5,9.5,64.0], [43.007,7.661,65.091], color="blue red", name="Arrows_15.7950000763_4")

cluster_dict["15.7950000763"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(39.5613346205), float(16.3112290136), float(61.3079713059), float(1.0)]


cluster_dict["15.7950000763"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(36.0), float(16.5), float(60.5), float(1.0)]

cluster_dict["15.7950000763_arrows"] += cgo_arrow([36.0,16.5,60.5], [37.483,15.325,57.746], color="red blue", name="Arrows_15.7950000763_5")

cluster_dict["15.7950000763"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(36.0), float(16.5), float(60.5), float(1.0)]

cluster_dict["15.7950000763_arrows"] += cgo_arrow([36.0,16.5,60.5], [37.483,15.325,57.746], color="red blue", name="Arrows_15.7950000763_6")

cluster_dict["15.7950000763"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(40.5), float(10.5), float(65.0), float(1.0)]

cluster_dict["15.7950000763_arrows"] += cgo_arrow([40.5,10.5,65.0], [37.79,10.544,67.284], color="red blue", name="Arrows_15.7950000763_7")

cluster_dict["15.7950000763"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(43.0), float(13.0), float(61.0), float(1.0)]

cluster_dict["15.7950000763_arrows"] += cgo_arrow([43.0,13.0,61.0], [43.012,15.508,59.647], color="red blue", name="Arrows_15.7950000763_8")

cmd.load_cgo(cluster_dict["15.7950000763"], "Features_15.7950000763", 1)
cmd.load_cgo(cluster_dict["15.7950000763_arrows"], "Arrows_15.7950000763")
cmd.set("transparency", 0.2,"Features_15.7950000763")
cmd.group("Pharmacophore_15.7950000763", members="Features_15.7950000763")
cmd.group("Pharmacophore_15.7950000763", members="Arrows_15.7950000763")

if dirpath:
    f = join(dirpath, "label_threshold_15.7950000763.mol2")
else:
    f = "label_threshold_15.7950000763.mol2"

cmd.load(f, 'label_threshold_15.7950000763')
cmd.hide('everything', 'label_threshold_15.7950000763')
cmd.label("label_threshold_15.7950000763", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.7950000763', members= 'label_threshold_15.7950000763')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
