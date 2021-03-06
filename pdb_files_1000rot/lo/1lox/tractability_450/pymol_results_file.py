
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


cluster_dict = {"15.4099998474":[], "15.4099998474_arrows":[]}

cluster_dict["15.4099998474"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-27.5), float(150.0), float(58.0), float(1.0)]

cluster_dict["15.4099998474_arrows"] += cgo_arrow([-27.5,150.0,58.0], [-28.645,151.11,60.419], color="blue red", name="Arrows_15.4099998474_1")

cluster_dict["15.4099998474"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-27.0), float(152.5), float(59.5), float(1.0)]

cluster_dict["15.4099998474_arrows"] += cgo_arrow([-27.0,152.5,59.5], [-28.645,151.11,60.419], color="blue red", name="Arrows_15.4099998474_2")

cluster_dict["15.4099998474"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-24.5), float(150.0), float(56.5), float(1.0)]

cluster_dict["15.4099998474_arrows"] += cgo_arrow([-24.5,150.0,56.5], [-21.963,148.062,56.622], color="blue red", name="Arrows_15.4099998474_3")

cluster_dict["15.4099998474"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-27.0549381753), float(152.149255437), float(55.1518974655), float(1.0)]


cluster_dict["15.4099998474"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-31.0), float(150.5), float(55.0), float(1.0)]

cluster_dict["15.4099998474_arrows"] += cgo_arrow([-31.0,150.5,55.0], [-34.61,150.654,54.431], color="red blue", name="Arrows_15.4099998474_4")

cluster_dict["15.4099998474"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-28.5), float(150.0), float(57.5), float(1.0)]

cluster_dict["15.4099998474_arrows"] += cgo_arrow([-28.5,150.0,57.5], [-30.336,148.482,58.963], color="red blue", name="Arrows_15.4099998474_5")

cluster_dict["15.4099998474"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-25.0), float(151.0), float(47.0), float(1.0)]

cluster_dict["15.4099998474_arrows"] += cgo_arrow([-25.0,151.0,47.0], [-22.506,150.484,44.302], color="red blue", name="Arrows_15.4099998474_6")

cluster_dict["15.4099998474"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-25.0), float(154.5), float(58.5), float(1.0)]

cluster_dict["15.4099998474_arrows"] += cgo_arrow([-25.0,154.5,58.5], [-23.239,151.396,61.984], color="red blue", name="Arrows_15.4099998474_7")

cluster_dict["15.4099998474"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-24.5), float(153.0), float(51.5), float(1.0)]

cluster_dict["15.4099998474_arrows"] += cgo_arrow([-24.5,153.0,51.5], [-21.467,154.331,51.348], color="red blue", name="Arrows_15.4099998474_8")

cluster_dict["15.4099998474"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-21.5), float(151.0), float(54.0), float(1.0)]

cluster_dict["15.4099998474_arrows"] += cgo_arrow([-21.5,151.0,54.0], [-18.5,150.526,55.499], color="red blue", name="Arrows_15.4099998474_9")

cmd.load_cgo(cluster_dict["15.4099998474"], "Features_15.4099998474", 1)
cmd.load_cgo(cluster_dict["15.4099998474_arrows"], "Arrows_15.4099998474")
cmd.set("transparency", 0.2,"Features_15.4099998474")
cmd.group("Pharmacophore_15.4099998474", members="Features_15.4099998474")
cmd.group("Pharmacophore_15.4099998474", members="Arrows_15.4099998474")

if dirpath:
    f = join(dirpath, "label_threshold_15.4099998474.mol2")
else:
    f = "label_threshold_15.4099998474.mol2"

cmd.load(f, 'label_threshold_15.4099998474')
cmd.hide('everything', 'label_threshold_15.4099998474')
cmd.label("label_threshold_15.4099998474", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.4099998474', members= 'label_threshold_15.4099998474')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
