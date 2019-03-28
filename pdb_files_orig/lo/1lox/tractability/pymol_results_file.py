
from os.path import join
import tempfile
import zipfile
from pymol import cmd
from pymol.cgo import *

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


cluster_dict = {"16.6609992981":[], "16.6609992981_arrows":[]}

cluster_dict["16.6609992981"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-12.5), float(139.5), float(46.0), float(1.0)]

cluster_dict["16.6609992981_arrows"] += cgo_arrow([-12.5,139.5,46.0], [-13.941,136.939,47.444], color="blue red", name="Arrows_16.6609992981_1")

cluster_dict["16.6609992981"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-12.5), float(139.5), float(46.0), float(1.0)]

cluster_dict["16.6609992981_arrows"] += cgo_arrow([-12.5,139.5,46.0], [-13.941,136.939,47.444], color="blue red", name="Arrows_16.6609992981_2")

cluster_dict["16.6609992981"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-10.0), float(138.5), float(49.5), float(1.0)]

cluster_dict["16.6609992981_arrows"] += cgo_arrow([-10.0,138.5,49.5], [-7.298,137.092,49.558], color="blue red", name="Arrows_16.6609992981_3")

cluster_dict["16.6609992981"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-6.5), float(139.0), float(49.5), float(1.0)]

cluster_dict["16.6609992981_arrows"] += cgo_arrow([-6.5,139.0,49.5], [-7.298,137.092,49.558], color="blue red", name="Arrows_16.6609992981_4")

cluster_dict["16.6609992981"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-2.5), float(143.5), float(46.5), float(1.0)]

cluster_dict["16.6609992981_arrows"] += cgo_arrow([-2.5,143.5,46.5], [-4.668,144.707,44.258], color="blue red", name="Arrows_16.6609992981_5")

cluster_dict["16.6609992981"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-7.94463434977), float(140.499283839), float(48.0233458023), float(1.0)]


cluster_dict["16.6609992981"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-2.36704956368), float(138.001197199), float(49.1765858838), float(1.0)]


cluster_dict["16.6609992981"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-1.49957541703), float(141.256299762), float(46.7562997623), float(1.0)]


cluster_dict["16.6609992981"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-0.5), float(139.250513625), float(47.1884939574), float(1.0)]


cluster_dict["16.6609992981"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-15.5), float(133.0), float(48.5), float(1.0)]

cluster_dict["16.6609992981_arrows"] += cgo_arrow([-15.5,133.0,48.5], [-13.181,131.247,49.727], color="red blue", name="Arrows_16.6609992981_6")

cluster_dict["16.6609992981"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-13.5), float(135.0), float(45.5), float(1.0)]

cluster_dict["16.6609992981_arrows"] += cgo_arrow([-13.5,135.0,45.5], [-13.941,136.939,47.444], color="red blue", name="Arrows_16.6609992981_7")

cluster_dict["16.6609992981"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-10.0), float(140.5), float(45.0), float(1.0)]


cluster_dict["16.6609992981"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-2.0), float(143.5), float(48.0), float(1.0)]

cluster_dict["16.6609992981_arrows"] += cgo_arrow([-2.0,143.5,48.0], [-0.44,144.992,50.646], color="red blue", name="Arrows_16.6609992981_8")

cluster_dict["16.6609992981"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-2.0), float(143.5), float(48.0), float(1.0)]

cluster_dict["16.6609992981_arrows"] += cgo_arrow([-2.0,143.5,48.0], [-0.44,144.992,50.646], color="red blue", name="Arrows_16.6609992981_9")

cmd.load_cgo(cluster_dict["16.6609992981"], "Features_16.6609992981", 1)
cmd.load_cgo(cluster_dict["16.6609992981_arrows"], "Arrows_16.6609992981")
cmd.set("transparency", 0.2,"Features_16.6609992981")
cmd.group("Pharmacophore_16.6609992981", members="Features_16.6609992981")
cmd.group("Pharmacophore_16.6609992981", members="Arrows_16.6609992981")

if dirpath:
    f = join(dirpath, "label_threshold_16.6609992981.mol2")
else:
    f = "label_threshold_16.6609992981.mol2"

cmd.load(f, 'label_threshold_16.6609992981')
cmd.hide('everything', 'label_threshold_16.6609992981')
cmd.label("label_threshold_16.6609992981", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.6609992981', members= 'label_threshold_16.6609992981')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
