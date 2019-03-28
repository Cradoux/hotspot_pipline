
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


cluster_dict = {"16.4039993286":[], "16.4039993286_arrows":[]}

cluster_dict["16.4039993286"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-17.5), float(39.0), float(159.0), float(1.0)]

cluster_dict["16.4039993286_arrows"] += cgo_arrow([-17.5,39.0,159.0], [-17.541,35.942,159.367], color="blue red", name="Arrows_16.4039993286_1")

cluster_dict["16.4039993286"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-14.5), float(43.5), float(154.0), float(1.0)]

cluster_dict["16.4039993286_arrows"] += cgo_arrow([-14.5,43.5,154.0], [-15.183,43.615,151.35], color="blue red", name="Arrows_16.4039993286_2")

cluster_dict["16.4039993286"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-11.5), float(44.5), float(156.0), float(1.0)]

cluster_dict["16.4039993286_arrows"] += cgo_arrow([-11.5,44.5,156.0], [-10.598,46.675,158.06], color="blue red", name="Arrows_16.4039993286_3")

cluster_dict["16.4039993286"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-8.5), float(46.5), float(156.5), float(1.0)]

cluster_dict["16.4039993286_arrows"] += cgo_arrow([-8.5,46.5,156.5], [-10.598,46.675,158.06], color="blue red", name="Arrows_16.4039993286_4")

cluster_dict["16.4039993286"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-7.5), float(41.5), float(154.0), float(1.0)]

cluster_dict["16.4039993286_arrows"] += cgo_arrow([-7.5,41.5,154.0], [-4.263,41.005,152.854], color="blue red", name="Arrows_16.4039993286_5")

cluster_dict["16.4039993286"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-9.48662565057), float(42.353400568), float(155.857394719), float(1.0)]


cluster_dict["16.4039993286"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-12.5), float(42.0), float(153.0), float(1.0)]


cluster_dict["16.4039993286"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-10.5), float(45.0), float(156.0), float(1.0)]

cluster_dict["16.4039993286_arrows"] += cgo_arrow([-10.5,45.0,156.0], [-10.598,46.675,158.06], color="red blue", name="Arrows_16.4039993286_6")

cluster_dict["16.4039993286"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-13.5), float(43.0), float(153.5), float(1.0)]

cluster_dict["16.4039993286_arrows"] += cgo_arrow([-13.5,43.0,153.5], [-15.183,43.615,151.35], color="red blue", name="Arrows_16.4039993286_7")

cluster_dict["16.4039993286"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-9.0), float(43.0), float(151.0), float(1.0)]

cluster_dict["16.4039993286_arrows"] += cgo_arrow([-9.0,43.0,151.0], [-8.173,39.85,150.66], color="red blue", name="Arrows_16.4039993286_8")

cluster_dict["16.4039993286"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-6.5), float(42.0), float(155.0), float(1.0)]

cluster_dict["16.4039993286_arrows"] += cgo_arrow([-6.5,42.0,155.0], [-4.263,41.005,152.854], color="red blue", name="Arrows_16.4039993286_9")

cmd.load_cgo(cluster_dict["16.4039993286"], "Features_16.4039993286", 1)
cmd.load_cgo(cluster_dict["16.4039993286_arrows"], "Arrows_16.4039993286")
cmd.set("transparency", 0.2,"Features_16.4039993286")
cmd.group("Pharmacophore_16.4039993286", members="Features_16.4039993286")
cmd.group("Pharmacophore_16.4039993286", members="Arrows_16.4039993286")

if dirpath:
    f = join(dirpath, "label_threshold_16.4039993286.mol2")
else:
    f = "label_threshold_16.4039993286.mol2"

cmd.load(f, 'label_threshold_16.4039993286')
cmd.hide('everything', 'label_threshold_16.4039993286')
cmd.label("label_threshold_16.4039993286", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.4039993286', members= 'label_threshold_16.4039993286')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
