
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


cluster_dict = {"18.0090007782":[], "18.0090007782_arrows":[]}

cluster_dict["18.0090007782"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-3.5), float(58.0), float(81.0), float(1.0)]

cluster_dict["18.0090007782_arrows"] += cgo_arrow([-3.5,58.0,81.0], [-3.493,55.183,80.502], color="blue red", name="Arrows_18.0090007782_1")

cluster_dict["18.0090007782"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-0.5), float(69.5), float(79.5), float(1.0)]

cluster_dict["18.0090007782_arrows"] += cgo_arrow([-0.5,69.5,79.5], [0.416,68.189,76.694], color="blue red", name="Arrows_18.0090007782_2")

cluster_dict["18.0090007782"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(0.0), float(61.0), float(79.5), float(1.0)]

cluster_dict["18.0090007782_arrows"] += cgo_arrow([0.0,61.0,79.5], [2.385,59.363,78.49], color="blue red", name="Arrows_18.0090007782_3")

cluster_dict["18.0090007782"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(0.0), float(65.0), float(86.0), float(1.0)]

cluster_dict["18.0090007782_arrows"] += cgo_arrow([0.0,65.0,86.0], [2.164,67.176,87.736], color="blue red", name="Arrows_18.0090007782_4")

cluster_dict["18.0090007782"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(0.417553846936), float(65.3811408994), float(82.4114182982), float(1.0)]


cluster_dict["18.0090007782"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-2.5), float(60.5), float(83.0), float(1.0)]

cluster_dict["18.0090007782_arrows"] += cgo_arrow([-2.5,60.5,83.0], [-1.978,58.881,85.643], color="red blue", name="Arrows_18.0090007782_5")

cluster_dict["18.0090007782"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-3.0), float(62.5), float(87.5), float(1.0)]

cluster_dict["18.0090007782_arrows"] += cgo_arrow([-3.0,62.5,87.5], [-4.02,62.974,90.745], color="red blue", name="Arrows_18.0090007782_6")

cluster_dict["18.0090007782"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-1.0), float(63.0), float(89.5), float(1.0)]

cluster_dict["18.0090007782_arrows"] += cgo_arrow([-1.0,63.0,89.5], [-0.829,63.238,92.424], color="red blue", name="Arrows_18.0090007782_7")

cluster_dict["18.0090007782"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-1.0), float(70.5), float(81.5), float(1.0)]

cluster_dict["18.0090007782_arrows"] += cgo_arrow([-1.0,70.5,81.5], [-0.351,72.07,85.578], color="red blue", name="Arrows_18.0090007782_8")

cluster_dict["18.0090007782"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(2.5), float(71.5), float(80.0), float(1.0)]

cluster_dict["18.0090007782_arrows"] += cgo_arrow([2.5,71.5,80.0], [2.986,72.519,77.302], color="red blue", name="Arrows_18.0090007782_9")

cmd.load_cgo(cluster_dict["18.0090007782"], "Features_18.0090007782", 1)
cmd.load_cgo(cluster_dict["18.0090007782_arrows"], "Arrows_18.0090007782")
cmd.set("transparency", 0.2,"Features_18.0090007782")
cmd.group("Pharmacophore_18.0090007782", members="Features_18.0090007782")
cmd.group("Pharmacophore_18.0090007782", members="Arrows_18.0090007782")

if dirpath:
    f = join(dirpath, "label_threshold_18.0090007782.mol2")
else:
    f = "label_threshold_18.0090007782.mol2"

cmd.load(f, 'label_threshold_18.0090007782')
cmd.hide('everything', 'label_threshold_18.0090007782')
cmd.label("label_threshold_18.0090007782", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_18.0090007782', members= 'label_threshold_18.0090007782')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
