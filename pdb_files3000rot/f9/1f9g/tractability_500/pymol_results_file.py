
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


cluster_dict = {"14.0360002518":[], "14.0360002518_arrows":[]}

cluster_dict["14.0360002518"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(2.5), float(-13.0), float(35.0), float(1.0)]

cluster_dict["14.0360002518_arrows"] += cgo_arrow([2.5,-13.0,35.0], [4.182,-11.8,37.801], color="blue red", name="Arrows_14.0360002518_1")

cluster_dict["14.0360002518"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(2.5), float(-13.0), float(35.0), float(1.0)]

cluster_dict["14.0360002518_arrows"] += cgo_arrow([2.5,-13.0,35.0], [4.182,-11.8,37.801], color="blue red", name="Arrows_14.0360002518_2")

cluster_dict["14.0360002518"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(11.5), float(-14.0), float(31.0), float(1.0)]

cluster_dict["14.0360002518_arrows"] += cgo_arrow([11.5,-14.0,31.0], [10.645,-16.255,29.051], color="blue red", name="Arrows_14.0360002518_3")

cluster_dict["14.0360002518"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(11.0), float(-11.5), float(32.5), float(1.0)]

cluster_dict["14.0360002518_arrows"] += cgo_arrow([11.0,-11.5,32.5], [10.675,-10.839,35.78], color="blue red", name="Arrows_14.0360002518_4")

cluster_dict["14.0360002518"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(11.0), float(-11.5), float(32.5), float(1.0)]

cluster_dict["14.0360002518_arrows"] += cgo_arrow([11.0,-11.5,32.5], [10.675,-10.839,35.78], color="blue red", name="Arrows_14.0360002518_5")

cluster_dict["14.0360002518"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(6.94072383761), float(-11.1170038627), float(32.3733264354), float(1.0)]


cluster_dict["14.0360002518"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(3.80694168047), float(-6.5), float(34.2959982679), float(1.0)]


cluster_dict["14.0360002518"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(2.0), float(-11.5), float(27.0), float(1.0)]

cluster_dict["14.0360002518_arrows"] += cgo_arrow([2.0,-11.5,27.0], [3.279,-10.274,24.325], color="red blue", name="Arrows_14.0360002518_6")

cluster_dict["14.0360002518"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(3.0), float(-14.0), float(35.5), float(1.0)]

cluster_dict["14.0360002518_arrows"] += cgo_arrow([3.0,-14.0,35.5], [4.182,-11.8,37.801], color="red blue", name="Arrows_14.0360002518_7")

cluster_dict["14.0360002518"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(2.0), float(-11.5), float(27.0), float(1.0)]

cluster_dict["14.0360002518_arrows"] += cgo_arrow([2.0,-11.5,27.0], [3.279,-10.274,24.325], color="red blue", name="Arrows_14.0360002518_8")

cluster_dict["14.0360002518"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(10.5), float(-13.5), float(32.0), float(1.0)]

cluster_dict["14.0360002518_arrows"] += cgo_arrow([10.5,-13.5,32.0], [10.136,-14.416,34.373], color="red blue", name="Arrows_14.0360002518_9")

cluster_dict["14.0360002518"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(13.0), float(-14.5), float(31.5), float(1.0)]

cluster_dict["14.0360002518_arrows"] += cgo_arrow([13.0,-14.5,31.5], [12.103,-17.645,30.042], color="red blue", name="Arrows_14.0360002518_10")

cmd.load_cgo(cluster_dict["14.0360002518"], "Features_14.0360002518", 1)
cmd.load_cgo(cluster_dict["14.0360002518_arrows"], "Arrows_14.0360002518")
cmd.set("transparency", 0.2,"Features_14.0360002518")
cmd.group("Pharmacophore_14.0360002518", members="Features_14.0360002518")
cmd.group("Pharmacophore_14.0360002518", members="Arrows_14.0360002518")

if dirpath:
    f = join(dirpath, "label_threshold_14.0360002518.mol2")
else:
    f = "label_threshold_14.0360002518.mol2"

cmd.load(f, 'label_threshold_14.0360002518')
cmd.hide('everything', 'label_threshold_14.0360002518')
cmd.label("label_threshold_14.0360002518", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.0360002518', members= 'label_threshold_14.0360002518')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")