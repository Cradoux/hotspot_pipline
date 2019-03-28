
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


cluster_dict = {"14.2899999619":[], "14.2899999619_arrows":[]}

cluster_dict["14.2899999619"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(59.5), float(32.0), float(117.0), float(1.0)]

cluster_dict["14.2899999619_arrows"] += cgo_arrow([59.5,32.0,117.0], [60.583,33.744,115.515], color="blue red", name="Arrows_14.2899999619_1")

cluster_dict["14.2899999619"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(64.0), float(30.0), float(113.0), float(1.0)]

cluster_dict["14.2899999619_arrows"] += cgo_arrow([64.0,30.0,113.0], [63.664,31.992,111.598], color="blue red", name="Arrows_14.2899999619_2")

cluster_dict["14.2899999619"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(64.0), float(26.5), float(117.5), float(1.0)]

cluster_dict["14.2899999619_arrows"] += cgo_arrow([64.0,26.5,117.5], [62.793,25.984,114.967], color="blue red", name="Arrows_14.2899999619_3")

cluster_dict["14.2899999619"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(67.0), float(31.5), float(115.5), float(1.0)]

cluster_dict["14.2899999619_arrows"] += cgo_arrow([67.0,31.5,115.5], [64.649,33.659,116.258], color="blue red", name="Arrows_14.2899999619_4")

cluster_dict["14.2899999619"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(62.9971122738), float(30.0506404945), float(116.567446074), float(1.0)]


cluster_dict["14.2899999619"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(60.4456041934), float(27.5), float(118.644761992), float(1.0)]


cluster_dict["14.2899999619"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(61.5), float(29.0), float(118.5), float(1.0)]

cluster_dict["14.2899999619_arrows"] += cgo_arrow([61.5,29.0,118.5], [62.626,33.117,120.114], color="red blue", name="Arrows_14.2899999619_5")

cluster_dict["14.2899999619"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(65.5), float(27.0), float(111.0), float(1.0)]

cluster_dict["14.2899999619_arrows"] += cgo_arrow([65.5,27.0,111.0], [65.367,23.887,110.637], color="red blue", name="Arrows_14.2899999619_6")

cluster_dict["14.2899999619"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(67.0), float(29.5), float(118.0), float(1.0)]

cluster_dict["14.2899999619_arrows"] += cgo_arrow([67.0,29.5,118.0], [68.071,26.881,119.097], color="red blue", name="Arrows_14.2899999619_7")

cmd.load_cgo(cluster_dict["14.2899999619"], "Features_14.2899999619", 1)
cmd.load_cgo(cluster_dict["14.2899999619_arrows"], "Arrows_14.2899999619")
cmd.set("transparency", 0.2,"Features_14.2899999619")
cmd.group("Pharmacophore_14.2899999619", members="Features_14.2899999619")
cmd.group("Pharmacophore_14.2899999619", members="Arrows_14.2899999619")

if dirpath:
    f = join(dirpath, "label_threshold_14.2899999619.mol2")
else:
    f = "label_threshold_14.2899999619.mol2"

cmd.load(f, 'label_threshold_14.2899999619')
cmd.hide('everything', 'label_threshold_14.2899999619')
cmd.label("label_threshold_14.2899999619", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.2899999619', members= 'label_threshold_14.2899999619')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
