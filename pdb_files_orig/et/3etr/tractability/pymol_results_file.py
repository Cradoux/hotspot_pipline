
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


cluster_dict = {"27.2840003967":[], "27.2840003967_arrows":[]}

cluster_dict["27.2840003967"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-41.5), float(-14.0), float(29.0), float(1.0)]

cluster_dict["27.2840003967_arrows"] += cgo_arrow([-41.5,-14.0,29.0], [-42.208,-16.006,27.081], color="blue red", name="Arrows_27.2840003967_1")

cluster_dict["27.2840003967"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-39.0), float(-17.5), float(30.5), float(1.0)]

cluster_dict["27.2840003967_arrows"] += cgo_arrow([-39.0,-17.5,30.5], [-39.859,-16.681,33.441], color="blue red", name="Arrows_27.2840003967_2")

cluster_dict["27.2840003967"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-37.5), float(-19.5), float(29.0), float(1.0)]

cluster_dict["27.2840003967_arrows"] += cgo_arrow([-37.5,-19.5,29.0], [-35.082,-20.727,28.166], color="blue red", name="Arrows_27.2840003967_3")

cluster_dict["27.2840003967"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-37.5), float(-19.5), float(29.0), float(1.0)]

cluster_dict["27.2840003967_arrows"] += cgo_arrow([-37.5,-19.5,29.0], [-35.082,-20.727,28.166], color="blue red", name="Arrows_27.2840003967_4")

cluster_dict["27.2840003967"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-38.789174964), float(-13.5281318463), float(30.4774433844), float(1.0)]


cluster_dict["27.2840003967"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-37.5), float(-12.5), float(33.0), float(1.0)]

cluster_dict["27.2840003967_arrows"] += cgo_arrow([-37.5,-12.5,33.0], [-34.828,-15.131,35.035], color="red blue", name="Arrows_27.2840003967_5")

cluster_dict["27.2840003967"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-39.5), float(-8.0), float(27.5), float(1.0)]

cluster_dict["27.2840003967_arrows"] += cgo_arrow([-39.5,-8.0,27.5], [-39.591,-6.845,29.882], color="red blue", name="Arrows_27.2840003967_6")

cluster_dict["27.2840003967"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-37.5), float(-12.5), float(33.0), float(1.0)]

cluster_dict["27.2840003967_arrows"] += cgo_arrow([-37.5,-12.5,33.0], [-34.828,-15.131,35.035], color="red blue", name="Arrows_27.2840003967_7")

cluster_dict["27.2840003967"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-38.0), float(-14.5), float(34.5), float(1.0)]

cluster_dict["27.2840003967_arrows"] += cgo_arrow([-38.0,-14.5,34.5], [-39.052,-16.878,35.506], color="red blue", name="Arrows_27.2840003967_8")

cluster_dict["27.2840003967"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(-12.5), float(27.5), float(1.0)]

cluster_dict["27.2840003967_arrows"] += cgo_arrow([-37.0,-12.5,27.5], [-36.488,-11.564,24.841], color="red blue", name="Arrows_27.2840003967_9")

cmd.load_cgo(cluster_dict["27.2840003967"], "Features_27.2840003967", 1)
cmd.load_cgo(cluster_dict["27.2840003967_arrows"], "Arrows_27.2840003967")
cmd.set("transparency", 0.2,"Features_27.2840003967")
cmd.group("Pharmacophore_27.2840003967", members="Features_27.2840003967")
cmd.group("Pharmacophore_27.2840003967", members="Arrows_27.2840003967")

if dirpath:
    f = join(dirpath, "label_threshold_27.2840003967.mol2")
else:
    f = "label_threshold_27.2840003967.mol2"

cmd.load(f, 'label_threshold_27.2840003967')
cmd.hide('everything', 'label_threshold_27.2840003967')
cmd.label("label_threshold_27.2840003967", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_27.2840003967', members= 'label_threshold_27.2840003967')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
