
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


cluster_dict = {"17.3129997253":[], "17.3129997253_arrows":[]}

cluster_dict["17.3129997253"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(75.0), float(24.5), float(1.0)]

cluster_dict["17.3129997253_arrows"] += cgo_arrow([18.0,75.0,24.5], [15.611,73.163,23.589], color="blue red", name="Arrows_17.3129997253_1")

cluster_dict["17.3129997253"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.0), float(67.5), float(17.0), float(1.0)]

cluster_dict["17.3129997253_arrows"] += cgo_arrow([19.0,67.5,17.0], [17.176,65.615,17.937], color="blue red", name="Arrows_17.3129997253_2")

cluster_dict["17.3129997253"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(19.2485320729), float(72.4462935829), float(19.8633768715), float(1.0)]


cluster_dict["17.3129997253"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(18.4027445746), float(67.6065001765), float(19.1065001765), float(1.0)]


cluster_dict["17.3129997253"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(73.0), float(25.5), float(1.0)]


cluster_dict["17.3129997253"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.5), float(76.5), float(24.5), float(1.0)]

cluster_dict["17.3129997253_arrows"] += cgo_arrow([16.5,76.5,24.5], [18.025,77.713,26.505], color="red blue", name="Arrows_17.3129997253_3")

cluster_dict["17.3129997253"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(68.0), float(15.5), float(1.0)]

cluster_dict["17.3129997253_arrows"] += cgo_arrow([17.5,68.0,15.5], [16.119,65.678,15.964], color="red blue", name="Arrows_17.3129997253_4")

cluster_dict["17.3129997253"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.0), float(70.0), float(22.5), float(1.0)]

cluster_dict["17.3129997253_arrows"] += cgo_arrow([17.0,70.0,22.5], [14.032,67.68,20.514], color="red blue", name="Arrows_17.3129997253_5")

cluster_dict["17.3129997253"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(71.5), float(24.0), float(1.0)]

cluster_dict["17.3129997253_arrows"] += cgo_arrow([17.5,71.5,24.0], [20.683,72.995,26.534], color="red blue", name="Arrows_17.3129997253_6")

cluster_dict["17.3129997253"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(74.0), float(16.5), float(1.0)]

cluster_dict["17.3129997253_arrows"] += cgo_arrow([18.5,74.0,16.5], [20.035,77.898,13.86], color="red blue", name="Arrows_17.3129997253_7")

cluster_dict["17.3129997253"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(19.5), float(72.0), float(17.0), float(1.0)]


cluster_dict["17.3129997253"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.5), float(73.0), float(19.0), float(1.0)]


cmd.load_cgo(cluster_dict["17.3129997253"], "Features_17.3129997253", 1)
cmd.load_cgo(cluster_dict["17.3129997253_arrows"], "Arrows_17.3129997253")
cmd.set("transparency", 0.2,"Features_17.3129997253")
cmd.group("Pharmacophore_17.3129997253", members="Features_17.3129997253")
cmd.group("Pharmacophore_17.3129997253", members="Arrows_17.3129997253")

if dirpath:
    f = join(dirpath, "label_threshold_17.3129997253.mol2")
else:
    f = "label_threshold_17.3129997253.mol2"

cmd.load(f, 'label_threshold_17.3129997253')
cmd.hide('everything', 'label_threshold_17.3129997253')
cmd.label("label_threshold_17.3129997253", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.3129997253', members= 'label_threshold_17.3129997253')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
