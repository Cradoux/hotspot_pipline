
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


cluster_dict = {"15.390999794":[], "15.390999794_arrows":[]}

cluster_dict["15.390999794"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(59.0), float(70.5), float(68.5), float(1.0)]

cluster_dict["15.390999794_arrows"] += cgo_arrow([59.0,70.5,68.5], [57.234,70.807,66.464], color="blue red", name="Arrows_15.390999794_1")

cluster_dict["15.390999794"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(63.5), float(71.0), float(68.5), float(1.0)]

cluster_dict["15.390999794_arrows"] += cgo_arrow([63.5,71.0,68.5], [62.225,71.597,67.806], color="blue red", name="Arrows_15.390999794_2")

cluster_dict["15.390999794"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(63.5), float(71.0), float(68.5), float(1.0)]

cluster_dict["15.390999794_arrows"] += cgo_arrow([63.5,71.0,68.5], [62.225,71.597,67.806], color="blue red", name="Arrows_15.390999794_3")

cluster_dict["15.390999794"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(67.5), float(69.0), float(70.5), float(1.0)]

cluster_dict["15.390999794_arrows"] += cgo_arrow([67.5,69.0,70.5], [70.684,70.92,69.926], color="blue red", name="Arrows_15.390999794_4")

cluster_dict["15.390999794"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(61.2913582521), float(67.5524655884), float(70.2293669851), float(1.0)]


cluster_dict["15.390999794"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(67.8602740577), float(72.8483444539), float(68.7133151742), float(1.0)]


cluster_dict["15.390999794"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(68.0), float(74.0), float(66.5), float(1.0)]

cluster_dict["15.390999794_arrows"] += cgo_arrow([68.0,74.0,66.5], [64.836,73.977,66.465], color="red blue", name="Arrows_15.390999794_5")

cluster_dict["15.390999794"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(67.0), float(76.0), float(75.0), float(1.0)]

cluster_dict["15.390999794_arrows"] += cgo_arrow([67.0,76.0,75.0], [65.846,79.346,75.307], color="red blue", name="Arrows_15.390999794_6")

cluster_dict["15.390999794"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(70.0), float(72.0), float(70.5), float(1.0)]

cluster_dict["15.390999794_arrows"] += cgo_arrow([70.0,72.0,70.5], [70.684,70.92,69.926], color="red blue", name="Arrows_15.390999794_7")

cluster_dict["15.390999794"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(67.0), float(77.0), float(70.5), float(1.0)]

cluster_dict["15.390999794_arrows"] += cgo_arrow([67.0,77.0,70.5], [64.782,77.81,70.101], color="red blue", name="Arrows_15.390999794_8")

cluster_dict["15.390999794"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(68.0), float(74.0), float(66.5), float(1.0)]

cluster_dict["15.390999794_arrows"] += cgo_arrow([68.0,74.0,66.5], [64.836,73.977,66.465], color="red blue", name="Arrows_15.390999794_9")

cluster_dict["15.390999794"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(70.0), float(72.0), float(70.5), float(1.0)]

cluster_dict["15.390999794_arrows"] += cgo_arrow([70.0,72.0,70.5], [70.684,70.92,69.926], color="red blue", name="Arrows_15.390999794_10")

cmd.load_cgo(cluster_dict["15.390999794"], "Features_15.390999794", 1)
cmd.load_cgo(cluster_dict["15.390999794_arrows"], "Arrows_15.390999794")
cmd.set("transparency", 0.2,"Features_15.390999794")
cmd.group("Pharmacophore_15.390999794", members="Features_15.390999794")
cmd.group("Pharmacophore_15.390999794", members="Arrows_15.390999794")

if dirpath:
    f = join(dirpath, "label_threshold_15.390999794.mol2")
else:
    f = "label_threshold_15.390999794.mol2"

cmd.load(f, 'label_threshold_15.390999794')
cmd.hide('everything', 'label_threshold_15.390999794')
cmd.label("label_threshold_15.390999794", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.390999794', members= 'label_threshold_15.390999794')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
