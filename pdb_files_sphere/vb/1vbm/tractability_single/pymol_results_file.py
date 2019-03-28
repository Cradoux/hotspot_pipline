
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


cluster_dict = {"14.2760000229":[], "14.2760000229_arrows":[]}

cluster_dict["14.2760000229"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(0.0), float(176.0), float(244.0), float(1.0)]

cluster_dict["14.2760000229_arrows"] += cgo_arrow([0.0,176.0,244.0], [-1.528,175.919,246.782], color="blue red", name="Arrows_14.2760000229_1")

cluster_dict["14.2760000229"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(0.5), float(169.5), float(248.5), float(1.0)]

cluster_dict["14.2760000229_arrows"] += cgo_arrow([0.5,169.5,248.5], [0.757,167.114,246.986], color="blue red", name="Arrows_14.2760000229_2")

cluster_dict["14.2760000229"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(6.5), float(172.0), float(243.0), float(1.0)]

cluster_dict["14.2760000229_arrows"] += cgo_arrow([6.5,172.0,243.0], [6.519,172.682,240.438], color="blue red", name="Arrows_14.2760000229_3")

cluster_dict["14.2760000229"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(7.0), float(169.5), float(245.0), float(1.0)]

cluster_dict["14.2760000229_arrows"] += cgo_arrow([7.0,169.5,245.0], [5.69,168.686,244.539], color="blue red", name="Arrows_14.2760000229_4")

cluster_dict["14.2760000229"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(10.5), float(173.0), float(244.0), float(1.0)]

cluster_dict["14.2760000229_arrows"] += cgo_arrow([10.5,173.0,244.0], [13.04,173.528,242.252], color="blue red", name="Arrows_14.2760000229_5")

cluster_dict["14.2760000229"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(3.09148812057), float(171.828665413), float(246.128708757), float(1.0)]


cluster_dict["14.2760000229"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(8.03858012344), float(168.561925606), float(246.228952095), float(1.0)]


cluster_dict["14.2760000229"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-0.5), float(171.5), float(250.5), float(1.0)]

cluster_dict["14.2760000229_arrows"] += cgo_arrow([-0.5,171.5,250.5], [-1.924,175.146,250.952], color="red blue", name="Arrows_14.2760000229_6")

cluster_dict["14.2760000229"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(8.0), float(173.0), float(242.5), float(1.0)]

cluster_dict["14.2760000229_arrows"] += cgo_arrow([8.0,173.0,242.5], [6.519,172.682,240.438], color="red blue", name="Arrows_14.2760000229_7")

cluster_dict["14.2760000229"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(8.0), float(171.0), float(247.5), float(1.0)]

cluster_dict["14.2760000229_arrows"] += cgo_arrow([8.0,171.0,247.5], [8.705,173.68,248.6], color="red blue", name="Arrows_14.2760000229_8")

cluster_dict["14.2760000229"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(11.0), float(172.0), float(244.5), float(1.0)]

cluster_dict["14.2760000229_arrows"] += cgo_arrow([11.0,172.0,244.5], [13.04,173.528,242.252], color="red blue", name="Arrows_14.2760000229_9")

cluster_dict["14.2760000229"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(11.0), float(172.0), float(244.5), float(1.0)]

cluster_dict["14.2760000229_arrows"] += cgo_arrow([11.0,172.0,244.5], [13.04,173.528,242.252], color="red blue", name="Arrows_14.2760000229_10")

cmd.load_cgo(cluster_dict["14.2760000229"], "Features_14.2760000229", 1)
cmd.load_cgo(cluster_dict["14.2760000229_arrows"], "Arrows_14.2760000229")
cmd.set("transparency", 0.2,"Features_14.2760000229")
cmd.group("Pharmacophore_14.2760000229", members="Features_14.2760000229")
cmd.group("Pharmacophore_14.2760000229", members="Arrows_14.2760000229")

if dirpath:
    f = join(dirpath, "label_threshold_14.2760000229.mol2")
else:
    f = "label_threshold_14.2760000229.mol2"

cmd.load(f, 'label_threshold_14.2760000229')
cmd.hide('everything', 'label_threshold_14.2760000229')
cmd.label("label_threshold_14.2760000229", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.2760000229', members= 'label_threshold_14.2760000229')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
