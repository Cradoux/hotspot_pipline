
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


cluster_dict = {"16.9120006561":[], "16.9120006561_arrows":[]}

cluster_dict["16.9120006561"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(14.5), float(72.0), float(25.0), float(1.0)]

cluster_dict["16.9120006561_arrows"] += cgo_arrow([14.5,72.0,25.0], [13.161,70.112,23.071], color="blue red", name="Arrows_16.9120006561_1")

cluster_dict["16.9120006561"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(70.5), float(26.5), float(1.0)]

cluster_dict["16.9120006561_arrows"] += cgo_arrow([17.5,70.5,26.5], [19.46,68.163,25.47], color="blue red", name="Arrows_16.9120006561_2")

cluster_dict["16.9120006561"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(76.5), float(28.0), float(1.0)]

cluster_dict["16.9120006561_arrows"] += cgo_arrow([18.0,76.5,28.0], [16.088,77.26,29.666], color="blue red", name="Arrows_16.9120006561_3")

cluster_dict["16.9120006561"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.5), float(72.5), float(24.0), float(1.0)]

cluster_dict["16.9120006561_arrows"] += cgo_arrow([19.5,72.5,24.0], [22.004,71.45,22.776], color="blue red", name="Arrows_16.9120006561_4")

cluster_dict["16.9120006561"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(76.5), float(18.0), float(1.0)]

cluster_dict["16.9120006561_arrows"] += cgo_arrow([22.5,76.5,18.0], [24.525,75.162,17.929], color="blue red", name="Arrows_16.9120006561_5")

cluster_dict["16.9120006561"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(18.5411248797), float(74.4497145988), float(25.472341983), float(1.0)]


cluster_dict["16.9120006561"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(17.290338783), float(69.5462763018), float(26.3262962662), float(1.0)]


cluster_dict["16.9120006561"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(14.5), float(74.5), float(26.0), float(1.0)]

cluster_dict["16.9120006561_arrows"] += cgo_arrow([14.5,74.5,26.0], [12.602,75.566,23.891], color="red blue", name="Arrows_16.9120006561_6")

cluster_dict["16.9120006561"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(73.5), float(27.0), float(1.0)]


cluster_dict["16.9120006561"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(77.0), float(23.5), float(1.0)]

cluster_dict["16.9120006561_arrows"] += cgo_arrow([18.0,77.0,23.5], [18.608,79.736,23.169], color="red blue", name="Arrows_16.9120006561_7")

cluster_dict["16.9120006561"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(76.5), float(27.5), float(1.0)]

cluster_dict["16.9120006561_arrows"] += cgo_arrow([21.0,76.5,27.5], [22.441,80.537,26.502], color="red blue", name="Arrows_16.9120006561_8")

cmd.load_cgo(cluster_dict["16.9120006561"], "Features_16.9120006561", 1)
cmd.load_cgo(cluster_dict["16.9120006561_arrows"], "Arrows_16.9120006561")
cmd.set("transparency", 0.2,"Features_16.9120006561")
cmd.group("Pharmacophore_16.9120006561", members="Features_16.9120006561")
cmd.group("Pharmacophore_16.9120006561", members="Arrows_16.9120006561")

if dirpath:
    f = join(dirpath, "label_threshold_16.9120006561.mol2")
else:
    f = "label_threshold_16.9120006561.mol2"

cmd.load(f, 'label_threshold_16.9120006561')
cmd.hide('everything', 'label_threshold_16.9120006561')
cmd.label("label_threshold_16.9120006561", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.9120006561', members= 'label_threshold_16.9120006561')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")