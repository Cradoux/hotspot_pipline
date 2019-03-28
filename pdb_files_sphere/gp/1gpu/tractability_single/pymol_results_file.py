
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


cluster_dict = {"18.3530006409":[], "18.3530006409_arrows":[]}

cluster_dict["18.3530006409"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(74.0), float(28.0), float(1.0)]

cluster_dict["18.3530006409_arrows"] += cgo_arrow([18.0,74.0,28.0], [16.088,77.26,29.666], color="blue red", name="Arrows_18.3530006409_1")

cluster_dict["18.3530006409"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(74.0), float(28.0), float(1.0)]

cluster_dict["18.3530006409_arrows"] += cgo_arrow([18.0,74.0,28.0], [16.088,77.26,29.666], color="blue red", name="Arrows_18.3530006409_2")

cluster_dict["18.3530006409"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.5), float(75.5), float(22.0), float(1.0)]

cluster_dict["18.3530006409_arrows"] += cgo_arrow([19.5,75.5,22.0], [20.418,79.948,22.04], color="blue red", name="Arrows_18.3530006409_3")

cluster_dict["18.3530006409"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.5), float(72.5), float(22.5), float(1.0)]

cluster_dict["18.3530006409_arrows"] += cgo_arrow([19.5,72.5,22.5], [22.004,71.45,22.776], color="blue red", name="Arrows_18.3530006409_4")

cluster_dict["18.3530006409"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(18.6944988438), float(74.3076920031), float(25.3197305358), float(1.0)]


cluster_dict["18.3530006409"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(17.4296459048), float(69.060698831), float(25.8899540611), float(1.0)]


cluster_dict["18.3530006409"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(18.990645823), float(69.6253939323), float(24.9559636482), float(1.0)]


cluster_dict["18.3530006409"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(14.0), float(74.5), float(26.0), float(1.0)]

cluster_dict["18.3530006409_arrows"] += cgo_arrow([14.0,74.5,26.0], [12.602,75.566,23.891], color="red blue", name="Arrows_18.3530006409_5")

cluster_dict["18.3530006409"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(77.0), float(23.5), float(1.0)]

cluster_dict["18.3530006409_arrows"] += cgo_arrow([18.5,77.0,23.5], [18.608,79.736,23.169], color="red blue", name="Arrows_18.3530006409_6")

cluster_dict["18.3530006409"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(76.0), float(27.0), float(1.0)]

cluster_dict["18.3530006409_arrows"] += cgo_arrow([21.0,76.0,27.0], [22.441,80.537,26.502], color="red blue", name="Arrows_18.3530006409_7")

cluster_dict["18.3530006409"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.5), float(73.5), float(22.5), float(1.0)]

cluster_dict["18.3530006409_arrows"] += cgo_arrow([21.5,73.5,22.5], [22.004,71.45,22.776], color="red blue", name="Arrows_18.3530006409_8")

cmd.load_cgo(cluster_dict["18.3530006409"], "Features_18.3530006409", 1)
cmd.load_cgo(cluster_dict["18.3530006409_arrows"], "Arrows_18.3530006409")
cmd.set("transparency", 0.2,"Features_18.3530006409")
cmd.group("Pharmacophore_18.3530006409", members="Features_18.3530006409")
cmd.group("Pharmacophore_18.3530006409", members="Arrows_18.3530006409")

if dirpath:
    f = join(dirpath, "label_threshold_18.3530006409.mol2")
else:
    f = "label_threshold_18.3530006409.mol2"

cmd.load(f, 'label_threshold_18.3530006409')
cmd.hide('everything', 'label_threshold_18.3530006409')
cmd.label("label_threshold_18.3530006409", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_18.3530006409', members= 'label_threshold_18.3530006409')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
