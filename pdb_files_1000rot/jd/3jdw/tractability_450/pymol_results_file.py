
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


cluster_dict = {"12.0010004044":[], "12.0010004044_arrows":[]}

cluster_dict["12.0010004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(44.5), float(73.5), float(15.0), float(1.0)]

cluster_dict["12.0010004044_arrows"] += cgo_arrow([44.5,73.5,15.0], [44.9,70.913,13.25], color="blue red", name="Arrows_12.0010004044_1")

cluster_dict["12.0010004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(44.0), float(73.0), float(19.0), float(1.0)]

cluster_dict["12.0010004044_arrows"] += cgo_arrow([44.0,73.0,19.0], [42.775,75.507,19.924], color="blue red", name="Arrows_12.0010004044_2")

cluster_dict["12.0010004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(45.0), float(73.0), float(15.0), float(1.0)]

cluster_dict["12.0010004044_arrows"] += cgo_arrow([45.0,73.0,15.0], [44.9,70.913,13.25], color="blue red", name="Arrows_12.0010004044_3")

cluster_dict["12.0010004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(45.0), float(75.0), float(13.5), float(1.0)]

cluster_dict["12.0010004044_arrows"] += cgo_arrow([45.0,75.0,13.5], [42.811,75.821,11.69], color="blue red", name="Arrows_12.0010004044_4")

cluster_dict["12.0010004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(45.5), float(72.5), float(15.0), float(1.0)]

cluster_dict["12.0010004044_arrows"] += cgo_arrow([45.5,72.5,15.0], [44.9,70.913,13.25], color="blue red", name="Arrows_12.0010004044_5")

cluster_dict["12.0010004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(48.0), float(70.5), float(16.5), float(1.0)]

cluster_dict["12.0010004044_arrows"] += cgo_arrow([48.0,70.5,16.5], [49.527,70.997,19.052], color="blue red", name="Arrows_12.0010004044_6")

cluster_dict["12.0010004044"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(42.2959182624), float(80.5126810668), float(14.5806047189), float(1.0)]


cluster_dict["12.0010004044"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(45.5042560249), float(72.3664096762), float(16.2554346136), float(1.0)]


cluster_dict["12.0010004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(42.0), float(72.5), float(13.5), float(1.0)]

cluster_dict["12.0010004044_arrows"] += cgo_arrow([42.0,72.5,13.5], [43.61,72.345,11.062], color="red blue", name="Arrows_12.0010004044_7")

cluster_dict["12.0010004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(45.0), float(73.5), float(17.5), float(1.0)]

cluster_dict["12.0010004044_arrows"] += cgo_arrow([45.0,73.5,17.5], [45.656,76.341,19.264], color="red blue", name="Arrows_12.0010004044_8")

cluster_dict["12.0010004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(44.0), float(77.5), float(16.0), float(1.0)]

cluster_dict["12.0010004044_arrows"] += cgo_arrow([44.0,77.5,16.0], [43.632,79.741,17.764], color="red blue", name="Arrows_12.0010004044_9")

cmd.load_cgo(cluster_dict["12.0010004044"], "Features_12.0010004044", 1)
cmd.load_cgo(cluster_dict["12.0010004044_arrows"], "Arrows_12.0010004044")
cmd.set("transparency", 0.2,"Features_12.0010004044")
cmd.group("Pharmacophore_12.0010004044", members="Features_12.0010004044")
cmd.group("Pharmacophore_12.0010004044", members="Arrows_12.0010004044")

if dirpath:
    f = join(dirpath, "label_threshold_12.0010004044.mol2")
else:
    f = "label_threshold_12.0010004044.mol2"

cmd.load(f, 'label_threshold_12.0010004044')
cmd.hide('everything', 'label_threshold_12.0010004044')
cmd.label("label_threshold_12.0010004044", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.0010004044', members= 'label_threshold_12.0010004044')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
