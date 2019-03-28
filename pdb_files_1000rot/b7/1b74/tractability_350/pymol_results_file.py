
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


cluster_dict = {"14.0010004044":[], "14.0010004044_arrows":[]}

cluster_dict["14.0010004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(37.0), float(31.5), float(1.0)]

cluster_dict["14.0010004044_arrows"] += cgo_arrow([21.0,37.0,31.5], [24.062,36.854,31.532], color="blue red", name="Arrows_14.0010004044_1")

cluster_dict["14.0010004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(27.0), float(24.5), float(1.0)]

cluster_dict["14.0010004044_arrows"] += cgo_arrow([21.0,27.0,24.5], [20.72,28.331,26.689], color="blue red", name="Arrows_14.0010004044_2")

cluster_dict["14.0010004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(34.5), float(27.5), float(1.0)]

cluster_dict["14.0010004044_arrows"] += cgo_arrow([24.5,34.5,27.5], [26.805,36.028,27.954], color="blue red", name="Arrows_14.0010004044_3")

cluster_dict["14.0010004044"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(19.521953446), float(39.3425147302), float(29.9665187275), float(1.0)]


cluster_dict["14.0010004044"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.6933600548), float(29.1925600959), float(23.1096163426), float(1.0)]


cluster_dict["14.0010004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(35.0), float(25.5), float(1.0)]

cluster_dict["14.0010004044_arrows"] += cgo_arrow([17.5,35.0,25.5], [15.614,33.4,25.626], color="red blue", name="Arrows_14.0010004044_4")

cluster_dict["14.0010004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(41.0), float(26.5), float(1.0)]

cluster_dict["14.0010004044_arrows"] += cgo_arrow([18.5,41.0,26.5], [20.521,43.214,26.409], color="red blue", name="Arrows_14.0010004044_5")

cluster_dict["14.0010004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.0), float(27.0), float(24.5), float(1.0)]

cluster_dict["14.0010004044_arrows"] += cgo_arrow([20.0,27.0,24.5], [20.72,28.331,26.689], color="red blue", name="Arrows_14.0010004044_6")

cluster_dict["14.0010004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(30.5), float(26.0), float(1.0)]

cluster_dict["14.0010004044_arrows"] += cgo_arrow([26.0,30.5,26.0], [27.697,33.03,25.727], color="red blue", name="Arrows_14.0010004044_7")

cluster_dict["14.0010004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(34.0), float(29.0), float(1.0)]

cluster_dict["14.0010004044_arrows"] += cgo_arrow([22.5,34.0,29.0], [23.801,33.08,31.81], color="red blue", name="Arrows_14.0010004044_8")

cluster_dict["14.0010004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(34.5), float(23.5), float(1.0)]

cluster_dict["14.0010004044_arrows"] += cgo_arrow([26.0,34.5,23.5], [28.493,36.168,24.374], color="red blue", name="Arrows_14.0010004044_9")

cmd.load_cgo(cluster_dict["14.0010004044"], "Features_14.0010004044", 1)
cmd.load_cgo(cluster_dict["14.0010004044_arrows"], "Arrows_14.0010004044")
cmd.set("transparency", 0.2,"Features_14.0010004044")
cmd.group("Pharmacophore_14.0010004044", members="Features_14.0010004044")
cmd.group("Pharmacophore_14.0010004044", members="Arrows_14.0010004044")

if dirpath:
    f = join(dirpath, "label_threshold_14.0010004044.mol2")
else:
    f = "label_threshold_14.0010004044.mol2"

cmd.load(f, 'label_threshold_14.0010004044')
cmd.hide('everything', 'label_threshold_14.0010004044')
cmd.label("label_threshold_14.0010004044", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.0010004044', members= 'label_threshold_14.0010004044')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
