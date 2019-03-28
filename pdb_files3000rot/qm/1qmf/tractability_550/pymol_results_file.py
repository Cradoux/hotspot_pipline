
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


cluster_dict = {"15.9259996414":[], "15.9259996414_arrows":[]}

cluster_dict["15.9259996414"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(97.0), float(58.0), float(51.5), float(1.0)]

cluster_dict["15.9259996414_arrows"] += cgo_arrow([97.0,58.0,51.5], [96.892,59.428,54.314], color="blue red", name="Arrows_15.9259996414_1")

cluster_dict["15.9259996414"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(98.5), float(64.5), float(52.5), float(1.0)]

cluster_dict["15.9259996414_arrows"] += cgo_arrow([98.5,64.5,52.5], [100.719,65.589,53.902], color="blue red", name="Arrows_15.9259996414_2")

cluster_dict["15.9259996414"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(102.5), float(53.0), float(53.0), float(1.0)]

cluster_dict["15.9259996414_arrows"] += cgo_arrow([102.5,53.0,53.0], [100.549,51.341,51.975], color="blue red", name="Arrows_15.9259996414_3")

cluster_dict["15.9259996414"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(103.0), float(57.0), float(53.0), float(1.0)]

cluster_dict["15.9259996414_arrows"] += cgo_arrow([103.0,57.0,53.0], [104.159,54.788,54.922], color="blue red", name="Arrows_15.9259996414_4")

cluster_dict["15.9259996414"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(100.860133035), float(58.9139409571), float(50.0001763024), float(1.0)]


cluster_dict["15.9259996414"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(102.75), float(53.0), float(46.0), float(1.0)]


cluster_dict["15.9259996414"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(97.0), float(65.0), float(55.5), float(1.0)]

cluster_dict["15.9259996414_arrows"] += cgo_arrow([97.0,65.0,55.5], [99.525,64.248,56.728], color="red blue", name="Arrows_15.9259996414_5")

cluster_dict["15.9259996414"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(98.5), float(54.0), float(48.0), float(1.0)]

cluster_dict["15.9259996414_arrows"] += cgo_arrow([98.5,54.0,48.0], [97.986,50.759,47.681], color="red blue", name="Arrows_15.9259996414_6")

cluster_dict["15.9259996414"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(98.5), float(65.5), float(49.0), float(1.0)]

cluster_dict["15.9259996414_arrows"] += cgo_arrow([98.5,65.5,49.0], [98.672,66.443,45.392], color="red blue", name="Arrows_15.9259996414_7")

cluster_dict["15.9259996414"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(102.0), float(66.0), float(49.0), float(1.0)]

cluster_dict["15.9259996414_arrows"] += cgo_arrow([102.0,66.0,49.0], [101.484,68.93,48.907], color="red blue", name="Arrows_15.9259996414_8")

cluster_dict["15.9259996414"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(102.5), float(60.5), float(49.0), float(1.0)]

cluster_dict["15.9259996414_arrows"] += cgo_arrow([102.5,60.5,49.0], [100.246,57.359,46.359], color="red blue", name="Arrows_15.9259996414_9")

cmd.load_cgo(cluster_dict["15.9259996414"], "Features_15.9259996414", 1)
cmd.load_cgo(cluster_dict["15.9259996414_arrows"], "Arrows_15.9259996414")
cmd.set("transparency", 0.2,"Features_15.9259996414")
cmd.group("Pharmacophore_15.9259996414", members="Features_15.9259996414")
cmd.group("Pharmacophore_15.9259996414", members="Arrows_15.9259996414")

if dirpath:
    f = join(dirpath, "label_threshold_15.9259996414.mol2")
else:
    f = "label_threshold_15.9259996414.mol2"

cmd.load(f, 'label_threshold_15.9259996414')
cmd.hide('everything', 'label_threshold_15.9259996414')
cmd.label("label_threshold_15.9259996414", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.9259996414', members= 'label_threshold_15.9259996414')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
