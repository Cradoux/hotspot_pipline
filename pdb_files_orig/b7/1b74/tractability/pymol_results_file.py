
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


cluster_dict = {"13.7729997635":[], "13.7729997635_arrows":[]}

cluster_dict["13.7729997635"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(36.5), float(31.5), float(1.0)]

cluster_dict["13.7729997635_arrows"] += cgo_arrow([21.0,36.5,31.5], [20.448,33.722,32.582], color="blue red", name="Arrows_13.7729997635_1")

cluster_dict["13.7729997635"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(28.0), float(24.5), float(1.0)]

cluster_dict["13.7729997635_arrows"] += cgo_arrow([21.0,28.0,24.5], [20.72,28.331,26.689], color="blue red", name="Arrows_13.7729997635_2")

cluster_dict["13.7729997635"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(34.5), float(27.5), float(1.0)]

cluster_dict["13.7729997635_arrows"] += cgo_arrow([24.5,34.5,27.5], [26.805,36.028,27.954], color="blue red", name="Arrows_13.7729997635_3")

cluster_dict["13.7729997635"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(20.0562328902), float(37.0862625828), float(27.9539665172), float(1.0)]


cluster_dict["13.7729997635"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.7458981567), float(26.8150628555), float(20.6390580393), float(1.0)]


cluster_dict["13.7729997635"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(41.0), float(26.5), float(1.0)]

cluster_dict["13.7729997635_arrows"] += cgo_arrow([18.5,41.0,26.5], [20.521,43.214,26.409], color="red blue", name="Arrows_13.7729997635_4")

cluster_dict["13.7729997635"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(30.5), float(26.0), float(1.0)]

cluster_dict["13.7729997635_arrows"] += cgo_arrow([26.0,30.5,26.0], [27.697,33.03,25.727], color="red blue", name="Arrows_13.7729997635_5")

cluster_dict["13.7729997635"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(34.0), float(29.0), float(1.0)]

cluster_dict["13.7729997635_arrows"] += cgo_arrow([22.5,34.0,29.0], [23.801,33.08,31.81], color="red blue", name="Arrows_13.7729997635_6")

cluster_dict["13.7729997635"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(34.5), float(23.5), float(1.0)]

cluster_dict["13.7729997635_arrows"] += cgo_arrow([26.5,34.5,23.5], [28.493,36.168,24.374], color="red blue", name="Arrows_13.7729997635_7")

cmd.load_cgo(cluster_dict["13.7729997635"], "Features_13.7729997635", 1)
cmd.load_cgo(cluster_dict["13.7729997635_arrows"], "Arrows_13.7729997635")
cmd.set("transparency", 0.2,"Features_13.7729997635")
cmd.group("Pharmacophore_13.7729997635", members="Features_13.7729997635")
cmd.group("Pharmacophore_13.7729997635", members="Arrows_13.7729997635")

if dirpath:
    f = join(dirpath, "label_threshold_13.7729997635.mol2")
else:
    f = "label_threshold_13.7729997635.mol2"

cmd.load(f, 'label_threshold_13.7729997635')
cmd.hide('everything', 'label_threshold_13.7729997635')
cmd.label("label_threshold_13.7729997635", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.7729997635', members= 'label_threshold_13.7729997635')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
