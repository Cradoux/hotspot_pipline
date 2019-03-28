
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


cluster_dict = {"15.3219995499":[], "15.3219995499_arrows":[]}

cluster_dict["15.3219995499"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(98.5), float(9.0), float(22.0), float(1.0)]

cluster_dict["15.3219995499_arrows"] += cgo_arrow([98.5,9.0,22.0], [96.849,6.806,22.678], color="blue red", name="Arrows_15.3219995499_1")

cluster_dict["15.3219995499"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(99.5), float(11.0), float(16.0), float(1.0)]

cluster_dict["15.3219995499_arrows"] += cgo_arrow([99.5,11.0,16.0], [101.793,12.269,14.709], color="blue red", name="Arrows_15.3219995499_2")

cluster_dict["15.3219995499"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(103.5), float(9.5), float(19.5), float(1.0)]

cluster_dict["15.3219995499_arrows"] += cgo_arrow([103.5,9.5,19.5], [104.648,6.085,21.622], color="blue red", name="Arrows_15.3219995499_3")

cluster_dict["15.3219995499"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(109.5), float(9.5), float(20.0), float(1.0)]

cluster_dict["15.3219995499_arrows"] += cgo_arrow([109.5,9.5,20.0], [109.895,6.591,19.309], color="blue red", name="Arrows_15.3219995499_4")

cluster_dict["15.3219995499"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(100.307883631), float(11.396720045), float(20.1708377161), float(1.0)]


cluster_dict["15.3219995499"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(105.195528138), float(15.9617164997), float(24.1082899644), float(1.0)]


cluster_dict["15.3219995499"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(97.0), float(12.5), float(21.0), float(1.0)]

cluster_dict["15.3219995499_arrows"] += cgo_arrow([97.0,12.5,21.0], [94.861,13.295,22.615], color="red blue", name="Arrows_15.3219995499_5")

cluster_dict["15.3219995499"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(98.5), float(14.0), float(19.0), float(1.0)]

cluster_dict["15.3219995499_arrows"] += cgo_arrow([98.5,14.0,19.0], [102.548,16.345,18.531], color="red blue", name="Arrows_15.3219995499_6")

cluster_dict["15.3219995499"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(99.0), float(5.5), float(20.0), float(1.0)]

cluster_dict["15.3219995499_arrows"] += cgo_arrow([99.0,5.5,20.0], [96.932,3.526,18.588], color="red blue", name="Arrows_15.3219995499_7")

cluster_dict["15.3219995499"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(101.0), float(10.5), float(18.5), float(1.0)]

cluster_dict["15.3219995499_arrows"] += cgo_arrow([101.0,10.5,18.5], [103.102,11.6,16.39], color="red blue", name="Arrows_15.3219995499_8")

cluster_dict["15.3219995499"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(105.5), float(11.0), float(18.0), float(1.0)]

cluster_dict["15.3219995499_arrows"] += cgo_arrow([105.5,11.0,18.0], [103.102,11.6,16.39], color="red blue", name="Arrows_15.3219995499_9")

cmd.load_cgo(cluster_dict["15.3219995499"], "Features_15.3219995499", 1)
cmd.load_cgo(cluster_dict["15.3219995499_arrows"], "Arrows_15.3219995499")
cmd.set("transparency", 0.2,"Features_15.3219995499")
cmd.group("Pharmacophore_15.3219995499", members="Features_15.3219995499")
cmd.group("Pharmacophore_15.3219995499", members="Arrows_15.3219995499")

if dirpath:
    f = join(dirpath, "label_threshold_15.3219995499.mol2")
else:
    f = "label_threshold_15.3219995499.mol2"

cmd.load(f, 'label_threshold_15.3219995499')
cmd.hide('everything', 'label_threshold_15.3219995499')
cmd.label("label_threshold_15.3219995499", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.3219995499', members= 'label_threshold_15.3219995499')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
