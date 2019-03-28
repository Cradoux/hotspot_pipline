
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


cluster_dict = {"16.9699993134":[], "16.9699993134_arrows":[]}

cluster_dict["16.9699993134"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(67.6358537645), float(11.0504850752), float(18.0994987027), float(1.0)]


cluster_dict["16.9699993134"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(69.1781314046), float(15.5675692403), float(9.3705601519), float(1.0)]


cluster_dict["16.9699993134"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(64.5), float(11.0), float(18.5), float(1.0)]

cluster_dict["16.9699993134_arrows"] += cgo_arrow([64.5,11.0,18.5], [62.285,12.705,18.43], color="red blue", name="Arrows_16.9699993134_1")

cluster_dict["16.9699993134"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(66.0), float(8.5), float(18.5), float(1.0)]


cluster_dict["16.9699993134"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(66.5), float(11.0), float(21.0), float(1.0)]

cluster_dict["16.9699993134_arrows"] += cgo_arrow([66.5,11.0,21.0], [63.534,10.951,22.284], color="red blue", name="Arrows_16.9699993134_2")

cluster_dict["16.9699993134"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(65.5), float(15.5), float(8.5), float(1.0)]

cluster_dict["16.9699993134_arrows"] += cgo_arrow([65.5,15.5,8.5], [65.199,17.927,6.197], color="red blue", name="Arrows_16.9699993134_3")

cluster_dict["16.9699993134"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(69.5), float(12.5), float(14.5), float(1.0)]

cluster_dict["16.9699993134_arrows"] += cgo_arrow([69.5,12.5,14.5], [71.288,11.805,14.01], color="red blue", name="Arrows_16.9699993134_4")

cluster_dict["16.9699993134"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(67.0), float(15.0), float(9.0), float(1.0)]

cluster_dict["16.9699993134_arrows"] += cgo_arrow([67.0,15.0,9.0], [65.199,17.927,6.197], color="red blue", name="Arrows_16.9699993134_5")

cluster_dict["16.9699993134"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(68.0), float(17.0), float(6.0), float(1.0)]

cluster_dict["16.9699993134_arrows"] += cgo_arrow([68.0,17.0,6.0], [65.199,17.927,6.197], color="red blue", name="Arrows_16.9699993134_6")

cluster_dict["16.9699993134"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(69.5), float(13.5), float(13.5), float(1.0)]

cluster_dict["16.9699993134_arrows"] += cgo_arrow([69.5,13.5,13.5], [71.288,11.805,14.01], color="red blue", name="Arrows_16.9699993134_7")

cluster_dict["16.9699993134"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(72.0), float(15.5), float(8.5), float(1.0)]

cluster_dict["16.9699993134_arrows"] += cgo_arrow([72.0,15.5,8.5], [73.123,11.343,7.892], color="red blue", name="Arrows_16.9699993134_8")

cmd.load_cgo(cluster_dict["16.9699993134"], "Features_16.9699993134", 1)
cmd.load_cgo(cluster_dict["16.9699993134_arrows"], "Arrows_16.9699993134")
cmd.set("transparency", 0.2,"Features_16.9699993134")
cmd.group("Pharmacophore_16.9699993134", members="Features_16.9699993134")
cmd.group("Pharmacophore_16.9699993134", members="Arrows_16.9699993134")

if dirpath:
    f = join(dirpath, "label_threshold_16.9699993134.mol2")
else:
    f = "label_threshold_16.9699993134.mol2"

cmd.load(f, 'label_threshold_16.9699993134')
cmd.hide('everything', 'label_threshold_16.9699993134')
cmd.label("label_threshold_16.9699993134", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.9699993134', members= 'label_threshold_16.9699993134')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
