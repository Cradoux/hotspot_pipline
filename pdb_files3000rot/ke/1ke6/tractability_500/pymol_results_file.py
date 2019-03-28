
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


cluster_dict = {"16.9829998016":[], "16.9829998016_arrows":[]}

cluster_dict["16.9829998016"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(12.0), float(10.5), float(1.0)]

cluster_dict["16.9829998016_arrows"] += cgo_arrow([26.5,12.0,10.5], [24.487,14.204,8.555], color="blue red", name="Arrows_16.9829998016_1")

cluster_dict["16.9829998016"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(13.0), float(10.5), float(1.0)]

cluster_dict["16.9829998016_arrows"] += cgo_arrow([26.0,13.0,10.5], [24.487,14.204,8.555], color="blue red", name="Arrows_16.9829998016_2")

cluster_dict["16.9829998016"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(16.0), float(14.0), float(1.0)]

cluster_dict["16.9829998016_arrows"] += cgo_arrow([27.0,16.0,14.0], [25.103,18.159,13.412], color="blue red", name="Arrows_16.9829998016_3")

cluster_dict["16.9829998016"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(12.0), float(10.5), float(1.0)]

cluster_dict["16.9829998016_arrows"] += cgo_arrow([26.5,12.0,10.5], [24.487,14.204,8.555], color="blue red", name="Arrows_16.9829998016_4")

cluster_dict["16.9829998016"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(10.5), float(16.5), float(1.0)]

cluster_dict["16.9829998016_arrows"] += cgo_arrow([30.5,10.5,16.5], [33.465,11.952,16.935], color="blue red", name="Arrows_16.9829998016_5")

cluster_dict["16.9829998016"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(29.0260278639), float(10.9995023525), float(14.0305309133), float(1.0)]


cluster_dict["16.9829998016"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(25.5), float(14.5), float(12.0), float(1.0)]

cluster_dict["16.9829998016_arrows"] += cgo_arrow([25.5,14.5,12.0], [24.278,16.63,9.644], color="red blue", name="Arrows_16.9829998016_6")

cluster_dict["16.9829998016"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(28.5), float(10.5), float(14.5), float(1.0)]


cluster_dict["16.9829998016"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(31.5), float(11.5), float(13.5), float(1.0)]

cluster_dict["16.9829998016_arrows"] += cgo_arrow([31.5,11.5,13.5], [34.329,10.687,14.006], color="red blue", name="Arrows_16.9829998016_7")

cluster_dict["16.9829998016"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(16.5), float(17.0), float(1.0)]

cluster_dict["16.9829998016_arrows"] += cgo_arrow([30.5,16.5,17.0], [33.505,15.926,16.422], color="red blue", name="Arrows_16.9829998016_8")

cmd.load_cgo(cluster_dict["16.9829998016"], "Features_16.9829998016", 1)
cmd.load_cgo(cluster_dict["16.9829998016_arrows"], "Arrows_16.9829998016")
cmd.set("transparency", 0.2,"Features_16.9829998016")
cmd.group("Pharmacophore_16.9829998016", members="Features_16.9829998016")
cmd.group("Pharmacophore_16.9829998016", members="Arrows_16.9829998016")

if dirpath:
    f = join(dirpath, "label_threshold_16.9829998016.mol2")
else:
    f = "label_threshold_16.9829998016.mol2"

cmd.load(f, 'label_threshold_16.9829998016')
cmd.hide('everything', 'label_threshold_16.9829998016')
cmd.label("label_threshold_16.9829998016", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.9829998016', members= 'label_threshold_16.9829998016')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
