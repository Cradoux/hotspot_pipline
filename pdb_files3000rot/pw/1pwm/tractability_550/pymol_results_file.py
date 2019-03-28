
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


cluster_dict = {"15.890999794":[], "15.890999794_arrows":[]}

cluster_dict["15.890999794"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(15.5), float(-4.5), float(19.5), float(1.0)]

cluster_dict["15.890999794_arrows"] += cgo_arrow([15.5,-4.5,19.5], [12.907,-4.46,18.205], color="blue red", name="Arrows_15.890999794_1")

cluster_dict["15.890999794"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(17.3174248665), float(-7.53904230351), float(17.71710574), float(1.0)]


cluster_dict["15.890999794"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(13.711530398), float(-3.07032697585), float(20.6429402961), float(1.0)]


cluster_dict["15.890999794"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(14.3553399628), float(-2.16012372109), float(18.4418551293), float(1.0)]


cluster_dict["15.890999794"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(-0.5), float(19.5), float(1.0)]


cluster_dict["15.890999794"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(14.0), float(-8.0), float(14.5), float(1.0)]

cluster_dict["15.890999794_arrows"] += cgo_arrow([14.0,-8.0,14.5], [14.389,-3.838,14.455], color="red blue", name="Arrows_15.890999794_2")

cluster_dict["15.890999794"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(-8.0), float(13.5), float(1.0)]

cluster_dict["15.890999794_arrows"] += cgo_arrow([18.5,-8.0,13.5], [20.712,-4.579,11.589], color="red blue", name="Arrows_15.890999794_3")

cluster_dict["15.890999794"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(-3.5), float(23.0), float(1.0)]

cluster_dict["15.890999794_arrows"] += cgo_arrow([15.0,-3.5,23.0], [13.814,-5.448,21.846], color="red blue", name="Arrows_15.890999794_4")

cluster_dict["15.890999794"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(-2.5), float(18.0), float(1.0)]

cluster_dict["15.890999794_arrows"] += cgo_arrow([16.0,-2.5,18.0], [15.247,-0.45,15.978], color="red blue", name="Arrows_15.890999794_5")

cluster_dict["15.890999794"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(-7.5), float(16.5), float(1.0)]

cluster_dict["15.890999794_arrows"] += cgo_arrow([18.0,-7.5,16.5], [20.681,-4.461,18.07], color="red blue", name="Arrows_15.890999794_6")

cluster_dict["15.890999794"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(19.0), float(-6.0), float(24.5), float(1.0)]

cluster_dict["15.890999794_arrows"] += cgo_arrow([19.0,-6.0,24.5], [20.616,-9.446,24.69], color="red blue", name="Arrows_15.890999794_7")

cluster_dict["15.890999794"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(-8.0), float(13.5), float(1.0)]

cluster_dict["15.890999794_arrows"] += cgo_arrow([18.5,-8.0,13.5], [20.712,-4.579,11.589], color="red blue", name="Arrows_15.890999794_8")

cluster_dict["15.890999794"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(-8.0), float(12.0), float(1.0)]

cluster_dict["15.890999794_arrows"] += cgo_arrow([21.0,-8.0,12.0], [21.313,-6.25,9.361], color="red blue", name="Arrows_15.890999794_9")

cmd.load_cgo(cluster_dict["15.890999794"], "Features_15.890999794", 1)
cmd.load_cgo(cluster_dict["15.890999794_arrows"], "Arrows_15.890999794")
cmd.set("transparency", 0.2,"Features_15.890999794")
cmd.group("Pharmacophore_15.890999794", members="Features_15.890999794")
cmd.group("Pharmacophore_15.890999794", members="Arrows_15.890999794")

if dirpath:
    f = join(dirpath, "label_threshold_15.890999794.mol2")
else:
    f = "label_threshold_15.890999794.mol2"

cmd.load(f, 'label_threshold_15.890999794')
cmd.hide('everything', 'label_threshold_15.890999794')
cmd.label("label_threshold_15.890999794", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.890999794', members= 'label_threshold_15.890999794')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
