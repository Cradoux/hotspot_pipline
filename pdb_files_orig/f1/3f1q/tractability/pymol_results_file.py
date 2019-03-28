
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
gfiles = ['apolar.grd', 'acceptor.grd']
grids = ['apolar', 'acceptor']
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


cluster_dict = {"15.4790000916":[], "15.4790000916_arrows":[]}

cluster_dict["15.4790000916"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(6.47108653988), float(-35.0555815055), float(-1.68457325995), float(1.0)]


cluster_dict["15.4790000916"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(9.4463452259), float(-41.3062615921), float(2.4127503898), float(1.0)]


cluster_dict["15.4790000916"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(0.0), float(-39.5), float(-2.5), float(1.0)]

cluster_dict["15.4790000916_arrows"] += cgo_arrow([0.0,-39.5,-2.5], [-0.513,-35.685,-3.003], color="red blue", name="Arrows_15.4790000916_1")

cluster_dict["15.4790000916"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(2.0), float(-37.5), float(-5.5), float(1.0)]

cluster_dict["15.4790000916_arrows"] += cgo_arrow([2.0,-37.5,-5.5], [1.834,-38.561,-7.841], color="red blue", name="Arrows_15.4790000916_2")

cluster_dict["15.4790000916"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(4.0), float(-37.5), float(0.0), float(1.0)]

cluster_dict["15.4790000916_arrows"] += cgo_arrow([4.0,-37.5,0.0], [2.326,-39.824,2.141], color="red blue", name="Arrows_15.4790000916_3")

cluster_dict["15.4790000916"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(4.0), float(-35.0), float(-1.0), float(1.0)]

cluster_dict["15.4790000916_arrows"] += cgo_arrow([4.0,-35.0,-1.0], [0.932,-33.395,-1.451], color="red blue", name="Arrows_15.4790000916_4")

cluster_dict["15.4790000916"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(4.5), float(-34.5), float(-3.0), float(1.0)]

cluster_dict["15.4790000916_arrows"] += cgo_arrow([4.5,-34.5,-3.0], [0.932,-33.395,-1.451], color="red blue", name="Arrows_15.4790000916_5")

cluster_dict["15.4790000916"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(9.0), float(-32.0), float(-5.0), float(1.0)]


cmd.load_cgo(cluster_dict["15.4790000916"], "Features_15.4790000916", 1)
cmd.load_cgo(cluster_dict["15.4790000916_arrows"], "Arrows_15.4790000916")
cmd.set("transparency", 0.2,"Features_15.4790000916")
cmd.group("Pharmacophore_15.4790000916", members="Features_15.4790000916")
cmd.group("Pharmacophore_15.4790000916", members="Arrows_15.4790000916")

if dirpath:
    f = join(dirpath, "label_threshold_15.4790000916.mol2")
else:
    f = "label_threshold_15.4790000916.mol2"

cmd.load(f, 'label_threshold_15.4790000916')
cmd.hide('everything', 'label_threshold_15.4790000916')
cmd.label("label_threshold_15.4790000916", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.4790000916', members= 'label_threshold_15.4790000916')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
