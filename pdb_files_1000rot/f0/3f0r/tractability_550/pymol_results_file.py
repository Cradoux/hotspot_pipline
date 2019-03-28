
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


cluster_dict = {"12.6599998474":[], "12.6599998474_arrows":[]}

cluster_dict["12.6599998474"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(33.5), float(-6.5), float(6.5), float(1.0)]

cluster_dict["12.6599998474_arrows"] += cgo_arrow([33.5,-6.5,6.5], [31.839,-6.19,4.213], color="blue red", name="Arrows_12.6599998474_1")

cluster_dict["12.6599998474"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(33.0), float(-4.5), float(2.5), float(1.0)]

cluster_dict["12.6599998474_arrows"] += cgo_arrow([33.0,-4.5,2.5], [31.597,-4.137,0.858], color="blue red", name="Arrows_12.6599998474_2")

cluster_dict["12.6599998474"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(-4.0), float(2.0), float(1.0)]

cluster_dict["12.6599998474_arrows"] += cgo_arrow([36.5,-4.0,2.0], [37.036,-1.381,0.495], color="blue red", name="Arrows_12.6599998474_3")

cluster_dict["12.6599998474"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(39.5), float(-4.5), float(3.0), float(1.0)]

cluster_dict["12.6599998474_arrows"] += cgo_arrow([39.5,-4.5,3.0], [38.748,-6.913,1.552], color="blue red", name="Arrows_12.6599998474_4")

cluster_dict["12.6599998474"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(36.7748465598), float(-3.40858737549), float(5.31558303584), float(1.0)]


cluster_dict["12.6599998474"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(39.9971658274), float(-4.24419957772), float(13.8564767021), float(1.0)]


cluster_dict["12.6599998474"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(37.0), float(-6.0), float(8.0), float(1.0)]

cluster_dict["12.6599998474_arrows"] += cgo_arrow([37.0,-6.0,8.0], [39.742,-6.281,9.454], color="red blue", name="Arrows_12.6599998474_5")

cluster_dict["12.6599998474"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(40.5), float(-8.0), float(11.0), float(1.0)]

cluster_dict["12.6599998474_arrows"] += cgo_arrow([40.5,-8.0,11.0], [39.742,-6.281,9.454], color="red blue", name="Arrows_12.6599998474_6")

cluster_dict["12.6599998474"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(41.0), float(-3.0), float(3.5), float(1.0)]

cluster_dict["12.6599998474_arrows"] += cgo_arrow([41.0,-3.0,3.5], [43.985,-5.458,4.808], color="red blue", name="Arrows_12.6599998474_7")

cluster_dict["12.6599998474"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(41.0), float(-3.0), float(3.5), float(1.0)]

cluster_dict["12.6599998474_arrows"] += cgo_arrow([41.0,-3.0,3.5], [43.985,-5.458,4.808], color="red blue", name="Arrows_12.6599998474_8")

cmd.load_cgo(cluster_dict["12.6599998474"], "Features_12.6599998474", 1)
cmd.load_cgo(cluster_dict["12.6599998474_arrows"], "Arrows_12.6599998474")
cmd.set("transparency", 0.2,"Features_12.6599998474")
cmd.group("Pharmacophore_12.6599998474", members="Features_12.6599998474")
cmd.group("Pharmacophore_12.6599998474", members="Arrows_12.6599998474")

if dirpath:
    f = join(dirpath, "label_threshold_12.6599998474.mol2")
else:
    f = "label_threshold_12.6599998474.mol2"

cmd.load(f, 'label_threshold_12.6599998474')
cmd.hide('everything', 'label_threshold_12.6599998474')
cmd.label("label_threshold_12.6599998474", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.6599998474', members= 'label_threshold_12.6599998474')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
