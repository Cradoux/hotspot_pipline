
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


cluster_dict = {"17.2399997711":[], "17.2399997711_arrows":[]}

cluster_dict["17.2399997711"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(43.0), float(2.5), float(-8.5), float(1.0)]

cluster_dict["17.2399997711_arrows"] += cgo_arrow([43.0,2.5,-8.5], [41.697,2.863,-5.936], color="blue red", name="Arrows_17.2399997711_1")

cluster_dict["17.2399997711"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(44.0), float(1.0), float(-10.5), float(1.0)]

cluster_dict["17.2399997711_arrows"] += cgo_arrow([44.0,1.0,-10.5], [41.901,0.681,-13.01], color="blue red", name="Arrows_17.2399997711_2")

cluster_dict["17.2399997711"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(49.0), float(4.0), float(-5.0), float(1.0)]

cluster_dict["17.2399997711_arrows"] += cgo_arrow([49.0,4.0,-5.0], [48.938,6.518,-5.861], color="blue red", name="Arrows_17.2399997711_3")

cluster_dict["17.2399997711"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(53.0), float(3.0), float(-15.5), float(1.0)]

cluster_dict["17.2399997711_arrows"] += cgo_arrow([53.0,3.0,-15.5], [53.133,0.151,-16.057], color="blue red", name="Arrows_17.2399997711_4")

cluster_dict["17.2399997711"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(52.5), float(4.5), float(-12.5), float(1.0)]

cluster_dict["17.2399997711_arrows"] += cgo_arrow([52.5,4.5,-12.5], [53.874,1.645,-12.802], color="blue red", name="Arrows_17.2399997711_5")

cluster_dict["17.2399997711"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(52.0), float(6.5), float(-12.5), float(1.0)]

cluster_dict["17.2399997711_arrows"] += cgo_arrow([52.0,6.5,-12.5], [53.758,8.584,-13.579], color="blue red", name="Arrows_17.2399997711_6")

cluster_dict["17.2399997711"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(47.8797352213), float(1.76940427729), float(-8.8529982529), float(1.0)]


cluster_dict["17.2399997711"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(51.8420978291), float(5.12592766064), float(-15.6567064325), float(1.0)]


cluster_dict["17.2399997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(44.0), float(0.0), float(-6.0), float(1.0)]

cluster_dict["17.2399997711_arrows"] += cgo_arrow([44.0,0.0,-6.0], [40.435,-0.517,-3.997], color="red blue", name="Arrows_17.2399997711_7")

cluster_dict["17.2399997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(49.5), float(4.0), float(-9.0), float(1.0)]

cluster_dict["17.2399997711_arrows"] += cgo_arrow([49.5,4.0,-9.0], [49.051,8.384,-6.905], color="red blue", name="Arrows_17.2399997711_8")

cluster_dict["17.2399997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(50.0), float(1.5), float(-15.0), float(1.0)]

cluster_dict["17.2399997711_arrows"] += cgo_arrow([50.0,1.5,-15.0], [52.853,-1.315,-14.363], color="red blue", name="Arrows_17.2399997711_9")

cmd.load_cgo(cluster_dict["17.2399997711"], "Features_17.2399997711", 1)
cmd.load_cgo(cluster_dict["17.2399997711_arrows"], "Arrows_17.2399997711")
cmd.set("transparency", 0.2,"Features_17.2399997711")
cmd.group("Pharmacophore_17.2399997711", members="Features_17.2399997711")
cmd.group("Pharmacophore_17.2399997711", members="Arrows_17.2399997711")

if dirpath:
    f = join(dirpath, "label_threshold_17.2399997711.mol2")
else:
    f = "label_threshold_17.2399997711.mol2"

cmd.load(f, 'label_threshold_17.2399997711')
cmd.hide('everything', 'label_threshold_17.2399997711')
cmd.label("label_threshold_17.2399997711", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.2399997711', members= 'label_threshold_17.2399997711')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
