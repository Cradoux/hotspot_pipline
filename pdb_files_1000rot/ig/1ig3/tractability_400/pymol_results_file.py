
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


cluster_dict = {"15.4750003815":[], "15.4750003815_arrows":[]}

cluster_dict["15.4750003815"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(42.0), float(20.0), float(1.0)]

cluster_dict["15.4750003815_arrows"] += cgo_arrow([23.5,42.0,20.0], [21.792,42.89,17.917], color="blue red", name="Arrows_15.4750003815_1")

cluster_dict["15.4750003815"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(43.0), float(23.0), float(1.0)]

cluster_dict["15.4750003815_arrows"] += cgo_arrow([23.0,43.0,23.0], [26.259,42.843,23.416], color="blue red", name="Arrows_15.4750003815_2")

cluster_dict["15.4750003815"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(45.0), float(14.0), float(1.0)]

cluster_dict["15.4750003815_arrows"] += cgo_arrow([27.0,45.0,14.0], [23.869,44.199,14.749], color="blue red", name="Arrows_15.4750003815_3")

cluster_dict["15.4750003815"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(46.5), float(16.0), float(1.0)]

cluster_dict["15.4750003815_arrows"] += cgo_arrow([27.0,46.5,16.0], [25.327,48.612,14.609], color="blue red", name="Arrows_15.4750003815_4")

cluster_dict["15.4750003815"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(42.0), float(16.5), float(1.0)]

cluster_dict["15.4750003815_arrows"] += cgo_arrow([30.5,42.0,16.5], [32.135,42.72,14.753], color="blue red", name="Arrows_15.4750003815_5")

cluster_dict["15.4750003815"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(23.953361701), float(40.9277963811), float(20.5891864046), float(1.0)]


cluster_dict["15.4750003815"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(44.0), float(20.5), float(1.0)]

cluster_dict["15.4750003815_arrows"] += cgo_arrow([22.5,44.0,20.5], [22.183,45.976,18.664], color="red blue", name="Arrows_15.4750003815_6")

cluster_dict["15.4750003815"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(44.5), float(24.0), float(1.0)]

cluster_dict["15.4750003815_arrows"] += cgo_arrow([22.5,44.5,24.0], [19.407,45.482,25.054], color="red blue", name="Arrows_15.4750003815_7")

cluster_dict["15.4750003815"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(45.0), float(15.5), float(1.0)]

cluster_dict["15.4750003815_arrows"] += cgo_arrow([27.0,45.0,15.5], [23.869,44.199,14.749], color="red blue", name="Arrows_15.4750003815_8")

cmd.load_cgo(cluster_dict["15.4750003815"], "Features_15.4750003815", 1)
cmd.load_cgo(cluster_dict["15.4750003815_arrows"], "Arrows_15.4750003815")
cmd.set("transparency", 0.2,"Features_15.4750003815")
cmd.group("Pharmacophore_15.4750003815", members="Features_15.4750003815")
cmd.group("Pharmacophore_15.4750003815", members="Arrows_15.4750003815")

if dirpath:
    f = join(dirpath, "label_threshold_15.4750003815.mol2")
else:
    f = "label_threshold_15.4750003815.mol2"

cmd.load(f, 'label_threshold_15.4750003815')
cmd.hide('everything', 'label_threshold_15.4750003815')
cmd.label("label_threshold_15.4750003815", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.4750003815', members= 'label_threshold_15.4750003815')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
