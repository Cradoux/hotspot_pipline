
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


cluster_dict = {"31.2439994812":[], "31.2439994812_arrows":[]}

cluster_dict["31.2439994812"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(-39.0), float(36.0), float(1.0)]

cluster_dict["31.2439994812_arrows"] += cgo_arrow([30.5,-39.0,36.0], [27.575,-41.176,36.727], color="blue red", name="Arrows_31.2439994812_1")

cluster_dict["31.2439994812"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(32.5), float(-36.0), float(36.0), float(1.0)]

cluster_dict["31.2439994812_arrows"] += cgo_arrow([32.5,-36.0,36.0], [34.071,-33.616,34.804], color="blue red", name="Arrows_31.2439994812_2")

cluster_dict["31.2439994812"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(34.5), float(-34.0), float(38.0), float(1.0)]

cluster_dict["31.2439994812_arrows"] += cgo_arrow([34.5,-34.0,38.0], [34.071,-33.616,34.804], color="blue red", name="Arrows_31.2439994812_3")

cluster_dict["31.2439994812"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(34.5), float(-38.5), float(34.5), float(1.0)]

cluster_dict["31.2439994812_arrows"] += cgo_arrow([34.5,-38.5,34.5], [33.968,-36.953,33.337], color="blue red", name="Arrows_31.2439994812_4")

cluster_dict["31.2439994812"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(32.7758529906), float(-40.3582192195), float(35.2253980264), float(1.0)]


cluster_dict["31.2439994812"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(30.0), float(-43.5), float(36.0), float(1.0)]

cluster_dict["31.2439994812_arrows"] += cgo_arrow([30.0,-43.5,36.0], [28.683,-42.257,38.348], color="red blue", name="Arrows_31.2439994812_5")

cluster_dict["31.2439994812"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(31.5), float(-41.0), float(32.0), float(1.0)]

cluster_dict["31.2439994812_arrows"] += cgo_arrow([31.5,-41.0,32.0], [30.124,-41.797,28.213], color="red blue", name="Arrows_31.2439994812_6")

cluster_dict["31.2439994812"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(32.5), float(-46.0), float(31.0), float(1.0)]

cluster_dict["31.2439994812_arrows"] += cgo_arrow([32.5,-46.0,31.0], [32.469,-46.804,26.319], color="red blue", name="Arrows_31.2439994812_7")

cluster_dict["31.2439994812"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(33.5), float(-39.0), float(38.0), float(1.0)]


cluster_dict["31.2439994812"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(34.0), float(-39.5), float(34.5), float(1.0)]

cluster_dict["31.2439994812_arrows"] += cgo_arrow([34.0,-39.5,34.5], [33.968,-36.953,33.337], color="red blue", name="Arrows_31.2439994812_8")

cmd.load_cgo(cluster_dict["31.2439994812"], "Features_31.2439994812", 1)
cmd.load_cgo(cluster_dict["31.2439994812_arrows"], "Arrows_31.2439994812")
cmd.set("transparency", 0.2,"Features_31.2439994812")
cmd.group("Pharmacophore_31.2439994812", members="Features_31.2439994812")
cmd.group("Pharmacophore_31.2439994812", members="Arrows_31.2439994812")

if dirpath:
    f = join(dirpath, "label_threshold_31.2439994812.mol2")
else:
    f = "label_threshold_31.2439994812.mol2"

cmd.load(f, 'label_threshold_31.2439994812')
cmd.hide('everything', 'label_threshold_31.2439994812')
cmd.label("label_threshold_31.2439994812", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_31.2439994812', members= 'label_threshold_31.2439994812')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
