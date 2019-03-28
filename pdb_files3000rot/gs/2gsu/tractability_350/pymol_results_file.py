
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


cluster_dict = {"14.8459997177":[], "14.8459997177_arrows":[]}

cluster_dict["14.8459997177"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(40.0), float(19.0), float(1.0)]

cluster_dict["14.8459997177_arrows"] += cgo_arrow([22.0,40.0,19.0], [19.46,38.609,18.635], color="blue red", name="Arrows_14.8459997177_1")

cluster_dict["14.8459997177"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(38.0), float(15.5), float(1.0)]

cluster_dict["14.8459997177_arrows"] += cgo_arrow([23.0,38.0,15.5], [21.742,36.194,13.802], color="blue red", name="Arrows_14.8459997177_2")

cluster_dict["14.8459997177"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(38.0), float(13.5), float(1.0)]

cluster_dict["14.8459997177_arrows"] += cgo_arrow([23.5,38.0,13.5], [21.742,36.194,13.802], color="blue red", name="Arrows_14.8459997177_3")

cluster_dict["14.8459997177"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(38.5), float(14.0), float(1.0)]

cluster_dict["14.8459997177_arrows"] += cgo_arrow([25.0,38.5,14.0], [28.084,38.378,13.672], color="blue red", name="Arrows_14.8459997177_4")

cluster_dict["14.8459997177"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(39.5), float(10.5), float(1.0)]

cluster_dict["14.8459997177_arrows"] += cgo_arrow([25.0,39.5,10.5], [22.076,40.437,10.904], color="blue red", name="Arrows_14.8459997177_5")

cluster_dict["14.8459997177"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(25.1270400898), float(41.1826912763), float(17.4134310624), float(1.0)]


cluster_dict["14.8459997177"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(45.0), float(18.5), float(1.0)]

cluster_dict["14.8459997177_arrows"] += cgo_arrow([23.5,45.0,18.5], [22.695,43.688,21.05], color="red blue", name="Arrows_14.8459997177_6")

cluster_dict["14.8459997177"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(37.0), float(15.0), float(1.0)]

cluster_dict["14.8459997177_arrows"] += cgo_arrow([24.0,37.0,15.0], [21.742,36.194,13.802], color="red blue", name="Arrows_14.8459997177_7")

cluster_dict["14.8459997177"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(42.0), float(13.0), float(1.0)]

cluster_dict["14.8459997177_arrows"] += cgo_arrow([26.5,42.0,13.0], [25.832,44.224,13.785], color="red blue", name="Arrows_14.8459997177_8")

cluster_dict["14.8459997177"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(43.0), float(16.0), float(1.0)]

cluster_dict["14.8459997177_arrows"] += cgo_arrow([26.5,43.0,16.0], [25.832,44.224,13.785], color="red blue", name="Arrows_14.8459997177_9")

cmd.load_cgo(cluster_dict["14.8459997177"], "Features_14.8459997177", 1)
cmd.load_cgo(cluster_dict["14.8459997177_arrows"], "Arrows_14.8459997177")
cmd.set("transparency", 0.2,"Features_14.8459997177")
cmd.group("Pharmacophore_14.8459997177", members="Features_14.8459997177")
cmd.group("Pharmacophore_14.8459997177", members="Arrows_14.8459997177")

if dirpath:
    f = join(dirpath, "label_threshold_14.8459997177.mol2")
else:
    f = "label_threshold_14.8459997177.mol2"

cmd.load(f, 'label_threshold_14.8459997177')
cmd.hide('everything', 'label_threshold_14.8459997177')
cmd.label("label_threshold_14.8459997177", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.8459997177', members= 'label_threshold_14.8459997177')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
