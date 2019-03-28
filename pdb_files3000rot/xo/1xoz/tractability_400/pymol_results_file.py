
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


cluster_dict = {"19.8369998932":[], "19.8369998932_arrows":[]}

cluster_dict["19.8369998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(41.5), float(36.5), float(15.5), float(1.0)]

cluster_dict["19.8369998932_arrows"] += cgo_arrow([41.5,36.5,15.5], [39.181,37.637,16.548], color="blue red", name="Arrows_19.8369998932_1")

cluster_dict["19.8369998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(43.0), float(39.0), float(19.5), float(1.0)]

cluster_dict["19.8369998932_arrows"] += cgo_arrow([43.0,39.0,19.5], [45.249,40.474,18.915], color="blue red", name="Arrows_19.8369998932_2")

cluster_dict["19.8369998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(45.0), float(36.5), float(10.0), float(1.0)]

cluster_dict["19.8369998932_arrows"] += cgo_arrow([45.0,36.5,10.0], [46.98,35.513,7.652], color="blue red", name="Arrows_19.8369998932_3")

cluster_dict["19.8369998932"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(47.2498642336), float(34.9207494515), float(12.8610107653), float(1.0)]


cluster_dict["19.8369998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(45.0), float(36.5), float(12.0), float(1.0)]

cluster_dict["19.8369998932_arrows"] += cgo_arrow([45.0,36.5,12.0], [45.729,33.684,7.995], color="red blue", name="Arrows_19.8369998932_4")

cluster_dict["19.8369998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(42.5), float(38.0), float(17.0), float(1.0)]

cluster_dict["19.8369998932_arrows"] += cgo_arrow([42.5,38.0,17.0], [41.879,39.772,14.796], color="red blue", name="Arrows_19.8369998932_5")

cluster_dict["19.8369998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(45.0), float(36.5), float(12.0), float(1.0)]

cluster_dict["19.8369998932_arrows"] += cgo_arrow([45.0,36.5,12.0], [45.729,33.684,7.995], color="red blue", name="Arrows_19.8369998932_6")

cluster_dict["19.8369998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(46.0), float(37.5), float(10.0), float(1.0)]

cluster_dict["19.8369998932_arrows"] += cgo_arrow([46.0,37.5,10.0], [45.729,33.684,7.995], color="red blue", name="Arrows_19.8369998932_7")

cluster_dict["19.8369998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(48.5), float(37.5), float(19.5), float(1.0)]

cluster_dict["19.8369998932_arrows"] += cgo_arrow([48.5,37.5,19.5], [46.853,41.904,18.83], color="red blue", name="Arrows_19.8369998932_8")

cluster_dict["19.8369998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(51.0), float(34.5), float(8.0), float(1.0)]

cluster_dict["19.8369998932_arrows"] += cgo_arrow([51.0,34.5,8.0], [53.535,33.389,4.555], color="red blue", name="Arrows_19.8369998932_9")

cmd.load_cgo(cluster_dict["19.8369998932"], "Features_19.8369998932", 1)
cmd.load_cgo(cluster_dict["19.8369998932_arrows"], "Arrows_19.8369998932")
cmd.set("transparency", 0.2,"Features_19.8369998932")
cmd.group("Pharmacophore_19.8369998932", members="Features_19.8369998932")
cmd.group("Pharmacophore_19.8369998932", members="Arrows_19.8369998932")

if dirpath:
    f = join(dirpath, "label_threshold_19.8369998932.mol2")
else:
    f = "label_threshold_19.8369998932.mol2"

cmd.load(f, 'label_threshold_19.8369998932')
cmd.hide('everything', 'label_threshold_19.8369998932')
cmd.label("label_threshold_19.8369998932", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_19.8369998932', members= 'label_threshold_19.8369998932')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")