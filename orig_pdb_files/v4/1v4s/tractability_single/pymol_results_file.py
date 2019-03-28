
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


cluster_dict = {"14.3999996185":[], "14.3999996185_arrows":[]}

cluster_dict["14.3999996185"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(19.0), float(57.0), float(1.0)]

cluster_dict["14.3999996185_arrows"] += cgo_arrow([36.5,19.0,57.0], [34.748,16.322,56.81], color="blue red", name="Arrows_14.3999996185_1")

cluster_dict["14.3999996185"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(39.0), float(17.5), float(61.0), float(1.0)]

cluster_dict["14.3999996185_arrows"] += cgo_arrow([39.0,17.5,61.0], [40.059,15.519,58.856], color="blue red", name="Arrows_14.3999996185_2")

cluster_dict["14.3999996185"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(40.0), float(24.0), float(60.0), float(1.0)]

cluster_dict["14.3999996185_arrows"] += cgo_arrow([40.0,24.0,60.0], [37.531,24.144,61.642], color="blue red", name="Arrows_14.3999996185_3")

cluster_dict["14.3999996185"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(39.5291315519), float(17.409878203), float(60.7613772988), float(1.0)]


cluster_dict["14.3999996185"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(36.7603060458), float(17.4921525954), float(63.7603060458), float(1.0)]


cluster_dict["14.3999996185"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(18.0), float(59.0), float(1.0)]

cluster_dict["14.3999996185_arrows"] += cgo_arrow([36.5,18.0,59.0], [37.483,15.325,57.746], color="red blue", name="Arrows_14.3999996185_4")

cluster_dict["14.3999996185"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(37.5), float(18.0), float(62.5), float(1.0)]

cluster_dict["14.3999996185_arrows"] += cgo_arrow([37.5,18.0,62.5], [35.617,18.465,66.87], color="red blue", name="Arrows_14.3999996185_5")

cluster_dict["14.3999996185"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(39.5), float(16.5), float(65.5), float(1.0)]

cluster_dict["14.3999996185_arrows"] += cgo_arrow([39.5,16.5,65.5], [42.04,17.748,66.666], color="red blue", name="Arrows_14.3999996185_6")

cluster_dict["14.3999996185"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(41.0), float(11.0), float(64.0), float(1.0)]

cluster_dict["14.3999996185_arrows"] += cgo_arrow([41.0,11.0,64.0], [37.79,10.544,67.284], color="red blue", name="Arrows_14.3999996185_7")

cluster_dict["14.3999996185"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(43.0), float(12.5), float(61.0), float(1.0)]

cluster_dict["14.3999996185_arrows"] += cgo_arrow([43.0,12.5,61.0], [43.012,15.508,59.647], color="red blue", name="Arrows_14.3999996185_8")

cluster_dict["14.3999996185"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(43.0), float(18.0), float(63.0), float(1.0)]

cluster_dict["14.3999996185_arrows"] += cgo_arrow([43.0,18.0,63.0], [45.018,17.321,62.219], color="red blue", name="Arrows_14.3999996185_9")

cmd.load_cgo(cluster_dict["14.3999996185"], "Features_14.3999996185", 1)
cmd.load_cgo(cluster_dict["14.3999996185_arrows"], "Arrows_14.3999996185")
cmd.set("transparency", 0.2,"Features_14.3999996185")
cmd.group("Pharmacophore_14.3999996185", members="Features_14.3999996185")
cmd.group("Pharmacophore_14.3999996185", members="Arrows_14.3999996185")

if dirpath:
    f = join(dirpath, "label_threshold_14.3999996185.mol2")
else:
    f = "label_threshold_14.3999996185.mol2"

cmd.load(f, 'label_threshold_14.3999996185')
cmd.hide('everything', 'label_threshold_14.3999996185')
cmd.label("label_threshold_14.3999996185", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.3999996185', members= 'label_threshold_14.3999996185')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
