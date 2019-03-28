
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


cluster_dict = {"16.3819999695":[], "16.3819999695_arrows":[]}

cluster_dict["16.3819999695"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(52.5), float(25.0), float(10.5), float(1.0)]

cluster_dict["16.3819999695_arrows"] += cgo_arrow([52.5,25.0,10.5], [54.742,24.374,11.722], color="blue red", name="Arrows_16.3819999695_1")

cluster_dict["16.3819999695"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(52.5), float(28.5), float(11.5), float(1.0)]

cluster_dict["16.3819999695_arrows"] += cgo_arrow([52.5,28.5,11.5], [54.302,27.83,13.922], color="blue red", name="Arrows_16.3819999695_2")

cluster_dict["16.3819999695"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(56.5), float(30.0), float(18.0), float(1.0)]

cluster_dict["16.3819999695_arrows"] += cgo_arrow([56.5,30.0,18.0], [58.208,32.181,18.361], color="blue red", name="Arrows_16.3819999695_3")

cluster_dict["16.3819999695"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(51.5016134034), float(24.1443093419), float(11.1375211954), float(1.0)]


cluster_dict["16.3819999695"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(47.0), float(27.5), float(15.5), float(1.0)]

cluster_dict["16.3819999695_arrows"] += cgo_arrow([47.0,27.5,15.5], [46.452,25.595,13.666], color="red blue", name="Arrows_16.3819999695_4")

cluster_dict["16.3819999695"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(48.5), float(28.0), float(10.5), float(1.0)]

cluster_dict["16.3819999695_arrows"] += cgo_arrow([48.5,28.0,10.5], [45.233,26.182,9.804], color="red blue", name="Arrows_16.3819999695_5")

cluster_dict["16.3819999695"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(50.5), float(28.5), float(10.0), float(1.0)]

cluster_dict["16.3819999695_arrows"] += cgo_arrow([50.5,28.5,10.0], [50.293,31.808,9.059], color="red blue", name="Arrows_16.3819999695_6")

cluster_dict["16.3819999695"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(51.0), float(30.0), float(7.5), float(1.0)]

cluster_dict["16.3819999695_arrows"] += cgo_arrow([51.0,30.0,7.5], [50.293,31.808,9.059], color="red blue", name="Arrows_16.3819999695_7")

cluster_dict["16.3819999695"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(53.0), float(28.0), float(10.5), float(1.0)]

cluster_dict["16.3819999695_arrows"] += cgo_arrow([53.0,28.0,10.5], [56.329,28.064,9.71], color="red blue", name="Arrows_16.3819999695_8")

cluster_dict["16.3819999695"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(55.5), float(21.0), float(8.0), float(1.0)]

cluster_dict["16.3819999695_arrows"] += cgo_arrow([55.5,21.0,8.0], [56.446,19.309,6.191], color="red blue", name="Arrows_16.3819999695_9")

cmd.load_cgo(cluster_dict["16.3819999695"], "Features_16.3819999695", 1)
cmd.load_cgo(cluster_dict["16.3819999695_arrows"], "Arrows_16.3819999695")
cmd.set("transparency", 0.2,"Features_16.3819999695")
cmd.group("Pharmacophore_16.3819999695", members="Features_16.3819999695")
cmd.group("Pharmacophore_16.3819999695", members="Arrows_16.3819999695")

if dirpath:
    f = join(dirpath, "label_threshold_16.3819999695.mol2")
else:
    f = "label_threshold_16.3819999695.mol2"

cmd.load(f, 'label_threshold_16.3819999695')
cmd.hide('everything', 'label_threshold_16.3819999695')
cmd.label("label_threshold_16.3819999695", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.3819999695', members= 'label_threshold_16.3819999695')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
