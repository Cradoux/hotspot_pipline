
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


cluster_dict = {"15.0839996338":[], "15.0839996338_arrows":[]}

cluster_dict["15.0839996338"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(32.0), float(57.0), float(30.0), float(1.0)]

cluster_dict["15.0839996338_arrows"] += cgo_arrow([32.0,57.0,30.0], [32.435,56.176,27.574], color="blue red", name="Arrows_15.0839996338_1")

cluster_dict["15.0839996338"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(35.0), float(56.0), float(37.5), float(1.0)]

cluster_dict["15.0839996338_arrows"] += cgo_arrow([35.0,56.0,37.5], [32.702,57.592,38.987], color="blue red", name="Arrows_15.0839996338_2")

cluster_dict["15.0839996338"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(39.0), float(56.5), float(39.5), float(1.0)]

cluster_dict["15.0839996338_arrows"] += cgo_arrow([39.0,56.5,39.5], [37.874,55.634,42.017], color="blue red", name="Arrows_15.0839996338_3")

cluster_dict["15.0839996338"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(31.8400246118), float(57.8400246118), float(30.5), float(1.0)]


cluster_dict["15.0839996338"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(38.2578742094), float(56.1095176685), float(34.1180738027), float(1.0)]


cluster_dict["15.0839996338"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(33.7654919872), float(57.7654919872), float(36.5), float(1.0)]


cluster_dict["15.0839996338"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(33.5), float(54.5), float(33.5), float(1.0)]

cluster_dict["15.0839996338_arrows"] += cgo_arrow([33.5,54.5,33.5], [34.74,51.905,31.389], color="red blue", name="Arrows_15.0839996338_4")

cluster_dict["15.0839996338"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(35.5), float(56.0), float(29.0), float(1.0)]

cluster_dict["15.0839996338_arrows"] += cgo_arrow([35.5,56.0,29.0], [34.661,58.006,26.634], color="red blue", name="Arrows_15.0839996338_5")

cluster_dict["15.0839996338"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(57.0), float(39.0), float(1.0)]

cluster_dict["15.0839996338_arrows"] += cgo_arrow([36.5,57.0,39.0], [35.202,56.434,41.571], color="red blue", name="Arrows_15.0839996338_6")

cluster_dict["15.0839996338"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(40.5), float(57.5), float(37.0), float(1.0)]


cluster_dict["15.0839996338"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(41.5), float(58.5), float(42.0), float(1.0)]

cluster_dict["15.0839996338_arrows"] += cgo_arrow([41.5,58.5,42.0], [39.664,59.531,44.321], color="red blue", name="Arrows_15.0839996338_7")

cmd.load_cgo(cluster_dict["15.0839996338"], "Features_15.0839996338", 1)
cmd.load_cgo(cluster_dict["15.0839996338_arrows"], "Arrows_15.0839996338")
cmd.set("transparency", 0.2,"Features_15.0839996338")
cmd.group("Pharmacophore_15.0839996338", members="Features_15.0839996338")
cmd.group("Pharmacophore_15.0839996338", members="Arrows_15.0839996338")

if dirpath:
    f = join(dirpath, "label_threshold_15.0839996338.mol2")
else:
    f = "label_threshold_15.0839996338.mol2"

cmd.load(f, 'label_threshold_15.0839996338')
cmd.hide('everything', 'label_threshold_15.0839996338')
cmd.label("label_threshold_15.0839996338", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.0839996338', members= 'label_threshold_15.0839996338')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
