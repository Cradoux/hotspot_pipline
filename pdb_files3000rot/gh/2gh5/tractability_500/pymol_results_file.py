
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


cluster_dict = {"13.8020000458":[], "13.8020000458_arrows":[]}

cluster_dict["13.8020000458"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(95.5), float(16.5), float(100.0), float(1.0)]

cluster_dict["13.8020000458_arrows"] += cgo_arrow([95.5,16.5,100.0], [95.246,15.863,102.643], color="blue red", name="Arrows_13.8020000458_1")

cluster_dict["13.8020000458"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(97.5), float(17.0), float(88.0), float(1.0)]

cluster_dict["13.8020000458_arrows"] += cgo_arrow([97.5,17.0,88.0], [95.688,16.84,85.897], color="blue red", name="Arrows_13.8020000458_2")

cluster_dict["13.8020000458"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(99.5), float(14.0), float(101.0), float(1.0)]

cluster_dict["13.8020000458_arrows"] += cgo_arrow([99.5,14.0,101.0], [99.4,10.661,101.72], color="blue red", name="Arrows_13.8020000458_3")

cluster_dict["13.8020000458"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(97.7424501654), float(19.9729687839), float(95.3940129435), float(1.0)]


cluster_dict["13.8020000458"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(101.339607178), float(12.2402945638), float(90.8835296507), float(1.0)]


cluster_dict["13.8020000458"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(96.0), float(18.0), float(99.0), float(1.0)]

cluster_dict["13.8020000458_arrows"] += cgo_arrow([96.0,18.0,99.0], [93.358,18.6,99.599], color="red blue", name="Arrows_13.8020000458_4")

cluster_dict["13.8020000458"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(99.5), float(24.0), float(92.0), float(1.0)]

cluster_dict["13.8020000458_arrows"] += cgo_arrow([99.5,24.0,92.0], [99.478,26.029,94.674], color="red blue", name="Arrows_13.8020000458_5")

cluster_dict["13.8020000458"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(101.0), float(21.5), float(100.0), float(1.0)]

cluster_dict["13.8020000458_arrows"] += cgo_arrow([101.0,21.5,100.0], [101.55,19.69,101.979], color="red blue", name="Arrows_13.8020000458_6")

cluster_dict["13.8020000458"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(101.5), float(10.0), float(89.5), float(1.0)]

cluster_dict["13.8020000458_arrows"] += cgo_arrow([101.5,10.0,89.5], [99.053,9.698,91.013], color="red blue", name="Arrows_13.8020000458_7")

cmd.load_cgo(cluster_dict["13.8020000458"], "Features_13.8020000458", 1)
cmd.load_cgo(cluster_dict["13.8020000458_arrows"], "Arrows_13.8020000458")
cmd.set("transparency", 0.2,"Features_13.8020000458")
cmd.group("Pharmacophore_13.8020000458", members="Features_13.8020000458")
cmd.group("Pharmacophore_13.8020000458", members="Arrows_13.8020000458")

if dirpath:
    f = join(dirpath, "label_threshold_13.8020000458.mol2")
else:
    f = "label_threshold_13.8020000458.mol2"

cmd.load(f, 'label_threshold_13.8020000458')
cmd.hide('everything', 'label_threshold_13.8020000458')
cmd.label("label_threshold_13.8020000458", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.8020000458', members= 'label_threshold_13.8020000458')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
