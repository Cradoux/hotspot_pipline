
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


cluster_dict = {"12.8680000305":[], "12.8680000305_arrows":[]}

cluster_dict["12.8680000305"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(39.0), float(54.5), float(0.5), float(1.0)]

cluster_dict["12.8680000305_arrows"] += cgo_arrow([39.0,54.5,0.5], [38.382,52.515,1.817], color="blue red", name="Arrows_12.8680000305_1")

cluster_dict["12.8680000305"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(43.5), float(50.5), float(1.0), float(1.0)]

cluster_dict["12.8680000305_arrows"] += cgo_arrow([43.5,50.5,1.0], [40.801,50.178,1.888], color="blue red", name="Arrows_12.8680000305_2")

cluster_dict["12.8680000305"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(42.0276292132), float(54.4410967158), float(-0.329018111108), float(1.0)]


cluster_dict["12.8680000305"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(43.1727521417), float(54.1408480037), float(3.05685061756), float(1.0)]


cluster_dict["12.8680000305"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(42.5), float(53.0), float(2.0), float(1.0)]

cluster_dict["12.8680000305_arrows"] += cgo_arrow([42.5,53.0,2.0], [44.208,54.939,3.202], color="red blue", name="Arrows_12.8680000305_3")

cluster_dict["12.8680000305"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(43.0), float(55.5), float(-2.5), float(1.0)]

cluster_dict["12.8680000305_arrows"] += cgo_arrow([43.0,55.5,-2.5], [46.103,56.982,-0.149], color="red blue", name="Arrows_12.8680000305_4")

cluster_dict["12.8680000305"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(42.5), float(49.5), float(-3.0), float(1.0)]

cluster_dict["12.8680000305_arrows"] += cgo_arrow([42.5,49.5,-3.0], [42.385,51.023,-5.187], color="red blue", name="Arrows_12.8680000305_5")

cmd.load_cgo(cluster_dict["12.8680000305"], "Features_12.8680000305", 1)
cmd.load_cgo(cluster_dict["12.8680000305_arrows"], "Arrows_12.8680000305")
cmd.set("transparency", 0.2,"Features_12.8680000305")
cmd.group("Pharmacophore_12.8680000305", members="Features_12.8680000305")
cmd.group("Pharmacophore_12.8680000305", members="Arrows_12.8680000305")

if dirpath:
    f = join(dirpath, "label_threshold_12.8680000305.mol2")
else:
    f = "label_threshold_12.8680000305.mol2"

cmd.load(f, 'label_threshold_12.8680000305')
cmd.hide('everything', 'label_threshold_12.8680000305')
cmd.label("label_threshold_12.8680000305", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.8680000305', members= 'label_threshold_12.8680000305')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
