
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


cluster_dict = {"14.5240001678":[], "14.5240001678_arrows":[]}

cluster_dict["14.5240001678"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-32.0), float(-10.0), float(20.0), float(1.0)]

cluster_dict["14.5240001678_arrows"] += cgo_arrow([-32.0,-10.0,20.0], [-32.693,-13.194,18.406], color="blue red", name="Arrows_14.5240001678_1")

cluster_dict["14.5240001678"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-31.0), float(-11.5), float(22.5), float(1.0)]

cluster_dict["14.5240001678_arrows"] += cgo_arrow([-31.0,-11.5,22.5], [-28.145,-12.556,19.953], color="blue red", name="Arrows_14.5240001678_2")

cluster_dict["14.5240001678"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-32.5), float(-6.0), float(23.5), float(1.0)]

cluster_dict["14.5240001678_arrows"] += cgo_arrow([-32.5,-6.0,23.5], [-33.765,-7.447,25.176], color="blue red", name="Arrows_14.5240001678_3")

cluster_dict["14.5240001678"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-31.0), float(-11.5), float(22.5), float(1.0)]

cluster_dict["14.5240001678_arrows"] += cgo_arrow([-31.0,-11.5,22.5], [-28.145,-12.556,19.953], color="blue red", name="Arrows_14.5240001678_4")

cluster_dict["14.5240001678"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-27.5), float(-17.5), float(27.5), float(1.0)]

cluster_dict["14.5240001678_arrows"] += cgo_arrow([-27.5,-17.5,27.5], [-29.467,-18.731,29.67], color="blue red", name="Arrows_14.5240001678_5")

cluster_dict["14.5240001678"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-30.8911130531), float(-12.3238134027), float(23.9837660693), float(1.0)]


cluster_dict["14.5240001678"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-30.5), float(-12.5), float(28.0), float(1.0)]

cluster_dict["14.5240001678_arrows"] += cgo_arrow([-30.5,-12.5,28.0], [-29.806,-8.212,28.81], color="red blue", name="Arrows_14.5240001678_6")

cluster_dict["14.5240001678"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-34.5), float(-9.0), float(18.0), float(1.0)]

cluster_dict["14.5240001678_arrows"] += cgo_arrow([-34.5,-9.0,18.0], [-35.149,-8.851,17.552], color="red blue", name="Arrows_14.5240001678_7")

cluster_dict["14.5240001678"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-30.5), float(-12.0), float(23.5), float(1.0)]

cluster_dict["14.5240001678_arrows"] += cgo_arrow([-30.5,-12.0,23.5], [-26.862,-11.631,22.737], color="red blue", name="Arrows_14.5240001678_8")

cluster_dict["14.5240001678"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-30.5), float(-12.5), float(28.0), float(1.0)]

cluster_dict["14.5240001678_arrows"] += cgo_arrow([-30.5,-12.5,28.0], [-29.806,-8.212,28.81], color="red blue", name="Arrows_14.5240001678_9")

cluster_dict["14.5240001678"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-27.0), float(-14.0), float(26.5), float(1.0)]

cluster_dict["14.5240001678_arrows"] += cgo_arrow([-27.0,-14.0,26.5], [-24.688,-13.685,28.563], color="red blue", name="Arrows_14.5240001678_10")

cmd.load_cgo(cluster_dict["14.5240001678"], "Features_14.5240001678", 1)
cmd.load_cgo(cluster_dict["14.5240001678_arrows"], "Arrows_14.5240001678")
cmd.set("transparency", 0.2,"Features_14.5240001678")
cmd.group("Pharmacophore_14.5240001678", members="Features_14.5240001678")
cmd.group("Pharmacophore_14.5240001678", members="Arrows_14.5240001678")

if dirpath:
    f = join(dirpath, "label_threshold_14.5240001678.mol2")
else:
    f = "label_threshold_14.5240001678.mol2"

cmd.load(f, 'label_threshold_14.5240001678')
cmd.hide('everything', 'label_threshold_14.5240001678')
cmd.label("label_threshold_14.5240001678", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.5240001678', members= 'label_threshold_14.5240001678')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
