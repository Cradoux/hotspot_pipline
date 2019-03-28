
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


cluster_dict = {"10.8440003395":[], "10.8440003395_arrows":[]}

cluster_dict["10.8440003395"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-36.1094658749), float(-46.2779890967), float(-7.4734367213), float(1.0)]


cluster_dict["10.8440003395"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-39.5), float(-46.0), float(-8.0), float(1.0)]

cluster_dict["10.8440003395_arrows"] += cgo_arrow([-39.5,-46.0,-8.0], [-41.624,-46.288,-5.26], color="red blue", name="Arrows_10.8440003395_1")

cluster_dict["10.8440003395"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-33.0), float(-40.5), float(-3.0), float(1.0)]

cluster_dict["10.8440003395_arrows"] += cgo_arrow([-33.0,-40.5,-3.0], [-30.776,-41.63,-2.197], color="red blue", name="Arrows_10.8440003395_2")

cluster_dict["10.8440003395"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-35.0), float(-40.5), float(-1.5), float(1.0)]

cluster_dict["10.8440003395_arrows"] += cgo_arrow([-35.0,-40.5,-1.5], [-36.047,-37.613,-0.7], color="red blue", name="Arrows_10.8440003395_3")

cluster_dict["10.8440003395"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-32.5), float(-50.0), float(-7.0), float(1.0)]

cluster_dict["10.8440003395_arrows"] += cgo_arrow([-32.5,-50.0,-7.0], [-33.632,-52.301,-6.73], color="red blue", name="Arrows_10.8440003395_4")

cluster_dict["10.8440003395"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-33.0), float(-42.5), float(-1.0), float(1.0)]

cluster_dict["10.8440003395_arrows"] += cgo_arrow([-33.0,-42.5,-1.0], [-30.776,-41.63,-2.197], color="red blue", name="Arrows_10.8440003395_5")

cluster_dict["10.8440003395"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-33.0), float(-40.5), float(-3.0), float(1.0)]

cluster_dict["10.8440003395_arrows"] += cgo_arrow([-33.0,-40.5,-3.0], [-30.776,-41.63,-2.197], color="red blue", name="Arrows_10.8440003395_6")

cmd.load_cgo(cluster_dict["10.8440003395"], "Features_10.8440003395", 1)
cmd.load_cgo(cluster_dict["10.8440003395_arrows"], "Arrows_10.8440003395")
cmd.set("transparency", 0.2,"Features_10.8440003395")
cmd.group("Pharmacophore_10.8440003395", members="Features_10.8440003395")
cmd.group("Pharmacophore_10.8440003395", members="Arrows_10.8440003395")

if dirpath:
    f = join(dirpath, "label_threshold_10.8440003395.mol2")
else:
    f = "label_threshold_10.8440003395.mol2"

cmd.load(f, 'label_threshold_10.8440003395')
cmd.hide('everything', 'label_threshold_10.8440003395')
cmd.label("label_threshold_10.8440003395", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.8440003395', members= 'label_threshold_10.8440003395')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")