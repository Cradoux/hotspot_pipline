
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


cluster_dict = {"13.1440000534":[], "13.1440000534_arrows":[]}

cluster_dict["13.1440000534"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(11.5), float(-14.0), float(1.0)]

cluster_dict["13.1440000534_arrows"] += cgo_arrow([25.0,11.5,-14.0], [22.615,11.131,-15.816], color="blue red", name="Arrows_13.1440000534_1")

cluster_dict["13.1440000534"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(14.5), float(-14.5), float(1.0)]

cluster_dict["13.1440000534_arrows"] += cgo_arrow([25.0,14.5,-14.5], [21.438,14.463,-13.601], color="blue red", name="Arrows_13.1440000534_2")

cluster_dict["13.1440000534"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(18.0), float(-16.5), float(1.0)]

cluster_dict["13.1440000534_arrows"] += cgo_arrow([24.5,18.0,-16.5], [21.965,15.819,-16.753], color="blue red", name="Arrows_13.1440000534_3")

cluster_dict["13.1440000534"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(11.5), float(-14.0), float(1.0)]

cluster_dict["13.1440000534_arrows"] += cgo_arrow([25.0,11.5,-14.0], [22.615,11.131,-15.816], color="blue red", name="Arrows_13.1440000534_4")

cluster_dict["13.1440000534"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(12.0), float(-23.0), float(1.0)]

cluster_dict["13.1440000534_arrows"] += cgo_arrow([29.5,12.0,-23.0], [28.974,12.038,-25.36], color="blue red", name="Arrows_13.1440000534_5")

cluster_dict["13.1440000534"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(10.0), float(-13.0), float(1.0)]

cluster_dict["13.1440000534_arrows"] += cgo_arrow([29.5,10.0,-13.0], [30.129,7.657,-14.523], color="blue red", name="Arrows_13.1440000534_6")

cluster_dict["13.1440000534"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.0), float(4.5), float(-13.0), float(1.0)]

cluster_dict["13.1440000534_arrows"] += cgo_arrow([30.0,4.5,-13.0], [30.724,3.981,-15.898], color="blue red", name="Arrows_13.1440000534_7")

cluster_dict["13.1440000534"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.0522414329), float(11.3736451211), float(-15.9760071374), float(1.0)]


cluster_dict["13.1440000534"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(10.0), float(-20.5), float(1.0)]

cluster_dict["13.1440000534_arrows"] += cgo_arrow([26.0,10.0,-20.5], [23.813,11.613,-19.373], color="red blue", name="Arrows_13.1440000534_8")

cluster_dict["13.1440000534"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(28.5), float(11.0), float(-22.0), float(1.0)]

cluster_dict["13.1440000534_arrows"] += cgo_arrow([28.5,11.0,-22.0], [26.191,13.098,-22.189], color="red blue", name="Arrows_13.1440000534_9")

cluster_dict["13.1440000534"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(10.5), float(-13.0), float(1.0)]

cluster_dict["13.1440000534_arrows"] += cgo_arrow([29.5,10.5,-13.0], [32.505,10.662,-12.636], color="red blue", name="Arrows_13.1440000534_10")

cmd.load_cgo(cluster_dict["13.1440000534"], "Features_13.1440000534", 1)
cmd.load_cgo(cluster_dict["13.1440000534_arrows"], "Arrows_13.1440000534")
cmd.set("transparency", 0.2,"Features_13.1440000534")
cmd.group("Pharmacophore_13.1440000534", members="Features_13.1440000534")
cmd.group("Pharmacophore_13.1440000534", members="Arrows_13.1440000534")

if dirpath:
    f = join(dirpath, "label_threshold_13.1440000534.mol2")
else:
    f = "label_threshold_13.1440000534.mol2"

cmd.load(f, 'label_threshold_13.1440000534')
cmd.hide('everything', 'label_threshold_13.1440000534')
cmd.label("label_threshold_13.1440000534", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.1440000534', members= 'label_threshold_13.1440000534')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
