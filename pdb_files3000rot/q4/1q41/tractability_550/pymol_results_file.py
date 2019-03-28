
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


cluster_dict = {"17.545999527":[], "17.545999527_arrows":[]}

cluster_dict["17.545999527"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.0), float(-15.5), float(13.5), float(1.0)]

cluster_dict["17.545999527_arrows"] += cgo_arrow([19.0,-15.5,13.5], [19.297,-16.053,16.327], color="blue red", name="Arrows_17.545999527_1")

cluster_dict["17.545999527"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.0), float(-18.5), float(5.5), float(1.0)]

cluster_dict["17.545999527_arrows"] += cgo_arrow([19.0,-18.5,5.5], [19.375,-20.634,3.366], color="blue red", name="Arrows_17.545999527_2")

cluster_dict["17.545999527"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.0), float(-15.5), float(13.5), float(1.0)]

cluster_dict["17.545999527_arrows"] += cgo_arrow([19.0,-15.5,13.5], [19.297,-16.053,16.327], color="blue red", name="Arrows_17.545999527_3")

cluster_dict["17.545999527"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(-20.0), float(7.5), float(1.0)]

cluster_dict["17.545999527_arrows"] += cgo_arrow([24.0,-20.0,7.5], [25.475,-20.386,4.935], color="blue red", name="Arrows_17.545999527_4")

cluster_dict["17.545999527"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(-14.0), float(14.0), float(1.0)]

cluster_dict["17.545999527_arrows"] += cgo_arrow([24.5,-14.0,14.0], [25.791,-11.385,14.706], color="blue red", name="Arrows_17.545999527_5")

cluster_dict["17.545999527"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(-19.0), float(11.5), float(1.0)]

cluster_dict["17.545999527_arrows"] += cgo_arrow([26.5,-19.0,11.5], [24.718,-20.334,13.2], color="blue red", name="Arrows_17.545999527_6")

cluster_dict["17.545999527"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.2722396541), float(-17.3924225738), float(9.72293026761), float(1.0)]


cluster_dict["17.545999527"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(-16.5), float(6.5), float(1.0)]

cluster_dict["17.545999527_arrows"] += cgo_arrow([17.5,-16.5,6.5], [17.502,-12.43,6.943], color="red blue", name="Arrows_17.545999527_7")

cluster_dict["17.545999527"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(-20.0), float(6.0), float(1.0)]

cluster_dict["17.545999527_arrows"] += cgo_arrow([22.0,-20.0,6.0], [23.476,-21.001,3.239], color="red blue", name="Arrows_17.545999527_8")

cluster_dict["17.545999527"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(-19.0), float(12.0), float(1.0)]

cluster_dict["17.545999527_arrows"] += cgo_arrow([22.0,-19.0,12.0], [22.758,-20.731,14.24], color="red blue", name="Arrows_17.545999527_9")

cluster_dict["17.545999527"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(-22.5), float(9.0), float(1.0)]

cluster_dict["17.545999527_arrows"] += cgo_arrow([26.5,-22.5,9.0], [26.733,-24.229,6.442], color="red blue", name="Arrows_17.545999527_10")

cmd.load_cgo(cluster_dict["17.545999527"], "Features_17.545999527", 1)
cmd.load_cgo(cluster_dict["17.545999527_arrows"], "Arrows_17.545999527")
cmd.set("transparency", 0.2,"Features_17.545999527")
cmd.group("Pharmacophore_17.545999527", members="Features_17.545999527")
cmd.group("Pharmacophore_17.545999527", members="Arrows_17.545999527")

if dirpath:
    f = join(dirpath, "label_threshold_17.545999527.mol2")
else:
    f = "label_threshold_17.545999527.mol2"

cmd.load(f, 'label_threshold_17.545999527')
cmd.hide('everything', 'label_threshold_17.545999527')
cmd.label("label_threshold_17.545999527", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.545999527', members= 'label_threshold_17.545999527')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
