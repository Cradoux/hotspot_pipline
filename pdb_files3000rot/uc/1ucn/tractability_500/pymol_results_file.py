
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


cluster_dict = {"17.7019996643":[], "17.7019996643_arrows":[]}

cluster_dict["17.7019996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-12.5), float(92.0), float(25.5), float(1.0)]

cluster_dict["17.7019996643_arrows"] += cgo_arrow([-12.5,92.0,25.5], [-11.32,91.433,27.948], color="blue red", name="Arrows_17.7019996643_1")

cluster_dict["17.7019996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-13.5), float(103.0), float(14.0), float(1.0)]

cluster_dict["17.7019996643_arrows"] += cgo_arrow([-13.5,103.0,14.0], [-11.668,104.855,13.029], color="blue red", name="Arrows_17.7019996643_2")

cluster_dict["17.7019996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-10.0), float(96.5), float(15.0), float(1.0)]

cluster_dict["17.7019996643_arrows"] += cgo_arrow([-10.0,96.5,15.0], [-8.732,93.236,15.72], color="blue red", name="Arrows_17.7019996643_3")

cluster_dict["17.7019996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-11.0), float(101.5), float(15.0), float(1.0)]

cluster_dict["17.7019996643_arrows"] += cgo_arrow([-11.0,101.5,15.0], [-8.147,102.532,14.595], color="blue red", name="Arrows_17.7019996643_4")

cluster_dict["17.7019996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-10.0), float(96.5), float(15.0), float(1.0)]

cluster_dict["17.7019996643_arrows"] += cgo_arrow([-10.0,96.5,15.0], [-8.732,93.236,15.72], color="blue red", name="Arrows_17.7019996643_5")

cluster_dict["17.7019996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-9.5), float(96.0), float(22.0), float(1.0)]

cluster_dict["17.7019996643_arrows"] += cgo_arrow([-9.5,96.0,22.0], [-8.429,98.737,22.114], color="blue red", name="Arrows_17.7019996643_6")

cluster_dict["17.7019996643"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-12.5368098043), float(95.2370322277), float(20.7519997994), float(1.0)]


cluster_dict["17.7019996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-13.5), float(98.0), float(15.0), float(1.0)]

cluster_dict["17.7019996643_arrows"] += cgo_arrow([-13.5,98.0,15.0], [-11.776,98.791,11.921], color="red blue", name="Arrows_17.7019996643_7")

cluster_dict["17.7019996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-10.5), float(95.5), float(15.5), float(1.0)]

cluster_dict["17.7019996643_arrows"] += cgo_arrow([-10.5,95.5,15.5], [-8.732,93.236,15.72], color="red blue", name="Arrows_17.7019996643_8")

cluster_dict["17.7019996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-10.5), float(95.5), float(15.5), float(1.0)]

cluster_dict["17.7019996643_arrows"] += cgo_arrow([-10.5,95.5,15.5], [-8.732,93.236,15.72], color="red blue", name="Arrows_17.7019996643_9")

cluster_dict["17.7019996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-9.5), float(99.5), float(19.0), float(1.0)]

cluster_dict["17.7019996643_arrows"] += cgo_arrow([-9.5,99.5,19.0], [-8.429,98.737,22.114], color="red blue", name="Arrows_17.7019996643_10")

cmd.load_cgo(cluster_dict["17.7019996643"], "Features_17.7019996643", 1)
cmd.load_cgo(cluster_dict["17.7019996643_arrows"], "Arrows_17.7019996643")
cmd.set("transparency", 0.2,"Features_17.7019996643")
cmd.group("Pharmacophore_17.7019996643", members="Features_17.7019996643")
cmd.group("Pharmacophore_17.7019996643", members="Arrows_17.7019996643")

if dirpath:
    f = join(dirpath, "label_threshold_17.7019996643.mol2")
else:
    f = "label_threshold_17.7019996643.mol2"

cmd.load(f, 'label_threshold_17.7019996643')
cmd.hide('everything', 'label_threshold_17.7019996643')
cmd.label("label_threshold_17.7019996643", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.7019996643', members= 'label_threshold_17.7019996643')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")