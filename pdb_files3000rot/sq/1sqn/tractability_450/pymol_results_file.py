
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


cluster_dict = {"16.3269996643":[], "16.3269996643_arrows":[]}

cluster_dict["16.3269996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-5.5), float(2.0), float(29.0), float(1.0)]

cluster_dict["16.3269996643_arrows"] += cgo_arrow([-5.5,2.0,29.0], [-7.046,3.455,26.953], color="blue red", name="Arrows_16.3269996643_1")

cluster_dict["16.3269996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-4.0), float(5.5), float(21.0), float(1.0)]

cluster_dict["16.3269996643_arrows"] += cgo_arrow([-4.0,5.5,21.0], [-5.708,5.87,23.416], color="blue red", name="Arrows_16.3269996643_2")

cluster_dict["16.3269996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-2.5), float(-0.5), float(30.5), float(1.0)]

cluster_dict["16.3269996643_arrows"] += cgo_arrow([-2.5,-0.5,30.5], [-0.137,-1.469,31.664], color="blue red", name="Arrows_16.3269996643_3")

cluster_dict["16.3269996643"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-2.46616675421), float(1.54456917018), float(25.1877161744), float(1.0)]


cluster_dict["16.3269996643"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-4.94494357133), float(4.61071967524), float(19.5), float(1.0)]


cluster_dict["16.3269996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-5.0), float(-2.0), float(30.5), float(1.0)]

cluster_dict["16.3269996643_arrows"] += cgo_arrow([-5.0,-2.0,30.5], [-4.293,-3.037,32.774], color="red blue", name="Arrows_16.3269996643_4")

cluster_dict["16.3269996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-4.5), float(0.0), float(24.5), float(1.0)]


cluster_dict["16.3269996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-2.5), float(5.5), float(20.5), float(1.0)]

cluster_dict["16.3269996643_arrows"] += cgo_arrow([-2.5,5.5,20.5], [-2.48,5.061,16.995], color="red blue", name="Arrows_16.3269996643_5")

cluster_dict["16.3269996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-0.5), float(0.0), float(23.0), float(1.0)]


cmd.load_cgo(cluster_dict["16.3269996643"], "Features_16.3269996643", 1)
cmd.load_cgo(cluster_dict["16.3269996643_arrows"], "Arrows_16.3269996643")
cmd.set("transparency", 0.2,"Features_16.3269996643")
cmd.group("Pharmacophore_16.3269996643", members="Features_16.3269996643")
cmd.group("Pharmacophore_16.3269996643", members="Arrows_16.3269996643")

if dirpath:
    f = join(dirpath, "label_threshold_16.3269996643.mol2")
else:
    f = "label_threshold_16.3269996643.mol2"

cmd.load(f, 'label_threshold_16.3269996643')
cmd.hide('everything', 'label_threshold_16.3269996643')
cmd.label("label_threshold_16.3269996643", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.3269996643', members= 'label_threshold_16.3269996643')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
