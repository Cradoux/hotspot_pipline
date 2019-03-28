
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


cluster_dict = {"15.0399999619":[], "15.0399999619_arrows":[]}

cluster_dict["15.0399999619"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-8.5), float(18.0), float(34.0), float(1.0)]

cluster_dict["15.0399999619_arrows"] += cgo_arrow([-8.5,18.0,34.0], [-9.569,19.744,36.273], color="blue red", name="Arrows_15.0399999619_1")

cluster_dict["15.0399999619"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-5.5), float(20.5), float(25.0), float(1.0)]

cluster_dict["15.0399999619_arrows"] += cgo_arrow([-5.5,20.5,25.0], [-8.164,21.467,24.804], color="blue red", name="Arrows_15.0399999619_2")

cluster_dict["15.0399999619"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-9.53177467605), float(15.3820218023), float(28.7148559129), float(1.0)]


cluster_dict["15.0399999619"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-14.0), float(12.0), float(29.5), float(1.0)]

cluster_dict["15.0399999619_arrows"] += cgo_arrow([-14.0,12.0,29.5], [-16.374,9.898,26.209], color="red blue", name="Arrows_15.0399999619_3")

cluster_dict["15.0399999619"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-12.0), float(15.5), float(33.5), float(1.0)]

cluster_dict["15.0399999619_arrows"] += cgo_arrow([-12.0,15.5,33.5], [-13.393,18.264,34.841], color="red blue", name="Arrows_15.0399999619_4")

cluster_dict["15.0399999619"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-11.0), float(12.0), float(35.5), float(1.0)]

cluster_dict["15.0399999619_arrows"] += cgo_arrow([-11.0,12.0,35.5], [-9.078,9.709,34.737], color="red blue", name="Arrows_15.0399999619_5")

cluster_dict["15.0399999619"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-10.5), float(16.5), float(23.5), float(1.0)]

cluster_dict["15.0399999619_arrows"] += cgo_arrow([-10.5,16.5,23.5], [-13.176,18.246,26.384], color="red blue", name="Arrows_15.0399999619_6")

cluster_dict["15.0399999619"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-9.5), float(17.5), float(30.0), float(1.0)]

cluster_dict["15.0399999619_arrows"] += cgo_arrow([-9.5,17.5,30.0], [-9.296,20.715,29.569], color="red blue", name="Arrows_15.0399999619_7")

cluster_dict["15.0399999619"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-8.0), float(14.0), float(28.5), float(1.0)]

cluster_dict["15.0399999619_arrows"] += cgo_arrow([-8.0,14.0,28.5], [-8.616,10.95,27.786], color="red blue", name="Arrows_15.0399999619_8")

cluster_dict["15.0399999619"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-3.0), float(19.0), float(24.5), float(1.0)]

cluster_dict["15.0399999619_arrows"] += cgo_arrow([-3.0,19.0,24.5], [-3.494,21.076,22.873], color="red blue", name="Arrows_15.0399999619_9")

cmd.load_cgo(cluster_dict["15.0399999619"], "Features_15.0399999619", 1)
cmd.load_cgo(cluster_dict["15.0399999619_arrows"], "Arrows_15.0399999619")
cmd.set("transparency", 0.2,"Features_15.0399999619")
cmd.group("Pharmacophore_15.0399999619", members="Features_15.0399999619")
cmd.group("Pharmacophore_15.0399999619", members="Arrows_15.0399999619")

if dirpath:
    f = join(dirpath, "label_threshold_15.0399999619.mol2")
else:
    f = "label_threshold_15.0399999619.mol2"

cmd.load(f, 'label_threshold_15.0399999619')
cmd.hide('everything', 'label_threshold_15.0399999619')
cmd.label("label_threshold_15.0399999619", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.0399999619', members= 'label_threshold_15.0399999619')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
