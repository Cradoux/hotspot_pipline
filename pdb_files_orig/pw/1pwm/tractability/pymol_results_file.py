
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


cluster_dict = {"15.968000412":[], "15.968000412_arrows":[]}

cluster_dict["15.968000412"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(15.5), float(-4.5), float(19.5), float(1.0)]

cluster_dict["15.968000412_arrows"] += cgo_arrow([15.5,-4.5,19.5], [12.907,-4.46,18.205], color="blue red", name="Arrows_15.968000412_1")

cluster_dict["15.968000412"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(-2.0), float(23.5), float(1.0)]

cluster_dict["15.968000412_arrows"] += cgo_arrow([16.0,-2.0,23.5], [14.255,-0.537,22.458], color="blue red", name="Arrows_15.968000412_2")

cluster_dict["15.968000412"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.0), float(-4.5), float(22.5), float(1.0)]

cluster_dict["15.968000412_arrows"] += cgo_arrow([19.0,-4.5,22.5], [21.795,-3.156,22.131], color="blue red", name="Arrows_15.968000412_3")

cluster_dict["15.968000412"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(17.3057757037), float(-7.30981903322), float(17.8615479205), float(1.0)]


cluster_dict["15.968000412"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(14.0784223155), float(-3.3170931798), float(20.6732607887), float(1.0)]


cluster_dict["15.968000412"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(14.3384988467), float(-2.28352298984), float(18.4242646564), float(1.0)]


cluster_dict["15.968000412"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(-7.5), float(16.5), float(1.0)]

cluster_dict["15.968000412_arrows"] += cgo_arrow([18.0,-7.5,16.5], [20.681,-4.461,18.07], color="red blue", name="Arrows_15.968000412_4")

cluster_dict["15.968000412"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(14.5), float(-4.0), float(23.0), float(1.0)]

cluster_dict["15.968000412_arrows"] += cgo_arrow([14.5,-4.0,23.0], [13.814,-5.448,21.846], color="red blue", name="Arrows_15.968000412_5")

cluster_dict["15.968000412"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(14.5), float(-5.5), float(19.5), float(1.0)]

cluster_dict["15.968000412_arrows"] += cgo_arrow([14.5,-5.5,19.5], [13.814,-5.448,21.846], color="red blue", name="Arrows_15.968000412_6")

cluster_dict["15.968000412"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(-2.5), float(18.0), float(1.0)]

cluster_dict["15.968000412_arrows"] += cgo_arrow([16.0,-2.5,18.0], [15.247,-0.45,15.978], color="red blue", name="Arrows_15.968000412_7")

cluster_dict["15.968000412"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(-7.5), float(16.5), float(1.0)]

cluster_dict["15.968000412_arrows"] += cgo_arrow([18.0,-7.5,16.5], [20.681,-4.461,18.07], color="red blue", name="Arrows_15.968000412_8")

cluster_dict["15.968000412"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(-8.5), float(13.5), float(1.0)]

cluster_dict["15.968000412_arrows"] += cgo_arrow([18.5,-8.5,13.5], [20.712,-4.579,11.589], color="red blue", name="Arrows_15.968000412_9")

cluster_dict["15.968000412"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.5), float(-6.0), float(14.0), float(1.0)]

cluster_dict["15.968000412_arrows"] += cgo_arrow([20.5,-6.0,14.0], [20.712,-4.579,11.589], color="red blue", name="Arrows_15.968000412_10")

cmd.load_cgo(cluster_dict["15.968000412"], "Features_15.968000412", 1)
cmd.load_cgo(cluster_dict["15.968000412_arrows"], "Arrows_15.968000412")
cmd.set("transparency", 0.2,"Features_15.968000412")
cmd.group("Pharmacophore_15.968000412", members="Features_15.968000412")
cmd.group("Pharmacophore_15.968000412", members="Arrows_15.968000412")

if dirpath:
    f = join(dirpath, "label_threshold_15.968000412.mol2")
else:
    f = "label_threshold_15.968000412.mol2"

cmd.load(f, 'label_threshold_15.968000412')
cmd.hide('everything', 'label_threshold_15.968000412')
cmd.label("label_threshold_15.968000412", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.968000412', members= 'label_threshold_15.968000412')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
