
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


cluster_dict = {"16.9950008392":[], "16.9950008392_arrows":[]}

cluster_dict["16.9950008392"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-4.5), float(4.5), float(34.5), float(1.0)]

cluster_dict["16.9950008392_arrows"] += cgo_arrow([-4.5,4.5,34.5], [-5.806,3.277,33.264], color="blue red", name="Arrows_16.9950008392_1")

cluster_dict["16.9950008392"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(4.0), float(5.5), float(32.5), float(1.0)]

cluster_dict["16.9950008392_arrows"] += cgo_arrow([4.0,5.5,32.5], [4.428,7.941,30.839], color="blue red", name="Arrows_16.9950008392_2")

cluster_dict["16.9950008392"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(3.2886519388), float(4.85401211308), float(33.8303639535), float(1.0)]


cluster_dict["16.9950008392"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-5.5), float(5.5), float(37.0), float(1.0)]

cluster_dict["16.9950008392_arrows"] += cgo_arrow([-5.5,5.5,37.0], [-7.684,3.077,36.685], color="red blue", name="Arrows_16.9950008392_3")

cluster_dict["16.9950008392"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(1.5), float(5.5), float(34.5), float(1.0)]

cluster_dict["16.9950008392_arrows"] += cgo_arrow([1.5,5.5,34.5], [4.099,4.885,36.872], color="red blue", name="Arrows_16.9950008392_4")

cluster_dict["16.9950008392"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(0.0), float(5.0), float(32.5), float(1.0)]

cluster_dict["16.9950008392_arrows"] += cgo_arrow([0.0,5.0,32.5], [-1.528,3.649,30.394], color="red blue", name="Arrows_16.9950008392_5")

cluster_dict["16.9950008392"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(2.0), float(9.0), float(33.5), float(1.0)]

cluster_dict["16.9950008392_arrows"] += cgo_arrow([2.0,9.0,33.5], [3.428,11.483,32.678], color="red blue", name="Arrows_16.9950008392_6")

cluster_dict["16.9950008392"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(5.5), float(0.5), float(36.0), float(1.0)]

cluster_dict["16.9950008392_arrows"] += cgo_arrow([5.5,0.5,36.0], [3.023,1.895,38.655], color="red blue", name="Arrows_16.9950008392_7")

cluster_dict["16.9950008392"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(8.0), float(-0.5), float(34.0), float(1.0)]


cluster_dict["16.9950008392"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(11.0), float(1.5), float(34.0), float(1.0)]

cluster_dict["16.9950008392_arrows"] += cgo_arrow([11.0,1.5,34.0], [13.026,0.93,37.598], color="red blue", name="Arrows_16.9950008392_8")

cmd.load_cgo(cluster_dict["16.9950008392"], "Features_16.9950008392", 1)
cmd.load_cgo(cluster_dict["16.9950008392_arrows"], "Arrows_16.9950008392")
cmd.set("transparency", 0.2,"Features_16.9950008392")
cmd.group("Pharmacophore_16.9950008392", members="Features_16.9950008392")
cmd.group("Pharmacophore_16.9950008392", members="Arrows_16.9950008392")

if dirpath:
    f = join(dirpath, "label_threshold_16.9950008392.mol2")
else:
    f = "label_threshold_16.9950008392.mol2"

cmd.load(f, 'label_threshold_16.9950008392')
cmd.hide('everything', 'label_threshold_16.9950008392')
cmd.label("label_threshold_16.9950008392", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.9950008392', members= 'label_threshold_16.9950008392')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
