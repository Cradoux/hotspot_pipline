
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


cluster_dict = {"17.263999939":[], "17.263999939_arrows":[]}

cluster_dict["17.263999939"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-1.5), float(36.0), float(14.0), float(1.0)]

cluster_dict["17.263999939_arrows"] += cgo_arrow([-1.5,36.0,14.0], [-3.802,37.699,13.646], color="blue red", name="Arrows_17.263999939_1")

cluster_dict["17.263999939"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(1.5), float(37.0), float(11.5), float(1.0)]

cluster_dict["17.263999939_arrows"] += cgo_arrow([1.5,37.0,11.5], [-0.422,38.866,12.703], color="blue red", name="Arrows_17.263999939_2")

cluster_dict["17.263999939"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(2.5), float(28.5), float(19.5), float(1.0)]

cluster_dict["17.263999939_arrows"] += cgo_arrow([2.5,28.5,19.5], [5.259,27.287,18.789], color="blue red", name="Arrows_17.263999939_3")

cluster_dict["17.263999939"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(2.0), float(40.0), float(14.0), float(1.0)]

cluster_dict["17.263999939_arrows"] += cgo_arrow([2.0,40.0,14.0], [-0.422,38.866,12.703], color="blue red", name="Arrows_17.263999939_4")

cluster_dict["17.263999939"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(4.0), float(32.5), float(15.5), float(1.0)]

cluster_dict["17.263999939_arrows"] += cgo_arrow([4.0,32.5,15.5], [6.051,31.184,13.663], color="blue red", name="Arrows_17.263999939_5")

cluster_dict["17.263999939"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(2.4253272141), float(43.5868161354), float(15.6378714123), float(1.0)]


cluster_dict["17.263999939"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(0.871109226087), float(32.9749310419), float(15.663272766), float(1.0)]


cluster_dict["17.263999939"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-1.5), float(32.0), float(16.0), float(1.0)]


cluster_dict["17.263999939"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(2.0), float(37.5), float(16.0), float(1.0)]

cluster_dict["17.263999939_arrows"] += cgo_arrow([2.0,37.5,16.0], [2.014,39.461,18.603], color="red blue", name="Arrows_17.263999939_6")

cluster_dict["17.263999939"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(2.5), float(35.5), float(16.0), float(1.0)]

cluster_dict["17.263999939_arrows"] += cgo_arrow([2.5,35.5,16.0], [2.757,36.169,20.371], color="red blue", name="Arrows_17.263999939_7")

cluster_dict["17.263999939"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(3.0), float(29.0), float(15.5), float(1.0)]

cluster_dict["17.263999939_arrows"] += cgo_arrow([3.0,29.0,15.5], [6.367,28.535,16.294], color="red blue", name="Arrows_17.263999939_8")

cmd.load_cgo(cluster_dict["17.263999939"], "Features_17.263999939", 1)
cmd.load_cgo(cluster_dict["17.263999939_arrows"], "Arrows_17.263999939")
cmd.set("transparency", 0.2,"Features_17.263999939")
cmd.group("Pharmacophore_17.263999939", members="Features_17.263999939")
cmd.group("Pharmacophore_17.263999939", members="Arrows_17.263999939")

if dirpath:
    f = join(dirpath, "label_threshold_17.263999939.mol2")
else:
    f = "label_threshold_17.263999939.mol2"

cmd.load(f, 'label_threshold_17.263999939')
cmd.hide('everything', 'label_threshold_17.263999939')
cmd.label("label_threshold_17.263999939", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.263999939', members= 'label_threshold_17.263999939')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
