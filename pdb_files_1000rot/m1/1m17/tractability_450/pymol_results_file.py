
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


cluster_dict = {"16.375":[], "16.375_arrows":[]}

cluster_dict["16.375"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(20.0), float(0.0), float(54.5), float(1.0)]

cluster_dict["16.375_arrows"] += cgo_arrow([20.0,0.0,54.5], [18.565,-1.027,56.748], color="blue red", name="Arrows_16.375_1")

cluster_dict["16.375"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(-1.0), float(56.5), float(1.0)]

cluster_dict["16.375_arrows"] += cgo_arrow([24.5,-1.0,56.5], [24.253,-3.657,58.017], color="blue red", name="Arrows_16.375_2")

cluster_dict["16.375"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(2.5), float(52.5), float(1.0)]

cluster_dict["16.375_arrows"] += cgo_arrow([26.0,2.5,52.5], [26.827,5.51,53.247], color="blue red", name="Arrows_16.375_3")

cluster_dict["16.375"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.5), float(2.5), float(53.5), float(1.0)]

cluster_dict["16.375_arrows"] += cgo_arrow([27.5,2.5,53.5], [26.827,5.51,53.247], color="blue red", name="Arrows_16.375_4")

cluster_dict["16.375"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(28.5), float(-2.5), float(53.0), float(1.0)]

cluster_dict["16.375_arrows"] += cgo_arrow([28.5,-2.5,53.0], [27.403,-3.109,55.889], color="blue red", name="Arrows_16.375_5")

cluster_dict["16.375"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(23.9302795016), float(0.817936798106), float(52.4842495676), float(1.0)]


cluster_dict["16.375"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(-2.5), float(55.0), float(1.0)]

cluster_dict["16.375_arrows"] += cgo_arrow([22.5,-2.5,55.0], [20.898,-2.866,57.744], color="red blue", name="Arrows_16.375_6")

cluster_dict["16.375"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(1.5), float(50.5), float(1.0)]


cluster_dict["16.375"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(25.5), float(3.0), float(54.0), float(1.0)]

cluster_dict["16.375_arrows"] += cgo_arrow([25.5,3.0,54.0], [26.827,5.51,53.247], color="red blue", name="Arrows_16.375_7")

cluster_dict["16.375"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(25.5), float(3.0), float(54.0), float(1.0)]

cluster_dict["16.375_arrows"] += cgo_arrow([25.5,3.0,54.0], [26.827,5.51,53.247], color="red blue", name="Arrows_16.375_8")

cmd.load_cgo(cluster_dict["16.375"], "Features_16.375", 1)
cmd.load_cgo(cluster_dict["16.375_arrows"], "Arrows_16.375")
cmd.set("transparency", 0.2,"Features_16.375")
cmd.group("Pharmacophore_16.375", members="Features_16.375")
cmd.group("Pharmacophore_16.375", members="Arrows_16.375")

if dirpath:
    f = join(dirpath, "label_threshold_16.375.mol2")
else:
    f = "label_threshold_16.375.mol2"

cmd.load(f, 'label_threshold_16.375')
cmd.hide('everything', 'label_threshold_16.375')
cmd.label("label_threshold_16.375", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.375', members= 'label_threshold_16.375')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
