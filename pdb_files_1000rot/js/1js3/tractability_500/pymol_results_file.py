
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


cluster_dict = {"13.531999588":[], "13.531999588_arrows":[]}

cluster_dict["13.531999588"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(54.0), float(14.5), float(86.5), float(1.0)]

cluster_dict["13.531999588_arrows"] += cgo_arrow([54.0,14.5,86.5], [56.522,13.369,87.687], color="blue red", name="Arrows_13.531999588_1")

cluster_dict["13.531999588"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(54.0), float(17.0), float(81.5), float(1.0)]

cluster_dict["13.531999588_arrows"] += cgo_arrow([54.0,17.0,81.5], [56.376,16.304,83.272], color="blue red", name="Arrows_13.531999588_2")

cluster_dict["13.531999588"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(55.0), float(18.5), float(86.5), float(1.0)]

cluster_dict["13.531999588_arrows"] += cgo_arrow([55.0,18.5,86.5], [56.042,20.207,88.973], color="blue red", name="Arrows_13.531999588_3")

cluster_dict["13.531999588"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(55.0), float(17.0), float(88.0), float(1.0)]

cluster_dict["13.531999588_arrows"] += cgo_arrow([55.0,17.0,88.0], [56.042,20.207,88.973], color="blue red", name="Arrows_13.531999588_4")

cluster_dict["13.531999588"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(53.1814502378), float(13.7075281988), float(82.9932178369), float(1.0)]


cluster_dict["13.531999588"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(50.5), float(11.0), float(83.5), float(1.0)]

cluster_dict["13.531999588_arrows"] += cgo_arrow([50.5,11.0,83.5], [50.505,8.932,85.935], color="red blue", name="Arrows_13.531999588_5")

cluster_dict["13.531999588"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(53.5), float(20.0), float(85.0), float(1.0)]

cluster_dict["13.531999588_arrows"] += cgo_arrow([53.5,20.0,85.0], [53.076,22.34,83.177], color="red blue", name="Arrows_13.531999588_6")

cluster_dict["13.531999588"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(55.5), float(9.0), float(84.5), float(1.0)]

cluster_dict["13.531999588_arrows"] += cgo_arrow([55.5,9.0,84.5], [54.183,6.069,87.338], color="red blue", name="Arrows_13.531999588_7")

cmd.load_cgo(cluster_dict["13.531999588"], "Features_13.531999588", 1)
cmd.load_cgo(cluster_dict["13.531999588_arrows"], "Arrows_13.531999588")
cmd.set("transparency", 0.2,"Features_13.531999588")
cmd.group("Pharmacophore_13.531999588", members="Features_13.531999588")
cmd.group("Pharmacophore_13.531999588", members="Arrows_13.531999588")

if dirpath:
    f = join(dirpath, "label_threshold_13.531999588.mol2")
else:
    f = "label_threshold_13.531999588.mol2"

cmd.load(f, 'label_threshold_13.531999588')
cmd.hide('everything', 'label_threshold_13.531999588')
cmd.label("label_threshold_13.531999588", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.531999588', members= 'label_threshold_13.531999588')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
