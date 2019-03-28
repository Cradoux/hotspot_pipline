
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


cluster_dict = {"14.045999527":[], "14.045999527_arrows":[]}

cluster_dict["14.045999527"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(57.0), float(24.0), float(1.0)]

cluster_dict["14.045999527_arrows"] += cgo_arrow([22.0,57.0,24.0], [22.472,54.511,23.015], color="blue red", name="Arrows_14.045999527_1")

cluster_dict["14.045999527"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(57.0), float(24.0), float(1.0)]

cluster_dict["14.045999527_arrows"] += cgo_arrow([22.0,57.0,24.0], [22.472,54.511,23.015], color="blue red", name="Arrows_14.045999527_2")

cluster_dict["14.045999527"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(56.0), float(28.0), float(1.0)]

cluster_dict["14.045999527_arrows"] += cgo_arrow([23.5,56.0,28.0], [22.078,55.253,29.415], color="blue red", name="Arrows_14.045999527_3")

cluster_dict["14.045999527"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(28.5), float(57.0), float(22.0), float(1.0)]

cluster_dict["14.045999527_arrows"] += cgo_arrow([28.5,57.0,22.0], [30.436,55.621,20.697], color="blue red", name="Arrows_14.045999527_4")

cluster_dict["14.045999527"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(14.6833207163), float(57.025024804), float(23.1583482444), float(1.0)]


cluster_dict["14.045999527"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.115352679), float(58.7387929468), float(24.8935365677), float(1.0)]


cluster_dict["14.045999527"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(57.0), float(22.5), float(1.0)]

cluster_dict["14.045999527_arrows"] += cgo_arrow([18.5,57.0,22.5], [18.518,57.133,19.453], color="red blue", name="Arrows_14.045999527_5")

cluster_dict["14.045999527"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(57.0), float(26.5), float(1.0)]

cluster_dict["14.045999527_arrows"] += cgo_arrow([26.5,57.0,26.5], [27.568,55.302,26.456], color="red blue", name="Arrows_14.045999527_6")

cluster_dict["14.045999527"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(25.5), float(61.0), float(27.0), float(1.0)]

cluster_dict["14.045999527_arrows"] += cgo_arrow([25.5,61.0,27.0], [25.159,63.146,28.81], color="red blue", name="Arrows_14.045999527_7")

cluster_dict["14.045999527"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(57.0), float(22.5), float(1.0)]

cluster_dict["14.045999527_arrows"] += cgo_arrow([26.5,57.0,22.5], [24.547,55.693,21.198], color="red blue", name="Arrows_14.045999527_8")

cluster_dict["14.045999527"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(58.5), float(21.5), float(1.0)]

cluster_dict["14.045999527_arrows"] += cgo_arrow([29.5,58.5,21.5], [30.942,58.212,18.848], color="red blue", name="Arrows_14.045999527_9")

cluster_dict["14.045999527"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(57.5), float(26.0), float(1.0)]

cluster_dict["14.045999527_arrows"] += cgo_arrow([29.5,57.5,26.0], [27.568,55.302,26.456], color="red blue", name="Arrows_14.045999527_10")

cmd.load_cgo(cluster_dict["14.045999527"], "Features_14.045999527", 1)
cmd.load_cgo(cluster_dict["14.045999527_arrows"], "Arrows_14.045999527")
cmd.set("transparency", 0.2,"Features_14.045999527")
cmd.group("Pharmacophore_14.045999527", members="Features_14.045999527")
cmd.group("Pharmacophore_14.045999527", members="Arrows_14.045999527")

if dirpath:
    f = join(dirpath, "label_threshold_14.045999527.mol2")
else:
    f = "label_threshold_14.045999527.mol2"

cmd.load(f, 'label_threshold_14.045999527')
cmd.hide('everything', 'label_threshold_14.045999527')
cmd.label("label_threshold_14.045999527", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.045999527', members= 'label_threshold_14.045999527')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
