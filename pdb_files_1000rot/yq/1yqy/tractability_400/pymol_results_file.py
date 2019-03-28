
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


cluster_dict = {"14.1499996185":[], "14.1499996185_arrows":[]}

cluster_dict["14.1499996185"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(7.0), float(12.0), float(27.0), float(1.0)]

cluster_dict["14.1499996185_arrows"] += cgo_arrow([7.0,12.0,27.0], [5.78,10.478,24.794], color="blue red", name="Arrows_14.1499996185_1")

cluster_dict["14.1499996185"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(13.5), float(18.0), float(25.5), float(1.0)]

cluster_dict["14.1499996185_arrows"] += cgo_arrow([13.5,18.0,25.5], [12.721,19.487,23.066], color="blue red", name="Arrows_14.1499996185_2")

cluster_dict["14.1499996185"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(13.0), float(18.5), float(29.5), float(1.0)]

cluster_dict["14.1499996185_arrows"] += cgo_arrow([13.0,18.5,29.5], [11.952,19.883,31.578], color="blue red", name="Arrows_14.1499996185_3")

cluster_dict["14.1499996185"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(12.5), float(11.0), float(24.5), float(1.0)]

cluster_dict["14.1499996185_arrows"] += cgo_arrow([12.5,11.0,24.5], [10.052,12.521,24.252], color="blue red", name="Arrows_14.1499996185_4")

cluster_dict["14.1499996185"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.0), float(15.5), float(22.5), float(1.0)]

cluster_dict["14.1499996185_arrows"] += cgo_arrow([19.0,15.5,22.5], [19.446,13.652,20.202], color="blue red", name="Arrows_14.1499996185_5")

cluster_dict["14.1499996185"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(7.26089857836), float(11.6483063088), float(28.5799135927), float(1.0)]


cluster_dict["14.1499996185"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(6.08839038419), float(13.3666586788), float(26.3490225513), float(1.0)]


cluster_dict["14.1499996185"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(12.7580758279), float(16.3682634974), float(26.4969574601), float(1.0)]


cluster_dict["14.1499996185"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(7.5), float(19.5), float(26.5), float(1.0)]

cluster_dict["14.1499996185_arrows"] += cgo_arrow([7.5,19.5,26.5], [4.514,19.97,26.919], color="red blue", name="Arrows_14.1499996185_6")

cluster_dict["14.1499996185"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(13.5), float(17.0), float(25.5), float(1.0)]

cluster_dict["14.1499996185_arrows"] += cgo_arrow([13.5,17.0,25.5], [12.721,19.487,23.066], color="red blue", name="Arrows_14.1499996185_7")

cluster_dict["14.1499996185"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.5), float(11.5), float(23.0), float(1.0)]

cluster_dict["14.1499996185_arrows"] += cgo_arrow([16.5,11.5,23.0], [16.579,8.495,22.858], color="red blue", name="Arrows_14.1499996185_8")

cluster_dict["14.1499996185"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(18.0), float(25.5), float(1.0)]

cluster_dict["14.1499996185_arrows"] += cgo_arrow([21.0,18.0,25.5], [21.277,21.713,27.522], color="red blue", name="Arrows_14.1499996185_9")

cmd.load_cgo(cluster_dict["14.1499996185"], "Features_14.1499996185", 1)
cmd.load_cgo(cluster_dict["14.1499996185_arrows"], "Arrows_14.1499996185")
cmd.set("transparency", 0.2,"Features_14.1499996185")
cmd.group("Pharmacophore_14.1499996185", members="Features_14.1499996185")
cmd.group("Pharmacophore_14.1499996185", members="Arrows_14.1499996185")

if dirpath:
    f = join(dirpath, "label_threshold_14.1499996185.mol2")
else:
    f = "label_threshold_14.1499996185.mol2"

cmd.load(f, 'label_threshold_14.1499996185')
cmd.hide('everything', 'label_threshold_14.1499996185')
cmd.label("label_threshold_14.1499996185", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.1499996185', members= 'label_threshold_14.1499996185')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
