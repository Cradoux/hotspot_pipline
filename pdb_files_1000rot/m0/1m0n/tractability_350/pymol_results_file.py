
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


cluster_dict = {"12.625":[], "12.625_arrows":[]}

cluster_dict["12.625"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(13.0), float(54.0), float(24.5), float(1.0)]

cluster_dict["12.625_arrows"] += cgo_arrow([13.0,54.0,24.5], [11.45,52.415,22.116], color="blue red", name="Arrows_12.625_1")

cluster_dict["12.625"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(14.5), float(53.5), float(22.5), float(1.0)]

cluster_dict["12.625_arrows"] += cgo_arrow([14.5,53.5,22.5], [14.391,54.985,20.064], color="blue red", name="Arrows_12.625_2")

cluster_dict["12.625"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(16.9137687618), float(48.9970161262), float(24.3236572521), float(1.0)]


cluster_dict["12.625"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(13.0), float(47.5), float(26.5), float(1.0)]

cluster_dict["12.625_arrows"] += cgo_arrow([13.0,47.5,26.5], [13.665,46.492,23.997], color="red blue", name="Arrows_12.625_3")

cluster_dict["12.625"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(14.5), float(52.0), float(19.0), float(1.0)]

cluster_dict["12.625_arrows"] += cgo_arrow([14.5,52.0,19.0], [12.54,51.945,16.891], color="red blue", name="Arrows_12.625_4")

cluster_dict["12.625"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(14.0), float(55.0), float(22.5), float(1.0)]

cluster_dict["12.625_arrows"] += cgo_arrow([14.0,55.0,22.5], [14.391,54.985,20.064], color="red blue", name="Arrows_12.625_5")

cluster_dict["12.625"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(14.5), float(52.0), float(19.0), float(1.0)]

cluster_dict["12.625_arrows"] += cgo_arrow([14.5,52.0,19.0], [12.54,51.945,16.891], color="red blue", name="Arrows_12.625_6")

cluster_dict["12.625"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.5), float(53.5), float(21.5), float(1.0)]

cluster_dict["12.625_arrows"] += cgo_arrow([15.5,53.5,21.5], [14.391,54.985,20.064], color="red blue", name="Arrows_12.625_7")

cluster_dict["12.625"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.0), float(45.5), float(27.5), float(1.0)]

cluster_dict["12.625_arrows"] += cgo_arrow([17.0,45.5,27.5], [19.435,44.038,28.278], color="red blue", name="Arrows_12.625_8")

cmd.load_cgo(cluster_dict["12.625"], "Features_12.625", 1)
cmd.load_cgo(cluster_dict["12.625_arrows"], "Arrows_12.625")
cmd.set("transparency", 0.2,"Features_12.625")
cmd.group("Pharmacophore_12.625", members="Features_12.625")
cmd.group("Pharmacophore_12.625", members="Arrows_12.625")

if dirpath:
    f = join(dirpath, "label_threshold_12.625.mol2")
else:
    f = "label_threshold_12.625.mol2"

cmd.load(f, 'label_threshold_12.625')
cmd.hide('everything', 'label_threshold_12.625')
cmd.label("label_threshold_12.625", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.625', members= 'label_threshold_12.625')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")