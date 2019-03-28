
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


cluster_dict = {"16.5569992065":[], "16.5569992065_arrows":[]}

cluster_dict["16.5569992065"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(15.5), float(24.5), float(28.5), float(1.0)]

cluster_dict["16.5569992065_arrows"] += cgo_arrow([15.5,24.5,28.5], [15.481,23.448,26.221], color="blue red", name="Arrows_16.5569992065_1")

cluster_dict["16.5569992065"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(23.5), float(27.5), float(1.0)]

cluster_dict["16.5569992065_arrows"] += cgo_arrow([25.0,23.5,27.5], [23.825,23.679,30.554], color="blue red", name="Arrows_16.5569992065_2")

cluster_dict["16.5569992065"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(29.5), float(25.0), float(1.0)]

cluster_dict["16.5569992065_arrows"] += cgo_arrow([25.0,29.5,25.0], [24.804,32.148,26.296], color="blue red", name="Arrows_16.5569992065_3")

cluster_dict["16.5569992065"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(21.9339281169), float(26.3346046445), float(26.6842481614), float(1.0)]


cluster_dict["16.5569992065"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(26.0), float(29.0), float(1.0)]

cluster_dict["16.5569992065_arrows"] += cgo_arrow([16.0,26.0,29.0], [17.178,26.794,32.704], color="red blue", name="Arrows_16.5569992065_4")

cluster_dict["16.5569992065"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(28.0), float(25.5), float(1.0)]

cluster_dict["16.5569992065_arrows"] += cgo_arrow([17.5,28.0,25.5], [14.267,27.323,25.392], color="red blue", name="Arrows_16.5569992065_5")

cluster_dict["16.5569992065"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(28.0), float(29.5), float(1.0)]

cluster_dict["16.5569992065_arrows"] += cgo_arrow([18.5,28.0,29.5], [17.178,26.794,32.704], color="red blue", name="Arrows_16.5569992065_6")

cluster_dict["16.5569992065"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(19.5), float(22.0), float(27.0), float(1.0)]

cluster_dict["16.5569992065_arrows"] += cgo_arrow([19.5,22.0,27.0], [19.763,22.288,29.47], color="red blue", name="Arrows_16.5569992065_7")

cluster_dict["16.5569992065"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.0), float(31.0), float(28.5), float(1.0)]

cluster_dict["16.5569992065_arrows"] += cgo_arrow([20.0,31.0,28.5], [22.1,33.173,28.772], color="red blue", name="Arrows_16.5569992065_8")

cluster_dict["16.5569992065"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(26.5), float(23.0), float(1.0)]

cluster_dict["16.5569992065_arrows"] += cgo_arrow([22.5,26.5,23.0], [21.419,24.54,21.082], color="red blue", name="Arrows_16.5569992065_9")

cmd.load_cgo(cluster_dict["16.5569992065"], "Features_16.5569992065", 1)
cmd.load_cgo(cluster_dict["16.5569992065_arrows"], "Arrows_16.5569992065")
cmd.set("transparency", 0.2,"Features_16.5569992065")
cmd.group("Pharmacophore_16.5569992065", members="Features_16.5569992065")
cmd.group("Pharmacophore_16.5569992065", members="Arrows_16.5569992065")

if dirpath:
    f = join(dirpath, "label_threshold_16.5569992065.mol2")
else:
    f = "label_threshold_16.5569992065.mol2"

cmd.load(f, 'label_threshold_16.5569992065')
cmd.hide('everything', 'label_threshold_16.5569992065')
cmd.label("label_threshold_16.5569992065", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.5569992065', members= 'label_threshold_16.5569992065')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
