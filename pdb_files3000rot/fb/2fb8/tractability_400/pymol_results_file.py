
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


cluster_dict = {"17.1990003586":[], "17.1990003586_arrows":[]}

cluster_dict["17.1990003586"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-44.0), float(28.5), float(-4.5), float(1.0)]

cluster_dict["17.1990003586_arrows"] += cgo_arrow([-44.0,28.5,-4.5], [-45.62,26.402,-4.591], color="blue red", name="Arrows_17.1990003586_1")

cluster_dict["17.1990003586"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(24.0), float(-0.5), float(1.0)]

cluster_dict["17.1990003586_arrows"] += cgo_arrow([-37.0,24.0,-0.5], [-35.435,21.826,0.786], color="blue red", name="Arrows_17.1990003586_2")

cluster_dict["17.1990003586"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-34.5), float(29.0), float(-1.0), float(1.0)]

cluster_dict["17.1990003586_arrows"] += cgo_arrow([-34.5,29.0,-1.0], [-32.426,27.812,-0.09], color="blue red", name="Arrows_17.1990003586_3")

cluster_dict["17.1990003586"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-39.3852470926), float(27.4683104269), float(-0.700966615353), float(1.0)]


cluster_dict["17.1990003586"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-44.5), float(23.0), float(-2.5), float(1.0)]

cluster_dict["17.1990003586_arrows"] += cgo_arrow([-44.5,23.0,-2.5], [-47.478,23.36,-2.426], color="red blue", name="Arrows_17.1990003586_4")

cluster_dict["17.1990003586"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-40.5), float(22.5), float(0.0), float(1.0)]

cluster_dict["17.1990003586_arrows"] += cgo_arrow([-40.5,22.5,0.0], [-38.239,21.157,0.384], color="red blue", name="Arrows_17.1990003586_5")

cluster_dict["17.1990003586"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-36.0), float(25.5), float(0.5), float(1.0)]

cluster_dict["17.1990003586_arrows"] += cgo_arrow([-36.0,25.5,0.5], [-32.936,24.781,0.433], color="red blue", name="Arrows_17.1990003586_6")

cmd.load_cgo(cluster_dict["17.1990003586"], "Features_17.1990003586", 1)
cmd.load_cgo(cluster_dict["17.1990003586_arrows"], "Arrows_17.1990003586")
cmd.set("transparency", 0.2,"Features_17.1990003586")
cmd.group("Pharmacophore_17.1990003586", members="Features_17.1990003586")
cmd.group("Pharmacophore_17.1990003586", members="Arrows_17.1990003586")

if dirpath:
    f = join(dirpath, "label_threshold_17.1990003586.mol2")
else:
    f = "label_threshold_17.1990003586.mol2"

cmd.load(f, 'label_threshold_17.1990003586')
cmd.hide('everything', 'label_threshold_17.1990003586')
cmd.label("label_threshold_17.1990003586", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.1990003586', members= 'label_threshold_17.1990003586')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
