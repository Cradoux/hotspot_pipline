
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


cluster_dict = {"13.8879995346":[], "13.8879995346_arrows":[]}

cluster_dict["13.8879995346"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(14.0), float(9.0), float(27.5), float(1.0)]

cluster_dict["13.8879995346_arrows"] += cgo_arrow([14.0,9.0,27.5], [15.449,10.37,29.635], color="blue red", name="Arrows_13.8879995346_1")

cluster_dict["13.8879995346"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(8.5), float(25.0), float(1.0)]

cluster_dict["13.8879995346_arrows"] += cgo_arrow([22.0,8.5,25.0], [21.083,5.903,23.79], color="blue red", name="Arrows_13.8879995346_2")

cluster_dict["13.8879995346"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(20.9631524649), float(11.0128669722), float(24.3907292261), float(1.0)]


cluster_dict["13.8879995346"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(14.0), float(8.0), float(28.5), float(1.0)]

cluster_dict["13.8879995346_arrows"] += cgo_arrow([14.0,8.0,28.5], [13.134,6.322,30.406], color="red blue", name="Arrows_13.8879995346_3")

cluster_dict["13.8879995346"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(8.5), float(22.5), float(1.0)]

cluster_dict["13.8879995346_arrows"] += cgo_arrow([17.5,8.5,22.5], [18.417,6.697,20.28], color="red blue", name="Arrows_13.8879995346_4")

cluster_dict["13.8879995346"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(6.5), float(25.5), float(1.0)]

cluster_dict["13.8879995346_arrows"] += cgo_arrow([24.5,6.5,25.5], [23.127,3.857,24.382], color="red blue", name="Arrows_13.8879995346_5")

cluster_dict["13.8879995346"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(9.0), float(22.0), float(1.0)]

cluster_dict["13.8879995346_arrows"] += cgo_arrow([24.5,9.0,22.0], [25.442,7.456,19.567], color="red blue", name="Arrows_13.8879995346_6")

cmd.load_cgo(cluster_dict["13.8879995346"], "Features_13.8879995346", 1)
cmd.load_cgo(cluster_dict["13.8879995346_arrows"], "Arrows_13.8879995346")
cmd.set("transparency", 0.2,"Features_13.8879995346")
cmd.group("Pharmacophore_13.8879995346", members="Features_13.8879995346")
cmd.group("Pharmacophore_13.8879995346", members="Arrows_13.8879995346")

if dirpath:
    f = join(dirpath, "label_threshold_13.8879995346.mol2")
else:
    f = "label_threshold_13.8879995346.mol2"

cmd.load(f, 'label_threshold_13.8879995346')
cmd.hide('everything', 'label_threshold_13.8879995346')
cmd.label("label_threshold_13.8879995346", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.8879995346', members= 'label_threshold_13.8879995346')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")