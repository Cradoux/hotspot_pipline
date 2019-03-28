
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


cluster_dict = {"13.5579996109":[], "13.5579996109_arrows":[]}

cluster_dict["13.5579996109"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(14.0), float(3.5), float(16.0), float(1.0)]

cluster_dict["13.5579996109_arrows"] += cgo_arrow([14.0,3.5,16.0], [13.465,1.988,13.526], color="blue red", name="Arrows_13.5579996109_1")

cluster_dict["13.5579996109"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(16.5), float(2.0), float(10.0), float(1.0)]

cluster_dict["13.5579996109_arrows"] += cgo_arrow([16.5,2.0,10.0], [14.173,3.186,8.735], color="blue red", name="Arrows_13.5579996109_2")

cluster_dict["13.5579996109"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(17.2102412213), float(4.97110437191), float(13.4032488438), float(1.0)]


cluster_dict["13.5579996109"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(4.0), float(16.0), float(1.0)]

cluster_dict["13.5579996109_arrows"] += cgo_arrow([18.0,4.0,16.0], [20.254,3.203,14.512], color="red blue", name="Arrows_13.5579996109_3")

cluster_dict["13.5579996109"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.5), float(3.5), float(12.0), float(1.0)]

cluster_dict["13.5579996109_arrows"] += cgo_arrow([15.5,3.5,12.0], [13.406,5.822,11.575], color="red blue", name="Arrows_13.5579996109_4")

cluster_dict["13.5579996109"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(8.0), float(15.0), float(1.0)]

cluster_dict["13.5579996109_arrows"] += cgo_arrow([16.0,8.0,15.0], [13.406,5.822,11.575], color="red blue", name="Arrows_13.5579996109_5")

cluster_dict["13.5579996109"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.0), float(1.0), float(15.0), float(1.0)]

cluster_dict["13.5579996109_arrows"] += cgo_arrow([17.0,1.0,15.0], [17.464,-0.224,17.579], color="red blue", name="Arrows_13.5579996109_6")

cluster_dict["13.5579996109"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(3.0), float(14.0), float(1.0)]

cluster_dict["13.5579996109_arrows"] += cgo_arrow([17.5,3.0,14.0], [20.254,3.203,14.512], color="red blue", name="Arrows_13.5579996109_7")

cluster_dict["13.5579996109"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(19.5), float(1.0), float(12.0), float(1.0)]

cluster_dict["13.5579996109_arrows"] += cgo_arrow([19.5,1.0,12.0], [20.254,3.203,14.512], color="red blue", name="Arrows_13.5579996109_8")

cluster_dict["13.5579996109"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(19.5), float(1.0), float(12.0), float(1.0)]

cluster_dict["13.5579996109_arrows"] += cgo_arrow([19.5,1.0,12.0], [20.254,3.203,14.512], color="red blue", name="Arrows_13.5579996109_9")

cmd.load_cgo(cluster_dict["13.5579996109"], "Features_13.5579996109", 1)
cmd.load_cgo(cluster_dict["13.5579996109_arrows"], "Arrows_13.5579996109")
cmd.set("transparency", 0.2,"Features_13.5579996109")
cmd.group("Pharmacophore_13.5579996109", members="Features_13.5579996109")
cmd.group("Pharmacophore_13.5579996109", members="Arrows_13.5579996109")

if dirpath:
    f = join(dirpath, "label_threshold_13.5579996109.mol2")
else:
    f = "label_threshold_13.5579996109.mol2"

cmd.load(f, 'label_threshold_13.5579996109')
cmd.hide('everything', 'label_threshold_13.5579996109')
cmd.label("label_threshold_13.5579996109", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.5579996109', members= 'label_threshold_13.5579996109')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
