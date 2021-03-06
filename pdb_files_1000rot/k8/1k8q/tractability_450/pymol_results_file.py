
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


cluster_dict = {"18.7619991302":[], "18.7619991302_arrows":[]}

cluster_dict["18.7619991302"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(33.5), float(80.5), float(22.5), float(1.0)]

cluster_dict["18.7619991302_arrows"] += cgo_arrow([33.5,80.5,22.5], [31.339,80.486,20.698], color="blue red", name="Arrows_18.7619991302_1")

cluster_dict["18.7619991302"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(38.0), float(87.0), float(25.5), float(1.0)]

cluster_dict["18.7619991302_arrows"] += cgo_arrow([38.0,87.0,25.5], [38.8,88.432,23.02], color="blue red", name="Arrows_18.7619991302_2")

cluster_dict["18.7619991302"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(33.3001115997), float(82.896159967), float(28.8653672922), float(1.0)]


cluster_dict["18.7619991302"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(79.5), float(25.0), float(1.0)]

cluster_dict["18.7619991302_arrows"] += cgo_arrow([26.0,79.5,25.0], [24.113,77.816,23.114], color="red blue", name="Arrows_18.7619991302_3")

cluster_dict["18.7619991302"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(32.5), float(79.0), float(22.5), float(1.0)]

cluster_dict["18.7619991302_arrows"] += cgo_arrow([32.5,79.0,22.5], [31.339,80.486,20.698], color="red blue", name="Arrows_18.7619991302_4")

cluster_dict["18.7619991302"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(34.0), float(86.0), float(27.0), float(1.0)]

cluster_dict["18.7619991302_arrows"] += cgo_arrow([34.0,86.0,27.0], [33.415,89.874,26.877], color="red blue", name="Arrows_18.7619991302_5")

cluster_dict["18.7619991302"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(36.0), float(82.5), float(33.0), float(1.0)]


cluster_dict["18.7619991302"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(38.0), float(86.0), float(25.5), float(1.0)]

cluster_dict["18.7619991302_arrows"] += cgo_arrow([38.0,86.0,25.5], [38.8,88.432,23.02], color="red blue", name="Arrows_18.7619991302_6")

cmd.load_cgo(cluster_dict["18.7619991302"], "Features_18.7619991302", 1)
cmd.load_cgo(cluster_dict["18.7619991302_arrows"], "Arrows_18.7619991302")
cmd.set("transparency", 0.2,"Features_18.7619991302")
cmd.group("Pharmacophore_18.7619991302", members="Features_18.7619991302")
cmd.group("Pharmacophore_18.7619991302", members="Arrows_18.7619991302")

if dirpath:
    f = join(dirpath, "label_threshold_18.7619991302.mol2")
else:
    f = "label_threshold_18.7619991302.mol2"

cmd.load(f, 'label_threshold_18.7619991302')
cmd.hide('everything', 'label_threshold_18.7619991302')
cmd.label("label_threshold_18.7619991302", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_18.7619991302', members= 'label_threshold_18.7619991302')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
