
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


cluster_dict = {"16.438999176":[], "16.438999176_arrows":[]}

cluster_dict["16.438999176"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(26.0), float(40.0), float(1.0)]

cluster_dict["16.438999176_arrows"] += cgo_arrow([23.5,26.0,40.0], [23.381,28.755,39.456], color="blue red", name="Arrows_16.438999176_1")

cluster_dict["16.438999176"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(25.5), float(25.5), float(46.0), float(1.0)]

cluster_dict["16.438999176_arrows"] += cgo_arrow([25.5,25.5,46.0], [24.992,22.671,45.425], color="blue red", name="Arrows_16.438999176_2")

cluster_dict["16.438999176"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.5), float(23.5), float(45.0), float(1.0)]

cluster_dict["16.438999176_arrows"] += cgo_arrow([27.5,23.5,45.0], [24.992,22.671,45.425], color="blue red", name="Arrows_16.438999176_3")

cluster_dict["16.438999176"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(25.5995540844), float(26.1197066298), float(39.5409655005), float(1.0)]


cluster_dict["16.438999176"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(26.0), float(32.0), float(1.0)]

cluster_dict["16.438999176_arrows"] += cgo_arrow([25.0,26.0,32.0], [26.885,27.06,28.588], color="red blue", name="Arrows_16.438999176_4")

cluster_dict["16.438999176"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(25.5), float(44.5), float(1.0)]

cluster_dict["16.438999176_arrows"] += cgo_arrow([25.0,25.5,44.5], [24.992,22.671,45.425], color="red blue", name="Arrows_16.438999176_5")

cluster_dict["16.438999176"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(25.5), float(24.5), float(38.0), float(1.0)]

cluster_dict["16.438999176_arrows"] += cgo_arrow([25.5,24.5,38.0], [28.183,26.981,38.218], color="red blue", name="Arrows_16.438999176_6")

cluster_dict["16.438999176"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(23.5), float(44.0), float(1.0)]

cluster_dict["16.438999176_arrows"] += cgo_arrow([26.0,23.5,44.0], [24.992,22.671,45.425], color="red blue", name="Arrows_16.438999176_7")

cluster_dict["16.438999176"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(26.5), float(41.0), float(1.0)]

cluster_dict["16.438999176_arrows"] += cgo_arrow([27.0,26.5,41.0], [28.183,26.981,38.218], color="red blue", name="Arrows_16.438999176_8")

cluster_dict["16.438999176"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.5), float(28.5), float(43.0), float(1.0)]

cluster_dict["16.438999176_arrows"] += cgo_arrow([27.5,28.5,43.0], [24.637,29.939,41.447], color="red blue", name="Arrows_16.438999176_9")

cmd.load_cgo(cluster_dict["16.438999176"], "Features_16.438999176", 1)
cmd.load_cgo(cluster_dict["16.438999176_arrows"], "Arrows_16.438999176")
cmd.set("transparency", 0.2,"Features_16.438999176")
cmd.group("Pharmacophore_16.438999176", members="Features_16.438999176")
cmd.group("Pharmacophore_16.438999176", members="Arrows_16.438999176")

if dirpath:
    f = join(dirpath, "label_threshold_16.438999176.mol2")
else:
    f = "label_threshold_16.438999176.mol2"

cmd.load(f, 'label_threshold_16.438999176')
cmd.hide('everything', 'label_threshold_16.438999176')
cmd.label("label_threshold_16.438999176", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.438999176', members= 'label_threshold_16.438999176')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
