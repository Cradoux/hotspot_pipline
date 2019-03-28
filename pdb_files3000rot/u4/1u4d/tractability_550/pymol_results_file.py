
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


cluster_dict = {"17.9084997177":[], "17.9084997177_arrows":[]}

cluster_dict["17.9084997177"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.5), float(15.5), float(56.5), float(1.0)]

cluster_dict["17.9084997177_arrows"] += cgo_arrow([19.5,15.5,56.5], [17.441,14.76,54.737], color="blue red", name="Arrows_17.9084997177_1")

cluster_dict["17.9084997177"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.5), float(13.5), float(54.0), float(1.0)]

cluster_dict["17.9084997177_arrows"] += cgo_arrow([19.5,13.5,54.0], [17.441,14.76,54.737], color="blue red", name="Arrows_17.9084997177_2")

cluster_dict["17.9084997177"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.5), float(12.0), float(57.5), float(1.0)]

cluster_dict["17.9084997177_arrows"] += cgo_arrow([27.5,12.0,57.5], [27.369,10.108,59.891], color="blue red", name="Arrows_17.9084997177_3")

cluster_dict["17.9084997177"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(28.5), float(18.5), float(56.5), float(1.0)]

cluster_dict["17.9084997177_arrows"] += cgo_arrow([28.5,18.5,56.5], [27.775,18.719,53.781], color="blue red", name="Arrows_17.9084997177_4")

cluster_dict["17.9084997177"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.6065631823), float(15.6713736453), float(56.3563111096), float(1.0)]


cluster_dict["17.9084997177"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(19.5), float(13.5), float(56.0), float(1.0)]

cluster_dict["17.9084997177_arrows"] += cgo_arrow([19.5,13.5,56.0], [17.441,14.76,54.737], color="red blue", name="Arrows_17.9084997177_5")

cluster_dict["17.9084997177"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(18.0), float(56.5), float(1.0)]

cluster_dict["17.9084997177_arrows"] += cgo_arrow([21.0,18.0,56.5], [19.197,20.094,58.5], color="red blue", name="Arrows_17.9084997177_6")

cluster_dict["17.9084997177"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(13.0), float(59.0), float(1.0)]

cluster_dict["17.9084997177_arrows"] += cgo_arrow([24.5,13.0,59.0], [21.913,9.676,59.261], color="red blue", name="Arrows_17.9084997177_7")

cluster_dict["17.9084997177"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(28.5), float(14.0), float(55.5), float(1.0)]

cluster_dict["17.9084997177_arrows"] += cgo_arrow([28.5,14.0,55.5], [28.508,13.873,52.827], color="red blue", name="Arrows_17.9084997177_8")

cmd.load_cgo(cluster_dict["17.9084997177"], "Features_17.9084997177", 1)
cmd.load_cgo(cluster_dict["17.9084997177_arrows"], "Arrows_17.9084997177")
cmd.set("transparency", 0.2,"Features_17.9084997177")
cmd.group("Pharmacophore_17.9084997177", members="Features_17.9084997177")
cmd.group("Pharmacophore_17.9084997177", members="Arrows_17.9084997177")

if dirpath:
    f = join(dirpath, "label_threshold_17.9084997177.mol2")
else:
    f = "label_threshold_17.9084997177.mol2"

cmd.load(f, 'label_threshold_17.9084997177')
cmd.hide('everything', 'label_threshold_17.9084997177')
cmd.label("label_threshold_17.9084997177", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.9084997177', members= 'label_threshold_17.9084997177')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")