
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


cluster_dict = {"13.1260004044":[], "13.1260004044_arrows":[]}

cluster_dict["13.1260004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-15.0), float(55.5), float(51.0), float(1.0)]

cluster_dict["13.1260004044_arrows"] += cgo_arrow([-15.0,55.5,51.0], [-16.296,53.603,52.731], color="blue red", name="Arrows_13.1260004044_1")

cluster_dict["13.1260004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-15.5), float(59.5), float(53.5), float(1.0)]

cluster_dict["13.1260004044_arrows"] += cgo_arrow([-15.5,59.5,53.5], [-13.483,59.635,55.763], color="blue red", name="Arrows_13.1260004044_2")

cluster_dict["13.1260004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-14.0), float(60.0), float(49.0), float(1.0)]

cluster_dict["13.1260004044_arrows"] += cgo_arrow([-14.0,60.0,49.0], [-12.485,60.798,46.613], color="blue red", name="Arrows_13.1260004044_3")

cluster_dict["13.1260004044"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-9.0), float(52.0), float(54.5), float(1.0)]

cluster_dict["13.1260004044_arrows"] += cgo_arrow([-9.0,52.0,54.5], [-9.093,53.902,56.945], color="blue red", name="Arrows_13.1260004044_4")

cluster_dict["13.1260004044"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-13.0327662492), float(56.8224640238), float(50.3254443824), float(1.0)]


cluster_dict["13.1260004044"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-8.43995867934), float(50.5637504855), float(45.2401652827), float(1.0)]


cluster_dict["13.1260004044"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-6.14227605471), float(50.5832710407), float(46.3120215082), float(1.0)]


cluster_dict["13.1260004044"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-5.5), float(54.0), float(50.75), float(1.0)]


cluster_dict["13.1260004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-15.0), float(56.0), float(53.0), float(1.0)]

cluster_dict["13.1260004044_arrows"] += cgo_arrow([-15.0,56.0,53.0], [-14.031,53.869,54.187], color="red blue", name="Arrows_13.1260004044_5")

cluster_dict["13.1260004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-16.0), float(55.5), float(50.0), float(1.0)]

cluster_dict["13.1260004044_arrows"] += cgo_arrow([-16.0,55.5,50.0], [-16.296,53.603,52.731], color="red blue", name="Arrows_13.1260004044_6")

cluster_dict["13.1260004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-13.5), float(59.0), float(48.5), float(1.0)]

cluster_dict["13.1260004044_arrows"] += cgo_arrow([-13.5,59.0,48.5], [-12.485,60.798,46.613], color="red blue", name="Arrows_13.1260004044_7")

cluster_dict["13.1260004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-12.0), float(52.5), float(53.5), float(1.0)]

cluster_dict["13.1260004044_arrows"] += cgo_arrow([-12.0,52.5,53.5], [-14.031,53.869,54.187], color="red blue", name="Arrows_13.1260004044_8")

cluster_dict["13.1260004044"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-7.0), float(50.5), float(45.0), float(1.0)]

cluster_dict["13.1260004044_arrows"] += cgo_arrow([-7.0,50.5,45.0], [-4.658,52.086,44.13], color="red blue", name="Arrows_13.1260004044_9")

cmd.load_cgo(cluster_dict["13.1260004044"], "Features_13.1260004044", 1)
cmd.load_cgo(cluster_dict["13.1260004044_arrows"], "Arrows_13.1260004044")
cmd.set("transparency", 0.2,"Features_13.1260004044")
cmd.group("Pharmacophore_13.1260004044", members="Features_13.1260004044")
cmd.group("Pharmacophore_13.1260004044", members="Arrows_13.1260004044")

if dirpath:
    f = join(dirpath, "label_threshold_13.1260004044.mol2")
else:
    f = "label_threshold_13.1260004044.mol2"

cmd.load(f, 'label_threshold_13.1260004044')
cmd.hide('everything', 'label_threshold_13.1260004044')
cmd.label("label_threshold_13.1260004044", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.1260004044', members= 'label_threshold_13.1260004044')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
