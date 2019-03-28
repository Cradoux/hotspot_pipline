
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
zip_dir = 'hotspot_boundaries.zip'
with zipfile.ZipFile(zip_dir) as hs_zip:
    hs_zip.extractall(dirpath)

cmd.load(join(dirpath,"protein.pdb"), "protein")
cmd.show("cartoon", "protein")

if dirpath:
    f = join(dirpath, "0/label_threshold_0.6.mol2")
else:
    f = "0/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
gfiles = ['0/donor.grd', '0/apolar.grd', '0/acceptor.grd']
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


cluster_dict = {"9.14599990845":[], "9.14599990845_arrows":[]}

cluster_dict["9.14599990845"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(9.0), float(36.5), float(15.0), float(1.0)]

cluster_dict["9.14599990845_arrows"] += cgo_arrow([9.0,36.5,15.0], [10.655,38.838,14.911], color="blue red", name="Arrows_9.14599990845_1")

cluster_dict["9.14599990845"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(12.0158887686), float(35.2643020793), float(16.2219633953), float(1.0)]


cluster_dict["9.14599990845"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(8.11989185836), float(36.3963460568), float(16.7524496191), float(1.0)]


cmd.load_cgo(cluster_dict["9.14599990845"], "Features_9.14599990845", 1)
cmd.load_cgo(cluster_dict["9.14599990845_arrows"], "Arrows_9.14599990845")
cmd.set("transparency", 0.2,"Features_9.14599990845")
cmd.group("Pharmacophore_9.14599990845", members="Features_9.14599990845")
cmd.group("Pharmacophore_9.14599990845", members="Arrows_9.14599990845")

if dirpath:
    f = join(dirpath, "0/label_threshold_9.14599990845.mol2")
else:
    f = "0/label_threshold_9.14599990845.mol2"

cmd.load(f, 'label_threshold_9.14599990845')
cmd.hide('everything', 'label_threshold_9.14599990845')
cmd.label("label_threshold_9.14599990845", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.14599990845', members= 'label_threshold_9.14599990845')


if dirpath:
    f = join(dirpath, '0/mesh.grd')
else:
    f = '0/mesh.grd'
cmd.load(f, 'mesh_0')
cmd.isomesh("isomesh_0", "mesh_0", 0.9)
cmd.color("grey80", "isomesh_0")
cmd.set('transparency', 0.4, "isomesh_0")

cmd.group('hotspot_0', "isomesh_0")
cmd.group('hotspot_0', "mesh_0")

if dirpath:
    f = join(dirpath, "1/label_threshold_3.2.mol2")
else:
    f = "1/label_threshold_3.2.mol2"

cmd.load(f, 'label_threshold_3.2')
cmd.hide('everything', 'label_threshold_3.2')
cmd.label("label_threshold_3.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [3.2]
gfiles = ['1/donor.grd', '1/apolar.grd', '1/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 1
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


cluster_dict = {"14.4499998093":[], "14.4499998093_arrows":[]}

cluster_dict["14.4499998093"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.0), float(57.5), float(17.0), float(1.0)]

cluster_dict["14.4499998093_arrows"] += cgo_arrow([30.0,57.5,17.0], [29.056,56.897,15.143], color="blue red", name="Arrows_14.4499998093_1")

cluster_dict["14.4499998093"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(28.2467382629), float(58.4668294205), float(20.3386778409), float(1.0)]


cluster_dict["14.4499998093"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(29.0), float(59.5), float(17.5), float(1.0)]

cluster_dict["14.4499998093_arrows"] += cgo_arrow([29.0,59.5,17.5], [29.461,59.893,15.268], color="red blue", name="Arrows_14.4499998093_2")

cmd.load_cgo(cluster_dict["14.4499998093"], "Features_14.4499998093", 1)
cmd.load_cgo(cluster_dict["14.4499998093_arrows"], "Arrows_14.4499998093")
cmd.set("transparency", 0.2,"Features_14.4499998093")
cmd.group("Pharmacophore_14.4499998093", members="Features_14.4499998093")
cmd.group("Pharmacophore_14.4499998093", members="Arrows_14.4499998093")

if dirpath:
    f = join(dirpath, "1/label_threshold_14.4499998093.mol2")
else:
    f = "1/label_threshold_14.4499998093.mol2"

cmd.load(f, 'label_threshold_14.4499998093')
cmd.hide('everything', 'label_threshold_14.4499998093')
cmd.label("label_threshold_14.4499998093", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.4499998093', members= 'label_threshold_14.4499998093')


if dirpath:
    f = join(dirpath, '1/mesh.grd')
else:
    f = '1/mesh.grd'
cmd.load(f, 'mesh_1')
cmd.isomesh("isomesh_1", "mesh_1", 0.9)
cmd.color("grey80", "isomesh_1")
cmd.set('transparency', 0.4, "isomesh_1")

cmd.group('hotspot_1', "isomesh_1")
cmd.group('hotspot_1', "mesh_1")

if dirpath:
    f = join(dirpath, "2/label_threshold_13.5.mol2")
else:
    f = "2/label_threshold_13.5.mol2"

cmd.load(f, 'label_threshold_13.5')
cmd.hide('everything', 'label_threshold_13.5')
cmd.label("label_threshold_13.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.5]
gfiles = ['2/donor.grd', '2/apolar.grd', '2/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 2
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


cluster_dict = {"14.1999998093":[], "14.1999998093_arrows":[]}

cluster_dict["14.1999998093"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(8.5), float(48.0), float(30.5), float(1.0)]

cluster_dict["14.1999998093_arrows"] += cgo_arrow([8.5,48.0,30.5], [6.903,49.386,28.3], color="blue red", name="Arrows_14.1999998093_1")

cluster_dict["14.1999998093"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(10.5466226769), float(50.9778157924), float(32.266763999), float(1.0)]


cluster_dict["14.1999998093"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(9.63553268743), float(42.3351705729), float(28.5430374556), float(1.0)]


cluster_dict["14.1999998093"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(11.8130540041), float(38.0209836191), float(18.2251644182), float(1.0)]


cluster_dict["14.1999998093"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(16.7834392587), float(49.5355418893), float(23.4477146961), float(1.0)]


cluster_dict["14.1999998093"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(25.1469282171), float(56.1167230082), float(21.7364278934), float(1.0)]


cluster_dict["14.1999998093"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(9.0), float(47.0), float(25.5), float(1.0)]

cluster_dict["14.1999998093_arrows"] += cgo_arrow([9.0,47.0,25.5], [7.61,46.423,23.631], color="red blue", name="Arrows_14.1999998093_2")

cluster_dict["14.1999998093"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(9.5), float(52.5), float(29.5), float(1.0)]

cluster_dict["14.1999998093_arrows"] += cgo_arrow([9.5,52.5,29.5], [7.675,53.976,27.475], color="red blue", name="Arrows_14.1999998093_3")

cluster_dict["14.1999998093"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(14.5), float(49.5), float(19.5), float(1.0)]

cluster_dict["14.1999998093_arrows"] += cgo_arrow([14.5,49.5,19.5], [13.87,48.102,17.509], color="red blue", name="Arrows_14.1999998093_4")

cluster_dict["14.1999998093"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(43.5), float(15.5), float(1.0)]

cluster_dict["14.1999998093_arrows"] += cgo_arrow([15.0,43.5,15.5], [12.884,44.157,16.767], color="red blue", name="Arrows_14.1999998093_5")

cluster_dict["14.1999998093"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(14.5), float(52.5), float(19.0), float(1.0)]

cluster_dict["14.1999998093_arrows"] += cgo_arrow([14.5,52.5,19.0], [14.391,54.985,20.064], color="red blue", name="Arrows_14.1999998093_6")

cluster_dict["14.1999998093"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.5), float(46.5), float(16.5), float(1.0)]

cluster_dict["14.1999998093_arrows"] += cgo_arrow([16.5,46.5,16.5], [13.87,48.102,17.509], color="red blue", name="Arrows_14.1999998093_7")

cluster_dict["14.1999998093"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(53.0), float(13.0), float(1.0)]

cluster_dict["14.1999998093_arrows"] += cgo_arrow([17.5,53.0,13.0], [17.796,54.621,11.551], color="red blue", name="Arrows_14.1999998093_8")

cluster_dict["14.1999998093"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(47.0), float(16.0), float(1.0)]

cluster_dict["14.1999998093_arrows"] += cgo_arrow([24.0,47.0,16.0], [23.609,46.125,13.689], color="red blue", name="Arrows_14.1999998093_9")

cmd.load_cgo(cluster_dict["14.1999998093"], "Features_14.1999998093", 1)
cmd.load_cgo(cluster_dict["14.1999998093_arrows"], "Arrows_14.1999998093")
cmd.set("transparency", 0.2,"Features_14.1999998093")
cmd.group("Pharmacophore_14.1999998093", members="Features_14.1999998093")
cmd.group("Pharmacophore_14.1999998093", members="Arrows_14.1999998093")

if dirpath:
    f = join(dirpath, "2/label_threshold_14.1999998093.mol2")
else:
    f = "2/label_threshold_14.1999998093.mol2"

cmd.load(f, 'label_threshold_14.1999998093')
cmd.hide('everything', 'label_threshold_14.1999998093')
cmd.label("label_threshold_14.1999998093", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.1999998093', members= 'label_threshold_14.1999998093')


if dirpath:
    f = join(dirpath, '2/mesh.grd')
else:
    f = '2/mesh.grd'
cmd.load(f, 'mesh_2')
cmd.isomesh("isomesh_2", "mesh_2", 0.9)
cmd.color("grey80", "isomesh_2")
cmd.set('transparency', 0.4, "isomesh_2")

cmd.group('hotspot_2', "isomesh_2")
cmd.group('hotspot_2', "mesh_2")

if dirpath:
    f = join(dirpath, "3/label_threshold_13.1.mol2")
else:
    f = "3/label_threshold_13.1.mol2"

cmd.load(f, 'label_threshold_13.1')
cmd.hide('everything', 'label_threshold_13.1')
cmd.label("label_threshold_13.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.1]
gfiles = ['3/donor.grd', '3/apolar.grd', '3/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 3
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


cluster_dict = {"14.0469999313":[], "14.0469999313_arrows":[]}

cluster_dict["14.0469999313"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(8.5), float(48.0), float(30.5), float(1.0)]

cluster_dict["14.0469999313_arrows"] += cgo_arrow([8.5,48.0,30.5], [6.903,49.386,28.3], color="blue red", name="Arrows_14.0469999313_1")

cluster_dict["14.0469999313"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(47.0), float(31.0), float(1.0)]

cluster_dict["14.0469999313_arrows"] += cgo_arrow([17.5,47.0,31.0], [19.253,49.535,30.637], color="blue red", name="Arrows_14.0469999313_2")

cluster_dict["14.0469999313"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(10.0179326225), float(50.5490932314), float(32.8043284009), float(1.0)]


cluster_dict["14.0469999313"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(9.54928674366), float(42.5005161364), float(28.4681332265), float(1.0)]


cluster_dict["14.0469999313"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(11.8882438181), float(37.324253857), float(18.7017274143), float(1.0)]


cluster_dict["14.0469999313"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(16.4966830712), float(48.045871464), float(25.3783710789), float(1.0)]


cluster_dict["14.0469999313"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(9.0), float(47.0), float(25.5), float(1.0)]

cluster_dict["14.0469999313_arrows"] += cgo_arrow([9.0,47.0,25.5], [7.61,46.423,23.631], color="red blue", name="Arrows_14.0469999313_3")

cluster_dict["14.0469999313"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(9.5), float(52.5), float(29.5), float(1.0)]

cluster_dict["14.0469999313_arrows"] += cgo_arrow([9.5,52.5,29.5], [7.675,53.976,27.475], color="red blue", name="Arrows_14.0469999313_4")

cluster_dict["14.0469999313"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(13.0), float(48.0), float(26.5), float(1.0)]

cluster_dict["14.0469999313_arrows"] += cgo_arrow([13.0,48.0,26.5], [13.665,46.492,23.997], color="red blue", name="Arrows_14.0469999313_5")

cluster_dict["14.0469999313"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(44.5), float(31.0), float(1.0)]

cluster_dict["14.0469999313_arrows"] += cgo_arrow([15.0,44.5,31.0], [13.536,44.497,33.984], color="red blue", name="Arrows_14.0469999313_6")

cluster_dict["14.0469999313"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(14.5), float(49.0), float(19.5), float(1.0)]

cluster_dict["14.0469999313_arrows"] += cgo_arrow([14.5,49.0,19.5], [13.87,48.102,17.509], color="red blue", name="Arrows_14.0469999313_7")

cluster_dict["14.0469999313"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(43.0), float(18.5), float(1.0)]

cluster_dict["14.0469999313_arrows"] += cgo_arrow([15.0,43.0,18.5], [12.884,44.157,16.767], color="red blue", name="Arrows_14.0469999313_8")

cluster_dict["14.0469999313"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.0), float(45.5), float(27.0), float(1.0)]

cluster_dict["14.0469999313_arrows"] += cgo_arrow([17.0,45.5,27.0], [19.435,44.038,28.278], color="red blue", name="Arrows_14.0469999313_9")

cmd.load_cgo(cluster_dict["14.0469999313"], "Features_14.0469999313", 1)
cmd.load_cgo(cluster_dict["14.0469999313_arrows"], "Arrows_14.0469999313")
cmd.set("transparency", 0.2,"Features_14.0469999313")
cmd.group("Pharmacophore_14.0469999313", members="Features_14.0469999313")
cmd.group("Pharmacophore_14.0469999313", members="Arrows_14.0469999313")

if dirpath:
    f = join(dirpath, "3/label_threshold_14.0469999313.mol2")
else:
    f = "3/label_threshold_14.0469999313.mol2"

cmd.load(f, 'label_threshold_14.0469999313')
cmd.hide('everything', 'label_threshold_14.0469999313')
cmd.label("label_threshold_14.0469999313", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.0469999313', members= 'label_threshold_14.0469999313')


if dirpath:
    f = join(dirpath, '3/mesh.grd')
else:
    f = '3/mesh.grd'
cmd.load(f, 'mesh_3')
cmd.isomesh("isomesh_3", "mesh_3", 0.9)
cmd.color("grey80", "isomesh_3")
cmd.set('transparency', 0.4, "isomesh_3")

cmd.group('hotspot_3', "isomesh_3")
cmd.group('hotspot_3', "mesh_3")

if dirpath:
    f = join(dirpath, "4/label_threshold_13.2.mol2")
else:
    f = "4/label_threshold_13.2.mol2"

cmd.load(f, 'label_threshold_13.2')
cmd.hide('everything', 'label_threshold_13.2')
cmd.label("label_threshold_13.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.2]
gfiles = ['4/donor.grd', '4/apolar.grd', '4/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 4
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


cluster_dict = {"14.029999733":[], "14.029999733_arrows":[]}

cluster_dict["14.029999733"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(8.5), float(48.0), float(31.0), float(1.0)]

cluster_dict["14.029999733_arrows"] += cgo_arrow([8.5,48.0,31.0], [10.691,46.223,32.917], color="blue red", name="Arrows_14.029999733_1")

cluster_dict["14.029999733"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(8.5), float(48.0), float(31.0), float(1.0)]

cluster_dict["14.029999733_arrows"] += cgo_arrow([8.5,48.0,31.0], [10.691,46.223,32.917], color="blue red", name="Arrows_14.029999733_2")

cluster_dict["14.029999733"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(47.0), float(31.0), float(1.0)]

cluster_dict["14.029999733_arrows"] += cgo_arrow([17.5,47.0,31.0], [19.253,49.535,30.637], color="blue red", name="Arrows_14.029999733_3")

cluster_dict["14.029999733"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(10.1082597901), float(51.0528100488), float(33.2701550528), float(1.0)]


cluster_dict["14.029999733"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(9.56101175825), float(42.5128776975), float(28.4222807974), float(1.0)]


cluster_dict["14.029999733"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(16.2538166464), float(48.0801179617), float(26.0585528729), float(1.0)]


cluster_dict["14.029999733"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(6.0), float(47.5), float(31.5), float(1.0)]

cluster_dict["14.029999733_arrows"] += cgo_arrow([6.0,47.5,31.5], [4.478,47.927,28.43], color="red blue", name="Arrows_14.029999733_4")

cluster_dict["14.029999733"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(9.0), float(47.0), float(25.5), float(1.0)]

cluster_dict["14.029999733_arrows"] += cgo_arrow([9.0,47.0,25.5], [7.61,46.423,23.631], color="red blue", name="Arrows_14.029999733_5")

cluster_dict["14.029999733"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(9.5), float(52.5), float(29.5), float(1.0)]

cluster_dict["14.029999733_arrows"] += cgo_arrow([9.5,52.5,29.5], [7.675,53.976,27.475], color="red blue", name="Arrows_14.029999733_6")

cluster_dict["14.029999733"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(13.0), float(48.0), float(26.5), float(1.0)]

cluster_dict["14.029999733_arrows"] += cgo_arrow([13.0,48.0,26.5], [13.665,46.492,23.997], color="red blue", name="Arrows_14.029999733_7")

cluster_dict["14.029999733"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(44.5), float(31.0), float(1.0)]

cluster_dict["14.029999733_arrows"] += cgo_arrow([15.0,44.5,31.0], [13.536,44.497,33.984], color="red blue", name="Arrows_14.029999733_8")

cluster_dict["14.029999733"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.0), float(45.5), float(27.0), float(1.0)]

cluster_dict["14.029999733_arrows"] += cgo_arrow([17.0,45.5,27.0], [19.435,44.038,28.278], color="red blue", name="Arrows_14.029999733_9")

cmd.load_cgo(cluster_dict["14.029999733"], "Features_14.029999733", 1)
cmd.load_cgo(cluster_dict["14.029999733_arrows"], "Arrows_14.029999733")
cmd.set("transparency", 0.2,"Features_14.029999733")
cmd.group("Pharmacophore_14.029999733", members="Features_14.029999733")
cmd.group("Pharmacophore_14.029999733", members="Arrows_14.029999733")

if dirpath:
    f = join(dirpath, "4/label_threshold_14.029999733.mol2")
else:
    f = "4/label_threshold_14.029999733.mol2"

cmd.load(f, 'label_threshold_14.029999733')
cmd.hide('everything', 'label_threshold_14.029999733')
cmd.label("label_threshold_14.029999733", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.029999733', members= 'label_threshold_14.029999733')


if dirpath:
    f = join(dirpath, '4/mesh.grd')
else:
    f = '4/mesh.grd'
cmd.load(f, 'mesh_4')
cmd.isomesh("isomesh_4", "mesh_4", 0.9)
cmd.color("grey80", "isomesh_4")
cmd.set('transparency', 0.4, "isomesh_4")

cmd.group('hotspot_4', "isomesh_4")
cmd.group('hotspot_4', "mesh_4")

if dirpath:
    f = join(dirpath, "5/label_threshold_13.4.mol2")
else:
    f = "5/label_threshold_13.4.mol2"

cmd.load(f, 'label_threshold_13.4')
cmd.hide('everything', 'label_threshold_13.4')
cmd.label("label_threshold_13.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.4]
gfiles = ['5/donor.grd', '5/apolar.grd', '5/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 5
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


cluster_dict = {"14.0069999695":[], "14.0069999695_arrows":[]}

cluster_dict["14.0069999695"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(8.5), float(48.0), float(30.5), float(1.0)]

cluster_dict["14.0069999695_arrows"] += cgo_arrow([8.5,48.0,30.5], [6.903,49.386,28.3], color="blue red", name="Arrows_14.0069999695_1")

cluster_dict["14.0069999695"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(47.0), float(31.0), float(1.0)]

cluster_dict["14.0069999695_arrows"] += cgo_arrow([17.5,47.0,31.0], [19.253,49.535,30.637], color="blue red", name="Arrows_14.0069999695_2")

cluster_dict["14.0069999695"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(10.0601115754), float(51.0281572296), float(33.2696142683), float(1.0)]


cluster_dict["14.0069999695"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(9.66812327841), float(42.3890340659), float(28.7282340608), float(1.0)]


cluster_dict["14.0069999695"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(16.7515374898), float(49.4054174547), float(23.5853131529), float(1.0)]


cluster_dict["14.0069999695"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(9.0), float(47.0), float(25.5), float(1.0)]

cluster_dict["14.0069999695_arrows"] += cgo_arrow([9.0,47.0,25.5], [7.61,46.423,23.631], color="red blue", name="Arrows_14.0069999695_3")

cluster_dict["14.0069999695"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(9.5), float(52.5), float(29.5), float(1.0)]

cluster_dict["14.0069999695_arrows"] += cgo_arrow([9.5,52.5,29.5], [7.675,53.976,27.475], color="red blue", name="Arrows_14.0069999695_4")

cluster_dict["14.0069999695"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(13.0), float(48.0), float(26.5), float(1.0)]

cluster_dict["14.0069999695_arrows"] += cgo_arrow([13.0,48.0,26.5], [13.665,46.492,23.997], color="red blue", name="Arrows_14.0069999695_5")

cluster_dict["14.0069999695"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(44.5), float(31.0), float(1.0)]

cluster_dict["14.0069999695_arrows"] += cgo_arrow([15.0,44.5,31.0], [13.536,44.497,33.984], color="red blue", name="Arrows_14.0069999695_6")

cluster_dict["14.0069999695"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(14.5), float(49.5), float(19.5), float(1.0)]

cluster_dict["14.0069999695_arrows"] += cgo_arrow([14.5,49.5,19.5], [13.87,48.102,17.509], color="red blue", name="Arrows_14.0069999695_7")

cluster_dict["14.0069999695"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(14.5), float(52.5), float(19.0), float(1.0)]

cluster_dict["14.0069999695_arrows"] += cgo_arrow([14.5,52.5,19.0], [14.391,54.985,20.064], color="red blue", name="Arrows_14.0069999695_8")

cluster_dict["14.0069999695"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.0), float(45.5), float(27.0), float(1.0)]

cluster_dict["14.0069999695_arrows"] += cgo_arrow([17.0,45.5,27.0], [19.435,44.038,28.278], color="red blue", name="Arrows_14.0069999695_9")

cmd.load_cgo(cluster_dict["14.0069999695"], "Features_14.0069999695", 1)
cmd.load_cgo(cluster_dict["14.0069999695_arrows"], "Arrows_14.0069999695")
cmd.set("transparency", 0.2,"Features_14.0069999695")
cmd.group("Pharmacophore_14.0069999695", members="Features_14.0069999695")
cmd.group("Pharmacophore_14.0069999695", members="Arrows_14.0069999695")

if dirpath:
    f = join(dirpath, "5/label_threshold_14.0069999695.mol2")
else:
    f = "5/label_threshold_14.0069999695.mol2"

cmd.load(f, 'label_threshold_14.0069999695')
cmd.hide('everything', 'label_threshold_14.0069999695')
cmd.label("label_threshold_14.0069999695", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.0069999695', members= 'label_threshold_14.0069999695')


if dirpath:
    f = join(dirpath, '5/mesh.grd')
else:
    f = '5/mesh.grd'
cmd.load(f, 'mesh_5')
cmd.isomesh("isomesh_5", "mesh_5", 0.9)
cmd.color("grey80", "isomesh_5")
cmd.set('transparency', 0.4, "isomesh_5")

cmd.group('hotspot_5', "isomesh_5")
cmd.group('hotspot_5', "mesh_5")

if dirpath:
    f = join(dirpath, "6/label_threshold_12.0.mol2")
else:
    f = "6/label_threshold_12.0.mol2"

cmd.load(f, 'label_threshold_12.0')
cmd.hide('everything', 'label_threshold_12.0')
cmd.label("label_threshold_12.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.0]
gfiles = ['6/donor.grd', '6/apolar.grd', '6/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 6
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


cluster_dict = {"13.8030004501":[], "13.8030004501_arrows":[]}

cluster_dict["13.8030004501"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(12.5), float(53.5), float(25.0), float(1.0)]

cluster_dict["13.8030004501_arrows"] += cgo_arrow([12.5,53.5,25.0], [9.808,52.921,25.117], color="blue red", name="Arrows_13.8030004501_1")

cluster_dict["13.8030004501"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(51.5), float(25.0), float(1.0)]

cluster_dict["13.8030004501_arrows"] += cgo_arrow([22.5,51.5,25.0], [22.456,50.136,27.967], color="blue red", name="Arrows_13.8030004501_2")

cluster_dict["13.8030004501"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(11.7081391514), float(52.4065538629), float(31.6261024367), float(1.0)]


cluster_dict["13.8030004501"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(16.868873535), float(49.8463037256), float(24.0535246589), float(1.0)]


cluster_dict["13.8030004501"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(25.4334689226), float(58.1746042479), float(20.8304589141), float(1.0)]


cluster_dict["13.8030004501"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(9.5), float(52.5), float(29.5), float(1.0)]

cluster_dict["13.8030004501_arrows"] += cgo_arrow([9.5,52.5,29.5], [7.675,53.976,27.475], color="red blue", name="Arrows_13.8030004501_3")

cluster_dict["13.8030004501"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(13.0), float(48.0), float(26.5), float(1.0)]

cluster_dict["13.8030004501_arrows"] += cgo_arrow([13.0,48.0,26.5], [13.665,46.492,23.997], color="red blue", name="Arrows_13.8030004501_4")

cluster_dict["13.8030004501"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(13.0), float(52.0), float(25.5), float(1.0)]

cluster_dict["13.8030004501_arrows"] += cgo_arrow([13.0,52.0,25.5], [10.627,51.048,26.078], color="red blue", name="Arrows_13.8030004501_5")

cluster_dict["13.8030004501"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(14.5), float(49.5), float(19.5), float(1.0)]

cluster_dict["13.8030004501_arrows"] += cgo_arrow([14.5,49.5,19.5], [13.87,48.102,17.509], color="red blue", name="Arrows_13.8030004501_6")

cluster_dict["13.8030004501"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(14.0), float(53.5), float(22.0), float(1.0)]

cluster_dict["13.8030004501_arrows"] += cgo_arrow([14.0,53.5,22.0], [14.391,54.985,20.064], color="red blue", name="Arrows_13.8030004501_7")

cluster_dict["13.8030004501"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.5), float(53.5), float(21.5), float(1.0)]

cluster_dict["13.8030004501_arrows"] += cgo_arrow([15.5,53.5,21.5], [14.391,54.985,20.064], color="red blue", name="Arrows_13.8030004501_8")

cluster_dict["13.8030004501"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(53.0), float(13.0), float(1.0)]

cluster_dict["13.8030004501_arrows"] += cgo_arrow([17.5,53.0,13.0], [17.796,54.621,11.551], color="red blue", name="Arrows_13.8030004501_9")

cmd.load_cgo(cluster_dict["13.8030004501"], "Features_13.8030004501", 1)
cmd.load_cgo(cluster_dict["13.8030004501_arrows"], "Arrows_13.8030004501")
cmd.set("transparency", 0.2,"Features_13.8030004501")
cmd.group("Pharmacophore_13.8030004501", members="Features_13.8030004501")
cmd.group("Pharmacophore_13.8030004501", members="Arrows_13.8030004501")

if dirpath:
    f = join(dirpath, "6/label_threshold_13.8030004501.mol2")
else:
    f = "6/label_threshold_13.8030004501.mol2"

cmd.load(f, 'label_threshold_13.8030004501')
cmd.hide('everything', 'label_threshold_13.8030004501')
cmd.label("label_threshold_13.8030004501", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.8030004501', members= 'label_threshold_13.8030004501')


if dirpath:
    f = join(dirpath, '6/mesh.grd')
else:
    f = '6/mesh.grd'
cmd.load(f, 'mesh_6')
cmd.isomesh("isomesh_6", "mesh_6", 0.9)
cmd.color("grey80", "isomesh_6")
cmd.set('transparency', 0.4, "isomesh_6")

cmd.group('hotspot_6', "isomesh_6")
cmd.group('hotspot_6', "mesh_6")

if dirpath:
    f = join(dirpath, "7/label_threshold_10.4.mol2")
else:
    f = "7/label_threshold_10.4.mol2"

cmd.load(f, 'label_threshold_10.4')
cmd.hide('everything', 'label_threshold_10.4')
cmd.label("label_threshold_10.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.4]
gfiles = ['7/donor.grd', '7/apolar.grd', '7/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 7
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


cluster_dict = {"13.6730003357":[], "13.6730003357_arrows":[]}

cluster_dict["13.6730003357"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(9.5), float(71.0), float(19.0), float(1.0)]

cluster_dict["13.6730003357_arrows"] += cgo_arrow([9.5,71.0,19.0], [7.119,71.943,19.278], color="blue red", name="Arrows_13.6730003357_1")

cluster_dict["13.6730003357"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(9.5), float(71.0), float(19.0), float(1.0)]

cluster_dict["13.6730003357_arrows"] += cgo_arrow([9.5,71.0,19.0], [7.119,71.943,19.278], color="blue red", name="Arrows_13.6730003357_2")

cluster_dict["13.6730003357"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(11.5), float(68.0), float(22.0), float(1.0)]

cluster_dict["13.6730003357_arrows"] += cgo_arrow([11.5,68.0,22.0], [12.568,67.69,25.092], color="blue red", name="Arrows_13.6730003357_3")

cluster_dict["13.6730003357"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(12.0), float(69.5), float(16.5), float(1.0)]

cluster_dict["13.6730003357_arrows"] += cgo_arrow([12.0,69.5,16.5], [10.812,69.187,14.423], color="blue red", name="Arrows_13.6730003357_4")

cluster_dict["13.6730003357"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(11.7131442811), float(71.6552138353), float(20.877953119), float(1.0)]


cluster_dict["13.6730003357"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(9.0), float(69.5), float(19.5), float(1.0)]

cluster_dict["13.6730003357_arrows"] += cgo_arrow([9.0,69.5,19.5], [4.371,70.166,19.381], color="red blue", name="Arrows_13.6730003357_5")

cluster_dict["13.6730003357"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(12.5), float(73.0), float(19.0), float(1.0)]

cluster_dict["13.6730003357_arrows"] += cgo_arrow([12.5,73.0,19.0], [12.12,74.968,20.641], color="red blue", name="Arrows_13.6730003357_6")

cluster_dict["13.6730003357"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(13.0), float(69.5), float(19.0), float(1.0)]

cluster_dict["13.6730003357_arrows"] += cgo_arrow([13.0,69.5,19.0], [12.834,66.793,18.647], color="red blue", name="Arrows_13.6730003357_7")

cmd.load_cgo(cluster_dict["13.6730003357"], "Features_13.6730003357", 1)
cmd.load_cgo(cluster_dict["13.6730003357_arrows"], "Arrows_13.6730003357")
cmd.set("transparency", 0.2,"Features_13.6730003357")
cmd.group("Pharmacophore_13.6730003357", members="Features_13.6730003357")
cmd.group("Pharmacophore_13.6730003357", members="Arrows_13.6730003357")

if dirpath:
    f = join(dirpath, "7/label_threshold_13.6730003357.mol2")
else:
    f = "7/label_threshold_13.6730003357.mol2"

cmd.load(f, 'label_threshold_13.6730003357')
cmd.hide('everything', 'label_threshold_13.6730003357')
cmd.label("label_threshold_13.6730003357", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.6730003357', members= 'label_threshold_13.6730003357')


if dirpath:
    f = join(dirpath, '7/mesh.grd')
else:
    f = '7/mesh.grd'
cmd.load(f, 'mesh_7')
cmd.isomesh("isomesh_7", "mesh_7", 0.9)
cmd.color("grey80", "isomesh_7")
cmd.set('transparency', 0.4, "isomesh_7")

cmd.group('hotspot_7', "isomesh_7")
cmd.group('hotspot_7', "mesh_7")

if dirpath:
    f = join(dirpath, "8/label_threshold_11.1.mol2")
else:
    f = "8/label_threshold_11.1.mol2"

cmd.load(f, 'label_threshold_11.1')
cmd.hide('everything', 'label_threshold_11.1')
cmd.label("label_threshold_11.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [11.1]
gfiles = ['8/donor.grd', '8/apolar.grd', '8/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 8
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


cluster_dict = {"13.2919998169":[], "13.2919998169_arrows":[]}

cluster_dict["13.2919998169"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(47.0), float(30.5), float(1.0)]

cluster_dict["13.2919998169_arrows"] += cgo_arrow([17.5,47.0,30.5], [19.253,49.535,30.637], color="blue red", name="Arrows_13.2919998169_1")

cluster_dict["13.2919998169"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(51.5), float(25.0), float(1.0)]

cluster_dict["13.2919998169_arrows"] += cgo_arrow([22.5,51.5,25.0], [22.456,50.136,27.967], color="blue red", name="Arrows_13.2919998169_2")

cluster_dict["13.2919998169"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(17.0364369555), float(49.0601345307), float(24.1689307181), float(1.0)]


cluster_dict["13.2919998169"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(26.7587366997), float(57.2602702291), float(20.8876712669), float(1.0)]


cluster_dict["13.2919998169"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(49.5), float(19.5), float(1.0)]

cluster_dict["13.2919998169_arrows"] += cgo_arrow([15.0,49.5,19.5], [13.87,48.102,17.509], color="red blue", name="Arrows_13.2919998169_3")

cluster_dict["13.2919998169"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.5), float(53.5), float(21.5), float(1.0)]

cluster_dict["13.2919998169_arrows"] += cgo_arrow([15.5,53.5,21.5], [14.391,54.985,20.064], color="red blue", name="Arrows_13.2919998169_4")

cluster_dict["13.2919998169"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.5), float(46.5), float(16.5), float(1.0)]

cluster_dict["13.2919998169_arrows"] += cgo_arrow([16.5,46.5,16.5], [13.87,48.102,17.509], color="red blue", name="Arrows_13.2919998169_5")

cluster_dict["13.2919998169"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.0), float(45.5), float(27.0), float(1.0)]

cluster_dict["13.2919998169_arrows"] += cgo_arrow([17.0,45.5,27.0], [19.435,44.038,28.278], color="red blue", name="Arrows_13.2919998169_6")

cluster_dict["13.2919998169"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(47.0), float(16.0), float(1.0)]

cluster_dict["13.2919998169_arrows"] += cgo_arrow([24.0,47.0,16.0], [23.609,46.125,13.689], color="red blue", name="Arrows_13.2919998169_7")

cmd.load_cgo(cluster_dict["13.2919998169"], "Features_13.2919998169", 1)
cmd.load_cgo(cluster_dict["13.2919998169_arrows"], "Arrows_13.2919998169")
cmd.set("transparency", 0.2,"Features_13.2919998169")
cmd.group("Pharmacophore_13.2919998169", members="Features_13.2919998169")
cmd.group("Pharmacophore_13.2919998169", members="Arrows_13.2919998169")

if dirpath:
    f = join(dirpath, "8/label_threshold_13.2919998169.mol2")
else:
    f = "8/label_threshold_13.2919998169.mol2"

cmd.load(f, 'label_threshold_13.2919998169')
cmd.hide('everything', 'label_threshold_13.2919998169')
cmd.label("label_threshold_13.2919998169", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.2919998169', members= 'label_threshold_13.2919998169')


if dirpath:
    f = join(dirpath, '8/mesh.grd')
else:
    f = '8/mesh.grd'
cmd.load(f, 'mesh_8')
cmd.isomesh("isomesh_8", "mesh_8", 0.9)
cmd.color("grey80", "isomesh_8")
cmd.set('transparency', 0.4, "isomesh_8")

cmd.group('hotspot_8', "isomesh_8")
cmd.group('hotspot_8', "mesh_8")

if dirpath:
    f = join(dirpath, "9/label_threshold_10.3.mol2")
else:
    f = "9/label_threshold_10.3.mol2"

cmd.load(f, 'label_threshold_10.3')
cmd.hide('everything', 'label_threshold_10.3')
cmd.label("label_threshold_10.3", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.3]
gfiles = ['9/donor.grd', '9/apolar.grd', '9/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 9
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


cluster_dict = {"12.9169998169":[], "12.9169998169_arrows":[]}

cluster_dict["12.9169998169"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(12.5), float(53.5), float(25.0), float(1.0)]

cluster_dict["12.9169998169_arrows"] += cgo_arrow([12.5,53.5,25.0], [9.808,52.921,25.117], color="blue red", name="Arrows_12.9169998169_1")

cluster_dict["12.9169998169"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(51.5), float(24.5), float(1.0)]

cluster_dict["12.9169998169_arrows"] += cgo_arrow([22.0,51.5,24.5], [23.453,52.763,22.084], color="blue red", name="Arrows_12.9169998169_2")

cluster_dict["12.9169998169"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(17.280066659), float(49.7251808966), float(23.0568923882), float(1.0)]


cluster_dict["12.9169998169"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(16.3921568627), float(52.4934640523), float(14.9869281046), float(1.0)]


cluster_dict["12.9169998169"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.8552973385), float(57.2985110025), float(20.639007015), float(1.0)]


cluster_dict["12.9169998169"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.5), float(53.5), float(21.5), float(1.0)]

cluster_dict["12.9169998169_arrows"] += cgo_arrow([15.5,53.5,21.5], [14.391,54.985,20.064], color="red blue", name="Arrows_12.9169998169_3")

cluster_dict["12.9169998169"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(14.0), float(53.5), float(22.0), float(1.0)]

cluster_dict["12.9169998169_arrows"] += cgo_arrow([14.0,53.5,22.0], [14.391,54.985,20.064], color="red blue", name="Arrows_12.9169998169_4")

cluster_dict["12.9169998169"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(43.5), float(15.5), float(1.0)]

cluster_dict["12.9169998169_arrows"] += cgo_arrow([15.0,43.5,15.5], [12.884,44.157,16.767], color="red blue", name="Arrows_12.9169998169_5")

cluster_dict["12.9169998169"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.5), float(46.5), float(16.5), float(1.0)]

cluster_dict["12.9169998169_arrows"] += cgo_arrow([16.5,46.5,16.5], [13.87,48.102,17.509], color="red blue", name="Arrows_12.9169998169_6")

cluster_dict["12.9169998169"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.5), float(53.5), float(21.5), float(1.0)]

cluster_dict["12.9169998169_arrows"] += cgo_arrow([15.5,53.5,21.5], [14.391,54.985,20.064], color="red blue", name="Arrows_12.9169998169_7")

cluster_dict["12.9169998169"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(53.0), float(13.0), float(1.0)]

cluster_dict["12.9169998169_arrows"] += cgo_arrow([17.5,53.0,13.0], [17.796,54.621,11.551], color="red blue", name="Arrows_12.9169998169_8")

cluster_dict["12.9169998169"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(47.0), float(16.0), float(1.0)]

cluster_dict["12.9169998169_arrows"] += cgo_arrow([24.0,47.0,16.0], [23.609,46.125,13.689], color="red blue", name="Arrows_12.9169998169_9")

cmd.load_cgo(cluster_dict["12.9169998169"], "Features_12.9169998169", 1)
cmd.load_cgo(cluster_dict["12.9169998169_arrows"], "Arrows_12.9169998169")
cmd.set("transparency", 0.2,"Features_12.9169998169")
cmd.group("Pharmacophore_12.9169998169", members="Features_12.9169998169")
cmd.group("Pharmacophore_12.9169998169", members="Arrows_12.9169998169")

if dirpath:
    f = join(dirpath, "9/label_threshold_12.9169998169.mol2")
else:
    f = "9/label_threshold_12.9169998169.mol2"

cmd.load(f, 'label_threshold_12.9169998169')
cmd.hide('everything', 'label_threshold_12.9169998169')
cmd.label("label_threshold_12.9169998169", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.9169998169', members= 'label_threshold_12.9169998169')


if dirpath:
    f = join(dirpath, '9/mesh.grd')
else:
    f = '9/mesh.grd'
cmd.load(f, 'mesh_9')
cmd.isomesh("isomesh_9", "mesh_9", 0.9)
cmd.color("grey80", "isomesh_9")
cmd.set('transparency', 0.4, "isomesh_9")

cmd.group('hotspot_9', "isomesh_9")
cmd.group('hotspot_9', "mesh_9")

if dirpath:
    f = join(dirpath, "10/label_threshold_1.8.mol2")
else:
    f = "10/label_threshold_1.8.mol2"

cmd.load(f, 'label_threshold_1.8')
cmd.hide('everything', 'label_threshold_1.8')
cmd.label("label_threshold_1.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [1.8]
gfiles = ['10/donor.grd', '10/apolar.grd', '10/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 10
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


cluster_dict = {"11.3680000305":[], "11.3680000305_arrows":[]}

cluster_dict["11.3680000305"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(39.5), float(-0.5), float(1.0)]

cluster_dict["11.3680000305_arrows"] += cgo_arrow([28.0,39.5,-0.5], [25.53,38.5,0.909], color="blue red", name="Arrows_11.3680000305_1")

cluster_dict["11.3680000305"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(25.8686260804), float(45.6866874305), float(-2.44401226525), float(1.0)]


cluster_dict["11.3680000305"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(41.0), float(0.0), float(1.0)]

cluster_dict["11.3680000305_arrows"] += cgo_arrow([23.5,41.0,0.0], [22.431,42.534,2.395], color="red blue", name="Arrows_11.3680000305_2")

cmd.load_cgo(cluster_dict["11.3680000305"], "Features_11.3680000305", 1)
cmd.load_cgo(cluster_dict["11.3680000305_arrows"], "Arrows_11.3680000305")
cmd.set("transparency", 0.2,"Features_11.3680000305")
cmd.group("Pharmacophore_11.3680000305", members="Features_11.3680000305")
cmd.group("Pharmacophore_11.3680000305", members="Arrows_11.3680000305")

if dirpath:
    f = join(dirpath, "10/label_threshold_11.3680000305.mol2")
else:
    f = "10/label_threshold_11.3680000305.mol2"

cmd.load(f, 'label_threshold_11.3680000305')
cmd.hide('everything', 'label_threshold_11.3680000305')
cmd.label("label_threshold_11.3680000305", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.3680000305', members= 'label_threshold_11.3680000305')


if dirpath:
    f = join(dirpath, '10/mesh.grd')
else:
    f = '10/mesh.grd'
cmd.load(f, 'mesh_10')
cmd.isomesh("isomesh_10", "mesh_10", 0.9)
cmd.color("grey80", "isomesh_10")
cmd.set('transparency', 0.4, "isomesh_10")

cmd.group('hotspot_10', "isomesh_10")
cmd.group('hotspot_10', "mesh_10")

if dirpath:
    f = join(dirpath, "11/label_threshold_0.6.mol2")
else:
    f = "11/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
gfiles = ['11/donor.grd', '11/apolar.grd', '11/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 11
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


cluster_dict = {"9.14599990845":[], "9.14599990845_arrows":[]}

cluster_dict["9.14599990845"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(9.0), float(36.5), float(15.0), float(1.0)]

cluster_dict["9.14599990845_arrows"] += cgo_arrow([9.0,36.5,15.0], [10.655,38.838,14.911], color="blue red", name="Arrows_9.14599990845_1")

cluster_dict["9.14599990845"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(12.0158887686), float(35.2643020793), float(16.2219633953), float(1.0)]


cluster_dict["9.14599990845"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(8.11989185836), float(36.3963460568), float(16.7524496191), float(1.0)]


cmd.load_cgo(cluster_dict["9.14599990845"], "Features_9.14599990845", 1)
cmd.load_cgo(cluster_dict["9.14599990845_arrows"], "Arrows_9.14599990845")
cmd.set("transparency", 0.2,"Features_9.14599990845")
cmd.group("Pharmacophore_9.14599990845", members="Features_9.14599990845")
cmd.group("Pharmacophore_9.14599990845", members="Arrows_9.14599990845")

if dirpath:
    f = join(dirpath, "11/label_threshold_9.14599990845.mol2")
else:
    f = "11/label_threshold_9.14599990845.mol2"

cmd.load(f, 'label_threshold_9.14599990845')
cmd.hide('everything', 'label_threshold_9.14599990845')
cmd.label("label_threshold_9.14599990845", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.14599990845', members= 'label_threshold_9.14599990845')


if dirpath:
    f = join(dirpath, '11/mesh.grd')
else:
    f = '11/mesh.grd'
cmd.load(f, 'mesh_11')
cmd.isomesh("isomesh_11", "mesh_11", 0.9)
cmd.color("grey80", "isomesh_11")
cmd.set('transparency', 0.4, "isomesh_11")

cmd.group('hotspot_11', "isomesh_11")
cmd.group('hotspot_11', "mesh_11")

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
