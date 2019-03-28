
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
    f = join(dirpath, "0/label_threshold_12.9.mol2")
else:
    f = "0/label_threshold_12.9.mol2"

cmd.load(f, 'label_threshold_12.9')
cmd.hide('everything', 'label_threshold_12.9')
cmd.label("label_threshold_12.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.9]
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


cluster_dict = {"15.8260002136":[], "15.8260002136_arrows":[]}

cluster_dict["15.8260002136"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(47.5), float(40.0), float(-10.5), float(1.0)]

cluster_dict["15.8260002136_arrows"] += cgo_arrow([47.5,40.0,-10.5], [49.122,39.952,-12.916], color="blue red", name="Arrows_15.8260002136_1")

cluster_dict["15.8260002136"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(47.5), float(39.0), float(-4.5), float(1.0)]

cluster_dict["15.8260002136_arrows"] += cgo_arrow([47.5,39.0,-4.5], [49.526,37.349,-2.499], color="blue red", name="Arrows_15.8260002136_2")

cluster_dict["15.8260002136"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(49.0), float(48.5), float(-4.5), float(1.0)]

cluster_dict["15.8260002136_arrows"] += cgo_arrow([49.0,48.5,-4.5], [50.269,47.32,-1.869], color="blue red", name="Arrows_15.8260002136_3")

cluster_dict["15.8260002136"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(50.0), float(45.5), float(-9.0), float(1.0)]

cluster_dict["15.8260002136_arrows"] += cgo_arrow([50.0,45.5,-9.0], [51.031,46.692,-11.391], color="blue red", name="Arrows_15.8260002136_4")

cluster_dict["15.8260002136"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(46.6231834463), float(42.3667636867), float(-6.61738139707), float(1.0)]


cluster_dict["15.8260002136"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(48.0), float(47.0), float(-5.0), float(1.0)]


cluster_dict["15.8260002136"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(48.0), float(42.5), float(-11.5), float(1.0)]

cluster_dict["15.8260002136_arrows"] += cgo_arrow([48.0,42.5,-11.5], [50.337,42.875,-12.964], color="red blue", name="Arrows_15.8260002136_5")

cluster_dict["15.8260002136"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(48.0), float(47.0), float(-5.0), float(1.0)]


cmd.load_cgo(cluster_dict["15.8260002136"], "Features_15.8260002136", 1)
cmd.load_cgo(cluster_dict["15.8260002136_arrows"], "Arrows_15.8260002136")
cmd.set("transparency", 0.2,"Features_15.8260002136")
cmd.group("Pharmacophore_15.8260002136", members="Features_15.8260002136")
cmd.group("Pharmacophore_15.8260002136", members="Arrows_15.8260002136")

if dirpath:
    f = join(dirpath, "0/label_threshold_15.8260002136.mol2")
else:
    f = "0/label_threshold_15.8260002136.mol2"

cmd.load(f, 'label_threshold_15.8260002136')
cmd.hide('everything', 'label_threshold_15.8260002136')
cmd.label("label_threshold_15.8260002136", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.8260002136', members= 'label_threshold_15.8260002136')


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
    f = join(dirpath, "1/label_threshold_8.8.mol2")
else:
    f = "1/label_threshold_8.8.mol2"

cmd.load(f, 'label_threshold_8.8')
cmd.hide('everything', 'label_threshold_8.8')
cmd.label("label_threshold_8.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [8.8]
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


cluster_dict = {"14.3319997787":[], "14.3319997787_arrows":[]}

cluster_dict["14.3319997787"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(69.0), float(48.5), float(6.0), float(1.0)]

cluster_dict["14.3319997787_arrows"] += cgo_arrow([69.0,48.5,6.0], [67.674,48.455,3.678], color="blue red", name="Arrows_14.3319997787_1")

cluster_dict["14.3319997787"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(70.0), float(46.0), float(0.5), float(1.0)]

cluster_dict["14.3319997787_arrows"] += cgo_arrow([70.0,46.0,0.5], [67.737,45.598,-0.45], color="blue red", name="Arrows_14.3319997787_2")

cluster_dict["14.3319997787"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(71.0), float(44.0), float(-0.5), float(1.0)]

cluster_dict["14.3319997787_arrows"] += cgo_arrow([71.0,44.0,-0.5], [72.426,42.145,-2.486], color="blue red", name="Arrows_14.3319997787_3")

cluster_dict["14.3319997787"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(72.3710610664), float(45.1717009121), float(3.83195332625), float(1.0)]


cluster_dict["14.3319997787"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(70.5), float(43.0), float(3.0), float(1.0)]

cluster_dict["14.3319997787_arrows"] += cgo_arrow([70.5,43.0,3.0], [70.844,40.247,1.133], color="red blue", name="Arrows_14.3319997787_4")

cluster_dict["14.3319997787"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(68.5), float(46.0), float(6.0), float(1.0)]


cluster_dict["14.3319997787"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(74.5), float(42.0), float(3.0), float(1.0)]

cluster_dict["14.3319997787_arrows"] += cgo_arrow([74.5,42.0,3.0], [72.687,40.241,4.233], color="red blue", name="Arrows_14.3319997787_5")

cmd.load_cgo(cluster_dict["14.3319997787"], "Features_14.3319997787", 1)
cmd.load_cgo(cluster_dict["14.3319997787_arrows"], "Arrows_14.3319997787")
cmd.set("transparency", 0.2,"Features_14.3319997787")
cmd.group("Pharmacophore_14.3319997787", members="Features_14.3319997787")
cmd.group("Pharmacophore_14.3319997787", members="Arrows_14.3319997787")

if dirpath:
    f = join(dirpath, "1/label_threshold_14.3319997787.mol2")
else:
    f = "1/label_threshold_14.3319997787.mol2"

cmd.load(f, 'label_threshold_14.3319997787')
cmd.hide('everything', 'label_threshold_14.3319997787')
cmd.label("label_threshold_14.3319997787", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.3319997787', members= 'label_threshold_14.3319997787')


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
    f = join(dirpath, "2/label_threshold_10.6.mol2")
else:
    f = "2/label_threshold_10.6.mol2"

cmd.load(f, 'label_threshold_10.6')
cmd.hide('everything', 'label_threshold_10.6')
cmd.label("label_threshold_10.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.6]
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


cluster_dict = {"14.1309995651":[], "14.1309995651_arrows":[]}

cluster_dict["14.1309995651"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(47.5), float(40.0), float(-10.5), float(1.0)]

cluster_dict["14.1309995651_arrows"] += cgo_arrow([47.5,40.0,-10.5], [49.122,39.952,-12.916], color="blue red", name="Arrows_14.1309995651_1")

cluster_dict["14.1309995651"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(47.5), float(39.0), float(-4.5), float(1.0)]

cluster_dict["14.1309995651_arrows"] += cgo_arrow([47.5,39.0,-4.5], [49.526,37.349,-2.499], color="blue red", name="Arrows_14.1309995651_2")

cluster_dict["14.1309995651"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(45.8135420546), float(40.166145912), float(-7.37082177534), float(1.0)]


cluster_dict["14.1309995651"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(46.5863896033), float(32.423670049), float(-11.8116214015), float(1.0)]


cluster_dict["14.1309995651"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(48.0), float(42.5), float(-11.5), float(1.0)]

cluster_dict["14.1309995651_arrows"] += cgo_arrow([48.0,42.5,-11.5], [50.337,42.875,-12.964], color="red blue", name="Arrows_14.1309995651_3")

cmd.load_cgo(cluster_dict["14.1309995651"], "Features_14.1309995651", 1)
cmd.load_cgo(cluster_dict["14.1309995651_arrows"], "Arrows_14.1309995651")
cmd.set("transparency", 0.2,"Features_14.1309995651")
cmd.group("Pharmacophore_14.1309995651", members="Features_14.1309995651")
cmd.group("Pharmacophore_14.1309995651", members="Arrows_14.1309995651")

if dirpath:
    f = join(dirpath, "2/label_threshold_14.1309995651.mol2")
else:
    f = "2/label_threshold_14.1309995651.mol2"

cmd.load(f, 'label_threshold_14.1309995651')
cmd.hide('everything', 'label_threshold_14.1309995651')
cmd.label("label_threshold_14.1309995651", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.1309995651', members= 'label_threshold_14.1309995651')


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
    f = join(dirpath, "3/label_threshold_0.5.mol2")
else:
    f = "3/label_threshold_0.5.mol2"

cmd.load(f, 'label_threshold_0.5')
cmd.hide('everything', 'label_threshold_0.5')
cmd.label("label_threshold_0.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.5]
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


cluster_dict = {"11.5069999695":[], "11.5069999695_arrows":[]}

cluster_dict["11.5069999695"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(57.0), float(63.0), float(1.5), float(1.0)]

cluster_dict["11.5069999695_arrows"] += cgo_arrow([57.0,63.0,1.5], [56.65,61.169,3.099], color="blue red", name="Arrows_11.5069999695_1")

cluster_dict["11.5069999695"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(59.6087650141), float(63.2259235886), float(-0.3881769583), float(1.0)]


cmd.load_cgo(cluster_dict["11.5069999695"], "Features_11.5069999695", 1)
cmd.load_cgo(cluster_dict["11.5069999695_arrows"], "Arrows_11.5069999695")
cmd.set("transparency", 0.2,"Features_11.5069999695")
cmd.group("Pharmacophore_11.5069999695", members="Features_11.5069999695")
cmd.group("Pharmacophore_11.5069999695", members="Arrows_11.5069999695")

if dirpath:
    f = join(dirpath, "3/label_threshold_11.5069999695.mol2")
else:
    f = "3/label_threshold_11.5069999695.mol2"

cmd.load(f, 'label_threshold_11.5069999695')
cmd.hide('everything', 'label_threshold_11.5069999695')
cmd.label("label_threshold_11.5069999695", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.5069999695', members= 'label_threshold_11.5069999695')


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
    f = join(dirpath, "4/label_threshold_1.9.mol2")
else:
    f = "4/label_threshold_1.9.mol2"

cmd.load(f, 'label_threshold_1.9')
cmd.hide('everything', 'label_threshold_1.9')
cmd.label("label_threshold_1.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [1.9]
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


cluster_dict = {"10.3339996338":[], "10.3339996338_arrows":[]}

cluster_dict["10.3339996338"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(45.0), float(48.5), float(10.5), float(1.0)]

cluster_dict["10.3339996338_arrows"] += cgo_arrow([45.0,48.5,10.5], [47.87,49.057,10.538], color="blue red", name="Arrows_10.3339996338_1")

cluster_dict["10.3339996338"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(43.4494100047), float(46.6783843293), float(10.8180872017), float(1.0)]


cluster_dict["10.3339996338"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(45.5), float(45.5), float(11.0), float(1.0)]

cluster_dict["10.3339996338_arrows"] += cgo_arrow([45.5,45.5,11.0], [44.739,42.649,13.009], color="red blue", name="Arrows_10.3339996338_2")

cluster_dict["10.3339996338"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(45.5), float(46.0), float(10.0), float(1.0)]

cluster_dict["10.3339996338_arrows"] += cgo_arrow([45.5,46.0,10.0], [48.134,47.667,8.79], color="red blue", name="Arrows_10.3339996338_3")

cmd.load_cgo(cluster_dict["10.3339996338"], "Features_10.3339996338", 1)
cmd.load_cgo(cluster_dict["10.3339996338_arrows"], "Arrows_10.3339996338")
cmd.set("transparency", 0.2,"Features_10.3339996338")
cmd.group("Pharmacophore_10.3339996338", members="Features_10.3339996338")
cmd.group("Pharmacophore_10.3339996338", members="Arrows_10.3339996338")

if dirpath:
    f = join(dirpath, "4/label_threshold_10.3339996338.mol2")
else:
    f = "4/label_threshold_10.3339996338.mol2"

cmd.load(f, 'label_threshold_10.3339996338')
cmd.hide('everything', 'label_threshold_10.3339996338')
cmd.label("label_threshold_10.3339996338", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.3339996338', members= 'label_threshold_10.3339996338')


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
    f = join(dirpath, "5/label_threshold_0.6.mol2")
else:
    f = "5/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
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


cluster_dict = {"10.1035003662":[], "10.1035003662_arrows":[]}

cluster_dict["10.1035003662"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(46.5), float(29.5), float(9.0), float(1.0)]

cluster_dict["10.1035003662_arrows"] += cgo_arrow([46.5,29.5,9.0], [46.213,31.336,6.479], color="blue red", name="Arrows_10.1035003662_1")

cluster_dict["10.1035003662"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(44.3249523871), float(29.5018454316), float(10.864817071), float(1.0)]


cluster_dict["10.1035003662"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(46.2952582098), float(28.1735753337), float(13.1880835216), float(1.0)]


cluster_dict["10.1035003662"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(46.68182728), float(28.7962660165), float(13.1933395645), float(1.0)]


cluster_dict["10.1035003662"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(48.2323977162), float(28.6516650793), float(8.99329210121), float(1.0)]


cluster_dict["10.1035003662"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(46.5), float(30.5), float(9.5), float(1.0)]

cluster_dict["10.1035003662_arrows"] += cgo_arrow([46.5,30.5,9.5], [43.976,33.707,7.915], color="red blue", name="Arrows_10.1035003662_2")

cluster_dict["10.1035003662"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(46.5), float(30.5), float(9.5), float(1.0)]

cluster_dict["10.1035003662_arrows"] += cgo_arrow([46.5,30.5,9.5], [43.976,33.707,7.915], color="red blue", name="Arrows_10.1035003662_3")

cmd.load_cgo(cluster_dict["10.1035003662"], "Features_10.1035003662", 1)
cmd.load_cgo(cluster_dict["10.1035003662_arrows"], "Arrows_10.1035003662")
cmd.set("transparency", 0.2,"Features_10.1035003662")
cmd.group("Pharmacophore_10.1035003662", members="Features_10.1035003662")
cmd.group("Pharmacophore_10.1035003662", members="Arrows_10.1035003662")

if dirpath:
    f = join(dirpath, "5/label_threshold_10.1035003662.mol2")
else:
    f = "5/label_threshold_10.1035003662.mol2"

cmd.load(f, 'label_threshold_10.1035003662')
cmd.hide('everything', 'label_threshold_10.1035003662')
cmd.label("label_threshold_10.1035003662", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.1035003662', members= 'label_threshold_10.1035003662')


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
    f = join(dirpath, "6/label_threshold_0.6.mol2")
else:
    f = "6/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
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


cluster_dict = {"9.30399990082":[], "9.30399990082_arrows":[]}

cluster_dict["9.30399990082"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(58.3371803824), float(29.3639712172), float(-14.5790803241), float(1.0)]


cluster_dict["9.30399990082"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(61.1980838774), float(34.158856875), float(-15.1197790311), float(1.0)]


cluster_dict["9.30399990082"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(55.5), float(27.0), float(-15.0), float(1.0)]

cluster_dict["9.30399990082_arrows"] += cgo_arrow([55.5,27.0,-15.0], [53.232,24.817,-14.54], color="red blue", name="Arrows_9.30399990082_1")

cmd.load_cgo(cluster_dict["9.30399990082"], "Features_9.30399990082", 1)
cmd.load_cgo(cluster_dict["9.30399990082_arrows"], "Arrows_9.30399990082")
cmd.set("transparency", 0.2,"Features_9.30399990082")
cmd.group("Pharmacophore_9.30399990082", members="Features_9.30399990082")
cmd.group("Pharmacophore_9.30399990082", members="Arrows_9.30399990082")

if dirpath:
    f = join(dirpath, "6/label_threshold_9.30399990082.mol2")
else:
    f = "6/label_threshold_9.30399990082.mol2"

cmd.load(f, 'label_threshold_9.30399990082')
cmd.hide('everything', 'label_threshold_9.30399990082')
cmd.label("label_threshold_9.30399990082", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.30399990082', members= 'label_threshold_9.30399990082')


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
    f = join(dirpath, "7/label_threshold_0.6.mol2")
else:
    f = "7/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
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


cluster_dict = {"8.3140001297":[], "8.3140001297_arrows":[]}

cluster_dict["8.3140001297"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(52.3222907861), float(53.6110647815), float(6.93610888053), float(1.0)]


cmd.load_cgo(cluster_dict["8.3140001297"], "Features_8.3140001297", 1)
cmd.load_cgo(cluster_dict["8.3140001297_arrows"], "Arrows_8.3140001297")
cmd.set("transparency", 0.2,"Features_8.3140001297")
cmd.group("Pharmacophore_8.3140001297", members="Features_8.3140001297")
cmd.group("Pharmacophore_8.3140001297", members="Arrows_8.3140001297")

if dirpath:
    f = join(dirpath, "7/label_threshold_8.3140001297.mol2")
else:
    f = "7/label_threshold_8.3140001297.mol2"

cmd.load(f, 'label_threshold_8.3140001297')
cmd.hide('everything', 'label_threshold_8.3140001297')
cmd.label("label_threshold_8.3140001297", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_8.3140001297', members= 'label_threshold_8.3140001297')


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
    f = join(dirpath, "8/label_threshold_30.0.mol2")
else:
    f = "8/label_threshold_30.0.mol2"

cmd.load(f, 'label_threshold_30.0')
cmd.hide('everything', 'label_threshold_30.0')
cmd.label("label_threshold_30.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [30.0]
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


cluster_dict = {"0":[], "0_arrows":[]}

cmd.load_cgo(cluster_dict["0"], "Features_0", 1)
cmd.load_cgo(cluster_dict["0_arrows"], "Arrows_0")
cmd.set("transparency", 0.2,"Features_0")
cmd.group("Pharmacophore_0", members="Features_0")
cmd.group("Pharmacophore_0", members="Arrows_0")

if dirpath:
    f = join(dirpath, "8/label_threshold_0.mol2")
else:
    f = "8/label_threshold_0.mol2"

cmd.load(f, 'label_threshold_0')
cmd.hide('everything', 'label_threshold_0')
cmd.label("label_threshold_0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_0', members= 'label_threshold_0')


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

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
