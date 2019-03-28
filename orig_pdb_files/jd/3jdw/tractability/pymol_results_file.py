
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
    f = join(dirpath, "0/label_threshold_1.5.mol2")
else:
    f = "0/label_threshold_1.5.mol2"

cmd.load(f, 'label_threshold_1.5')
cmd.hide('everything', 'label_threshold_1.5')
cmd.label("label_threshold_1.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [1.5]
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


cluster_dict = {"10.1459999084":[], "10.1459999084_arrows":[]}

cluster_dict["10.1459999084"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(35.5), float(74.0), float(7.0), float(1.0)]

cluster_dict["10.1459999084_arrows"] += cgo_arrow([35.5,74.0,7.0], [34.116,73.871,9.597], color="blue red", name="Arrows_10.1459999084_1")

cluster_dict["10.1459999084"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(34.6444879604), float(72.2906312222), float(4.6447164033), float(1.0)]


cluster_dict["10.1459999084"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(35.6666268606), float(74.3868666615), float(8.71691188344), float(1.0)]


cmd.load_cgo(cluster_dict["10.1459999084"], "Features_10.1459999084", 1)
cmd.load_cgo(cluster_dict["10.1459999084_arrows"], "Arrows_10.1459999084")
cmd.set("transparency", 0.2,"Features_10.1459999084")
cmd.group("Pharmacophore_10.1459999084", members="Features_10.1459999084")
cmd.group("Pharmacophore_10.1459999084", members="Arrows_10.1459999084")

if dirpath:
    f = join(dirpath, "0/label_threshold_10.1459999084.mol2")
else:
    f = "0/label_threshold_10.1459999084.mol2"

cmd.load(f, 'label_threshold_10.1459999084')
cmd.hide('everything', 'label_threshold_10.1459999084')
cmd.label("label_threshold_10.1459999084", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.1459999084', members= 'label_threshold_10.1459999084')


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
    f = join(dirpath, "1/label_threshold_11.5.mol2")
else:
    f = "1/label_threshold_11.5.mol2"

cmd.load(f, 'label_threshold_11.5')
cmd.hide('everything', 'label_threshold_11.5')
cmd.label("label_threshold_11.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [11.5]
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


cluster_dict = {"13.6309995651":[], "13.6309995651_arrows":[]}

cluster_dict["13.6309995651"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(46.0), float(73.5), float(14.5), float(1.0)]

cluster_dict["13.6309995651_arrows"] += cgo_arrow([46.0,73.5,14.5], [44.9,70.913,13.25], color="blue red", name="Arrows_13.6309995651_1")

cluster_dict["13.6309995651"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(44.0), float(73.0), float(19.0), float(1.0)]

cluster_dict["13.6309995651_arrows"] += cgo_arrow([44.0,73.0,19.0], [42.775,75.507,19.924], color="blue red", name="Arrows_13.6309995651_2")

cluster_dict["13.6309995651"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(46.0), float(73.5), float(14.5), float(1.0)]

cluster_dict["13.6309995651_arrows"] += cgo_arrow([46.0,73.5,14.5], [44.9,70.913,13.25], color="blue red", name="Arrows_13.6309995651_3")

cluster_dict["13.6309995651"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(47.5), float(70.0), float(17.0), float(1.0)]

cluster_dict["13.6309995651_arrows"] += cgo_arrow([47.5,70.0,17.0], [49.527,70.997,19.052], color="blue red", name="Arrows_13.6309995651_4")

cluster_dict["13.6309995651"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(47.0), float(80.0), float(21.0), float(1.0)]

cluster_dict["13.6309995651_arrows"] += cgo_arrow([47.0,80.0,21.0], [47.337,77.828,19.543], color="blue red", name="Arrows_13.6309995651_5")

cluster_dict["13.6309995651"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(47.5), float(81.0), float(24.5), float(1.0)]

cluster_dict["13.6309995651_arrows"] += cgo_arrow([47.5,81.0,24.5], [50.222,81.088,26.088], color="blue red", name="Arrows_13.6309995651_6")

cluster_dict["13.6309995651"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(38.1679747), float(72.893707006), float(8.81416075632), float(1.0)]


cluster_dict["13.6309995651"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(43.5494153858), float(81.9082845005), float(13.3695984721), float(1.0)]


cluster_dict["13.6309995651"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(45.4836703814), float(72.0240099425), float(16.2928234305), float(1.0)]


cluster_dict["13.6309995651"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(45.6736520546), float(81.6290606566), float(21.0710972221), float(1.0)]


cluster_dict["13.6309995651"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(45.0), float(73.5), float(17.5), float(1.0)]

cluster_dict["13.6309995651_arrows"] += cgo_arrow([45.0,73.5,17.5], [45.656,76.341,19.264], color="red blue", name="Arrows_13.6309995651_7")

cluster_dict["13.6309995651"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(44.0), float(77.5), float(16.0), float(1.0)]

cluster_dict["13.6309995651_arrows"] += cgo_arrow([44.0,77.5,16.0], [43.632,79.741,17.764], color="red blue", name="Arrows_13.6309995651_8")

cluster_dict["13.6309995651"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(48.0), float(81.0), float(21.5), float(1.0)]

cluster_dict["13.6309995651_arrows"] += cgo_arrow([48.0,81.0,21.5], [50.487,82.373,22.334], color="red blue", name="Arrows_13.6309995651_9")

cmd.load_cgo(cluster_dict["13.6309995651"], "Features_13.6309995651", 1)
cmd.load_cgo(cluster_dict["13.6309995651_arrows"], "Arrows_13.6309995651")
cmd.set("transparency", 0.2,"Features_13.6309995651")
cmd.group("Pharmacophore_13.6309995651", members="Features_13.6309995651")
cmd.group("Pharmacophore_13.6309995651", members="Arrows_13.6309995651")

if dirpath:
    f = join(dirpath, "1/label_threshold_13.6309995651.mol2")
else:
    f = "1/label_threshold_13.6309995651.mol2"

cmd.load(f, 'label_threshold_13.6309995651')
cmd.hide('everything', 'label_threshold_13.6309995651')
cmd.label("label_threshold_13.6309995651", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.6309995651', members= 'label_threshold_13.6309995651')


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
    f = join(dirpath, "2/label_threshold_11.5.mol2")
else:
    f = "2/label_threshold_11.5.mol2"

cmd.load(f, 'label_threshold_11.5')
cmd.hide('everything', 'label_threshold_11.5')
cmd.label("label_threshold_11.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [11.5]
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


cluster_dict = {"13.3924999237":[], "13.3924999237_arrows":[]}

cluster_dict["13.3924999237"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(46.0), float(73.5), float(14.5), float(1.0)]

cluster_dict["13.3924999237_arrows"] += cgo_arrow([46.0,73.5,14.5], [44.9,70.913,13.25], color="blue red", name="Arrows_13.3924999237_1")

cluster_dict["13.3924999237"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(44.0), float(73.0), float(19.0), float(1.0)]

cluster_dict["13.3924999237_arrows"] += cgo_arrow([44.0,73.0,19.0], [42.775,75.507,19.924], color="blue red", name="Arrows_13.3924999237_2")

cluster_dict["13.3924999237"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(44.0), float(83.5), float(24.0), float(1.0)]

cluster_dict["13.3924999237_arrows"] += cgo_arrow([44.0,83.5,24.0], [43.187,84.409,26.679], color="blue red", name="Arrows_13.3924999237_3")

cluster_dict["13.3924999237"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(46.0), float(73.5), float(14.5), float(1.0)]

cluster_dict["13.3924999237_arrows"] += cgo_arrow([46.0,73.5,14.5], [44.9,70.913,13.25], color="blue red", name="Arrows_13.3924999237_4")

cluster_dict["13.3924999237"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(47.0), float(70.0), float(16.5), float(1.0)]

cluster_dict["13.3924999237_arrows"] += cgo_arrow([47.0,70.0,16.5], [49.527,70.997,19.052], color="blue red", name="Arrows_13.3924999237_5")

cluster_dict["13.3924999237"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(47.0), float(80.0), float(21.0), float(1.0)]

cluster_dict["13.3924999237_arrows"] += cgo_arrow([47.0,80.0,21.0], [47.337,77.828,19.543], color="blue red", name="Arrows_13.3924999237_6")

cluster_dict["13.3924999237"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(38.1827902512), float(73.0490973433), float(8.48261487004), float(1.0)]


cluster_dict["13.3924999237"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(43.5520006274), float(81.9061198547), float(13.3627140045), float(1.0)]


cluster_dict["13.3924999237"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(45.4593222205), float(72.4587546543), float(16.0491081695), float(1.0)]


cluster_dict["13.3924999237"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(45.6611283489), float(81.8700532854), float(21.2386512082), float(1.0)]


cluster_dict["13.3924999237"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(45.0), float(73.5), float(17.5), float(1.0)]

cluster_dict["13.3924999237_arrows"] += cgo_arrow([45.0,73.5,17.5], [45.656,76.341,19.264], color="red blue", name="Arrows_13.3924999237_7")

cluster_dict["13.3924999237"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(44.0), float(77.5), float(16.0), float(1.0)]

cluster_dict["13.3924999237_arrows"] += cgo_arrow([44.0,77.5,16.0], [43.632,79.741,17.764], color="red blue", name="Arrows_13.3924999237_8")

cluster_dict["13.3924999237"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(48.0), float(81.0), float(21.5), float(1.0)]

cluster_dict["13.3924999237_arrows"] += cgo_arrow([48.0,81.0,21.5], [50.487,82.373,22.334], color="red blue", name="Arrows_13.3924999237_9")

cmd.load_cgo(cluster_dict["13.3924999237"], "Features_13.3924999237", 1)
cmd.load_cgo(cluster_dict["13.3924999237_arrows"], "Arrows_13.3924999237")
cmd.set("transparency", 0.2,"Features_13.3924999237")
cmd.group("Pharmacophore_13.3924999237", members="Features_13.3924999237")
cmd.group("Pharmacophore_13.3924999237", members="Arrows_13.3924999237")

if dirpath:
    f = join(dirpath, "2/label_threshold_13.3924999237.mol2")
else:
    f = "2/label_threshold_13.3924999237.mol2"

cmd.load(f, 'label_threshold_13.3924999237')
cmd.hide('everything', 'label_threshold_13.3924999237')
cmd.label("label_threshold_13.3924999237", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.3924999237', members= 'label_threshold_13.3924999237')


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
    f = join(dirpath, "3/label_threshold_10.8.mol2")
else:
    f = "3/label_threshold_10.8.mol2"

cmd.load(f, 'label_threshold_10.8')
cmd.hide('everything', 'label_threshold_10.8')
cmd.label("label_threshold_10.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.8]
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


cluster_dict = {"12.9770002365":[], "12.9770002365_arrows":[]}

cluster_dict["12.9770002365"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(45.0), float(73.5), float(15.0), float(1.0)]

cluster_dict["12.9770002365_arrows"] += cgo_arrow([45.0,73.5,15.0], [44.9,70.913,13.25], color="blue red", name="Arrows_12.9770002365_1")

cluster_dict["12.9770002365"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(44.0), float(73.0), float(19.0), float(1.0)]

cluster_dict["12.9770002365_arrows"] += cgo_arrow([44.0,73.0,19.0], [42.775,75.507,19.924], color="blue red", name="Arrows_12.9770002365_2")

cluster_dict["12.9770002365"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(44.0), float(83.5), float(24.0), float(1.0)]

cluster_dict["12.9770002365_arrows"] += cgo_arrow([44.0,83.5,24.0], [43.187,84.409,26.679], color="blue red", name="Arrows_12.9770002365_3")

cluster_dict["12.9770002365"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(45.0), float(73.5), float(15.0), float(1.0)]

cluster_dict["12.9770002365_arrows"] += cgo_arrow([45.0,73.5,15.0], [44.9,70.913,13.25], color="blue red", name="Arrows_12.9770002365_4")

cluster_dict["12.9770002365"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(47.5), float(70.5), float(17.5), float(1.0)]

cluster_dict["12.9770002365_arrows"] += cgo_arrow([47.5,70.5,17.5], [49.527,70.997,19.052], color="blue red", name="Arrows_12.9770002365_5")

cluster_dict["12.9770002365"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(47.0), float(80.0), float(21.0), float(1.0)]

cluster_dict["12.9770002365_arrows"] += cgo_arrow([47.0,80.0,21.0], [47.337,77.828,19.543], color="blue red", name="Arrows_12.9770002365_6")

cluster_dict["12.9770002365"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(43.2229013098), float(81.7089766431), float(13.8166532408), float(1.0)]


cluster_dict["12.9770002365"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(45.5347237744), float(73.0362818397), float(16.3699406417), float(1.0)]


cluster_dict["12.9770002365"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(45.6533621079), float(81.8812827172), float(21.4244892422), float(1.0)]


cluster_dict["12.9770002365"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(46.4995907516), float(77.4543780348), float(15.6820786126), float(1.0)]


cluster_dict["12.9770002365"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(45.0), float(73.5), float(17.5), float(1.0)]

cluster_dict["12.9770002365_arrows"] += cgo_arrow([45.0,73.5,17.5], [45.656,76.341,19.264], color="red blue", name="Arrows_12.9770002365_7")

cluster_dict["12.9770002365"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(44.0), float(77.5), float(16.0), float(1.0)]

cluster_dict["12.9770002365_arrows"] += cgo_arrow([44.0,77.5,16.0], [43.632,79.741,17.764], color="red blue", name="Arrows_12.9770002365_8")

cluster_dict["12.9770002365"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(48.0), float(81.0), float(21.5), float(1.0)]

cluster_dict["12.9770002365_arrows"] += cgo_arrow([48.0,81.0,21.5], [50.487,82.373,22.334], color="red blue", name="Arrows_12.9770002365_9")

cmd.load_cgo(cluster_dict["12.9770002365"], "Features_12.9770002365", 1)
cmd.load_cgo(cluster_dict["12.9770002365_arrows"], "Arrows_12.9770002365")
cmd.set("transparency", 0.2,"Features_12.9770002365")
cmd.group("Pharmacophore_12.9770002365", members="Features_12.9770002365")
cmd.group("Pharmacophore_12.9770002365", members="Arrows_12.9770002365")

if dirpath:
    f = join(dirpath, "3/label_threshold_12.9770002365.mol2")
else:
    f = "3/label_threshold_12.9770002365.mol2"

cmd.load(f, 'label_threshold_12.9770002365')
cmd.hide('everything', 'label_threshold_12.9770002365')
cmd.label("label_threshold_12.9770002365", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.9770002365', members= 'label_threshold_12.9770002365')


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
    f = join(dirpath, "4/label_threshold_0.6.mol2")
else:
    f = "4/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
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


cluster_dict = {"12.2390003204":[], "12.2390003204_arrows":[]}

cluster_dict["12.2390003204"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(44.0), float(83.5), float(24.0), float(1.0)]

cluster_dict["12.2390003204_arrows"] += cgo_arrow([44.0,83.5,24.0], [43.187,84.409,26.679], color="blue red", name="Arrows_12.2390003204_1")

cluster_dict["12.2390003204"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(47.0), float(80.0), float(21.0), float(1.0)]

cluster_dict["12.2390003204_arrows"] += cgo_arrow([47.0,80.0,21.0], [47.337,77.828,19.543], color="blue red", name="Arrows_12.2390003204_2")

cluster_dict["12.2390003204"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(47.5), float(81.5), float(25.0), float(1.0)]

cluster_dict["12.2390003204_arrows"] += cgo_arrow([47.5,81.5,25.0], [50.222,81.088,26.088], color="blue red", name="Arrows_12.2390003204_3")

cluster_dict["12.2390003204"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(41.2029636825), float(79.4179057008), float(18.5), float(1.0)]


cluster_dict["12.2390003204"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(45.5465613132), float(81.9045568936), float(21.9240051034), float(1.0)]


cluster_dict["12.2390003204"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(48.0), float(81.0), float(21.5), float(1.0)]

cluster_dict["12.2390003204_arrows"] += cgo_arrow([48.0,81.0,21.5], [50.487,82.373,22.334], color="red blue", name="Arrows_12.2390003204_4")

cmd.load_cgo(cluster_dict["12.2390003204"], "Features_12.2390003204", 1)
cmd.load_cgo(cluster_dict["12.2390003204_arrows"], "Arrows_12.2390003204")
cmd.set("transparency", 0.2,"Features_12.2390003204")
cmd.group("Pharmacophore_12.2390003204", members="Features_12.2390003204")
cmd.group("Pharmacophore_12.2390003204", members="Arrows_12.2390003204")

if dirpath:
    f = join(dirpath, "4/label_threshold_12.2390003204.mol2")
else:
    f = "4/label_threshold_12.2390003204.mol2"

cmd.load(f, 'label_threshold_12.2390003204')
cmd.hide('everything', 'label_threshold_12.2390003204')
cmd.label("label_threshold_12.2390003204", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.2390003204', members= 'label_threshold_12.2390003204')


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


cluster_dict = {"10.8260002136":[], "10.8260002136_arrows":[]}

cluster_dict["10.8260002136"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(67.0861822634), float(66.574957458), float(8.50935379786), float(1.0)]


cluster_dict["10.8260002136"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(65.0), float(69.0), float(10.0), float(1.0)]

cluster_dict["10.8260002136_arrows"] += cgo_arrow([65.0,69.0,10.0], [64.317,71.994,9.601], color="red blue", name="Arrows_10.8260002136_1")

cmd.load_cgo(cluster_dict["10.8260002136"], "Features_10.8260002136", 1)
cmd.load_cgo(cluster_dict["10.8260002136_arrows"], "Arrows_10.8260002136")
cmd.set("transparency", 0.2,"Features_10.8260002136")
cmd.group("Pharmacophore_10.8260002136", members="Features_10.8260002136")
cmd.group("Pharmacophore_10.8260002136", members="Arrows_10.8260002136")

if dirpath:
    f = join(dirpath, "5/label_threshold_10.8260002136.mol2")
else:
    f = "5/label_threshold_10.8260002136.mol2"

cmd.load(f, 'label_threshold_10.8260002136')
cmd.hide('everything', 'label_threshold_10.8260002136')
cmd.label("label_threshold_10.8260002136", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.8260002136', members= 'label_threshold_10.8260002136')


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
    f = join(dirpath, "6/label_threshold_1.5.mol2")
else:
    f = "6/label_threshold_1.5.mol2"

cmd.load(f, 'label_threshold_1.5')
cmd.hide('everything', 'label_threshold_1.5')
cmd.label("label_threshold_1.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [1.5]
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


cluster_dict = {"10.1459999084":[], "10.1459999084_arrows":[]}

cluster_dict["10.1459999084"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(35.5), float(74.0), float(7.0), float(1.0)]

cluster_dict["10.1459999084_arrows"] += cgo_arrow([35.5,74.0,7.0], [34.116,73.871,9.597], color="blue red", name="Arrows_10.1459999084_1")

cluster_dict["10.1459999084"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(34.6444879604), float(72.2906312222), float(4.6447164033), float(1.0)]


cluster_dict["10.1459999084"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(35.6666268606), float(74.3868666615), float(8.71691188344), float(1.0)]


cmd.load_cgo(cluster_dict["10.1459999084"], "Features_10.1459999084", 1)
cmd.load_cgo(cluster_dict["10.1459999084_arrows"], "Arrows_10.1459999084")
cmd.set("transparency", 0.2,"Features_10.1459999084")
cmd.group("Pharmacophore_10.1459999084", members="Features_10.1459999084")
cmd.group("Pharmacophore_10.1459999084", members="Arrows_10.1459999084")

if dirpath:
    f = join(dirpath, "6/label_threshold_10.1459999084.mol2")
else:
    f = "6/label_threshold_10.1459999084.mol2"

cmd.load(f, 'label_threshold_10.1459999084')
cmd.hide('everything', 'label_threshold_10.1459999084')
cmd.label("label_threshold_10.1459999084", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.1459999084', members= 'label_threshold_10.1459999084')


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

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
