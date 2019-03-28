
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
    f = join(dirpath, "0/label_threshold_8.0.mol2")
else:
    f = "0/label_threshold_8.0.mol2"

cmd.load(f, 'label_threshold_8.0')
cmd.hide('everything', 'label_threshold_8.0')
cmd.label("label_threshold_8.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [8.0]
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


cluster_dict = {"14.0049996376":[], "14.0049996376_arrows":[]}

cluster_dict["14.0049996376"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(14.0), float(3.5), float(15.5), float(1.0)]

cluster_dict["14.0049996376_arrows"] += cgo_arrow([14.0,3.5,15.5], [13.465,1.988,13.526], color="blue red", name="Arrows_14.0049996376_1")

cluster_dict["14.0049996376"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(16.5), float(2.0), float(10.0), float(1.0)]

cluster_dict["14.0049996376_arrows"] += cgo_arrow([16.5,2.0,10.0], [14.173,3.186,8.735], color="blue red", name="Arrows_14.0049996376_2")

cluster_dict["14.0049996376"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(17.1755463246), float(4.9450631864), float(13.4456497661), float(1.0)]


cluster_dict["14.0049996376"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(2.5), float(17.0), float(1.0)]

cluster_dict["14.0049996376_arrows"] += cgo_arrow([16.0,2.5,17.0], [17.464,-0.224,17.579], color="red blue", name="Arrows_14.0049996376_3")

cluster_dict["14.0049996376"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.5), float(3.5), float(12.0), float(1.0)]

cluster_dict["14.0049996376_arrows"] += cgo_arrow([15.5,3.5,12.0], [13.406,5.822,11.575], color="red blue", name="Arrows_14.0049996376_4")

cluster_dict["14.0049996376"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.0), float(1.0), float(15.0), float(1.0)]

cluster_dict["14.0049996376_arrows"] += cgo_arrow([17.0,1.0,15.0], [17.464,-0.224,17.579], color="red blue", name="Arrows_14.0049996376_5")

cluster_dict["14.0049996376"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(3.0), float(13.5), float(1.0)]

cluster_dict["14.0049996376_arrows"] += cgo_arrow([17.5,3.0,13.5], [20.254,3.203,14.512], color="red blue", name="Arrows_14.0049996376_6")

cluster_dict["14.0049996376"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(2.5), float(10.5), float(1.0)]

cluster_dict["14.0049996376_arrows"] += cgo_arrow([18.0,2.5,10.5], [18.378,3.68,7.501], color="red blue", name="Arrows_14.0049996376_7")

cmd.load_cgo(cluster_dict["14.0049996376"], "Features_14.0049996376", 1)
cmd.load_cgo(cluster_dict["14.0049996376_arrows"], "Arrows_14.0049996376")
cmd.set("transparency", 0.2,"Features_14.0049996376")
cmd.group("Pharmacophore_14.0049996376", members="Features_14.0049996376")
cmd.group("Pharmacophore_14.0049996376", members="Arrows_14.0049996376")

if dirpath:
    f = join(dirpath, "0/label_threshold_14.0049996376.mol2")
else:
    f = "0/label_threshold_14.0049996376.mol2"

cmd.load(f, 'label_threshold_14.0049996376')
cmd.hide('everything', 'label_threshold_14.0049996376')
cmd.label("label_threshold_14.0049996376", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.0049996376', members= 'label_threshold_14.0049996376')


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
    f = join(dirpath, "1/label_threshold_8.0.mol2")
else:
    f = "1/label_threshold_8.0.mol2"

cmd.load(f, 'label_threshold_8.0')
cmd.hide('everything', 'label_threshold_8.0')
cmd.label("label_threshold_8.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [8.0]
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


cluster_dict = {"14.0049996376":[], "14.0049996376_arrows":[]}

cluster_dict["14.0049996376"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(14.0), float(3.5), float(15.5), float(1.0)]

cluster_dict["14.0049996376_arrows"] += cgo_arrow([14.0,3.5,15.5], [13.465,1.988,13.526], color="blue red", name="Arrows_14.0049996376_1")

cluster_dict["14.0049996376"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(16.5), float(2.0), float(10.0), float(1.0)]

cluster_dict["14.0049996376_arrows"] += cgo_arrow([16.5,2.0,10.0], [14.173,3.186,8.735], color="blue red", name="Arrows_14.0049996376_2")

cluster_dict["14.0049996376"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(17.1755463246), float(4.9450631864), float(13.4456497661), float(1.0)]


cluster_dict["14.0049996376"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(2.5), float(17.0), float(1.0)]

cluster_dict["14.0049996376_arrows"] += cgo_arrow([16.0,2.5,17.0], [17.464,-0.224,17.579], color="red blue", name="Arrows_14.0049996376_3")

cluster_dict["14.0049996376"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.5), float(3.5), float(12.0), float(1.0)]

cluster_dict["14.0049996376_arrows"] += cgo_arrow([15.5,3.5,12.0], [13.406,5.822,11.575], color="red blue", name="Arrows_14.0049996376_4")

cluster_dict["14.0049996376"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.0), float(1.0), float(15.0), float(1.0)]

cluster_dict["14.0049996376_arrows"] += cgo_arrow([17.0,1.0,15.0], [17.464,-0.224,17.579], color="red blue", name="Arrows_14.0049996376_5")

cluster_dict["14.0049996376"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(3.0), float(13.5), float(1.0)]

cluster_dict["14.0049996376_arrows"] += cgo_arrow([17.5,3.0,13.5], [20.254,3.203,14.512], color="red blue", name="Arrows_14.0049996376_6")

cluster_dict["14.0049996376"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(2.5), float(10.5), float(1.0)]

cluster_dict["14.0049996376_arrows"] += cgo_arrow([18.0,2.5,10.5], [18.378,3.68,7.501], color="red blue", name="Arrows_14.0049996376_7")

cmd.load_cgo(cluster_dict["14.0049996376"], "Features_14.0049996376", 1)
cmd.load_cgo(cluster_dict["14.0049996376_arrows"], "Arrows_14.0049996376")
cmd.set("transparency", 0.2,"Features_14.0049996376")
cmd.group("Pharmacophore_14.0049996376", members="Features_14.0049996376")
cmd.group("Pharmacophore_14.0049996376", members="Arrows_14.0049996376")

if dirpath:
    f = join(dirpath, "1/label_threshold_14.0049996376.mol2")
else:
    f = "1/label_threshold_14.0049996376.mol2"

cmd.load(f, 'label_threshold_14.0049996376')
cmd.hide('everything', 'label_threshold_14.0049996376')
cmd.label("label_threshold_14.0049996376", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.0049996376', members= 'label_threshold_14.0049996376')


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
    f = join(dirpath, "2/label_threshold_0.6.mol2")
else:
    f = "2/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
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


cluster_dict = {"11.935500145":[], "11.935500145_arrows":[]}

cluster_dict["11.935500145"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(9.0), float(-7.5), float(-2.5), float(1.0)]

cluster_dict["11.935500145_arrows"] += cgo_arrow([9.0,-7.5,-2.5], [7.021,-6.273,-0.862], color="blue red", name="Arrows_11.935500145_1")

cluster_dict["11.935500145"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(12.5), float(-7.5), float(-1.0), float(1.0)]

cluster_dict["11.935500145_arrows"] += cgo_arrow([12.5,-7.5,-1.0], [11.875,-9.22,1.403], color="blue red", name="Arrows_11.935500145_2")

cluster_dict["11.935500145"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(11.6144974652), float(-8.70469606357), float(-2.95732415569), float(1.0)]


cluster_dict["11.935500145"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(14.0), float(-7.5), float(-4.0), float(1.0)]

cluster_dict["11.935500145_arrows"] += cgo_arrow([14.0,-7.5,-4.0], [17.434,-9.25,-5.548], color="red blue", name="Arrows_11.935500145_3")

cmd.load_cgo(cluster_dict["11.935500145"], "Features_11.935500145", 1)
cmd.load_cgo(cluster_dict["11.935500145_arrows"], "Arrows_11.935500145")
cmd.set("transparency", 0.2,"Features_11.935500145")
cmd.group("Pharmacophore_11.935500145", members="Features_11.935500145")
cmd.group("Pharmacophore_11.935500145", members="Arrows_11.935500145")

if dirpath:
    f = join(dirpath, "2/label_threshold_11.935500145.mol2")
else:
    f = "2/label_threshold_11.935500145.mol2"

cmd.load(f, 'label_threshold_11.935500145')
cmd.hide('everything', 'label_threshold_11.935500145')
cmd.label("label_threshold_11.935500145", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.935500145', members= 'label_threshold_11.935500145')


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
    f = join(dirpath, "3/label_threshold_0.6.mol2")
else:
    f = "3/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
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


cluster_dict = {"9.58699989319":[], "9.58699989319_arrows":[]}

cluster_dict["9.58699989319"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(33.2038423552), float(1.02053464859), float(13.9506391644), float(1.0)]


cluster_dict["9.58699989319"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(36.6947238335), float(-2.38671762764), float(13.7107017629), float(1.0)]


cluster_dict["9.58699989319"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(37.1131731058), float(-0.787603346957), float(12.7159031208), float(1.0)]


cluster_dict["9.58699989319"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(35.0), float(-1.5), float(14.0), float(1.0)]

cluster_dict["9.58699989319_arrows"] += cgo_arrow([35.0,-1.5,14.0], [36.771,-1.884,16.386], color="red blue", name="Arrows_9.58699989319_1")

cmd.load_cgo(cluster_dict["9.58699989319"], "Features_9.58699989319", 1)
cmd.load_cgo(cluster_dict["9.58699989319_arrows"], "Arrows_9.58699989319")
cmd.set("transparency", 0.2,"Features_9.58699989319")
cmd.group("Pharmacophore_9.58699989319", members="Features_9.58699989319")
cmd.group("Pharmacophore_9.58699989319", members="Arrows_9.58699989319")

if dirpath:
    f = join(dirpath, "3/label_threshold_9.58699989319.mol2")
else:
    f = "3/label_threshold_9.58699989319.mol2"

cmd.load(f, 'label_threshold_9.58699989319')
cmd.hide('everything', 'label_threshold_9.58699989319')
cmd.label("label_threshold_9.58699989319", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.58699989319', members= 'label_threshold_9.58699989319')


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
    f = join(dirpath, "4/label_threshold_1.3.mol2")
else:
    f = "4/label_threshold_1.3.mol2"

cmd.load(f, 'label_threshold_1.3')
cmd.hide('everything', 'label_threshold_1.3')
cmd.label("label_threshold_1.3", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [1.3]
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


cluster_dict = {"8.43700027466":[], "8.43700027466_arrows":[]}

cluster_dict["8.43700027466"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.0558883447), float(-4.35318747273), float(5.05928305907), float(1.0)]


cluster_dict["8.43700027466"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.5), float(-4.0), float(7.0), float(1.0)]

cluster_dict["8.43700027466_arrows"] += cgo_arrow([27.5,-4.0,7.0], [29.655,-3.121,8.829], color="red blue", name="Arrows_8.43700027466_1")

cmd.load_cgo(cluster_dict["8.43700027466"], "Features_8.43700027466", 1)
cmd.load_cgo(cluster_dict["8.43700027466_arrows"], "Arrows_8.43700027466")
cmd.set("transparency", 0.2,"Features_8.43700027466")
cmd.group("Pharmacophore_8.43700027466", members="Features_8.43700027466")
cmd.group("Pharmacophore_8.43700027466", members="Arrows_8.43700027466")

if dirpath:
    f = join(dirpath, "4/label_threshold_8.43700027466.mol2")
else:
    f = "4/label_threshold_8.43700027466.mol2"

cmd.load(f, 'label_threshold_8.43700027466')
cmd.hide('everything', 'label_threshold_8.43700027466')
cmd.label("label_threshold_8.43700027466", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_8.43700027466', members= 'label_threshold_8.43700027466')


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


cluster_dict = {"7.56699991226":[], "7.56699991226_arrows":[]}

cluster_dict["7.56699991226"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(28.5767207638), float(-3.27311461058), float(24.3704578587), float(1.0)]


cmd.load_cgo(cluster_dict["7.56699991226"], "Features_7.56699991226", 1)
cmd.load_cgo(cluster_dict["7.56699991226_arrows"], "Arrows_7.56699991226")
cmd.set("transparency", 0.2,"Features_7.56699991226")
cmd.group("Pharmacophore_7.56699991226", members="Features_7.56699991226")
cmd.group("Pharmacophore_7.56699991226", members="Arrows_7.56699991226")

if dirpath:
    f = join(dirpath, "5/label_threshold_7.56699991226.mol2")
else:
    f = "5/label_threshold_7.56699991226.mol2"

cmd.load(f, 'label_threshold_7.56699991226')
cmd.hide('everything', 'label_threshold_7.56699991226')
cmd.label("label_threshold_7.56699991226", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_7.56699991226', members= 'label_threshold_7.56699991226')


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

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
