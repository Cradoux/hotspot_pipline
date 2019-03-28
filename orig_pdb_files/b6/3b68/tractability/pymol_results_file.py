
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
    f = join(dirpath, "0/label_threshold_14.6.mol2")
else:
    f = "0/label_threshold_14.6.mol2"

cmd.load(f, 'label_threshold_14.6')
cmd.hide('everything', 'label_threshold_14.6')
cmd.label("label_threshold_14.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [14.6]
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


cluster_dict = {"17.8740005493":[], "17.8740005493_arrows":[]}

cluster_dict["17.8740005493"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(31.5), float(-2.0), float(3.5), float(1.0)]

cluster_dict["17.8740005493_arrows"] += cgo_arrow([31.5,-2.0,3.5], [33.071,0.429,3.28], color="blue red", name="Arrows_17.8740005493_1")

cluster_dict["17.8740005493"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.0), float(4.0), float(4.5), float(1.0)]

cluster_dict["17.8740005493_arrows"] += cgo_arrow([29.0,4.0,4.5], [31.0,4.435,2.491], color="blue red", name="Arrows_17.8740005493_2")

cluster_dict["17.8740005493"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.3065143859), float(2.4693032185), float(5.38323686298), float(1.0)]


cluster_dict["17.8740005493"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.8462707522), float(11.4922844082), float(8.73771758043), float(1.0)]


cluster_dict["17.8740005493"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(32.0), float(0.0), float(5.0), float(1.0)]


cluster_dict["17.8740005493"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(10.5), float(6.5), float(1.0)]

cluster_dict["17.8740005493_arrows"] += cgo_arrow([23.5,10.5,6.5], [26.091,9.42,6.004], color="red blue", name="Arrows_17.8740005493_3")

cluster_dict["17.8740005493"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(2.0), float(3.0), float(1.0)]


cluster_dict["17.8740005493"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(8.0), float(4.5), float(1.0)]

cluster_dict["17.8740005493_arrows"] += cgo_arrow([24.0,8.0,4.5], [26.091,9.42,6.004], color="red blue", name="Arrows_17.8740005493_4")

cluster_dict["17.8740005493"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(25.5), float(10.0), float(9.0), float(1.0)]

cluster_dict["17.8740005493_arrows"] += cgo_arrow([25.5,10.0,9.0], [26.091,9.42,6.004], color="red blue", name="Arrows_17.8740005493_5")

cluster_dict["17.8740005493"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(1.0), float(4.5), float(1.0)]

cluster_dict["17.8740005493_arrows"] += cgo_arrow([26.0,1.0,4.5], [28.704,-2.798,6.142], color="red blue", name="Arrows_17.8740005493_6")

cluster_dict["17.8740005493"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(2.0), float(2.0), float(1.0)]


cluster_dict["17.8740005493"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(0.0), float(5.5), float(1.0)]

cluster_dict["17.8740005493_arrows"] += cgo_arrow([29.5,0.0,5.5], [28.704,-2.798,6.142], color="red blue", name="Arrows_17.8740005493_7")

cmd.load_cgo(cluster_dict["17.8740005493"], "Features_17.8740005493", 1)
cmd.load_cgo(cluster_dict["17.8740005493_arrows"], "Arrows_17.8740005493")
cmd.set("transparency", 0.2,"Features_17.8740005493")
cmd.group("Pharmacophore_17.8740005493", members="Features_17.8740005493")
cmd.group("Pharmacophore_17.8740005493", members="Arrows_17.8740005493")

if dirpath:
    f = join(dirpath, "0/label_threshold_17.8740005493.mol2")
else:
    f = "0/label_threshold_17.8740005493.mol2"

cmd.load(f, 'label_threshold_17.8740005493')
cmd.hide('everything', 'label_threshold_17.8740005493')
cmd.label("label_threshold_17.8740005493", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.8740005493', members= 'label_threshold_17.8740005493')


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
    f = join(dirpath, "1/label_threshold_12.9.mol2")
else:
    f = "1/label_threshold_12.9.mol2"

cmd.load(f, 'label_threshold_12.9')
cmd.hide('everything', 'label_threshold_12.9')
cmd.label("label_threshold_12.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.9]
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


cluster_dict = {"17.5869998932":[], "17.5869998932_arrows":[]}

cluster_dict["17.5869998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(7.5), float(2.5), float(1.0)]

cluster_dict["17.5869998932_arrows"] += cgo_arrow([25.0,7.5,2.5], [26.569,10.052,3.904], color="blue red", name="Arrows_17.5869998932_1")

cluster_dict["17.5869998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(31.5), float(-2.0), float(3.5), float(1.0)]

cluster_dict["17.5869998932_arrows"] += cgo_arrow([31.5,-2.0,3.5], [33.071,0.429,3.28], color="blue red", name="Arrows_17.5869998932_2")

cluster_dict["17.5869998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.0), float(4.0), float(4.5), float(1.0)]

cluster_dict["17.5869998932_arrows"] += cgo_arrow([29.0,4.0,4.5], [31.0,4.435,2.491], color="blue red", name="Arrows_17.5869998932_3")

cluster_dict["17.5869998932"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.4022574756), float(2.37731946865), float(5.47814800412), float(1.0)]


cluster_dict["17.5869998932"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(30.6950287519), float(-1.53321002707), float(5.36748687294), float(1.0)]


cluster_dict["17.5869998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.5), float(4.5), float(3.0), float(1.0)]

cluster_dict["17.5869998932_arrows"] += cgo_arrow([21.5,4.5,3.0], [21.282,5.516,7.331], color="red blue", name="Arrows_17.5869998932_4")

cluster_dict["17.5869998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(2.0), float(3.0), float(1.0)]


cluster_dict["17.5869998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(8.0), float(4.5), float(1.0)]

cluster_dict["17.5869998932_arrows"] += cgo_arrow([24.0,8.0,4.5], [26.091,9.42,6.004], color="red blue", name="Arrows_17.5869998932_5")

cluster_dict["17.5869998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(1.0), float(4.5), float(1.0)]

cluster_dict["17.5869998932_arrows"] += cgo_arrow([26.0,1.0,4.5], [28.704,-2.798,6.142], color="red blue", name="Arrows_17.5869998932_6")

cluster_dict["17.5869998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(2.0), float(2.0), float(1.0)]


cluster_dict["17.5869998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(0.0), float(5.5), float(1.0)]

cluster_dict["17.5869998932_arrows"] += cgo_arrow([29.5,0.0,5.5], [28.704,-2.798,6.142], color="red blue", name="Arrows_17.5869998932_7")

cmd.load_cgo(cluster_dict["17.5869998932"], "Features_17.5869998932", 1)
cmd.load_cgo(cluster_dict["17.5869998932_arrows"], "Arrows_17.5869998932")
cmd.set("transparency", 0.2,"Features_17.5869998932")
cmd.group("Pharmacophore_17.5869998932", members="Features_17.5869998932")
cmd.group("Pharmacophore_17.5869998932", members="Arrows_17.5869998932")

if dirpath:
    f = join(dirpath, "1/label_threshold_17.5869998932.mol2")
else:
    f = "1/label_threshold_17.5869998932.mol2"

cmd.load(f, 'label_threshold_17.5869998932')
cmd.hide('everything', 'label_threshold_17.5869998932')
cmd.label("label_threshold_17.5869998932", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.5869998932', members= 'label_threshold_17.5869998932')


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
    f = join(dirpath, "2/label_threshold_7.1.mol2")
else:
    f = "2/label_threshold_7.1.mol2"

cmd.load(f, 'label_threshold_7.1')
cmd.hide('everything', 'label_threshold_7.1')
cmd.label("label_threshold_7.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [7.1]
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


cluster_dict = {"13.7324995995":[], "13.7324995995_arrows":[]}

cluster_dict["13.7324995995"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(20.5), float(16.5), float(6.5), float(1.0)]

cluster_dict["13.7324995995_arrows"] += cgo_arrow([20.5,16.5,6.5], [20.086,16.052,8.963], color="blue red", name="Arrows_13.7324995995_1")

cluster_dict["13.7324995995"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(14.5), float(8.0), float(1.0)]

cluster_dict["13.7324995995_arrows"] += cgo_arrow([22.0,14.5,8.0], [20.086,16.052,8.963], color="blue red", name="Arrows_13.7324995995_2")

cluster_dict["13.7324995995"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(25.5), float(14.5), float(5.5), float(1.0)]

cluster_dict["13.7324995995_arrows"] += cgo_arrow([25.5,14.5,5.5], [25.609,17.502,4.964], color="blue red", name="Arrows_13.7324995995_3")

cluster_dict["13.7324995995"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(25.5), float(8.5), float(2.5), float(1.0)]

cluster_dict["13.7324995995_arrows"] += cgo_arrow([25.5,8.5,2.5], [26.569,10.052,3.904], color="blue red", name="Arrows_13.7324995995_4")

cluster_dict["13.7324995995"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(23.8466757961), float(13.2285446627), float(7.17272046796), float(1.0)]


cluster_dict["13.7324995995"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(25.6909076901), float(8.70889710719), float(3.549359656), float(1.0)]


cluster_dict["13.7324995995"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(10.5), float(6.5), float(1.0)]

cluster_dict["13.7324995995_arrows"] += cgo_arrow([23.5,10.5,6.5], [26.091,9.42,6.004], color="red blue", name="Arrows_13.7324995995_5")

cluster_dict["13.7324995995"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(8.5), float(4.5), float(1.0)]

cluster_dict["13.7324995995_arrows"] += cgo_arrow([24.0,8.5,4.5], [26.091,9.42,6.004], color="red blue", name="Arrows_13.7324995995_6")

cluster_dict["13.7324995995"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(12.0), float(11.0), float(1.0)]

cluster_dict["13.7324995995_arrows"] += cgo_arrow([25.0,12.0,11.0], [28.53,12.868,12.24], color="red blue", name="Arrows_13.7324995995_7")

cluster_dict["13.7324995995"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(15.0), float(2.0), float(1.0)]

cluster_dict["13.7324995995_arrows"] += cgo_arrow([24.5,15.0,2.0], [27.494,16.157,2.344], color="red blue", name="Arrows_13.7324995995_8")

cluster_dict["13.7324995995"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(25.5), float(10.0), float(9.0), float(1.0)]

cluster_dict["13.7324995995_arrows"] += cgo_arrow([25.5,10.0,9.0], [26.091,9.42,6.004], color="red blue", name="Arrows_13.7324995995_9")

cmd.load_cgo(cluster_dict["13.7324995995"], "Features_13.7324995995", 1)
cmd.load_cgo(cluster_dict["13.7324995995_arrows"], "Arrows_13.7324995995")
cmd.set("transparency", 0.2,"Features_13.7324995995")
cmd.group("Pharmacophore_13.7324995995", members="Features_13.7324995995")
cmd.group("Pharmacophore_13.7324995995", members="Arrows_13.7324995995")

if dirpath:
    f = join(dirpath, "2/label_threshold_13.7324995995.mol2")
else:
    f = "2/label_threshold_13.7324995995.mol2"

cmd.load(f, 'label_threshold_13.7324995995')
cmd.hide('everything', 'label_threshold_13.7324995995')
cmd.label("label_threshold_13.7324995995", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.7324995995', members= 'label_threshold_13.7324995995')


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
    f = join(dirpath, "3/label_threshold_7.8.mol2")
else:
    f = "3/label_threshold_7.8.mol2"

cmd.load(f, 'label_threshold_7.8')
cmd.hide('everything', 'label_threshold_7.8')
cmd.label("label_threshold_7.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [7.8]
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


cluster_dict = {"11.0229997635":[], "11.0229997635_arrows":[]}

cluster_dict["11.0229997635"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(3.0), float(19.0), float(14.0), float(1.0)]

cluster_dict["11.0229997635_arrows"] += cgo_arrow([3.0,19.0,14.0], [1.104,19.888,15.522], color="blue red", name="Arrows_11.0229997635_1")

cluster_dict["11.0229997635"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(6.0), float(20.0), float(12.0), float(1.0)]

cluster_dict["11.0229997635_arrows"] += cgo_arrow([6.0,20.0,12.0], [8.67,20.706,10.175], color="blue red", name="Arrows_11.0229997635_2")

cluster_dict["11.0229997635"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(4.5389380239), float(20.2690179095), float(13.9011426552), float(1.0)]


cluster_dict["11.0229997635"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(2.5), float(15.5), float(11.0), float(1.0)]

cluster_dict["11.0229997635_arrows"] += cgo_arrow([2.5,15.5,11.0], [0.21,14.277,12.128], color="red blue", name="Arrows_11.0229997635_3")

cluster_dict["11.0229997635"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(6.0), float(18.5), float(12.5), float(1.0)]

cluster_dict["11.0229997635_arrows"] += cgo_arrow([6.0,18.5,12.5], [7.108,17.247,10.439], color="red blue", name="Arrows_11.0229997635_4")

cmd.load_cgo(cluster_dict["11.0229997635"], "Features_11.0229997635", 1)
cmd.load_cgo(cluster_dict["11.0229997635_arrows"], "Arrows_11.0229997635")
cmd.set("transparency", 0.2,"Features_11.0229997635")
cmd.group("Pharmacophore_11.0229997635", members="Features_11.0229997635")
cmd.group("Pharmacophore_11.0229997635", members="Arrows_11.0229997635")

if dirpath:
    f = join(dirpath, "3/label_threshold_11.0229997635.mol2")
else:
    f = "3/label_threshold_11.0229997635.mol2"

cmd.load(f, 'label_threshold_11.0229997635')
cmd.hide('everything', 'label_threshold_11.0229997635')
cmd.label("label_threshold_11.0229997635", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.0229997635', members= 'label_threshold_11.0229997635')


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
    f = join(dirpath, "4/label_threshold_6.9.mol2")
else:
    f = "4/label_threshold_6.9.mol2"

cmd.load(f, 'label_threshold_6.9')
cmd.hide('everything', 'label_threshold_6.9')
cmd.label("label_threshold_6.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [6.9]
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


cluster_dict = {"10.4720001221":[], "10.4720001221_arrows":[]}

cluster_dict["10.4720001221"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(3.0), float(19.0), float(14.0), float(1.0)]

cluster_dict["10.4720001221_arrows"] += cgo_arrow([3.0,19.0,14.0], [1.104,19.888,15.522], color="blue red", name="Arrows_10.4720001221_1")

cluster_dict["10.4720001221"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(6.0), float(20.0), float(12.0), float(1.0)]

cluster_dict["10.4720001221_arrows"] += cgo_arrow([6.0,20.0,12.0], [8.67,20.706,10.175], color="blue red", name="Arrows_10.4720001221_2")

cluster_dict["10.4720001221"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(5.11819011694), float(20.642907456), float(14.5760854424), float(1.0)]


cluster_dict["10.4720001221"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(2.5), float(15.5), float(11.0), float(1.0)]

cluster_dict["10.4720001221_arrows"] += cgo_arrow([2.5,15.5,11.0], [0.21,14.277,12.128], color="red blue", name="Arrows_10.4720001221_3")

cluster_dict["10.4720001221"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(6.0), float(18.5), float(12.5), float(1.0)]

cluster_dict["10.4720001221_arrows"] += cgo_arrow([6.0,18.5,12.5], [7.108,17.247,10.439], color="red blue", name="Arrows_10.4720001221_4")

cmd.load_cgo(cluster_dict["10.4720001221"], "Features_10.4720001221", 1)
cmd.load_cgo(cluster_dict["10.4720001221_arrows"], "Arrows_10.4720001221")
cmd.set("transparency", 0.2,"Features_10.4720001221")
cmd.group("Pharmacophore_10.4720001221", members="Features_10.4720001221")
cmd.group("Pharmacophore_10.4720001221", members="Arrows_10.4720001221")

if dirpath:
    f = join(dirpath, "4/label_threshold_10.4720001221.mol2")
else:
    f = "4/label_threshold_10.4720001221.mol2"

cmd.load(f, 'label_threshold_10.4720001221')
cmd.hide('everything', 'label_threshold_10.4720001221')
cmd.label("label_threshold_10.4720001221", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.4720001221', members= 'label_threshold_10.4720001221')


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


cluster_dict = {"10.2440004349":[], "10.2440004349_arrows":[]}

cluster_dict["10.2440004349"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(31.391403842), float(0.211211098522), float(26.6810261055), float(1.0)]


cluster_dict["10.2440004349"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(28.5), float(-1.0), float(25.0), float(1.0)]

cluster_dict["10.2440004349_arrows"] += cgo_arrow([28.5,-1.0,25.0], [29.026,-1.43,22.61], color="red blue", name="Arrows_10.2440004349_1")

cmd.load_cgo(cluster_dict["10.2440004349"], "Features_10.2440004349", 1)
cmd.load_cgo(cluster_dict["10.2440004349_arrows"], "Arrows_10.2440004349")
cmd.set("transparency", 0.2,"Features_10.2440004349")
cmd.group("Pharmacophore_10.2440004349", members="Features_10.2440004349")
cmd.group("Pharmacophore_10.2440004349", members="Arrows_10.2440004349")

if dirpath:
    f = join(dirpath, "5/label_threshold_10.2440004349.mol2")
else:
    f = "5/label_threshold_10.2440004349.mol2"

cmd.load(f, 'label_threshold_10.2440004349')
cmd.hide('everything', 'label_threshold_10.2440004349')
cmd.label("label_threshold_10.2440004349", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.2440004349', members= 'label_threshold_10.2440004349')


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
