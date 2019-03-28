
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
    f = join(dirpath, "0/label_threshold_15.4.mol2")
else:
    f = "0/label_threshold_15.4.mol2"

cmd.load(f, 'label_threshold_15.4')
cmd.hide('everything', 'label_threshold_15.4')
cmd.label("label_threshold_15.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.4]
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


cluster_dict = {"16.4419994354":[], "16.4419994354_arrows":[]}

cluster_dict["16.4419994354"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(8.5), float(26.5), float(1.0)]

cluster_dict["16.4419994354_arrows"] += cgo_arrow([21.0,8.5,26.5], [20.263,11.401,25.376], color="blue red", name="Arrows_16.4419994354_1")

cluster_dict["16.4419994354"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(6.5), float(22.5), float(1.0)]

cluster_dict["16.4419994354_arrows"] += cgo_arrow([26.5,6.5,22.5], [24.383,5.1,20.493], color="blue red", name="Arrows_16.4419994354_2")

cluster_dict["16.4419994354"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(10.5), float(24.0), float(1.0)]

cluster_dict["16.4419994354_arrows"] += cgo_arrow([26.0,10.5,24.0], [24.828,13.48,24.46], color="blue red", name="Arrows_16.4419994354_3")

cluster_dict["16.4419994354"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(7.5), float(16.0), float(1.0)]

cluster_dict["16.4419994354_arrows"] += cgo_arrow([28.0,7.5,16.0], [25.678,8.126,17.598], color="blue red", name="Arrows_16.4419994354_4")

cluster_dict["16.4419994354"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(4.0), float(24.5), float(1.0)]

cluster_dict["16.4419994354_arrows"] += cgo_arrow([29.5,4.0,24.5], [31.281,4.479,26.704], color="blue red", name="Arrows_16.4419994354_5")

cluster_dict["16.4419994354"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.0), float(5.0), float(20.0), float(1.0)]

cluster_dict["16.4419994354_arrows"] += cgo_arrow([30.0,5.0,20.0], [31.989,4.998,17.834], color="blue red", name="Arrows_16.4419994354_6")

cluster_dict["16.4419994354"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.2203010927), float(7.56959229961), float(22.7641545129), float(1.0)]


cluster_dict["16.4419994354"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(35.3245614035), float(9.51169590643), float(17.7368421053), float(1.0)]


cluster_dict["16.4419994354"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(7.5), float(31.0), float(1.0)]

cluster_dict["16.4419994354_arrows"] += cgo_arrow([23.0,7.5,31.0], [24.327,5.898,29.592], color="red blue", name="Arrows_16.4419994354_7")

cluster_dict["16.4419994354"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(5.5), float(21.0), float(1.0)]

cluster_dict["16.4419994354_arrows"] += cgo_arrow([27.0,5.5,21.0], [24.383,5.1,20.493], color="red blue", name="Arrows_16.4419994354_8")

cluster_dict["16.4419994354"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(28.5), float(5.5), float(16.5), float(1.0)]

cluster_dict["16.4419994354_arrows"] += cgo_arrow([28.5,5.5,16.5], [30.035,4.383,14.053], color="red blue", name="Arrows_16.4419994354_9")

cmd.load_cgo(cluster_dict["16.4419994354"], "Features_16.4419994354", 1)
cmd.load_cgo(cluster_dict["16.4419994354_arrows"], "Arrows_16.4419994354")
cmd.set("transparency", 0.2,"Features_16.4419994354")
cmd.group("Pharmacophore_16.4419994354", members="Features_16.4419994354")
cmd.group("Pharmacophore_16.4419994354", members="Arrows_16.4419994354")

if dirpath:
    f = join(dirpath, "0/label_threshold_16.4419994354.mol2")
else:
    f = "0/label_threshold_16.4419994354.mol2"

cmd.load(f, 'label_threshold_16.4419994354')
cmd.hide('everything', 'label_threshold_16.4419994354')
cmd.label("label_threshold_16.4419994354", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.4419994354', members= 'label_threshold_16.4419994354')


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
    f = join(dirpath, "1/label_threshold_29.8.mol2")
else:
    f = "1/label_threshold_29.8.mol2"

cmd.load(f, 'label_threshold_29.8')
cmd.hide('everything', 'label_threshold_29.8')
cmd.label("label_threshold_29.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [29.8]
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


cluster_dict = {"27.6810007095":[], "27.6810007095_arrows":[]}

cluster_dict["27.6810007095"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-38.8041929116), float(-13.6382885499), float(30.4825558877), float(1.0)]


cluster_dict["27.6810007095"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-39.5), float(-13.0), float(27.0), float(1.0)]

cluster_dict["27.6810007095_arrows"] += cgo_arrow([-39.5,-13.0,27.0], [-40.302,-14.538,24.96], color="red blue", name="Arrows_27.6810007095_1")

cmd.load_cgo(cluster_dict["27.6810007095"], "Features_27.6810007095", 1)
cmd.load_cgo(cluster_dict["27.6810007095_arrows"], "Arrows_27.6810007095")
cmd.set("transparency", 0.2,"Features_27.6810007095")
cmd.group("Pharmacophore_27.6810007095", members="Features_27.6810007095")
cmd.group("Pharmacophore_27.6810007095", members="Arrows_27.6810007095")

if dirpath:
    f = join(dirpath, "1/label_threshold_27.6810007095.mol2")
else:
    f = "1/label_threshold_27.6810007095.mol2"

cmd.load(f, 'label_threshold_27.6810007095')
cmd.hide('everything', 'label_threshold_27.6810007095')
cmd.label("label_threshold_27.6810007095", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_27.6810007095', members= 'label_threshold_27.6810007095')


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
    f = join(dirpath, "2/label_threshold_29.9.mol2")
else:
    f = "2/label_threshold_29.9.mol2"

cmd.load(f, 'label_threshold_29.9')
cmd.hide('everything', 'label_threshold_29.9')
cmd.label("label_threshold_29.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [29.9]
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


cluster_dict = {"26.2000007629":[], "26.2000007629_arrows":[]}

cluster_dict["26.2000007629"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(8.08044336007), float(-8.94127604978), float(17.1242150732), float(1.0)]


cmd.load_cgo(cluster_dict["26.2000007629"], "Features_26.2000007629", 1)
cmd.load_cgo(cluster_dict["26.2000007629_arrows"], "Arrows_26.2000007629")
cmd.set("transparency", 0.2,"Features_26.2000007629")
cmd.group("Pharmacophore_26.2000007629", members="Features_26.2000007629")
cmd.group("Pharmacophore_26.2000007629", members="Arrows_26.2000007629")

if dirpath:
    f = join(dirpath, "2/label_threshold_26.2000007629.mol2")
else:
    f = "2/label_threshold_26.2000007629.mol2"

cmd.load(f, 'label_threshold_26.2000007629')
cmd.hide('everything', 'label_threshold_26.2000007629')
cmd.label("label_threshold_26.2000007629", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_26.2000007629', members= 'label_threshold_26.2000007629')


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
    f = join(dirpath, "3/label_threshold_30.0.mol2")
else:
    f = "3/label_threshold_30.0.mol2"

cmd.load(f, 'label_threshold_30.0')
cmd.hide('everything', 'label_threshold_30.0')
cmd.label("label_threshold_30.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [30.0]
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


cluster_dict = {"25.7539997101":[], "25.7539997101_arrows":[]}

cluster_dict["25.7539997101"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-6.29390903054), float(1.37945973546), float(14.6673279813), float(1.0)]


cmd.load_cgo(cluster_dict["25.7539997101"], "Features_25.7539997101", 1)
cmd.load_cgo(cluster_dict["25.7539997101_arrows"], "Arrows_25.7539997101")
cmd.set("transparency", 0.2,"Features_25.7539997101")
cmd.group("Pharmacophore_25.7539997101", members="Features_25.7539997101")
cmd.group("Pharmacophore_25.7539997101", members="Arrows_25.7539997101")

if dirpath:
    f = join(dirpath, "3/label_threshold_25.7539997101.mol2")
else:
    f = "3/label_threshold_25.7539997101.mol2"

cmd.load(f, 'label_threshold_25.7539997101')
cmd.hide('everything', 'label_threshold_25.7539997101')
cmd.label("label_threshold_25.7539997101", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_25.7539997101', members= 'label_threshold_25.7539997101')


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
    f = join(dirpath, "4/label_threshold_28.9.mol2")
else:
    f = "4/label_threshold_28.9.mol2"

cmd.load(f, 'label_threshold_28.9')
cmd.hide('everything', 'label_threshold_28.9')
cmd.label("label_threshold_28.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [28.9]
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


cluster_dict = {"25.7060003281":[], "25.7060003281_arrows":[]}

cluster_dict["25.7060003281"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(3.02585518737), float(-25.1932771678), float(27.5191384089), float(1.0)]


cmd.load_cgo(cluster_dict["25.7060003281"], "Features_25.7060003281", 1)
cmd.load_cgo(cluster_dict["25.7060003281_arrows"], "Arrows_25.7060003281")
cmd.set("transparency", 0.2,"Features_25.7060003281")
cmd.group("Pharmacophore_25.7060003281", members="Features_25.7060003281")
cmd.group("Pharmacophore_25.7060003281", members="Arrows_25.7060003281")

if dirpath:
    f = join(dirpath, "4/label_threshold_25.7060003281.mol2")
else:
    f = "4/label_threshold_25.7060003281.mol2"

cmd.load(f, 'label_threshold_25.7060003281')
cmd.hide('everything', 'label_threshold_25.7060003281')
cmd.label("label_threshold_25.7060003281", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_25.7060003281', members= 'label_threshold_25.7060003281')


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
    f = join(dirpath, "5/label_threshold_29.1.mol2")
else:
    f = "5/label_threshold_29.1.mol2"

cmd.load(f, 'label_threshold_29.1')
cmd.hide('everything', 'label_threshold_29.1')
cmd.label("label_threshold_29.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [29.1]
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


cluster_dict = {"25.3819999695":[], "25.3819999695_arrows":[]}

cluster_dict["25.3819999695"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-23.5269558924), float(-19.9248243685), float(37.5352392294), float(1.0)]


cmd.load_cgo(cluster_dict["25.3819999695"], "Features_25.3819999695", 1)
cmd.load_cgo(cluster_dict["25.3819999695_arrows"], "Arrows_25.3819999695")
cmd.set("transparency", 0.2,"Features_25.3819999695")
cmd.group("Pharmacophore_25.3819999695", members="Features_25.3819999695")
cmd.group("Pharmacophore_25.3819999695", members="Arrows_25.3819999695")

if dirpath:
    f = join(dirpath, "5/label_threshold_25.3819999695.mol2")
else:
    f = "5/label_threshold_25.3819999695.mol2"

cmd.load(f, 'label_threshold_25.3819999695')
cmd.hide('everything', 'label_threshold_25.3819999695')
cmd.label("label_threshold_25.3819999695", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_25.3819999695', members= 'label_threshold_25.3819999695')


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
    f = join(dirpath, "6/label_threshold_18.4.mol2")
else:
    f = "6/label_threshold_18.4.mol2"

cmd.load(f, 'label_threshold_18.4')
cmd.hide('everything', 'label_threshold_18.4')
cmd.label("label_threshold_18.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [18.4]
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


cluster_dict = {"25.3159999847":[], "25.3159999847_arrows":[]}

cluster_dict["25.3159999847"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-38.0403299442), float(-4.23557956062), float(15.4502288612), float(1.0)]


cmd.load_cgo(cluster_dict["25.3159999847"], "Features_25.3159999847", 1)
cmd.load_cgo(cluster_dict["25.3159999847_arrows"], "Arrows_25.3159999847")
cmd.set("transparency", 0.2,"Features_25.3159999847")
cmd.group("Pharmacophore_25.3159999847", members="Features_25.3159999847")
cmd.group("Pharmacophore_25.3159999847", members="Arrows_25.3159999847")

if dirpath:
    f = join(dirpath, "6/label_threshold_25.3159999847.mol2")
else:
    f = "6/label_threshold_25.3159999847.mol2"

cmd.load(f, 'label_threshold_25.3159999847')
cmd.hide('everything', 'label_threshold_25.3159999847')
cmd.label("label_threshold_25.3159999847", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_25.3159999847', members= 'label_threshold_25.3159999847')


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
    f = join(dirpath, "7/label_threshold_30.0.mol2")
else:
    f = "7/label_threshold_30.0.mol2"

cmd.load(f, 'label_threshold_30.0')
cmd.hide('everything', 'label_threshold_30.0')
cmd.label("label_threshold_30.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [30.0]
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


cluster_dict = {"25.1299991608":[], "25.1299991608_arrows":[]}

cluster_dict["25.1299991608"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(3.91237316383), float(-24.8918291732), float(26.754673353), float(1.0)]


cmd.load_cgo(cluster_dict["25.1299991608"], "Features_25.1299991608", 1)
cmd.load_cgo(cluster_dict["25.1299991608_arrows"], "Arrows_25.1299991608")
cmd.set("transparency", 0.2,"Features_25.1299991608")
cmd.group("Pharmacophore_25.1299991608", members="Features_25.1299991608")
cmd.group("Pharmacophore_25.1299991608", members="Arrows_25.1299991608")

if dirpath:
    f = join(dirpath, "7/label_threshold_25.1299991608.mol2")
else:
    f = "7/label_threshold_25.1299991608.mol2"

cmd.load(f, 'label_threshold_25.1299991608')
cmd.hide('everything', 'label_threshold_25.1299991608')
cmd.label("label_threshold_25.1299991608", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_25.1299991608', members= 'label_threshold_25.1299991608')


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
    f = join(dirpath, "8/label_threshold_27.5.mol2")
else:
    f = "8/label_threshold_27.5.mol2"

cmd.load(f, 'label_threshold_27.5')
cmd.hide('everything', 'label_threshold_27.5')
cmd.label("label_threshold_27.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [27.5]
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


cluster_dict = {"24.2740001678":[], "24.2740001678_arrows":[]}

cluster_dict["24.2740001678"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-50.715439082), float(5.90364767065), float(31.7928813755), float(1.0)]


cmd.load_cgo(cluster_dict["24.2740001678"], "Features_24.2740001678", 1)
cmd.load_cgo(cluster_dict["24.2740001678_arrows"], "Arrows_24.2740001678")
cmd.set("transparency", 0.2,"Features_24.2740001678")
cmd.group("Pharmacophore_24.2740001678", members="Features_24.2740001678")
cmd.group("Pharmacophore_24.2740001678", members="Arrows_24.2740001678")

if dirpath:
    f = join(dirpath, "8/label_threshold_24.2740001678.mol2")
else:
    f = "8/label_threshold_24.2740001678.mol2"

cmd.load(f, 'label_threshold_24.2740001678')
cmd.hide('everything', 'label_threshold_24.2740001678')
cmd.label("label_threshold_24.2740001678", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_24.2740001678', members= 'label_threshold_24.2740001678')


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
    f = join(dirpath, "9/label_threshold_29.9.mol2")
else:
    f = "9/label_threshold_29.9.mol2"

cmd.load(f, 'label_threshold_29.9')
cmd.hide('everything', 'label_threshold_29.9')
cmd.label("label_threshold_29.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [29.9]
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


cluster_dict = {"24.267999649":[], "24.267999649_arrows":[]}

cluster_dict["24.267999649"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.6529325721), float(-12.0638237975), float(33.6578377989), float(1.0)]


cmd.load_cgo(cluster_dict["24.267999649"], "Features_24.267999649", 1)
cmd.load_cgo(cluster_dict["24.267999649_arrows"], "Arrows_24.267999649")
cmd.set("transparency", 0.2,"Features_24.267999649")
cmd.group("Pharmacophore_24.267999649", members="Features_24.267999649")
cmd.group("Pharmacophore_24.267999649", members="Arrows_24.267999649")

if dirpath:
    f = join(dirpath, "9/label_threshold_24.267999649.mol2")
else:
    f = "9/label_threshold_24.267999649.mol2"

cmd.load(f, 'label_threshold_24.267999649')
cmd.hide('everything', 'label_threshold_24.267999649')
cmd.label("label_threshold_24.267999649", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_24.267999649', members= 'label_threshold_24.267999649')


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
    f = join(dirpath, "10/label_threshold_29.9.mol2")
else:
    f = "10/label_threshold_29.9.mol2"

cmd.load(f, 'label_threshold_29.9')
cmd.hide('everything', 'label_threshold_29.9')
cmd.label("label_threshold_29.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [29.9]
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


cluster_dict = {"23.9659996033":[], "23.9659996033_arrows":[]}

cluster_dict["23.9659996033"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-51.6075353891), float(6.13103055353), float(31.8186906486), float(1.0)]


cmd.load_cgo(cluster_dict["23.9659996033"], "Features_23.9659996033", 1)
cmd.load_cgo(cluster_dict["23.9659996033_arrows"], "Arrows_23.9659996033")
cmd.set("transparency", 0.2,"Features_23.9659996033")
cmd.group("Pharmacophore_23.9659996033", members="Features_23.9659996033")
cmd.group("Pharmacophore_23.9659996033", members="Arrows_23.9659996033")

if dirpath:
    f = join(dirpath, "10/label_threshold_23.9659996033.mol2")
else:
    f = "10/label_threshold_23.9659996033.mol2"

cmd.load(f, 'label_threshold_23.9659996033')
cmd.hide('everything', 'label_threshold_23.9659996033')
cmd.label("label_threshold_23.9659996033", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_23.9659996033', members= 'label_threshold_23.9659996033')


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
    f = join(dirpath, "11/label_threshold_29.8.mol2")
else:
    f = "11/label_threshold_29.8.mol2"

cmd.load(f, 'label_threshold_29.8')
cmd.hide('everything', 'label_threshold_29.8')
cmd.label("label_threshold_29.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [29.8]
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


cluster_dict = {"23.9580001831":[], "23.9580001831_arrows":[]}

cluster_dict["23.9580001831"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.8391521381), float(-12.0110479997), float(33.547187864), float(1.0)]


cmd.load_cgo(cluster_dict["23.9580001831"], "Features_23.9580001831", 1)
cmd.load_cgo(cluster_dict["23.9580001831_arrows"], "Arrows_23.9580001831")
cmd.set("transparency", 0.2,"Features_23.9580001831")
cmd.group("Pharmacophore_23.9580001831", members="Features_23.9580001831")
cmd.group("Pharmacophore_23.9580001831", members="Arrows_23.9580001831")

if dirpath:
    f = join(dirpath, "11/label_threshold_23.9580001831.mol2")
else:
    f = "11/label_threshold_23.9580001831.mol2"

cmd.load(f, 'label_threshold_23.9580001831')
cmd.hide('everything', 'label_threshold_23.9580001831')
cmd.label("label_threshold_23.9580001831", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_23.9580001831', members= 'label_threshold_23.9580001831')


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

if dirpath:
    f = join(dirpath, "12/label_threshold_29.9.mol2")
else:
    f = "12/label_threshold_29.9.mol2"

cmd.load(f, 'label_threshold_29.9')
cmd.hide('everything', 'label_threshold_29.9')
cmd.label("label_threshold_29.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [29.9]
gfiles = ['12/donor.grd', '12/apolar.grd', '12/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 12
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


cluster_dict = {"23.841999054":[], "23.841999054_arrows":[]}

cluster_dict["23.841999054"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-51.7043814304), float(6.09349160577), float(31.8614441553), float(1.0)]


cmd.load_cgo(cluster_dict["23.841999054"], "Features_23.841999054", 1)
cmd.load_cgo(cluster_dict["23.841999054_arrows"], "Arrows_23.841999054")
cmd.set("transparency", 0.2,"Features_23.841999054")
cmd.group("Pharmacophore_23.841999054", members="Features_23.841999054")
cmd.group("Pharmacophore_23.841999054", members="Arrows_23.841999054")

if dirpath:
    f = join(dirpath, "12/label_threshold_23.841999054.mol2")
else:
    f = "12/label_threshold_23.841999054.mol2"

cmd.load(f, 'label_threshold_23.841999054')
cmd.hide('everything', 'label_threshold_23.841999054')
cmd.label("label_threshold_23.841999054", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_23.841999054', members= 'label_threshold_23.841999054')


if dirpath:
    f = join(dirpath, '12/mesh.grd')
else:
    f = '12/mesh.grd'
cmd.load(f, 'mesh_12')
cmd.isomesh("isomesh_12", "mesh_12", 0.9)
cmd.color("grey80", "isomesh_12")
cmd.set('transparency', 0.4, "isomesh_12")

cmd.group('hotspot_12', "isomesh_12")
cmd.group('hotspot_12', "mesh_12")

if dirpath:
    f = join(dirpath, "13/label_threshold_30.0.mol2")
else:
    f = "13/label_threshold_30.0.mol2"

cmd.load(f, 'label_threshold_30.0')
cmd.hide('everything', 'label_threshold_30.0')
cmd.label("label_threshold_30.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [30.0]
gfiles = ['13/donor.grd', '13/apolar.grd', '13/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 13
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


cluster_dict = {"23.7779998779":[], "23.7779998779_arrows":[]}

cluster_dict["23.7779998779"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-36.1355390256), float(-3.13896158601), float(14.03676015), float(1.0)]


cmd.load_cgo(cluster_dict["23.7779998779"], "Features_23.7779998779", 1)
cmd.load_cgo(cluster_dict["23.7779998779_arrows"], "Arrows_23.7779998779")
cmd.set("transparency", 0.2,"Features_23.7779998779")
cmd.group("Pharmacophore_23.7779998779", members="Features_23.7779998779")
cmd.group("Pharmacophore_23.7779998779", members="Arrows_23.7779998779")

if dirpath:
    f = join(dirpath, "13/label_threshold_23.7779998779.mol2")
else:
    f = "13/label_threshold_23.7779998779.mol2"

cmd.load(f, 'label_threshold_23.7779998779')
cmd.hide('everything', 'label_threshold_23.7779998779')
cmd.label("label_threshold_23.7779998779", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_23.7779998779', members= 'label_threshold_23.7779998779')


if dirpath:
    f = join(dirpath, '13/mesh.grd')
else:
    f = '13/mesh.grd'
cmd.load(f, 'mesh_13')
cmd.isomesh("isomesh_13", "mesh_13", 0.9)
cmd.color("grey80", "isomesh_13")
cmd.set('transparency', 0.4, "isomesh_13")

cmd.group('hotspot_13', "isomesh_13")
cmd.group('hotspot_13', "mesh_13")

if dirpath:
    f = join(dirpath, "14/label_threshold_29.9.mol2")
else:
    f = "14/label_threshold_29.9.mol2"

cmd.load(f, 'label_threshold_29.9')
cmd.hide('everything', 'label_threshold_29.9')
cmd.label("label_threshold_29.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [29.9]
gfiles = ['14/donor.grd', '14/apolar.grd', '14/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 14
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


cluster_dict = {"23.6700000763":[], "23.6700000763_arrows":[]}

cluster_dict["23.6700000763"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-5.98236849507), float(2.72436410511), float(13.450961029), float(1.0)]


cluster_dict["23.6700000763"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-11.281484589), float(10.5632662501), float(17.730728245), float(1.0)]


cmd.load_cgo(cluster_dict["23.6700000763"], "Features_23.6700000763", 1)
cmd.load_cgo(cluster_dict["23.6700000763_arrows"], "Arrows_23.6700000763")
cmd.set("transparency", 0.2,"Features_23.6700000763")
cmd.group("Pharmacophore_23.6700000763", members="Features_23.6700000763")
cmd.group("Pharmacophore_23.6700000763", members="Arrows_23.6700000763")

if dirpath:
    f = join(dirpath, "14/label_threshold_23.6700000763.mol2")
else:
    f = "14/label_threshold_23.6700000763.mol2"

cmd.load(f, 'label_threshold_23.6700000763')
cmd.hide('everything', 'label_threshold_23.6700000763')
cmd.label("label_threshold_23.6700000763", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_23.6700000763', members= 'label_threshold_23.6700000763')


if dirpath:
    f = join(dirpath, '14/mesh.grd')
else:
    f = '14/mesh.grd'
cmd.load(f, 'mesh_14')
cmd.isomesh("isomesh_14", "mesh_14", 0.9)
cmd.color("grey80", "isomesh_14")
cmd.set('transparency', 0.4, "isomesh_14")

cmd.group('hotspot_14', "isomesh_14")
cmd.group('hotspot_14', "mesh_14")

if dirpath:
    f = join(dirpath, "15/label_threshold_29.8.mol2")
else:
    f = "15/label_threshold_29.8.mol2"

cmd.load(f, 'label_threshold_29.8')
cmd.hide('everything', 'label_threshold_29.8')
cmd.label("label_threshold_29.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [29.8]
gfiles = ['15/donor.grd', '15/apolar.grd', '15/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 15
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


cluster_dict = {"23.6420001984":[], "23.6420001984_arrows":[]}

cluster_dict["23.6420001984"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-36.0556639107), float(-3.11042957812), float(13.9973178806), float(1.0)]


cmd.load_cgo(cluster_dict["23.6420001984"], "Features_23.6420001984", 1)
cmd.load_cgo(cluster_dict["23.6420001984_arrows"], "Arrows_23.6420001984")
cmd.set("transparency", 0.2,"Features_23.6420001984")
cmd.group("Pharmacophore_23.6420001984", members="Features_23.6420001984")
cmd.group("Pharmacophore_23.6420001984", members="Arrows_23.6420001984")

if dirpath:
    f = join(dirpath, "15/label_threshold_23.6420001984.mol2")
else:
    f = "15/label_threshold_23.6420001984.mol2"

cmd.load(f, 'label_threshold_23.6420001984')
cmd.hide('everything', 'label_threshold_23.6420001984')
cmd.label("label_threshold_23.6420001984", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_23.6420001984', members= 'label_threshold_23.6420001984')


if dirpath:
    f = join(dirpath, '15/mesh.grd')
else:
    f = '15/mesh.grd'
cmd.load(f, 'mesh_15')
cmd.isomesh("isomesh_15", "mesh_15", 0.9)
cmd.color("grey80", "isomesh_15")
cmd.set('transparency', 0.4, "isomesh_15")

cmd.group('hotspot_15', "isomesh_15")
cmd.group('hotspot_15', "mesh_15")

if dirpath:
    f = join(dirpath, "16/label_threshold_17.9.mol2")
else:
    f = "16/label_threshold_17.9.mol2"

cmd.load(f, 'label_threshold_17.9')
cmd.hide('everything', 'label_threshold_17.9')
cmd.label("label_threshold_17.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [17.9]
gfiles = ['16/donor.grd', '16/apolar.grd', '16/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 16
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


cluster_dict = {"21.7159996033":[], "21.7159996033_arrows":[]}

cluster_dict["21.7159996033"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-61.5590480947), float(8.17223258053), float(78.8173595514), float(1.0)]


cluster_dict["21.7159996033"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-59.5), float(9.0), float(81.0), float(1.0)]

cluster_dict["21.7159996033_arrows"] += cgo_arrow([-59.5,9.0,81.0], [-58.393,11.66,81.384], color="red blue", name="Arrows_21.7159996033_1")

cmd.load_cgo(cluster_dict["21.7159996033"], "Features_21.7159996033", 1)
cmd.load_cgo(cluster_dict["21.7159996033_arrows"], "Arrows_21.7159996033")
cmd.set("transparency", 0.2,"Features_21.7159996033")
cmd.group("Pharmacophore_21.7159996033", members="Features_21.7159996033")
cmd.group("Pharmacophore_21.7159996033", members="Arrows_21.7159996033")

if dirpath:
    f = join(dirpath, "16/label_threshold_21.7159996033.mol2")
else:
    f = "16/label_threshold_21.7159996033.mol2"

cmd.load(f, 'label_threshold_21.7159996033')
cmd.hide('everything', 'label_threshold_21.7159996033')
cmd.label("label_threshold_21.7159996033", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_21.7159996033', members= 'label_threshold_21.7159996033')


if dirpath:
    f = join(dirpath, '16/mesh.grd')
else:
    f = '16/mesh.grd'
cmd.load(f, 'mesh_16')
cmd.isomesh("isomesh_16", "mesh_16", 0.9)
cmd.color("grey80", "isomesh_16")
cmd.set('transparency', 0.4, "isomesh_16")

cmd.group('hotspot_16', "isomesh_16")
cmd.group('hotspot_16', "mesh_16")

if dirpath:
    f = join(dirpath, "17/label_threshold_14.6.mol2")
else:
    f = "17/label_threshold_14.6.mol2"

cmd.load(f, 'label_threshold_14.6')
cmd.hide('everything', 'label_threshold_14.6')
cmd.label("label_threshold_14.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [14.6]
gfiles = ['17/donor.grd', '17/apolar.grd', '17/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 17
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


cluster_dict = {"20.9440002441":[], "20.9440002441_arrows":[]}

cluster_dict["20.9440002441"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-8.5), float(2.0), float(19.5), float(1.0)]

cluster_dict["20.9440002441_arrows"] += cgo_arrow([-8.5,2.0,19.5], [-10.06,1.446,18.913], color="blue red", name="Arrows_20.9440002441_1")

cluster_dict["20.9440002441"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-11.764043916), float(10.4708266929), float(18.9786579928), float(1.0)]


cluster_dict["20.9440002441"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-7.59676543083), float(4.03882327262), float(14.097390648), float(1.0)]


cmd.load_cgo(cluster_dict["20.9440002441"], "Features_20.9440002441", 1)
cmd.load_cgo(cluster_dict["20.9440002441_arrows"], "Arrows_20.9440002441")
cmd.set("transparency", 0.2,"Features_20.9440002441")
cmd.group("Pharmacophore_20.9440002441", members="Features_20.9440002441")
cmd.group("Pharmacophore_20.9440002441", members="Arrows_20.9440002441")

if dirpath:
    f = join(dirpath, "17/label_threshold_20.9440002441.mol2")
else:
    f = "17/label_threshold_20.9440002441.mol2"

cmd.load(f, 'label_threshold_20.9440002441')
cmd.hide('everything', 'label_threshold_20.9440002441')
cmd.label("label_threshold_20.9440002441", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_20.9440002441', members= 'label_threshold_20.9440002441')


if dirpath:
    f = join(dirpath, '17/mesh.grd')
else:
    f = '17/mesh.grd'
cmd.load(f, 'mesh_17')
cmd.isomesh("isomesh_17", "mesh_17", 0.9)
cmd.color("grey80", "isomesh_17")
cmd.set('transparency', 0.4, "isomesh_17")

cmd.group('hotspot_17', "isomesh_17")
cmd.group('hotspot_17', "mesh_17")

if dirpath:
    f = join(dirpath, "18/label_threshold_23.7.mol2")
else:
    f = "18/label_threshold_23.7.mol2"

cmd.load(f, 'label_threshold_23.7')
cmd.hide('everything', 'label_threshold_23.7')
cmd.label("label_threshold_23.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [23.7]
gfiles = ['18/donor.grd', '18/apolar.grd', '18/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 18
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


cluster_dict = {"20.1259994507":[], "20.1259994507_arrows":[]}

cluster_dict["20.1259994507"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-38.0), float(10.0), float(23.5), float(1.0)]

cluster_dict["20.1259994507_arrows"] += cgo_arrow([-38.0,10.0,23.5], [-37.669,7.405,22.176], color="blue red", name="Arrows_20.1259994507_1")

cluster_dict["20.1259994507"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-52.0895817334), float(8.7866029574), float(29.0882561775), float(1.0)]


cluster_dict["20.1259994507"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-45.2890497982), float(13.3383916983), float(21.6172547527), float(1.0)]


cluster_dict["20.1259994507"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-51.328246033), float(18.1805289489), float(28.6744370379), float(1.0)]


cluster_dict["20.1259994507"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-47.8112324518), float(8.25576313172), float(29.1472126879), float(1.0)]


cmd.load_cgo(cluster_dict["20.1259994507"], "Features_20.1259994507", 1)
cmd.load_cgo(cluster_dict["20.1259994507_arrows"], "Arrows_20.1259994507")
cmd.set("transparency", 0.2,"Features_20.1259994507")
cmd.group("Pharmacophore_20.1259994507", members="Features_20.1259994507")
cmd.group("Pharmacophore_20.1259994507", members="Arrows_20.1259994507")

if dirpath:
    f = join(dirpath, "18/label_threshold_20.1259994507.mol2")
else:
    f = "18/label_threshold_20.1259994507.mol2"

cmd.load(f, 'label_threshold_20.1259994507')
cmd.hide('everything', 'label_threshold_20.1259994507')
cmd.label("label_threshold_20.1259994507", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_20.1259994507', members= 'label_threshold_20.1259994507')


if dirpath:
    f = join(dirpath, '18/mesh.grd')
else:
    f = '18/mesh.grd'
cmd.load(f, 'mesh_18')
cmd.isomesh("isomesh_18", "mesh_18", 0.9)
cmd.color("grey80", "isomesh_18")
cmd.set('transparency', 0.4, "isomesh_18")

cmd.group('hotspot_18', "isomesh_18")
cmd.group('hotspot_18', "mesh_18")

if dirpath:
    f = join(dirpath, "19/label_threshold_29.6.mol2")
else:
    f = "19/label_threshold_29.6.mol2"

cmd.load(f, 'label_threshold_29.6')
cmd.hide('everything', 'label_threshold_29.6')
cmd.label("label_threshold_29.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [29.6]
gfiles = ['19/donor.grd', '19/apolar.grd', '19/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 19
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


cluster_dict = {"20.0699996948":[], "20.0699996948_arrows":[]}

cluster_dict["20.0699996948"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(39.0935796651), float(-8.93337253144), float(8.0563313493), float(1.0)]


cmd.load_cgo(cluster_dict["20.0699996948"], "Features_20.0699996948", 1)
cmd.load_cgo(cluster_dict["20.0699996948_arrows"], "Arrows_20.0699996948")
cmd.set("transparency", 0.2,"Features_20.0699996948")
cmd.group("Pharmacophore_20.0699996948", members="Features_20.0699996948")
cmd.group("Pharmacophore_20.0699996948", members="Arrows_20.0699996948")

if dirpath:
    f = join(dirpath, "19/label_threshold_20.0699996948.mol2")
else:
    f = "19/label_threshold_20.0699996948.mol2"

cmd.load(f, 'label_threshold_20.0699996948')
cmd.hide('everything', 'label_threshold_20.0699996948')
cmd.label("label_threshold_20.0699996948", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_20.0699996948', members= 'label_threshold_20.0699996948')


if dirpath:
    f = join(dirpath, '19/mesh.grd')
else:
    f = '19/mesh.grd'
cmd.load(f, 'mesh_19')
cmd.isomesh("isomesh_19", "mesh_19", 0.9)
cmd.color("grey80", "isomesh_19")
cmd.set('transparency', 0.4, "isomesh_19")

cmd.group('hotspot_19', "isomesh_19")
cmd.group('hotspot_19', "mesh_19")

if dirpath:
    f = join(dirpath, "20/label_threshold_2.7.mol2")
else:
    f = "20/label_threshold_2.7.mol2"

cmd.load(f, 'label_threshold_2.7')
cmd.hide('everything', 'label_threshold_2.7')
cmd.label("label_threshold_2.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [2.7]
gfiles = ['20/donor.grd', '20/apolar.grd', '20/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 20
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


cluster_dict = {"18.7759990692":[], "18.7759990692_arrows":[]}

cluster_dict["18.7759990692"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(38.4999324163), float(2.7092025764), float(44.3757999714), float(1.0)]


cluster_dict["18.7759990692"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(38.5), float(2.0), float(40.5), float(1.0)]

cluster_dict["18.7759990692_arrows"] += cgo_arrow([38.5,2.0,40.5], [38.765,2.462,38.467], color="red blue", name="Arrows_18.7759990692_1")

cmd.load_cgo(cluster_dict["18.7759990692"], "Features_18.7759990692", 1)
cmd.load_cgo(cluster_dict["18.7759990692_arrows"], "Arrows_18.7759990692")
cmd.set("transparency", 0.2,"Features_18.7759990692")
cmd.group("Pharmacophore_18.7759990692", members="Features_18.7759990692")
cmd.group("Pharmacophore_18.7759990692", members="Arrows_18.7759990692")

if dirpath:
    f = join(dirpath, "20/label_threshold_18.7759990692.mol2")
else:
    f = "20/label_threshold_18.7759990692.mol2"

cmd.load(f, 'label_threshold_18.7759990692')
cmd.hide('everything', 'label_threshold_18.7759990692')
cmd.label("label_threshold_18.7759990692", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_18.7759990692', members= 'label_threshold_18.7759990692')


if dirpath:
    f = join(dirpath, '20/mesh.grd')
else:
    f = '20/mesh.grd'
cmd.load(f, 'mesh_20')
cmd.isomesh("isomesh_20", "mesh_20", 0.9)
cmd.color("grey80", "isomesh_20")
cmd.set('transparency', 0.4, "isomesh_20")

cmd.group('hotspot_20', "isomesh_20")
cmd.group('hotspot_20', "mesh_20")

if dirpath:
    f = join(dirpath, "21/label_threshold_21.5.mol2")
else:
    f = "21/label_threshold_21.5.mol2"

cmd.load(f, 'label_threshold_21.5')
cmd.hide('everything', 'label_threshold_21.5')
cmd.label("label_threshold_21.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [21.5]
gfiles = ['21/donor.grd', '21/apolar.grd', '21/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 21
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


cluster_dict = {"18.2880001068":[], "18.2880001068_arrows":[]}

cluster_dict["18.2880001068"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-38.0), float(10.0), float(23.5), float(1.0)]

cluster_dict["18.2880001068_arrows"] += cgo_arrow([-38.0,10.0,23.5], [-37.669,7.405,22.176], color="blue red", name="Arrows_18.2880001068_1")

cluster_dict["18.2880001068"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-51.0), float(8.35606546412), float(29.3088706584), float(1.0)]


cluster_dict["18.2880001068"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-43.7529163408), float(12.6797107898), float(21.6574395737), float(1.0)]


cluster_dict["18.2880001068"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-47.8067722619), float(8.2589823406), float(29.1420368087), float(1.0)]


cmd.load_cgo(cluster_dict["18.2880001068"], "Features_18.2880001068", 1)
cmd.load_cgo(cluster_dict["18.2880001068_arrows"], "Arrows_18.2880001068")
cmd.set("transparency", 0.2,"Features_18.2880001068")
cmd.group("Pharmacophore_18.2880001068", members="Features_18.2880001068")
cmd.group("Pharmacophore_18.2880001068", members="Arrows_18.2880001068")

if dirpath:
    f = join(dirpath, "21/label_threshold_18.2880001068.mol2")
else:
    f = "21/label_threshold_18.2880001068.mol2"

cmd.load(f, 'label_threshold_18.2880001068')
cmd.hide('everything', 'label_threshold_18.2880001068')
cmd.label("label_threshold_18.2880001068", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_18.2880001068', members= 'label_threshold_18.2880001068')


if dirpath:
    f = join(dirpath, '21/mesh.grd')
else:
    f = '21/mesh.grd'
cmd.load(f, 'mesh_21')
cmd.isomesh("isomesh_21", "mesh_21", 0.9)
cmd.color("grey80", "isomesh_21")
cmd.set('transparency', 0.4, "isomesh_21")

cmd.group('hotspot_21', "isomesh_21")
cmd.group('hotspot_21', "mesh_21")

if dirpath:
    f = join(dirpath, "22/label_threshold_2.7.mol2")
else:
    f = "22/label_threshold_2.7.mol2"

cmd.load(f, 'label_threshold_2.7')
cmd.hide('everything', 'label_threshold_2.7')
cmd.label("label_threshold_2.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [2.7]
gfiles = ['22/donor.grd', '22/apolar.grd', '22/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 22
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


cluster_dict = {"18.1180000305":[], "18.1180000305_arrows":[]}

cluster_dict["18.1180000305"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-60.1620573066), float(17.7286031751), float(50.3805654434), float(1.0)]


cluster_dict["18.1180000305"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-61.5), float(14.5), float(50.0), float(1.0)]

cluster_dict["18.1180000305_arrows"] += cgo_arrow([-61.5,14.5,50.0], [-61.651,12.379,49.581], color="red blue", name="Arrows_18.1180000305_1")

cmd.load_cgo(cluster_dict["18.1180000305"], "Features_18.1180000305", 1)
cmd.load_cgo(cluster_dict["18.1180000305_arrows"], "Arrows_18.1180000305")
cmd.set("transparency", 0.2,"Features_18.1180000305")
cmd.group("Pharmacophore_18.1180000305", members="Features_18.1180000305")
cmd.group("Pharmacophore_18.1180000305", members="Arrows_18.1180000305")

if dirpath:
    f = join(dirpath, "22/label_threshold_18.1180000305.mol2")
else:
    f = "22/label_threshold_18.1180000305.mol2"

cmd.load(f, 'label_threshold_18.1180000305')
cmd.hide('everything', 'label_threshold_18.1180000305')
cmd.label("label_threshold_18.1180000305", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_18.1180000305', members= 'label_threshold_18.1180000305')


if dirpath:
    f = join(dirpath, '22/mesh.grd')
else:
    f = '22/mesh.grd'
cmd.load(f, 'mesh_22')
cmd.isomesh("isomesh_22", "mesh_22", 0.9)
cmd.color("grey80", "isomesh_22")
cmd.set('transparency', 0.4, "isomesh_22")

cmd.group('hotspot_22', "isomesh_22")
cmd.group('hotspot_22', "mesh_22")

if dirpath:
    f = join(dirpath, "23/label_threshold_16.3.mol2")
else:
    f = "23/label_threshold_16.3.mol2"

cmd.load(f, 'label_threshold_16.3')
cmd.hide('everything', 'label_threshold_16.3')
cmd.label("label_threshold_16.3", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [16.3]
gfiles = ['23/donor.grd', '23/apolar.grd', '23/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 23
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


cluster_dict = {"18.0340003967":[], "18.0340003967_arrows":[]}

cluster_dict["18.0340003967"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-47.0), float(8.0), float(59.5), float(1.0)]

cluster_dict["18.0340003967_arrows"] += cgo_arrow([-47.0,8.0,59.5], [-49.244,7.68,61.514], color="blue red", name="Arrows_18.0340003967_1")

cluster_dict["18.0340003967"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-46.5), float(11.5), float(53.5), float(1.0)]

cluster_dict["18.0340003967_arrows"] += cgo_arrow([-46.5,11.5,53.5], [-45.246,12.959,55.728], color="blue red", name="Arrows_18.0340003967_2")

cluster_dict["18.0340003967"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-43.5), float(21.0), float(67.5), float(1.0)]

cluster_dict["18.0340003967_arrows"] += cgo_arrow([-43.5,21.0,67.5], [-43.345,21.71,70.432], color="blue red", name="Arrows_18.0340003967_3")

cluster_dict["18.0340003967"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-41.0), float(14.5), float(71.0), float(1.0)]

cluster_dict["18.0340003967_arrows"] += cgo_arrow([-41.0,14.5,71.0], [-43.662,14.908,72.712], color="blue red", name="Arrows_18.0340003967_4")

cluster_dict["18.0340003967"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-41.5), float(15.5), float(62.0), float(1.0)]

cluster_dict["18.0340003967_arrows"] += cgo_arrow([-41.5,15.5,62.0], [-38.82,15.592,63.892], color="blue red", name="Arrows_18.0340003967_5")

cluster_dict["18.0340003967"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-45.8126910745), float(13.8802324496), float(52.3255942073), float(1.0)]


cluster_dict["18.0340003967"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-42.8227072733), float(16.4520139675), float(64.8106791333), float(1.0)]


cluster_dict["18.0340003967"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-45.0), float(17.0), float(61.0), float(1.0)]

cluster_dict["18.0340003967_arrows"] += cgo_arrow([-45.0,17.0,61.0], [-47.29,14.896,62.214], color="red blue", name="Arrows_18.0340003967_6")

cluster_dict["18.0340003967"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-43.0), float(17.0), float(69.0), float(1.0)]

cluster_dict["18.0340003967_arrows"] += cgo_arrow([-43.0,17.0,69.0], [-44.825,17.307,71.758], color="red blue", name="Arrows_18.0340003967_7")

cluster_dict["18.0340003967"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-41.5), float(11.5), float(68.0), float(1.0)]

cluster_dict["18.0340003967_arrows"] += cgo_arrow([-41.5,11.5,68.0], [-42.298,9.428,66.293], color="red blue", name="Arrows_18.0340003967_8")

cluster_dict["18.0340003967"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-39.5), float(12.0), float(62.0), float(1.0)]

cluster_dict["18.0340003967_arrows"] += cgo_arrow([-39.5,12.0,62.0], [-37.759,9.32,60.96], color="red blue", name="Arrows_18.0340003967_9")

cmd.load_cgo(cluster_dict["18.0340003967"], "Features_18.0340003967", 1)
cmd.load_cgo(cluster_dict["18.0340003967_arrows"], "Arrows_18.0340003967")
cmd.set("transparency", 0.2,"Features_18.0340003967")
cmd.group("Pharmacophore_18.0340003967", members="Features_18.0340003967")
cmd.group("Pharmacophore_18.0340003967", members="Arrows_18.0340003967")

if dirpath:
    f = join(dirpath, "23/label_threshold_18.0340003967.mol2")
else:
    f = "23/label_threshold_18.0340003967.mol2"

cmd.load(f, 'label_threshold_18.0340003967')
cmd.hide('everything', 'label_threshold_18.0340003967')
cmd.label("label_threshold_18.0340003967", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_18.0340003967', members= 'label_threshold_18.0340003967')


if dirpath:
    f = join(dirpath, '23/mesh.grd')
else:
    f = '23/mesh.grd'
cmd.load(f, 'mesh_23')
cmd.isomesh("isomesh_23", "mesh_23", 0.9)
cmd.color("grey80", "isomesh_23")
cmd.set('transparency', 0.4, "isomesh_23")

cmd.group('hotspot_23', "isomesh_23")
cmd.group('hotspot_23', "mesh_23")

if dirpath:
    f = join(dirpath, "24/label_threshold_16.0.mol2")
else:
    f = "24/label_threshold_16.0.mol2"

cmd.load(f, 'label_threshold_16.0')
cmd.hide('everything', 'label_threshold_16.0')
cmd.label("label_threshold_16.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [16.0]
gfiles = ['24/donor.grd', '24/apolar.grd', '24/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 24
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


cluster_dict = {"17.6499996185":[], "17.6499996185_arrows":[]}

cluster_dict["17.6499996185"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-45.5), float(9.0), float(48.0), float(1.0)]

cluster_dict["17.6499996185_arrows"] += cgo_arrow([-45.5,9.0,48.0], [-46.33,6.271,48.222], color="blue red", name="Arrows_17.6499996185_1")

cluster_dict["17.6499996185"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-36.5), float(15.5), float(44.0), float(1.0)]

cluster_dict["17.6499996185_arrows"] += cgo_arrow([-36.5,15.5,44.0], [-36.562,17.765,45.674], color="blue red", name="Arrows_17.6499996185_2")

cluster_dict["17.6499996185"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-33.5), float(9.5), float(52.0), float(1.0)]

cluster_dict["17.6499996185_arrows"] += cgo_arrow([-33.5,9.5,52.0], [-33.337,6.587,52.943], color="blue red", name="Arrows_17.6499996185_3")

cluster_dict["17.6499996185"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-33.0), float(12.5), float(49.0), float(1.0)]

cluster_dict["17.6499996185_arrows"] += cgo_arrow([-33.0,12.5,49.0], [-31.181,13.279,46.18], color="blue red", name="Arrows_17.6499996185_4")

cluster_dict["17.6499996185"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-31.0), float(11.0), float(48.0), float(1.0)]

cluster_dict["17.6499996185_arrows"] += cgo_arrow([-31.0,11.0,48.0], [-31.181,13.279,46.18], color="blue red", name="Arrows_17.6499996185_5")

cluster_dict["17.6499996185"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-44.07090425), float(11.1326111792), float(49.7462545888), float(1.0)]


cluster_dict["17.6499996185"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-43.1877422504), float(4.97951180796), float(54.7711602397), float(1.0)]


cluster_dict["17.6499996185"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-31.9226720588), float(10.9838266686), float(50.7295345256), float(1.0)]


cluster_dict["17.6499996185"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-28.1597633136), float(6.68934911243), float(43.0), float(1.0)]


cluster_dict["17.6499996185"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-45.5), float(8.5), float(49.5), float(1.0)]

cluster_dict["17.6499996185_arrows"] += cgo_arrow([-45.5,8.5,49.5], [-46.33,6.271,48.222], color="red blue", name="Arrows_17.6499996185_6")

cluster_dict["17.6499996185"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-32.5), float(11.5), float(55.0), float(1.0)]

cluster_dict["17.6499996185_arrows"] += cgo_arrow([-32.5,11.5,55.0], [-34.887,10.154,57.128], color="red blue", name="Arrows_17.6499996185_7")

cluster_dict["17.6499996185"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-32.5), float(12.0), float(48.0), float(1.0)]

cluster_dict["17.6499996185_arrows"] += cgo_arrow([-32.5,12.0,48.0], [-31.181,13.279,46.18], color="red blue", name="Arrows_17.6499996185_8")

cluster_dict["17.6499996185"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-27.0), float(12.5), float(51.5), float(1.0)]

cluster_dict["17.6499996185_arrows"] += cgo_arrow([-27.0,12.5,51.5], [-27.105,12.341,55.227], color="red blue", name="Arrows_17.6499996185_9")

cmd.load_cgo(cluster_dict["17.6499996185"], "Features_17.6499996185", 1)
cmd.load_cgo(cluster_dict["17.6499996185_arrows"], "Arrows_17.6499996185")
cmd.set("transparency", 0.2,"Features_17.6499996185")
cmd.group("Pharmacophore_17.6499996185", members="Features_17.6499996185")
cmd.group("Pharmacophore_17.6499996185", members="Arrows_17.6499996185")

if dirpath:
    f = join(dirpath, "24/label_threshold_17.6499996185.mol2")
else:
    f = "24/label_threshold_17.6499996185.mol2"

cmd.load(f, 'label_threshold_17.6499996185')
cmd.hide('everything', 'label_threshold_17.6499996185')
cmd.label("label_threshold_17.6499996185", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.6499996185', members= 'label_threshold_17.6499996185')


if dirpath:
    f = join(dirpath, '24/mesh.grd')
else:
    f = '24/mesh.grd'
cmd.load(f, 'mesh_24')
cmd.isomesh("isomesh_24", "mesh_24", 0.9)
cmd.color("grey80", "isomesh_24")
cmd.set('transparency', 0.4, "isomesh_24")

cmd.group('hotspot_24', "isomesh_24")
cmd.group('hotspot_24', "mesh_24")

if dirpath:
    f = join(dirpath, "25/label_threshold_14.9.mol2")
else:
    f = "25/label_threshold_14.9.mol2"

cmd.load(f, 'label_threshold_14.9')
cmd.hide('everything', 'label_threshold_14.9')
cmd.label("label_threshold_14.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [14.9]
gfiles = ['25/donor.grd', '25/apolar.grd', '25/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 25
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


cluster_dict = {"17.5769996643":[], "17.5769996643_arrows":[]}

cluster_dict["17.5769996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-47.0), float(8.0), float(59.5), float(1.0)]

cluster_dict["17.5769996643_arrows"] += cgo_arrow([-47.0,8.0,59.5], [-49.244,7.68,61.514], color="blue red", name="Arrows_17.5769996643_1")

cluster_dict["17.5769996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-43.5), float(17.5), float(59.5), float(1.0)]

cluster_dict["17.5769996643_arrows"] += cgo_arrow([-43.5,17.5,59.5], [-43.033,16.649,56.4], color="blue red", name="Arrows_17.5769996643_2")

cluster_dict["17.5769996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-43.5), float(21.0), float(67.5), float(1.0)]

cluster_dict["17.5769996643_arrows"] += cgo_arrow([-43.5,21.0,67.5], [-43.345,21.71,70.432], color="blue red", name="Arrows_17.5769996643_3")

cluster_dict["17.5769996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-41.0), float(14.5), float(71.0), float(1.0)]

cluster_dict["17.5769996643_arrows"] += cgo_arrow([-41.0,14.5,71.0], [-43.662,14.908,72.712], color="blue red", name="Arrows_17.5769996643_4")

cluster_dict["17.5769996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-41.5), float(15.5), float(62.0), float(1.0)]

cluster_dict["17.5769996643_arrows"] += cgo_arrow([-41.5,15.5,62.0], [-38.82,15.592,63.892], color="blue red", name="Arrows_17.5769996643_5")

cluster_dict["17.5769996643"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-45.6529850902), float(7.00620345603), float(59.8791551805), float(1.0)]


cluster_dict["17.5769996643"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-42.8005981666), float(16.4014651316), float(64.8189306991), float(1.0)]


cluster_dict["17.5769996643"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-38.0), float(5.5), float(67.0), float(1.0)]


cluster_dict["17.5769996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-43.0), float(17.0), float(69.0), float(1.0)]

cluster_dict["17.5769996643_arrows"] += cgo_arrow([-43.0,17.0,69.0], [-44.825,17.307,71.758], color="red blue", name="Arrows_17.5769996643_6")

cluster_dict["17.5769996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-45.0), float(17.0), float(61.0), float(1.0)]

cluster_dict["17.5769996643_arrows"] += cgo_arrow([-45.0,17.0,61.0], [-47.29,14.896,62.214], color="red blue", name="Arrows_17.5769996643_7")

cluster_dict["17.5769996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-41.5), float(11.5), float(68.0), float(1.0)]

cluster_dict["17.5769996643_arrows"] += cgo_arrow([-41.5,11.5,68.0], [-42.298,9.428,66.293], color="red blue", name="Arrows_17.5769996643_8")

cluster_dict["17.5769996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-41.0), float(18.5), float(63.5), float(1.0)]

cluster_dict["17.5769996643_arrows"] += cgo_arrow([-41.0,18.5,63.5], [-38.16,18.205,64.636], color="red blue", name="Arrows_17.5769996643_9")

cmd.load_cgo(cluster_dict["17.5769996643"], "Features_17.5769996643", 1)
cmd.load_cgo(cluster_dict["17.5769996643_arrows"], "Arrows_17.5769996643")
cmd.set("transparency", 0.2,"Features_17.5769996643")
cmd.group("Pharmacophore_17.5769996643", members="Features_17.5769996643")
cmd.group("Pharmacophore_17.5769996643", members="Arrows_17.5769996643")

if dirpath:
    f = join(dirpath, "25/label_threshold_17.5769996643.mol2")
else:
    f = "25/label_threshold_17.5769996643.mol2"

cmd.load(f, 'label_threshold_17.5769996643')
cmd.hide('everything', 'label_threshold_17.5769996643')
cmd.label("label_threshold_17.5769996643", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.5769996643', members= 'label_threshold_17.5769996643')


if dirpath:
    f = join(dirpath, '25/mesh.grd')
else:
    f = '25/mesh.grd'
cmd.load(f, 'mesh_25')
cmd.isomesh("isomesh_25", "mesh_25", 0.9)
cmd.color("grey80", "isomesh_25")
cmd.set('transparency', 0.4, "isomesh_25")

cmd.group('hotspot_25', "isomesh_25")
cmd.group('hotspot_25', "mesh_25")

if dirpath:
    f = join(dirpath, "26/label_threshold_15.5.mol2")
else:
    f = "26/label_threshold_15.5.mol2"

cmd.load(f, 'label_threshold_15.5')
cmd.hide('everything', 'label_threshold_15.5')
cmd.label("label_threshold_15.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.5]
gfiles = ['26/donor.grd', '26/apolar.grd', '26/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 26
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


cluster_dict = {"17.5580005646":[], "17.5580005646_arrows":[]}

cluster_dict["17.5580005646"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-48.5), float(6.5), float(59.0), float(1.0)]

cluster_dict["17.5580005646_arrows"] += cgo_arrow([-48.5,6.5,59.0], [-49.244,7.68,61.514], color="blue red", name="Arrows_17.5580005646_1")

cluster_dict["17.5580005646"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-47.0), float(8.0), float(59.5), float(1.0)]

cluster_dict["17.5580005646_arrows"] += cgo_arrow([-47.0,8.0,59.5], [-49.244,7.68,61.514], color="blue red", name="Arrows_17.5580005646_2")

cluster_dict["17.5580005646"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-46.0), float(9.0), float(48.0), float(1.0)]

cluster_dict["17.5580005646_arrows"] += cgo_arrow([-46.0,9.0,48.0], [-47.285,10.878,46.617], color="blue red", name="Arrows_17.5580005646_3")

cluster_dict["17.5580005646"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-46.0), float(9.5), float(45.0), float(1.0)]

cluster_dict["17.5580005646_arrows"] += cgo_arrow([-46.0,9.5,45.0], [-47.285,10.878,46.617], color="blue red", name="Arrows_17.5580005646_4")

cluster_dict["17.5580005646"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-46.5), float(11.5), float(53.5), float(1.0)]

cluster_dict["17.5580005646_arrows"] += cgo_arrow([-46.5,11.5,53.5], [-45.246,12.959,55.728], color="blue red", name="Arrows_17.5580005646_5")

cluster_dict["17.5580005646"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-47.0), float(18.0), float(51.0), float(1.0)]

cluster_dict["17.5580005646_arrows"] += cgo_arrow([-47.0,18.0,51.0], [-50.267,18.042,52.078], color="blue red", name="Arrows_17.5580005646_6")

cluster_dict["17.5580005646"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-45.0919376864), float(12.059570209), float(50.3808698713), float(1.0)]


cluster_dict["17.5580005646"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-45.1571951698), float(4.91333020259), float(55.3754649263), float(1.0)]


cluster_dict["17.5580005646"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-44.1856812337), float(17.0815697325), float(59.3428128632), float(1.0)]


cluster_dict["17.5580005646"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-48.0), float(19.0), float(51.0), float(1.0)]

cluster_dict["17.5580005646_arrows"] += cgo_arrow([-48.0,19.0,51.0], [-50.267,18.042,52.078], color="red blue", name="Arrows_17.5580005646_7")

cluster_dict["17.5580005646"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-45.5), float(8.5), float(49.5), float(1.0)]

cluster_dict["17.5580005646_arrows"] += cgo_arrow([-45.5,8.5,49.5], [-46.33,6.271,48.222], color="red blue", name="Arrows_17.5580005646_8")

cluster_dict["17.5580005646"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-44.0), float(15.5), float(59.5), float(1.0)]

cluster_dict["17.5580005646_arrows"] += cgo_arrow([-44.0,15.5,59.5], [-45.765,13.142,60.591], color="red blue", name="Arrows_17.5580005646_9")

cmd.load_cgo(cluster_dict["17.5580005646"], "Features_17.5580005646", 1)
cmd.load_cgo(cluster_dict["17.5580005646_arrows"], "Arrows_17.5580005646")
cmd.set("transparency", 0.2,"Features_17.5580005646")
cmd.group("Pharmacophore_17.5580005646", members="Features_17.5580005646")
cmd.group("Pharmacophore_17.5580005646", members="Arrows_17.5580005646")

if dirpath:
    f = join(dirpath, "26/label_threshold_17.5580005646.mol2")
else:
    f = "26/label_threshold_17.5580005646.mol2"

cmd.load(f, 'label_threshold_17.5580005646')
cmd.hide('everything', 'label_threshold_17.5580005646')
cmd.label("label_threshold_17.5580005646", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.5580005646', members= 'label_threshold_17.5580005646')


if dirpath:
    f = join(dirpath, '26/mesh.grd')
else:
    f = '26/mesh.grd'
cmd.load(f, 'mesh_26')
cmd.isomesh("isomesh_26", "mesh_26", 0.9)
cmd.color("grey80", "isomesh_26")
cmd.set('transparency', 0.4, "isomesh_26")

cmd.group('hotspot_26', "isomesh_26")
cmd.group('hotspot_26', "mesh_26")

if dirpath:
    f = join(dirpath, "27/label_threshold_16.1.mol2")
else:
    f = "27/label_threshold_16.1.mol2"

cmd.load(f, 'label_threshold_16.1')
cmd.hide('everything', 'label_threshold_16.1')
cmd.label("label_threshold_16.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [16.1]
gfiles = ['27/donor.grd', '27/apolar.grd', '27/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 27
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


cluster_dict = {"17.4335002899":[], "17.4335002899_arrows":[]}

cluster_dict["17.4335002899"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-55.0), float(-3.5), float(49.0), float(1.0)]

cluster_dict["17.4335002899_arrows"] += cgo_arrow([-55.0,-3.5,49.0], [-56.059,-0.837,49.615], color="blue red", name="Arrows_17.4335002899_1")

cluster_dict["17.4335002899"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-48.0), float(-4.0), float(48.5), float(1.0)]

cluster_dict["17.4335002899_arrows"] += cgo_arrow([-48.0,-4.0,48.5], [-48.278,-5.299,45.404], color="blue red", name="Arrows_17.4335002899_2")

cluster_dict["17.4335002899"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-47.0), float(8.0), float(59.5), float(1.0)]

cluster_dict["17.4335002899_arrows"] += cgo_arrow([-47.0,8.0,59.5], [-49.244,7.68,61.514], color="blue red", name="Arrows_17.4335002899_3")

cluster_dict["17.4335002899"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-46.0), float(9.0), float(48.0), float(1.0)]

cluster_dict["17.4335002899_arrows"] += cgo_arrow([-46.0,9.0,48.0], [-47.285,10.878,46.617], color="blue red", name="Arrows_17.4335002899_4")

cluster_dict["17.4335002899"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-46.0), float(9.5), float(45.0), float(1.0)]

cluster_dict["17.4335002899_arrows"] += cgo_arrow([-46.0,9.5,45.0], [-47.285,10.878,46.617], color="blue red", name="Arrows_17.4335002899_5")

cluster_dict["17.4335002899"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-46.5), float(11.5), float(53.5), float(1.0)]

cluster_dict["17.4335002899_arrows"] += cgo_arrow([-46.5,11.5,53.5], [-45.246,12.959,55.728], color="blue red", name="Arrows_17.4335002899_6")

cluster_dict["17.4335002899"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-51.0146379384), float(-4.8305837183), float(51.3245883589), float(1.0)]


cluster_dict["17.4335002899"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-44.8882262892), float(10.9565859382), float(50.5267943719), float(1.0)]


cluster_dict["17.4335002899"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-42.1093900997), float(12.8874540851), float(61.5179952277), float(1.0)]


cluster_dict["17.4335002899"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-45.5), float(8.5), float(49.5), float(1.0)]

cluster_dict["17.4335002899_arrows"] += cgo_arrow([-45.5,8.5,49.5], [-46.33,6.271,48.222], color="red blue", name="Arrows_17.4335002899_7")

cluster_dict["17.4335002899"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-44.0), float(15.0), float(59.5), float(1.0)]

cluster_dict["17.4335002899_arrows"] += cgo_arrow([-44.0,15.0,59.5], [-45.765,13.142,60.591], color="red blue", name="Arrows_17.4335002899_8")

cluster_dict["17.4335002899"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-40.0), float(11.5), float(61.5), float(1.0)]

cluster_dict["17.4335002899_arrows"] += cgo_arrow([-40.0,11.5,61.5], [-37.759,9.32,60.96], color="red blue", name="Arrows_17.4335002899_9")

cmd.load_cgo(cluster_dict["17.4335002899"], "Features_17.4335002899", 1)
cmd.load_cgo(cluster_dict["17.4335002899_arrows"], "Arrows_17.4335002899")
cmd.set("transparency", 0.2,"Features_17.4335002899")
cmd.group("Pharmacophore_17.4335002899", members="Features_17.4335002899")
cmd.group("Pharmacophore_17.4335002899", members="Arrows_17.4335002899")

if dirpath:
    f = join(dirpath, "27/label_threshold_17.4335002899.mol2")
else:
    f = "27/label_threshold_17.4335002899.mol2"

cmd.load(f, 'label_threshold_17.4335002899')
cmd.hide('everything', 'label_threshold_17.4335002899')
cmd.label("label_threshold_17.4335002899", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.4335002899', members= 'label_threshold_17.4335002899')


if dirpath:
    f = join(dirpath, '27/mesh.grd')
else:
    f = '27/mesh.grd'
cmd.load(f, 'mesh_27')
cmd.isomesh("isomesh_27", "mesh_27", 0.9)
cmd.color("grey80", "isomesh_27")
cmd.set('transparency', 0.4, "isomesh_27")

cmd.group('hotspot_27', "isomesh_27")
cmd.group('hotspot_27', "mesh_27")

if dirpath:
    f = join(dirpath, "28/label_threshold_16.0.mol2")
else:
    f = "28/label_threshold_16.0.mol2"

cmd.load(f, 'label_threshold_16.0')
cmd.hide('everything', 'label_threshold_16.0')
cmd.label("label_threshold_16.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [16.0]
gfiles = ['28/donor.grd', '28/apolar.grd', '28/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 28
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


cluster_dict = {"17.2999992371":[], "17.2999992371_arrows":[]}

cluster_dict["17.2999992371"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(12.0), float(33.5), float(39.0), float(1.0)]

cluster_dict["17.2999992371_arrows"] += cgo_arrow([12.0,33.5,39.0], [11.415,33.832,35.865], color="blue red", name="Arrows_17.2999992371_1")

cluster_dict["17.2999992371"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(16.5), float(31.5), float(27.5), float(1.0)]

cluster_dict["17.2999992371_arrows"] += cgo_arrow([16.5,31.5,27.5], [18.584,34.583,29.301], color="blue red", name="Arrows_17.2999992371_2")

cluster_dict["17.2999992371"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(31.5), float(26.0), float(1.0)]

cluster_dict["17.2999992371_arrows"] += cgo_arrow([16.0,31.5,26.0], [14.264,28.282,24.68], color="blue red", name="Arrows_17.2999992371_3")

cluster_dict["17.2999992371"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(16.5), float(31.5), float(27.5), float(1.0)]

cluster_dict["17.2999992371_arrows"] += cgo_arrow([16.5,31.5,27.5], [18.584,34.583,29.301], color="blue red", name="Arrows_17.2999992371_4")

cluster_dict["17.2999992371"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(28.0), float(40.0), float(1.0)]

cluster_dict["17.2999992371_arrows"] += cgo_arrow([18.0,28.0,40.0], [19.6,30.549,39.07], color="blue red", name="Arrows_17.2999992371_5")

cluster_dict["17.2999992371"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(17.0), float(37.5), float(40.5), float(1.0)]

cluster_dict["17.2999992371_arrows"] += cgo_arrow([17.0,37.5,40.5], [16.019,37.942,37.737], color="blue red", name="Arrows_17.2999992371_6")

cluster_dict["17.2999992371"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.5), float(33.0), float(34.0), float(1.0)]

cluster_dict["17.2999992371_arrows"] += cgo_arrow([19.5,33.0,34.0], [16.762,32.791,33.752], color="blue red", name="Arrows_17.2999992371_7")

cluster_dict["17.2999992371"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(13.9507069323), float(35.4969849971), float(40.8339925715), float(1.0)]


cluster_dict["17.2999992371"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(14.023030533), float(23.6688922864), float(37.9297947025), float(1.0)]


cluster_dict["17.2999992371"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(19.0600451992), float(29.7991617379), float(32.1416507133), float(1.0)]


cluster_dict["17.2999992371"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(21.0938360018), float(28.1228606777), float(24.2721536102), float(1.0)]


cluster_dict["17.2999992371"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(26.842647015), float(35.4332821446), float(41.3000581966), float(1.0)]


cluster_dict["17.2999992371"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(27.0), float(29.5), float(1.0)]

cluster_dict["17.2999992371_arrows"] += cgo_arrow([16.0,27.0,29.5], [12.214,25.278,28.555], color="red blue", name="Arrows_17.2999992371_8")

cluster_dict["17.2999992371"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(27.5), float(30.0), float(1.0)]

cluster_dict["17.2999992371_arrows"] += cgo_arrow([18.5,27.5,30.0], [20.516,24.849,30.866], color="red blue", name="Arrows_17.2999992371_9")

cluster_dict["17.2999992371"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(31.5), float(29.0), float(1.0)]

cluster_dict["17.2999992371_arrows"] += cgo_arrow([21.0,31.5,29.0], [20.757,32.781,26.406], color="red blue", name="Arrows_17.2999992371_10")

cmd.load_cgo(cluster_dict["17.2999992371"], "Features_17.2999992371", 1)
cmd.load_cgo(cluster_dict["17.2999992371_arrows"], "Arrows_17.2999992371")
cmd.set("transparency", 0.2,"Features_17.2999992371")
cmd.group("Pharmacophore_17.2999992371", members="Features_17.2999992371")
cmd.group("Pharmacophore_17.2999992371", members="Arrows_17.2999992371")

if dirpath:
    f = join(dirpath, "28/label_threshold_17.2999992371.mol2")
else:
    f = "28/label_threshold_17.2999992371.mol2"

cmd.load(f, 'label_threshold_17.2999992371')
cmd.hide('everything', 'label_threshold_17.2999992371')
cmd.label("label_threshold_17.2999992371", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.2999992371', members= 'label_threshold_17.2999992371')


if dirpath:
    f = join(dirpath, '28/mesh.grd')
else:
    f = '28/mesh.grd'
cmd.load(f, 'mesh_28')
cmd.isomesh("isomesh_28", "mesh_28", 0.9)
cmd.color("grey80", "isomesh_28")
cmd.set('transparency', 0.4, "isomesh_28")

cmd.group('hotspot_28', "isomesh_28")
cmd.group('hotspot_28', "mesh_28")

if dirpath:
    f = join(dirpath, "29/label_threshold_15.6.mol2")
else:
    f = "29/label_threshold_15.6.mol2"

cmd.load(f, 'label_threshold_15.6')
cmd.hide('everything', 'label_threshold_15.6')
cmd.label("label_threshold_15.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.6]
gfiles = ['29/donor.grd', '29/apolar.grd', '29/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 29
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


cluster_dict = {"17.2539997101":[], "17.2539997101_arrows":[]}

cluster_dict["17.2539997101"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(17.5), float(48.0), float(1.0)]

cluster_dict["17.2539997101_arrows"] += cgo_arrow([23.0,17.5,48.0], [20.183,19.86,48.328], color="blue red", name="Arrows_17.2539997101_1")

cluster_dict["17.2539997101"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(14.5), float(48.0), float(1.0)]

cluster_dict["17.2539997101_arrows"] += cgo_arrow([24.5,14.5,48.0], [23.499,12.057,46.806], color="blue red", name="Arrows_17.2539997101_2")

cluster_dict["17.2539997101"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(27.0), float(46.0), float(1.0)]

cluster_dict["17.2539997101_arrows"] += cgo_arrow([24.5,27.0,46.0], [27.036,27.958,46.663], color="blue red", name="Arrows_17.2539997101_3")

cluster_dict["17.2539997101"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(23.5886597462), float(7.40804774466), float(41.8167044071), float(1.0)]


cluster_dict["17.2539997101"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.736926168), float(20.159460061), float(48.0382985965), float(1.0)]


cluster_dict["17.2539997101"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(19.5), float(51.0), float(1.0)]

cluster_dict["17.2539997101_arrows"] += cgo_arrow([23.5,19.5,51.0], [20.57,20.868,50.662], color="red blue", name="Arrows_17.2539997101_4")

cluster_dict["17.2539997101"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(15.0), float(46.5), float(1.0)]

cluster_dict["17.2539997101_arrows"] += cgo_arrow([25.0,15.0,46.5], [25.971,15.865,43.79], color="red blue", name="Arrows_17.2539997101_5")

cluster_dict["17.2539997101"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(19.5), float(46.0), float(1.0)]

cluster_dict["17.2539997101_arrows"] += cgo_arrow([25.0,19.5,46.0], [28.358,17.279,45.014], color="red blue", name="Arrows_17.2539997101_6")

cluster_dict["17.2539997101"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(24.0), float(48.5), float(1.0)]

cluster_dict["17.2539997101_arrows"] += cgo_arrow([26.5,24.0,48.5], [28.337,26.463,48.684], color="red blue", name="Arrows_17.2539997101_7")

cluster_dict["17.2539997101"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(16.0), float(48.5), float(1.0)]

cluster_dict["17.2539997101_arrows"] += cgo_arrow([26.5,16.0,48.5], [28.358,17.279,45.014], color="red blue", name="Arrows_17.2539997101_8")

cluster_dict["17.2539997101"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(5.5), float(47.5), float(1.0)]

cluster_dict["17.2539997101_arrows"] += cgo_arrow([27.0,5.5,47.5], [29.655,6.397,46.414], color="red blue", name="Arrows_17.2539997101_9")

cmd.load_cgo(cluster_dict["17.2539997101"], "Features_17.2539997101", 1)
cmd.load_cgo(cluster_dict["17.2539997101_arrows"], "Arrows_17.2539997101")
cmd.set("transparency", 0.2,"Features_17.2539997101")
cmd.group("Pharmacophore_17.2539997101", members="Features_17.2539997101")
cmd.group("Pharmacophore_17.2539997101", members="Arrows_17.2539997101")

if dirpath:
    f = join(dirpath, "29/label_threshold_17.2539997101.mol2")
else:
    f = "29/label_threshold_17.2539997101.mol2"

cmd.load(f, 'label_threshold_17.2539997101')
cmd.hide('everything', 'label_threshold_17.2539997101')
cmd.label("label_threshold_17.2539997101", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.2539997101', members= 'label_threshold_17.2539997101')


if dirpath:
    f = join(dirpath, '29/mesh.grd')
else:
    f = '29/mesh.grd'
cmd.load(f, 'mesh_29')
cmd.isomesh("isomesh_29", "mesh_29", 0.9)
cmd.color("grey80", "isomesh_29")
cmd.set('transparency', 0.4, "isomesh_29")

cmd.group('hotspot_29', "isomesh_29")
cmd.group('hotspot_29', "mesh_29")

if dirpath:
    f = join(dirpath, "30/label_threshold_16.0.mol2")
else:
    f = "30/label_threshold_16.0.mol2"

cmd.load(f, 'label_threshold_16.0')
cmd.hide('everything', 'label_threshold_16.0')
cmd.label("label_threshold_16.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [16.0]
gfiles = ['30/donor.grd', '30/apolar.grd', '30/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 30
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


cluster_dict = {"17.1770000458":[], "17.1770000458_arrows":[]}

cluster_dict["17.1770000458"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(16.5), float(31.5), float(27.5), float(1.0)]

cluster_dict["17.1770000458_arrows"] += cgo_arrow([16.5,31.5,27.5], [18.584,34.583,29.301], color="blue red", name="Arrows_17.1770000458_1")

cluster_dict["17.1770000458"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(31.5), float(26.0), float(1.0)]

cluster_dict["17.1770000458_arrows"] += cgo_arrow([16.0,31.5,26.0], [14.264,28.282,24.68], color="blue red", name="Arrows_17.1770000458_2")

cluster_dict["17.1770000458"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(16.5), float(31.5), float(27.5), float(1.0)]

cluster_dict["17.1770000458_arrows"] += cgo_arrow([16.5,31.5,27.5], [18.584,34.583,29.301], color="blue red", name="Arrows_17.1770000458_3")

cluster_dict["17.1770000458"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(28.0), float(40.0), float(1.0)]

cluster_dict["17.1770000458_arrows"] += cgo_arrow([18.0,28.0,40.0], [19.6,30.549,39.07], color="blue red", name="Arrows_17.1770000458_4")

cluster_dict["17.1770000458"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.5), float(33.0), float(34.0), float(1.0)]

cluster_dict["17.1770000458_arrows"] += cgo_arrow([19.5,33.0,34.0], [16.762,32.791,33.752], color="blue red", name="Arrows_17.1770000458_5")

cluster_dict["17.1770000458"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(13.5686560093), float(34.5240737596), float(40.5060185822), float(1.0)]


cluster_dict["17.1770000458"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(13.545177045), float(22.9346666888), float(38.0961426547), float(1.0)]


cluster_dict["17.1770000458"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(19.0796077114), float(29.8161042763), float(32.1359048255), float(1.0)]


cluster_dict["17.1770000458"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(20.4102658895), float(27.8250703292), float(24.1446848278), float(1.0)]


cluster_dict["17.1770000458"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(11.5), float(22.5), float(38.5), float(1.0)]

cluster_dict["17.1770000458_arrows"] += cgo_arrow([11.5,22.5,38.5], [10.91,21.842,41.086], color="red blue", name="Arrows_17.1770000458_6")

cluster_dict["17.1770000458"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(27.0), float(29.5), float(1.0)]

cluster_dict["17.1770000458_arrows"] += cgo_arrow([16.0,27.0,29.5], [12.214,25.278,28.555], color="red blue", name="Arrows_17.1770000458_7")

cluster_dict["17.1770000458"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(27.5), float(30.0), float(1.0)]

cluster_dict["17.1770000458_arrows"] += cgo_arrow([18.5,27.5,30.0], [20.516,24.849,30.866], color="red blue", name="Arrows_17.1770000458_8")

cluster_dict["17.1770000458"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(31.5), float(29.0), float(1.0)]

cluster_dict["17.1770000458_arrows"] += cgo_arrow([21.0,31.5,29.0], [20.757,32.781,26.406], color="red blue", name="Arrows_17.1770000458_9")

cluster_dict["17.1770000458"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(23.0), float(22.5), float(1.0)]

cluster_dict["17.1770000458_arrows"] += cgo_arrow([22.0,23.0,22.5], [20.807,22.261,23.987], color="red blue", name="Arrows_17.1770000458_10")

cmd.load_cgo(cluster_dict["17.1770000458"], "Features_17.1770000458", 1)
cmd.load_cgo(cluster_dict["17.1770000458_arrows"], "Arrows_17.1770000458")
cmd.set("transparency", 0.2,"Features_17.1770000458")
cmd.group("Pharmacophore_17.1770000458", members="Features_17.1770000458")
cmd.group("Pharmacophore_17.1770000458", members="Arrows_17.1770000458")

if dirpath:
    f = join(dirpath, "30/label_threshold_17.1770000458.mol2")
else:
    f = "30/label_threshold_17.1770000458.mol2"

cmd.load(f, 'label_threshold_17.1770000458')
cmd.hide('everything', 'label_threshold_17.1770000458')
cmd.label("label_threshold_17.1770000458", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.1770000458', members= 'label_threshold_17.1770000458')


if dirpath:
    f = join(dirpath, '30/mesh.grd')
else:
    f = '30/mesh.grd'
cmd.load(f, 'mesh_30')
cmd.isomesh("isomesh_30", "mesh_30", 0.9)
cmd.color("grey80", "isomesh_30")
cmd.set('transparency', 0.4, "isomesh_30")

cmd.group('hotspot_30', "isomesh_30")
cmd.group('hotspot_30', "mesh_30")

if dirpath:
    f = join(dirpath, "31/label_threshold_15.2.mol2")
else:
    f = "31/label_threshold_15.2.mol2"

cmd.load(f, 'label_threshold_15.2')
cmd.hide('everything', 'label_threshold_15.2')
cmd.label("label_threshold_15.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.2]
gfiles = ['31/donor.grd', '31/apolar.grd', '31/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 31
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


cluster_dict = {"17.1639995575":[], "17.1639995575_arrows":[]}

cluster_dict["17.1639995575"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(26.0), float(34.0), float(1.0)]

cluster_dict["17.1639995575_arrows"] += cgo_arrow([16.0,26.0,34.0], [17.727,23.511,33.873], color="blue red", name="Arrows_17.1639995575_1")

cluster_dict["17.1639995575"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(16.5), float(31.5), float(27.5), float(1.0)]

cluster_dict["17.1639995575_arrows"] += cgo_arrow([16.5,31.5,27.5], [18.584,34.583,29.301], color="blue red", name="Arrows_17.1639995575_2")

cluster_dict["17.1639995575"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(31.5), float(26.0), float(1.0)]

cluster_dict["17.1639995575_arrows"] += cgo_arrow([16.0,31.5,26.0], [14.264,28.282,24.68], color="blue red", name="Arrows_17.1639995575_3")

cluster_dict["17.1639995575"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(16.5), float(31.5), float(27.5), float(1.0)]

cluster_dict["17.1639995575_arrows"] += cgo_arrow([16.5,31.5,27.5], [18.584,34.583,29.301], color="blue red", name="Arrows_17.1639995575_4")

cluster_dict["17.1639995575"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(20.0), float(29.5), float(34.5), float(1.0)]

cluster_dict["17.1639995575_arrows"] += cgo_arrow([20.0,29.5,34.5], [21.495,28.702,32.346], color="blue red", name="Arrows_17.1639995575_5")

cluster_dict["17.1639995575"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.5), float(27.0), float(30.5), float(1.0)]

cluster_dict["17.1639995575_arrows"] += cgo_arrow([19.5,27.0,30.5], [20.516,24.849,30.866], color="blue red", name="Arrows_17.1639995575_6")

cluster_dict["17.1639995575"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.5), float(33.0), float(34.0), float(1.0)]

cluster_dict["17.1639995575_arrows"] += cgo_arrow([19.5,33.0,34.0], [16.762,32.791,33.752], color="blue red", name="Arrows_17.1639995575_7")

cluster_dict["17.1639995575"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(19.1574828273), float(29.8387693495), float(31.853603929), float(1.0)]


cluster_dict["17.1639995575"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(19.9076767473), float(27.4253208617), float(24.3383195277), float(1.0)]


cluster_dict["17.1639995575"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(27.0), float(29.5), float(1.0)]

cluster_dict["17.1639995575_arrows"] += cgo_arrow([16.0,27.0,29.5], [12.214,25.278,28.555], color="red blue", name="Arrows_17.1639995575_8")

cluster_dict["17.1639995575"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(31.5), float(29.0), float(1.0)]

cluster_dict["17.1639995575_arrows"] += cgo_arrow([21.0,31.5,29.0], [20.757,32.781,26.406], color="red blue", name="Arrows_17.1639995575_9")

cluster_dict["17.1639995575"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(23.0), float(22.5), float(1.0)]

cluster_dict["17.1639995575_arrows"] += cgo_arrow([22.0,23.0,22.5], [20.807,22.261,23.987], color="red blue", name="Arrows_17.1639995575_10")

cmd.load_cgo(cluster_dict["17.1639995575"], "Features_17.1639995575", 1)
cmd.load_cgo(cluster_dict["17.1639995575_arrows"], "Arrows_17.1639995575")
cmd.set("transparency", 0.2,"Features_17.1639995575")
cmd.group("Pharmacophore_17.1639995575", members="Features_17.1639995575")
cmd.group("Pharmacophore_17.1639995575", members="Arrows_17.1639995575")

if dirpath:
    f = join(dirpath, "31/label_threshold_17.1639995575.mol2")
else:
    f = "31/label_threshold_17.1639995575.mol2"

cmd.load(f, 'label_threshold_17.1639995575')
cmd.hide('everything', 'label_threshold_17.1639995575')
cmd.label("label_threshold_17.1639995575", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.1639995575', members= 'label_threshold_17.1639995575')


if dirpath:
    f = join(dirpath, '31/mesh.grd')
else:
    f = '31/mesh.grd'
cmd.load(f, 'mesh_31')
cmd.isomesh("isomesh_31", "mesh_31", 0.9)
cmd.color("grey80", "isomesh_31")
cmd.set('transparency', 0.4, "isomesh_31")

cmd.group('hotspot_31', "isomesh_31")
cmd.group('hotspot_31', "mesh_31")

if dirpath:
    f = join(dirpath, "32/label_threshold_14.6.mol2")
else:
    f = "32/label_threshold_14.6.mol2"

cmd.load(f, 'label_threshold_14.6')
cmd.hide('everything', 'label_threshold_14.6')
cmd.label("label_threshold_14.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [14.6]
gfiles = ['32/donor.grd', '32/apolar.grd', '32/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 32
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


cluster_dict = {"17.0900001526":[], "17.0900001526_arrows":[]}

cluster_dict["17.0900001526"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(17.5), float(48.0), float(1.0)]

cluster_dict["17.0900001526_arrows"] += cgo_arrow([23.0,17.5,48.0], [20.183,19.86,48.328], color="blue red", name="Arrows_17.0900001526_1")

cluster_dict["17.0900001526"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(14.5), float(48.0), float(1.0)]

cluster_dict["17.0900001526_arrows"] += cgo_arrow([24.5,14.5,48.0], [23.499,12.057,46.806], color="blue red", name="Arrows_17.0900001526_2")

cluster_dict["17.0900001526"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(17.5), float(46.0), float(1.0)]

cluster_dict["17.0900001526_arrows"] += cgo_arrow([23.5,17.5,46.0], [25.971,15.865,43.79], color="blue red", name="Arrows_17.0900001526_3")

cluster_dict["17.0900001526"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(27.0), float(46.0), float(1.0)]

cluster_dict["17.0900001526_arrows"] += cgo_arrow([24.5,27.0,46.0], [27.036,27.958,46.663], color="blue red", name="Arrows_17.0900001526_4")

cluster_dict["17.0900001526"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.798750439), float(20.1430864569), float(48.0809410831), float(1.0)]


cluster_dict["17.0900001526"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(28.5), float(47.0), float(1.0)]


cluster_dict["17.0900001526"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(19.5), float(51.0), float(1.0)]

cluster_dict["17.0900001526_arrows"] += cgo_arrow([23.5,19.5,51.0], [20.57,20.868,50.662], color="red blue", name="Arrows_17.0900001526_5")

cluster_dict["17.0900001526"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(15.0), float(46.5), float(1.0)]

cluster_dict["17.0900001526_arrows"] += cgo_arrow([25.0,15.0,46.5], [25.971,15.865,43.79], color="red blue", name="Arrows_17.0900001526_6")

cluster_dict["17.0900001526"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(19.5), float(46.0), float(1.0)]

cluster_dict["17.0900001526_arrows"] += cgo_arrow([25.0,19.5,46.0], [28.358,17.279,45.014], color="red blue", name="Arrows_17.0900001526_7")

cluster_dict["17.0900001526"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(24.0), float(48.5), float(1.0)]

cluster_dict["17.0900001526_arrows"] += cgo_arrow([26.5,24.0,48.5], [28.337,26.463,48.684], color="red blue", name="Arrows_17.0900001526_8")

cluster_dict["17.0900001526"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(16.0), float(48.5), float(1.0)]

cluster_dict["17.0900001526_arrows"] += cgo_arrow([26.5,16.0,48.5], [28.358,17.279,45.014], color="red blue", name="Arrows_17.0900001526_9")

cmd.load_cgo(cluster_dict["17.0900001526"], "Features_17.0900001526", 1)
cmd.load_cgo(cluster_dict["17.0900001526_arrows"], "Arrows_17.0900001526")
cmd.set("transparency", 0.2,"Features_17.0900001526")
cmd.group("Pharmacophore_17.0900001526", members="Features_17.0900001526")
cmd.group("Pharmacophore_17.0900001526", members="Arrows_17.0900001526")

if dirpath:
    f = join(dirpath, "32/label_threshold_17.0900001526.mol2")
else:
    f = "32/label_threshold_17.0900001526.mol2"

cmd.load(f, 'label_threshold_17.0900001526')
cmd.hide('everything', 'label_threshold_17.0900001526')
cmd.label("label_threshold_17.0900001526", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.0900001526', members= 'label_threshold_17.0900001526')


if dirpath:
    f = join(dirpath, '32/mesh.grd')
else:
    f = '32/mesh.grd'
cmd.load(f, 'mesh_32')
cmd.isomesh("isomesh_32", "mesh_32", 0.9)
cmd.color("grey80", "isomesh_32")
cmd.set('transparency', 0.4, "isomesh_32")

cmd.group('hotspot_32', "isomesh_32")
cmd.group('hotspot_32', "mesh_32")

if dirpath:
    f = join(dirpath, "33/label_threshold_16.4.mol2")
else:
    f = "33/label_threshold_16.4.mol2"

cmd.load(f, 'label_threshold_16.4')
cmd.hide('everything', 'label_threshold_16.4')
cmd.label("label_threshold_16.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [16.4]
gfiles = ['33/donor.grd', '33/apolar.grd', '33/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 33
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


cluster_dict = {"17.0769996643":[], "17.0769996643_arrows":[]}

cluster_dict["17.0769996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(0.0), float(68.0), float(1.0)]

cluster_dict["17.0769996643_arrows"] += cgo_arrow([-37.0,0.0,68.0], [-38.407,0.945,65.456], color="blue red", name="Arrows_17.0769996643_1")

cluster_dict["17.0769996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(-2.5), float(73.5), float(1.0)]

cluster_dict["17.0769996643_arrows"] += cgo_arrow([-37.0,-2.5,73.5], [-35.226,-0.918,75.907], color="blue red", name="Arrows_17.0769996643_2")

cluster_dict["17.0769996643"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-43.4815195696), float(9.20018945833), float(80.3474430582), float(1.0)]


cluster_dict["17.0769996643"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-39.3397609155), float(-2.07954776609), float(72.8075442059), float(1.0)]


cluster_dict["17.0769996643"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-32.127757363), float(4.21270204321), float(74.5956944949), float(1.0)]


cluster_dict["17.0769996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-43.5), float(9.5), float(78.5), float(1.0)]

cluster_dict["17.0769996643_arrows"] += cgo_arrow([-43.5,9.5,78.5], [-44.539,7.475,75.686], color="red blue", name="Arrows_17.0769996643_3")

cluster_dict["17.0769996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-41.5), float(11.0), float(68.0), float(1.0)]

cluster_dict["17.0769996643_arrows"] += cgo_arrow([-41.5,11.0,68.0], [-42.298,9.428,66.293], color="red blue", name="Arrows_17.0769996643_4")

cluster_dict["17.0769996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-38.5), float(-6.5), float(74.0), float(1.0)]

cluster_dict["17.0769996643_arrows"] += cgo_arrow([-38.5,-6.5,74.0], [-36.711,-4.673,75.752], color="red blue", name="Arrows_17.0769996643_5")

cluster_dict["17.0769996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-35.5), float(4.0), float(70.5), float(1.0)]

cluster_dict["17.0769996643_arrows"] += cgo_arrow([-35.5,4.0,70.5], [-33.782,2.549,68.962], color="red blue", name="Arrows_17.0769996643_6")

cluster_dict["17.0769996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-35.0), float(4.5), float(67.0), float(1.0)]

cluster_dict["17.0769996643_arrows"] += cgo_arrow([-35.0,4.5,67.0], [-33.782,2.549,68.962], color="red blue", name="Arrows_17.0769996643_7")

cluster_dict["17.0769996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-36.0), float(6.0), float(65.0), float(1.0)]

cluster_dict["17.0769996643_arrows"] += cgo_arrow([-36.0,6.0,65.0], [-37.903,8.797,63.64], color="red blue", name="Arrows_17.0769996643_8")

cluster_dict["17.0769996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-32.5), float(6.0), float(70.5), float(1.0)]

cluster_dict["17.0769996643_arrows"] += cgo_arrow([-32.5,6.0,70.5], [-33.782,2.549,68.962], color="red blue", name="Arrows_17.0769996643_9")

cmd.load_cgo(cluster_dict["17.0769996643"], "Features_17.0769996643", 1)
cmd.load_cgo(cluster_dict["17.0769996643_arrows"], "Arrows_17.0769996643")
cmd.set("transparency", 0.2,"Features_17.0769996643")
cmd.group("Pharmacophore_17.0769996643", members="Features_17.0769996643")
cmd.group("Pharmacophore_17.0769996643", members="Arrows_17.0769996643")

if dirpath:
    f = join(dirpath, "33/label_threshold_17.0769996643.mol2")
else:
    f = "33/label_threshold_17.0769996643.mol2"

cmd.load(f, 'label_threshold_17.0769996643')
cmd.hide('everything', 'label_threshold_17.0769996643')
cmd.label("label_threshold_17.0769996643", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.0769996643', members= 'label_threshold_17.0769996643')


if dirpath:
    f = join(dirpath, '33/mesh.grd')
else:
    f = '33/mesh.grd'
cmd.load(f, 'mesh_33')
cmd.isomesh("isomesh_33", "mesh_33", 0.9)
cmd.color("grey80", "isomesh_33")
cmd.set('transparency', 0.4, "isomesh_33")

cmd.group('hotspot_33', "isomesh_33")
cmd.group('hotspot_33', "mesh_33")

if dirpath:
    f = join(dirpath, "34/label_threshold_16.2.mol2")
else:
    f = "34/label_threshold_16.2.mol2"

cmd.load(f, 'label_threshold_16.2')
cmd.hide('everything', 'label_threshold_16.2')
cmd.label("label_threshold_16.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [16.2]
gfiles = ['34/donor.grd', '34/apolar.grd', '34/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 34
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


cluster_dict = {"17.0659999847":[], "17.0659999847_arrows":[]}

cluster_dict["17.0659999847"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(0.0), float(68.0), float(1.0)]

cluster_dict["17.0659999847_arrows"] += cgo_arrow([-37.0,0.0,68.0], [-38.407,0.945,65.456], color="blue red", name="Arrows_17.0659999847_1")

cluster_dict["17.0659999847"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(-2.5), float(73.5), float(1.0)]

cluster_dict["17.0659999847_arrows"] += cgo_arrow([-37.0,-2.5,73.5], [-35.226,-0.918,75.907], color="blue red", name="Arrows_17.0659999847_2")

cluster_dict["17.0659999847"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-41.1913122579), float(-9.44466663169), float(67.9655497924), float(1.0)]


cluster_dict["17.0659999847"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-43.9557996207), float(7.29268507496), float(78.6890920318), float(1.0)]


cluster_dict["17.0659999847"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-39.2510902394), float(-1.98905385467), float(72.5806423876), float(1.0)]


cluster_dict["17.0659999847"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-33.798312702), float(4.05451260087), float(74.9820527454), float(1.0)]


cluster_dict["17.0659999847"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-40.5), float(-3.0), float(71.5), float(1.0)]

cluster_dict["17.0659999847_arrows"] += cgo_arrow([-40.5,-3.0,71.5], [-41.252,-0.556,70.849], color="red blue", name="Arrows_17.0659999847_3")

cluster_dict["17.0659999847"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-38.5), float(-6.5), float(74.0), float(1.0)]

cluster_dict["17.0659999847_arrows"] += cgo_arrow([-38.5,-6.5,74.0], [-36.711,-4.673,75.752], color="red blue", name="Arrows_17.0659999847_4")

cluster_dict["17.0659999847"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-38.0), float(-5.5), float(63.0), float(1.0)]

cluster_dict["17.0659999847_arrows"] += cgo_arrow([-38.0,-5.5,63.0], [-35.514,-4.848,62.445], color="red blue", name="Arrows_17.0659999847_5")

cluster_dict["17.0659999847"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-35.5), float(4.0), float(70.5), float(1.0)]

cluster_dict["17.0659999847_arrows"] += cgo_arrow([-35.5,4.0,70.5], [-33.782,2.549,68.962], color="red blue", name="Arrows_17.0659999847_6")

cluster_dict["17.0659999847"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-35.0), float(4.5), float(67.0), float(1.0)]

cluster_dict["17.0659999847_arrows"] += cgo_arrow([-35.0,4.5,67.0], [-33.782,2.549,68.962], color="red blue", name="Arrows_17.0659999847_7")

cluster_dict["17.0659999847"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-36.0), float(5.5), float(65.0), float(1.0)]

cluster_dict["17.0659999847_arrows"] += cgo_arrow([-36.0,5.5,65.0], [-37.903,8.797,63.64], color="red blue", name="Arrows_17.0659999847_8")

cluster_dict["17.0659999847"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-32.5), float(6.0), float(70.5), float(1.0)]

cluster_dict["17.0659999847_arrows"] += cgo_arrow([-32.5,6.0,70.5], [-33.782,2.549,68.962], color="red blue", name="Arrows_17.0659999847_9")

cmd.load_cgo(cluster_dict["17.0659999847"], "Features_17.0659999847", 1)
cmd.load_cgo(cluster_dict["17.0659999847_arrows"], "Arrows_17.0659999847")
cmd.set("transparency", 0.2,"Features_17.0659999847")
cmd.group("Pharmacophore_17.0659999847", members="Features_17.0659999847")
cmd.group("Pharmacophore_17.0659999847", members="Arrows_17.0659999847")

if dirpath:
    f = join(dirpath, "34/label_threshold_17.0659999847.mol2")
else:
    f = "34/label_threshold_17.0659999847.mol2"

cmd.load(f, 'label_threshold_17.0659999847')
cmd.hide('everything', 'label_threshold_17.0659999847')
cmd.label("label_threshold_17.0659999847", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.0659999847', members= 'label_threshold_17.0659999847')


if dirpath:
    f = join(dirpath, '34/mesh.grd')
else:
    f = '34/mesh.grd'
cmd.load(f, 'mesh_34')
cmd.isomesh("isomesh_34", "mesh_34", 0.9)
cmd.color("grey80", "isomesh_34")
cmd.set('transparency', 0.4, "isomesh_34")

cmd.group('hotspot_34', "isomesh_34")
cmd.group('hotspot_34', "mesh_34")

if dirpath:
    f = join(dirpath, "35/label_threshold_16.1.mol2")
else:
    f = "35/label_threshold_16.1.mol2"

cmd.load(f, 'label_threshold_16.1')
cmd.hide('everything', 'label_threshold_16.1')
cmd.label("label_threshold_16.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [16.1]
gfiles = ['35/donor.grd', '35/apolar.grd', '35/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 35
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


cluster_dict = {"17.0480003357":[], "17.0480003357_arrows":[]}

cluster_dict["17.0480003357"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(0.0), float(68.0), float(1.0)]

cluster_dict["17.0480003357_arrows"] += cgo_arrow([-37.0,0.0,68.0], [-38.407,0.945,65.456], color="blue red", name="Arrows_17.0480003357_1")

cluster_dict["17.0480003357"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(-2.5), float(73.5), float(1.0)]

cluster_dict["17.0480003357_arrows"] += cgo_arrow([-37.0,-2.5,73.5], [-35.226,-0.918,75.907], color="blue red", name="Arrows_17.0480003357_2")

cluster_dict["17.0480003357"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-41.1466778425), float(-9.44558431507), float(67.9329518783), float(1.0)]


cluster_dict["17.0480003357"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-39.2393125754), float(-2.00428445017), float(72.5662398405), float(1.0)]


cluster_dict["17.0480003357"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-33.3473291798), float(3.86170797046), float(74.5368041915), float(1.0)]


cluster_dict["17.0480003357"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-40.5), float(-3.0), float(71.5), float(1.0)]

cluster_dict["17.0480003357_arrows"] += cgo_arrow([-40.5,-3.0,71.5], [-41.252,-0.556,70.849], color="red blue", name="Arrows_17.0480003357_3")

cluster_dict["17.0480003357"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-38.5), float(-6.5), float(74.0), float(1.0)]

cluster_dict["17.0480003357_arrows"] += cgo_arrow([-38.5,-6.5,74.0], [-36.711,-4.673,75.752], color="red blue", name="Arrows_17.0480003357_4")

cluster_dict["17.0480003357"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-38.0), float(-5.5), float(63.0), float(1.0)]

cluster_dict["17.0480003357_arrows"] += cgo_arrow([-38.0,-5.5,63.0], [-35.514,-4.848,62.445], color="red blue", name="Arrows_17.0480003357_5")

cluster_dict["17.0480003357"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-35.5), float(4.0), float(70.5), float(1.0)]

cluster_dict["17.0480003357_arrows"] += cgo_arrow([-35.5,4.0,70.5], [-33.782,2.549,68.962], color="red blue", name="Arrows_17.0480003357_6")

cluster_dict["17.0480003357"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-35.0), float(4.5), float(67.0), float(1.0)]

cluster_dict["17.0480003357_arrows"] += cgo_arrow([-35.0,4.5,67.0], [-33.782,2.549,68.962], color="red blue", name="Arrows_17.0480003357_7")

cluster_dict["17.0480003357"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-36.0), float(6.0), float(65.0), float(1.0)]

cluster_dict["17.0480003357_arrows"] += cgo_arrow([-36.0,6.0,65.0], [-37.903,8.797,63.64], color="red blue", name="Arrows_17.0480003357_8")

cluster_dict["17.0480003357"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-32.5), float(6.0), float(70.5), float(1.0)]

cluster_dict["17.0480003357_arrows"] += cgo_arrow([-32.5,6.0,70.5], [-33.782,2.549,68.962], color="red blue", name="Arrows_17.0480003357_9")

cmd.load_cgo(cluster_dict["17.0480003357"], "Features_17.0480003357", 1)
cmd.load_cgo(cluster_dict["17.0480003357_arrows"], "Arrows_17.0480003357")
cmd.set("transparency", 0.2,"Features_17.0480003357")
cmd.group("Pharmacophore_17.0480003357", members="Features_17.0480003357")
cmd.group("Pharmacophore_17.0480003357", members="Arrows_17.0480003357")

if dirpath:
    f = join(dirpath, "35/label_threshold_17.0480003357.mol2")
else:
    f = "35/label_threshold_17.0480003357.mol2"

cmd.load(f, 'label_threshold_17.0480003357')
cmd.hide('everything', 'label_threshold_17.0480003357')
cmd.label("label_threshold_17.0480003357", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.0480003357', members= 'label_threshold_17.0480003357')


if dirpath:
    f = join(dirpath, '35/mesh.grd')
else:
    f = '35/mesh.grd'
cmd.load(f, 'mesh_35')
cmd.isomesh("isomesh_35", "mesh_35", 0.9)
cmd.color("grey80", "isomesh_35")
cmd.set('transparency', 0.4, "isomesh_35")

cmd.group('hotspot_35', "isomesh_35")
cmd.group('hotspot_35', "mesh_35")

if dirpath:
    f = join(dirpath, "36/label_threshold_15.5.mol2")
else:
    f = "36/label_threshold_15.5.mol2"

cmd.load(f, 'label_threshold_15.5')
cmd.hide('everything', 'label_threshold_15.5')
cmd.label("label_threshold_15.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.5]
gfiles = ['36/donor.grd', '36/apolar.grd', '36/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 36
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


cluster_dict = {"16.9969997406":[], "16.9969997406_arrows":[]}

cluster_dict["16.9969997406"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(9.0), float(6.5), float(43.0), float(1.0)]

cluster_dict["16.9969997406_arrows"] += cgo_arrow([9.0,6.5,43.0], [10.755,4.658,42.464], color="blue red", name="Arrows_16.9969997406_1")

cluster_dict["16.9969997406"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(10.0), float(7.0), float(44.5), float(1.0)]

cluster_dict["16.9969997406_arrows"] += cgo_arrow([10.0,7.0,44.5], [9.073,4.46,45.522], color="blue red", name="Arrows_16.9969997406_2")

cluster_dict["16.9969997406"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(20.5), float(7.5), float(36.5), float(1.0)]

cluster_dict["16.9969997406_arrows"] += cgo_arrow([20.5,7.5,36.5], [19.802,5.432,35.84], color="blue red", name="Arrows_16.9969997406_3")

cluster_dict["16.9969997406"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(9.0), float(41.5), float(1.0)]

cluster_dict["16.9969997406_arrows"] += cgo_arrow([23.5,9.0,41.5], [24.082,12.057,42.288], color="blue red", name="Arrows_16.9969997406_4")

cluster_dict["16.9969997406"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(3.54361820039), float(10.5189474952), float(47.8947576848), float(1.0)]


cluster_dict["16.9969997406"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(8.25849106806), float(8.56164233022), float(44.0539478586), float(1.0)]


cluster_dict["16.9969997406"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(21.0818277686), float(5.28850115063), float(41.0475434735), float(1.0)]


cluster_dict["16.9969997406"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(1.5), float(8.5), float(48.0), float(1.0)]

cluster_dict["16.9969997406_arrows"] += cgo_arrow([1.5,8.5,48.0], [1.122,5.617,48.391], color="red blue", name="Arrows_16.9969997406_5")

cluster_dict["16.9969997406"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(5.5), float(8.0), float(43.0), float(1.0)]

cluster_dict["16.9969997406_arrows"] += cgo_arrow([5.5,8.0,43.0], [3.715,8.369,43.968], color="red blue", name="Arrows_16.9969997406_6")

cluster_dict["16.9969997406"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(10.0), float(6.5), float(50.5), float(1.0)]

cluster_dict["16.9969997406_arrows"] += cgo_arrow([10.0,6.5,50.5], [7.899,6.991,52.135], color="red blue", name="Arrows_16.9969997406_7")

cluster_dict["16.9969997406"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(13.5), float(13.0), float(45.0), float(1.0)]

cluster_dict["16.9969997406_arrows"] += cgo_arrow([13.5,13.0,45.0], [14.379,15.215,43.186], color="red blue", name="Arrows_16.9969997406_8")

cluster_dict["16.9969997406"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(14.5), float(-2.5), float(48.0), float(1.0)]

cluster_dict["16.9969997406_arrows"] += cgo_arrow([14.5,-2.5,48.0], [15.22,-0.595,51.059], color="red blue", name="Arrows_16.9969997406_9")

cmd.load_cgo(cluster_dict["16.9969997406"], "Features_16.9969997406", 1)
cmd.load_cgo(cluster_dict["16.9969997406_arrows"], "Arrows_16.9969997406")
cmd.set("transparency", 0.2,"Features_16.9969997406")
cmd.group("Pharmacophore_16.9969997406", members="Features_16.9969997406")
cmd.group("Pharmacophore_16.9969997406", members="Arrows_16.9969997406")

if dirpath:
    f = join(dirpath, "36/label_threshold_16.9969997406.mol2")
else:
    f = "36/label_threshold_16.9969997406.mol2"

cmd.load(f, 'label_threshold_16.9969997406')
cmd.hide('everything', 'label_threshold_16.9969997406')
cmd.label("label_threshold_16.9969997406", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.9969997406', members= 'label_threshold_16.9969997406')


if dirpath:
    f = join(dirpath, '36/mesh.grd')
else:
    f = '36/mesh.grd'
cmd.load(f, 'mesh_36')
cmd.isomesh("isomesh_36", "mesh_36", 0.9)
cmd.color("grey80", "isomesh_36")
cmd.set('transparency', 0.4, "isomesh_36")

cmd.group('hotspot_36', "isomesh_36")
cmd.group('hotspot_36', "mesh_36")

if dirpath:
    f = join(dirpath, "37/label_threshold_5.8.mol2")
else:
    f = "37/label_threshold_5.8.mol2"

cmd.load(f, 'label_threshold_5.8')
cmd.hide('everything', 'label_threshold_5.8')
cmd.label("label_threshold_5.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [5.8]
gfiles = ['37/donor.grd', '37/apolar.grd', '37/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 37
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


cluster_dict = {"16.8869991302":[], "16.8869991302_arrows":[]}

cluster_dict["16.8869991302"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-45.3418011853), float(10.6605967389), float(80.8328514563), float(1.0)]


cluster_dict["16.8869991302"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-40.0), float(6.0), float(80.9637209972), float(1.0)]


cluster_dict["16.8869991302"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-45.0), float(13.5), float(82.5), float(1.0)]

cluster_dict["16.8869991302_arrows"] += cgo_arrow([-45.0,13.5,82.5], [-45.229,15.633,82.483], color="red blue", name="Arrows_16.8869991302_1")

cluster_dict["16.8869991302"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-43.5), float(9.5), float(78.5), float(1.0)]

cluster_dict["16.8869991302_arrows"] += cgo_arrow([-43.5,9.5,78.5], [-44.539,7.475,75.686], color="red blue", name="Arrows_16.8869991302_2")

cmd.load_cgo(cluster_dict["16.8869991302"], "Features_16.8869991302", 1)
cmd.load_cgo(cluster_dict["16.8869991302_arrows"], "Arrows_16.8869991302")
cmd.set("transparency", 0.2,"Features_16.8869991302")
cmd.group("Pharmacophore_16.8869991302", members="Features_16.8869991302")
cmd.group("Pharmacophore_16.8869991302", members="Arrows_16.8869991302")

if dirpath:
    f = join(dirpath, "37/label_threshold_16.8869991302.mol2")
else:
    f = "37/label_threshold_16.8869991302.mol2"

cmd.load(f, 'label_threshold_16.8869991302')
cmd.hide('everything', 'label_threshold_16.8869991302')
cmd.label("label_threshold_16.8869991302", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.8869991302', members= 'label_threshold_16.8869991302')


if dirpath:
    f = join(dirpath, '37/mesh.grd')
else:
    f = '37/mesh.grd'
cmd.load(f, 'mesh_37')
cmd.isomesh("isomesh_37", "mesh_37", 0.9)
cmd.color("grey80", "isomesh_37")
cmd.set('transparency', 0.4, "isomesh_37")

cmd.group('hotspot_37', "isomesh_37")
cmd.group('hotspot_37', "mesh_37")

if dirpath:
    f = join(dirpath, "38/label_threshold_15.7.mol2")
else:
    f = "38/label_threshold_15.7.mol2"

cmd.load(f, 'label_threshold_15.7')
cmd.hide('everything', 'label_threshold_15.7')
cmd.label("label_threshold_15.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.7]
gfiles = ['38/donor.grd', '38/apolar.grd', '38/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 38
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


cluster_dict = {"16.7999992371":[], "16.7999992371_arrows":[]}

cluster_dict["16.7999992371"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-56.0), float(-8.5), float(49.5), float(1.0)]

cluster_dict["16.7999992371_arrows"] += cgo_arrow([-56.0,-8.5,49.5], [-58.938,-9.438,49.082], color="blue red", name="Arrows_16.7999992371_1")

cluster_dict["16.7999992371"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-55.0), float(-12.5), float(51.5), float(1.0)]

cluster_dict["16.7999992371_arrows"] += cgo_arrow([-55.0,-12.5,51.5], [-52.315,-11.683,50.879], color="blue red", name="Arrows_16.7999992371_2")

cluster_dict["16.7999992371"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-55.0), float(-3.5), float(49.0), float(1.0)]

cluster_dict["16.7999992371_arrows"] += cgo_arrow([-55.0,-3.5,49.0], [-56.059,-0.837,49.615], color="blue red", name="Arrows_16.7999992371_3")

cluster_dict["16.7999992371"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-53.5), float(3.0), float(42.5), float(1.0)]

cluster_dict["16.7999992371_arrows"] += cgo_arrow([-53.5,3.0,42.5], [-52.061,3.063,39.683], color="blue red", name="Arrows_16.7999992371_4")

cluster_dict["16.7999992371"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-52.5), float(-6.5), float(48.0), float(1.0)]

cluster_dict["16.7999992371_arrows"] += cgo_arrow([-52.5,-6.5,48.0], [-50.9,-8.682,47.613], color="blue red", name="Arrows_16.7999992371_5")

cluster_dict["16.7999992371"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-48.0), float(-4.0), float(48.5), float(1.0)]

cluster_dict["16.7999992371_arrows"] += cgo_arrow([-48.0,-4.0,48.5], [-48.278,-5.299,45.404], color="blue red", name="Arrows_16.7999992371_6")

cluster_dict["16.7999992371"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-45.5), float(-3.0), float(50.5), float(1.0)]

cluster_dict["16.7999992371_arrows"] += cgo_arrow([-45.5,-3.0,50.5], [-42.563,-1.737,50.027], color="blue red", name="Arrows_16.7999992371_7")

cluster_dict["16.7999992371"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-57.2563957242), float(-10.2919670621), float(53.5139227967), float(1.0)]


cluster_dict["16.7999992371"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-53.1055301446), float(-5.80910728451), float(50.5315252189), float(1.0)]


cluster_dict["16.7999992371"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-56.5), float(-4.5), float(50.5), float(1.0)]

cluster_dict["16.7999992371_arrows"] += cgo_arrow([-56.5,-4.5,50.5], [-59.288,-5.348,50.967], color="red blue", name="Arrows_16.7999992371_8")

cluster_dict["16.7999992371"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-53.0), float(-7.5), float(49.0), float(1.0)]

cluster_dict["16.7999992371_arrows"] += cgo_arrow([-53.0,-7.5,49.0], [-50.9,-8.682,47.613], color="red blue", name="Arrows_16.7999992371_9")

cmd.load_cgo(cluster_dict["16.7999992371"], "Features_16.7999992371", 1)
cmd.load_cgo(cluster_dict["16.7999992371_arrows"], "Arrows_16.7999992371")
cmd.set("transparency", 0.2,"Features_16.7999992371")
cmd.group("Pharmacophore_16.7999992371", members="Features_16.7999992371")
cmd.group("Pharmacophore_16.7999992371", members="Arrows_16.7999992371")

if dirpath:
    f = join(dirpath, "38/label_threshold_16.7999992371.mol2")
else:
    f = "38/label_threshold_16.7999992371.mol2"

cmd.load(f, 'label_threshold_16.7999992371')
cmd.hide('everything', 'label_threshold_16.7999992371')
cmd.label("label_threshold_16.7999992371", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.7999992371', members= 'label_threshold_16.7999992371')


if dirpath:
    f = join(dirpath, '38/mesh.grd')
else:
    f = '38/mesh.grd'
cmd.load(f, 'mesh_38')
cmd.isomesh("isomesh_38", "mesh_38", 0.9)
cmd.color("grey80", "isomesh_38")
cmd.set('transparency', 0.4, "isomesh_38")

cmd.group('hotspot_38', "isomesh_38")
cmd.group('hotspot_38', "mesh_38")

if dirpath:
    f = join(dirpath, "39/label_threshold_15.7.mol2")
else:
    f = "39/label_threshold_15.7.mol2"

cmd.load(f, 'label_threshold_15.7')
cmd.hide('everything', 'label_threshold_15.7')
cmd.label("label_threshold_15.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.7]
gfiles = ['39/donor.grd', '39/apolar.grd', '39/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 39
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


cluster_dict = {"16.7269992828":[], "16.7269992828_arrows":[]}

cluster_dict["16.7269992828"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(20.5), float(7.5), float(36.5), float(1.0)]

cluster_dict["16.7269992828_arrows"] += cgo_arrow([20.5,7.5,36.5], [19.802,5.432,35.84], color="blue red", name="Arrows_16.7269992828_1")

cluster_dict["16.7269992828"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(8.5), float(26.5), float(1.0)]

cluster_dict["16.7269992828_arrows"] += cgo_arrow([21.0,8.5,26.5], [20.263,11.401,25.376], color="blue red", name="Arrows_16.7269992828_2")

cluster_dict["16.7269992828"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(9.0), float(41.0), float(1.0)]

cluster_dict["16.7269992828_arrows"] += cgo_arrow([23.5,9.0,41.0], [24.082,12.057,42.288], color="blue red", name="Arrows_16.7269992828_3")

cluster_dict["16.7269992828"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(6.5), float(22.5), float(1.0)]

cluster_dict["16.7269992828_arrows"] += cgo_arrow([26.5,6.5,22.5], [24.383,5.1,20.493], color="blue red", name="Arrows_16.7269992828_4")

cluster_dict["16.7269992828"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(10.5), float(24.0), float(1.0)]

cluster_dict["16.7269992828_arrows"] += cgo_arrow([26.0,10.5,24.0], [24.828,13.48,24.46], color="blue red", name="Arrows_16.7269992828_5")

cluster_dict["16.7269992828"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.0), float(4.0), float(24.5), float(1.0)]

cluster_dict["16.7269992828_arrows"] += cgo_arrow([29.0,4.0,24.5], [31.281,4.479,26.704], color="blue red", name="Arrows_16.7269992828_6")

cluster_dict["16.7269992828"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.1160899266), float(5.98590208424), float(40.761845009), float(1.0)]


cluster_dict["16.7269992828"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(26.1291690455), float(7.97590920165), float(24.0448020184), float(1.0)]


cluster_dict["16.7269992828"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(7.0), float(38.0), float(1.0)]

cluster_dict["16.7269992828_arrows"] += cgo_arrow([22.5,7.0,38.0], [22.648,5.212,35.929], color="red blue", name="Arrows_16.7269992828_7")

cluster_dict["16.7269992828"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(8.0), float(31.5), float(1.0)]

cluster_dict["16.7269992828_arrows"] += cgo_arrow([23.0,8.0,31.5], [24.327,5.898,29.592], color="red blue", name="Arrows_16.7269992828_8")

cluster_dict["16.7269992828"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(5.5), float(21.5), float(1.0)]

cluster_dict["16.7269992828_arrows"] += cgo_arrow([26.5,5.5,21.5], [24.383,5.1,20.493], color="red blue", name="Arrows_16.7269992828_9")

cmd.load_cgo(cluster_dict["16.7269992828"], "Features_16.7269992828", 1)
cmd.load_cgo(cluster_dict["16.7269992828_arrows"], "Arrows_16.7269992828")
cmd.set("transparency", 0.2,"Features_16.7269992828")
cmd.group("Pharmacophore_16.7269992828", members="Features_16.7269992828")
cmd.group("Pharmacophore_16.7269992828", members="Arrows_16.7269992828")

if dirpath:
    f = join(dirpath, "39/label_threshold_16.7269992828.mol2")
else:
    f = "39/label_threshold_16.7269992828.mol2"

cmd.load(f, 'label_threshold_16.7269992828')
cmd.hide('everything', 'label_threshold_16.7269992828')
cmd.label("label_threshold_16.7269992828", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.7269992828', members= 'label_threshold_16.7269992828')


if dirpath:
    f = join(dirpath, '39/mesh.grd')
else:
    f = '39/mesh.grd'
cmd.load(f, 'mesh_39')
cmd.isomesh("isomesh_39", "mesh_39", 0.9)
cmd.color("grey80", "isomesh_39")
cmd.set('transparency', 0.4, "isomesh_39")

cmd.group('hotspot_39', "isomesh_39")
cmd.group('hotspot_39', "mesh_39")

if dirpath:
    f = join(dirpath, "40/label_threshold_15.5.mol2")
else:
    f = "40/label_threshold_15.5.mol2"

cmd.load(f, 'label_threshold_15.5')
cmd.hide('everything', 'label_threshold_15.5')
cmd.label("label_threshold_15.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.5]
gfiles = ['40/donor.grd', '40/apolar.grd', '40/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 40
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


cluster_dict = {"16.7250003815":[], "16.7250003815_arrows":[]}

cluster_dict["16.7250003815"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-56.0), float(-8.5), float(49.5), float(1.0)]

cluster_dict["16.7250003815_arrows"] += cgo_arrow([-56.0,-8.5,49.5], [-58.938,-9.438,49.082], color="blue red", name="Arrows_16.7250003815_1")

cluster_dict["16.7250003815"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-55.0), float(-12.5), float(51.5), float(1.0)]

cluster_dict["16.7250003815_arrows"] += cgo_arrow([-55.0,-12.5,51.5], [-52.315,-11.683,50.879], color="blue red", name="Arrows_16.7250003815_2")

cluster_dict["16.7250003815"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-55.0), float(-3.5), float(49.0), float(1.0)]

cluster_dict["16.7250003815_arrows"] += cgo_arrow([-55.0,-3.5,49.0], [-56.059,-0.837,49.615], color="blue red", name="Arrows_16.7250003815_3")

cluster_dict["16.7250003815"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-52.5), float(-1.0), float(53.5), float(1.0)]

cluster_dict["16.7250003815_arrows"] += cgo_arrow([-52.5,-1.0,53.5], [-51.425,-1.643,56.198], color="blue red", name="Arrows_16.7250003815_4")

cluster_dict["16.7250003815"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-52.5), float(-6.5), float(48.0), float(1.0)]

cluster_dict["16.7250003815_arrows"] += cgo_arrow([-52.5,-6.5,48.0], [-50.9,-8.682,47.613], color="blue red", name="Arrows_16.7250003815_5")

cluster_dict["16.7250003815"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-49.0), float(-4.5), float(48.5), float(1.0)]

cluster_dict["16.7250003815_arrows"] += cgo_arrow([-49.0,-4.5,48.5], [-48.278,-5.299,45.404], color="blue red", name="Arrows_16.7250003815_6")

cluster_dict["16.7250003815"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-48.0), float(-5.5), float(51.5), float(1.0)]

cluster_dict["16.7250003815_arrows"] += cgo_arrow([-48.0,-5.5,51.5], [-44.987,-5.308,53.292], color="blue red", name="Arrows_16.7250003815_7")

cluster_dict["16.7250003815"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-57.3867758936), float(-10.2244478619), float(53.5507579016), float(1.0)]


cluster_dict["16.7250003815"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-53.4252402571), float(-5.92656648334), float(50.5875151974), float(1.0)]


cluster_dict["16.7250003815"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-49.1889270934), float(-18.3946351804), float(52.3740779799), float(1.0)]


cluster_dict["16.7250003815"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-56.5), float(-4.5), float(50.5), float(1.0)]

cluster_dict["16.7250003815_arrows"] += cgo_arrow([-56.5,-4.5,50.5], [-59.288,-5.348,50.967], color="red blue", name="Arrows_16.7250003815_8")

cluster_dict["16.7250003815"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-53.0), float(-7.5), float(49.0), float(1.0)]

cluster_dict["16.7250003815_arrows"] += cgo_arrow([-53.0,-7.5,49.0], [-50.9,-8.682,47.613], color="red blue", name="Arrows_16.7250003815_9")

cmd.load_cgo(cluster_dict["16.7250003815"], "Features_16.7250003815", 1)
cmd.load_cgo(cluster_dict["16.7250003815_arrows"], "Arrows_16.7250003815")
cmd.set("transparency", 0.2,"Features_16.7250003815")
cmd.group("Pharmacophore_16.7250003815", members="Features_16.7250003815")
cmd.group("Pharmacophore_16.7250003815", members="Arrows_16.7250003815")

if dirpath:
    f = join(dirpath, "40/label_threshold_16.7250003815.mol2")
else:
    f = "40/label_threshold_16.7250003815.mol2"

cmd.load(f, 'label_threshold_16.7250003815')
cmd.hide('everything', 'label_threshold_16.7250003815')
cmd.label("label_threshold_16.7250003815", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.7250003815', members= 'label_threshold_16.7250003815')


if dirpath:
    f = join(dirpath, '40/mesh.grd')
else:
    f = '40/mesh.grd'
cmd.load(f, 'mesh_40')
cmd.isomesh("isomesh_40", "mesh_40", 0.9)
cmd.color("grey80", "isomesh_40")
cmd.set('transparency', 0.4, "isomesh_40")

cmd.group('hotspot_40', "isomesh_40")
cmd.group('hotspot_40', "mesh_40")

if dirpath:
    f = join(dirpath, "41/label_threshold_15.7.mol2")
else:
    f = "41/label_threshold_15.7.mol2"

cmd.load(f, 'label_threshold_15.7')
cmd.hide('everything', 'label_threshold_15.7')
cmd.label("label_threshold_15.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.7]
gfiles = ['41/donor.grd', '41/apolar.grd', '41/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 41
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


cluster_dict = {"16.7129993439":[], "16.7129993439_arrows":[]}

cluster_dict["16.7129993439"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(0.0), float(68.0), float(1.0)]

cluster_dict["16.7129993439_arrows"] += cgo_arrow([-37.0,0.0,68.0], [-38.407,0.945,65.456], color="blue red", name="Arrows_16.7129993439_1")

cluster_dict["16.7129993439"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(6.0), float(70.5), float(1.0)]

cluster_dict["16.7129993439_arrows"] += cgo_arrow([-37.0,6.0,70.5], [-37.849,5.602,73.03], color="blue red", name="Arrows_16.7129993439_2")

cluster_dict["16.7129993439"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(-2.5), float(73.5), float(1.0)]

cluster_dict["16.7129993439_arrows"] += cgo_arrow([-37.0,-2.5,73.5], [-35.226,-0.918,75.907], color="blue red", name="Arrows_16.7129993439_3")

cluster_dict["16.7129993439"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-38.5127664534), float(-1.91608545832), float(71.5094264979), float(1.0)]


cluster_dict["16.7129993439"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-38.8981054908), float(-7.0), float(67.6623143099), float(1.0)]


cluster_dict["16.7129993439"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-31.972308732), float(4.05823210041), float(74.1853747661), float(1.0)]


cluster_dict["16.7129993439"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-40.5), float(-3.0), float(71.5), float(1.0)]

cluster_dict["16.7129993439_arrows"] += cgo_arrow([-40.5,-3.0,71.5], [-41.252,-0.556,70.849], color="red blue", name="Arrows_16.7129993439_4")

cluster_dict["16.7129993439"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-38.0), float(-5.5), float(63.0), float(1.0)]

cluster_dict["16.7129993439_arrows"] += cgo_arrow([-38.0,-5.5,63.0], [-35.514,-4.848,62.445], color="red blue", name="Arrows_16.7129993439_5")

cluster_dict["16.7129993439"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-35.5), float(4.0), float(70.5), float(1.0)]

cluster_dict["16.7129993439_arrows"] += cgo_arrow([-35.5,4.0,70.5], [-33.782,2.549,68.962], color="red blue", name="Arrows_16.7129993439_6")

cluster_dict["16.7129993439"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-35.0), float(4.5), float(67.0), float(1.0)]

cluster_dict["16.7129993439_arrows"] += cgo_arrow([-35.0,4.5,67.0], [-33.782,2.549,68.962], color="red blue", name="Arrows_16.7129993439_7")

cluster_dict["16.7129993439"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-36.0), float(6.0), float(65.0), float(1.0)]

cluster_dict["16.7129993439_arrows"] += cgo_arrow([-36.0,6.0,65.0], [-37.903,8.797,63.64], color="red blue", name="Arrows_16.7129993439_8")

cluster_dict["16.7129993439"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-32.5), float(6.0), float(70.5), float(1.0)]

cluster_dict["16.7129993439_arrows"] += cgo_arrow([-32.5,6.0,70.5], [-33.782,2.549,68.962], color="red blue", name="Arrows_16.7129993439_9")

cmd.load_cgo(cluster_dict["16.7129993439"], "Features_16.7129993439", 1)
cmd.load_cgo(cluster_dict["16.7129993439_arrows"], "Arrows_16.7129993439")
cmd.set("transparency", 0.2,"Features_16.7129993439")
cmd.group("Pharmacophore_16.7129993439", members="Features_16.7129993439")
cmd.group("Pharmacophore_16.7129993439", members="Arrows_16.7129993439")

if dirpath:
    f = join(dirpath, "41/label_threshold_16.7129993439.mol2")
else:
    f = "41/label_threshold_16.7129993439.mol2"

cmd.load(f, 'label_threshold_16.7129993439')
cmd.hide('everything', 'label_threshold_16.7129993439')
cmd.label("label_threshold_16.7129993439", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.7129993439', members= 'label_threshold_16.7129993439')


if dirpath:
    f = join(dirpath, '41/mesh.grd')
else:
    f = '41/mesh.grd'
cmd.load(f, 'mesh_41')
cmd.isomesh("isomesh_41", "mesh_41", 0.9)
cmd.color("grey80", "isomesh_41")
cmd.set('transparency', 0.4, "isomesh_41")

cmd.group('hotspot_41', "isomesh_41")
cmd.group('hotspot_41', "mesh_41")

if dirpath:
    f = join(dirpath, "42/label_threshold_15.0.mol2")
else:
    f = "42/label_threshold_15.0.mol2"

cmd.load(f, 'label_threshold_15.0')
cmd.hide('everything', 'label_threshold_15.0')
cmd.label("label_threshold_15.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.0]
gfiles = ['42/donor.grd', '42/apolar.grd', '42/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 42
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


cluster_dict = {"16.6275005341":[], "16.6275005341_arrows":[]}

cluster_dict["16.6275005341"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-56.0), float(-8.5), float(49.0), float(1.0)]

cluster_dict["16.6275005341_arrows"] += cgo_arrow([-56.0,-8.5,49.0], [-58.938,-9.438,49.082], color="blue red", name="Arrows_16.6275005341_1")

cluster_dict["16.6275005341"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-55.0), float(-3.5), float(49.0), float(1.0)]

cluster_dict["16.6275005341_arrows"] += cgo_arrow([-55.0,-3.5,49.0], [-56.059,-0.837,49.615], color="blue red", name="Arrows_16.6275005341_2")

cluster_dict["16.6275005341"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-52.5), float(-6.5), float(48.0), float(1.0)]

cluster_dict["16.6275005341_arrows"] += cgo_arrow([-52.5,-6.5,48.0], [-50.9,-8.682,47.613], color="blue red", name="Arrows_16.6275005341_3")

cluster_dict["16.6275005341"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-48.0), float(-4.0), float(48.5), float(1.0)]

cluster_dict["16.6275005341_arrows"] += cgo_arrow([-48.0,-4.0,48.5], [-48.278,-5.299,45.404], color="blue red", name="Arrows_16.6275005341_4")

cluster_dict["16.6275005341"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-46.0), float(9.0), float(48.0), float(1.0)]

cluster_dict["16.6275005341_arrows"] += cgo_arrow([-46.0,9.0,48.0], [-47.285,10.878,46.617], color="blue red", name="Arrows_16.6275005341_5")

cluster_dict["16.6275005341"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-46.0), float(9.5), float(45.0), float(1.0)]

cluster_dict["16.6275005341_arrows"] += cgo_arrow([-46.0,9.5,45.0], [-47.285,10.878,46.617], color="blue red", name="Arrows_16.6275005341_6")

cluster_dict["16.6275005341"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-52.514160347), float(-3.74128449495), float(49.3404289655), float(1.0)]


cluster_dict["16.6275005341"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-46.6156435545), float(2.11819235224), float(52.969012571), float(1.0)]


cluster_dict["16.6275005341"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-45.7319106788), float(8.6161192276), float(49.5247088625), float(1.0)]


cluster_dict["16.6275005341"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-44.3021726387), float(7.5), float(48.4173811099), float(1.0)]


cluster_dict["16.6275005341"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-56.5), float(-4.5), float(50.5), float(1.0)]

cluster_dict["16.6275005341_arrows"] += cgo_arrow([-56.5,-4.5,50.5], [-59.288,-5.348,50.967], color="red blue", name="Arrows_16.6275005341_7")

cluster_dict["16.6275005341"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-53.0), float(-7.5), float(49.0), float(1.0)]

cluster_dict["16.6275005341_arrows"] += cgo_arrow([-53.0,-7.5,49.0], [-50.9,-8.682,47.613], color="red blue", name="Arrows_16.6275005341_8")

cluster_dict["16.6275005341"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-46.0), float(8.0), float(49.0), float(1.0)]

cluster_dict["16.6275005341_arrows"] += cgo_arrow([-46.0,8.0,49.0], [-46.33,6.271,48.222], color="red blue", name="Arrows_16.6275005341_9")

cmd.load_cgo(cluster_dict["16.6275005341"], "Features_16.6275005341", 1)
cmd.load_cgo(cluster_dict["16.6275005341_arrows"], "Arrows_16.6275005341")
cmd.set("transparency", 0.2,"Features_16.6275005341")
cmd.group("Pharmacophore_16.6275005341", members="Features_16.6275005341")
cmd.group("Pharmacophore_16.6275005341", members="Arrows_16.6275005341")

if dirpath:
    f = join(dirpath, "42/label_threshold_16.6275005341.mol2")
else:
    f = "42/label_threshold_16.6275005341.mol2"

cmd.load(f, 'label_threshold_16.6275005341')
cmd.hide('everything', 'label_threshold_16.6275005341')
cmd.label("label_threshold_16.6275005341", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.6275005341', members= 'label_threshold_16.6275005341')


if dirpath:
    f = join(dirpath, '42/mesh.grd')
else:
    f = '42/mesh.grd'
cmd.load(f, 'mesh_42')
cmd.isomesh("isomesh_42", "mesh_42", 0.9)
cmd.color("grey80", "isomesh_42")
cmd.set('transparency', 0.4, "isomesh_42")

cmd.group('hotspot_42', "isomesh_42")
cmd.group('hotspot_42', "mesh_42")

if dirpath:
    f = join(dirpath, "43/label_threshold_10.8.mol2")
else:
    f = "43/label_threshold_10.8.mol2"

cmd.load(f, 'label_threshold_10.8')
cmd.hide('everything', 'label_threshold_10.8')
cmd.label("label_threshold_10.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.8]
gfiles = ['43/donor.grd', '43/apolar.grd', '43/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 43
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


cluster_dict = {"16.468000412":[], "16.468000412_arrows":[]}

cluster_dict["16.468000412"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(44.3084486314), float(31.568911294), float(36.4662167112), float(1.0)]


cmd.load_cgo(cluster_dict["16.468000412"], "Features_16.468000412", 1)
cmd.load_cgo(cluster_dict["16.468000412_arrows"], "Arrows_16.468000412")
cmd.set("transparency", 0.2,"Features_16.468000412")
cmd.group("Pharmacophore_16.468000412", members="Features_16.468000412")
cmd.group("Pharmacophore_16.468000412", members="Arrows_16.468000412")

if dirpath:
    f = join(dirpath, "43/label_threshold_16.468000412.mol2")
else:
    f = "43/label_threshold_16.468000412.mol2"

cmd.load(f, 'label_threshold_16.468000412')
cmd.hide('everything', 'label_threshold_16.468000412')
cmd.label("label_threshold_16.468000412", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.468000412', members= 'label_threshold_16.468000412')


if dirpath:
    f = join(dirpath, '43/mesh.grd')
else:
    f = '43/mesh.grd'
cmd.load(f, 'mesh_43')
cmd.isomesh("isomesh_43", "mesh_43", 0.9)
cmd.color("grey80", "isomesh_43")
cmd.set('transparency', 0.4, "isomesh_43")

cmd.group('hotspot_43', "isomesh_43")
cmd.group('hotspot_43', "mesh_43")

if dirpath:
    f = join(dirpath, "44/label_threshold_15.6.mol2")
else:
    f = "44/label_threshold_15.6.mol2"

cmd.load(f, 'label_threshold_15.6')
cmd.hide('everything', 'label_threshold_15.6')
cmd.label("label_threshold_15.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.6]
gfiles = ['44/donor.grd', '44/apolar.grd', '44/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 44
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


cluster_dict = {"16.466999054":[], "16.466999054_arrows":[]}

cluster_dict["16.466999054"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(20.5), float(7.5), float(36.5), float(1.0)]

cluster_dict["16.466999054_arrows"] += cgo_arrow([20.5,7.5,36.5], [19.802,5.432,35.84], color="blue red", name="Arrows_16.466999054_1")

cluster_dict["16.466999054"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(6.5), float(22.5), float(1.0)]

cluster_dict["16.466999054_arrows"] += cgo_arrow([26.5,6.5,22.5], [24.383,5.1,20.493], color="blue red", name="Arrows_16.466999054_2")

cluster_dict["16.466999054"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(10.5), float(24.0), float(1.0)]

cluster_dict["16.466999054_arrows"] += cgo_arrow([26.0,10.5,24.0], [24.828,13.48,24.46], color="blue red", name="Arrows_16.466999054_3")

cluster_dict["16.466999054"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(14.5), float(36.0), float(1.0)]

cluster_dict["16.466999054_arrows"] += cgo_arrow([26.5,14.5,36.0], [28.472,16.776,37.727], color="blue red", name="Arrows_16.466999054_4")

cluster_dict["16.466999054"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(7.5), float(16.0), float(1.0)]

cluster_dict["16.466999054_arrows"] += cgo_arrow([28.0,7.5,16.0], [25.678,8.126,17.598], color="blue red", name="Arrows_16.466999054_5")

cluster_dict["16.466999054"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(4.0), float(24.5), float(1.0)]

cluster_dict["16.466999054_arrows"] += cgo_arrow([29.5,4.0,24.5], [31.281,4.479,26.704], color="blue red", name="Arrows_16.466999054_6")

cluster_dict["16.466999054"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.4718438726), float(8.53107317106), float(37.7238137515), float(1.0)]


cluster_dict["16.466999054"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.3418079096), float(11.5706214689), float(34.2175141243), float(1.0)]


cluster_dict["16.466999054"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.1156081298), float(7.47006251869), float(22.9209994333), float(1.0)]


cluster_dict["16.466999054"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(7.0), float(37.5), float(1.0)]

cluster_dict["16.466999054_arrows"] += cgo_arrow([22.5,7.0,37.5], [22.648,5.212,35.929], color="red blue", name="Arrows_16.466999054_7")

cluster_dict["16.466999054"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(5.5), float(21.0), float(1.0)]

cluster_dict["16.466999054_arrows"] += cgo_arrow([27.0,5.5,21.0], [24.383,5.1,20.493], color="red blue", name="Arrows_16.466999054_8")

cluster_dict["16.466999054"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(28.5), float(5.5), float(16.5), float(1.0)]

cluster_dict["16.466999054_arrows"] += cgo_arrow([28.5,5.5,16.5], [30.035,4.383,14.053], color="red blue", name="Arrows_16.466999054_9")

cmd.load_cgo(cluster_dict["16.466999054"], "Features_16.466999054", 1)
cmd.load_cgo(cluster_dict["16.466999054_arrows"], "Arrows_16.466999054")
cmd.set("transparency", 0.2,"Features_16.466999054")
cmd.group("Pharmacophore_16.466999054", members="Features_16.466999054")
cmd.group("Pharmacophore_16.466999054", members="Arrows_16.466999054")

if dirpath:
    f = join(dirpath, "44/label_threshold_16.466999054.mol2")
else:
    f = "44/label_threshold_16.466999054.mol2"

cmd.load(f, 'label_threshold_16.466999054')
cmd.hide('everything', 'label_threshold_16.466999054')
cmd.label("label_threshold_16.466999054", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.466999054', members= 'label_threshold_16.466999054')


if dirpath:
    f = join(dirpath, '44/mesh.grd')
else:
    f = '44/mesh.grd'
cmd.load(f, 'mesh_44')
cmd.isomesh("isomesh_44", "mesh_44", 0.9)
cmd.color("grey80", "isomesh_44")
cmd.set('transparency', 0.4, "isomesh_44")

cmd.group('hotspot_44', "isomesh_44")
cmd.group('hotspot_44', "mesh_44")

if dirpath:
    f = join(dirpath, "45/label_threshold_15.7.mol2")
else:
    f = "45/label_threshold_15.7.mol2"

cmd.load(f, 'label_threshold_15.7')
cmd.hide('everything', 'label_threshold_15.7')
cmd.label("label_threshold_15.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.7]
gfiles = ['45/donor.grd', '45/apolar.grd', '45/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 45
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


cluster_dict = {"16.4610004425":[], "16.4610004425_arrows":[]}

cluster_dict["16.4610004425"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(8.5), float(26.5), float(1.0)]

cluster_dict["16.4610004425_arrows"] += cgo_arrow([21.0,8.5,26.5], [20.263,11.401,25.376], color="blue red", name="Arrows_16.4610004425_1")

cluster_dict["16.4610004425"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(7.5), float(10.0), float(1.0)]

cluster_dict["16.4610004425_arrows"] += cgo_arrow([21.0,7.5,10.0], [23.733,6.22,10.733], color="blue red", name="Arrows_16.4610004425_2")

cluster_dict["16.4610004425"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(6.5), float(22.5), float(1.0)]

cluster_dict["16.4610004425_arrows"] += cgo_arrow([26.5,6.5,22.5], [24.383,5.1,20.493], color="blue red", name="Arrows_16.4610004425_3")

cluster_dict["16.4610004425"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(10.5), float(24.0), float(1.0)]

cluster_dict["16.4610004425_arrows"] += cgo_arrow([26.0,10.5,24.0], [24.828,13.48,24.46], color="blue red", name="Arrows_16.4610004425_4")

cluster_dict["16.4610004425"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(7.5), float(16.0), float(1.0)]

cluster_dict["16.4610004425_arrows"] += cgo_arrow([28.0,7.5,16.0], [25.678,8.126,17.598], color="blue red", name="Arrows_16.4610004425_5")

cluster_dict["16.4610004425"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(4.0), float(24.5), float(1.0)]

cluster_dict["16.4610004425_arrows"] += cgo_arrow([29.5,4.0,24.5], [31.281,4.479,26.704], color="blue red", name="Arrows_16.4610004425_6")

cluster_dict["16.4610004425"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(20.1399369534), float(12.2013626473), float(10.8732126943), float(1.0)]


cluster_dict["16.4610004425"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.2416666667), float(10.7083333333), float(32.9166666667), float(1.0)]


cluster_dict["16.4610004425"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.0924161791), float(7.46689519444), float(22.8420218614), float(1.0)]


cluster_dict["16.4610004425"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(8.0), float(31.5), float(1.0)]

cluster_dict["16.4610004425_arrows"] += cgo_arrow([23.0,8.0,31.5], [24.327,5.898,29.592], color="red blue", name="Arrows_16.4610004425_7")

cluster_dict["16.4610004425"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(5.5), float(21.0), float(1.0)]

cluster_dict["16.4610004425_arrows"] += cgo_arrow([27.0,5.5,21.0], [24.383,5.1,20.493], color="red blue", name="Arrows_16.4610004425_8")

cluster_dict["16.4610004425"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(28.5), float(5.5), float(16.5), float(1.0)]

cluster_dict["16.4610004425_arrows"] += cgo_arrow([28.5,5.5,16.5], [30.035,4.383,14.053], color="red blue", name="Arrows_16.4610004425_9")

cmd.load_cgo(cluster_dict["16.4610004425"], "Features_16.4610004425", 1)
cmd.load_cgo(cluster_dict["16.4610004425_arrows"], "Arrows_16.4610004425")
cmd.set("transparency", 0.2,"Features_16.4610004425")
cmd.group("Pharmacophore_16.4610004425", members="Features_16.4610004425")
cmd.group("Pharmacophore_16.4610004425", members="Arrows_16.4610004425")

if dirpath:
    f = join(dirpath, "45/label_threshold_16.4610004425.mol2")
else:
    f = "45/label_threshold_16.4610004425.mol2"

cmd.load(f, 'label_threshold_16.4610004425')
cmd.hide('everything', 'label_threshold_16.4610004425')
cmd.label("label_threshold_16.4610004425", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.4610004425', members= 'label_threshold_16.4610004425')


if dirpath:
    f = join(dirpath, '45/mesh.grd')
else:
    f = '45/mesh.grd'
cmd.load(f, 'mesh_45')
cmd.isomesh("isomesh_45", "mesh_45", 0.9)
cmd.color("grey80", "isomesh_45")
cmd.set('transparency', 0.4, "isomesh_45")

cmd.group('hotspot_45', "isomesh_45")
cmd.group('hotspot_45', "mesh_45")

if dirpath:
    f = join(dirpath, "46/label_threshold_15.4.mol2")
else:
    f = "46/label_threshold_15.4.mol2"

cmd.load(f, 'label_threshold_15.4')
cmd.hide('everything', 'label_threshold_15.4')
cmd.label("label_threshold_15.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.4]
gfiles = ['46/donor.grd', '46/apolar.grd', '46/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 46
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


cluster_dict = {"16.4419994354":[], "16.4419994354_arrows":[]}

cluster_dict["16.4419994354"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(9.0), float(26.0), float(1.0)]

cluster_dict["16.4419994354_arrows"] += cgo_arrow([22.0,9.0,26.0], [20.263,11.401,25.376], color="blue red", name="Arrows_16.4419994354_1")

cluster_dict["16.4419994354"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(5.0), float(25.5), float(1.0)]

cluster_dict["16.4419994354_arrows"] += cgo_arrow([23.0,5.0,25.5], [21.848,2.901,23.953], color="blue red", name="Arrows_16.4419994354_2")

cluster_dict["16.4419994354"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(6.5), float(22.5), float(1.0)]

cluster_dict["16.4419994354_arrows"] += cgo_arrow([26.5,6.5,22.5], [24.383,5.1,20.493], color="blue red", name="Arrows_16.4419994354_3")

cluster_dict["16.4419994354"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(10.5), float(24.0), float(1.0)]

cluster_dict["16.4419994354_arrows"] += cgo_arrow([26.0,10.5,24.0], [24.828,13.48,24.46], color="blue red", name="Arrows_16.4419994354_4")

cluster_dict["16.4419994354"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(7.5), float(16.0), float(1.0)]

cluster_dict["16.4419994354_arrows"] += cgo_arrow([28.0,7.5,16.0], [25.678,8.126,17.598], color="blue red", name="Arrows_16.4419994354_5")

cluster_dict["16.4419994354"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(4.0), float(24.5), float(1.0)]

cluster_dict["16.4419994354_arrows"] += cgo_arrow([29.5,4.0,24.5], [31.281,4.479,26.704], color="blue red", name="Arrows_16.4419994354_6")

cluster_dict["16.4419994354"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.0), float(5.0), float(20.0), float(1.0)]

cluster_dict["16.4419994354_arrows"] += cgo_arrow([30.0,5.0,20.0], [31.989,4.998,17.834], color="blue red", name="Arrows_16.4419994354_7")

cluster_dict["16.4419994354"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.2585955207), float(7.54896635327), float(22.8028688937), float(1.0)]


cluster_dict["16.4419994354"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(1.5), float(25.5), float(1.0)]


cluster_dict["16.4419994354"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(35.3245614035), float(9.51169590643), float(17.7368421053), float(1.0)]


cluster_dict["16.4419994354"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(5.5), float(21.0), float(1.0)]

cluster_dict["16.4419994354_arrows"] += cgo_arrow([27.0,5.5,21.0], [24.383,5.1,20.493], color="red blue", name="Arrows_16.4419994354_8")

cluster_dict["16.4419994354"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(28.5), float(5.5), float(16.5), float(1.0)]

cluster_dict["16.4419994354_arrows"] += cgo_arrow([28.5,5.5,16.5], [30.035,4.383,14.053], color="red blue", name="Arrows_16.4419994354_9")

cmd.load_cgo(cluster_dict["16.4419994354"], "Features_16.4419994354", 1)
cmd.load_cgo(cluster_dict["16.4419994354_arrows"], "Arrows_16.4419994354")
cmd.set("transparency", 0.2,"Features_16.4419994354")
cmd.group("Pharmacophore_16.4419994354", members="Features_16.4419994354")
cmd.group("Pharmacophore_16.4419994354", members="Arrows_16.4419994354")

if dirpath:
    f = join(dirpath, "46/label_threshold_16.4419994354.mol2")
else:
    f = "46/label_threshold_16.4419994354.mol2"

cmd.load(f, 'label_threshold_16.4419994354')
cmd.hide('everything', 'label_threshold_16.4419994354')
cmd.label("label_threshold_16.4419994354", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.4419994354', members= 'label_threshold_16.4419994354')


if dirpath:
    f = join(dirpath, '46/mesh.grd')
else:
    f = '46/mesh.grd'
cmd.load(f, 'mesh_46')
cmd.isomesh("isomesh_46", "mesh_46", 0.9)
cmd.color("grey80", "isomesh_46")
cmd.set('transparency', 0.4, "isomesh_46")

cmd.group('hotspot_46', "isomesh_46")
cmd.group('hotspot_46', "mesh_46")

if dirpath:
    f = join(dirpath, "47/label_threshold_14.9.mol2")
else:
    f = "47/label_threshold_14.9.mol2"

cmd.load(f, 'label_threshold_14.9')
cmd.hide('everything', 'label_threshold_14.9')
cmd.label("label_threshold_14.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [14.9]
gfiles = ['47/donor.grd', '47/apolar.grd', '47/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 47
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


cluster_dict = {"16.4009990692":[], "16.4009990692_arrows":[]}

cluster_dict["16.4009990692"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-41.5), float(1.0), float(75.0), float(1.0)]

cluster_dict["16.4009990692_arrows"] += cgo_arrow([-41.5,1.0,75.0], [-40.56,4.447,76.446], color="blue red", name="Arrows_16.4009990692_1")

cluster_dict["16.4009990692"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-40.5), float(-3.0), float(69.5), float(1.0)]

cluster_dict["16.4009990692_arrows"] += cgo_arrow([-40.5,-3.0,69.5], [-41.252,-0.556,70.849], color="blue red", name="Arrows_16.4009990692_2")

cluster_dict["16.4009990692"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-39.5), float(-2.0), float(71.5), float(1.0)]

cluster_dict["16.4009990692_arrows"] += cgo_arrow([-39.5,-2.0,71.5], [-41.252,-0.556,70.849], color="blue red", name="Arrows_16.4009990692_3")

cluster_dict["16.4009990692"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(0.0), float(68.0), float(1.0)]

cluster_dict["16.4009990692_arrows"] += cgo_arrow([-37.0,0.0,68.0], [-38.407,0.945,65.456], color="blue red", name="Arrows_16.4009990692_4")

cluster_dict["16.4009990692"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(-6.5), float(69.0), float(1.0)]

cluster_dict["16.4009990692_arrows"] += cgo_arrow([-37.0,-6.5,69.0], [-35.002,-6.704,67.036], color="blue red", name="Arrows_16.4009990692_5")

cluster_dict["16.4009990692"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(-2.5), float(73.5), float(1.0)]

cluster_dict["16.4009990692_arrows"] += cgo_arrow([-37.0,-2.5,73.5], [-35.226,-0.918,75.907], color="blue red", name="Arrows_16.4009990692_6")

cluster_dict["16.4009990692"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-41.9040966894), float(-9.19860713302), float(68.3046610712), float(1.0)]


cluster_dict["16.4009990692"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-39.3892303926), float(-2.51203380601), float(72.3042708787), float(1.0)]


cluster_dict["16.4009990692"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-40.5), float(-3.0), float(71.5), float(1.0)]

cluster_dict["16.4009990692_arrows"] += cgo_arrow([-40.5,-3.0,71.5], [-41.252,-0.556,70.849], color="red blue", name="Arrows_16.4009990692_7")

cluster_dict["16.4009990692"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-38.5), float(-6.5), float(74.0), float(1.0)]

cluster_dict["16.4009990692_arrows"] += cgo_arrow([-38.5,-6.5,74.0], [-36.711,-4.673,75.752], color="red blue", name="Arrows_16.4009990692_8")

cluster_dict["16.4009990692"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-38.0), float(-5.5), float(63.0), float(1.0)]

cluster_dict["16.4009990692_arrows"] += cgo_arrow([-38.0,-5.5,63.0], [-35.514,-4.848,62.445], color="red blue", name="Arrows_16.4009990692_9")

cmd.load_cgo(cluster_dict["16.4009990692"], "Features_16.4009990692", 1)
cmd.load_cgo(cluster_dict["16.4009990692_arrows"], "Arrows_16.4009990692")
cmd.set("transparency", 0.2,"Features_16.4009990692")
cmd.group("Pharmacophore_16.4009990692", members="Features_16.4009990692")
cmd.group("Pharmacophore_16.4009990692", members="Arrows_16.4009990692")

if dirpath:
    f = join(dirpath, "47/label_threshold_16.4009990692.mol2")
else:
    f = "47/label_threshold_16.4009990692.mol2"

cmd.load(f, 'label_threshold_16.4009990692')
cmd.hide('everything', 'label_threshold_16.4009990692')
cmd.label("label_threshold_16.4009990692", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.4009990692', members= 'label_threshold_16.4009990692')


if dirpath:
    f = join(dirpath, '47/mesh.grd')
else:
    f = '47/mesh.grd'
cmd.load(f, 'mesh_47')
cmd.isomesh("isomesh_47", "mesh_47", 0.9)
cmd.color("grey80", "isomesh_47")
cmd.set('transparency', 0.4, "isomesh_47")

cmd.group('hotspot_47', "isomesh_47")
cmd.group('hotspot_47', "mesh_47")

if dirpath:
    f = join(dirpath, "48/label_threshold_14.9.mol2")
else:
    f = "48/label_threshold_14.9.mol2"

cmd.load(f, 'label_threshold_14.9')
cmd.hide('everything', 'label_threshold_14.9')
cmd.label("label_threshold_14.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [14.9]
gfiles = ['48/donor.grd', '48/apolar.grd', '48/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 48
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


cluster_dict = {"16.3659992218":[], "16.3659992218_arrows":[]}

cluster_dict["16.3659992218"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-40.5), float(-3.0), float(69.5), float(1.0)]

cluster_dict["16.3659992218_arrows"] += cgo_arrow([-40.5,-3.0,69.5], [-41.252,-0.556,70.849], color="blue red", name="Arrows_16.3659992218_1")

cluster_dict["16.3659992218"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-39.5), float(-2.0), float(71.5), float(1.0)]

cluster_dict["16.3659992218_arrows"] += cgo_arrow([-39.5,-2.0,71.5], [-41.252,-0.556,70.849], color="blue red", name="Arrows_16.3659992218_2")

cluster_dict["16.3659992218"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(0.0), float(68.0), float(1.0)]

cluster_dict["16.3659992218_arrows"] += cgo_arrow([-37.0,0.0,68.0], [-38.407,0.945,65.456], color="blue red", name="Arrows_16.3659992218_3")

cluster_dict["16.3659992218"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(-6.5), float(69.0), float(1.0)]

cluster_dict["16.3659992218_arrows"] += cgo_arrow([-37.0,-6.5,69.0], [-35.002,-6.704,67.036], color="blue red", name="Arrows_16.3659992218_4")

cluster_dict["16.3659992218"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(-2.5), float(73.5), float(1.0)]

cluster_dict["16.3659992218_arrows"] += cgo_arrow([-37.0,-2.5,73.5], [-35.226,-0.918,75.907], color="blue red", name="Arrows_16.3659992218_5")

cluster_dict["16.3659992218"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-39.9769157481), float(-4.83180071927), float(70.6990508832), float(1.0)]


cluster_dict["16.3659992218"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-44.5), float(-3.5), float(76.0), float(1.0)]

cluster_dict["16.3659992218_arrows"] += cgo_arrow([-44.5,-3.5,76.0], [-46.114,-1.645,74.864], color="red blue", name="Arrows_16.3659992218_6")

cluster_dict["16.3659992218"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-40.5), float(-3.0), float(71.5), float(1.0)]

cluster_dict["16.3659992218_arrows"] += cgo_arrow([-40.5,-3.0,71.5], [-41.252,-0.556,70.849], color="red blue", name="Arrows_16.3659992218_7")

cluster_dict["16.3659992218"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-38.5), float(-6.5), float(74.0), float(1.0)]

cluster_dict["16.3659992218_arrows"] += cgo_arrow([-38.5,-6.5,74.0], [-36.711,-4.673,75.752], color="red blue", name="Arrows_16.3659992218_8")

cluster_dict["16.3659992218"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-38.0), float(-5.5), float(63.0), float(1.0)]

cluster_dict["16.3659992218_arrows"] += cgo_arrow([-38.0,-5.5,63.0], [-35.514,-4.848,62.445], color="red blue", name="Arrows_16.3659992218_9")

cmd.load_cgo(cluster_dict["16.3659992218"], "Features_16.3659992218", 1)
cmd.load_cgo(cluster_dict["16.3659992218_arrows"], "Arrows_16.3659992218")
cmd.set("transparency", 0.2,"Features_16.3659992218")
cmd.group("Pharmacophore_16.3659992218", members="Features_16.3659992218")
cmd.group("Pharmacophore_16.3659992218", members="Arrows_16.3659992218")

if dirpath:
    f = join(dirpath, "48/label_threshold_16.3659992218.mol2")
else:
    f = "48/label_threshold_16.3659992218.mol2"

cmd.load(f, 'label_threshold_16.3659992218')
cmd.hide('everything', 'label_threshold_16.3659992218')
cmd.label("label_threshold_16.3659992218", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.3659992218', members= 'label_threshold_16.3659992218')


if dirpath:
    f = join(dirpath, '48/mesh.grd')
else:
    f = '48/mesh.grd'
cmd.load(f, 'mesh_48')
cmd.isomesh("isomesh_48", "mesh_48", 0.9)
cmd.color("grey80", "isomesh_48")
cmd.set('transparency', 0.4, "isomesh_48")

cmd.group('hotspot_48', "isomesh_48")
cmd.group('hotspot_48', "mesh_48")

if dirpath:
    f = join(dirpath, "49/label_threshold_15.0.mol2")
else:
    f = "49/label_threshold_15.0.mol2"

cmd.load(f, 'label_threshold_15.0')
cmd.hide('everything', 'label_threshold_15.0')
cmd.label("label_threshold_15.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.0]
gfiles = ['49/donor.grd', '49/apolar.grd', '49/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 49
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


cluster_dict = {"16.3349990845":[], "16.3349990845_arrows":[]}

cluster_dict["16.3349990845"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(6.0), float(9.0), float(42.5), float(1.0)]

cluster_dict["16.3349990845_arrows"] += cgo_arrow([6.0,9.0,42.5], [3.715,8.369,43.968], color="blue red", name="Arrows_16.3349990845_1")

cluster_dict["16.3349990845"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(9.0), float(6.5), float(43.0), float(1.0)]

cluster_dict["16.3349990845_arrows"] += cgo_arrow([9.0,6.5,43.0], [10.755,4.658,42.464], color="blue red", name="Arrows_16.3349990845_2")

cluster_dict["16.3349990845"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(10.0), float(7.0), float(44.5), float(1.0)]

cluster_dict["16.3349990845_arrows"] += cgo_arrow([10.0,7.0,44.5], [9.073,4.46,45.522], color="blue red", name="Arrows_16.3349990845_3")

cluster_dict["16.3349990845"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(11.0), float(9.5), float(41.5), float(1.0)]

cluster_dict["16.3349990845_arrows"] += cgo_arrow([11.0,9.5,41.5], [11.409,11.758,39.462], color="blue red", name="Arrows_16.3349990845_4")

cluster_dict["16.3349990845"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(3.64165303242), float(10.4715019872), float(47.8968019702), float(1.0)]


cluster_dict["16.3349990845"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(7.95145132269), float(7.89890354036), float(44.2028458091), float(1.0)]


cluster_dict["16.3349990845"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(3.5), float(0.0), float(41.0), float(1.0)]

cluster_dict["16.3349990845_arrows"] += cgo_arrow([3.5,0.0,41.0], [3.075,-1.763,38.461], color="red blue", name="Arrows_16.3349990845_5")

cluster_dict["16.3349990845"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(4.0), float(11.5), float(46.0), float(1.0)]

cluster_dict["16.3349990845_arrows"] += cgo_arrow([4.0,11.5,46.0], [4.114,14.533,46.591], color="red blue", name="Arrows_16.3349990845_6")

cluster_dict["16.3349990845"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(5.5), float(8.0), float(43.0), float(1.0)]

cluster_dict["16.3349990845_arrows"] += cgo_arrow([5.5,8.0,43.0], [3.715,8.369,43.968], color="red blue", name="Arrows_16.3349990845_7")

cluster_dict["16.3349990845"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(10.0), float(6.5), float(50.5), float(1.0)]

cluster_dict["16.3349990845_arrows"] += cgo_arrow([10.0,6.5,50.5], [7.899,6.991,52.135], color="red blue", name="Arrows_16.3349990845_8")

cluster_dict["16.3349990845"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(13.5), float(13.0), float(45.0), float(1.0)]

cluster_dict["16.3349990845_arrows"] += cgo_arrow([13.5,13.0,45.0], [14.379,15.215,43.186], color="red blue", name="Arrows_16.3349990845_9")

cmd.load_cgo(cluster_dict["16.3349990845"], "Features_16.3349990845", 1)
cmd.load_cgo(cluster_dict["16.3349990845_arrows"], "Arrows_16.3349990845")
cmd.set("transparency", 0.2,"Features_16.3349990845")
cmd.group("Pharmacophore_16.3349990845", members="Features_16.3349990845")
cmd.group("Pharmacophore_16.3349990845", members="Arrows_16.3349990845")

if dirpath:
    f = join(dirpath, "49/label_threshold_16.3349990845.mol2")
else:
    f = "49/label_threshold_16.3349990845.mol2"

cmd.load(f, 'label_threshold_16.3349990845')
cmd.hide('everything', 'label_threshold_16.3349990845')
cmd.label("label_threshold_16.3349990845", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.3349990845', members= 'label_threshold_16.3349990845')


if dirpath:
    f = join(dirpath, '49/mesh.grd')
else:
    f = '49/mesh.grd'
cmd.load(f, 'mesh_49')
cmd.isomesh("isomesh_49", "mesh_49", 0.9)
cmd.color("grey80", "isomesh_49")
cmd.set('transparency', 0.4, "isomesh_49")

cmd.group('hotspot_49', "isomesh_49")
cmd.group('hotspot_49', "mesh_49")

if dirpath:
    f = join(dirpath, "50/label_threshold_15.0.mol2")
else:
    f = "50/label_threshold_15.0.mol2"

cmd.load(f, 'label_threshold_15.0')
cmd.hide('everything', 'label_threshold_15.0')
cmd.label("label_threshold_15.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.0]
gfiles = ['50/donor.grd', '50/apolar.grd', '50/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 50
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


cluster_dict = {"16.3320007324":[], "16.3320007324_arrows":[]}

cluster_dict["16.3320007324"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(6.0), float(9.0), float(42.5), float(1.0)]

cluster_dict["16.3320007324_arrows"] += cgo_arrow([6.0,9.0,42.5], [3.715,8.369,43.968], color="blue red", name="Arrows_16.3320007324_1")

cluster_dict["16.3320007324"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(9.0), float(6.5), float(43.0), float(1.0)]

cluster_dict["16.3320007324_arrows"] += cgo_arrow([9.0,6.5,43.0], [10.755,4.658,42.464], color="blue red", name="Arrows_16.3320007324_2")

cluster_dict["16.3320007324"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(11.0), float(9.5), float(41.5), float(1.0)]

cluster_dict["16.3320007324_arrows"] += cgo_arrow([11.0,9.5,41.5], [11.409,11.758,39.462], color="blue red", name="Arrows_16.3320007324_3")

cluster_dict["16.3320007324"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(3.64165303242), float(10.4715019872), float(47.8968019702), float(1.0)]


cluster_dict["16.3320007324"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(7.95509658616), float(7.87673549401), float(44.2074561156), float(1.0)]


cluster_dict["16.3320007324"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(3.5), float(0.0), float(41.0), float(1.0)]

cluster_dict["16.3320007324_arrows"] += cgo_arrow([3.5,0.0,41.0], [3.075,-1.763,38.461], color="red blue", name="Arrows_16.3320007324_4")

cluster_dict["16.3320007324"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(4.0), float(11.5), float(46.0), float(1.0)]

cluster_dict["16.3320007324_arrows"] += cgo_arrow([4.0,11.5,46.0], [4.114,14.533,46.591], color="red blue", name="Arrows_16.3320007324_5")

cluster_dict["16.3320007324"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(5.5), float(8.0), float(43.0), float(1.0)]

cluster_dict["16.3320007324_arrows"] += cgo_arrow([5.5,8.0,43.0], [3.715,8.369,43.968], color="red blue", name="Arrows_16.3320007324_6")

cluster_dict["16.3320007324"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(10.0), float(6.5), float(50.5), float(1.0)]

cluster_dict["16.3320007324_arrows"] += cgo_arrow([10.0,6.5,50.5], [7.899,6.991,52.135], color="red blue", name="Arrows_16.3320007324_7")

cluster_dict["16.3320007324"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(12.0), float(5.5), float(52.0), float(1.0)]

cluster_dict["16.3320007324_arrows"] += cgo_arrow([12.0,5.5,52.0], [11.081,7.326,54.016], color="red blue", name="Arrows_16.3320007324_8")

cluster_dict["16.3320007324"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(13.5), float(13.0), float(45.0), float(1.0)]

cluster_dict["16.3320007324_arrows"] += cgo_arrow([13.5,13.0,45.0], [14.379,15.215,43.186], color="red blue", name="Arrows_16.3320007324_9")

cmd.load_cgo(cluster_dict["16.3320007324"], "Features_16.3320007324", 1)
cmd.load_cgo(cluster_dict["16.3320007324_arrows"], "Arrows_16.3320007324")
cmd.set("transparency", 0.2,"Features_16.3320007324")
cmd.group("Pharmacophore_16.3320007324", members="Features_16.3320007324")
cmd.group("Pharmacophore_16.3320007324", members="Arrows_16.3320007324")

if dirpath:
    f = join(dirpath, "50/label_threshold_16.3320007324.mol2")
else:
    f = "50/label_threshold_16.3320007324.mol2"

cmd.load(f, 'label_threshold_16.3320007324')
cmd.hide('everything', 'label_threshold_16.3320007324')
cmd.label("label_threshold_16.3320007324", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.3320007324', members= 'label_threshold_16.3320007324')


if dirpath:
    f = join(dirpath, '50/mesh.grd')
else:
    f = '50/mesh.grd'
cmd.load(f, 'mesh_50')
cmd.isomesh("isomesh_50", "mesh_50", 0.9)
cmd.color("grey80", "isomesh_50")
cmd.set('transparency', 0.4, "isomesh_50")

cmd.group('hotspot_50', "isomesh_50")
cmd.group('hotspot_50', "mesh_50")

if dirpath:
    f = join(dirpath, "51/label_threshold_14.9.mol2")
else:
    f = "51/label_threshold_14.9.mol2"

cmd.load(f, 'label_threshold_14.9')
cmd.hide('everything', 'label_threshold_14.9')
cmd.label("label_threshold_14.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [14.9]
gfiles = ['51/donor.grd', '51/apolar.grd', '51/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 51
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


cluster_dict = {"16.329000473":[], "16.329000473_arrows":[]}

cluster_dict["16.329000473"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(0.0), float(68.0), float(1.0)]

cluster_dict["16.329000473_arrows"] += cgo_arrow([-37.0,0.0,68.0], [-38.407,0.945,65.456], color="blue red", name="Arrows_16.329000473_1")

cluster_dict["16.329000473"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(6.0), float(70.5), float(1.0)]

cluster_dict["16.329000473_arrows"] += cgo_arrow([-37.0,6.0,70.5], [-37.849,5.602,73.03], color="blue red", name="Arrows_16.329000473_2")

cluster_dict["16.329000473"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-37.0), float(-2.5), float(73.5), float(1.0)]

cluster_dict["16.329000473_arrows"] += cgo_arrow([-37.0,-2.5,73.5], [-35.226,-0.918,75.907], color="blue red", name="Arrows_16.329000473_3")

cluster_dict["16.329000473"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-28.5), float(3.0), float(74.5), float(1.0)]

cluster_dict["16.329000473_arrows"] += cgo_arrow([-28.5,3.0,74.5], [-29.451,0.249,74.504], color="blue red", name="Arrows_16.329000473_4")

cluster_dict["16.329000473"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-32.3854228263), float(3.67970921751), float(73.9991577272), float(1.0)]


cluster_dict["16.329000473"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-34.1079985211), float(5.0214856723), float(65.132123156), float(1.0)]


cluster_dict["16.329000473"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-35.5), float(4.0), float(70.5), float(1.0)]

cluster_dict["16.329000473_arrows"] += cgo_arrow([-35.5,4.0,70.5], [-33.782,2.549,68.962], color="red blue", name="Arrows_16.329000473_5")

cluster_dict["16.329000473"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-35.0), float(4.5), float(67.0), float(1.0)]

cluster_dict["16.329000473_arrows"] += cgo_arrow([-35.0,4.5,67.0], [-33.782,2.549,68.962], color="red blue", name="Arrows_16.329000473_6")

cluster_dict["16.329000473"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-36.0), float(6.0), float(65.0), float(1.0)]

cluster_dict["16.329000473_arrows"] += cgo_arrow([-36.0,6.0,65.0], [-37.903,8.797,63.64], color="red blue", name="Arrows_16.329000473_7")

cluster_dict["16.329000473"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-32.5), float(6.0), float(70.5), float(1.0)]

cluster_dict["16.329000473_arrows"] += cgo_arrow([-32.5,6.0,70.5], [-33.782,2.549,68.962], color="red blue", name="Arrows_16.329000473_8")

cluster_dict["16.329000473"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-32.0), float(2.5), float(75.5), float(1.0)]

cluster_dict["16.329000473_arrows"] += cgo_arrow([-32.0,2.5,75.5], [-31.095,0.678,77.69], color="red blue", name="Arrows_16.329000473_9")

cmd.load_cgo(cluster_dict["16.329000473"], "Features_16.329000473", 1)
cmd.load_cgo(cluster_dict["16.329000473_arrows"], "Arrows_16.329000473")
cmd.set("transparency", 0.2,"Features_16.329000473")
cmd.group("Pharmacophore_16.329000473", members="Features_16.329000473")
cmd.group("Pharmacophore_16.329000473", members="Arrows_16.329000473")

if dirpath:
    f = join(dirpath, "51/label_threshold_16.329000473.mol2")
else:
    f = "51/label_threshold_16.329000473.mol2"

cmd.load(f, 'label_threshold_16.329000473')
cmd.hide('everything', 'label_threshold_16.329000473')
cmd.label("label_threshold_16.329000473", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.329000473', members= 'label_threshold_16.329000473')


if dirpath:
    f = join(dirpath, '51/mesh.grd')
else:
    f = '51/mesh.grd'
cmd.load(f, 'mesh_51')
cmd.isomesh("isomesh_51", "mesh_51", 0.9)
cmd.color("grey80", "isomesh_51")
cmd.set('transparency', 0.4, "isomesh_51")

cmd.group('hotspot_51', "isomesh_51")
cmd.group('hotspot_51', "mesh_51")

if dirpath:
    f = join(dirpath, "52/label_threshold_15.1.mol2")
else:
    f = "52/label_threshold_15.1.mol2"

cmd.load(f, 'label_threshold_15.1')
cmd.hide('everything', 'label_threshold_15.1')
cmd.label("label_threshold_15.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.1]
gfiles = ['52/donor.grd', '52/apolar.grd', '52/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 52
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


cluster_dict = {"16.3229999542":[], "16.3229999542_arrows":[]}

cluster_dict["16.3229999542"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(12.0), float(33.5), float(39.0), float(1.0)]

cluster_dict["16.3229999542_arrows"] += cgo_arrow([12.0,33.5,39.0], [11.415,33.832,35.865], color="blue red", name="Arrows_16.3229999542_1")

cluster_dict["16.3229999542"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(16.5), float(26.0), float(34.0), float(1.0)]

cluster_dict["16.3229999542_arrows"] += cgo_arrow([16.5,26.0,34.0], [17.727,23.511,33.873], color="blue red", name="Arrows_16.3229999542_2")

cluster_dict["16.3229999542"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(28.0), float(27.0), float(1.0)]

cluster_dict["16.3229999542_arrows"] += cgo_arrow([15.0,28.0,27.0], [13.275,26.581,27.115], color="blue red", name="Arrows_16.3229999542_3")

cluster_dict["16.3229999542"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(16.5), float(26.0), float(34.0), float(1.0)]

cluster_dict["16.3229999542_arrows"] += cgo_arrow([16.5,26.0,34.0], [17.727,23.511,33.873], color="blue red", name="Arrows_16.3229999542_4")

cluster_dict["16.3229999542"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.0), float(32.0), float(34.0), float(1.0)]

cluster_dict["16.3229999542_arrows"] += cgo_arrow([19.0,32.0,34.0], [16.762,32.791,33.752], color="blue red", name="Arrows_16.3229999542_5")

cluster_dict["16.3229999542"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(11.1822025241), float(32.4927801981), float(40.6361080125), float(1.0)]


cluster_dict["16.3229999542"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(16.4634111463), float(26.8794766978), float(34.9644706115), float(1.0)]


cluster_dict["16.3229999542"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(16.4862847443), float(25.4130947578), float(26.5946235073), float(1.0)]


cluster_dict["16.3229999542"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.4633132989), float(21.1564002118), float(44.4303078464), float(1.0)]


cluster_dict["16.3229999542"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(11.5), float(22.5), float(38.5), float(1.0)]

cluster_dict["16.3229999542_arrows"] += cgo_arrow([11.5,22.5,38.5], [10.91,21.842,41.086], color="red blue", name="Arrows_16.3229999542_6")

cluster_dict["16.3229999542"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(12.0), float(30.5), float(39.0), float(1.0)]

cluster_dict["16.3229999542_arrows"] += cgo_arrow([12.0,30.5,39.0], [13.35,28.307,37.593], color="red blue", name="Arrows_16.3229999542_7")

cluster_dict["16.3229999542"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(13.0), float(20.5), float(39.5), float(1.0)]

cluster_dict["16.3229999542_arrows"] += cgo_arrow([13.0,20.5,39.5], [10.91,21.842,41.086], color="red blue", name="Arrows_16.3229999542_8")

cluster_dict["16.3229999542"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(27.0), float(29.5), float(1.0)]

cluster_dict["16.3229999542_arrows"] += cgo_arrow([16.0,27.0,29.5], [12.214,25.278,28.555], color="red blue", name="Arrows_16.3229999542_9")

cmd.load_cgo(cluster_dict["16.3229999542"], "Features_16.3229999542", 1)
cmd.load_cgo(cluster_dict["16.3229999542_arrows"], "Arrows_16.3229999542")
cmd.set("transparency", 0.2,"Features_16.3229999542")
cmd.group("Pharmacophore_16.3229999542", members="Features_16.3229999542")
cmd.group("Pharmacophore_16.3229999542", members="Arrows_16.3229999542")

if dirpath:
    f = join(dirpath, "52/label_threshold_16.3229999542.mol2")
else:
    f = "52/label_threshold_16.3229999542.mol2"

cmd.load(f, 'label_threshold_16.3229999542')
cmd.hide('everything', 'label_threshold_16.3229999542')
cmd.label("label_threshold_16.3229999542", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.3229999542', members= 'label_threshold_16.3229999542')


if dirpath:
    f = join(dirpath, '52/mesh.grd')
else:
    f = '52/mesh.grd'
cmd.load(f, 'mesh_52')
cmd.isomesh("isomesh_52", "mesh_52", 0.9)
cmd.color("grey80", "isomesh_52")
cmd.set('transparency', 0.4, "isomesh_52")

cmd.group('hotspot_52', "isomesh_52")
cmd.group('hotspot_52', "mesh_52")

if dirpath:
    f = join(dirpath, "53/label_threshold_14.2.mol2")
else:
    f = "53/label_threshold_14.2.mol2"

cmd.load(f, 'label_threshold_14.2')
cmd.hide('everything', 'label_threshold_14.2')
cmd.label("label_threshold_14.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [14.2]
gfiles = ['53/donor.grd', '53/apolar.grd', '53/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 53
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


cluster_dict = {"16.2280006409":[], "16.2280006409_arrows":[]}

cluster_dict["16.2280006409"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(12.0), float(33.5), float(39.0), float(1.0)]

cluster_dict["16.2280006409_arrows"] += cgo_arrow([12.0,33.5,39.0], [11.415,33.832,35.865], color="blue red", name="Arrows_16.2280006409_1")

cluster_dict["16.2280006409"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(36.5), float(43.0), float(1.0)]

cluster_dict["16.2280006409_arrows"] += cgo_arrow([15.0,36.5,43.0], [16.087,33.906,43.704], color="blue red", name="Arrows_16.2280006409_2")

cluster_dict["16.2280006409"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(17.0), float(37.5), float(40.5), float(1.0)]

cluster_dict["16.2280006409_arrows"] += cgo_arrow([17.0,37.5,40.5], [16.019,37.942,37.737], color="blue red", name="Arrows_16.2280006409_3")

cluster_dict["16.2280006409"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(29.5), float(39.5), float(1.0)]

cluster_dict["16.2280006409_arrows"] += cgo_arrow([17.5,29.5,39.5], [19.6,30.549,39.07], color="blue red", name="Arrows_16.2280006409_4")

cluster_dict["16.2280006409"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.0), float(29.5), float(36.5), float(1.0)]

cluster_dict["16.2280006409_arrows"] += cgo_arrow([19.0,29.5,36.5], [19.6,30.549,39.07], color="blue red", name="Arrows_16.2280006409_5")

cluster_dict["16.2280006409"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(11.4237813717), float(35.2931543491), float(41.391286656), float(1.0)]


cluster_dict["16.2280006409"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(16.7028350381), float(29.8766050452), float(36.6731262164), float(1.0)]


cluster_dict["16.2280006409"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(7.5), float(32.0), float(44.5), float(1.0)]

cluster_dict["16.2280006409_arrows"] += cgo_arrow([7.5,32.0,44.5], [11.491,29.94,43.112], color="red blue", name="Arrows_16.2280006409_6")

cluster_dict["16.2280006409"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(8.5), float(38.5), float(39.5), float(1.0)]

cluster_dict["16.2280006409_arrows"] += cgo_arrow([8.5,38.5,39.5], [10.654,40.531,39.037], color="red blue", name="Arrows_16.2280006409_7")

cluster_dict["16.2280006409"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(12.0), float(30.5), float(39.0), float(1.0)]

cluster_dict["16.2280006409_arrows"] += cgo_arrow([12.0,30.5,39.0], [13.35,28.307,37.593], color="red blue", name="Arrows_16.2280006409_8")

cluster_dict["16.2280006409"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(14.5), float(31.0), float(40.0), float(1.0)]

cluster_dict["16.2280006409_arrows"] += cgo_arrow([14.5,31.0,40.0], [13.35,28.307,37.593], color="red blue", name="Arrows_16.2280006409_9")

cmd.load_cgo(cluster_dict["16.2280006409"], "Features_16.2280006409", 1)
cmd.load_cgo(cluster_dict["16.2280006409_arrows"], "Arrows_16.2280006409")
cmd.set("transparency", 0.2,"Features_16.2280006409")
cmd.group("Pharmacophore_16.2280006409", members="Features_16.2280006409")
cmd.group("Pharmacophore_16.2280006409", members="Arrows_16.2280006409")

if dirpath:
    f = join(dirpath, "53/label_threshold_16.2280006409.mol2")
else:
    f = "53/label_threshold_16.2280006409.mol2"

cmd.load(f, 'label_threshold_16.2280006409')
cmd.hide('everything', 'label_threshold_16.2280006409')
cmd.label("label_threshold_16.2280006409", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.2280006409', members= 'label_threshold_16.2280006409')


if dirpath:
    f = join(dirpath, '53/mesh.grd')
else:
    f = '53/mesh.grd'
cmd.load(f, 'mesh_53')
cmd.isomesh("isomesh_53", "mesh_53", 0.9)
cmd.color("grey80", "isomesh_53")
cmd.set('transparency', 0.4, "isomesh_53")

cmd.group('hotspot_53', "isomesh_53")
cmd.group('hotspot_53', "mesh_53")

if dirpath:
    f = join(dirpath, "54/label_threshold_15.0.mol2")
else:
    f = "54/label_threshold_15.0.mol2"

cmd.load(f, 'label_threshold_15.0')
cmd.hide('everything', 'label_threshold_15.0')
cmd.label("label_threshold_15.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.0]
gfiles = ['54/donor.grd', '54/apolar.grd', '54/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 54
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


cluster_dict = {"16.2019996643":[], "16.2019996643_arrows":[]}

cluster_dict["16.2019996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(6.5), float(22.5), float(1.0)]

cluster_dict["16.2019996643_arrows"] += cgo_arrow([26.5,6.5,22.5], [24.383,5.1,20.493], color="blue red", name="Arrows_16.2019996643_1")

cluster_dict["16.2019996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(10.5), float(24.0), float(1.0)]

cluster_dict["16.2019996643_arrows"] += cgo_arrow([26.0,10.5,24.0], [24.828,13.48,24.46], color="blue red", name="Arrows_16.2019996643_2")

cluster_dict["16.2019996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(7.5), float(16.0), float(1.0)]

cluster_dict["16.2019996643_arrows"] += cgo_arrow([28.0,7.5,16.0], [25.678,8.126,17.598], color="blue red", name="Arrows_16.2019996643_3")

cluster_dict["16.2019996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(4.0), float(24.5), float(1.0)]

cluster_dict["16.2019996643_arrows"] += cgo_arrow([29.5,4.0,24.5], [31.281,4.479,26.704], color="blue red", name="Arrows_16.2019996643_4")

cluster_dict["16.2019996643"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.0), float(5.0), float(20.0), float(1.0)]

cluster_dict["16.2019996643_arrows"] += cgo_arrow([30.0,5.0,20.0], [31.989,4.998,17.834], color="blue red", name="Arrows_16.2019996643_5")

cluster_dict["16.2019996643"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(28.6887177831), float(7.88665139408), float(21.5447611563), float(1.0)]


cluster_dict["16.2019996643"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(29.6666666667), float(1.5), float(25.3333333333), float(1.0)]


cluster_dict["16.2019996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(5.5), float(21.0), float(1.0)]

cluster_dict["16.2019996643_arrows"] += cgo_arrow([27.0,5.5,21.0], [24.383,5.1,20.493], color="red blue", name="Arrows_16.2019996643_6")

cluster_dict["16.2019996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(28.5), float(5.5), float(16.5), float(1.0)]

cluster_dict["16.2019996643_arrows"] += cgo_arrow([28.5,5.5,16.5], [30.035,4.383,14.053], color="red blue", name="Arrows_16.2019996643_7")

cluster_dict["16.2019996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(31.5), float(4.5), float(20.5), float(1.0)]

cluster_dict["16.2019996643_arrows"] += cgo_arrow([31.5,4.5,20.5], [33.26,6.403,21.589], color="red blue", name="Arrows_16.2019996643_8")

cluster_dict["16.2019996643"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(43.5), float(12.0), float(18.0), float(1.0)]

cluster_dict["16.2019996643_arrows"] += cgo_arrow([43.5,12.0,18.0], [41.724,10.524,16.554], color="red blue", name="Arrows_16.2019996643_9")

cmd.load_cgo(cluster_dict["16.2019996643"], "Features_16.2019996643", 1)
cmd.load_cgo(cluster_dict["16.2019996643_arrows"], "Arrows_16.2019996643")
cmd.set("transparency", 0.2,"Features_16.2019996643")
cmd.group("Pharmacophore_16.2019996643", members="Features_16.2019996643")
cmd.group("Pharmacophore_16.2019996643", members="Arrows_16.2019996643")

if dirpath:
    f = join(dirpath, "54/label_threshold_16.2019996643.mol2")
else:
    f = "54/label_threshold_16.2019996643.mol2"

cmd.load(f, 'label_threshold_16.2019996643')
cmd.hide('everything', 'label_threshold_16.2019996643')
cmd.label("label_threshold_16.2019996643", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.2019996643', members= 'label_threshold_16.2019996643')


if dirpath:
    f = join(dirpath, '54/mesh.grd')
else:
    f = '54/mesh.grd'
cmd.load(f, 'mesh_54')
cmd.isomesh("isomesh_54", "mesh_54", 0.9)
cmd.color("grey80", "isomesh_54")
cmd.set('transparency', 0.4, "isomesh_54")

cmd.group('hotspot_54', "isomesh_54")
cmd.group('hotspot_54', "mesh_54")

if dirpath:
    f = join(dirpath, "55/label_threshold_14.6.mol2")
else:
    f = "55/label_threshold_14.6.mol2"

cmd.load(f, 'label_threshold_14.6')
cmd.hide('everything', 'label_threshold_14.6')
cmd.label("label_threshold_14.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [14.6]
gfiles = ['55/donor.grd', '55/apolar.grd', '55/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 55
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


cluster_dict = {"16.1790008545":[], "16.1790008545_arrows":[]}

cluster_dict["16.1790008545"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-33.5), float(9.5), float(52.0), float(1.0)]

cluster_dict["16.1790008545_arrows"] += cgo_arrow([-33.5,9.5,52.0], [-33.337,6.587,52.943], color="blue red", name="Arrows_16.1790008545_1")

cluster_dict["16.1790008545"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-33.0), float(12.5), float(49.0), float(1.0)]

cluster_dict["16.1790008545_arrows"] += cgo_arrow([-33.0,12.5,49.0], [-31.181,13.279,46.18], color="blue red", name="Arrows_16.1790008545_2")

cluster_dict["16.1790008545"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-29.0), float(9.0), float(50.5), float(1.0)]

cluster_dict["16.1790008545_arrows"] += cgo_arrow([-29.0,9.0,50.5], [-29.403,7.109,52.477], color="blue red", name="Arrows_16.1790008545_3")

cluster_dict["16.1790008545"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-30.4675880721), float(10.3187859129), float(48.59071801), float(1.0)]


cluster_dict["16.1790008545"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-34.833344282), float(7.83337712809), float(50.333311436), float(1.0)]


cluster_dict["16.1790008545"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-25.1975127975), float(14.0680886707), float(51.2770480854), float(1.0)]


cluster_dict["16.1790008545"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-32.5), float(11.5), float(55.0), float(1.0)]

cluster_dict["16.1790008545_arrows"] += cgo_arrow([-32.5,11.5,55.0], [-34.887,10.154,57.128], color="red blue", name="Arrows_16.1790008545_4")

cluster_dict["16.1790008545"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-32.5), float(12.0), float(48.0), float(1.0)]

cluster_dict["16.1790008545_arrows"] += cgo_arrow([-32.5,12.0,48.0], [-31.181,13.279,46.18], color="red blue", name="Arrows_16.1790008545_5")

cluster_dict["16.1790008545"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-28.5), float(9.0), float(51.0), float(1.0)]

cluster_dict["16.1790008545_arrows"] += cgo_arrow([-28.5,9.0,51.0], [-29.403,7.109,52.477], color="red blue", name="Arrows_16.1790008545_6")

cluster_dict["16.1790008545"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-31.0), float(10.5), float(45.0), float(1.0)]

cluster_dict["16.1790008545_arrows"] += cgo_arrow([-31.0,10.5,45.0], [-33.454,10.054,46.148], color="red blue", name="Arrows_16.1790008545_7")

cluster_dict["16.1790008545"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-28.5), float(9.0), float(51.0), float(1.0)]

cluster_dict["16.1790008545_arrows"] += cgo_arrow([-28.5,9.0,51.0], [-29.403,7.109,52.477], color="red blue", name="Arrows_16.1790008545_8")

cluster_dict["16.1790008545"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-27.0), float(12.5), float(51.5), float(1.0)]

cluster_dict["16.1790008545_arrows"] += cgo_arrow([-27.0,12.5,51.5], [-27.105,12.341,55.227], color="red blue", name="Arrows_16.1790008545_9")

cmd.load_cgo(cluster_dict["16.1790008545"], "Features_16.1790008545", 1)
cmd.load_cgo(cluster_dict["16.1790008545_arrows"], "Arrows_16.1790008545")
cmd.set("transparency", 0.2,"Features_16.1790008545")
cmd.group("Pharmacophore_16.1790008545", members="Features_16.1790008545")
cmd.group("Pharmacophore_16.1790008545", members="Arrows_16.1790008545")

if dirpath:
    f = join(dirpath, "55/label_threshold_16.1790008545.mol2")
else:
    f = "55/label_threshold_16.1790008545.mol2"

cmd.load(f, 'label_threshold_16.1790008545')
cmd.hide('everything', 'label_threshold_16.1790008545')
cmd.label("label_threshold_16.1790008545", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.1790008545', members= 'label_threshold_16.1790008545')


if dirpath:
    f = join(dirpath, '55/mesh.grd')
else:
    f = '55/mesh.grd'
cmd.load(f, 'mesh_55')
cmd.isomesh("isomesh_55", "mesh_55", 0.9)
cmd.color("grey80", "isomesh_55")
cmd.set('transparency', 0.4, "isomesh_55")

cmd.group('hotspot_55', "isomesh_55")
cmd.group('hotspot_55', "mesh_55")

if dirpath:
    f = join(dirpath, "56/label_threshold_14.6.mol2")
else:
    f = "56/label_threshold_14.6.mol2"

cmd.load(f, 'label_threshold_14.6')
cmd.hide('everything', 'label_threshold_14.6')
cmd.label("label_threshold_14.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [14.6]
gfiles = ['56/donor.grd', '56/apolar.grd', '56/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 56
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


cluster_dict = {"16.0854997635":[], "16.0854997635_arrows":[]}

cluster_dict["16.0854997635"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(0.5), float(3.0), float(39.5), float(1.0)]

cluster_dict["16.0854997635_arrows"] += cgo_arrow([0.5,3.0,39.5], [0.089,2.153,36.754], color="blue red", name="Arrows_16.0854997635_1")

cluster_dict["16.0854997635"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(6.0), float(9.0), float(42.5), float(1.0)]

cluster_dict["16.0854997635_arrows"] += cgo_arrow([6.0,9.0,42.5], [3.715,8.369,43.968], color="blue red", name="Arrows_16.0854997635_2")

cluster_dict["16.0854997635"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(9.0), float(6.5), float(43.0), float(1.0)]

cluster_dict["16.0854997635_arrows"] += cgo_arrow([9.0,6.5,43.0], [10.755,4.658,42.464], color="blue red", name="Arrows_16.0854997635_3")

cluster_dict["16.0854997635"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(9.0), float(10.5), float(44.0), float(1.0)]

cluster_dict["16.0854997635_arrows"] += cgo_arrow([9.0,10.5,44.0], [8.374,12.84,46.167], color="blue red", name="Arrows_16.0854997635_4")

cluster_dict["16.0854997635"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(11.0), float(9.5), float(41.5), float(1.0)]

cluster_dict["16.0854997635_arrows"] += cgo_arrow([11.0,9.5,41.5], [11.409,11.758,39.462], color="blue red", name="Arrows_16.0854997635_5")

cluster_dict["16.0854997635"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(3.42527865445), float(10.3555159003), float(47.7502698465), float(1.0)]


cluster_dict["16.0854997635"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(7.25004312143), float(7.66316177769), float(43.9679580642), float(1.0)]


cluster_dict["16.0854997635"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(12.1621799459), float(11.2499862928), float(41.5851197973), float(1.0)]


cluster_dict["16.0854997635"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(4.0), float(11.5), float(46.0), float(1.0)]

cluster_dict["16.0854997635_arrows"] += cgo_arrow([4.0,11.5,46.0], [4.114,14.533,46.591], color="red blue", name="Arrows_16.0854997635_6")

cluster_dict["16.0854997635"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(5.5), float(8.0), float(43.0), float(1.0)]

cluster_dict["16.0854997635_arrows"] += cgo_arrow([5.5,8.0,43.0], [3.715,8.369,43.968], color="red blue", name="Arrows_16.0854997635_7")

cluster_dict["16.0854997635"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(8.5), float(10.0), float(39.5), float(1.0)]

cluster_dict["16.0854997635_arrows"] += cgo_arrow([8.5,10.0,39.5], [7.5,12.038,40.855], color="red blue", name="Arrows_16.0854997635_8")

cluster_dict["16.0854997635"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(9.5), float(6.0), float(50.0), float(1.0)]

cluster_dict["16.0854997635_arrows"] += cgo_arrow([9.5,6.0,50.0], [7.899,6.991,52.135], color="red blue", name="Arrows_16.0854997635_9")

cmd.load_cgo(cluster_dict["16.0854997635"], "Features_16.0854997635", 1)
cmd.load_cgo(cluster_dict["16.0854997635_arrows"], "Arrows_16.0854997635")
cmd.set("transparency", 0.2,"Features_16.0854997635")
cmd.group("Pharmacophore_16.0854997635", members="Features_16.0854997635")
cmd.group("Pharmacophore_16.0854997635", members="Arrows_16.0854997635")

if dirpath:
    f = join(dirpath, "56/label_threshold_16.0854997635.mol2")
else:
    f = "56/label_threshold_16.0854997635.mol2"

cmd.load(f, 'label_threshold_16.0854997635')
cmd.hide('everything', 'label_threshold_16.0854997635')
cmd.label("label_threshold_16.0854997635", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.0854997635', members= 'label_threshold_16.0854997635')


if dirpath:
    f = join(dirpath, '56/mesh.grd')
else:
    f = '56/mesh.grd'
cmd.load(f, 'mesh_56')
cmd.isomesh("isomesh_56", "mesh_56", 0.9)
cmd.color("grey80", "isomesh_56")
cmd.set('transparency', 0.4, "isomesh_56")

cmd.group('hotspot_56', "isomesh_56")
cmd.group('hotspot_56', "mesh_56")

if dirpath:
    f = join(dirpath, "57/label_threshold_7.5.mol2")
else:
    f = "57/label_threshold_7.5.mol2"

cmd.load(f, 'label_threshold_7.5')
cmd.hide('everything', 'label_threshold_7.5')
cmd.label("label_threshold_7.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [7.5]
gfiles = ['57/donor.grd', '57/apolar.grd', '57/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 57
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


cluster_dict = {"16.0419998169":[], "16.0419998169_arrows":[]}

cluster_dict["16.0419998169"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(45.3923891182), float(26.2478267365), float(38.7407470164), float(1.0)]


cmd.load_cgo(cluster_dict["16.0419998169"], "Features_16.0419998169", 1)
cmd.load_cgo(cluster_dict["16.0419998169_arrows"], "Arrows_16.0419998169")
cmd.set("transparency", 0.2,"Features_16.0419998169")
cmd.group("Pharmacophore_16.0419998169", members="Features_16.0419998169")
cmd.group("Pharmacophore_16.0419998169", members="Arrows_16.0419998169")

if dirpath:
    f = join(dirpath, "57/label_threshold_16.0419998169.mol2")
else:
    f = "57/label_threshold_16.0419998169.mol2"

cmd.load(f, 'label_threshold_16.0419998169')
cmd.hide('everything', 'label_threshold_16.0419998169')
cmd.label("label_threshold_16.0419998169", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.0419998169', members= 'label_threshold_16.0419998169')


if dirpath:
    f = join(dirpath, '57/mesh.grd')
else:
    f = '57/mesh.grd'
cmd.load(f, 'mesh_57')
cmd.isomesh("isomesh_57", "mesh_57", 0.9)
cmd.color("grey80", "isomesh_57")
cmd.set('transparency', 0.4, "isomesh_57")

cmd.group('hotspot_57', "isomesh_57")
cmd.group('hotspot_57', "mesh_57")

if dirpath:
    f = join(dirpath, "58/label_threshold_14.8.mol2")
else:
    f = "58/label_threshold_14.8.mol2"

cmd.load(f, 'label_threshold_14.8')
cmd.hide('everything', 'label_threshold_14.8')
cmd.label("label_threshold_14.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [14.8]
gfiles = ['58/donor.grd', '58/apolar.grd', '58/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 58
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


cluster_dict = {"16.0340003967":[], "16.0340003967_arrows":[]}

cluster_dict["16.0340003967"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(21.5), float(5.5), float(38.0), float(1.0)]

cluster_dict["16.0340003967_arrows"] += cgo_arrow([21.5,5.5,38.0], [22.648,5.212,35.929], color="blue red", name="Arrows_16.0340003967_1")

cluster_dict["16.0340003967"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(9.0), float(41.0), float(1.0)]

cluster_dict["16.0340003967_arrows"] += cgo_arrow([23.5,9.0,41.0], [24.082,12.057,42.288], color="blue red", name="Arrows_16.0340003967_2")

cluster_dict["16.0340003967"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(16.0), float(45.5), float(1.0)]

cluster_dict["16.0340003967_arrows"] += cgo_arrow([23.5,16.0,45.5], [25.971,15.865,43.79], color="blue red", name="Arrows_16.0340003967_3")

cluster_dict["16.0340003967"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(15.0), float(36.5), float(1.0)]

cluster_dict["16.0340003967_arrows"] += cgo_arrow([26.5,15.0,36.5], [28.472,16.776,37.727], color="blue red", name="Arrows_16.0340003967_4")

cluster_dict["16.0340003967"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.2206122954), float(5.46950127292), float(41.3765765547), float(1.0)]


cluster_dict["16.0340003967"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.7750251756), float(11.4969795874), float(34.0825699141), float(1.0)]


cluster_dict["16.0340003967"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(23.2748932444), float(16.1824154588), float(45.7678372386), float(1.0)]


cluster_dict["16.0340003967"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(19.0), float(-0.5), float(42.5), float(1.0)]

cluster_dict["16.0340003967_arrows"] += cgo_arrow([19.0,-0.5,42.5], [15.944,-0.291,41.945], color="red blue", name="Arrows_16.0340003967_5")

cluster_dict["16.0340003967"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(7.0), float(38.0), float(1.0)]

cluster_dict["16.0340003967_arrows"] += cgo_arrow([22.5,7.0,38.0], [22.648,5.212,35.929], color="red blue", name="Arrows_16.0340003967_6")

cluster_dict["16.0340003967"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(8.0), float(31.5), float(1.0)]

cluster_dict["16.0340003967_arrows"] += cgo_arrow([23.0,8.0,31.5], [24.327,5.898,29.592], color="red blue", name="Arrows_16.0340003967_7")

cluster_dict["16.0340003967"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(15.0), float(37.0), float(1.0)]

cluster_dict["16.0340003967_arrows"] += cgo_arrow([23.0,15.0,37.0], [20.949,16.359,38.69], color="red blue", name="Arrows_16.0340003967_8")

cluster_dict["16.0340003967"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(25.0), float(15.0), float(46.0), float(1.0)]

cluster_dict["16.0340003967_arrows"] += cgo_arrow([25.0,15.0,46.0], [25.971,15.865,43.79], color="red blue", name="Arrows_16.0340003967_9")

cmd.load_cgo(cluster_dict["16.0340003967"], "Features_16.0340003967", 1)
cmd.load_cgo(cluster_dict["16.0340003967_arrows"], "Arrows_16.0340003967")
cmd.set("transparency", 0.2,"Features_16.0340003967")
cmd.group("Pharmacophore_16.0340003967", members="Features_16.0340003967")
cmd.group("Pharmacophore_16.0340003967", members="Arrows_16.0340003967")

if dirpath:
    f = join(dirpath, "58/label_threshold_16.0340003967.mol2")
else:
    f = "58/label_threshold_16.0340003967.mol2"

cmd.load(f, 'label_threshold_16.0340003967')
cmd.hide('everything', 'label_threshold_16.0340003967')
cmd.label("label_threshold_16.0340003967", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_16.0340003967', members= 'label_threshold_16.0340003967')


if dirpath:
    f = join(dirpath, '58/mesh.grd')
else:
    f = '58/mesh.grd'
cmd.load(f, 'mesh_58')
cmd.isomesh("isomesh_58", "mesh_58", 0.9)
cmd.color("grey80", "isomesh_58")
cmd.set('transparency', 0.4, "isomesh_58")

cmd.group('hotspot_58', "isomesh_58")
cmd.group('hotspot_58', "mesh_58")

if dirpath:
    f = join(dirpath, "59/label_threshold_14.1.mol2")
else:
    f = "59/label_threshold_14.1.mol2"

cmd.load(f, 'label_threshold_14.1')
cmd.hide('everything', 'label_threshold_14.1')
cmd.label("label_threshold_14.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [14.1]
gfiles = ['59/donor.grd', '59/apolar.grd', '59/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 59
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


cluster_dict = {"15.9499998093":[], "15.9499998093_arrows":[]}

cluster_dict["15.9499998093"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(8.0), float(26.5), float(1.0)]

cluster_dict["15.9499998093_arrows"] += cgo_arrow([21.0,8.0,26.5], [20.263,11.401,25.376], color="blue red", name="Arrows_15.9499998093_1")

cluster_dict["15.9499998093"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(5.0), float(26.0), float(1.0)]

cluster_dict["15.9499998093_arrows"] += cgo_arrow([23.0,5.0,26.0], [21.848,2.901,23.953], color="blue red", name="Arrows_15.9499998093_2")

cluster_dict["15.9499998093"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(21.5), float(5.5), float(38.0), float(1.0)]

cluster_dict["15.9499998093_arrows"] += cgo_arrow([21.5,5.5,38.0], [22.648,5.212,35.929], color="blue red", name="Arrows_15.9499998093_3")

cluster_dict["15.9499998093"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.5), float(6.5), float(22.5), float(1.0)]

cluster_dict["15.9499998093_arrows"] += cgo_arrow([26.5,6.5,22.5], [24.383,5.1,20.493], color="blue red", name="Arrows_15.9499998093_4")

cluster_dict["15.9499998093"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(4.0), float(24.5), float(1.0)]

cluster_dict["15.9499998093_arrows"] += cgo_arrow([29.5,4.0,24.5], [31.281,4.479,26.704], color="blue red", name="Arrows_15.9499998093_5")

cluster_dict["15.9499998093"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(30.0), float(4.5), float(20.0), float(1.0)]

cluster_dict["15.9499998093_arrows"] += cgo_arrow([30.0,4.5,20.0], [31.989,4.998,17.834], color="blue red", name="Arrows_15.9499998093_6")

cluster_dict["15.9499998093"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.2254589387), float(3.39278000187), float(39.8753166486), float(1.0)]


cluster_dict["15.9499998093"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.3062343537), float(6.40181565209), float(24.9081882365), float(1.0)]


cluster_dict["15.9499998093"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.4986250672), float(9.14089035022), float(31.5037872321), float(1.0)]


cluster_dict["15.9499998093"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(21.5094222256), float(6.92887273134), float(36.8676098123), float(1.0)]


cluster_dict["15.9499998093"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(1.5), float(25.5), float(1.0)]


cluster_dict["15.9499998093"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(6.5), float(37.5), float(1.0)]

cluster_dict["15.9499998093_arrows"] += cgo_arrow([22.5,6.5,37.5], [22.648,5.212,35.929], color="red blue", name="Arrows_15.9499998093_7")

cluster_dict["15.9499998093"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(8.0), float(31.5), float(1.0)]

cluster_dict["15.9499998093_arrows"] += cgo_arrow([23.0,8.0,31.5], [24.327,5.898,29.592], color="red blue", name="Arrows_15.9499998093_8")

cluster_dict["15.9499998093"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(27.0), float(5.5), float(21.0), float(1.0)]

cluster_dict["15.9499998093_arrows"] += cgo_arrow([27.0,5.5,21.0], [24.383,5.1,20.493], color="red blue", name="Arrows_15.9499998093_9")

cmd.load_cgo(cluster_dict["15.9499998093"], "Features_15.9499998093", 1)
cmd.load_cgo(cluster_dict["15.9499998093_arrows"], "Arrows_15.9499998093")
cmd.set("transparency", 0.2,"Features_15.9499998093")
cmd.group("Pharmacophore_15.9499998093", members="Features_15.9499998093")
cmd.group("Pharmacophore_15.9499998093", members="Arrows_15.9499998093")

if dirpath:
    f = join(dirpath, "59/label_threshold_15.9499998093.mol2")
else:
    f = "59/label_threshold_15.9499998093.mol2"

cmd.load(f, 'label_threshold_15.9499998093')
cmd.hide('everything', 'label_threshold_15.9499998093')
cmd.label("label_threshold_15.9499998093", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.9499998093', members= 'label_threshold_15.9499998093')


if dirpath:
    f = join(dirpath, '59/mesh.grd')
else:
    f = '59/mesh.grd'
cmd.load(f, 'mesh_59')
cmd.isomesh("isomesh_59", "mesh_59", 0.9)
cmd.color("grey80", "isomesh_59")
cmd.set('transparency', 0.4, "isomesh_59")

cmd.group('hotspot_59', "isomesh_59")
cmd.group('hotspot_59', "mesh_59")

if dirpath:
    f = join(dirpath, "60/label_threshold_14.3.mol2")
else:
    f = "60/label_threshold_14.3.mol2"

cmd.load(f, 'label_threshold_14.3')
cmd.hide('everything', 'label_threshold_14.3')
cmd.label("label_threshold_14.3", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [14.3]
gfiles = ['60/donor.grd', '60/apolar.grd', '60/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 60
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


cluster_dict = {"15.906999588":[], "15.906999588_arrows":[]}

cluster_dict["15.906999588"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-33.5), float(9.5), float(52.0), float(1.0)]

cluster_dict["15.906999588_arrows"] += cgo_arrow([-33.5,9.5,52.0], [-33.337,6.587,52.943], color="blue red", name="Arrows_15.906999588_1")

cluster_dict["15.906999588"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-33.0), float(12.5), float(49.0), float(1.0)]

cluster_dict["15.906999588_arrows"] += cgo_arrow([-33.0,12.5,49.0], [-31.181,13.279,46.18], color="blue red", name="Arrows_15.906999588_2")

cluster_dict["15.906999588"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-31.0), float(11.5), float(52.0), float(1.0)]

cluster_dict["15.906999588_arrows"] += cgo_arrow([-31.0,11.5,52.0], [-29.015,12.348,54.074], color="blue red", name="Arrows_15.906999588_3")

cluster_dict["15.906999588"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-29.0), float(9.0), float(50.5), float(1.0)]

cluster_dict["15.906999588_arrows"] += cgo_arrow([-29.0,9.0,50.5], [-29.403,7.109,52.477], color="blue red", name="Arrows_15.906999588_4")

cluster_dict["15.906999588"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-29.7859250617), float(9.85005256231), float(47.7996637871), float(1.0)]


cluster_dict["15.906999588"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-34.833344282), float(7.83337712809), float(50.333311436), float(1.0)]


cluster_dict["15.906999588"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-32.5), float(12.0), float(48.0), float(1.0)]

cluster_dict["15.906999588_arrows"] += cgo_arrow([-32.5,12.0,48.0], [-31.181,13.279,46.18], color="red blue", name="Arrows_15.906999588_5")

cluster_dict["15.906999588"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-28.5), float(9.0), float(51.0), float(1.0)]

cluster_dict["15.906999588_arrows"] += cgo_arrow([-28.5,9.0,51.0], [-29.403,7.109,52.477], color="red blue", name="Arrows_15.906999588_6")

cluster_dict["15.906999588"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-31.0), float(10.5), float(45.0), float(1.0)]

cluster_dict["15.906999588_arrows"] += cgo_arrow([-31.0,10.5,45.0], [-33.454,10.054,46.148], color="red blue", name="Arrows_15.906999588_7")

cluster_dict["15.906999588"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-28.5), float(9.0), float(51.0), float(1.0)]

cluster_dict["15.906999588_arrows"] += cgo_arrow([-28.5,9.0,51.0], [-29.403,7.109,52.477], color="red blue", name="Arrows_15.906999588_8")

cluster_dict["15.906999588"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-27.0), float(12.5), float(51.5), float(1.0)]

cluster_dict["15.906999588_arrows"] += cgo_arrow([-27.0,12.5,51.5], [-27.105,12.341,55.227], color="red blue", name="Arrows_15.906999588_9")

cmd.load_cgo(cluster_dict["15.906999588"], "Features_15.906999588", 1)
cmd.load_cgo(cluster_dict["15.906999588_arrows"], "Arrows_15.906999588")
cmd.set("transparency", 0.2,"Features_15.906999588")
cmd.group("Pharmacophore_15.906999588", members="Features_15.906999588")
cmd.group("Pharmacophore_15.906999588", members="Arrows_15.906999588")

if dirpath:
    f = join(dirpath, "60/label_threshold_15.906999588.mol2")
else:
    f = "60/label_threshold_15.906999588.mol2"

cmd.load(f, 'label_threshold_15.906999588')
cmd.hide('everything', 'label_threshold_15.906999588')
cmd.label("label_threshold_15.906999588", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.906999588', members= 'label_threshold_15.906999588')


if dirpath:
    f = join(dirpath, '60/mesh.grd')
else:
    f = '60/mesh.grd'
cmd.load(f, 'mesh_60')
cmd.isomesh("isomesh_60", "mesh_60", 0.9)
cmd.color("grey80", "isomesh_60")
cmd.set('transparency', 0.4, "isomesh_60")

cmd.group('hotspot_60', "isomesh_60")
cmd.group('hotspot_60', "mesh_60")

if dirpath:
    f = join(dirpath, "61/label_threshold_13.0.mol2")
else:
    f = "61/label_threshold_13.0.mol2"

cmd.load(f, 'label_threshold_13.0')
cmd.hide('everything', 'label_threshold_13.0')
cmd.label("label_threshold_13.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.0]
gfiles = ['61/donor.grd', '61/apolar.grd', '61/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 61
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


cluster_dict = {"15.8819999695":[], "15.8819999695_arrows":[]}

cluster_dict["15.8819999695"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-19.0), float(-8.0), float(33.5), float(1.0)]

cluster_dict["15.8819999695_arrows"] += cgo_arrow([-19.0,-8.0,33.5], [-20.183,-8.769,31.779], color="blue red", name="Arrows_15.8819999695_1")

cluster_dict["15.8819999695"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-17.5), float(-10.5), float(28.0), float(1.0)]

cluster_dict["15.8819999695_arrows"] += cgo_arrow([-17.5,-10.5,28.0], [-19.572,-9.949,26.024], color="blue red", name="Arrows_15.8819999695_2")

cluster_dict["15.8819999695"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-15.5), float(-13.0), float(33.5), float(1.0)]

cluster_dict["15.8819999695_arrows"] += cgo_arrow([-15.5,-13.0,33.5], [-18.251,-13.39,33.288], color="blue red", name="Arrows_15.8819999695_3")

cluster_dict["15.8819999695"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-23.1857651346), float(-17.2200114681), float(37.0519769517), float(1.0)]


cluster_dict["15.8819999695"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-18.8844205294), float(-7.51376315249), float(37.0975887935), float(1.0)]


cluster_dict["15.8819999695"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-15.5298795438), float(-13.3566427357), float(31.427262498), float(1.0)]


cluster_dict["15.8819999695"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-12.2572565123), float(-13.9725546503), float(43.256809705), float(1.0)]


cluster_dict["15.8819999695"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-17.5), float(-8.0), float(36.5), float(1.0)]

cluster_dict["15.8819999695_arrows"] += cgo_arrow([-17.5,-8.0,36.5], [-14.978,-7.526,35.757], color="red blue", name="Arrows_15.8819999695_4")

cluster_dict["15.8819999695"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-16.5), float(-10.0), float(27.0), float(1.0)]

cluster_dict["15.8819999695_arrows"] += cgo_arrow([-16.5,-10.0,27.0], [-19.572,-9.949,26.024], color="red blue", name="Arrows_15.8819999695_5")

cluster_dict["15.8819999695"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-14.5), float(-14.0), float(43.5), float(1.0)]

cluster_dict["15.8819999695_arrows"] += cgo_arrow([-14.5,-14.0,43.5], [-15.04,-16.828,41.999], color="red blue", name="Arrows_15.8819999695_6")

cluster_dict["15.8819999695"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-15.5), float(-4.5), float(35.0), float(1.0)]

cluster_dict["15.8819999695_arrows"] += cgo_arrow([-15.5,-4.5,35.0], [-14.978,-7.526,35.757], color="red blue", name="Arrows_15.8819999695_7")

cluster_dict["15.8819999695"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-12.5), float(-13.0), float(42.5), float(1.0)]

cluster_dict["15.8819999695_arrows"] += cgo_arrow([-12.5,-13.0,42.5], [-13.853,-9.961,40.813], color="red blue", name="Arrows_15.8819999695_8")

cmd.load_cgo(cluster_dict["15.8819999695"], "Features_15.8819999695", 1)
cmd.load_cgo(cluster_dict["15.8819999695_arrows"], "Arrows_15.8819999695")
cmd.set("transparency", 0.2,"Features_15.8819999695")
cmd.group("Pharmacophore_15.8819999695", members="Features_15.8819999695")
cmd.group("Pharmacophore_15.8819999695", members="Arrows_15.8819999695")

if dirpath:
    f = join(dirpath, "61/label_threshold_15.8819999695.mol2")
else:
    f = "61/label_threshold_15.8819999695.mol2"

cmd.load(f, 'label_threshold_15.8819999695')
cmd.hide('everything', 'label_threshold_15.8819999695')
cmd.label("label_threshold_15.8819999695", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.8819999695', members= 'label_threshold_15.8819999695')


if dirpath:
    f = join(dirpath, '61/mesh.grd')
else:
    f = '61/mesh.grd'
cmd.load(f, 'mesh_61')
cmd.isomesh("isomesh_61", "mesh_61", 0.9)
cmd.color("grey80", "isomesh_61")
cmd.set('transparency', 0.4, "isomesh_61")

cmd.group('hotspot_61', "isomesh_61")
cmd.group('hotspot_61', "mesh_61")

if dirpath:
    f = join(dirpath, "62/label_threshold_10.4.mol2")
else:
    f = "62/label_threshold_10.4.mol2"

cmd.load(f, 'label_threshold_10.4')
cmd.hide('everything', 'label_threshold_10.4')
cmd.label("label_threshold_10.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.4]
gfiles = ['62/donor.grd', '62/apolar.grd', '62/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 62
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


cluster_dict = {"15.5719995499":[], "15.5719995499_arrows":[]}

cluster_dict["15.5719995499"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-70.5), float(14.5), float(69.0), float(1.0)]

cluster_dict["15.5719995499_arrows"] += cgo_arrow([-70.5,14.5,69.0], [-71.191,13.36,66.435], color="blue red", name="Arrows_15.5719995499_1")

cluster_dict["15.5719995499"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-69.3048645571), float(13.6232113038), float(73.3737298246), float(1.0)]


cluster_dict["15.5719995499"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-62.5418070196), float(9.25298450452), float(77.4937774538), float(1.0)]


cluster_dict["15.5719995499"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-67.5), float(12.5), float(72.0), float(1.0)]

cluster_dict["15.5719995499_arrows"] += cgo_arrow([-67.5,12.5,72.0], [-67.149,9.709,73.764], color="red blue", name="Arrows_15.5719995499_2")

cmd.load_cgo(cluster_dict["15.5719995499"], "Features_15.5719995499", 1)
cmd.load_cgo(cluster_dict["15.5719995499_arrows"], "Arrows_15.5719995499")
cmd.set("transparency", 0.2,"Features_15.5719995499")
cmd.group("Pharmacophore_15.5719995499", members="Features_15.5719995499")
cmd.group("Pharmacophore_15.5719995499", members="Arrows_15.5719995499")

if dirpath:
    f = join(dirpath, "62/label_threshold_15.5719995499.mol2")
else:
    f = "62/label_threshold_15.5719995499.mol2"

cmd.load(f, 'label_threshold_15.5719995499')
cmd.hide('everything', 'label_threshold_15.5719995499')
cmd.label("label_threshold_15.5719995499", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.5719995499', members= 'label_threshold_15.5719995499')


if dirpath:
    f = join(dirpath, '62/mesh.grd')
else:
    f = '62/mesh.grd'
cmd.load(f, 'mesh_62')
cmd.isomesh("isomesh_62", "mesh_62", 0.9)
cmd.color("grey80", "isomesh_62")
cmd.set('transparency', 0.4, "isomesh_62")

cmd.group('hotspot_62', "isomesh_62")
cmd.group('hotspot_62', "mesh_62")

if dirpath:
    f = join(dirpath, "63/label_threshold_13.8.mol2")
else:
    f = "63/label_threshold_13.8.mol2"

cmd.load(f, 'label_threshold_13.8')
cmd.hide('everything', 'label_threshold_13.8')
cmd.label("label_threshold_13.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.8]
gfiles = ['63/donor.grd', '63/apolar.grd', '63/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 63
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


cluster_dict = {"15.4250001907":[], "15.4250001907_arrows":[]}

cluster_dict["15.4250001907"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(15.0), float(13.0), float(1.0)]

cluster_dict["15.4250001907_arrows"] += cgo_arrow([22.0,15.0,13.0], [22.028,12.401,14.823], color="blue red", name="Arrows_15.4250001907_1")

cluster_dict["15.4250001907"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(13.5), float(20.5), float(1.0)]

cluster_dict["15.4250001907_arrows"] += cgo_arrow([22.5,13.5,20.5], [23.814,16.273,20.013], color="blue red", name="Arrows_15.4250001907_2")

cluster_dict["15.4250001907"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(10.5), float(24.0), float(1.0)]

cluster_dict["15.4250001907_arrows"] += cgo_arrow([26.0,10.5,24.0], [24.828,13.48,24.46], color="blue red", name="Arrows_15.4250001907_3")

cluster_dict["15.4250001907"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(9.5), float(16.5), float(1.0)]

cluster_dict["15.4250001907_arrows"] += cgo_arrow([28.0,9.5,16.5], [25.678,8.126,17.598], color="blue red", name="Arrows_15.4250001907_4")

cluster_dict["15.4250001907"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(29.0), float(10.0), float(27.5), float(1.0)]

cluster_dict["15.4250001907_arrows"] += cgo_arrow([29.0,10.0,27.5], [27.929,12.58,27.784], color="blue red", name="Arrows_15.4250001907_5")

cluster_dict["15.4250001907"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(32.5), float(14.5), float(24.5), float(1.0)]

cluster_dict["15.4250001907_arrows"] += cgo_arrow([32.5,14.5,24.5], [32.987,16.799,26.549], color="blue red", name="Arrows_15.4250001907_6")

cluster_dict["15.4250001907"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.8976240499), float(26.0754316458), float(23.3001049547), float(1.0)]


cluster_dict["15.4250001907"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(24.9), float(26.7), float(1.0)]


cluster_dict["15.4250001907"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(23.4435881104), float(15.6366603509), float(12.5624097365), float(1.0)]


cluster_dict["15.4250001907"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.8543112464), float(10.6931334121), float(20.6761831407), float(1.0)]


cluster_dict["15.4250001907"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(23.5), float(25.0), float(1.0)]

cluster_dict["15.4250001907_arrows"] += cgo_arrow([18.0,23.5,25.0], [20.807,22.261,23.987], color="red blue", name="Arrows_15.4250001907_7")

cluster_dict["15.4250001907"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(23.0), float(22.5), float(1.0)]

cluster_dict["15.4250001907_arrows"] += cgo_arrow([22.0,23.0,22.5], [20.807,22.261,23.987], color="red blue", name="Arrows_15.4250001907_8")

cluster_dict["15.4250001907"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(9.0), float(26.0), float(1.0)]

cluster_dict["15.4250001907_arrows"] += cgo_arrow([30.5,9.0,26.0], [34.056,6.45,27.561], color="red blue", name="Arrows_15.4250001907_9")

cmd.load_cgo(cluster_dict["15.4250001907"], "Features_15.4250001907", 1)
cmd.load_cgo(cluster_dict["15.4250001907_arrows"], "Arrows_15.4250001907")
cmd.set("transparency", 0.2,"Features_15.4250001907")
cmd.group("Pharmacophore_15.4250001907", members="Features_15.4250001907")
cmd.group("Pharmacophore_15.4250001907", members="Arrows_15.4250001907")

if dirpath:
    f = join(dirpath, "63/label_threshold_15.4250001907.mol2")
else:
    f = "63/label_threshold_15.4250001907.mol2"

cmd.load(f, 'label_threshold_15.4250001907')
cmd.hide('everything', 'label_threshold_15.4250001907')
cmd.label("label_threshold_15.4250001907", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.4250001907', members= 'label_threshold_15.4250001907')


if dirpath:
    f = join(dirpath, '63/mesh.grd')
else:
    f = '63/mesh.grd'
cmd.load(f, 'mesh_63')
cmd.isomesh("isomesh_63", "mesh_63", 0.9)
cmd.color("grey80", "isomesh_63")
cmd.set('transparency', 0.4, "isomesh_63")

cmd.group('hotspot_63', "isomesh_63")
cmd.group('hotspot_63', "mesh_63")

if dirpath:
    f = join(dirpath, "64/label_threshold_13.5.mol2")
else:
    f = "64/label_threshold_13.5.mol2"

cmd.load(f, 'label_threshold_13.5')
cmd.hide('everything', 'label_threshold_13.5')
cmd.label("label_threshold_13.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.5]
gfiles = ['64/donor.grd', '64/apolar.grd', '64/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 64
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


cluster_dict = {"15.4040002823":[], "15.4040002823_arrows":[]}

cluster_dict["15.4040002823"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-19.0), float(-8.0), float(33.5), float(1.0)]

cluster_dict["15.4040002823_arrows"] += cgo_arrow([-19.0,-8.0,33.5], [-20.183,-8.769,31.779], color="blue red", name="Arrows_15.4040002823_1")

cluster_dict["15.4040002823"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-17.5), float(-10.5), float(28.0), float(1.0)]

cluster_dict["15.4040002823_arrows"] += cgo_arrow([-17.5,-10.5,28.0], [-19.572,-9.949,26.024], color="blue red", name="Arrows_15.4040002823_2")

cluster_dict["15.4040002823"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-15.5), float(-13.0), float(33.5), float(1.0)]

cluster_dict["15.4040002823_arrows"] += cgo_arrow([-15.5,-13.0,33.5], [-18.251,-13.39,33.288], color="blue red", name="Arrows_15.4040002823_3")

cluster_dict["15.4040002823"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-14.5), float(-4.5), float(46.0), float(1.0)]

cluster_dict["15.4040002823_arrows"] += cgo_arrow([-14.5,-4.5,46.0], [-12.281,-6.363,46.478], color="blue red", name="Arrows_15.4040002823_4")

cluster_dict["15.4040002823"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-23.2625743437), float(-16.6871735313), float(37.0445950279), float(1.0)]


cluster_dict["15.4040002823"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-19.2210815991), float(-7.05031784267), float(45.0524312044), float(1.0)]


cluster_dict["15.4040002823"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-19.4304403912), float(-7.91517236404), float(36.5054969484), float(1.0)]


cluster_dict["15.4040002823"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-18.3357520372), float(-15.5800382872), float(45.595405281), float(1.0)]


cluster_dict["15.4040002823"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-15.434186382), float(-13.4233635096), float(31.5802422629), float(1.0)]


cluster_dict["15.4040002823"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-13.0478973301), float(-12.5695111455), float(44.2829302378), float(1.0)]


cluster_dict["15.4040002823"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-21.0), float(-6.0), float(45.5), float(1.0)]

cluster_dict["15.4040002823_arrows"] += cgo_arrow([-21.0,-6.0,45.5], [-22.337,-3.396,48.966], color="red blue", name="Arrows_15.4040002823_5")

cluster_dict["15.4040002823"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-17.5), float(-8.0), float(36.5), float(1.0)]

cluster_dict["15.4040002823_arrows"] += cgo_arrow([-17.5,-8.0,36.5], [-14.978,-7.526,35.757], color="red blue", name="Arrows_15.4040002823_6")

cluster_dict["15.4040002823"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-16.5), float(-10.0), float(27.0), float(1.0)]

cluster_dict["15.4040002823_arrows"] += cgo_arrow([-16.5,-10.0,27.0], [-19.572,-9.949,26.024], color="red blue", name="Arrows_15.4040002823_7")

cluster_dict["15.4040002823"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-15.5), float(-4.5), float(35.0), float(1.0)]

cluster_dict["15.4040002823_arrows"] += cgo_arrow([-15.5,-4.5,35.0], [-14.978,-7.526,35.757], color="red blue", name="Arrows_15.4040002823_8")

cluster_dict["15.4040002823"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-12.5), float(-13.0), float(42.5), float(1.0)]

cluster_dict["15.4040002823_arrows"] += cgo_arrow([-12.5,-13.0,42.5], [-13.853,-9.961,40.813], color="red blue", name="Arrows_15.4040002823_9")

cmd.load_cgo(cluster_dict["15.4040002823"], "Features_15.4040002823", 1)
cmd.load_cgo(cluster_dict["15.4040002823_arrows"], "Arrows_15.4040002823")
cmd.set("transparency", 0.2,"Features_15.4040002823")
cmd.group("Pharmacophore_15.4040002823", members="Features_15.4040002823")
cmd.group("Pharmacophore_15.4040002823", members="Arrows_15.4040002823")

if dirpath:
    f = join(dirpath, "64/label_threshold_15.4040002823.mol2")
else:
    f = "64/label_threshold_15.4040002823.mol2"

cmd.load(f, 'label_threshold_15.4040002823')
cmd.hide('everything', 'label_threshold_15.4040002823')
cmd.label("label_threshold_15.4040002823", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.4040002823', members= 'label_threshold_15.4040002823')


if dirpath:
    f = join(dirpath, '64/mesh.grd')
else:
    f = '64/mesh.grd'
cmd.load(f, 'mesh_64')
cmd.isomesh("isomesh_64", "mesh_64", 0.9)
cmd.color("grey80", "isomesh_64")
cmd.set('transparency', 0.4, "isomesh_64")

cmd.group('hotspot_64', "isomesh_64")
cmd.group('hotspot_64', "mesh_64")

if dirpath:
    f = join(dirpath, "65/label_threshold_13.8.mol2")
else:
    f = "65/label_threshold_13.8.mol2"

cmd.load(f, 'label_threshold_13.8')
cmd.hide('everything', 'label_threshold_13.8')
cmd.label("label_threshold_13.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.8]
gfiles = ['65/donor.grd', '65/apolar.grd', '65/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 65
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


cluster_dict = {"15.2840003967":[], "15.2840003967_arrows":[]}

cluster_dict["15.2840003967"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-19.0), float(-8.0), float(33.5), float(1.0)]

cluster_dict["15.2840003967_arrows"] += cgo_arrow([-19.0,-8.0,33.5], [-20.183,-8.769,31.779], color="blue red", name="Arrows_15.2840003967_1")

cluster_dict["15.2840003967"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-16.0), float(-13.0), float(33.5), float(1.0)]

cluster_dict["15.2840003967_arrows"] += cgo_arrow([-16.0,-13.0,33.5], [-18.251,-13.39,33.288], color="blue red", name="Arrows_15.2840003967_2")

cluster_dict["15.2840003967"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-16.0), float(-3.0), float(54.0), float(1.0)]

cluster_dict["15.2840003967_arrows"] += cgo_arrow([-16.0,-3.0,54.0], [-17.363,-0.587,52.277], color="blue red", name="Arrows_15.2840003967_3")

cluster_dict["15.2840003967"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-16.0), float(-4.5), float(48.5), float(1.0)]

cluster_dict["15.2840003967_arrows"] += cgo_arrow([-16.0,-4.5,48.5], [-18.731,-2.721,49.33], color="blue red", name="Arrows_15.2840003967_4")

cluster_dict["15.2840003967"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-14.0), float(-4.0), float(46.5), float(1.0)]

cluster_dict["15.2840003967_arrows"] += cgo_arrow([-14.0,-4.0,46.5], [-13.027,-1.485,45.573], color="blue red", name="Arrows_15.2840003967_5")

cluster_dict["15.2840003967"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-12.0), float(-13.0), float(48.0), float(1.0)]

cluster_dict["15.2840003967_arrows"] += cgo_arrow([-12.0,-13.0,48.0], [-9.173,-12.246,47.746], color="blue red", name="Arrows_15.2840003967_6")

cluster_dict["15.2840003967"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-22.9024657152), float(-16.1591007107), float(38.0320013085), float(1.0)]


cluster_dict["15.2840003967"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-19.5291297208), float(-7.74296762811), float(36.6158643012), float(1.0)]


cluster_dict["15.2840003967"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-19.0490797546), float(-6.86503067485), float(44.6993865031), float(1.0)]


cluster_dict["15.2840003967"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-17.9671021997), float(-13.0626775495), float(55.0400113954), float(1.0)]


cluster_dict["15.2840003967"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-18.1102794022), float(-16.3819907902), float(45.9523278989), float(1.0)]


cluster_dict["15.2840003967"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-16.3508429422), float(-3.55202164871), float(55.8010868519), float(1.0)]


cluster_dict["15.2840003967"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-14.9906902375), float(-11.6074083108), float(34.6482648242), float(1.0)]


cluster_dict["15.2840003967"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-14.984883501), float(-0.424593773066), float(49.3412767608), float(1.0)]


cluster_dict["15.2840003967"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-12.6250515948), float(-14.3372037449), float(46.0130440081), float(1.0)]


cluster_dict["15.2840003967"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-19.0), float(-4.0), float(56.5), float(1.0)]

cluster_dict["15.2840003967_arrows"] += cgo_arrow([-19.0,-4.0,56.5], [-21.532,-6.13,56.958], color="red blue", name="Arrows_15.2840003967_7")

cluster_dict["15.2840003967"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-17.0), float(-3.5), float(55.5), float(1.0)]


cluster_dict["15.2840003967"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-15.5), float(-4.5), float(35.0), float(1.0)]

cluster_dict["15.2840003967_arrows"] += cgo_arrow([-15.5,-4.5,35.0], [-14.978,-7.526,35.757], color="red blue", name="Arrows_15.2840003967_8")

cmd.load_cgo(cluster_dict["15.2840003967"], "Features_15.2840003967", 1)
cmd.load_cgo(cluster_dict["15.2840003967_arrows"], "Arrows_15.2840003967")
cmd.set("transparency", 0.2,"Features_15.2840003967")
cmd.group("Pharmacophore_15.2840003967", members="Features_15.2840003967")
cmd.group("Pharmacophore_15.2840003967", members="Arrows_15.2840003967")

if dirpath:
    f = join(dirpath, "65/label_threshold_15.2840003967.mol2")
else:
    f = "65/label_threshold_15.2840003967.mol2"

cmd.load(f, 'label_threshold_15.2840003967')
cmd.hide('everything', 'label_threshold_15.2840003967')
cmd.label("label_threshold_15.2840003967", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.2840003967', members= 'label_threshold_15.2840003967')


if dirpath:
    f = join(dirpath, '65/mesh.grd')
else:
    f = '65/mesh.grd'
cmd.load(f, 'mesh_65')
cmd.isomesh("isomesh_65", "mesh_65", 0.9)
cmd.color("grey80", "isomesh_65")
cmd.set('transparency', 0.4, "isomesh_65")

cmd.group('hotspot_65', "isomesh_65")
cmd.group('hotspot_65', "mesh_65")

if dirpath:
    f = join(dirpath, "66/label_threshold_12.7.mol2")
else:
    f = "66/label_threshold_12.7.mol2"

cmd.load(f, 'label_threshold_12.7')
cmd.hide('everything', 'label_threshold_12.7')
cmd.label("label_threshold_12.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.7]
gfiles = ['66/donor.grd', '66/apolar.grd', '66/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 66
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


cluster_dict = {"14.9829998016":[], "14.9829998016_arrows":[]}

cluster_dict["14.9829998016"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(18.0), float(9.0), float(6.5), float(1.0)]

cluster_dict["14.9829998016_arrows"] += cgo_arrow([18.0,9.0,6.5], [14.787,8.784,6.131], color="blue red", name="Arrows_14.9829998016_1")

cluster_dict["14.9829998016"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(7.5), float(10.0), float(1.0)]

cluster_dict["14.9829998016_arrows"] += cgo_arrow([21.0,7.5,10.0], [23.733,6.22,10.733], color="blue red", name="Arrows_14.9829998016_2")

cluster_dict["14.9829998016"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(21.5), float(12.0), float(12.5), float(1.0)]

cluster_dict["14.9829998016_arrows"] += cgo_arrow([21.5,12.0,12.5], [22.028,12.401,14.823], color="blue red", name="Arrows_14.9829998016_3")

cluster_dict["14.9829998016"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(13.5), float(20.5), float(1.0)]

cluster_dict["14.9829998016_arrows"] += cgo_arrow([22.5,13.5,20.5], [23.814,16.273,20.013], color="blue red", name="Arrows_14.9829998016_4")

cluster_dict["14.9829998016"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.5), float(7.5), float(15.5), float(1.0)]

cluster_dict["14.9829998016_arrows"] += cgo_arrow([27.5,7.5,15.5], [25.678,8.126,17.598], color="blue red", name="Arrows_14.9829998016_5")

cluster_dict["14.9829998016"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(20.4914244746), float(12.1937267005), float(10.1115693938), float(1.0)]


cluster_dict["14.9829998016"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(25.9452667449), float(11.7863555349), float(19.1131130858), float(1.0)]


cluster_dict["14.9829998016"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(19.5), float(14.5), float(12.5), float(1.0)]

cluster_dict["14.9829998016_arrows"] += cgo_arrow([19.5,14.5,12.5], [22.028,12.401,14.823], color="red blue", name="Arrows_14.9829998016_6")

cluster_dict["14.9829998016"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(21.0), float(8.0), float(7.5), float(1.0)]

cluster_dict["14.9829998016_arrows"] += cgo_arrow([21.0,8.0,7.5], [22.516,5.137,7.048], color="red blue", name="Arrows_14.9829998016_7")

cluster_dict["14.9829998016"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(10.0), float(21.0), float(1.0)]

cluster_dict["14.9829998016_arrows"] += cgo_arrow([22.0,10.0,21.0], [22.403,9.214,18.681], color="red blue", name="Arrows_14.9829998016_8")

cmd.load_cgo(cluster_dict["14.9829998016"], "Features_14.9829998016", 1)
cmd.load_cgo(cluster_dict["14.9829998016_arrows"], "Arrows_14.9829998016")
cmd.set("transparency", 0.2,"Features_14.9829998016")
cmd.group("Pharmacophore_14.9829998016", members="Features_14.9829998016")
cmd.group("Pharmacophore_14.9829998016", members="Arrows_14.9829998016")

if dirpath:
    f = join(dirpath, "66/label_threshold_14.9829998016.mol2")
else:
    f = "66/label_threshold_14.9829998016.mol2"

cmd.load(f, 'label_threshold_14.9829998016')
cmd.hide('everything', 'label_threshold_14.9829998016')
cmd.label("label_threshold_14.9829998016", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.9829998016', members= 'label_threshold_14.9829998016')


if dirpath:
    f = join(dirpath, '66/mesh.grd')
else:
    f = '66/mesh.grd'
cmd.load(f, 'mesh_66')
cmd.isomesh("isomesh_66", "mesh_66", 0.9)
cmd.color("grey80", "isomesh_66")
cmd.set('transparency', 0.4, "isomesh_66")

cmd.group('hotspot_66', "isomesh_66")
cmd.group('hotspot_66', "mesh_66")

if dirpath:
    f = join(dirpath, "67/label_threshold_12.6.mol2")
else:
    f = "67/label_threshold_12.6.mol2"

cmd.load(f, 'label_threshold_12.6')
cmd.hide('everything', 'label_threshold_12.6')
cmd.label("label_threshold_12.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.6]
gfiles = ['67/donor.grd', '67/apolar.grd', '67/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 67
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


cluster_dict = {"14.7410001755":[], "14.7410001755_arrows":[]}

cluster_dict["14.7410001755"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-54.5), float(-11.0), float(52.5), float(1.0)]

cluster_dict["14.7410001755_arrows"] += cgo_arrow([-54.5,-11.0,52.5], [-52.315,-11.683,50.879], color="blue red", name="Arrows_14.7410001755_1")

cluster_dict["14.7410001755"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-50.0), float(-17.5), float(53.5), float(1.0)]

cluster_dict["14.7410001755_arrows"] += cgo_arrow([-50.0,-17.5,53.5], [-48.745,-15.756,55.271], color="blue red", name="Arrows_14.7410001755_2")

cluster_dict["14.7410001755"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-49.5), float(-18.5), float(49.5), float(1.0)]

cluster_dict["14.7410001755_arrows"] += cgo_arrow([-49.5,-18.5,49.5], [-52.311,-18.251,48.129], color="blue red", name="Arrows_14.7410001755_3")

cluster_dict["14.7410001755"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-48.0), float(-18.5), float(56.0), float(1.0)]

cluster_dict["14.7410001755_arrows"] += cgo_arrow([-48.0,-18.5,56.0], [-48.745,-15.756,55.271], color="blue red", name="Arrows_14.7410001755_4")

cluster_dict["14.7410001755"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-55.6991480055), float(-12.9561640993), float(52.7671634441), float(1.0)]


cluster_dict["14.7410001755"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-48.2057900668), float(-20.514560702), float(52.5696291034), float(1.0)]


cluster_dict["14.7410001755"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-49.8642123418), float(-10.2142181319), float(54.2067451282), float(1.0)]


cluster_dict["14.7410001755"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-56.0), float(-12.0), float(48.5), float(1.0)]

cluster_dict["14.7410001755_arrows"] += cgo_arrow([-56.0,-12.0,48.5], [-58.036,-13.564,47.78], color="red blue", name="Arrows_14.7410001755_5")

cluster_dict["14.7410001755"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-52.5), float(-16.5), float(52.5), float(1.0)]

cluster_dict["14.7410001755_arrows"] += cgo_arrow([-52.5,-16.5,52.5], [-55.225,-17.308,51.745], color="red blue", name="Arrows_14.7410001755_6")

cluster_dict["14.7410001755"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-50.0), float(-22.5), float(49.0), float(1.0)]

cluster_dict["14.7410001755_arrows"] += cgo_arrow([-50.0,-22.5,49.0], [-52.292,-22.084,46.561], color="red blue", name="Arrows_14.7410001755_7")

cluster_dict["14.7410001755"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-45.5), float(-24.0), float(52.5), float(1.0)]

cluster_dict["14.7410001755_arrows"] += cgo_arrow([-45.5,-24.0,52.5], [-46.667,-24.409,55.215], color="red blue", name="Arrows_14.7410001755_8")

cmd.load_cgo(cluster_dict["14.7410001755"], "Features_14.7410001755", 1)
cmd.load_cgo(cluster_dict["14.7410001755_arrows"], "Arrows_14.7410001755")
cmd.set("transparency", 0.2,"Features_14.7410001755")
cmd.group("Pharmacophore_14.7410001755", members="Features_14.7410001755")
cmd.group("Pharmacophore_14.7410001755", members="Arrows_14.7410001755")

if dirpath:
    f = join(dirpath, "67/label_threshold_14.7410001755.mol2")
else:
    f = "67/label_threshold_14.7410001755.mol2"

cmd.load(f, 'label_threshold_14.7410001755')
cmd.hide('everything', 'label_threshold_14.7410001755')
cmd.label("label_threshold_14.7410001755", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.7410001755', members= 'label_threshold_14.7410001755')


if dirpath:
    f = join(dirpath, '67/mesh.grd')
else:
    f = '67/mesh.grd'
cmd.load(f, 'mesh_67')
cmd.isomesh("isomesh_67", "mesh_67", 0.9)
cmd.color("grey80", "isomesh_67")
cmd.set('transparency', 0.4, "isomesh_67")

cmd.group('hotspot_67', "isomesh_67")
cmd.group('hotspot_67', "mesh_67")

if dirpath:
    f = join(dirpath, "68/label_threshold_13.5.mol2")
else:
    f = "68/label_threshold_13.5.mol2"

cmd.load(f, 'label_threshold_13.5')
cmd.hide('everything', 'label_threshold_13.5')
cmd.label("label_threshold_13.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.5]
gfiles = ['68/donor.grd', '68/apolar.grd', '68/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 68
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


cluster_dict = {"14.6719999313":[], "14.6719999313_arrows":[]}

cluster_dict["14.6719999313"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-16.0), float(-3.0), float(54.0), float(1.0)]

cluster_dict["14.6719999313_arrows"] += cgo_arrow([-16.0,-3.0,54.0], [-17.363,-0.587,52.277], color="blue red", name="Arrows_14.6719999313_1")

cluster_dict["14.6719999313"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-16.0), float(-4.5), float(48.5), float(1.0)]

cluster_dict["14.6719999313_arrows"] += cgo_arrow([-16.0,-4.5,48.5], [-18.731,-2.721,49.33], color="blue red", name="Arrows_14.6719999313_2")

cluster_dict["14.6719999313"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-14.0), float(-4.0), float(46.5), float(1.0)]

cluster_dict["14.6719999313_arrows"] += cgo_arrow([-14.0,-4.0,46.5], [-13.027,-1.485,45.573], color="blue red", name="Arrows_14.6719999313_3")

cluster_dict["14.6719999313"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-12.0), float(-13.0), float(48.0), float(1.0)]

cluster_dict["14.6719999313_arrows"] += cgo_arrow([-12.0,-13.0,48.0], [-9.173,-12.246,47.746], color="blue red", name="Arrows_14.6719999313_4")

cluster_dict["14.6719999313"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-19.0710749899), float(-7.09927869789), float(45.0060518138), float(1.0)]


cluster_dict["14.6719999313"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-18.1608580596), float(-13.9680461792), float(56.3264427109), float(1.0)]


cluster_dict["14.6719999313"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-18.0925404926), float(-16.3691893375), float(45.9427934622), float(1.0)]


cluster_dict["14.6719999313"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-16.2796885207), float(-2.4488654791), float(56.9594561523), float(1.0)]


cluster_dict["14.6719999313"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-15.030006436), float(-0.833798249699), float(49.3738831333), float(1.0)]


cluster_dict["14.6719999313"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-12.6555050534), float(-15.0718888402), float(46.6962526756), float(1.0)]


cluster_dict["14.6719999313"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-21.0), float(-6.0), float(45.5), float(1.0)]

cluster_dict["14.6719999313_arrows"] += cgo_arrow([-21.0,-6.0,45.5], [-22.337,-3.396,48.966], color="red blue", name="Arrows_14.6719999313_5")

cluster_dict["14.6719999313"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-19.0), float(-4.0), float(56.5), float(1.0)]

cluster_dict["14.6719999313_arrows"] += cgo_arrow([-19.0,-4.0,56.5], [-21.532,-6.13,56.958], color="red blue", name="Arrows_14.6719999313_6")

cluster_dict["14.6719999313"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-18.0), float(-2.0), float(57.5), float(1.0)]


cluster_dict["14.6719999313"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-17.0), float(-3.5), float(55.5), float(1.0)]


cluster_dict["14.6719999313"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-11.0), float(-16.5), float(44.5), float(1.0)]

cluster_dict["14.6719999313_arrows"] += cgo_arrow([-11.0,-16.5,44.5], [-15.04,-16.828,41.999], color="red blue", name="Arrows_14.6719999313_7")

cmd.load_cgo(cluster_dict["14.6719999313"], "Features_14.6719999313", 1)
cmd.load_cgo(cluster_dict["14.6719999313_arrows"], "Arrows_14.6719999313")
cmd.set("transparency", 0.2,"Features_14.6719999313")
cmd.group("Pharmacophore_14.6719999313", members="Features_14.6719999313")
cmd.group("Pharmacophore_14.6719999313", members="Arrows_14.6719999313")

if dirpath:
    f = join(dirpath, "68/label_threshold_14.6719999313.mol2")
else:
    f = "68/label_threshold_14.6719999313.mol2"

cmd.load(f, 'label_threshold_14.6719999313')
cmd.hide('everything', 'label_threshold_14.6719999313')
cmd.label("label_threshold_14.6719999313", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.6719999313', members= 'label_threshold_14.6719999313')


if dirpath:
    f = join(dirpath, '68/mesh.grd')
else:
    f = '68/mesh.grd'
cmd.load(f, 'mesh_68')
cmd.isomesh("isomesh_68", "mesh_68", 0.9)
cmd.color("grey80", "isomesh_68")
cmd.set('transparency', 0.4, "isomesh_68")

cmd.group('hotspot_68', "isomesh_68")
cmd.group('hotspot_68', "mesh_68")

if dirpath:
    f = join(dirpath, "69/label_threshold_2.7.mol2")
else:
    f = "69/label_threshold_2.7.mol2"

cmd.load(f, 'label_threshold_2.7')
cmd.hide('everything', 'label_threshold_2.7')
cmd.label("label_threshold_2.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [2.7]
gfiles = ['69/donor.grd', '69/apolar.grd', '69/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 69
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


cluster_dict = {"14.6440000534":[], "14.6440000534_arrows":[]}

cluster_dict["14.6440000534"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(25.9974980177), float(-13.5446274677), float(48.521400508), float(1.0)]


cmd.load_cgo(cluster_dict["14.6440000534"], "Features_14.6440000534", 1)
cmd.load_cgo(cluster_dict["14.6440000534_arrows"], "Arrows_14.6440000534")
cmd.set("transparency", 0.2,"Features_14.6440000534")
cmd.group("Pharmacophore_14.6440000534", members="Features_14.6440000534")
cmd.group("Pharmacophore_14.6440000534", members="Arrows_14.6440000534")

if dirpath:
    f = join(dirpath, "69/label_threshold_14.6440000534.mol2")
else:
    f = "69/label_threshold_14.6440000534.mol2"

cmd.load(f, 'label_threshold_14.6440000534')
cmd.hide('everything', 'label_threshold_14.6440000534')
cmd.label("label_threshold_14.6440000534", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.6440000534', members= 'label_threshold_14.6440000534')


if dirpath:
    f = join(dirpath, '69/mesh.grd')
else:
    f = '69/mesh.grd'
cmd.load(f, 'mesh_69')
cmd.isomesh("isomesh_69", "mesh_69", 0.9)
cmd.color("grey80", "isomesh_69")
cmd.set('transparency', 0.4, "isomesh_69")

cmd.group('hotspot_69', "isomesh_69")
cmd.group('hotspot_69', "mesh_69")

if dirpath:
    f = join(dirpath, "70/label_threshold_13.4.mol2")
else:
    f = "70/label_threshold_13.4.mol2"

cmd.load(f, 'label_threshold_13.4')
cmd.hide('everything', 'label_threshold_13.4')
cmd.label("label_threshold_13.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.4]
gfiles = ['70/donor.grd', '70/apolar.grd', '70/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 70
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


cluster_dict = {"14.6140003204":[], "14.6140003204_arrows":[]}

cluster_dict["14.6140003204"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-13.5), float(13.0), float(24.5), float(1.0)]

cluster_dict["14.6140003204_arrows"] += cgo_arrow([-13.5,13.0,24.5], [-16.326,13.302,26.053], color="blue red", name="Arrows_14.6140003204_1")

cluster_dict["14.6140003204"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-11.0), float(8.0), float(30.5), float(1.0)]

cluster_dict["14.6140003204_arrows"] += cgo_arrow([-11.0,8.0,30.5], [-12.682,10.491,30.534], color="blue red", name="Arrows_14.6140003204_2")

cluster_dict["14.6140003204"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-11.0), float(11.0), float(33.0), float(1.0)]

cluster_dict["14.6140003204_arrows"] += cgo_arrow([-11.0,11.0,33.0], [-11.513,9.439,35.357], color="blue red", name="Arrows_14.6140003204_3")

cluster_dict["14.6140003204"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-8.5), float(2.0), float(19.5), float(1.0)]

cluster_dict["14.6140003204_arrows"] += cgo_arrow([-8.5,2.0,19.5], [-10.06,1.446,18.913], color="blue red", name="Arrows_14.6140003204_4")

cluster_dict["14.6140003204"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-14.9787234043), float(-1.42553191489), float(22.5106382979), float(1.0)]


cluster_dict["14.6140003204"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-13.3140770583), float(11.6276170016), float(22.386631294), float(1.0)]


cluster_dict["14.6140003204"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-8.44971239398), float(6.87123375698), float(29.2691769446), float(1.0)]


cluster_dict["14.6140003204"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-10.5), float(8.0), float(32.5), float(1.0)]


cluster_dict["14.6140003204"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-7.539630314), float(2.9654565324), float(18.3783623799), float(1.0)]


cluster_dict["14.6140003204"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-9.8), float(12.4), float(37.5), float(1.0)]


cluster_dict["14.6140003204"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-8.53571428571), float(12.0), float(35.9285714286), float(1.0)]


cluster_dict["14.6140003204"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-6.91531311447), float(18.3794851021), float(32.3316862645), float(1.0)]


cluster_dict["14.6140003204"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-16.5), float(8.5), float(23.0), float(1.0)]

cluster_dict["14.6140003204_arrows"] += cgo_arrow([-16.5,8.5,23.0], [-13.826,4.855,21.317], color="red blue", name="Arrows_14.6140003204_5")

cluster_dict["14.6140003204"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-10.5), float(-1.5), float(27.5), float(1.0)]

cluster_dict["14.6140003204_arrows"] += cgo_arrow([-10.5,-1.5,27.5], [-12.063,-0.184,28.703], color="red blue", name="Arrows_14.6140003204_6")

cluster_dict["14.6140003204"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-9.5), float(1.5), float(30.0), float(1.0)]

cluster_dict["14.6140003204_arrows"] += cgo_arrow([-9.5,1.5,30.0], [-12.063,-0.184,28.703], color="red blue", name="Arrows_14.6140003204_7")

cluster_dict["14.6140003204"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-4.5), float(7.5), float(29.5), float(1.0)]

cluster_dict["14.6140003204_arrows"] += cgo_arrow([-4.5,7.5,29.5], [-6.661,9.212,32.894], color="red blue", name="Arrows_14.6140003204_8")

cmd.load_cgo(cluster_dict["14.6140003204"], "Features_14.6140003204", 1)
cmd.load_cgo(cluster_dict["14.6140003204_arrows"], "Arrows_14.6140003204")
cmd.set("transparency", 0.2,"Features_14.6140003204")
cmd.group("Pharmacophore_14.6140003204", members="Features_14.6140003204")
cmd.group("Pharmacophore_14.6140003204", members="Arrows_14.6140003204")

if dirpath:
    f = join(dirpath, "70/label_threshold_14.6140003204.mol2")
else:
    f = "70/label_threshold_14.6140003204.mol2"

cmd.load(f, 'label_threshold_14.6140003204')
cmd.hide('everything', 'label_threshold_14.6140003204')
cmd.label("label_threshold_14.6140003204", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.6140003204', members= 'label_threshold_14.6140003204')


if dirpath:
    f = join(dirpath, '70/mesh.grd')
else:
    f = '70/mesh.grd'
cmd.load(f, 'mesh_70')
cmd.isomesh("isomesh_70", "mesh_70", 0.9)
cmd.color("grey80", "isomesh_70")
cmd.set('transparency', 0.4, "isomesh_70")

cmd.group('hotspot_70', "isomesh_70")
cmd.group('hotspot_70', "mesh_70")

if dirpath:
    f = join(dirpath, "71/label_threshold_17.6.mol2")
else:
    f = "71/label_threshold_17.6.mol2"

cmd.load(f, 'label_threshold_17.6')
cmd.hide('everything', 'label_threshold_17.6')
cmd.label("label_threshold_17.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [17.6]
gfiles = ['71/donor.grd', '71/apolar.grd', '71/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 71
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


cluster_dict = {"14.545999527":[], "14.545999527_arrows":[]}

cluster_dict["14.545999527"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(16.6337497023), float(-21.3062797053), float(41.7640959164), float(1.0)]


cluster_dict["14.545999527"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(10.5), float(-19.0), float(41.5), float(1.0)]

cluster_dict["14.545999527_arrows"] += cgo_arrow([10.5,-19.0,41.5], [8.555,-21.148,42.25], color="red blue", name="Arrows_14.545999527_1")

cmd.load_cgo(cluster_dict["14.545999527"], "Features_14.545999527", 1)
cmd.load_cgo(cluster_dict["14.545999527_arrows"], "Arrows_14.545999527")
cmd.set("transparency", 0.2,"Features_14.545999527")
cmd.group("Pharmacophore_14.545999527", members="Features_14.545999527")
cmd.group("Pharmacophore_14.545999527", members="Arrows_14.545999527")

if dirpath:
    f = join(dirpath, "71/label_threshold_14.545999527.mol2")
else:
    f = "71/label_threshold_14.545999527.mol2"

cmd.load(f, 'label_threshold_14.545999527')
cmd.hide('everything', 'label_threshold_14.545999527')
cmd.label("label_threshold_14.545999527", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.545999527', members= 'label_threshold_14.545999527')


if dirpath:
    f = join(dirpath, '71/mesh.grd')
else:
    f = '71/mesh.grd'
cmd.load(f, 'mesh_71')
cmd.isomesh("isomesh_71", "mesh_71", 0.9)
cmd.color("grey80", "isomesh_71")
cmd.set('transparency', 0.4, "isomesh_71")

cmd.group('hotspot_71', "isomesh_71")
cmd.group('hotspot_71', "mesh_71")

if dirpath:
    f = join(dirpath, "72/label_threshold_13.3.mol2")
else:
    f = "72/label_threshold_13.3.mol2"

cmd.load(f, 'label_threshold_13.3')
cmd.hide('everything', 'label_threshold_13.3')
cmd.label("label_threshold_13.3", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.3]
gfiles = ['72/donor.grd', '72/apolar.grd', '72/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 72
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


cluster_dict = {"14.5109996796":[], "14.5109996796_arrows":[]}

cluster_dict["14.5109996796"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-13.5), float(13.0), float(24.5), float(1.0)]

cluster_dict["14.5109996796_arrows"] += cgo_arrow([-13.5,13.0,24.5], [-16.326,13.302,26.053], color="blue red", name="Arrows_14.5109996796_1")

cluster_dict["14.5109996796"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-11.0), float(8.0), float(30.5), float(1.0)]

cluster_dict["14.5109996796_arrows"] += cgo_arrow([-11.0,8.0,30.5], [-12.682,10.491,30.534], color="blue red", name="Arrows_14.5109996796_2")

cluster_dict["14.5109996796"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-8.5), float(2.0), float(19.5), float(1.0)]

cluster_dict["14.5109996796_arrows"] += cgo_arrow([-8.5,2.0,19.5], [-10.06,1.446,18.913], color="blue red", name="Arrows_14.5109996796_3")

cluster_dict["14.5109996796"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-13.1131504989), float(10.9899565554), float(23.4589195903), float(1.0)]


cluster_dict["14.5109996796"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-14.6231884058), float(-1.94927536232), float(22.731884058), float(1.0)]


cluster_dict["14.5109996796"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-12.6329791075), float(-2.38977050276), float(21.7899892836), float(1.0)]


cluster_dict["14.5109996796"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-8.47693727892), float(6.94507360384), float(29.2159172503), float(1.0)]


cluster_dict["14.5109996796"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-10.5), float(8.0), float(32.5), float(1.0)]


cluster_dict["14.5109996796"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-8.45614035088), float(11.8274853801), float(37.9707602339), float(1.0)]


cluster_dict["14.5109996796"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-6.66426459502), float(1.54031732798), float(19.6215703431), float(1.0)]


cluster_dict["14.5109996796"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-6.72051632938), float(17.4152363341), float(34.4432035902), float(1.0)]


cluster_dict["14.5109996796"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-5.56025389861), float(-2.34575101793), float(28.3225062061), float(1.0)]


cluster_dict["14.5109996796"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-16.5), float(8.5), float(23.0), float(1.0)]

cluster_dict["14.5109996796_arrows"] += cgo_arrow([-16.5,8.5,23.0], [-13.826,4.855,21.317], color="red blue", name="Arrows_14.5109996796_4")

cluster_dict["14.5109996796"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-15.5), float(-4.5), float(35.0), float(1.0)]

cluster_dict["14.5109996796_arrows"] += cgo_arrow([-15.5,-4.5,35.0], [-14.978,-7.526,35.757], color="red blue", name="Arrows_14.5109996796_5")

cluster_dict["14.5109996796"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-10.5), float(-2.0), float(27.5), float(1.0)]

cluster_dict["14.5109996796_arrows"] += cgo_arrow([-10.5,-2.0,27.5], [-12.063,-0.184,28.703], color="red blue", name="Arrows_14.5109996796_6")

cluster_dict["14.5109996796"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-9.5), float(1.5), float(30.0), float(1.0)]

cluster_dict["14.5109996796_arrows"] += cgo_arrow([-9.5,1.5,30.0], [-12.063,-0.184,28.703], color="red blue", name="Arrows_14.5109996796_7")

cluster_dict["14.5109996796"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-8.0), float(-2.5), float(28.0), float(1.0)]

cluster_dict["14.5109996796_arrows"] += cgo_arrow([-8.0,-2.5,28.0], [-9.153,-5.211,25.883], color="red blue", name="Arrows_14.5109996796_8")

cluster_dict["14.5109996796"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-6.5), float(18.0), float(32.5), float(1.0)]


cmd.load_cgo(cluster_dict["14.5109996796"], "Features_14.5109996796", 1)
cmd.load_cgo(cluster_dict["14.5109996796_arrows"], "Arrows_14.5109996796")
cmd.set("transparency", 0.2,"Features_14.5109996796")
cmd.group("Pharmacophore_14.5109996796", members="Features_14.5109996796")
cmd.group("Pharmacophore_14.5109996796", members="Arrows_14.5109996796")

if dirpath:
    f = join(dirpath, "72/label_threshold_14.5109996796.mol2")
else:
    f = "72/label_threshold_14.5109996796.mol2"

cmd.load(f, 'label_threshold_14.5109996796')
cmd.hide('everything', 'label_threshold_14.5109996796')
cmd.label("label_threshold_14.5109996796", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.5109996796', members= 'label_threshold_14.5109996796')


if dirpath:
    f = join(dirpath, '72/mesh.grd')
else:
    f = '72/mesh.grd'
cmd.load(f, 'mesh_72')
cmd.isomesh("isomesh_72", "mesh_72", 0.9)
cmd.color("grey80", "isomesh_72")
cmd.set('transparency', 0.4, "isomesh_72")

cmd.group('hotspot_72', "isomesh_72")
cmd.group('hotspot_72', "mesh_72")

if dirpath:
    f = join(dirpath, "73/label_threshold_12.5.mol2")
else:
    f = "73/label_threshold_12.5.mol2"

cmd.load(f, 'label_threshold_12.5')
cmd.hide('everything', 'label_threshold_12.5')
cmd.label("label_threshold_12.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.5]
gfiles = ['73/donor.grd', '73/apolar.grd', '73/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 73
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


cluster_dict = {"14.4350004196":[], "14.4350004196_arrows":[]}

cluster_dict["14.4350004196"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-13.5), float(13.0), float(24.5), float(1.0)]

cluster_dict["14.4350004196_arrows"] += cgo_arrow([-13.5,13.0,24.5], [-16.326,13.302,26.053], color="blue red", name="Arrows_14.4350004196_1")

cluster_dict["14.4350004196"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-11.0), float(10.5), float(33.0), float(1.0)]

cluster_dict["14.4350004196_arrows"] += cgo_arrow([-11.0,10.5,33.0], [-11.513,9.439,35.357], color="blue red", name="Arrows_14.4350004196_2")

cluster_dict["14.4350004196"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-11.0), float(10.5), float(33.0), float(1.0)]

cluster_dict["14.4350004196_arrows"] += cgo_arrow([-11.0,10.5,33.0], [-11.513,9.439,35.357], color="blue red", name="Arrows_14.4350004196_3")

cluster_dict["14.4350004196"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-12.4015534206), float(11.9905860998), float(25.8257580363), float(1.0)]


cluster_dict["14.4350004196"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-13.0), float(7.0), float(31.5), float(1.0)]


cluster_dict["14.4350004196"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-8.67473428879), float(12.1814458432), float(38.4033376962), float(1.0)]


cluster_dict["14.4350004196"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-5.95211700933), float(19.9774944511), float(35.1770540576), float(1.0)]


cluster_dict["14.4350004196"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-6.5), float(18.0), float(32.5), float(1.0)]


cmd.load_cgo(cluster_dict["14.4350004196"], "Features_14.4350004196", 1)
cmd.load_cgo(cluster_dict["14.4350004196_arrows"], "Arrows_14.4350004196")
cmd.set("transparency", 0.2,"Features_14.4350004196")
cmd.group("Pharmacophore_14.4350004196", members="Features_14.4350004196")
cmd.group("Pharmacophore_14.4350004196", members="Arrows_14.4350004196")

if dirpath:
    f = join(dirpath, "73/label_threshold_14.4350004196.mol2")
else:
    f = "73/label_threshold_14.4350004196.mol2"

cmd.load(f, 'label_threshold_14.4350004196')
cmd.hide('everything', 'label_threshold_14.4350004196')
cmd.label("label_threshold_14.4350004196", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.4350004196', members= 'label_threshold_14.4350004196')


if dirpath:
    f = join(dirpath, '73/mesh.grd')
else:
    f = '73/mesh.grd'
cmd.load(f, 'mesh_73')
cmd.isomesh("isomesh_73", "mesh_73", 0.9)
cmd.color("grey80", "isomesh_73")
cmd.set('transparency', 0.4, "isomesh_73")

cmd.group('hotspot_73', "isomesh_73")
cmd.group('hotspot_73', "mesh_73")

if dirpath:
    f = join(dirpath, "74/label_threshold_13.0.mol2")
else:
    f = "74/label_threshold_13.0.mol2"

cmd.load(f, 'label_threshold_13.0')
cmd.hide('everything', 'label_threshold_13.0')
cmd.label("label_threshold_13.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.0]
gfiles = ['74/donor.grd', '74/apolar.grd', '74/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 74
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


cluster_dict = {"14.3979997635":[], "14.3979997635_arrows":[]}

cluster_dict["14.3979997635"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-16.0), float(-4.5), float(48.5), float(1.0)]

cluster_dict["14.3979997635_arrows"] += cgo_arrow([-16.0,-4.5,48.5], [-18.731,-2.721,49.33], color="blue red", name="Arrows_14.3979997635_1")

cluster_dict["14.3979997635"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-14.0), float(-5.5), float(43.0), float(1.0)]

cluster_dict["14.3979997635_arrows"] += cgo_arrow([-14.0,-5.5,43.0], [-13.446,-7.97,41.855], color="blue red", name="Arrows_14.3979997635_2")

cluster_dict["14.3979997635"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-14.0), float(-4.0), float(46.5), float(1.0)]

cluster_dict["14.3979997635_arrows"] += cgo_arrow([-14.0,-4.0,46.5], [-13.027,-1.485,45.573], color="blue red", name="Arrows_14.3979997635_3")

cluster_dict["14.3979997635"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-12.0), float(-13.0), float(48.0), float(1.0)]

cluster_dict["14.3979997635_arrows"] += cgo_arrow([-12.0,-13.0,48.0], [-9.173,-12.246,47.746], color="blue red", name="Arrows_14.3979997635_4")

cluster_dict["14.3979997635"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-20.8679002332), float(-19.6480078668), float(40.0102412462), float(1.0)]


cluster_dict["14.3979997635"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-14.1056849424), float(-13.6304995703), float(46.1431455873), float(1.0)]


cluster_dict["14.3979997635"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-18.0837023639), float(-13.8142208138), float(55.8847244256), float(1.0)]


cluster_dict["14.3979997635"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-18.085071121), float(-16.3107865497), float(46.0202898246), float(1.0)]


cluster_dict["14.3979997635"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-18.4335852691), float(-10.0338372049), float(38.3916036827), float(1.0)]


cluster_dict["14.3979997635"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-14.7558259818), float(-3.5), float(49.5827162496), float(1.0)]


cluster_dict["14.3979997635"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-14.0), float(-11.75), float(36.0), float(1.0)]


cluster_dict["14.3979997635"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-17.0), float(-9.5), float(47.5), float(1.0)]

cluster_dict["14.3979997635_arrows"] += cgo_arrow([-17.0,-9.5,47.5], [-17.067,-12.548,47.742], color="red blue", name="Arrows_14.3979997635_5")

cluster_dict["14.3979997635"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-14.5), float(-14.0), float(43.5), float(1.0)]

cluster_dict["14.3979997635_arrows"] += cgo_arrow([-14.5,-14.0,43.5], [-15.04,-16.828,41.999], color="red blue", name="Arrows_14.3979997635_6")

cluster_dict["14.3979997635"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-12.5), float(-13.0), float(42.5), float(1.0)]

cluster_dict["14.3979997635_arrows"] += cgo_arrow([-12.5,-13.0,42.5], [-13.853,-9.961,40.813], color="red blue", name="Arrows_14.3979997635_7")

cluster_dict["14.3979997635"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-11.0), float(-16.5), float(44.5), float(1.0)]

cluster_dict["14.3979997635_arrows"] += cgo_arrow([-11.0,-16.5,44.5], [-15.04,-16.828,41.999], color="red blue", name="Arrows_14.3979997635_8")

cmd.load_cgo(cluster_dict["14.3979997635"], "Features_14.3979997635", 1)
cmd.load_cgo(cluster_dict["14.3979997635_arrows"], "Arrows_14.3979997635")
cmd.set("transparency", 0.2,"Features_14.3979997635")
cmd.group("Pharmacophore_14.3979997635", members="Features_14.3979997635")
cmd.group("Pharmacophore_14.3979997635", members="Arrows_14.3979997635")

if dirpath:
    f = join(dirpath, "74/label_threshold_14.3979997635.mol2")
else:
    f = "74/label_threshold_14.3979997635.mol2"

cmd.load(f, 'label_threshold_14.3979997635')
cmd.hide('everything', 'label_threshold_14.3979997635')
cmd.label("label_threshold_14.3979997635", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.3979997635', members= 'label_threshold_14.3979997635')


if dirpath:
    f = join(dirpath, '74/mesh.grd')
else:
    f = '74/mesh.grd'
cmd.load(f, 'mesh_74')
cmd.isomesh("isomesh_74", "mesh_74", 0.9)
cmd.color("grey80", "isomesh_74")
cmd.set('transparency', 0.4, "isomesh_74")

cmd.group('hotspot_74', "isomesh_74")
cmd.group('hotspot_74', "mesh_74")

if dirpath:
    f = join(dirpath, "75/label_threshold_13.1.mol2")
else:
    f = "75/label_threshold_13.1.mol2"

cmd.load(f, 'label_threshold_13.1')
cmd.hide('everything', 'label_threshold_13.1')
cmd.label("label_threshold_13.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.1]
gfiles = ['75/donor.grd', '75/apolar.grd', '75/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 75
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


cluster_dict = {"14.3870000839":[], "14.3870000839_arrows":[]}

cluster_dict["14.3870000839"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-46.0), float(9.5), float(45.0), float(1.0)]

cluster_dict["14.3870000839_arrows"] += cgo_arrow([-46.0,9.5,45.0], [-47.285,10.878,46.617], color="blue red", name="Arrows_14.3870000839_1")

cluster_dict["14.3870000839"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-45.0), float(13.0), float(45.0), float(1.0)]

cluster_dict["14.3870000839_arrows"] += cgo_arrow([-45.0,13.0,45.0], [-47.048,13.086,46.553], color="blue red", name="Arrows_14.3870000839_2")

cluster_dict["14.3870000839"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-41.0), float(13.5), float(43.0), float(1.0)]

cluster_dict["14.3870000839_arrows"] += cgo_arrow([-41.0,13.5,43.0], [-39.32,11.233,41.968], color="blue red", name="Arrows_14.3870000839_3")

cluster_dict["14.3870000839"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-38.0), float(20.0), float(35.5), float(1.0)]

cluster_dict["14.3870000839_arrows"] += cgo_arrow([-38.0,20.0,35.5], [-38.645,17.942,33.279], color="blue red", name="Arrows_14.3870000839_4")

cluster_dict["14.3870000839"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-36.5), float(15.5), float(44.0), float(1.0)]

cluster_dict["14.3870000839_arrows"] += cgo_arrow([-36.5,15.5,44.0], [-36.562,17.765,45.674], color="blue red", name="Arrows_14.3870000839_5")

cluster_dict["14.3870000839"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-48.8734280503), float(20.5159214739), float(34.3380811109), float(1.0)]


cluster_dict["14.3870000839"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-45.7238769663), float(8.60103568042), float(33.6354710348), float(1.0)]


cluster_dict["14.3870000839"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-43.7543323322), float(12.8467936817), float(45.970359973), float(1.0)]


cluster_dict["14.3870000839"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-39.6161817866), float(18.7727111802), float(37.5148123677), float(1.0)]


cluster_dict["14.3870000839"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-32.9891746288), float(17.5294014846), float(47.0941986814), float(1.0)]


cluster_dict["14.3870000839"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-42.0), float(15.5), float(37.5), float(1.0)]

cluster_dict["14.3870000839_arrows"] += cgo_arrow([-42.0,15.5,37.5], [-43.081,12.971,38.706], color="red blue", name="Arrows_14.3870000839_6")

cluster_dict["14.3870000839"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-40.0), float(22.0), float(35.0), float(1.0)]

cluster_dict["14.3870000839_arrows"] += cgo_arrow([-40.0,22.0,35.0], [-42.566,21.223,34.185], color="red blue", name="Arrows_14.3870000839_7")

cluster_dict["14.3870000839"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-37.5), float(17.0), float(42.0), float(1.0)]

cluster_dict["14.3870000839_arrows"] += cgo_arrow([-37.5,17.0,42.0], [-38.333,18.672,44.584], color="red blue", name="Arrows_14.3870000839_8")

cluster_dict["14.3870000839"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-34.0), float(15.0), float(39.0), float(1.0)]

cluster_dict["14.3870000839_arrows"] += cgo_arrow([-34.0,15.0,39.0], [-32.359,12.909,40.197], color="red blue", name="Arrows_14.3870000839_9")

cmd.load_cgo(cluster_dict["14.3870000839"], "Features_14.3870000839", 1)
cmd.load_cgo(cluster_dict["14.3870000839_arrows"], "Arrows_14.3870000839")
cmd.set("transparency", 0.2,"Features_14.3870000839")
cmd.group("Pharmacophore_14.3870000839", members="Features_14.3870000839")
cmd.group("Pharmacophore_14.3870000839", members="Arrows_14.3870000839")

if dirpath:
    f = join(dirpath, "75/label_threshold_14.3870000839.mol2")
else:
    f = "75/label_threshold_14.3870000839.mol2"

cmd.load(f, 'label_threshold_14.3870000839')
cmd.hide('everything', 'label_threshold_14.3870000839')
cmd.label("label_threshold_14.3870000839", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.3870000839', members= 'label_threshold_14.3870000839')


if dirpath:
    f = join(dirpath, '75/mesh.grd')
else:
    f = '75/mesh.grd'
cmd.load(f, 'mesh_75')
cmd.isomesh("isomesh_75", "mesh_75", 0.9)
cmd.color("grey80", "isomesh_75")
cmd.set('transparency', 0.4, "isomesh_75")

cmd.group('hotspot_75', "isomesh_75")
cmd.group('hotspot_75', "mesh_75")

if dirpath:
    f = join(dirpath, "76/label_threshold_16.3.mol2")
else:
    f = "76/label_threshold_16.3.mol2"

cmd.load(f, 'label_threshold_16.3')
cmd.hide('everything', 'label_threshold_16.3')
cmd.label("label_threshold_16.3", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [16.3]
gfiles = ['76/donor.grd', '76/apolar.grd', '76/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 76
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


cluster_dict = {"14.3559999466":[], "14.3559999466_arrows":[]}

cluster_dict["14.3559999466"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(17.3721111742), float(-21.7180996811), float(41.9880169312), float(1.0)]


cluster_dict["14.3559999466"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(20.5569643483), float(-12.7021606166), float(38.2249543884), float(1.0)]


cluster_dict["14.3559999466"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.5437621188), float(-14.634295755), float(45.8543634652), float(1.0)]


cluster_dict["14.3559999466"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(10.5), float(-19.0), float(41.5), float(1.0)]

cluster_dict["14.3559999466_arrows"] += cgo_arrow([10.5,-19.0,41.5], [8.555,-21.148,42.25], color="red blue", name="Arrows_14.3559999466_1")

cmd.load_cgo(cluster_dict["14.3559999466"], "Features_14.3559999466", 1)
cmd.load_cgo(cluster_dict["14.3559999466_arrows"], "Arrows_14.3559999466")
cmd.set("transparency", 0.2,"Features_14.3559999466")
cmd.group("Pharmacophore_14.3559999466", members="Features_14.3559999466")
cmd.group("Pharmacophore_14.3559999466", members="Arrows_14.3559999466")

if dirpath:
    f = join(dirpath, "76/label_threshold_14.3559999466.mol2")
else:
    f = "76/label_threshold_14.3559999466.mol2"

cmd.load(f, 'label_threshold_14.3559999466')
cmd.hide('everything', 'label_threshold_14.3559999466')
cmd.label("label_threshold_14.3559999466", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.3559999466', members= 'label_threshold_14.3559999466')


if dirpath:
    f = join(dirpath, '76/mesh.grd')
else:
    f = '76/mesh.grd'
cmd.load(f, 'mesh_76')
cmd.isomesh("isomesh_76", "mesh_76", 0.9)
cmd.color("grey80", "isomesh_76")
cmd.set('transparency', 0.4, "isomesh_76")

cmd.group('hotspot_76', "isomesh_76")
cmd.group('hotspot_76', "mesh_76")

if dirpath:
    f = join(dirpath, "77/label_threshold_11.1.mol2")
else:
    f = "77/label_threshold_11.1.mol2"

cmd.load(f, 'label_threshold_11.1')
cmd.hide('everything', 'label_threshold_11.1')
cmd.label("label_threshold_11.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [11.1]
gfiles = ['77/donor.grd', '77/apolar.grd', '77/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 77
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


cluster_dict = {"14.2330002785":[], "14.2330002785_arrows":[]}

cluster_dict["14.2330002785"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-49.0), float(-24.5), float(9.5), float(1.0)]

cluster_dict["14.2330002785_arrows"] += cgo_arrow([-49.0,-24.5,9.5], [-47.904,-21.861,10.188], color="blue red", name="Arrows_14.2330002785_1")

cluster_dict["14.2330002785"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-46.8987996308), float(-26.2550980033), float(14.909405836), float(1.0)]


cmd.load_cgo(cluster_dict["14.2330002785"], "Features_14.2330002785", 1)
cmd.load_cgo(cluster_dict["14.2330002785_arrows"], "Arrows_14.2330002785")
cmd.set("transparency", 0.2,"Features_14.2330002785")
cmd.group("Pharmacophore_14.2330002785", members="Features_14.2330002785")
cmd.group("Pharmacophore_14.2330002785", members="Arrows_14.2330002785")

if dirpath:
    f = join(dirpath, "77/label_threshold_14.2330002785.mol2")
else:
    f = "77/label_threshold_14.2330002785.mol2"

cmd.load(f, 'label_threshold_14.2330002785')
cmd.hide('everything', 'label_threshold_14.2330002785')
cmd.label("label_threshold_14.2330002785", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.2330002785', members= 'label_threshold_14.2330002785')


if dirpath:
    f = join(dirpath, '77/mesh.grd')
else:
    f = '77/mesh.grd'
cmd.load(f, 'mesh_77')
cmd.isomesh("isomesh_77", "mesh_77", 0.9)
cmd.color("grey80", "isomesh_77")
cmd.set('transparency', 0.4, "isomesh_77")

cmd.group('hotspot_77', "isomesh_77")
cmd.group('hotspot_77', "mesh_77")

if dirpath:
    f = join(dirpath, "78/label_threshold_10.9.mol2")
else:
    f = "78/label_threshold_10.9.mol2"

cmd.load(f, 'label_threshold_10.9')
cmd.hide('everything', 'label_threshold_10.9')
cmd.label("label_threshold_10.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.9]
gfiles = ['78/donor.grd', '78/apolar.grd', '78/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 78
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


cluster_dict = {"14.2189998627":[], "14.2189998627_arrows":[]}

cluster_dict["14.2189998627"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-49.0), float(-24.5), float(9.5), float(1.0)]

cluster_dict["14.2189998627_arrows"] += cgo_arrow([-49.0,-24.5,9.5], [-47.904,-21.861,10.188], color="blue red", name="Arrows_14.2189998627_1")

cluster_dict["14.2189998627"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-46.9469413476), float(-26.3513894029), float(15.0132774893), float(1.0)]


cmd.load_cgo(cluster_dict["14.2189998627"], "Features_14.2189998627", 1)
cmd.load_cgo(cluster_dict["14.2189998627_arrows"], "Arrows_14.2189998627")
cmd.set("transparency", 0.2,"Features_14.2189998627")
cmd.group("Pharmacophore_14.2189998627", members="Features_14.2189998627")
cmd.group("Pharmacophore_14.2189998627", members="Arrows_14.2189998627")

if dirpath:
    f = join(dirpath, "78/label_threshold_14.2189998627.mol2")
else:
    f = "78/label_threshold_14.2189998627.mol2"

cmd.load(f, 'label_threshold_14.2189998627')
cmd.hide('everything', 'label_threshold_14.2189998627')
cmd.label("label_threshold_14.2189998627", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.2189998627', members= 'label_threshold_14.2189998627')


if dirpath:
    f = join(dirpath, '78/mesh.grd')
else:
    f = '78/mesh.grd'
cmd.load(f, 'mesh_78')
cmd.isomesh("isomesh_78", "mesh_78", 0.9)
cmd.color("grey80", "isomesh_78")
cmd.set('transparency', 0.4, "isomesh_78")

cmd.group('hotspot_78', "isomesh_78")
cmd.group('hotspot_78', "mesh_78")

if dirpath:
    f = join(dirpath, "79/label_threshold_10.5.mol2")
else:
    f = "79/label_threshold_10.5.mol2"

cmd.load(f, 'label_threshold_10.5')
cmd.hide('everything', 'label_threshold_10.5')
cmd.label("label_threshold_10.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.5]
gfiles = ['79/donor.grd', '79/apolar.grd', '79/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 79
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


cluster_dict = {"14.2080001831":[], "14.2080001831_arrows":[]}

cluster_dict["14.2080001831"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-49.0), float(-24.5), float(9.5), float(1.0)]

cluster_dict["14.2080001831_arrows"] += cgo_arrow([-49.0,-24.5,9.5], [-47.904,-21.861,10.188], color="blue red", name="Arrows_14.2080001831_1")

cluster_dict["14.2080001831"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-46.8692818971), float(-26.1992167387), float(14.7767451744), float(1.0)]


cmd.load_cgo(cluster_dict["14.2080001831"], "Features_14.2080001831", 1)
cmd.load_cgo(cluster_dict["14.2080001831_arrows"], "Arrows_14.2080001831")
cmd.set("transparency", 0.2,"Features_14.2080001831")
cmd.group("Pharmacophore_14.2080001831", members="Features_14.2080001831")
cmd.group("Pharmacophore_14.2080001831", members="Arrows_14.2080001831")

if dirpath:
    f = join(dirpath, "79/label_threshold_14.2080001831.mol2")
else:
    f = "79/label_threshold_14.2080001831.mol2"

cmd.load(f, 'label_threshold_14.2080001831')
cmd.hide('everything', 'label_threshold_14.2080001831')
cmd.label("label_threshold_14.2080001831", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.2080001831', members= 'label_threshold_14.2080001831')


if dirpath:
    f = join(dirpath, '79/mesh.grd')
else:
    f = '79/mesh.grd'
cmd.load(f, 'mesh_79')
cmd.isomesh("isomesh_79", "mesh_79", 0.9)
cmd.color("grey80", "isomesh_79")
cmd.set('transparency', 0.4, "isomesh_79")

cmd.group('hotspot_79', "isomesh_79")
cmd.group('hotspot_79', "mesh_79")

if dirpath:
    f = join(dirpath, "80/label_threshold_13.0.mol2")
else:
    f = "80/label_threshold_13.0.mol2"

cmd.load(f, 'label_threshold_13.0')
cmd.hide('everything', 'label_threshold_13.0')
cmd.label("label_threshold_13.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [13.0]
gfiles = ['80/donor.grd', '80/apolar.grd', '80/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 80
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


cluster_dict = {"14.1920003891":[], "14.1920003891_arrows":[]}

cluster_dict["14.1920003891"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-13.0), float(13.5), float(26.0), float(1.0)]

cluster_dict["14.1920003891_arrows"] += cgo_arrow([-13.0,13.5,26.0], [-16.326,13.302,26.053], color="blue red", name="Arrows_14.1920003891_1")

cluster_dict["14.1920003891"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-11.0), float(8.0), float(30.5), float(1.0)]

cluster_dict["14.1920003891_arrows"] += cgo_arrow([-11.0,8.0,30.5], [-12.682,10.491,30.534], color="blue red", name="Arrows_14.1920003891_2")

cluster_dict["14.1920003891"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-11.0), float(11.0), float(33.0), float(1.0)]

cluster_dict["14.1920003891_arrows"] += cgo_arrow([-11.0,11.0,33.0], [-11.513,9.439,35.357], color="blue red", name="Arrows_14.1920003891_3")

cluster_dict["14.1920003891"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(0.0), float(3.0), float(39.5), float(1.0)]

cluster_dict["14.1920003891_arrows"] += cgo_arrow([0.0,3.0,39.5], [0.089,2.153,36.754], color="blue red", name="Arrows_14.1920003891_4")

cluster_dict["14.1920003891"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(2.0), float(4.5), float(38.0), float(1.0)]

cluster_dict["14.1920003891_arrows"] += cgo_arrow([2.0,4.5,38.0], [1.181,3.849,35.951], color="blue red", name="Arrows_14.1920003891_5")

cluster_dict["14.1920003891"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-11.3995771982), float(11.7982399038), float(27.1448450638), float(1.0)]


cluster_dict["14.1920003891"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-13.0), float(7.0), float(31.5), float(1.0)]


cluster_dict["14.1920003891"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-8.44192529802), float(7.0648532135), float(29.2370534073), float(1.0)]


cluster_dict["14.1920003891"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-10.5), float(8.0), float(32.5), float(1.0)]


cluster_dict["14.1920003891"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-8.36712058005), float(11.903422477), float(38.6293595598), float(1.0)]


cluster_dict["14.1920003891"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-5.909877518), float(20.0318244416), float(35.2771432446), float(1.0)]


cluster_dict["14.1920003891"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(0.571253071286), float(11.2453121261), float(44.9148796318), float(1.0)]


cluster_dict["14.1920003891"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(1.96690911927), float(4.49531258052), float(39.6723968907), float(1.0)]


cluster_dict["14.1920003891"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-9.5), float(1.5), float(30.0), float(1.0)]

cluster_dict["14.1920003891_arrows"] += cgo_arrow([-9.5,1.5,30.0], [-12.063,-0.184,28.703], color="red blue", name="Arrows_14.1920003891_6")

cluster_dict["14.1920003891"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-6.5), float(18.0), float(32.5), float(1.0)]


cluster_dict["14.1920003891"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-4.5), float(7.5), float(29.5), float(1.0)]

cluster_dict["14.1920003891_arrows"] += cgo_arrow([-4.5,7.5,29.5], [-6.661,9.212,32.894], color="red blue", name="Arrows_14.1920003891_7")

cluster_dict["14.1920003891"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(3.5), float(6.5), float(38.0), float(1.0)]

cluster_dict["14.1920003891_arrows"] += cgo_arrow([3.5,6.5,38.0], [3.554,7.603,40.919], color="red blue", name="Arrows_14.1920003891_8")

cmd.load_cgo(cluster_dict["14.1920003891"], "Features_14.1920003891", 1)
cmd.load_cgo(cluster_dict["14.1920003891_arrows"], "Arrows_14.1920003891")
cmd.set("transparency", 0.2,"Features_14.1920003891")
cmd.group("Pharmacophore_14.1920003891", members="Features_14.1920003891")
cmd.group("Pharmacophore_14.1920003891", members="Arrows_14.1920003891")

if dirpath:
    f = join(dirpath, "80/label_threshold_14.1920003891.mol2")
else:
    f = "80/label_threshold_14.1920003891.mol2"

cmd.load(f, 'label_threshold_14.1920003891')
cmd.hide('everything', 'label_threshold_14.1920003891')
cmd.label("label_threshold_14.1920003891", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_14.1920003891', members= 'label_threshold_14.1920003891')


if dirpath:
    f = join(dirpath, '80/mesh.grd')
else:
    f = '80/mesh.grd'
cmd.load(f, 'mesh_80')
cmd.isomesh("isomesh_80", "mesh_80", 0.9)
cmd.color("grey80", "isomesh_80")
cmd.set('transparency', 0.4, "isomesh_80")

cmd.group('hotspot_80', "isomesh_80")
cmd.group('hotspot_80', "mesh_80")

if dirpath:
    f = join(dirpath, "81/label_threshold_12.4.mol2")
else:
    f = "81/label_threshold_12.4.mol2"

cmd.load(f, 'label_threshold_12.4')
cmd.hide('everything', 'label_threshold_12.4')
cmd.label("label_threshold_12.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.4]
gfiles = ['81/donor.grd', '81/apolar.grd', '81/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 81
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


cluster_dict = {"13.9469995499":[], "13.9469995499_arrows":[]}

cluster_dict["13.9469995499"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-16.0), float(-9.5), float(28.5), float(1.0)]

cluster_dict["13.9469995499_arrows"] += cgo_arrow([-16.0,-9.5,28.5], [-14.66,-11.484,26.814], color="blue red", name="Arrows_13.9469995499_1")

cluster_dict["13.9469995499"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-13.5), float(-6.0), float(23.0), float(1.0)]

cluster_dict["13.9469995499_arrows"] += cgo_arrow([-13.5,-6.0,23.0], [-15.302,-8.206,23.93], color="blue red", name="Arrows_13.9469995499_2")

cluster_dict["13.9469995499"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-11.0), float(8.0), float(30.5), float(1.0)]

cluster_dict["13.9469995499_arrows"] += cgo_arrow([-11.0,8.0,30.5], [-12.682,10.491,30.534], color="blue red", name="Arrows_13.9469995499_3")

cluster_dict["13.9469995499"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-8.5), float(2.0), float(19.5), float(1.0)]

cluster_dict["13.9469995499_arrows"] += cgo_arrow([-8.5,2.0,19.5], [-10.06,1.446,18.913], color="blue red", name="Arrows_13.9469995499_4")

cluster_dict["13.9469995499"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-14.4433137159), float(-2.64585542932), float(22.4807206617), float(1.0)]


cluster_dict["13.9469995499"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-13.0), float(-5.4), float(21.751667117), float(1.0)]


cluster_dict["13.9469995499"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-13.0), float(7.0), float(31.5), float(1.0)]


cluster_dict["13.9469995499"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-8.21415240199), float(6.51388732341), float(29.0773565518), float(1.0)]


cluster_dict["13.9469995499"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-4.94423986875), float(-0.346612196086), float(19.2273797962), float(1.0)]


cluster_dict["13.9469995499"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-5.52432432048), float(-2.42051252927), float(28.3113777882), float(1.0)]


cluster_dict["13.9469995499"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-15.5), float(-4.5), float(35.0), float(1.0)]

cluster_dict["13.9469995499_arrows"] += cgo_arrow([-15.5,-4.5,35.0], [-14.978,-7.526,35.757], color="red blue", name="Arrows_13.9469995499_5")

cluster_dict["13.9469995499"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-10.5), float(-2.0), float(27.5), float(1.0)]

cluster_dict["13.9469995499_arrows"] += cgo_arrow([-10.5,-2.0,27.5], [-12.063,-0.184,28.703], color="red blue", name="Arrows_13.9469995499_6")

cluster_dict["13.9469995499"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-9.5), float(1.5), float(30.0), float(1.0)]

cluster_dict["13.9469995499_arrows"] += cgo_arrow([-9.5,1.5,30.0], [-12.063,-0.184,28.703], color="red blue", name="Arrows_13.9469995499_7")

cluster_dict["13.9469995499"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-8.0), float(-2.5), float(28.0), float(1.0)]

cluster_dict["13.9469995499_arrows"] += cgo_arrow([-8.0,-2.5,28.0], [-9.153,-5.211,25.883], color="red blue", name="Arrows_13.9469995499_8")

cluster_dict["13.9469995499"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-4.5), float(7.5), float(29.5), float(1.0)]

cluster_dict["13.9469995499_arrows"] += cgo_arrow([-4.5,7.5,29.5], [-6.661,9.212,32.894], color="red blue", name="Arrows_13.9469995499_9")

cmd.load_cgo(cluster_dict["13.9469995499"], "Features_13.9469995499", 1)
cmd.load_cgo(cluster_dict["13.9469995499_arrows"], "Arrows_13.9469995499")
cmd.set("transparency", 0.2,"Features_13.9469995499")
cmd.group("Pharmacophore_13.9469995499", members="Features_13.9469995499")
cmd.group("Pharmacophore_13.9469995499", members="Arrows_13.9469995499")

if dirpath:
    f = join(dirpath, "81/label_threshold_13.9469995499.mol2")
else:
    f = "81/label_threshold_13.9469995499.mol2"

cmd.load(f, 'label_threshold_13.9469995499')
cmd.hide('everything', 'label_threshold_13.9469995499')
cmd.label("label_threshold_13.9469995499", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.9469995499', members= 'label_threshold_13.9469995499')


if dirpath:
    f = join(dirpath, '81/mesh.grd')
else:
    f = '81/mesh.grd'
cmd.load(f, 'mesh_81')
cmd.isomesh("isomesh_81", "mesh_81", 0.9)
cmd.color("grey80", "isomesh_81")
cmd.set('transparency', 0.4, "isomesh_81")

cmd.group('hotspot_81', "isomesh_81")
cmd.group('hotspot_81', "mesh_81")

if dirpath:
    f = join(dirpath, "82/label_threshold_12.1.mol2")
else:
    f = "82/label_threshold_12.1.mol2"

cmd.load(f, 'label_threshold_12.1')
cmd.hide('everything', 'label_threshold_12.1')
cmd.label("label_threshold_12.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.1]
gfiles = ['82/donor.grd', '82/apolar.grd', '82/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 82
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


cluster_dict = {"13.9250001907":[], "13.9250001907_arrows":[]}

cluster_dict["13.9250001907"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-56.5), float(-8.5), float(50.0), float(1.0)]

cluster_dict["13.9250001907_arrows"] += cgo_arrow([-56.5,-8.5,50.0], [-58.938,-9.438,49.082], color="blue red", name="Arrows_13.9250001907_1")

cluster_dict["13.9250001907"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-57.5), float(-4.5), float(50.0), float(1.0)]

cluster_dict["13.9250001907_arrows"] += cgo_arrow([-57.5,-4.5,50.0], [-59.288,-5.348,50.967], color="blue red", name="Arrows_13.9250001907_2")

cluster_dict["13.9250001907"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-55.0), float(-3.5), float(50.0), float(1.0)]

cluster_dict["13.9250001907_arrows"] += cgo_arrow([-55.0,-3.5,50.0], [-56.059,-0.837,49.615], color="blue red", name="Arrows_13.9250001907_3")

cluster_dict["13.9250001907"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-53.0), float(-1.5), float(53.5), float(1.0)]

cluster_dict["13.9250001907_arrows"] += cgo_arrow([-53.0,-1.5,53.5], [-51.425,-1.643,56.198], color="blue red", name="Arrows_13.9250001907_4")

cluster_dict["13.9250001907"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-52.0), float(-5.5), float(65.5), float(1.0)]

cluster_dict["13.9250001907_arrows"] += cgo_arrow([-52.0,-5.5,65.5], [-49.259,-5.105,64.737], color="blue red", name="Arrows_13.9250001907_5")

cluster_dict["13.9250001907"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-57.1377796136), float(-7.66490791315), float(53.9435505648), float(1.0)]


cluster_dict["13.9250001907"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-52.944787615), float(-8.08664816046), float(65.2375923453), float(1.0)]


cluster_dict["13.9250001907"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-58.5), float(-6.5), float(49.0), float(1.0)]

cluster_dict["13.9250001907_arrows"] += cgo_arrow([-58.5,-6.5,49.0], [-59.288,-5.348,50.967], color="red blue", name="Arrows_13.9250001907_6")

cluster_dict["13.9250001907"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-56.0), float(-7.0), float(58.0), float(1.0)]

cluster_dict["13.9250001907_arrows"] += cgo_arrow([-56.0,-7.0,58.0], [-54.596,-8.072,60.438], color="red blue", name="Arrows_13.9250001907_7")

cluster_dict["13.9250001907"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-56.5), float(-4.5), float(50.5), float(1.0)]

cluster_dict["13.9250001907_arrows"] += cgo_arrow([-56.5,-4.5,50.5], [-59.288,-5.348,50.967], color="red blue", name="Arrows_13.9250001907_8")

cluster_dict["13.9250001907"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-56.0), float(-10.0), float(51.0), float(1.0)]

cluster_dict["13.9250001907_arrows"] += cgo_arrow([-56.0,-10.0,51.0], [-60.269,-11.261,49.226], color="red blue", name="Arrows_13.9250001907_9")

cmd.load_cgo(cluster_dict["13.9250001907"], "Features_13.9250001907", 1)
cmd.load_cgo(cluster_dict["13.9250001907_arrows"], "Arrows_13.9250001907")
cmd.set("transparency", 0.2,"Features_13.9250001907")
cmd.group("Pharmacophore_13.9250001907", members="Features_13.9250001907")
cmd.group("Pharmacophore_13.9250001907", members="Arrows_13.9250001907")

if dirpath:
    f = join(dirpath, "82/label_threshold_13.9250001907.mol2")
else:
    f = "82/label_threshold_13.9250001907.mol2"

cmd.load(f, 'label_threshold_13.9250001907')
cmd.hide('everything', 'label_threshold_13.9250001907')
cmd.label("label_threshold_13.9250001907", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.9250001907', members= 'label_threshold_13.9250001907')


if dirpath:
    f = join(dirpath, '82/mesh.grd')
else:
    f = '82/mesh.grd'
cmd.load(f, 'mesh_82')
cmd.isomesh("isomesh_82", "mesh_82", 0.9)
cmd.color("grey80", "isomesh_82")
cmd.set('transparency', 0.4, "isomesh_82")

cmd.group('hotspot_82', "isomesh_82")
cmd.group('hotspot_82', "mesh_82")

if dirpath:
    f = join(dirpath, "83/label_threshold_11.2.mol2")
else:
    f = "83/label_threshold_11.2.mol2"

cmd.load(f, 'label_threshold_11.2')
cmd.hide('everything', 'label_threshold_11.2')
cmd.label("label_threshold_11.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [11.2]
gfiles = ['83/donor.grd', '83/apolar.grd', '83/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 83
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


cluster_dict = {"13.8090000153":[], "13.8090000153_arrows":[]}

cluster_dict["13.8090000153"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(9.5), float(9.5), float(46.5), float(1.0)]

cluster_dict["13.8090000153_arrows"] += cgo_arrow([9.5,9.5,46.5], [8.374,12.84,46.167], color="blue red", name="Arrows_13.8090000153_1")

cluster_dict["13.8090000153"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(12.0), float(10.0), float(48.5), float(1.0)]

cluster_dict["13.8090000153_arrows"] += cgo_arrow([12.0,10.0,48.5], [9.932,10.567,50.575], color="blue red", name="Arrows_13.8090000153_2")

cluster_dict["13.8090000153"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(13.0), float(3.5), float(47.5), float(1.0)]

cluster_dict["13.8090000153_arrows"] += cgo_arrow([13.0,3.5,47.5], [14.765,4.449,46.227], color="blue red", name="Arrows_13.8090000153_3")

cluster_dict["13.8090000153"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(14.5), float(-7.5), float(52.5), float(1.0)]

cluster_dict["13.8090000153_arrows"] += cgo_arrow([14.5,-7.5,52.5], [13.149,-6.193,54.717], color="blue red", name="Arrows_13.8090000153_4")

cluster_dict["13.8090000153"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(10.2685917642), float(6.57644047103), float(49.6366144794), float(1.0)]


cluster_dict["13.8090000153"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(15.758455154), float(-5.04849426), float(51.1604033959), float(1.0)]


cluster_dict["13.8090000153"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(13.75), float(7.5), float(51.0), float(1.0)]


cluster_dict["13.8090000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(2.0), float(7.5), float(49.5), float(1.0)]

cluster_dict["13.8090000153_arrows"] += cgo_arrow([2.0,7.5,49.5], [1.122,5.617,48.391], color="red blue", name="Arrows_13.8090000153_5")

cluster_dict["13.8090000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(5.5), float(11.0), float(50.0), float(1.0)]

cluster_dict["13.8090000153_arrows"] += cgo_arrow([5.5,11.0,50.0], [7.497,8.965,50.988], color="red blue", name="Arrows_13.8090000153_6")

cluster_dict["13.8090000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(10.0), float(6.5), float(50.5), float(1.0)]

cluster_dict["13.8090000153_arrows"] += cgo_arrow([10.0,6.5,50.5], [7.899,6.991,52.135], color="red blue", name="Arrows_13.8090000153_7")

cluster_dict["13.8090000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(11.0), float(6.5), float(45.0), float(1.0)]

cluster_dict["13.8090000153_arrows"] += cgo_arrow([11.0,6.5,45.0], [11.032,3.68,45.103], color="red blue", name="Arrows_13.8090000153_8")

cluster_dict["13.8090000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(-2.5), float(48.0), float(1.0)]

cluster_dict["13.8090000153_arrows"] += cgo_arrow([15.0,-2.5,48.0], [15.22,-0.595,51.059], color="red blue", name="Arrows_13.8090000153_9")

cluster_dict["13.8090000153"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(16.5), float(-1.0), float(48.0), float(1.0)]

cluster_dict["13.8090000153_arrows"] += cgo_arrow([16.5,-1.0,48.0], [17.164,1.471,50.049], color="red blue", name="Arrows_13.8090000153_10")

cmd.load_cgo(cluster_dict["13.8090000153"], "Features_13.8090000153", 1)
cmd.load_cgo(cluster_dict["13.8090000153_arrows"], "Arrows_13.8090000153")
cmd.set("transparency", 0.2,"Features_13.8090000153")
cmd.group("Pharmacophore_13.8090000153", members="Features_13.8090000153")
cmd.group("Pharmacophore_13.8090000153", members="Arrows_13.8090000153")

if dirpath:
    f = join(dirpath, "83/label_threshold_13.8090000153.mol2")
else:
    f = "83/label_threshold_13.8090000153.mol2"

cmd.load(f, 'label_threshold_13.8090000153')
cmd.hide('everything', 'label_threshold_13.8090000153')
cmd.label("label_threshold_13.8090000153", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.8090000153', members= 'label_threshold_13.8090000153')


if dirpath:
    f = join(dirpath, '83/mesh.grd')
else:
    f = '83/mesh.grd'
cmd.load(f, 'mesh_83')
cmd.isomesh("isomesh_83", "mesh_83", 0.9)
cmd.color("grey80", "isomesh_83")
cmd.set('transparency', 0.4, "isomesh_83")

cmd.group('hotspot_83', "isomesh_83")
cmd.group('hotspot_83', "mesh_83")

if dirpath:
    f = join(dirpath, "84/label_threshold_11.1.mol2")
else:
    f = "84/label_threshold_11.1.mol2"

cmd.load(f, 'label_threshold_11.1')
cmd.hide('everything', 'label_threshold_11.1')
cmd.label("label_threshold_11.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [11.1]
gfiles = ['84/donor.grd', '84/apolar.grd', '84/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 84
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


cluster_dict = {"13.7740001678":[], "13.7740001678_arrows":[]}

cluster_dict["13.7740001678"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-17.5), float(-10.5), float(28.0), float(1.0)]

cluster_dict["13.7740001678_arrows"] += cgo_arrow([-17.5,-10.5,28.0], [-19.572,-9.949,26.024], color="blue red", name="Arrows_13.7740001678_1")

cluster_dict["13.7740001678"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-16.5), float(-11.5), float(32.0), float(1.0)]

cluster_dict["13.7740001678_arrows"] += cgo_arrow([-16.5,-11.5,32.0], [-18.259,-9.567,30.976], color="blue red", name="Arrows_13.7740001678_2")

cluster_dict["13.7740001678"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-13.5), float(-6.0), float(23.0), float(1.0)]

cluster_dict["13.7740001678_arrows"] += cgo_arrow([-13.5,-6.0,23.0], [-15.302,-8.206,23.93], color="blue red", name="Arrows_13.7740001678_3")

cluster_dict["13.7740001678"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-14.0), float(-2.0), float(25.0), float(1.0)]

cluster_dict["13.7740001678_arrows"] += cgo_arrow([-14.0,-2.0,25.0], [-14.141,-2.336,27.942], color="blue red", name="Arrows_13.7740001678_4")

cluster_dict["13.7740001678"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-8.5), float(2.0), float(19.5), float(1.0)]

cluster_dict["13.7740001678_arrows"] += cgo_arrow([-8.5,2.0,19.5], [-10.06,1.446,18.913], color="blue red", name="Arrows_13.7740001678_5")

cluster_dict["13.7740001678"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-16.6001710987), float(-12.5644473837), float(29.0537832513), float(1.0)]


cluster_dict["13.7740001678"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-14.3716340986), float(-3.36595784527), float(22.2961981316), float(1.0)]


cluster_dict["13.7740001678"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-14.876422277), float(7.93809100776), float(23.0695410219), float(1.0)]


cluster_dict["13.7740001678"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-9.26183906667), float(3.39759338228), float(29.2837605687), float(1.0)]


cluster_dict["13.7740001678"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-5.85602411237), float(0.0967710881865), float(18.878155802), float(1.0)]


cluster_dict["13.7740001678"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-5.4564138048), float(-2.39918215313), float(28.2879158581), float(1.0)]


cluster_dict["13.7740001678"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-16.5), float(-10.0), float(27.0), float(1.0)]

cluster_dict["13.7740001678_arrows"] += cgo_arrow([-16.5,-10.0,27.0], [-19.572,-9.949,26.024], color="red blue", name="Arrows_13.7740001678_6")

cluster_dict["13.7740001678"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-16.0), float(8.5), float(22.5), float(1.0)]

cluster_dict["13.7740001678_arrows"] += cgo_arrow([-16.0,8.5,22.5], [-12.952,6.981,19.841], color="red blue", name="Arrows_13.7740001678_7")

cluster_dict["13.7740001678"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-10.5), float(-2.0), float(27.5), float(1.0)]

cluster_dict["13.7740001678_arrows"] += cgo_arrow([-10.5,-2.0,27.5], [-12.063,-0.184,28.703], color="red blue", name="Arrows_13.7740001678_8")

cluster_dict["13.7740001678"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-9.5), float(1.5), float(30.0), float(1.0)]

cluster_dict["13.7740001678_arrows"] += cgo_arrow([-9.5,1.5,30.0], [-12.063,-0.184,28.703], color="red blue", name="Arrows_13.7740001678_9")

cluster_dict["13.7740001678"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-8.0), float(-2.5), float(28.0), float(1.0)]

cluster_dict["13.7740001678_arrows"] += cgo_arrow([-8.0,-2.5,28.0], [-9.153,-5.211,25.883], color="red blue", name="Arrows_13.7740001678_10")

cmd.load_cgo(cluster_dict["13.7740001678"], "Features_13.7740001678", 1)
cmd.load_cgo(cluster_dict["13.7740001678_arrows"], "Arrows_13.7740001678")
cmd.set("transparency", 0.2,"Features_13.7740001678")
cmd.group("Pharmacophore_13.7740001678", members="Features_13.7740001678")
cmd.group("Pharmacophore_13.7740001678", members="Arrows_13.7740001678")

if dirpath:
    f = join(dirpath, "84/label_threshold_13.7740001678.mol2")
else:
    f = "84/label_threshold_13.7740001678.mol2"

cmd.load(f, 'label_threshold_13.7740001678')
cmd.hide('everything', 'label_threshold_13.7740001678')
cmd.label("label_threshold_13.7740001678", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.7740001678', members= 'label_threshold_13.7740001678')


if dirpath:
    f = join(dirpath, '84/mesh.grd')
else:
    f = '84/mesh.grd'
cmd.load(f, 'mesh_84')
cmd.isomesh("isomesh_84", "mesh_84", 0.9)
cmd.color("grey80", "isomesh_84")
cmd.set('transparency', 0.4, "isomesh_84")

cmd.group('hotspot_84', "isomesh_84")
cmd.group('hotspot_84', "mesh_84")

if dirpath:
    f = join(dirpath, "85/label_threshold_11.4.mol2")
else:
    f = "85/label_threshold_11.4.mol2"

cmd.load(f, 'label_threshold_11.4')
cmd.hide('everything', 'label_threshold_11.4')
cmd.label("label_threshold_11.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [11.4]
gfiles = ['85/donor.grd', '85/apolar.grd', '85/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 85
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


cluster_dict = {"13.7399997711":[], "13.7399997711_arrows":[]}

cluster_dict["13.7399997711"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-13.5), float(13.0), float(24.5), float(1.0)]

cluster_dict["13.7399997711_arrows"] += cgo_arrow([-13.5,13.0,24.5], [-16.326,13.302,26.053], color="blue red", name="Arrows_13.7399997711_1")

cluster_dict["13.7399997711"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-11.5), float(8.5), float(29.0), float(1.0)]

cluster_dict["13.7399997711_arrows"] += cgo_arrow([-11.5,8.5,29.0], [-12.709,6.338,28.608], color="blue red", name="Arrows_13.7399997711_2")

cluster_dict["13.7399997711"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-12.5993399479), float(12.0289663623), float(23.1913352335), float(1.0)]


cluster_dict["13.7399997711"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-16.5), float(8.5), float(23.0), float(1.0)]

cluster_dict["13.7399997711_arrows"] += cgo_arrow([-16.5,8.5,23.0], [-13.826,4.855,21.317], color="red blue", name="Arrows_13.7399997711_3")

cmd.load_cgo(cluster_dict["13.7399997711"], "Features_13.7399997711", 1)
cmd.load_cgo(cluster_dict["13.7399997711_arrows"], "Arrows_13.7399997711")
cmd.set("transparency", 0.2,"Features_13.7399997711")
cmd.group("Pharmacophore_13.7399997711", members="Features_13.7399997711")
cmd.group("Pharmacophore_13.7399997711", members="Arrows_13.7399997711")

if dirpath:
    f = join(dirpath, "85/label_threshold_13.7399997711.mol2")
else:
    f = "85/label_threshold_13.7399997711.mol2"

cmd.load(f, 'label_threshold_13.7399997711')
cmd.hide('everything', 'label_threshold_13.7399997711')
cmd.label("label_threshold_13.7399997711", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.7399997711', members= 'label_threshold_13.7399997711')


if dirpath:
    f = join(dirpath, '85/mesh.grd')
else:
    f = '85/mesh.grd'
cmd.load(f, 'mesh_85')
cmd.isomesh("isomesh_85", "mesh_85", 0.9)
cmd.color("grey80", "isomesh_85")
cmd.set('transparency', 0.4, "isomesh_85")

cmd.group('hotspot_85', "isomesh_85")
cmd.group('hotspot_85', "mesh_85")

if dirpath:
    f = join(dirpath, "86/label_threshold_12.0.mol2")
else:
    f = "86/label_threshold_12.0.mol2"

cmd.load(f, 'label_threshold_12.0')
cmd.hide('everything', 'label_threshold_12.0')
cmd.label("label_threshold_12.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.0]
gfiles = ['86/donor.grd', '86/apolar.grd', '86/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 86
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


cluster_dict = {"13.7110004425":[], "13.7110004425_arrows":[]}

cluster_dict["13.7110004425"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(13.0), float(3.5), float(47.5), float(1.0)]

cluster_dict["13.7110004425_arrows"] += cgo_arrow([13.0,3.5,47.5], [14.765,4.449,46.227], color="blue red", name="Arrows_13.7110004425_1")

cluster_dict["13.7110004425"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(-7.5), float(52.5), float(1.0)]

cluster_dict["13.7110004425_arrows"] += cgo_arrow([15.0,-7.5,52.5], [13.149,-6.193,54.717], color="blue red", name="Arrows_13.7110004425_2")

cluster_dict["13.7110004425"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(14.5), float(-8.5), float(50.0), float(1.0)]

cluster_dict["13.7110004425_arrows"] += cgo_arrow([14.5,-8.5,50.0], [14.608,-9.503,48.087], color="blue red", name="Arrows_13.7110004425_3")

cluster_dict["13.7110004425"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(12.2444643293), float(3.95395410119), float(50.1112878361), float(1.0)]


cluster_dict["13.7110004425"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(16.0426836234), float(-5.48987016669), float(48.8048528235), float(1.0)]


cluster_dict["13.7110004425"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(20.7505034629), float(1.02424384858), float(42.9517870854), float(1.0)]


cluster_dict["13.7110004425"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(20.0326935465), float(-9.60817703608), float(39.2354955571), float(1.0)]


cluster_dict["13.7110004425"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.7512504406), float(-9.27946762598), float(48.221140557), float(1.0)]


cluster_dict["13.7110004425"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(12.0), float(5.0), float(51.5), float(1.0)]

cluster_dict["13.7110004425_arrows"] += cgo_arrow([12.0,5.0,51.5], [11.081,7.326,54.016], color="red blue", name="Arrows_13.7110004425_4")

cluster_dict["13.7110004425"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(13.5), float(-7.0), float(48.0), float(1.0)]

cluster_dict["13.7110004425_arrows"] += cgo_arrow([13.5,-7.0,48.0], [11.961,-9.207,47.387], color="red blue", name="Arrows_13.7110004425_5")

cluster_dict["13.7110004425"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(-2.5), float(48.0), float(1.0)]

cluster_dict["13.7110004425_arrows"] += cgo_arrow([15.0,-2.5,48.0], [15.22,-0.595,51.059], color="red blue", name="Arrows_13.7110004425_6")

cluster_dict["13.7110004425"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.0), float(-8.5), float(52.0), float(1.0)]

cluster_dict["13.7110004425_arrows"] += cgo_arrow([17.0,-8.5,52.0], [19.182,-9.441,50.638], color="red blue", name="Arrows_13.7110004425_7")

cluster_dict["13.7110004425"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(-5.5), float(45.5), float(1.0)]

cluster_dict["13.7110004425_arrows"] += cgo_arrow([18.5,-5.5,45.5], [19.124,-4.683,42.679], color="red blue", name="Arrows_13.7110004425_8")

cluster_dict["13.7110004425"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(19.0), float(-0.5), float(42.5), float(1.0)]

cluster_dict["13.7110004425_arrows"] += cgo_arrow([19.0,-0.5,42.5], [15.944,-0.291,41.945], color="red blue", name="Arrows_13.7110004425_9")

cmd.load_cgo(cluster_dict["13.7110004425"], "Features_13.7110004425", 1)
cmd.load_cgo(cluster_dict["13.7110004425_arrows"], "Arrows_13.7110004425")
cmd.set("transparency", 0.2,"Features_13.7110004425")
cmd.group("Pharmacophore_13.7110004425", members="Features_13.7110004425")
cmd.group("Pharmacophore_13.7110004425", members="Arrows_13.7110004425")

if dirpath:
    f = join(dirpath, "86/label_threshold_13.7110004425.mol2")
else:
    f = "86/label_threshold_13.7110004425.mol2"

cmd.load(f, 'label_threshold_13.7110004425')
cmd.hide('everything', 'label_threshold_13.7110004425')
cmd.label("label_threshold_13.7110004425", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.7110004425', members= 'label_threshold_13.7110004425')


if dirpath:
    f = join(dirpath, '86/mesh.grd')
else:
    f = '86/mesh.grd'
cmd.load(f, 'mesh_86')
cmd.isomesh("isomesh_86", "mesh_86", 0.9)
cmd.color("grey80", "isomesh_86")
cmd.set('transparency', 0.4, "isomesh_86")

cmd.group('hotspot_86', "isomesh_86")
cmd.group('hotspot_86', "mesh_86")

if dirpath:
    f = join(dirpath, "87/label_threshold_12.8.mol2")
else:
    f = "87/label_threshold_12.8.mol2"

cmd.load(f, 'label_threshold_12.8')
cmd.hide('everything', 'label_threshold_12.8')
cmd.label("label_threshold_12.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.8]
gfiles = ['87/donor.grd', '87/apolar.grd', '87/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 87
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


cluster_dict = {"13.6499996185":[], "13.6499996185_arrows":[]}

cluster_dict["13.6499996185"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-16.0), float(-3.0), float(54.0), float(1.0)]

cluster_dict["13.6499996185_arrows"] += cgo_arrow([-16.0,-3.0,54.0], [-17.363,-0.587,52.277], color="blue red", name="Arrows_13.6499996185_1")

cluster_dict["13.6499996185"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-16.0), float(-4.5), float(48.5), float(1.0)]

cluster_dict["13.6499996185_arrows"] += cgo_arrow([-16.0,-4.5,48.5], [-18.731,-2.721,49.33], color="blue red", name="Arrows_13.6499996185_2")

cluster_dict["13.6499996185"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-14.0), float(-4.0), float(46.5), float(1.0)]

cluster_dict["13.6499996185_arrows"] += cgo_arrow([-14.0,-4.0,46.5], [-13.027,-1.485,45.573], color="blue red", name="Arrows_13.6499996185_3")

cluster_dict["13.6499996185"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-12.5), float(-13.0), float(48.0), float(1.0)]

cluster_dict["13.6499996185_arrows"] += cgo_arrow([-12.5,-13.0,48.0], [-9.173,-12.246,47.746], color="blue red", name="Arrows_13.6499996185_4")

cluster_dict["13.6499996185"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-16.0649186988), float(-8.23737237506), float(45.202876122), float(1.0)]


cluster_dict["13.6499996185"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-22.1971454389), float(9.5), float(51.2066606425), float(1.0)]


cluster_dict["13.6499996185"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-16.177163156), float(-2.44419296213), float(56.9850837406), float(1.0)]


cluster_dict["13.6499996185"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-17.1196181405), float(-12.2370910542), float(54.9576158399), float(1.0)]


cluster_dict["13.6499996185"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-14.5931309941), float(-0.837609457242), float(49.3634838222), float(1.0)]


cluster_dict["13.6499996185"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-21.0), float(-6.0), float(45.5), float(1.0)]

cluster_dict["13.6499996185_arrows"] += cgo_arrow([-21.0,-6.0,45.5], [-22.337,-3.396,48.966], color="red blue", name="Arrows_13.6499996185_5")

cluster_dict["13.6499996185"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-19.0), float(-4.0), float(56.5), float(1.0)]

cluster_dict["13.6499996185_arrows"] += cgo_arrow([-19.0,-4.0,56.5], [-21.532,-6.13,56.958], color="red blue", name="Arrows_13.6499996185_6")

cluster_dict["13.6499996185"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-18.0), float(-2.0), float(57.5), float(1.0)]


cluster_dict["13.6499996185"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-17.0), float(-3.5), float(55.5), float(1.0)]


cluster_dict["13.6499996185"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-15.0), float(-0.5), float(46.5), float(1.0)]

cluster_dict["13.6499996185_arrows"] += cgo_arrow([-15.0,-0.5,46.5], [-17.71,0.005,45.762], color="red blue", name="Arrows_13.6499996185_7")

cmd.load_cgo(cluster_dict["13.6499996185"], "Features_13.6499996185", 1)
cmd.load_cgo(cluster_dict["13.6499996185_arrows"], "Arrows_13.6499996185")
cmd.set("transparency", 0.2,"Features_13.6499996185")
cmd.group("Pharmacophore_13.6499996185", members="Features_13.6499996185")
cmd.group("Pharmacophore_13.6499996185", members="Arrows_13.6499996185")

if dirpath:
    f = join(dirpath, "87/label_threshold_13.6499996185.mol2")
else:
    f = "87/label_threshold_13.6499996185.mol2"

cmd.load(f, 'label_threshold_13.6499996185')
cmd.hide('everything', 'label_threshold_13.6499996185')
cmd.label("label_threshold_13.6499996185", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.6499996185', members= 'label_threshold_13.6499996185')


if dirpath:
    f = join(dirpath, '87/mesh.grd')
else:
    f = '87/mesh.grd'
cmd.load(f, 'mesh_87')
cmd.isomesh("isomesh_87", "mesh_87", 0.9)
cmd.color("grey80", "isomesh_87")
cmd.set('transparency', 0.4, "isomesh_87")

cmd.group('hotspot_87', "isomesh_87")
cmd.group('hotspot_87', "mesh_87")

if dirpath:
    f = join(dirpath, "88/label_threshold_15.0.mol2")
else:
    f = "88/label_threshold_15.0.mol2"

cmd.load(f, 'label_threshold_15.0')
cmd.hide('everything', 'label_threshold_15.0')
cmd.label("label_threshold_15.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.0]
gfiles = ['88/donor.grd', '88/apolar.grd', '88/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 88
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


cluster_dict = {"13.5979995728":[], "13.5979995728_arrows":[]}

cluster_dict["13.5979995728"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(15.5220735046), float(-20.6223570596), float(41.4891414415), float(1.0)]


cluster_dict["13.5979995728"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(20.6028925057), float(-12.7138436582), float(38.2258240654), float(1.0)]


cluster_dict["13.5979995728"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(10.5), float(-19.0), float(41.5), float(1.0)]

cluster_dict["13.5979995728_arrows"] += cgo_arrow([10.5,-19.0,41.5], [8.555,-21.148,42.25], color="red blue", name="Arrows_13.5979995728_1")

cmd.load_cgo(cluster_dict["13.5979995728"], "Features_13.5979995728", 1)
cmd.load_cgo(cluster_dict["13.5979995728_arrows"], "Arrows_13.5979995728")
cmd.set("transparency", 0.2,"Features_13.5979995728")
cmd.group("Pharmacophore_13.5979995728", members="Features_13.5979995728")
cmd.group("Pharmacophore_13.5979995728", members="Arrows_13.5979995728")

if dirpath:
    f = join(dirpath, "88/label_threshold_13.5979995728.mol2")
else:
    f = "88/label_threshold_13.5979995728.mol2"

cmd.load(f, 'label_threshold_13.5979995728')
cmd.hide('everything', 'label_threshold_13.5979995728')
cmd.label("label_threshold_13.5979995728", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.5979995728', members= 'label_threshold_13.5979995728')


if dirpath:
    f = join(dirpath, '88/mesh.grd')
else:
    f = '88/mesh.grd'
cmd.load(f, 'mesh_88')
cmd.isomesh("isomesh_88", "mesh_88", 0.9)
cmd.color("grey80", "isomesh_88")
cmd.set('transparency', 0.4, "isomesh_88")

cmd.group('hotspot_88', "isomesh_88")
cmd.group('hotspot_88', "mesh_88")

if dirpath:
    f = join(dirpath, "89/label_threshold_11.9.mol2")
else:
    f = "89/label_threshold_11.9.mol2"

cmd.load(f, 'label_threshold_11.9')
cmd.hide('everything', 'label_threshold_11.9')
cmd.label("label_threshold_11.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [11.9]
gfiles = ['89/donor.grd', '89/apolar.grd', '89/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 89
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


cluster_dict = {"13.5570001602":[], "13.5570001602_arrows":[]}

cluster_dict["13.5570001602"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-11.0), float(8.5), float(30.5), float(1.0)]

cluster_dict["13.5570001602_arrows"] += cgo_arrow([-11.0,8.5,30.5], [-12.682,10.491,30.534], color="blue red", name="Arrows_13.5570001602_1")

cluster_dict["13.5570001602"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-8.50413955772), float(12.2589652309), float(28.6219068985), float(1.0)]


cluster_dict["13.5570001602"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-10.3910607962), float(8.31459573776), float(32.6452522717), float(1.0)]


cluster_dict["13.5570001602"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-8.32572708714), float(12.0406848078), float(38.4137080344), float(1.0)]


cluster_dict["13.5570001602"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-5.93172994162), float(19.9786139874), float(35.0098352585), float(1.0)]


cluster_dict["13.5570001602"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-6.88535031847), float(14.4872611465), float(29.9808917197), float(1.0)]


cluster_dict["13.5570001602"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-0.0151718870713), float(12.9404145944), float(43.999253171), float(1.0)]


cluster_dict["13.5570001602"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-6.5), float(18.0), float(32.5), float(1.0)]


cluster_dict["13.5570001602"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-4.5), float(8.0), float(30.0), float(1.0)]

cluster_dict["13.5570001602_arrows"] += cgo_arrow([-4.5,8.0,30.0], [-6.661,9.212,32.894], color="red blue", name="Arrows_13.5570001602_2")

cmd.load_cgo(cluster_dict["13.5570001602"], "Features_13.5570001602", 1)
cmd.load_cgo(cluster_dict["13.5570001602_arrows"], "Arrows_13.5570001602")
cmd.set("transparency", 0.2,"Features_13.5570001602")
cmd.group("Pharmacophore_13.5570001602", members="Features_13.5570001602")
cmd.group("Pharmacophore_13.5570001602", members="Arrows_13.5570001602")

if dirpath:
    f = join(dirpath, "89/label_threshold_13.5570001602.mol2")
else:
    f = "89/label_threshold_13.5570001602.mol2"

cmd.load(f, 'label_threshold_13.5570001602')
cmd.hide('everything', 'label_threshold_13.5570001602')
cmd.label("label_threshold_13.5570001602", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.5570001602', members= 'label_threshold_13.5570001602')


if dirpath:
    f = join(dirpath, '89/mesh.grd')
else:
    f = '89/mesh.grd'
cmd.load(f, 'mesh_89')
cmd.isomesh("isomesh_89", "mesh_89", 0.9)
cmd.color("grey80", "isomesh_89")
cmd.set('transparency', 0.4, "isomesh_89")

cmd.group('hotspot_89', "isomesh_89")
cmd.group('hotspot_89', "mesh_89")

if dirpath:
    f = join(dirpath, "90/label_threshold_10.2.mol2")
else:
    f = "90/label_threshold_10.2.mol2"

cmd.load(f, 'label_threshold_10.2')
cmd.hide('everything', 'label_threshold_10.2')
cmd.label("label_threshold_10.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10.2]
gfiles = ['90/donor.grd', '90/apolar.grd', '90/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 90
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


cluster_dict = {"13.4659996033":[], "13.4659996033_arrows":[]}

cluster_dict["13.4659996033"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-19.5), float(-6.5), float(33.0), float(1.0)]

cluster_dict["13.4659996033_arrows"] += cgo_arrow([-19.5,-6.5,33.0], [-20.183,-8.769,31.779], color="blue red", name="Arrows_13.4659996033_1")

cluster_dict["13.4659996033"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-17.5), float(-10.0), float(28.0), float(1.0)]

cluster_dict["13.4659996033_arrows"] += cgo_arrow([-17.5,-10.0,28.0], [-19.572,-9.949,26.024], color="blue red", name="Arrows_13.4659996033_2")

cluster_dict["13.4659996033"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-13.5), float(-6.0), float(23.0), float(1.0)]

cluster_dict["13.4659996033_arrows"] += cgo_arrow([-13.5,-6.0,23.0], [-15.302,-8.206,23.93], color="blue red", name="Arrows_13.4659996033_3")

cluster_dict["13.4659996033"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-14.0), float(-2.0), float(25.0), float(1.0)]

cluster_dict["13.4659996033_arrows"] += cgo_arrow([-14.0,-2.0,25.0], [-14.141,-2.336,27.942], color="blue red", name="Arrows_13.4659996033_4")

cluster_dict["13.4659996033"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-8.5), float(2.0), float(19.5), float(1.0)]

cluster_dict["13.4659996033_arrows"] += cgo_arrow([-8.5,2.0,19.5], [-10.06,1.446,18.913], color="blue red", name="Arrows_13.4659996033_5")

cluster_dict["13.4659996033"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-14.7717843313), float(8.93850836527), float(21.8753066901), float(1.0)]


cluster_dict["13.4659996033"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-14.1872328545), float(-3.53099565771), float(22.3341834313), float(1.0)]


cluster_dict["13.4659996033"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-14.1142857143), float(-2.15714285714), float(33.3428571429), float(1.0)]


cluster_dict["13.4659996033"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-9.17845819463), float(2.03565045652), float(28.4798427785), float(1.0)]


cluster_dict["13.4659996033"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-6.43308383127), float(0.343236179665), float(18.2977611508), float(1.0)]


cluster_dict["13.4659996033"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-5.87822968794), float(-2.41174750151), float(27.6207561883), float(1.0)]


cluster_dict["13.4659996033"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-16.5), float(-10.0), float(27.0), float(1.0)]

cluster_dict["13.4659996033_arrows"] += cgo_arrow([-16.5,-10.0,27.0], [-19.572,-9.949,26.024], color="red blue", name="Arrows_13.4659996033_6")

cluster_dict["13.4659996033"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-16.5), float(8.5), float(23.0), float(1.0)]

cluster_dict["13.4659996033_arrows"] += cgo_arrow([-16.5,8.5,23.0], [-13.826,4.855,21.317], color="red blue", name="Arrows_13.4659996033_7")

cluster_dict["13.4659996033"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-10.5), float(-2.0), float(27.5), float(1.0)]

cluster_dict["13.4659996033_arrows"] += cgo_arrow([-10.5,-2.0,27.5], [-12.063,-0.184,28.703], color="red blue", name="Arrows_13.4659996033_8")

cluster_dict["13.4659996033"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-9.5), float(1.5), float(30.0), float(1.0)]

cluster_dict["13.4659996033_arrows"] += cgo_arrow([-9.5,1.5,30.0], [-12.063,-0.184,28.703], color="red blue", name="Arrows_13.4659996033_9")

cluster_dict["13.4659996033"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-8.0), float(-2.5), float(28.0), float(1.0)]

cluster_dict["13.4659996033_arrows"] += cgo_arrow([-8.0,-2.5,28.0], [-9.153,-5.211,25.883], color="red blue", name="Arrows_13.4659996033_10")

cmd.load_cgo(cluster_dict["13.4659996033"], "Features_13.4659996033", 1)
cmd.load_cgo(cluster_dict["13.4659996033_arrows"], "Arrows_13.4659996033")
cmd.set("transparency", 0.2,"Features_13.4659996033")
cmd.group("Pharmacophore_13.4659996033", members="Features_13.4659996033")
cmd.group("Pharmacophore_13.4659996033", members="Arrows_13.4659996033")

if dirpath:
    f = join(dirpath, "90/label_threshold_13.4659996033.mol2")
else:
    f = "90/label_threshold_13.4659996033.mol2"

cmd.load(f, 'label_threshold_13.4659996033')
cmd.hide('everything', 'label_threshold_13.4659996033')
cmd.label("label_threshold_13.4659996033", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.4659996033', members= 'label_threshold_13.4659996033')


if dirpath:
    f = join(dirpath, '90/mesh.grd')
else:
    f = '90/mesh.grd'
cmd.load(f, 'mesh_90')
cmd.isomesh("isomesh_90", "mesh_90", 0.9)
cmd.color("grey80", "isomesh_90")
cmd.set('transparency', 0.4, "isomesh_90")

cmd.group('hotspot_90', "isomesh_90")
cmd.group('hotspot_90', "mesh_90")

if dirpath:
    f = join(dirpath, "91/label_threshold_11.9.mol2")
else:
    f = "91/label_threshold_11.9.mol2"

cmd.load(f, 'label_threshold_11.9')
cmd.hide('everything', 'label_threshold_11.9')
cmd.label("label_threshold_11.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [11.9]
gfiles = ['91/donor.grd', '91/apolar.grd', '91/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 91
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


cluster_dict = {"13.2349996567":[], "13.2349996567_arrows":[]}

cluster_dict["13.2349996567"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-16.0), float(-3.0), float(54.0), float(1.0)]

cluster_dict["13.2349996567_arrows"] += cgo_arrow([-16.0,-3.0,54.0], [-17.363,-0.587,52.277], color="blue red", name="Arrows_13.2349996567_1")

cluster_dict["13.2349996567"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-16.0), float(-4.5), float(48.5), float(1.0)]

cluster_dict["13.2349996567_arrows"] += cgo_arrow([-16.0,-4.5,48.5], [-18.731,-2.721,49.33], color="blue red", name="Arrows_13.2349996567_2")

cluster_dict["13.2349996567"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-14.0), float(-4.0), float(46.5), float(1.0)]

cluster_dict["13.2349996567_arrows"] += cgo_arrow([-14.0,-4.0,46.5], [-13.027,-1.485,45.573], color="blue red", name="Arrows_13.2349996567_3")

cluster_dict["13.2349996567"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-12.5), float(-11.5), float(49.0), float(1.0)]

cluster_dict["13.2349996567_arrows"] += cgo_arrow([-12.5,-11.5,49.0], [-9.173,-12.246,47.746], color="blue red", name="Arrows_13.2349996567_4")

cluster_dict["13.2349996567"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-16.440752504), float(-8.06355100462), float(47.5568299782), float(1.0)]


cluster_dict["13.2349996567"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-17.6646354123), float(-13.2677194032), float(55.9536277695), float(1.0)]


cluster_dict["13.2349996567"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-16.0752768125), float(-2.40509635209), float(56.9283773099), float(1.0)]


cluster_dict["13.2349996567"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-14.4861385615), float(-0.639112063669), float(49.530595123), float(1.0)]


cluster_dict["13.2349996567"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-21.0), float(-6.0), float(45.5), float(1.0)]

cluster_dict["13.2349996567_arrows"] += cgo_arrow([-21.0,-6.0,45.5], [-22.337,-3.396,48.966], color="red blue", name="Arrows_13.2349996567_5")

cluster_dict["13.2349996567"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-19.0), float(-4.0), float(56.5), float(1.0)]

cluster_dict["13.2349996567_arrows"] += cgo_arrow([-19.0,-4.0,56.5], [-21.532,-6.13,56.958], color="red blue", name="Arrows_13.2349996567_6")

cluster_dict["13.2349996567"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-18.0), float(-2.0), float(57.5), float(1.0)]


cluster_dict["13.2349996567"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-17.0), float(-3.5), float(55.5), float(1.0)]


cluster_dict["13.2349996567"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-15.0), float(-0.5), float(46.5), float(1.0)]

cluster_dict["13.2349996567_arrows"] += cgo_arrow([-15.0,-0.5,46.5], [-17.71,0.005,45.762], color="red blue", name="Arrows_13.2349996567_7")

cmd.load_cgo(cluster_dict["13.2349996567"], "Features_13.2349996567", 1)
cmd.load_cgo(cluster_dict["13.2349996567_arrows"], "Arrows_13.2349996567")
cmd.set("transparency", 0.2,"Features_13.2349996567")
cmd.group("Pharmacophore_13.2349996567", members="Features_13.2349996567")
cmd.group("Pharmacophore_13.2349996567", members="Arrows_13.2349996567")

if dirpath:
    f = join(dirpath, "91/label_threshold_13.2349996567.mol2")
else:
    f = "91/label_threshold_13.2349996567.mol2"

cmd.load(f, 'label_threshold_13.2349996567')
cmd.hide('everything', 'label_threshold_13.2349996567')
cmd.label("label_threshold_13.2349996567", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.2349996567', members= 'label_threshold_13.2349996567')


if dirpath:
    f = join(dirpath, '91/mesh.grd')
else:
    f = '91/mesh.grd'
cmd.load(f, 'mesh_91')
cmd.isomesh("isomesh_91", "mesh_91", 0.9)
cmd.color("grey80", "isomesh_91")
cmd.set('transparency', 0.4, "isomesh_91")

cmd.group('hotspot_91', "isomesh_91")
cmd.group('hotspot_91', "mesh_91")

if dirpath:
    f = join(dirpath, "92/label_threshold_2.7.mol2")
else:
    f = "92/label_threshold_2.7.mol2"

cmd.load(f, 'label_threshold_2.7')
cmd.hide('everything', 'label_threshold_2.7')
cmd.label("label_threshold_2.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [2.7]
gfiles = ['92/donor.grd', '92/apolar.grd', '92/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 92
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


cluster_dict = {"12.998000145":[], "12.998000145_arrows":[]}

cluster_dict["12.998000145"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-50.733689033), float(20.8710615314), float(30.7644406774), float(1.0)]


cmd.load_cgo(cluster_dict["12.998000145"], "Features_12.998000145", 1)
cmd.load_cgo(cluster_dict["12.998000145_arrows"], "Arrows_12.998000145")
cmd.set("transparency", 0.2,"Features_12.998000145")
cmd.group("Pharmacophore_12.998000145", members="Features_12.998000145")
cmd.group("Pharmacophore_12.998000145", members="Arrows_12.998000145")

if dirpath:
    f = join(dirpath, "92/label_threshold_12.998000145.mol2")
else:
    f = "92/label_threshold_12.998000145.mol2"

cmd.load(f, 'label_threshold_12.998000145')
cmd.hide('everything', 'label_threshold_12.998000145')
cmd.label("label_threshold_12.998000145", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.998000145', members= 'label_threshold_12.998000145')


if dirpath:
    f = join(dirpath, '92/mesh.grd')
else:
    f = '92/mesh.grd'
cmd.load(f, 'mesh_92')
cmd.isomesh("isomesh_92", "mesh_92", 0.9)
cmd.color("grey80", "isomesh_92")
cmd.set('transparency', 0.4, "isomesh_92")

cmd.group('hotspot_92', "isomesh_92")
cmd.group('hotspot_92', "mesh_92")

if dirpath:
    f = join(dirpath, "93/label_threshold_2.7.mol2")
else:
    f = "93/label_threshold_2.7.mol2"

cmd.load(f, 'label_threshold_2.7')
cmd.hide('everything', 'label_threshold_2.7')
cmd.label("label_threshold_2.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [2.7]
gfiles = ['93/donor.grd', '93/apolar.grd', '93/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 93
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


cluster_dict = {"12.7760000229":[], "12.7760000229_arrows":[]}

cluster_dict["12.7760000229"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(28.4871425416), float(36.2365798492), float(42.9958944135), float(1.0)]


cluster_dict["12.7760000229"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(34.0), float(43.0), float(1.0)]

cluster_dict["12.7760000229_arrows"] += cgo_arrow([28.0,34.0,43.0], [30.922,33.091,43.322], color="red blue", name="Arrows_12.7760000229_1")

cmd.load_cgo(cluster_dict["12.7760000229"], "Features_12.7760000229", 1)
cmd.load_cgo(cluster_dict["12.7760000229_arrows"], "Arrows_12.7760000229")
cmd.set("transparency", 0.2,"Features_12.7760000229")
cmd.group("Pharmacophore_12.7760000229", members="Features_12.7760000229")
cmd.group("Pharmacophore_12.7760000229", members="Arrows_12.7760000229")

if dirpath:
    f = join(dirpath, "93/label_threshold_12.7760000229.mol2")
else:
    f = "93/label_threshold_12.7760000229.mol2"

cmd.load(f, 'label_threshold_12.7760000229')
cmd.hide('everything', 'label_threshold_12.7760000229')
cmd.label("label_threshold_12.7760000229", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.7760000229', members= 'label_threshold_12.7760000229')


if dirpath:
    f = join(dirpath, '93/mesh.grd')
else:
    f = '93/mesh.grd'
cmd.load(f, 'mesh_93')
cmd.isomesh("isomesh_93", "mesh_93", 0.9)
cmd.color("grey80", "isomesh_93")
cmd.set('transparency', 0.4, "isomesh_93")

cmd.group('hotspot_93', "isomesh_93")
cmd.group('hotspot_93', "mesh_93")

if dirpath:
    f = join(dirpath, "94/label_threshold_9.4.mol2")
else:
    f = "94/label_threshold_9.4.mol2"

cmd.load(f, 'label_threshold_9.4')
cmd.hide('everything', 'label_threshold_9.4')
cmd.label("label_threshold_9.4", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [9.4]
gfiles = ['94/donor.grd', '94/apolar.grd', '94/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 94
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


cluster_dict = {"12.2770004272":[], "12.2770004272_arrows":[]}

cluster_dict["12.2770004272"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(34.5), float(9.5), float(24.5), float(1.0)]

cluster_dict["12.2770004272_arrows"] += cgo_arrow([34.5,9.5,24.5], [34.779,7.271,23.0], color="blue red", name="Arrows_12.2770004272_1")

cluster_dict["12.2770004272"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(35.1193397295), float(10.6958269607), float(19.7163648091), float(1.0)]


cluster_dict["12.2770004272"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(44.8725482561), float(13.8419708239), float(19.4295354797), float(1.0)]


cluster_dict["12.2770004272"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(43.5), float(12.0), float(18.0), float(1.0)]

cluster_dict["12.2770004272_arrows"] += cgo_arrow([43.5,12.0,18.0], [41.724,10.524,16.554], color="red blue", name="Arrows_12.2770004272_2")

cmd.load_cgo(cluster_dict["12.2770004272"], "Features_12.2770004272", 1)
cmd.load_cgo(cluster_dict["12.2770004272_arrows"], "Arrows_12.2770004272")
cmd.set("transparency", 0.2,"Features_12.2770004272")
cmd.group("Pharmacophore_12.2770004272", members="Features_12.2770004272")
cmd.group("Pharmacophore_12.2770004272", members="Arrows_12.2770004272")

if dirpath:
    f = join(dirpath, "94/label_threshold_12.2770004272.mol2")
else:
    f = "94/label_threshold_12.2770004272.mol2"

cmd.load(f, 'label_threshold_12.2770004272')
cmd.hide('everything', 'label_threshold_12.2770004272')
cmd.label("label_threshold_12.2770004272", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.2770004272', members= 'label_threshold_12.2770004272')


if dirpath:
    f = join(dirpath, '94/mesh.grd')
else:
    f = '94/mesh.grd'
cmd.load(f, 'mesh_94')
cmd.isomesh("isomesh_94", "mesh_94", 0.9)
cmd.color("grey80", "isomesh_94")
cmd.set('transparency', 0.4, "isomesh_94")

cmd.group('hotspot_94', "isomesh_94")
cmd.group('hotspot_94', "mesh_94")

if dirpath:
    f = join(dirpath, "95/label_threshold_8.1.mol2")
else:
    f = "95/label_threshold_8.1.mol2"

cmd.load(f, 'label_threshold_8.1')
cmd.hide('everything', 'label_threshold_8.1')
cmd.label("label_threshold_8.1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [8.1]
gfiles = ['95/donor.grd', '95/apolar.grd', '95/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 95
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


cluster_dict = {"12.140999794":[], "12.140999794_arrows":[]}

cluster_dict["12.140999794"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(9.5), float(-22.5), float(-1.5), float(1.0)]

cluster_dict["12.140999794_arrows"] += cgo_arrow([9.5,-22.5,-1.5], [10.264,-25.83,-1.724], color="blue red", name="Arrows_12.140999794_1")

cluster_dict["12.140999794"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(9.5), float(-20.0), float(0.0), float(1.0)]

cluster_dict["12.140999794_arrows"] += cgo_arrow([9.5,-20.0,0.0], [7.395,-17.828,1.205], color="blue red", name="Arrows_12.140999794_2")

cluster_dict["12.140999794"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(9.5), float(-20.0), float(1.5), float(1.0)]

cluster_dict["12.140999794_arrows"] += cgo_arrow([9.5,-20.0,1.5], [7.395,-17.828,1.205], color="blue red", name="Arrows_12.140999794_3")

cluster_dict["12.140999794"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(12.0), float(-15.5), float(-1.0), float(1.0)]

cluster_dict["12.140999794_arrows"] += cgo_arrow([12.0,-15.5,-1.0], [13.927,-14.08,0.735], color="blue red", name="Arrows_12.140999794_4")

cluster_dict["12.140999794"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(10.5427120886), float(-20.6071672345), float(-0.65659465427), float(1.0)]


cluster_dict["12.140999794"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(8.5), float(-18.5), float(-2.5), float(1.0)]

cluster_dict["12.140999794_arrows"] += cgo_arrow([8.5,-18.5,-2.5], [5.591,-17.952,-1.26], color="red blue", name="Arrows_12.140999794_5")

cmd.load_cgo(cluster_dict["12.140999794"], "Features_12.140999794", 1)
cmd.load_cgo(cluster_dict["12.140999794_arrows"], "Arrows_12.140999794")
cmd.set("transparency", 0.2,"Features_12.140999794")
cmd.group("Pharmacophore_12.140999794", members="Features_12.140999794")
cmd.group("Pharmacophore_12.140999794", members="Arrows_12.140999794")

if dirpath:
    f = join(dirpath, "95/label_threshold_12.140999794.mol2")
else:
    f = "95/label_threshold_12.140999794.mol2"

cmd.load(f, 'label_threshold_12.140999794')
cmd.hide('everything', 'label_threshold_12.140999794')
cmd.label("label_threshold_12.140999794", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.140999794', members= 'label_threshold_12.140999794')


if dirpath:
    f = join(dirpath, '95/mesh.grd')
else:
    f = '95/mesh.grd'
cmd.load(f, 'mesh_95')
cmd.isomesh("isomesh_95", "mesh_95", 0.9)
cmd.color("grey80", "isomesh_95")
cmd.set('transparency', 0.4, "isomesh_95")

cmd.group('hotspot_95', "isomesh_95")
cmd.group('hotspot_95', "mesh_95")

if dirpath:
    f = join(dirpath, "96/label_threshold_2.0.mol2")
else:
    f = "96/label_threshold_2.0.mol2"

cmd.load(f, 'label_threshold_2.0')
cmd.hide('everything', 'label_threshold_2.0')
cmd.label("label_threshold_2.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [2.0]
gfiles = ['96/donor.grd', '96/apolar.grd', '96/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 96
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


cluster_dict = {"12.0440001488":[], "12.0440001488_arrows":[]}

cluster_dict["12.0440001488"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(39.2327483685), float(-9.12699809355), float(25.6581805215), float(1.0)]


cmd.load_cgo(cluster_dict["12.0440001488"], "Features_12.0440001488", 1)
cmd.load_cgo(cluster_dict["12.0440001488_arrows"], "Arrows_12.0440001488")
cmd.set("transparency", 0.2,"Features_12.0440001488")
cmd.group("Pharmacophore_12.0440001488", members="Features_12.0440001488")
cmd.group("Pharmacophore_12.0440001488", members="Arrows_12.0440001488")

if dirpath:
    f = join(dirpath, "96/label_threshold_12.0440001488.mol2")
else:
    f = "96/label_threshold_12.0440001488.mol2"

cmd.load(f, 'label_threshold_12.0440001488')
cmd.hide('everything', 'label_threshold_12.0440001488')
cmd.label("label_threshold_12.0440001488", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_12.0440001488', members= 'label_threshold_12.0440001488')


if dirpath:
    f = join(dirpath, '96/mesh.grd')
else:
    f = '96/mesh.grd'
cmd.load(f, 'mesh_96')
cmd.isomesh("isomesh_96", "mesh_96", 0.9)
cmd.color("grey80", "isomesh_96")
cmd.set('transparency', 0.4, "isomesh_96")

cmd.group('hotspot_96', "isomesh_96")
cmd.group('hotspot_96', "mesh_96")

if dirpath:
    f = join(dirpath, "97/label_threshold_9.0.mol2")
else:
    f = "97/label_threshold_9.0.mol2"

cmd.load(f, 'label_threshold_9.0')
cmd.hide('everything', 'label_threshold_9.0')
cmd.label("label_threshold_9.0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [9.0]
gfiles = ['97/donor.grd', '97/apolar.grd', '97/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 97
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


cluster_dict = {"11.8369998932":[], "11.8369998932_arrows":[]}

cluster_dict["11.8369998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(32.5), float(14.5), float(24.5), float(1.0)]

cluster_dict["11.8369998932_arrows"] += cgo_arrow([32.5,14.5,24.5], [32.987,16.799,26.549], color="blue red", name="Arrows_11.8369998932_1")

cluster_dict["11.8369998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(35.0), float(10.0), float(25.0), float(1.0)]

cluster_dict["11.8369998932_arrows"] += cgo_arrow([35.0,10.0,25.0], [36.151,9.316,27.38], color="blue red", name="Arrows_11.8369998932_2")

cluster_dict["11.8369998932"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(31.1583439847), float(21.3713173324), float(21.7117958105), float(1.0)]


cluster_dict["11.8369998932"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(34.6546866306), float(12.5042200797), float(21.8145304399), float(1.0)]


cluster_dict["11.8369998932"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(44.6922539754), float(14.0369648614), float(19.5048898485), float(1.0)]


cluster_dict["11.8369998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(43.5), float(12.0), float(18.0), float(1.0)]

cluster_dict["11.8369998932_arrows"] += cgo_arrow([43.5,12.0,18.0], [41.724,10.524,16.554], color="red blue", name="Arrows_11.8369998932_3")

cmd.load_cgo(cluster_dict["11.8369998932"], "Features_11.8369998932", 1)
cmd.load_cgo(cluster_dict["11.8369998932_arrows"], "Arrows_11.8369998932")
cmd.set("transparency", 0.2,"Features_11.8369998932")
cmd.group("Pharmacophore_11.8369998932", members="Features_11.8369998932")
cmd.group("Pharmacophore_11.8369998932", members="Arrows_11.8369998932")

if dirpath:
    f = join(dirpath, "97/label_threshold_11.8369998932.mol2")
else:
    f = "97/label_threshold_11.8369998932.mol2"

cmd.load(f, 'label_threshold_11.8369998932')
cmd.hide('everything', 'label_threshold_11.8369998932')
cmd.label("label_threshold_11.8369998932", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.8369998932', members= 'label_threshold_11.8369998932')


if dirpath:
    f = join(dirpath, '97/mesh.grd')
else:
    f = '97/mesh.grd'
cmd.load(f, 'mesh_97')
cmd.isomesh("isomesh_97", "mesh_97", 0.9)
cmd.color("grey80", "isomesh_97")
cmd.set('transparency', 0.4, "isomesh_97")

cmd.group('hotspot_97', "isomesh_97")
cmd.group('hotspot_97', "mesh_97")

if dirpath:
    f = join(dirpath, "98/label_threshold_6.7.mol2")
else:
    f = "98/label_threshold_6.7.mol2"

cmd.load(f, 'label_threshold_6.7')
cmd.hide('everything', 'label_threshold_6.7')
cmd.label("label_threshold_6.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [6.7]
gfiles = ['98/donor.grd', '98/apolar.grd', '98/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 98
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


cluster_dict = {"11.6494998932":[], "11.6494998932_arrows":[]}

cluster_dict["11.6494998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(9.5), float(-22.5), float(-1.5), float(1.0)]

cluster_dict["11.6494998932_arrows"] += cgo_arrow([9.5,-22.5,-1.5], [10.264,-25.83,-1.724], color="blue red", name="Arrows_11.6494998932_1")

cluster_dict["11.6494998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(9.5), float(-20.0), float(0.0), float(1.0)]

cluster_dict["11.6494998932_arrows"] += cgo_arrow([9.5,-20.0,0.0], [7.395,-17.828,1.205], color="blue red", name="Arrows_11.6494998932_2")

cluster_dict["11.6494998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(9.5), float(-20.0), float(1.5), float(1.0)]

cluster_dict["11.6494998932_arrows"] += cgo_arrow([9.5,-20.0,1.5], [7.395,-17.828,1.205], color="blue red", name="Arrows_11.6494998932_3")

cluster_dict["11.6494998932"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(11.5), float(-16.0), float(-1.0), float(1.0)]

cluster_dict["11.6494998932_arrows"] += cgo_arrow([11.5,-16.0,-1.0], [12.41,-16.337,2.313], color="blue red", name="Arrows_11.6494998932_4")

cluster_dict["11.6494998932"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(10.1570070428), float(-21.7182986006), float(-0.0211001501216), float(1.0)]


cluster_dict["11.6494998932"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(8.5), float(-18.5), float(-2.5), float(1.0)]

cluster_dict["11.6494998932_arrows"] += cgo_arrow([8.5,-18.5,-2.5], [5.591,-17.952,-1.26], color="red blue", name="Arrows_11.6494998932_5")

cmd.load_cgo(cluster_dict["11.6494998932"], "Features_11.6494998932", 1)
cmd.load_cgo(cluster_dict["11.6494998932_arrows"], "Arrows_11.6494998932")
cmd.set("transparency", 0.2,"Features_11.6494998932")
cmd.group("Pharmacophore_11.6494998932", members="Features_11.6494998932")
cmd.group("Pharmacophore_11.6494998932", members="Arrows_11.6494998932")

if dirpath:
    f = join(dirpath, "98/label_threshold_11.6494998932.mol2")
else:
    f = "98/label_threshold_11.6494998932.mol2"

cmd.load(f, 'label_threshold_11.6494998932')
cmd.hide('everything', 'label_threshold_11.6494998932')
cmd.label("label_threshold_11.6494998932", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.6494998932', members= 'label_threshold_11.6494998932')


if dirpath:
    f = join(dirpath, '98/mesh.grd')
else:
    f = '98/mesh.grd'
cmd.load(f, 'mesh_98')
cmd.isomesh("isomesh_98", "mesh_98", 0.9)
cmd.color("grey80", "isomesh_98")
cmd.set('transparency', 0.4, "isomesh_98")

cmd.group('hotspot_98', "isomesh_98")
cmd.group('hotspot_98', "mesh_98")

if dirpath:
    f = join(dirpath, "99/label_threshold_15.8.mol2")
else:
    f = "99/label_threshold_15.8.mol2"

cmd.load(f, 'label_threshold_15.8')
cmd.hide('everything', 'label_threshold_15.8')
cmd.label("label_threshold_15.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.8]
gfiles = ['99/donor.grd', '99/apolar.grd', '99/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 99
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


cluster_dict = {"11.4879999161":[], "11.4879999161_arrows":[]}

cluster_dict["11.4879999161"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-70.5242772296), float(-16.6360862755), float(36.9788212062), float(1.0)]


cmd.load_cgo(cluster_dict["11.4879999161"], "Features_11.4879999161", 1)
cmd.load_cgo(cluster_dict["11.4879999161_arrows"], "Arrows_11.4879999161")
cmd.set("transparency", 0.2,"Features_11.4879999161")
cmd.group("Pharmacophore_11.4879999161", members="Features_11.4879999161")
cmd.group("Pharmacophore_11.4879999161", members="Arrows_11.4879999161")

if dirpath:
    f = join(dirpath, "99/label_threshold_11.4879999161.mol2")
else:
    f = "99/label_threshold_11.4879999161.mol2"

cmd.load(f, 'label_threshold_11.4879999161')
cmd.hide('everything', 'label_threshold_11.4879999161')
cmd.label("label_threshold_11.4879999161", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.4879999161', members= 'label_threshold_11.4879999161')


if dirpath:
    f = join(dirpath, '99/mesh.grd')
else:
    f = '99/mesh.grd'
cmd.load(f, 'mesh_99')
cmd.isomesh("isomesh_99", "mesh_99", 0.9)
cmd.color("grey80", "isomesh_99")
cmd.set('transparency', 0.4, "isomesh_99")

cmd.group('hotspot_99', "isomesh_99")
cmd.group('hotspot_99', "mesh_99")

if dirpath:
    f = join(dirpath, "100/label_threshold_15.6.mol2")
else:
    f = "100/label_threshold_15.6.mol2"

cmd.load(f, 'label_threshold_15.6')
cmd.hide('everything', 'label_threshold_15.6')
cmd.label("label_threshold_15.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.6]
gfiles = ['100/donor.grd', '100/apolar.grd', '100/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 100
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


cluster_dict = {"11.4610004425":[], "11.4610004425_arrows":[]}

cluster_dict["11.4610004425"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-70.5590119742), float(-16.6506665606), float(37.0017047379), float(1.0)]


cmd.load_cgo(cluster_dict["11.4610004425"], "Features_11.4610004425", 1)
cmd.load_cgo(cluster_dict["11.4610004425_arrows"], "Arrows_11.4610004425")
cmd.set("transparency", 0.2,"Features_11.4610004425")
cmd.group("Pharmacophore_11.4610004425", members="Features_11.4610004425")
cmd.group("Pharmacophore_11.4610004425", members="Arrows_11.4610004425")

if dirpath:
    f = join(dirpath, "100/label_threshold_11.4610004425.mol2")
else:
    f = "100/label_threshold_11.4610004425.mol2"

cmd.load(f, 'label_threshold_11.4610004425')
cmd.hide('everything', 'label_threshold_11.4610004425')
cmd.label("label_threshold_11.4610004425", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.4610004425', members= 'label_threshold_11.4610004425')


if dirpath:
    f = join(dirpath, '100/mesh.grd')
else:
    f = '100/mesh.grd'
cmd.load(f, 'mesh_100')
cmd.isomesh("isomesh_100", "mesh_100", 0.9)
cmd.color("grey80", "isomesh_100")
cmd.set('transparency', 0.4, "isomesh_100")

cmd.group('hotspot_100', "isomesh_100")
cmd.group('hotspot_100', "mesh_100")

if dirpath:
    f = join(dirpath, "101/label_threshold_0.6.mol2")
else:
    f = "101/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
gfiles = ['101/donor.grd', '101/apolar.grd', '101/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 101
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


cluster_dict = {"11.3070001602":[], "11.3070001602_arrows":[]}

cluster_dict["11.3070001602"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(28.5), float(-21.0), float(-3.0), float(1.0)]

cluster_dict["11.3070001602_arrows"] += cgo_arrow([28.5,-21.0,-3.0], [29.875,-18.848,-1.33], color="blue red", name="Arrows_11.3070001602_1")

cluster_dict["11.3070001602"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(27.0523116478), float(-19.9085537021), float(-4.53149552309), float(1.0)]


cmd.load_cgo(cluster_dict["11.3070001602"], "Features_11.3070001602", 1)
cmd.load_cgo(cluster_dict["11.3070001602_arrows"], "Arrows_11.3070001602")
cmd.set("transparency", 0.2,"Features_11.3070001602")
cmd.group("Pharmacophore_11.3070001602", members="Features_11.3070001602")
cmd.group("Pharmacophore_11.3070001602", members="Arrows_11.3070001602")

if dirpath:
    f = join(dirpath, "101/label_threshold_11.3070001602.mol2")
else:
    f = "101/label_threshold_11.3070001602.mol2"

cmd.load(f, 'label_threshold_11.3070001602')
cmd.hide('everything', 'label_threshold_11.3070001602')
cmd.label("label_threshold_11.3070001602", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_11.3070001602', members= 'label_threshold_11.3070001602')


if dirpath:
    f = join(dirpath, '101/mesh.grd')
else:
    f = '101/mesh.grd'
cmd.load(f, 'mesh_101')
cmd.isomesh("isomesh_101", "mesh_101", 0.9)
cmd.color("grey80", "isomesh_101")
cmd.set('transparency', 0.4, "isomesh_101")

cmd.group('hotspot_101', "isomesh_101")
cmd.group('hotspot_101', "mesh_101")

if dirpath:
    f = join(dirpath, "102/label_threshold_0.6.mol2")
else:
    f = "102/label_threshold_0.6.mol2"

cmd.load(f, 'label_threshold_0.6')
cmd.hide('everything', 'label_threshold_0.6')
cmd.label("label_threshold_0.6", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.6]
gfiles = ['102/donor.grd', '102/apolar.grd', '102/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 102
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


cluster_dict = {"10.4029998779":[], "10.4029998779_arrows":[]}

cluster_dict["10.4029998779"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-17.5), float(-35.5), float(27.5), float(1.0)]

cluster_dict["10.4029998779_arrows"] += cgo_arrow([-17.5,-35.5,27.5], [-19.935,-34.025,26.449], color="blue red", name="Arrows_10.4029998779_1")

cluster_dict["10.4029998779"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-18.1721581031), float(-37.5913142944), float(25.6828604984), float(1.0)]


cluster_dict["10.4029998779"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-19.0), float(-38.5), float(23.5), float(1.0)]

cluster_dict["10.4029998779_arrows"] += cgo_arrow([-19.0,-38.5,23.5], [-17.985,-40.3,21.035], color="red blue", name="Arrows_10.4029998779_2")

cmd.load_cgo(cluster_dict["10.4029998779"], "Features_10.4029998779", 1)
cmd.load_cgo(cluster_dict["10.4029998779_arrows"], "Arrows_10.4029998779")
cmd.set("transparency", 0.2,"Features_10.4029998779")
cmd.group("Pharmacophore_10.4029998779", members="Features_10.4029998779")
cmd.group("Pharmacophore_10.4029998779", members="Arrows_10.4029998779")

if dirpath:
    f = join(dirpath, "102/label_threshold_10.4029998779.mol2")
else:
    f = "102/label_threshold_10.4029998779.mol2"

cmd.load(f, 'label_threshold_10.4029998779')
cmd.hide('everything', 'label_threshold_10.4029998779')
cmd.label("label_threshold_10.4029998779", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_10.4029998779', members= 'label_threshold_10.4029998779')


if dirpath:
    f = join(dirpath, '102/mesh.grd')
else:
    f = '102/mesh.grd'
cmd.load(f, 'mesh_102')
cmd.isomesh("isomesh_102", "mesh_102", 0.9)
cmd.color("grey80", "isomesh_102")
cmd.set('transparency', 0.4, "isomesh_102")

cmd.group('hotspot_102', "isomesh_102")
cmd.group('hotspot_102', "mesh_102")

if dirpath:
    f = join(dirpath, "103/label_threshold_7.7.mol2")
else:
    f = "103/label_threshold_7.7.mol2"

cmd.load(f, 'label_threshold_7.7')
cmd.hide('everything', 'label_threshold_7.7')
cmd.label("label_threshold_7.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [7.7]
gfiles = ['103/donor.grd', '103/apolar.grd', '103/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 103
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


cluster_dict = {"9.93200016022":[], "9.93200016022_arrows":[]}

cluster_dict["9.93200016022"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-69.8305556866), float(-6.78822804206), float(60.2008412157), float(1.0)]


cluster_dict["9.93200016022"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-65.0), float(-6.0322755023), float(61.7949729637), float(1.0)]


cluster_dict["9.93200016022"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-60.3516282388), float(-7.00834808042), float(56.2576151527), float(1.0)]


cmd.load_cgo(cluster_dict["9.93200016022"], "Features_9.93200016022", 1)
cmd.load_cgo(cluster_dict["9.93200016022_arrows"], "Arrows_9.93200016022")
cmd.set("transparency", 0.2,"Features_9.93200016022")
cmd.group("Pharmacophore_9.93200016022", members="Features_9.93200016022")
cmd.group("Pharmacophore_9.93200016022", members="Arrows_9.93200016022")

if dirpath:
    f = join(dirpath, "103/label_threshold_9.93200016022.mol2")
else:
    f = "103/label_threshold_9.93200016022.mol2"

cmd.load(f, 'label_threshold_9.93200016022')
cmd.hide('everything', 'label_threshold_9.93200016022')
cmd.label("label_threshold_9.93200016022", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.93200016022', members= 'label_threshold_9.93200016022')


if dirpath:
    f = join(dirpath, '103/mesh.grd')
else:
    f = '103/mesh.grd'
cmd.load(f, 'mesh_103')
cmd.isomesh("isomesh_103", "mesh_103", 0.9)
cmd.color("grey80", "isomesh_103")
cmd.set('transparency', 0.4, "isomesh_103")

cmd.group('hotspot_103', "isomesh_103")
cmd.group('hotspot_103', "mesh_103")

if dirpath:
    f = join(dirpath, "104/label_threshold_2.7.mol2")
else:
    f = "104/label_threshold_2.7.mol2"

cmd.load(f, 'label_threshold_2.7')
cmd.hide('everything', 'label_threshold_2.7')
cmd.label("label_threshold_2.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [2.7]
gfiles = ['104/donor.grd', '104/apolar.grd', '104/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 104
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


cluster_dict = {"9.77200031281":[], "9.77200031281_arrows":[]}

cluster_dict["9.77200031281"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-18.2554216795), float(-13.9861642008), float(56.3348503788), float(1.0)]


cmd.load_cgo(cluster_dict["9.77200031281"], "Features_9.77200031281", 1)
cmd.load_cgo(cluster_dict["9.77200031281_arrows"], "Arrows_9.77200031281")
cmd.set("transparency", 0.2,"Features_9.77200031281")
cmd.group("Pharmacophore_9.77200031281", members="Features_9.77200031281")
cmd.group("Pharmacophore_9.77200031281", members="Arrows_9.77200031281")

if dirpath:
    f = join(dirpath, "104/label_threshold_9.77200031281.mol2")
else:
    f = "104/label_threshold_9.77200031281.mol2"

cmd.load(f, 'label_threshold_9.77200031281')
cmd.hide('everything', 'label_threshold_9.77200031281')
cmd.label("label_threshold_9.77200031281", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.77200031281', members= 'label_threshold_9.77200031281')


if dirpath:
    f = join(dirpath, '104/mesh.grd')
else:
    f = '104/mesh.grd'
cmd.load(f, 'mesh_104')
cmd.isomesh("isomesh_104", "mesh_104", 0.9)
cmd.color("grey80", "isomesh_104")
cmd.set('transparency', 0.4, "isomesh_104")

cmd.group('hotspot_104', "isomesh_104")
cmd.group('hotspot_104', "mesh_104")

if dirpath:
    f = join(dirpath, "105/label_threshold_0.5.mol2")
else:
    f = "105/label_threshold_0.5.mol2"

cmd.load(f, 'label_threshold_0.5')
cmd.hide('everything', 'label_threshold_0.5')
cmd.label("label_threshold_0.5", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [0.5]
gfiles = ['105/donor.grd', '105/apolar.grd', '105/acceptor.grd']
grids = ['donor', 'apolar', 'acceptor']
num = 105
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


cluster_dict = {"9.74199962616":[], "9.74199962616_arrows":[]}

cluster_dict["9.74199962616"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(8.98726170156), float(-2.76113114428), float(-1.6512616655), float(1.0)]


cmd.load_cgo(cluster_dict["9.74199962616"], "Features_9.74199962616", 1)
cmd.load_cgo(cluster_dict["9.74199962616_arrows"], "Arrows_9.74199962616")
cmd.set("transparency", 0.2,"Features_9.74199962616")
cmd.group("Pharmacophore_9.74199962616", members="Features_9.74199962616")
cmd.group("Pharmacophore_9.74199962616", members="Arrows_9.74199962616")

if dirpath:
    f = join(dirpath, "105/label_threshold_9.74199962616.mol2")
else:
    f = "105/label_threshold_9.74199962616.mol2"

cmd.load(f, 'label_threshold_9.74199962616')
cmd.hide('everything', 'label_threshold_9.74199962616')
cmd.label("label_threshold_9.74199962616", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.74199962616', members= 'label_threshold_9.74199962616')


if dirpath:
    f = join(dirpath, '105/mesh.grd')
else:
    f = '105/mesh.grd'
cmd.load(f, 'mesh_105')
cmd.isomesh("isomesh_105", "mesh_105", 0.9)
cmd.color("grey80", "isomesh_105")
cmd.set('transparency', 0.4, "isomesh_105")

cmd.group('hotspot_105', "isomesh_105")
cmd.group('hotspot_105', "mesh_105")

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
