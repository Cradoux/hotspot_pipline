
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
    f = join(dirpath, "0/label_threshold_15.9.mol2")
else:
    f = "0/label_threshold_15.9.mol2"

cmd.load(f, 'label_threshold_15.9')
cmd.hide('everything', 'label_threshold_15.9')
cmd.label("label_threshold_15.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.9]
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


cluster_dict = {"18.063999176":[], "18.063999176_arrows":[]}

cluster_dict["18.063999176"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(16.0), float(33.0), float(1.0)]

cluster_dict["18.063999176_arrows"] += cgo_arrow([16.0,16.0,33.0], [15.577,15.193,30.195], color="blue red", name="Arrows_18.063999176_1")

cluster_dict["18.063999176"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(19.0), float(26.5), float(1.0)]

cluster_dict["18.063999176_arrows"] += cgo_arrow([22.0,19.0,26.5], [20.586,17.072,25.369], color="blue red", name="Arrows_18.063999176_2")

cluster_dict["18.063999176"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(11.5), float(29.0), float(1.0)]

cluster_dict["18.063999176_arrows"] += cgo_arrow([22.5,11.5,29.0], [22.706,12.303,26.002], color="blue red", name="Arrows_18.063999176_3")

cluster_dict["18.063999176"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(4.0), float(32.5), float(1.0)]

cluster_dict["18.063999176_arrows"] += cgo_arrow([23.5,4.0,32.5], [21.545,2.876,32.024], color="blue red", name="Arrows_18.063999176_4")

cluster_dict["18.063999176"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(4.5), float(29.0), float(1.0)]

cluster_dict["18.063999176_arrows"] += cgo_arrow([26.0,4.5,29.0], [26.451,4.735,26.144], color="blue red", name="Arrows_18.063999176_5")

cluster_dict["18.063999176"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(21.8889839175), float(10.4796883837), float(32.4819449606), float(1.0)]


cluster_dict["18.063999176"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(16.0), float(33.5), float(1.0)]

cluster_dict["18.063999176_arrows"] += cgo_arrow([18.5,16.0,33.5], [20.616,17.371,33.653], color="red blue", name="Arrows_18.063999176_6")

cluster_dict["18.063999176"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.5), float(6.5), float(33.5), float(1.0)]

cluster_dict["18.063999176_arrows"] += cgo_arrow([20.5,6.5,33.5], [19.211,3.764,32.937], color="red blue", name="Arrows_18.063999176_7")

cluster_dict["18.063999176"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(11.0), float(32.0), float(1.0)]


cluster_dict["18.063999176"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(25.5), float(10.5), float(27.5), float(1.0)]

cluster_dict["18.063999176_arrows"] += cgo_arrow([25.5,10.5,27.5], [25.881,9.008,25.363], color="red blue", name="Arrows_18.063999176_8")

cmd.load_cgo(cluster_dict["18.063999176"], "Features_18.063999176", 1)
cmd.load_cgo(cluster_dict["18.063999176_arrows"], "Arrows_18.063999176")
cmd.set("transparency", 0.2,"Features_18.063999176")
cmd.group("Pharmacophore_18.063999176", members="Features_18.063999176")
cmd.group("Pharmacophore_18.063999176", members="Arrows_18.063999176")

if dirpath:
    f = join(dirpath, "0/label_threshold_18.063999176.mol2")
else:
    f = "0/label_threshold_18.063999176.mol2"

cmd.load(f, 'label_threshold_18.063999176')
cmd.hide('everything', 'label_threshold_18.063999176')
cmd.label("label_threshold_18.063999176", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_18.063999176', members= 'label_threshold_18.063999176')


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
    f = join(dirpath, "1/label_threshold_15.8.mol2")
else:
    f = "1/label_threshold_15.8.mol2"

cmd.load(f, 'label_threshold_15.8')
cmd.hide('everything', 'label_threshold_15.8')
cmd.label("label_threshold_15.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.8]
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


cluster_dict = {"18.0429992676":[], "18.0429992676_arrows":[]}

cluster_dict["18.0429992676"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(16.5), float(33.0), float(1.0)]

cluster_dict["18.0429992676_arrows"] += cgo_arrow([16.0,16.5,33.0], [15.577,15.193,30.195], color="blue red", name="Arrows_18.0429992676_1")

cluster_dict["18.0429992676"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(11.5), float(29.0), float(1.0)]

cluster_dict["18.0429992676_arrows"] += cgo_arrow([22.5,11.5,29.0], [22.706,12.303,26.002], color="blue red", name="Arrows_18.0429992676_2")

cluster_dict["18.0429992676"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(4.0), float(32.5), float(1.0)]

cluster_dict["18.0429992676_arrows"] += cgo_arrow([23.5,4.0,32.5], [21.545,2.876,32.024], color="blue red", name="Arrows_18.0429992676_3")

cluster_dict["18.0429992676"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(25.5), float(4.5), float(29.0), float(1.0)]

cluster_dict["18.0429992676_arrows"] += cgo_arrow([25.5,4.5,29.0], [26.451,4.735,26.144], color="blue red", name="Arrows_18.0429992676_4")

cluster_dict["18.0429992676"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(21.9154976093), float(10.511688745), float(32.4605916389), float(1.0)]


cluster_dict["18.0429992676"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(16.0), float(33.5), float(1.0)]

cluster_dict["18.0429992676_arrows"] += cgo_arrow([18.5,16.0,33.5], [20.616,17.371,33.653], color="red blue", name="Arrows_18.0429992676_5")

cluster_dict["18.0429992676"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.5), float(6.5), float(33.5), float(1.0)]

cluster_dict["18.0429992676_arrows"] += cgo_arrow([20.5,6.5,33.5], [19.211,3.764,32.937], color="red blue", name="Arrows_18.0429992676_6")

cluster_dict["18.0429992676"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(11.0), float(32.0), float(1.0)]


cluster_dict["18.0429992676"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(25.5), float(10.5), float(27.5), float(1.0)]

cluster_dict["18.0429992676_arrows"] += cgo_arrow([25.5,10.5,27.5], [25.881,9.008,25.363], color="red blue", name="Arrows_18.0429992676_7")

cmd.load_cgo(cluster_dict["18.0429992676"], "Features_18.0429992676", 1)
cmd.load_cgo(cluster_dict["18.0429992676_arrows"], "Arrows_18.0429992676")
cmd.set("transparency", 0.2,"Features_18.0429992676")
cmd.group("Pharmacophore_18.0429992676", members="Features_18.0429992676")
cmd.group("Pharmacophore_18.0429992676", members="Arrows_18.0429992676")

if dirpath:
    f = join(dirpath, "1/label_threshold_18.0429992676.mol2")
else:
    f = "1/label_threshold_18.0429992676.mol2"

cmd.load(f, 'label_threshold_18.0429992676')
cmd.hide('everything', 'label_threshold_18.0429992676')
cmd.label("label_threshold_18.0429992676", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_18.0429992676', members= 'label_threshold_18.0429992676')


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
    f = join(dirpath, "2/label_threshold_15.7.mol2")
else:
    f = "2/label_threshold_15.7.mol2"

cmd.load(f, 'label_threshold_15.7')
cmd.hide('everything', 'label_threshold_15.7')
cmd.label("label_threshold_15.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.7]
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


cluster_dict = {"17.3010005951":[], "17.3010005951_arrows":[]}

cluster_dict["17.3010005951"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(16.0), float(16.5), float(33.0), float(1.0)]

cluster_dict["17.3010005951_arrows"] += cgo_arrow([16.0,16.5,33.0], [15.577,15.193,30.195], color="blue red", name="Arrows_17.3010005951_1")

cluster_dict["17.3010005951"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(19.5), float(26.5), float(1.0)]

cluster_dict["17.3010005951_arrows"] += cgo_arrow([22.0,19.5,26.5], [20.586,17.072,25.369], color="blue red", name="Arrows_17.3010005951_2")

cluster_dict["17.3010005951"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(23.5), float(27.5), float(1.0)]

cluster_dict["17.3010005951_arrows"] += cgo_arrow([24.0,23.5,27.5], [23.825,23.679,30.554], color="blue red", name="Arrows_17.3010005951_3")

cluster_dict["17.3010005951"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(25.5), float(29.5), float(25.0), float(1.0)]

cluster_dict["17.3010005951_arrows"] += cgo_arrow([25.5,29.5,25.0], [24.804,32.148,26.296], color="blue red", name="Arrows_17.3010005951_4")

cluster_dict["17.3010005951"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(8.69957956514), float(26.9739733793), float(31.3327548703), float(1.0)]


cluster_dict["17.3010005951"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(21.5137090679), float(26.4217042324), float(26.8704549212), float(1.0)]


cluster_dict["17.3010005951"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(16.9833959195), float(15.617151671), float(33.7582549799), float(1.0)]


cluster_dict["17.3010005951"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.7584864267), float(14.9632433317), float(31.3236766459), float(1.0)]


cluster_dict["17.3010005951"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(26.5), float(29.0), float(1.0)]

cluster_dict["17.3010005951_arrows"] += cgo_arrow([15.0,26.5,29.0], [14.267,27.323,25.392], color="red blue", name="Arrows_17.3010005951_5")

cluster_dict["17.3010005951"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(16.0), float(33.5), float(1.0)]

cluster_dict["17.3010005951_arrows"] += cgo_arrow([18.5,16.0,33.5], [20.616,17.371,33.653], color="red blue", name="Arrows_17.3010005951_6")

cluster_dict["17.3010005951"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(19.0), float(14.5), float(32.5), float(1.0)]

cluster_dict["17.3010005951_arrows"] += cgo_arrow([19.0,14.5,32.5], [20.616,17.371,33.653], color="red blue", name="Arrows_17.3010005951_7")

cluster_dict["17.3010005951"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.5), float(31.0), float(28.5), float(1.0)]

cluster_dict["17.3010005951_arrows"] += cgo_arrow([20.5,31.0,28.5], [22.1,33.173,28.772], color="red blue", name="Arrows_17.3010005951_8")

cluster_dict["17.3010005951"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(26.5), float(23.0), float(1.0)]

cluster_dict["17.3010005951_arrows"] += cgo_arrow([22.5,26.5,23.0], [21.419,24.54,21.082], color="red blue", name="Arrows_17.3010005951_9")

cmd.load_cgo(cluster_dict["17.3010005951"], "Features_17.3010005951", 1)
cmd.load_cgo(cluster_dict["17.3010005951_arrows"], "Arrows_17.3010005951")
cmd.set("transparency", 0.2,"Features_17.3010005951")
cmd.group("Pharmacophore_17.3010005951", members="Features_17.3010005951")
cmd.group("Pharmacophore_17.3010005951", members="Arrows_17.3010005951")

if dirpath:
    f = join(dirpath, "2/label_threshold_17.3010005951.mol2")
else:
    f = "2/label_threshold_17.3010005951.mol2"

cmd.load(f, 'label_threshold_17.3010005951')
cmd.hide('everything', 'label_threshold_17.3010005951')
cmd.label("label_threshold_17.3010005951", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.3010005951', members= 'label_threshold_17.3010005951')


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
    f = join(dirpath, "3/label_threshold_14.8.mol2")
else:
    f = "3/label_threshold_14.8.mol2"

cmd.load(f, 'label_threshold_14.8')
cmd.hide('everything', 'label_threshold_14.8')
cmd.label("label_threshold_14.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [14.8]
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


cluster_dict = {"17.0909996033":[], "17.0909996033_arrows":[]}

cluster_dict["17.0909996033"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(15.5), float(24.5), float(28.5), float(1.0)]

cluster_dict["17.0909996033_arrows"] += cgo_arrow([15.5,24.5,28.5], [15.481,23.448,26.221], color="blue red", name="Arrows_17.0909996033_1")

cluster_dict["17.0909996033"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(17.0), float(24.5), float(28.0), float(1.0)]

cluster_dict["17.0909996033_arrows"] += cgo_arrow([17.0,24.5,28.0], [15.481,23.448,26.221], color="blue red", name="Arrows_17.0909996033_2")

cluster_dict["17.0909996033"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(19.5), float(26.5), float(1.0)]

cluster_dict["17.0909996033_arrows"] += cgo_arrow([22.0,19.5,26.5], [20.586,17.072,25.369], color="blue red", name="Arrows_17.0909996033_3")

cluster_dict["17.0909996033"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(23.5), float(27.5), float(1.0)]

cluster_dict["17.0909996033_arrows"] += cgo_arrow([24.0,23.5,27.5], [23.825,23.679,30.554], color="blue red", name="Arrows_17.0909996033_4")

cluster_dict["17.0909996033"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(25.5), float(29.5), float(25.0), float(1.0)]

cluster_dict["17.0909996033_arrows"] += cgo_arrow([25.5,29.5,25.0], [24.804,32.148,26.296], color="blue red", name="Arrows_17.0909996033_5")

cluster_dict["17.0909996033"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(21.7150717189), float(26.4788001982), float(26.7294506445), float(1.0)]


cluster_dict["17.0909996033"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(26.5), float(29.0), float(1.0)]

cluster_dict["17.0909996033_arrows"] += cgo_arrow([15.0,26.5,29.0], [14.267,27.323,25.392], color="red blue", name="Arrows_17.0909996033_6")

cluster_dict["17.0909996033"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.0), float(21.0), float(27.0), float(1.0)]

cluster_dict["17.0909996033_arrows"] += cgo_arrow([20.0,21.0,27.0], [19.863,20.032,28.99], color="red blue", name="Arrows_17.0909996033_7")

cluster_dict["17.0909996033"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.5), float(31.0), float(28.5), float(1.0)]

cluster_dict["17.0909996033_arrows"] += cgo_arrow([20.5,31.0,28.5], [22.1,33.173,28.772], color="red blue", name="Arrows_17.0909996033_8")

cluster_dict["17.0909996033"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(26.5), float(23.0), float(1.0)]

cluster_dict["17.0909996033_arrows"] += cgo_arrow([22.5,26.5,23.0], [21.419,24.54,21.082], color="red blue", name="Arrows_17.0909996033_9")

cmd.load_cgo(cluster_dict["17.0909996033"], "Features_17.0909996033", 1)
cmd.load_cgo(cluster_dict["17.0909996033_arrows"], "Arrows_17.0909996033")
cmd.set("transparency", 0.2,"Features_17.0909996033")
cmd.group("Pharmacophore_17.0909996033", members="Features_17.0909996033")
cmd.group("Pharmacophore_17.0909996033", members="Arrows_17.0909996033")

if dirpath:
    f = join(dirpath, "3/label_threshold_17.0909996033.mol2")
else:
    f = "3/label_threshold_17.0909996033.mol2"

cmd.load(f, 'label_threshold_17.0909996033')
cmd.hide('everything', 'label_threshold_17.0909996033')
cmd.label("label_threshold_17.0909996033", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.0909996033', members= 'label_threshold_17.0909996033')


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
    f = join(dirpath, "4/label_threshold_12.9.mol2")
else:
    f = "4/label_threshold_12.9.mol2"

cmd.load(f, 'label_threshold_12.9')
cmd.hide('everything', 'label_threshold_12.9')
cmd.label("label_threshold_12.9", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.9]
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


cluster_dict = {"15.7220001221":[], "15.7220001221_arrows":[]}

cluster_dict["15.7220001221"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(17.0), float(24.5), float(28.0), float(1.0)]

cluster_dict["15.7220001221_arrows"] += cgo_arrow([17.0,24.5,28.0], [15.481,23.448,26.221], color="blue red", name="Arrows_15.7220001221_1")

cluster_dict["15.7220001221"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.0), float(23.0), float(25.0), float(1.0)]

cluster_dict["15.7220001221_arrows"] += cgo_arrow([19.0,23.0,25.0], [20.129,22.061,22.222], color="blue red", name="Arrows_15.7220001221_2")

cluster_dict["15.7220001221"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(7.65729572299), float(26.5177653269), float(30.740019596), float(1.0)]


cluster_dict["15.7220001221"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(18.6485577843), float(27.3641257771), float(27.6354274012), float(1.0)]


cluster_dict["15.7220001221"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(10.0), float(28.0), float(28.5), float(1.0)]

cluster_dict["15.7220001221_arrows"] += cgo_arrow([10.0,28.0,28.5], [11.2,25.378,28.335], color="red blue", name="Arrows_15.7220001221_3")

cluster_dict["15.7220001221"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(10.5), float(29.0), float(30.0), float(1.0)]

cluster_dict["15.7220001221_arrows"] += cgo_arrow([10.5,29.0,30.0], [13.101,30.805,29.788], color="red blue", name="Arrows_15.7220001221_4")

cluster_dict["15.7220001221"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(15.0), float(26.5), float(29.0), float(1.0)]

cluster_dict["15.7220001221_arrows"] += cgo_arrow([15.0,26.5,29.0], [14.267,27.323,25.392], color="red blue", name="Arrows_15.7220001221_5")

cluster_dict["15.7220001221"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(17.5), float(28.0), float(25.5), float(1.0)]

cluster_dict["15.7220001221_arrows"] += cgo_arrow([17.5,28.0,25.5], [14.267,27.323,25.392], color="red blue", name="Arrows_15.7220001221_6")

cluster_dict["15.7220001221"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(19.5), float(21.0), float(27.0), float(1.0)]

cluster_dict["15.7220001221_arrows"] += cgo_arrow([19.5,21.0,27.0], [19.863,20.032,28.99], color="red blue", name="Arrows_15.7220001221_7")

cluster_dict["15.7220001221"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.0), float(31.5), float(29.0), float(1.0)]

cluster_dict["15.7220001221_arrows"] += cgo_arrow([20.0,31.5,29.0], [22.1,33.173,28.772], color="red blue", name="Arrows_15.7220001221_8")

cluster_dict["15.7220001221"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(20.0), float(31.5), float(29.0), float(1.0)]

cluster_dict["15.7220001221_arrows"] += cgo_arrow([20.0,31.5,29.0], [22.1,33.173,28.772], color="red blue", name="Arrows_15.7220001221_9")

cmd.load_cgo(cluster_dict["15.7220001221"], "Features_15.7220001221", 1)
cmd.load_cgo(cluster_dict["15.7220001221_arrows"], "Arrows_15.7220001221")
cmd.set("transparency", 0.2,"Features_15.7220001221")
cmd.group("Pharmacophore_15.7220001221", members="Features_15.7220001221")
cmd.group("Pharmacophore_15.7220001221", members="Arrows_15.7220001221")

if dirpath:
    f = join(dirpath, "4/label_threshold_15.7220001221.mol2")
else:
    f = "4/label_threshold_15.7220001221.mol2"

cmd.load(f, 'label_threshold_15.7220001221')
cmd.hide('everything', 'label_threshold_15.7220001221')
cmd.label("label_threshold_15.7220001221", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.7220001221', members= 'label_threshold_15.7220001221')


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
    f = join(dirpath, "5/label_threshold_12.2.mol2")
else:
    f = "5/label_threshold_12.2.mol2"

cmd.load(f, 'label_threshold_12.2')
cmd.hide('everything', 'label_threshold_12.2')
cmd.label("label_threshold_12.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [12.2]
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


cluster_dict = {"15.6330003738":[], "15.6330003738_arrows":[]}

cluster_dict["15.6330003738"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(19.5), float(27.0), float(1.0)]

cluster_dict["15.6330003738_arrows"] += cgo_arrow([22.5,19.5,27.0], [19.863,20.032,28.99], color="blue red", name="Arrows_15.6330003738_1")

cluster_dict["15.6330003738"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(15.5), float(30.5), float(1.0)]

cluster_dict["15.6330003738_arrows"] += cgo_arrow([22.5,15.5,30.5], [20.616,17.371,33.653], color="blue red", name="Arrows_15.6330003738_2")

cluster_dict["15.6330003738"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(23.5), float(27.5), float(1.0)]

cluster_dict["15.6330003738_arrows"] += cgo_arrow([24.0,23.5,27.5], [23.825,23.679,30.554], color="blue red", name="Arrows_15.6330003738_3")

cluster_dict["15.6330003738"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(12.0), float(29.5), float(1.0)]

cluster_dict["15.6330003738_arrows"] += cgo_arrow([23.5,12.0,29.5], [22.706,12.303,26.002], color="blue red", name="Arrows_15.6330003738_4")

cluster_dict["15.6330003738"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(25.0), float(28.0), float(1.0)]

cluster_dict["15.6330003738_arrows"] += cgo_arrow([28.0,25.0,28.0], [30.569,24.211,26.718], color="blue red", name="Arrows_15.6330003738_5")

cluster_dict["15.6330003738"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(23.9119642669), float(16.0139146476), float(30.9798658032), float(1.0)]


cluster_dict["15.6330003738"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(18.5), float(16.0), float(33.5), float(1.0)]

cluster_dict["15.6330003738_arrows"] += cgo_arrow([18.5,16.0,33.5], [20.616,17.371,33.653], color="red blue", name="Arrows_15.6330003738_6")

cluster_dict["15.6330003738"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.5), float(18.0), float(28.5), float(1.0)]

cluster_dict["15.6330003738_arrows"] += cgo_arrow([22.5,18.0,28.5], [21.002,15.7,27.106], color="red blue", name="Arrows_15.6330003738_7")

cluster_dict["15.6330003738"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.5), float(11.0), float(32.0), float(1.0)]


cluster_dict["15.6330003738"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(26.0), float(10.5), float(32.5), float(1.0)]

cluster_dict["15.6330003738_arrows"] += cgo_arrow([26.0,10.5,32.5], [28.772,9.721,33.856], color="red blue", name="Arrows_15.6330003738_8")

cmd.load_cgo(cluster_dict["15.6330003738"], "Features_15.6330003738", 1)
cmd.load_cgo(cluster_dict["15.6330003738_arrows"], "Arrows_15.6330003738")
cmd.set("transparency", 0.2,"Features_15.6330003738")
cmd.group("Pharmacophore_15.6330003738", members="Features_15.6330003738")
cmd.group("Pharmacophore_15.6330003738", members="Arrows_15.6330003738")

if dirpath:
    f = join(dirpath, "5/label_threshold_15.6330003738.mol2")
else:
    f = "5/label_threshold_15.6330003738.mol2"

cmd.load(f, 'label_threshold_15.6330003738')
cmd.hide('everything', 'label_threshold_15.6330003738')
cmd.label("label_threshold_15.6330003738", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_15.6330003738', members= 'label_threshold_15.6330003738')


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


cluster_dict = {"9.11299991608":[], "9.11299991608_arrows":[]}

cluster_dict["9.11299991608"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(19.5), float(-1.5), float(8.5), float(1.0)]

cluster_dict["9.11299991608_arrows"] += cgo_arrow([19.5,-1.5,8.5], [20.287,-2.168,11.435], color="blue red", name="Arrows_9.11299991608_1")

cluster_dict["9.11299991608"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(24.0), float(0.5), float(9.0), float(1.0)]

cluster_dict["9.11299991608_arrows"] += cgo_arrow([24.0,0.5,9.0], [23.276,0.393,5.779], color="blue red", name="Arrows_9.11299991608_2")

cluster_dict["9.11299991608"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(22.0364966367), float(-0.439892640061), float(8.51595527091), float(1.0)]


cluster_dict["9.11299991608"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(25.6572399031), float(2.19792574623), float(8.79137605577), float(1.0)]


cmd.load_cgo(cluster_dict["9.11299991608"], "Features_9.11299991608", 1)
cmd.load_cgo(cluster_dict["9.11299991608_arrows"], "Arrows_9.11299991608")
cmd.set("transparency", 0.2,"Features_9.11299991608")
cmd.group("Pharmacophore_9.11299991608", members="Features_9.11299991608")
cmd.group("Pharmacophore_9.11299991608", members="Arrows_9.11299991608")

if dirpath:
    f = join(dirpath, "6/label_threshold_9.11299991608.mol2")
else:
    f = "6/label_threshold_9.11299991608.mol2"

cmd.load(f, 'label_threshold_9.11299991608')
cmd.hide('everything', 'label_threshold_9.11299991608')
cmd.label("label_threshold_9.11299991608", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.11299991608', members= 'label_threshold_9.11299991608')


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


cluster_dict = {"0":[], "0_arrows":[]}

cmd.load_cgo(cluster_dict["0"], "Features_0", 1)
cmd.load_cgo(cluster_dict["0_arrows"], "Arrows_0")
cmd.set("transparency", 0.2,"Features_0")
cmd.group("Pharmacophore_0", members="Features_0")
cmd.group("Pharmacophore_0", members="Arrows_0")

if dirpath:
    f = join(dirpath, "7/label_threshold_0.mol2")
else:
    f = "7/label_threshold_0.mol2"

cmd.load(f, 'label_threshold_0')
cmd.hide('everything', 'label_threshold_0')
cmd.label("label_threshold_0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_0', members= 'label_threshold_0')


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

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
