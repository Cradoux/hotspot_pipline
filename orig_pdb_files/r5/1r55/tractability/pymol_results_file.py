
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
    f = join(dirpath, "0/label_threshold_3.7.mol2")
else:
    f = "0/label_threshold_3.7.mol2"

cmd.load(f, 'label_threshold_3.7')
cmd.hide('everything', 'label_threshold_3.7')
cmd.label("label_threshold_3.7", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [3.7]
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


cluster_dict = {"9.98200035095":[], "9.98200035095_arrows":[]}

cluster_dict["9.98200035095"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(9.77746449518), float(15.838065172), float(25.8169610088), float(1.0)]


cluster_dict["9.98200035095"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(10.0), float(13.5), float(26.5), float(1.0)]

cluster_dict["9.98200035095_arrows"] += cgo_arrow([10.0,13.5,26.5], [10.834,11.11,28.613], color="red blue", name="Arrows_9.98200035095_1")

cluster_dict["9.98200035095"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(13.5), float(16.5), float(26.0), float(1.0)]

cluster_dict["9.98200035095_arrows"] += cgo_arrow([13.5,16.5,26.0], [15.637,15.324,27.84], color="red blue", name="Arrows_9.98200035095_2")

cmd.load_cgo(cluster_dict["9.98200035095"], "Features_9.98200035095", 1)
cmd.load_cgo(cluster_dict["9.98200035095_arrows"], "Arrows_9.98200035095")
cmd.set("transparency", 0.2,"Features_9.98200035095")
cmd.group("Pharmacophore_9.98200035095", members="Features_9.98200035095")
cmd.group("Pharmacophore_9.98200035095", members="Arrows_9.98200035095")

if dirpath:
    f = join(dirpath, "0/label_threshold_9.98200035095.mol2")
else:
    f = "0/label_threshold_9.98200035095.mol2"

cmd.load(f, 'label_threshold_9.98200035095')
cmd.hide('everything', 'label_threshold_9.98200035095')
cmd.label("label_threshold_9.98200035095", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_9.98200035095', members= 'label_threshold_9.98200035095')


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

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")