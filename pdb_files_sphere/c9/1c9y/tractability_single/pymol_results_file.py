
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


cluster_dict = {"13.6379995346":[], "13.6379995346_arrows":[]}

cluster_dict["13.6379995346"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(0.0), float(2.0), float(-25.0), float(1.0)]

cluster_dict["13.6379995346_arrows"] += cgo_arrow([0.0,2.0,-25.0], [-1.832,0.903,-27.0], color="blue red", name="Arrows_13.6379995346_1")

cluster_dict["13.6379995346"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(6.5), float(3.5), float(-21.0), float(1.0)]

cluster_dict["13.6379995346_arrows"] += cgo_arrow([6.5,3.5,-21.0], [7.18,6.139,-21.718], color="blue red", name="Arrows_13.6379995346_2")

cluster_dict["13.6379995346"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(6.5), float(3.5), float(-21.0), float(1.0)]

cluster_dict["13.6379995346_arrows"] += cgo_arrow([6.5,3.5,-21.0], [7.18,6.139,-21.718], color="blue red", name="Arrows_13.6379995346_3")

cluster_dict["13.6379995346"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(10.5), float(3.5), float(-21.5), float(1.0)]

cluster_dict["13.6379995346_arrows"] += cgo_arrow([10.5,3.5,-21.5], [10.948,0.411,-23.095], color="blue red", name="Arrows_13.6379995346_4")

cluster_dict["13.6379995346"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(0.216889139932), float(4.21901718155), float(-19.2995737819), float(1.0)]


cluster_dict["13.6379995346"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(10.7623762376), float(4.57920792079), float(-21.9158415842), float(1.0)]


cluster_dict["13.6379995346"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-0.5), float(1.0), float(-22.5), float(1.0)]

cluster_dict["13.6379995346_arrows"] += cgo_arrow([-0.5,1.0,-22.5], [0.066,-1.533,-21.14], color="red blue", name="Arrows_13.6379995346_5")

cluster_dict["13.6379995346"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(4.5), float(2.5), float(-21.5), float(1.0)]

cluster_dict["13.6379995346_arrows"] += cgo_arrow([4.5,2.5,-21.5], [6.712,1.082,-23.059], color="red blue", name="Arrows_13.6379995346_6")

cluster_dict["13.6379995346"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(11.0), float(3.5), float(-20.5), float(1.0)]

cluster_dict["13.6379995346_arrows"] += cgo_arrow([11.0,3.5,-20.5], [13.023,5.316,-18.731], color="red blue", name="Arrows_13.6379995346_7")

cmd.load_cgo(cluster_dict["13.6379995346"], "Features_13.6379995346", 1)
cmd.load_cgo(cluster_dict["13.6379995346_arrows"], "Arrows_13.6379995346")
cmd.set("transparency", 0.2,"Features_13.6379995346")
cmd.group("Pharmacophore_13.6379995346", members="Features_13.6379995346")
cmd.group("Pharmacophore_13.6379995346", members="Arrows_13.6379995346")

if dirpath:
    f = join(dirpath, "label_threshold_13.6379995346.mol2")
else:
    f = "label_threshold_13.6379995346.mol2"

cmd.load(f, 'label_threshold_13.6379995346')
cmd.hide('everything', 'label_threshold_13.6379995346')
cmd.label("label_threshold_13.6379995346", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_13.6379995346', members= 'label_threshold_13.6379995346')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
