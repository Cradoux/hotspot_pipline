
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


cluster_dict = {"19.0729999542":[], "19.0729999542_arrows":[]}

cluster_dict["19.0729999542"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(27.5), float(10.0), float(-34.0), float(1.0)]

cluster_dict["19.0729999542_arrows"] += cgo_arrow([27.5,10.0,-34.0], [26.257,12.933,-32.801], color="blue red", name="Arrows_19.0729999542_1")

cluster_dict["19.0729999542"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(10.5), float(-39.0), float(1.0)]

cluster_dict["19.0729999542_arrows"] += cgo_arrow([28.0,10.5,-39.0], [27.961,13.592,-39.166], color="blue red", name="Arrows_19.0729999542_2")

cluster_dict["19.0729999542"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(31.0), float(10.0), float(-30.5), float(1.0)]

cluster_dict["19.0729999542_arrows"] += cgo_arrow([31.0,10.0,-30.5], [27.36,9.822,-28.982], color="blue red", name="Arrows_19.0729999542_3")

cluster_dict["19.0729999542"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(32.0), float(3.5), float(-37.0), float(1.0)]

cluster_dict["19.0729999542_arrows"] += cgo_arrow([32.0,3.5,-37.0], [32.104,1.695,-34.881], color="blue red", name="Arrows_19.0729999542_4")

cluster_dict["19.0729999542"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(33.5), float(3.5), float(-42.0), float(1.0)]

cluster_dict["19.0729999542_arrows"] += cgo_arrow([33.5,3.5,-42.0], [30.388,4.44,-42.566], color="blue red", name="Arrows_19.0729999542_5")

cluster_dict["19.0729999542"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(36.5), float(2.5), float(-43.5), float(1.0)]

cluster_dict["19.0729999542_arrows"] += cgo_arrow([36.5,2.5,-43.5], [38.722,3.854,-44.392], color="blue red", name="Arrows_19.0729999542_6")

cluster_dict["19.0729999542"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(32.4641646881), float(6.28149811423), float(-38.7841485447), float(1.0)]


cluster_dict["19.0729999542"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(30.5), float(7.5), float(-34.5), float(1.0)]

cluster_dict["19.0729999542_arrows"] += cgo_arrow([30.5,7.5,-34.5], [28.503,5.78,-35.679], color="red blue", name="Arrows_19.0729999542_7")

cluster_dict["19.0729999542"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(33.5), float(3.0), float(-44.0), float(1.0)]

cluster_dict["19.0729999542_arrows"] += cgo_arrow([33.5,3.0,-44.0], [30.965,4.679,-44.712], color="red blue", name="Arrows_19.0729999542_8")

cluster_dict["19.0729999542"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(33.0), float(13.0), float(-28.5), float(1.0)]

cluster_dict["19.0729999542_arrows"] += cgo_arrow([33.0,13.0,-28.5], [32.581,16.518,-27.932], color="red blue", name="Arrows_19.0729999542_9")

cmd.load_cgo(cluster_dict["19.0729999542"], "Features_19.0729999542", 1)
cmd.load_cgo(cluster_dict["19.0729999542_arrows"], "Arrows_19.0729999542")
cmd.set("transparency", 0.2,"Features_19.0729999542")
cmd.group("Pharmacophore_19.0729999542", members="Features_19.0729999542")
cmd.group("Pharmacophore_19.0729999542", members="Arrows_19.0729999542")

if dirpath:
    f = join(dirpath, "label_threshold_19.0729999542.mol2")
else:
    f = "label_threshold_19.0729999542.mol2"

cmd.load(f, 'label_threshold_19.0729999542')
cmd.hide('everything', 'label_threshold_19.0729999542')
cmd.label("label_threshold_19.0729999542", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_19.0729999542', members= 'label_threshold_19.0729999542')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
