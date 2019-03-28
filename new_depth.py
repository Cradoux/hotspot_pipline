from ccdc.protein import Protein
from hotspots.grid_extension import Grid
import numpy as np
import scipy.spatial as spatial


def multi_probes(mol, scaling=1):
    probe_sizes = [10,5,4,3]

    coords = [a.coordinates for a in mol.atoms]
    g = Grid.initalise_grid(coords=coords, padding=2)
    for probe in probe_sizes:
        for a in mol.heavy_atoms:
            g.set_sphere(point=a.coordinates,
                         radius=probe * scaling,
                         value=probe,
                         scaling='None')
    return g

def split_array(g, points_tree):
    probe_sizes = [3, 4,5,10]

    nx, ny, nz = g.nsteps

    solv_coords = []
    g2 = g.copy_and_clear()

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if g.value(i, j, k) == 0:
                    solv_coords.append((i, j, k))

    for probe in probe_sizes:
        print probe
        remaining = solv_coords

        for c in remaining:
            i, j, k = c
            coord = g2.indices_to_point(i, j, k)
            if len(points_tree.query_ball_point(coord, probe)) ==0:
                g2.set_value(i, j, k, 10 - probe)

                #solv_coords.remove((i,j,k))

    return g2.min_value_of_neighbours()


prot = Protein.from_file('pdb_files/aa/2aa2/protonated_no_wat.pdb')
prot_coords = [a.coordinates for a in prot.heavy_atoms]
prot_grid = multi_probes(prot)

points_tree = spatial.cKDTree(prot_coords)

#out_g = split_array(prot_grid, points_tree)
prot_grid.write(('prot_g.grd'))
#out_g.write('newdepth.grd')
