import luigi
from urllib2 import urlopen as urlopen
from ccdc.protein import Protein
from ccdc.cavity import Cavity
from ccdc import io
import subprocess
from hotspots.calculation import Runner
from hotspots.hs_io import HotspotWriter, HotspotReader
import sys
import os
import csv
import pandas as pd





class GetPDB(luigi.Task):
    pdb = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget("pdb_files/{}/{}/biological_assembly.pdb".format(self.pdb[1:3], self.pdb))

    def requires(self):
        pass

    def run(self):
        url = "https://files.rcsb.org/download/{}.pdb1".format(self.pdb)
        print(url)
        pdbfile = urlopen(url).read()

        with self.output().open('w') as out_file:
            out_file.write(pdbfile)


class Protonate(luigi.Task):
    pdb = luigi.Parameter()

    def requires(self):
        return GetPDB(self.pdb)

    def output(self):
        return luigi.LocalTarget("pdb_files/{}/{}/protonated.pdb".format(self.pdb[1:3], self.pdb))

    def run(self):
        cmd = "/home/cradoux/Research/pdb2pqr/pdb2pqr-linux-bin64-2.1.0/pdb2pqr --ff=amber --chain {} {}".format(self.input().path,
                                                                                           self.output().path)

        print cmd

        subprocess.call(cmd, shell=sys.platform != 'win32')


class CSDProteinPrep(luigi.Task):
    pdb = luigi.Parameter()

    def requires(self):
        return Protonate(self.pdb)

    def output(self):
        return luigi.LocalTarget("pdb_files/{}/{}/protonated_no_wat.pdb".format(self.pdb[1:3], self.pdb))

    def run(self):

        prot = Protein.from_file(self.input().path)
        prot.remove_all_waters()
        prot.detect_ligand_bonds()
        for l in prot.ligands:
            if 'HEM' not in l.identifier:
                prot.remove_ligand(l.identifier)

        prot.add_hydrogens()
        with io.MoleculeWriter(self.output().path) as w:
            w.write(prot)


class CSDLigandPrep(luigi.Task):
    pdb = luigi.Parameter()

    def requires(self):
        return GetPDB(self.pdb)

    def output(self):
        return luigi.LocalTarget("pdb_files/{}/{}/ligands.sdf".format(self.pdb[1:3], self.pdb))

    def run(self):
        prot = Protein.from_file(self.input().path)
        prot.detect_ligand_bonds()
        prot.add_hydrogens()
        with io.MoleculeWriter(self.output().path) as w:
            for l in prot.ligands:
                if 'HEM' not in l.identifier:
                    w.write(l)


class GetHotspots(luigi.Task):
    pdb = luigi.Parameter()

    def requires(self):
        return CSDProteinPrep(self.pdb)

    def output(self):
        return luigi.LocalTarget("pdb_files/{}/{}/hotspot/out.zip".format(self.pdb[1:3], self.pdb))

    def run(self):
        prot = Protein.from_file(self.input().path)
        #mol = io.MoleculeReader('ligands/{}.sdf'.format(self.pdb))[0]

        h = Runner()
        s = h.Settings()
        s.apolar_translation_threshold = 15
        s.polar_translation_threshold = 15
        s.polar_contributions = False
        s.nrotations = 3000
        hr = h.from_protein(prot, buriedness_method='ghecom', nprocesses=1, settings=s)

        out_settings = HotspotWriter.Settings()
        out_settings.charged = False
        w = HotspotWriter(os.path.dirname(self.output().path), grid_extension=".grd", zip_results=True,
                          settings=out_settings)

        w.write(hr)


class ScoreLigands(luigi.Task):
    pdb = luigi.Parameter()

    def requires(self):
        return {'hs_result': GetHotspots(self.pdb), 'ligands': CSDLigandPrep(self.pdb)}

    def output(self):
        return luigi.LocalTarget("pdb_files/{}/{}/ligand_scores.csv".format(self.pdb[1:3], self.pdb))

    def run(self):
        mols = io.MoleculeReader(self.input()['ligands'].path)
        hr = HotspotReader(self.input()['hs_result'].path).read()

        with open(self.output().path, 'w') as csv_file:
            csv_file.write("mol_id,atom_id,score\n")
            for mol in mols:
                scored_mol = hr.score(mol)
                for a in scored_mol.heavy_atoms:
                    out_str = "{},{},{}\n".format(mol.identifier, a.label, a.partial_charge)
                    csv_file.write(out_str)


class ScoreProtein(luigi.Task):
    pdb = luigi.Parameter()

    def requires(self):
        return {'protein_scores': ScoreProtein(self.pdb), 'hs_result': GetHotspots(self.pdb)}

    def output(self):
        return luigi.LocalTarget("pdb_files/{}/{}/protein_scores.csv".format(self.pdb[1:3], self.pdb))

    def run(self):
        prot = Protein.from_file(self.input()['protein'].path)
        hr = HotspotReader(self.input()['hs_result'].path).read()

        scored_prot = hr.score(prot)

        with open(self.output().path, 'w') as csv_file:
            csv_file.write("mol_id,atom_id,score\n")
            for a in scored_prot.heavy_atoms:
                out_str = "protein,{},{}\n".format(a.label, a.partial_charge)
                csv_file.write(out_str)


class ToJSON(luigi.Task):
    pdb = luigi.Parameter()

    def requires(self):
        return {'protein_scores': ScoreProtein(self.pdb), 'ligand_scores': ScoreLigands(self.pdb)}

    def output(self):
        return luigi.LocalTarget("pdb_files/{}/{}/scores.json".format(self.pdb[1:3], self.pdb))

    def run(self):
        prot_df = pd.read_csv(self.input()['protein_scores'].path)
        lig_df = pd.read_csv(self.input()['ligand_scores'].path)
        df = pd.concat([prot_df, lig_df])
        df = df.groupby('mol_id', as_index=False).apply(lambda x: x[['atom_id', 'score']].to_dict('r')).reset_index()
        df.to_json(self.output().path, orient='records')


class SaveTractabilityMap(luigi.Task):
    pdb = luigi.Parameter()
    volume = luigi.Parameter()

    def requires(self):
        return GetHotspots(self.pdb)

    def output(self):
        return luigi.LocalTarget("pdb_files/{}/{}/tractability_{}/out.zip".format(self.pdb[1:3], self.pdb, self.volume))

    def run(self):
        hr = HotspotReader(self.input().path).read()
        bcv = hr.tractability_map(volume=self.volume)
        out_settings = HotspotWriter.Settings()
        out_settings.charged = False
        w = HotspotWriter(os.path.dirname(self.output().path), grid_extension=".grd", zip_results=True,
                          settings=out_settings)

        w.write(bcv)


class WriteTractabilityScores(luigi.Task):
    pdb = luigi.Parameter()
    volume = luigi.Parameter()

    def requires(self):
        return SaveTractabilityMap(self.pdb, self.volume)

    def output(self):
        return luigi.LocalTarget("pdb_files/{}/{}/tractability_{}/scores.csv".format(self.pdb[1:3], self.pdb,self.volume))

    def run(self):
        hr = [HotspotReader(self.input().path).read()]
        i = 0
        all_cavs = []

        for cav in hr:
            i += 1
            hist = cav.map_values()

            all_points = []
            for x in hist.values():
                all_points += x.flatten().tolist()

            df = pd.DataFrame({'scores': all_points})
            df = df[df['scores'] != 0]
            df['cavity'] = i
            all_cavs.append(df)
        all_df = pd.concat(all_cavs)
        all_df.to_csv(self.output().path)


class LotsOTasks(luigi.WrapperTask):

    def requires(self):
        pdb_li = []
        volumes = [350,400,450,500,550]

        with open('info.csv') as csv_file:
            pdb_info = csv.reader(csv_file, delimiter=',')

            for line in pdb_info:
                pdb = line[0].lower()
                if line[2] == 'v':
                    pdb_li.append(pdb)

        for pdb in pdb_li:
            for volume in volumes:
                yield WriteTractabilityScores(pdb=pdb,volume=volume)



if __name__ == '__main__':
    luigi.build([LotsOTasks()], workers=12)
    # luigi.run()
