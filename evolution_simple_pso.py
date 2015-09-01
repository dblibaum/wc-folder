from random import Random
from time import time
from ecspy import swarm
from ecspy import ec
from ecspy import terminators
from ecspy import topologies
import rosetta
from toolbox import pose_from_rcsb
import rosetta.core.scoring
from extract_coords_pose import get_radius
import subprocess
import numpy as np
import itertools
import Bio.PDB
import PeptideBuilder
import Geometry
import pybel
from os import walk
from ramachandran_biopython import get_phi_psi
import random


class Evolution:

    def __init__(self, pdb, identifier=0, path=""):
        """ :param pdb: :type string: pdb ID of the protein to be folded.
        """
        self.gen_last = 0                                 # Stores generation for which energy score was last calculated
        self.lowest = 100000000                           # Lowest score observed

        # Rosetta inits
        # This has ended up only using rosetta to get the sequence of the pdb file, which is of course not necessary
        rosetta.init()                                    # Initialize rosetta libraries
        pose_native = pose_from_rcsb(pdb)                 # Create rosetta pose of natively folded protein from pdb file

        self.sequence = pose_native.sequence()            # Get sequence of protein
        self.id = identifier                              # Id of process
        self.rot_iter = 200                               # Number of iterations to try to resolve side chain clashes
        self.rot_mover_size = 5                           # Size of rotamer mover
        self.new_conf = False                             # Switch to build with new rotamer conformations
        self.mod_dict = {}                                # Dictionary of modified rotamers
        self.path = path

        self.c_size = len(self.sequence)*2                     # Number of residues * 2 (phi and psi for each residue)

        # Ecspy inits
        rand = Random()
        rand.seed(int(time()))
        self.es = swarm.PSO(rand)                         # Create ecspy evolution strategy seeded with current time
        self.bounder = ec.Bounder(-180, 180)
        self.es.topology = topologies.star_topology

    def generator(self, random, args):
        c_size = args.get('num_inputs', 10)
        return [random.randint(-180, 180) for x in range(c_size)]

    # Novelty criterion
    def eval_func(self, candidates, args):
        """ :param candidates: Current population being evaluated.
        """

        scores_list = []

        scores = []

        print "Calculating scores... "
        for conformation in candidates:

            score = 0.0
            surf_area = 0.0
            phobic_area = 0.0
            philic_area = 0.0
            e_score = 0.0
            
            # Create PeptideBuilder structure for conformation

            geo = Geometry.geometry(self.sequence[0])
            geo.phi = conformation[0]
            geo.psi_im1 = conformation[1]
            if self.sequence[0] != "G" and self.sequence[0] != "P" and self.sequence[0] != "A":
                if 0 in self.mod_dict:
                    geo.inputRotamers(self.mod_dict[0])
            structure = PeptideBuilder.initialize_res(geo)

            i = 2
            j = 1
            for aa in self.sequence[1:]:
                geo = Geometry.geometry(aa)
                geo.phi = conformation[i]
                geo.psi_im1 = conformation[i + 1]
                if aa != "G" and aa != "P" and aa != "A":
                    if j in self.mod_dict:
                        geo.inputRotamers(self.mod_dict[j])
                j += 1
                structure = PeptideBuilder.add_residue(structure, geo)
                i += 2

            out = Bio.PDB.PDBIO()
            out.set_structure(structure)
            out.save(self.path + "/msms/prot" + str(self.id) + ".pdb")
            
            # Add hydrogens with pybel

            mol = pybel.readfile("pdb", self.path + "/msms/prot" + str(self.id) + ".pdb").next()

            mol.OBMol.AddHydrogens()

            pybelmol = pybel.Molecule(mol)
            pybelmol.write("pdb", self.path + "/msms/prot" + str(self.id) + ".pdb",
                           overwrite=True)

            # Convert to xyzrn with msms executable

            subprocess.call('cd ' + self.path + '/msms; ./pdb_to_xyzrn prot' + str(self.id) +
                            '.pdb > prot' + str(self.id) + '.xyzrn', shell=True)
                            
            # Enforce constraints (atoms can't come closer than van-der-waals radii)

            check = []

            with open(self.path + "/msms/prot" + str(self.id) + ".xyzrn") as f:
                check = f.readlines()

            atoms_check = []

            unk = False
            for atom in check:
                entries = atom.split()
                name = entries[5]
                ent = name.split("_")
                resid = ent[2]
                atname = ent[0]
                res = ent[1]
                try:
                    r = get_radius(atname, res)[0]
                except KeyError:
                    unk = True
                    score = 1000000
                    break
                atoms_check.append([resid, entries[0], entries[1], entries[2], r, atname, res])

            if not unk:
                clash = False
                for atom in atoms_check:
                    id1 = int(atom[0])
                    atname1 = atom[5]
                    x = float(atom[1])
                    y = float(atom[2])
                    z = float(atom[3])
                    r = float(atom[4])
                    for atom2 in atoms_check:
                        id2 = int(atom2[0])
                        atname2 = atom2[5]
                        if (id1 != id2) and not (((id2 == (id1 + 1)) or (id2 == (id1 - 1)))
                                                 and (((atname1 == "CA") or (atname1 == "C") or (atname1 == "N")) and
                                                      ((atname2 == "CA") or (atname2 == "C") or (atname2 == "N")))):
                            x2 = float(atom2[1])
                            y2 = float(atom2[2])
                            z2 = float(atom2[3])
                            r2 = float(atom2[4])
                            distance = np.sqrt((x2 - x)**2 + (y2 - y)**2 + (z2 - z)**2)
                            if distance < r + r2:
                                solved = self.rebuild(id1, id2, self.rot_iter, conformation)
                                if not solved:
                                    score = 1000000
                                    clash = True
                                    print "CLASH"
                                    print "> " + str(atname1) + " " + str(atname2) + " " + str(distance) + " " + str(id1) + " " \
                                          + str(id2)
                                    break
                                else:
                                    print "Solved."
                    if clash:
                        break

                if not clash:
                    
                    # Compute ses for all atoms with msms executable
                    subprocess.call('cd ' + self.path + '/msms; ' +
                                    './msms.x86_64Linux2.2.6.1 -if prot' + str(self.id) + '.xyzrn -af ses_atoms' +
                                    str(self.id) + ' > /dev/null', shell=True)

                    with open(self.path + "/msms/ses_atoms" + str(self.id) + ".area") as f:
                        atoms = f.readlines()

                    ses_type = []
                    for i in range(len(atoms)):
                        if i != 0:
                            entries = atoms[i].split()
                            ses = entries[1]
                            atom_desc = entries[3]
                            atmname = atom_desc.split("_")
                            rad_hydro = get_radius(atmname[0].strip(), atmname[1].strip())
                            res_num = atmname[2].strip()
                            atm_type = rad_hydro[1]
                            ses_type.append([ses, atm_type, res_num])

                    for atom in ses_type:
                        surf_area += float(atom[0])
                        if atom[1] == 1:
                            phobic_area += float(atom[0])
                        if atom[1] == 2:
                            philic_area -= float(atom[0])

                    print "SA: " + str(surf_area/2)

                    print "SPHOBE: " + str(phobic_area)

                    print "SPHIL: " + str(philic_area)

                    # print "ES: " + str(e_score)

                    # Score the new conformation
                    score = surf_area/2 + phobic_area + philic_area

            scores.append([score, conformation])

            # For fitness
            scores_list.append(score)

        print "Done. \n"

        scores = sorted(scores)

        # Save lowest to pdb
        if scores[0][0] < self.lowest:
            conformation = scores[0][1]

            geo = Geometry.geometry(self.sequence[0])
            geo.phi = conformation[0]
            geo.psi_im1 = conformation[1]
            if self.sequence[0] != "G" and self.sequence[0] != "P" and self.sequence[0] != "A":
                if 0 in self.mod_dict:
                    geo.inputRotamers(self.mod_dict[0])
            structure = PeptideBuilder.initialize_res(geo)

            i = 2
            j = 1
            for aa in self.sequence[1:]:
                geo = Geometry.geometry(aa)
                geo.phi = conformation[i]
                geo.psi_im1 = conformation[i + 1]
                if aa != "G" and aa != "P" and aa != "A":
                    if j in self.mod_dict:
                        geo.inputRotamers(self.mod_dict[j])
                j += 1
                structure = PeptideBuilder.add_residue(structure, geo)
                i += 2

            out = Bio.PDB.PDBIO()
            out.set_structure(structure)
            out.save("lowest" + str(self.id) + ".pdb")

            mol = pybel.readfile("pdb", "lowest" + str(self.id) + ".pdb").next()

            mol.OBMol.AddHydrogens()

            pybelmol = pybel.Molecule(mol)
            pybelmol.write("pdb", "lowest" + str(self.id) + ".pdb", overwrite=True)
            pybelmol.write("pdb", "lowest_backup" + str(self.id) + ".pdb", overwrite=True)

        if scores[0][0] < self.lowest:
            self.lowest = scores[0][0]

        print "Lowest score: " + str(self.lowest)

        # return novelty_scores
        return scores_list

    def rebuild(self, id1, id2, it, conformation):
        """Rebuilds clashing residues.
        :param id1: First residue id.
        :param id2: Second residue id.
        :param it: Number of iterations.
        :param conformation: Current conformation.
        :return: :type boolean: True if successfully rebuilt, False if not.
        """
        print "Adjusting rotamers..."
        new_rot1 = []
        new_rot2 = []
        aa1 = self.sequence[id1 - 1]
        aa2 = self.sequence[id2 - 1]
        solved = False
        if (aa1 != "G" and aa1 != "P" and aa1 != "A") or (aa2 != "G" and aa2 != "P" and aa2 != "A"):
            for n in range(it):
                if solved:
                    break
                if n == 0:
                    new_rot1 = self.make_rot_list(aa1)
                    new_rot2 = self.make_rot_list(aa2)
                else:
                    moves1 = []
                    for j in range(len(new_rot1)):
                        rand = random.random()
                        if rand < 0.33:
                            moves1.append(-self.rot_mover_size)
                        elif rand < 0.66:
                            moves1.append(0)
                        else:
                            moves1.append(self.rot_mover_size)
                    new_rot1 = [new_rot1[i] + moves1[i] for i in range(len(new_rot1))]

                    moves2 = []
                    for j in range(len(new_rot2)):
                        rand = random.random()
                        if rand < 0.33:
                            moves2.append(-self.rot_mover_size)
                        elif rand < 0.66:
                            moves2.append(0)
                        else:
                            moves2.append(self.rot_mover_size)
                    new_rot2 = [new_rot2[i] + moves2[i] for i in range(len(new_rot2))]

                geo = Geometry.geometry(self.sequence[0])
                geo.phi = conformation[0]
                geo.psi_im1 = conformation[1]
                if id1 - 1 == 0:
                    if aa1 != "G" and aa1 != "P" and aa1 != "A":
                        geo.inputRotamers(new_rot1)
                elif id2 - 1 == 0:
                    if aa2 != "G" and aa2 != "P" and aa2 != "A":
                        geo.inputRotamers(new_rot2)
                elif 0 in self.mod_dict:
                    if self.sequence[0] != "G" and self.sequence[0] != "P" and self.sequence[0] != "A":
                        geo.inputRotamers(self.mod_dict[0])
                structure = PeptideBuilder.initialize_res(geo)

                i = 2
                j = 1
                for aa in self.sequence[1:]:
                    geo = Geometry.geometry(aa)
                    geo.phi = conformation[i]
                    geo.psi_im1 = conformation[i + 1]
                    if id1 == j + 1:
                        if aa1 != "G" and aa1 != "P" and aa1 != "A":
                            geo.inputRotamers(new_rot1)
                    elif id2 == j + 1:
                        if aa2 != "G" and aa2 != "P" and aa2 != "A":
                            geo.inputRotamers(new_rot2)
                    elif j in self.mod_dict:
                        geo.inputRotamers(self.mod_dict[j])
                    j += 1
                    structure = PeptideBuilder.add_residue(structure, geo)
                    i += 2

                out = Bio.PDB.PDBIO()
                out.set_structure(structure)
                out.save(self.path + "/msms/prot" + str(self.id) + ".pdb")

                mol = pybel.readfile("pdb", self.path + "/msms/prot" + str(self.id) + ".pdb").next()

                mol.OBMol.AddHydrogens()

                pybelmol = pybel.Molecule(mol)
                pybelmol.write("pdb", self.path + "/msms/prot" + str(self.id) + ".pdb",
                               overwrite=True)

                subprocess.call('cd ' + self.path + '/msms; ./pdb_to_xyzrn prot' + str(self.id) +
                                '.pdb > prot' + str(self.id) + '.xyzrn', shell=True)

                check = []

                with open(self.path + "/msms/prot" + str(self.id) + ".xyzrn") as f:
                    check = f.readlines()

                atoms_check = []

                clash = False
                unk = False
                for atom in check:
                    entries = atom.split()
                    name = entries[5]
                    ent = name.split("_")
                    resid = ent[2]
                    atname = ent[0]
                    res = ent[1]
                    try:
                        r = get_radius(atname, res)[0]
                    except KeyError:
                        unk = True
                        break
                    atoms_check.append([resid, entries[0], entries[1], entries[2], r, atname, res])

                if not unk:
                    for atom in atoms_check:
                        aid1 = int(atom[0])
                        atname1 = atom[5]
                        x = float(atom[1])
                        y = float(atom[2])
                        z = float(atom[3])
                        r = float(atom[4])
                        for atom2 in atoms_check:
                            aid2 = int(atom2[0])
                            atname2 = atom2[5]
                            if (aid1 != aid2) and not (((aid2 == (aid1 + 1)) or (aid2 == (aid1 - 1)))
                                                    and (((atname1 == "CA") or (atname1 == "C") or (atname1 == "N")) and
                                                         ((atname2 == "CA") or (atname2 == "C") or (atname2 == "N")))):
                                x2 = float(atom2[1])
                                y2 = float(atom2[2])
                                z2 = float(atom2[3])
                                r2 = float(atom2[4])
                                distance = np.sqrt((x2 - x)**2 + (y2 - y)**2 + (z2 - z)**2)
                                if distance < r + r2:
                                    clash = True
                                    break
                        if clash:
                            break
                    if not clash:
                        solved = True
        else:
            print str(aa1) + " " + str(aa2)
        if solved:
            if aa1 != "G" and aa1 != "P" and aa1 != "A":
                self.mod_dict[id1 - 1] = new_rot1
            if aa2 != "G" and aa2 != "P" and aa2 != "A":
                self.mod_dict[id2 - 1] = new_rot2
            return True
        else:
            return False

    def make_rot_list(self, aa):
        rotamers = []
        if aa == "R":
            for i in range(6):
                rand = random.randint(-180, 180)
                rotamers.append(rand)
        elif aa == "N":
            for i in range(2):
                rand = random.randint(-180, 180)
                rotamers.append(rand)
        elif aa == "D":
            for i in range(2):
                rand = random.randint(-180, 180)
                rotamers.append(rand)
        elif aa == "C":
            for i in range(1):
                rand = random.randint(-180, 180)
                rotamers.append(rand)
        elif aa == "Q":
            for i in range(3):
                rand = random.randint(-180, 180)
                rotamers.append(rand)
        elif aa == "E":
            for i in range(3):
                rand = random.randint(-180, 180)
                rotamers.append(rand)
        elif aa == "H":
            for i in range(2):
                rand = random.randint(-180, 180)
                rotamers.append(rand)
        elif aa == "I":
            for i in range(3):
                rand = random.randint(-180, 180)
                rotamers.append(rand)
        elif aa == "L":
            for i in range(3):
                rand = random.randint(-180, 180)
                rotamers.append(rand)
        elif aa == "K":
            for i in range(4):
                rand = random.randint(-180, 180)
                rotamers.append(rand)
        elif aa == "M":
            for i in range(3):
                rand = random.randint(-180, 180)
                rotamers.append(rand)
        elif aa == "F":
            for i in range(2):
                rand = random.randint(-180, 180)
                rotamers.append(rand)
        elif aa == "P":
            for i in range(1):
                rand = random.randint(-180, 180)
                rotamers.append(rand)
        elif aa == "S":
            for i in range(1):
                rand = random.randint(-180, 180)
                rotamers.append(rand)
        elif aa == "T":
            for i in range(2):
                rand = random.randint(-180, 180)
                rotamers.append(rand)
        elif aa == "W":
            for i in range(2):
                rand = random.randint(-180, 180)
                rotamers.append(rand)
        elif aa == "Y":
            for i in range(2):
                rand = random.randint(-180, 180)
                rotamers.append(rand)
        elif aa == "V":
            for i in range(2):
                rand = random.randint(-180, 180)
                rotamers.append(rand)
        else:
            pass
        return rotamers

    def create_seeds(self, directory):
        f = []
        seeds = []
        for (dirpath, dirnames, filenames) in walk(directory):
            f.extend(filenames)
            break
        for name in f:
            pdb = name.split(".")[0]
            seeds.append(get_phi_psi(pdb))
        seeds = [seed[:self.c_size] for seed in seeds if len(seed) >= self.c_size]
        return seeds

    def evolve(self):
        """ Run evolution.
        """
        pop_size = 100

        seeds = []

        seed = []
        for aa in self.sequence:
            geo = Geometry.geometry(aa)
            seed.append(geo.phi)
            seed.append(geo.psi_im1)

        template_seeds = self.create_seeds(self.path + "/pdbs2")

        seeds += template_seeds

        seeds += [seed for x in range(100 - len(template_seeds))]

        self.es.terminator = terminators.evaluation_termination
        self.es.evolve(generator=self.generator,
                       evaluator=self.eval_func,
                       pop_size=pop_size,
                       maximize=False,
                       max_evaluations=2000000000,
                       bounder=self.bounder,
                       seeds=seeds,
                       neighborhood_size=10,
                       cognitive_rate=1,
                       social_rate=1,
                       num_inputs=self.c_size)
