import rosetta
from toolbox import pose_from_rcsb
import rosetta.core.scoring
from extract_coords_pose import get_radius
import subprocess
import numpy as np
import Bio.PDB
import PeptideBuilder
import Geometry
import pybel
import random
import math


class MonteCarlo:

    def __init__(self, pdb, mc_temperature=1, identifier=0, local=True):
        """
        :param pdb: PDB id of protein to be folded.
        :param mc_temperature: Temperature for MC simulation.
        :param identifier: ID of process.
        :param local: :type boolean: If True, test moves per residue.
        """
        self.gen_last = 0                                 # Stores generation for which energy score was last calculated
        self.lowest = 100000000                           # Lowest score observed
        self.conformation = []                            # Current conformation of protein
        self.mover_size = 20                              # Degrees to move at each step
        self.current_score = 0                            # The current energy at any state
        self.steps = 0                                    # Number of steps elapsed
        self.max_steps = 1000000                          # Maximum number of steps to run simulation
        self.temperature = mc_temperature                 # Monte carlo temperature
        self.accepted = 1.0                               # Number of MC steps accepted
        self.rejected = 1.0                               # Number of MC steps rejected
        self.id = identifier                              # Id of process, for running multiple simultaneously
        self.threshold = 1                                # Threshold counter for increasing mover range
        self.score_10 = 10000000                          # The lowest score 10 steps ago
        self.local = local                                # Whether to use local mover (global otherwise)
        self.rot_conformation = []                        # Container for rotamer conformations
        self.rot_conf_local = []                          # Container for rotamer conformation within loop
        self.rot_iter = 100                               # Number of iterations to try to resolve side chain clashes
        self.rot_mover_size = 5                           # Size of rotamer mover
        self.new_conf = False                             # Switch to build with new rotamer conformations
        self.mod_dict = {}                                # Dicitonary of modified rotamers

        # Rosetta inits
        # This has ended up only using rosetta for getting the sequence from pdb, which is of course not necessary
        rosetta.init()                                    # Initialize rosetta libraries
        pose_native = pose_from_rcsb(pdb)                 # Create rosetta pose of natively folded protein from pdb file
        self.sequence = pose_native.sequence()            # Get sequence of protein

        self.c_size = len(self.sequence)*2                     # Number of residues * 2 (phi and psi for each residue)

    def local_mover(self, n):
        """ :param candidates: Current population being evaluated.
        """
        conformation = self.conformation
        move = 0

        score = 0.0
        surf_area = 0.0
        phobic_area = 0.0
        philic_area = 0.0
        e_score = 0.0
        
        if self.steps != 0:
            rand = random.random()
            if rand < 0.33:
                move = -self.mover_size
            elif rand < 0.66:
                move = 0
            else:
                move = self.mover_size

            conformation[n] += move

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
        out.save("/home/dan/solus/solus_design/PyRosetta/msms/prot" + str(self.id) + ".pdb")

        mol = pybel.readfile("pdb", "/home/dan/solus/solus_design/PyRosetta/msms/prot" + str(self.id) + ".pdb").next()

        mol.OBMol.AddHydrogens()

        pybelmol = pybel.Molecule(mol)
        pybelmol.write("pdb", "/home/dan/solus/solus_design/PyRosetta/msms/prot" + str(self.id) + ".pdb",
                       overwrite=True)

        subprocess.call('cd ~/solus/solus_design/PyRosetta/msms; ./pdb_to_xyzrn prot' + str(self.id) +
                        '.pdb > prot' + str(self.id) + '.xyzrn', shell=True)

        check = []

        with open("/home/dan/solus/solus_design/PyRosetta/msms/prot" + str(self.id) + ".xyzrn") as f:
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
                self.threshold = 0
                subprocess.call('cd ~/solus/solus_design/PyRosetta/msms; ' +
                                './msms.x86_64Linux2.2.6.1 -if prot' + str(self.id) + '.xyzrn -af ses_atoms' +
                                str(self.id) + ' > /dev/null', shell=True)

                with open("/home/dan/solus/solus_design/PyRosetta/msms/ses_atoms" + str(self.id) + ".area") as f:
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

                print "SA: " + str(surf_area)

                print "SPHOBE: " + str(phobic_area)

                # print "SPHIL: " + str(philic_area)

                # print "ES: " + str(e_score)

                # Score the new conformation
                score = surf_area + phobic_area

                if self.steps == 0:
                    self.current_score = score
                    self.accepted += 1
                else:
                    if score >= self.current_score:
                        try:
                            p_acc = math.exp(-(score - self.current_score)/self.temperature)
                            print "PACC " + str(p_acc)
                        except OverflowError:
                            p_acc = 0
                            print "OVERFLOW"
                        rand = random.random()
                        if p_acc < rand:
                            self.rejected += 1
                            pass
                        else:
                            self.conformation = conformation
                            self.current_score = score
                            self.accepted += 1
                    else:
                        self.conformation = conformation
                        self.current_score = score
                        self.accepted += 1
                        if score < self.lowest:
                            self.lowest = score
                            pybelmol.write("pdb", "lowest" + str(self.id) + ".pdb", overwrite=True)
                            pybelmol.write("pdb", "lowest_backup" + str(self.id) + ".pdb", overwrite=True)

        # For messing with dynamic parameter adjustment
        
        # if self.threshold % 10 == 0:
        #     if self.mover_size > 1:
        #         self.mover_size -= 1
        #         print "MS: " + str(self.mover_size)
        #     self.threshold = 0

        # if self.steps % 20 == 0:
        #     if self.lowest - self.score_10 < 10:
        #         self.mover_size += 1
        #     self.score_10 = self.lowest

        print "Done. \n"

        print "################################"
        print "Step " + str(self.steps) + " | MC Ratio: " + str(self.rejected/self.accepted) + \
              " | Lowest score: " + str(self.lowest)
        print "################################"

        self.steps += 1

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
                out.save("/home/dan/solus/solus_design/PyRosetta/msms/prot" + str(self.id) + ".pdb")

                mol = pybel.readfile("pdb", "/home/dan/solus/solus_design/PyRosetta/msms/prot" + str(self.id) + ".pdb").next()

                mol.OBMol.AddHydrogens()

                pybelmol = pybel.Molecule(mol)
                pybelmol.write("pdb", "/home/dan/solus/solus_design/PyRosetta/msms/prot" + str(self.id) + ".pdb",
                               overwrite=True)

                subprocess.call('cd ~/solus/solus_design/PyRosetta/msms; ./pdb_to_xyzrn prot' + str(self.id) +
                                '.pdb > prot' + str(self.id) + '.xyzrn', shell=True)

                check = []

                with open("/home/dan/solus/solus_design/PyRosetta/msms/prot" + str(self.id) + ".xyzrn") as f:
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
                        self.threshold += 1
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

    def rotamer_mover(self, conformation, n):
        if n == 0:
            self.rot_conf_local = self.rot_conformation
        for i in range(len(self.rot_conf_local)):
            aa = self.rot_conf_local[i]
            moves = []
            for j in range(len(aa)):
                rand = random.random()
                if rand < 0.33:
                    moves.append(-self.rot_mover_size)
                elif rand < 0.66:
                    moves.append(0)
                else:
                    moves.append(self.rot_mover_size)
            aa = [aa[i] + moves[i] for i in range(len(aa))]
            self.rot_conf_local[i] = aa

        geo = Geometry.geometry(self.sequence[0])
        geo.phi = conformation[0]
        geo.psi_im1 = conformation[1]
        geo.inputRotamers(self.rot_conf_local[0])
        structure = PeptideBuilder.initialize_res(geo)

        i = 2
        j = 1
        for aa in self.sequence[1:]:
            geo = Geometry.geometry(aa)
            geo.phi = conformation[i]
            geo.psi_im1 = conformation[i + 1]
            if aa != "G" and aa != "P" and aa != "A":
                geo.inputRotamers(self.rot_conf_local[j])
                j += 1
            structure = PeptideBuilder.add_residue(structure, geo)
            i += 2

        out = Bio.PDB.PDBIO()
        out.set_structure(structure)
        out.save("/home/dan/solus/solus_design/PyRosetta/msms/prot" + str(self.id) + ".pdb")

        mol = pybel.readfile("pdb", "/home/dan/solus/solus_design/PyRosetta/msms/prot" + str(self.id) + ".pdb").next()

        mol.OBMol.AddHydrogens()

        pybelmol = pybel.Molecule(mol)
        pybelmol.write("pdb", "/home/dan/solus/solus_design/PyRosetta/msms/prot" + str(self.id) + ".pdb",
                       overwrite=True)

        subprocess.call('cd ~/solus/solus_design/PyRosetta/msms; ./pdb_to_xyzrn prot' + str(self.id) +
                        '.pdb > prot' + str(self.id) + '.xyzrn', shell=True)

        check = []

        with open("/home/dan/solus/solus_design/PyRosetta/msms/prot" + str(self.id) + ".xyzrn") as f:
            check = f.readlines()

        atoms_check = []

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
                self.threshold += 1
                return False
            atoms_check.append([resid, entries[0], entries[1], entries[2], r, atname, res])

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
                        self.threshold += 1
        self.rot_conformation = self.rot_conf_local
        self.new_conf = True
        print "Success."
        return True

    def rotamer_cycle(self, n_iter, conformation):
        print "Adjusting rotamers..."
        for n in range(n_iter):
            if self.rotamer_mover(conformation, n):
                return True
        return False

    def mover(self):
        """ :param candidates: Current population being evaluated.
        """
        conformation = self.conformation
        moves = []

        score = 100000
        surf_area = 0.0
        phobic_area = 0.0
        philic_area = 0.0
        e_score = 0.0
        
        if self.steps != 0:
            for i in range(self.c_size):
                rand = random.random()
                if rand < 0.5:
                    moves.append(-self.mover_size)
                else:
                    moves.append(self.mover_size)

            conformation = [conformation[i] + moves[i] for i in range(self.c_size)]
            
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
        out.save("/home/dan/solus/solus_design/PyRosetta/msms/prot" + str(self.id) + ".pdb")

        mol = pybel.readfile("pdb", "/home/dan/solus/solus_design/PyRosetta/msms/prot" + str(self.id) + ".pdb").next()

        mol.OBMol.AddHydrogens()

        pybelmol = pybel.Molecule(mol)
        pybelmol.write("pdb", "/home/dan/solus/solus_design/PyRosetta/msms/prot" + str(self.id) + ".pdb",
                       overwrite=True)

        subprocess.call('cd ~/solus/solus_design/PyRosetta/msms; ./pdb_to_xyzrn prot' + str(self.id) +
                        '.pdb > prot' + str(self.id) + '.xyzrn', shell=True)

        check = []

        with open("/home/dan/solus/solus_design/PyRosetta/msms/prot" + str(self.id) + ".xyzrn") as f:
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
                self.threshold += 1
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
                                self.threshold += 1
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
                self.threshold = 0
                subprocess.call('cd ~/solus/solus_design/PyRosetta/msms; ' +
                                './msms.x86_64Linux2.2.6.1 -if prot' + str(self.id) + '.xyzrn -af ses_atoms' +
                                str(self.id) + ' > /dev/null', shell=True)

                with open("/home/dan/solus/solus_design/PyRosetta/msms/ses_atoms" + str(self.id) + ".area") as f:
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

                print "SA: " + str(surf_area)

                print "SPHOBE: " + str(phobic_area)

                # print "SPHIL: " + str(philic_area)

                # print "ES: " + str(e_score)

                # Score the new conformation
                score = surf_area + phobic_area

                if self.steps == 0:
                    self.current_score = score
                    self.accepted += 1
                else:
                    if score > self.current_score:
                        try:
                            p_acc = math.exp(-(self.current_score - score)/self.temperature)
                        except OverflowError:
                            p_acc = 0
                        rand = random.random()
                        if p_acc > rand:
                            self.rejected += 1
                            pass
                        else:
                            self.conformation = conformation
                            self.current_score = score
                            self.accepted += 1
                    else:
                        self.conformation = conformation
                        self.current_score = score
                        self.accepted += 1
                        if score < self.lowest:
                            self.lowest = score
                            pybelmol.write("pdb", "lowest" + str(self.id) + ".pdb", overwrite=True)
                            pybelmol.write("pdb", "lowest_backup" + str(self.id) + ".pdb", overwrite=True)

        if self.threshold % 10 == 0:
            if self.mover_size > 1:
                self.mover_size -= 1
            self.threshold = 0

        if self.steps % 10 == 0:
            if self.lowest - self.score_10 < 10:
                self.mover_size += 1
            self.score_10 = self.lowest

        print "Done. \n"

        print "################################"
        print "Step " + str(self.steps) + " | MC Ratio: " + str(self.rejected/self.accepted) + \
              " | Lowest score: " + str(self.lowest)
        print "################################"

        self.steps += 1

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

    def run(self):
        """ Run MC simulation.
        """

        seed = []

        for aa in self.sequence:
            geo = Geometry.geometry(aa)
            seed.append(geo.phi)
            seed.append(geo.psi_im1)

        self.conformation = seed

        if self.local:
            while self.steps < self.max_steps:
                for n in range(self.c_size):
                    self.local_mover(n)
        else:
            while self.steps < self.max_steps:
                self.mover()
