rd = {"C": .77, "CR": .72, "CD": .67, "O": .67, "OD": .60, "N": .70, "ND": .62, "H": .37, "S": 1.04}
hd = {"ALA": 1, "VAL": 1, "ILE": 1, "LEU": 1, "MET": 1, "PHE": 1, "TYR": 1, "TRP": 1, "CYS": 1, "LYS": 2, "ARG": 2,
      "HIS": 0, "ASP": 2, "SER": 0, "THR": 0, "GLU": 2, "ASN": 2, "GLN": 0, "PRO": 2, "GLY": 0}


def get_radius(atom, res):

    if res == "ASN":
        if atom == "OD1":
            return [rd["OD"], hd[res]]
        if atom == "CG":
            return [rd["CD"], hd[res]]
        if atom == "C":
            return [rd["CD"], hd[res]]
        if atom == "O":
            return [rd["OD"], hd[res]]
        if atom[0] == "1" or "2":
            return [rd["H"], hd[res]]
        return [rd[atom[0]], hd[res]]
    elif res == "GLN":
        if atom == "OE1":
            return [rd["OD"], hd[res]]
        if atom == "CD":
            return [rd["CD"], hd[res]]
        if atom == "C":
            return [rd["CD"], hd[res]]
        if atom == "O":
            return [rd["OD"], hd[res]]
        if atom[0] == "1" or "2":
            return [rd["H"], hd[res]]
        return [rd[atom[0]], hd[res]]
    elif res == "ARG":
        if atom == "CZ":
            return [rd["CD"], hd[res]]
        if atom == "NH1":
            return [rd["ND"], hd[res]]
        if atom == "C":
            return [rd["CD"], hd[res]]
        if atom == "O":
            return [rd["OD"], hd[res]]
        if atom[0] == "1" or "2":
            return [rd["H"], hd[res]]
        return [rd[atom[0]], hd[res]]
    elif res == "HIS":
        if atom == "NE2":
            return [rd["ND"], hd[res]]
        if atom == ("CD2" or "CE1" or "CG"):
            return [rd["CD"], hd[res]]
        if atom == "C":
            return [rd["CD"], hd[res]]
        if atom == "O":
            return [rd["OD"], hd[res]]
        if atom[0] == "1" or "2":
            return [rd["H"], hd[res]]
        return [rd[atom[0]], hd[res]]
    elif res == "ASP":
        if atom == "OD1":
            return [rd["OD"], hd[res]]
        if atom == "CG":
            return [rd["CD"], hd[res]]
        if atom == "C":
            return [rd["CD"], hd[res]]
        if atom == "O":
            return [rd["OD"], hd[res]]
        if atom[0] == "1" or "2":
            return [rd["H"], hd[res]]
        return [rd[atom[0]], hd[res]]
    elif res == "GLU":
        if atom == "OE1":
            return [rd["OD"], hd[res]]
        if atom == "CD":
            return [rd["CD"], hd[res]]
        if atom == "C":
            return [rd["CD"], hd[res]]
        if atom == "O":
            return [rd["OD"], hd[res]]
        if atom[0] == "1" or "2":
            return [rd["H"], hd[res]]
        return [rd[atom[0]], hd[res]]
    elif res == ("PHE" or "TYR"):
        if atom == ("CG" or "CE1" or "CE2"):
            return [rd["CD"], hd[res]]
        if atom == ("CD1" or "CD2" or "CZ"):
            return [rd["CR"], hd[res]]
        if atom == "C":
            return [rd["CD"], hd[res]]
        if atom == "O":
            return [rd["OD"], hd[res]]
        if atom[0] == "1" or "2":
            return [rd["H"], hd[res]]
        return [rd[atom[0]], hd[res]]
    elif res == "TRP":
        if atom == ("CD1" or "CE2" or "CG" or "CE3" or "CH2"):
            return [rd["CD"], hd[res]]
        if atom == ("CD2" or "CZ3" or "CZ2"):
            return [rd["CR"], hd[res]]
        if atom == "C":
            return [rd["CD"], hd[res]]
        if atom == "O":
            return [rd["OD"], hd[res]]
        if atom[0] == "1" or "2":
            return [rd["H"], hd[res]]
        return [rd[atom[0]], hd[res]]
    else:
        if atom == "C":
            return [rd["CD"], hd[res]]
        if atom == "O":
            return [rd["OD"], hd[res]]
        if atom[0] == "1" or "2":
            return [rd["H"], hd[res]]
        return [rd[atom[0]], hd[res]]
