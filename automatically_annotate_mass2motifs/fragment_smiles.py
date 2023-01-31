from rdkit import Chem
import numpy
import sys


class Constants:
    def __init__(self):
        self.typew = {Chem.rdchem.BondType.names["AROMATIC"]: 3.0,
                      Chem.rdchem.BondType.names["DOUBLE"]: 2.0,
                      Chem.rdchem.BondType.names["TRIPLE"]: 3.0,
                      Chem.rdchem.BondType.names["SINGLE"]: 1.0}
        self.missingfragmentpenalty = 10
        self.heterow = {False: 2, True: 1}
        self.mims = {'H': 1.0078250321,
                     'C': 12.0000000,
                     'N': 14.0030740052,
                     'O': 15.9949146221,
                     'F': 18.99840320,
                     'Na': 22.9897692809,
                     'P': 30.97376151,
                     'S': 31.97207069,
                     'Cl': 34.96885271,
                     'K': 38.96370668,
                     'Br': 78.9183376,
                     'I': 126.904468}
        self.Hmass = self.mims["H"]
        self.elmass = 0.0005486
        self.ionmasses = {1: {'+H': self.mims['H'],
                              '+NH4': self.mims['N'] + 4 * self.mims['H'],
                              '+Na': self.mims['Na'],
                              '+K': self.mims['K']},
                          -1: {'-H': -self.mims['H'],
                               '+Cl': self.mims['Cl']},
                          }


class Molecule:
    """Stores all information for a molecule"""
    def __init__(self, smile, constants: Constants):
        try:
            self.mol = Chem.MolFromSmiles(smile)
            self.accept = True
        except:
            self.accept = False
            return
        self.natoms = self.mol.GetNumAtoms()
        self.atom_masses = []
        self.atomHs = []
        self.neutral_loss_atoms = []
        self.bonded_atoms = []  # [[list of atom numbers]]
        self.bonds = set([])
        self.bondscore = {}
        self.constants = constants
        for x in range(self.natoms):
            # saves information about
            self.bonded_atoms.append([])
            atom = self.mol.GetAtomWithIdx(x)
            self.atomHs.append(atom.GetNumImplicitHs() + atom.GetNumExplicitHs())
            self.atom_masses.append(self.constants.mims[atom.GetSymbol()] + self.constants.Hmass * (self.atomHs[x]))
            if atom.GetSymbol() == 'O' and self.atomHs[x] == 1 and len(atom.GetBonds()) == 1:
                self.neutral_loss_atoms.append(x)
            if atom.GetSymbol() == 'N' and self.atomHs[x] == 2 and len(atom.GetBonds()) == 1:
                self.neutral_loss_atoms.append(x)
        for bond in self.mol.GetBonds():
            a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            self.bonded_atoms[a1].append(a2)
            self.bonded_atoms[a2].append(a1)
            bondbits = 1 << a1 | 1 << a2
            bondscore = self.constants.typew[bond.GetBondType()] * \
                        self.constants.heterow[bond.GetBeginAtom().GetSymbol() != 'C' or bond.GetEndAtom().GetSymbol() != 'C']
            self.bonds.add(bondbits)
            self.bondscore[bondbits] = bondscore


class FragmentEngine:
    """Code based on magma: https://github.com/NLeSC/MAGMa"""
    def __init__(self, smile: str, max_broken_bonds, max_water_losses, ionisation_mode, molcharge):
        self.constants = Constants()
        self.molecule = Molecule(smile, self.constants)

        self.max_broken_bonds = max_broken_bonds
        self.max_water_losses = max_water_losses
        self.ionisation_mode = ionisation_mode
        self.molcharge = molcharge
        self.new_fragment = 0
        self.template_fragment = 0
        # Predefine a list with zeroes, that is slowly filled with each generated fragment. I am not sure why this has length max_broken_bronds
        self.fragment_masses = ((max_broken_bonds + max_water_losses) * 2 + 1) * [0]
        self.fragment_info = [[0, 0, 0]]

    def extend(self, atom):
        for a in self.molecule.bonded_atoms[atom]:
            atombit = 1 << a
            if atombit & self.template_fragment and not atombit & self.new_fragment:
                self.new_fragment = self.new_fragment | atombit
                self.extend(a)

    def generate_fragments(self):
        """Creates all fragments in a molecule

        The method of calculating the fragments is by creating fragment numbers which are a binary number
        1000, would mean a molecule """
        # In binary frag is the representation of the molecule. The entire molecule for a molecule of 5 atoms is 11111.
        # Adds the fragment representing the entire molecule
        frag = (1 << self.molecule.natoms) - 1
        self.add_fragment(frag, self.calc_fragment_mass(frag), 0, 0)
        all_fragments = {frag}
        total_fragments = {frag}
        current_fragments = {frag}
        new_fragments = {frag}

        # generate fragments for max_broken_bond steps
        for step in range(self.max_broken_bonds):
            # loop over all fragments to be fragmented
            for fragment in current_fragments:
                # loop over all atoms
                for atom in range(self.molecule.natoms):
                    # check if atom in the fragment
                    if (1 << atom) & fragment:
                        # remove the atom
                        self.template_fragment = fragment ^ (1 << atom)
                        list_ext_atoms = set([])
                        extended_fragments = set([])
                        # find all its neighbor atoms
                        for a in self.molecule.bonded_atoms[atom]:
                            # present in the fragment
                            if (1 << a) & self.template_fragment:
                                list_ext_atoms.add(a)
                        # in case of one bonded atom, the new fragment is the remainder of the old fragment
                        if len(list_ext_atoms) == 1:
                            extended_fragments.add(self.template_fragment)
                        else:
                            # otherwise extend each neighbor atom to a complete fragment
                            for a in list_ext_atoms:
                                # except when deleted atom is in a ring and a previous extended
                                # fragment already contains this neighbor atom, then
                                # calculate fragment only once
                                for frag in extended_fragments:
                                    if (1 << a) & frag:
                                        break
                                else:
                                    # extend atom to complete fragment
                                    self.new_fragment = 1 << a
                                    self.extend(a)
                                    extended_fragments.add(self.new_fragment)
                        for frag in extended_fragments:
                            # add extended fragments, if not yet present, to the collection
                            if frag not in all_fragments:
                                all_fragments.add(frag)
                                bondbreaks, score = self.score_fragment(frag)
                                if bondbreaks <= self.max_broken_bonds and score < (self.constants.missingfragmentpenalty + 5):
                                    new_fragments.add(frag)
                                    total_fragments.add(frag)
                                    self.add_fragment(
                                        frag, self.calc_fragment_mass(frag), score, bondbreaks)
            current_fragments = new_fragments
            new_fragments = set([])
        # number of OH losses
        for step in range(self.max_water_losses):
            # loop of all fragments
            for fi in self.fragment_info:
                # on which to apply neutral loss rules
                if fi[2] == self.max_broken_bonds + step:
                    fragment = fi[0]
                    # loop over all atoms in the fragment
                    for atom in self.molecule.neutral_loss_atoms:
                        if (1 << atom) & fragment:
                            frag = fragment ^ (1 << atom)
                            # add extended fragments, if not yet present, to the collection
                            if frag not in total_fragments:
                                total_fragments.add(frag)
                                bondbreaks, score = self.score_fragment(frag)
                                if score < (self.constants.missingfragmentpenalty + 5):
                                    self.add_fragment(
                                        frag, self.calc_fragment_mass(frag), score, bondbreaks)
        self.convert_fragments_table()
        return len(self.fragment_info)

    def score_fragment(self, fragment):
        score = 0
        bondbreaks = 0
        for bond in self.molecule.bonds:
            if 0 < (fragment & bond) < bond:
                score += self.molecule.bondscore[bond]
                bondbreaks += 1
        if score == 0:
            print("score=0: ", fragment, bondbreaks)
        return bondbreaks, score

    def calc_fragment_mass(self, fragment):
        fragment_mass = 0.0
        for atom in range(self.molecule.natoms):
            if fragment & (1 << atom):
                fragment_mass += self.molecule.atom_masses[atom]
        return fragment_mass

    def add_fragment(self, fragment, fragmentmass, score, bondbreaks):
        mass_range = (
                (self.max_broken_bonds + self.max_water_losses - bondbreaks) * [0] +
                list(numpy.arange(
                    -bondbreaks + self.ionisation_mode * (1 - self.molcharge),
                    bondbreaks + self.ionisation_mode * (1 - self.molcharge) + 1) * self.constants.Hmass + fragmentmass) +
                      (self.max_broken_bonds + self.max_water_losses - bondbreaks) * [0])
        if bondbreaks == 0:
            # make sure that fragmentmass is included
            mass_range[self.max_broken_bonds + self.max_water_losses -
                       self.ionisation_mode] = fragmentmass
        # Adds the new mass fragments to fragment masses
        self.fragment_masses += mass_range
        self.fragment_info.append([fragment, score, bondbreaks])

    def convert_fragments_table(self):
        """Not sure what happens here a numpy 2d array is created with a specific, it contains all the masses of the
        fragments and probably helps with searching later, however it is currently not clear to me why it is full of 0's """
        self.fragment_masses_np = \
            numpy.array(self.fragment_masses).reshape(len(self.fragment_info),
                                                      (self.max_broken_bonds + self.max_water_losses) * 2 + 1)

    def find_fragments(self, mass, precision, mz_precision_abs):
        result = numpy.where(numpy.where(self.fragment_masses_np < max(mass * precision, mass + mz_precision_abs),
                                         self.fragment_masses_np, 0) > min(mass / precision, mass - mz_precision_abs))
        fragment_set = []
        for i in range(len(result[0])):
            fragment_id = result[0][i]
            fragment_set.append(self.fragment_info[fragment_id] +
                                [self.fragment_masses_np[fragment_id][
                                     self.max_broken_bonds + self.max_water_losses - self.ionisation_mode * (
                                                 1 - self.molcharge)]] +
                                [self.ionisation_mode * (1 - self.molcharge) + result[1][
                                    i] - self.max_broken_bonds - self.max_water_losses])
        return fragment_set

    def get_fragment_info(self, fragment):
        atomlist = []
        elements = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'F': 0,
                    'P': 0, 'S': 0, 'Cl': 0, 'Br': 0, 'I': 0}
        for atom in range(self.molecule.natoms):
            if ((1 << atom) & fragment):
                atomlist.append(atom)
                elements[self.molecule.mol.GetAtomWithIdx(atom).GetSymbol()] += 1
                elements['H'] += self.molecule.atomHs[atom]
        formula = ''
        for el in ('C', 'H', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I'):
            nel = elements[el]
            if nel > 0:
                formula += el
            if nel > 1:
                formula += str(nel)
        atomstring = ','.join(str(a) for a in atomlist)
        return atomstring, atomlist, formula, fragment2inchikey(self.molecule.mol, atomlist)

    def get_natoms(self):
        return self.molecule.natoms

    def accepted(self):
        return self.molecule.accept

def fragment2inchikey(mol, atomlist):
    emol = Chem.EditableMol(mol)
    for atom in reversed(range(mol.GetNumAtoms())):
        if atom not in atomlist:
            emol.RemoveAtom(atom)
    frag = emol.GetMol()
    return Chem.MolToSmiles(frag)

def fragment_molecule(smile: str):
    """Fragments a molecule in all possible fragments

    :returns
    A dictionary with the smiles and the substructures
    """
    fragmenter = FragmentEngine(smile,
                                max_broken_bonds=100,
                                max_water_losses=1,
                                ionisation_mode=1,
                                molcharge=1)
    print(fragmenter.accepted())
    print(fragmenter.generate_fragments())
    print(fragmenter.fragment_info)
    # print(fragmenter.fragment_masses)
    print(fragmenter.find_fragments(58, 0.1, 0.1))
    # print(fragmenter.fragment_masses_np)
    for fragment, score, bondbreaks in fragmenter.fragment_info:
    #     print("{0:b}".format(fragment))
    #     print(fragmenter.get_fragment_info(fragment))
        print(fragmenter.calc_fragment_mass(fragment))


if __name__ == "__main__":
    # fragment_molecule("O1-C(-C-O)-C(-O)-C(-O)-C(-O)-C-1-O-C-C-C-C-C(-C-O)-C(-O)-C(-O)-C(-O)-C-O-C")
    fragment_molecule("CC(-C)C")
    # print(Molecule("CC(-C)C").__dict__)