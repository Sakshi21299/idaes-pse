from abc import ABC

from .bblock import BBlock


class Atom(BBlock, ABC):
    """A class for representing Atoms."""

    NumberToSymbol = {None: '', 1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N',
                      8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si',
                      15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc',
                      22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni',
                      29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br',
                      36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo',
                      43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In',
                      50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba',
                      57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu',
                      64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
                      71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir',
                      78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po',
                      85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa',
                      92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf',
                      99: 'Es', 100: 'Fm', 101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf',
                      105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Uun',
                      111: 'Uuu', 112: 'Uub', 113: 'Uut', 114: 'Uuq', 115: 'Uup',
                      116: 'Uuh', 117: 'Uus', 118: 'Uuo'}

    SymbolToNumber = {'': None, 'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7,
                      'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14,
                      'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21,
                      'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28,
                      'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35,
                      'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42,
                      'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49,
                      'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56,
                      'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63,
                      'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70,
                      'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77,
                      'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84,
                      'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91,
                      'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98,
                      'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104,
                      'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Uun': 110,
                      'Uuu': 111, 'Uub': 112, 'Uut': 113, 'Uuq': 114, 'Uup': 115,
                      'Uuh': 116, 'Uus': 117, 'Uuo': 118}

    NumberToMass = {0: 0.0, 1: 1.00794, 2: 4.002602, 3: 6.941, 4: 9.01218, 5: 10.811,
                    6: 12.011, 7: 14.00674, 8: 15.9994, 9: 18.998403, 10: 20.1797,
                    11: 22.989768, 12: 24.305, 13: 26.981539, 14: 28.0855,
                    15: 30.973762, 16: 32.066, 17: 35.4527, 18: 39.948, 19: 39.0983,
                    20: 40.078, 21: 44.95591, 22: 47.88, 23: 50.9415, 24: 51.9961,
                    25: 54.93805, 26: 55.847, 27: 58.9332, 28: 58.6934, 29: 63.546,
                    30: 65.39, 31: 69.723, 32: 72.61, 33: 74.92159, 34: 78.96,
                    35: 79.904, 36: 83.8, 37: 85.4678, 38: 87.62, 39: 88.90585,
                    40: 91.224, 41: 92.90638, 42: 95.94, 43: 97.9072, 44: 101.07,
                    45: 102.9055, 46: 106.42, 47: 107.8682, 48: 112.411, 49: 114.818,
                    50: 118.71, 51: 121.76, 52: 127.6, 53: 126.90447, 54: 131.29,
                    55: 132.90543, 56: 137.327, 57: 138.9055, 58: 140.115,
                    59: 140.90765, 60: 144.24, 61: 144.9127, 62: 150.36, 63: 151.965,
                    64: 157.25, 65: 158.92534, 66: 162.5, 67: 164.93032, 68: 167.26,
                    69: 168.93421, 70: 173.04, 71: 174.967, 72: 178.49, 73: 180.9479,
                    74: 183.84, 75: 186.207, 76: 190.23, 77: 192.22, 78: 195.08,
                    79: 196.96654, 80: 200.59, 81: 204.3833, 82: 207.2, 83: 208.98037,
                    84: 208.9824, 85: 209.9871, 86: 222.0176, 87: 223.0197,
                    88: 226.0254, 89: 227.0278, 90: 232.0381, 91: 231.03588,
                    92: 238.0289, 93: 237.048, 94: 244.0642, 95: 243.0614, 96: 247.0703,
                    97: 247.0703, 98: 251.0796, 99: 252.083, 100: 257.0951, 101: 258.1,
                    102: 259.1009, 103: 262.11}

    # === STANDARD CONSTRUCTOR
    def __init__(self, arg=None):
        """Initialize an Atom object from an atomic number or symbol.

        Args:
            arg(int/str): if int, interpreted as an atomic number. If str,
        interpreted as an element symbol. In None, it serves as a placeholder. 

        """
        if arg is None:
            self._Number = None
        elif isinstance(arg, int):
            self._Number = arg
        elif isinstance(arg, str):
            self._Number = Atom.SymbolToNumber[arg]
        else:
            raise ValueError('Cannot make an Atom from this argument: {}'.format(arg))

    # === PROPERTY EVALUATION METHODS
    def __eq__(self, other):
        """Atom equality operator."""
        return other is not None and self.Number == other.Number

    def __ne__(self, other):
        """Atom not-equal operator."""
        return not self == other

    def __lt__(self, other):
        """Atom less-than operator."""
        return self.Symbol < other.Symbol

    def __hash__(self):
        """Atom hash operator."""
        return self.Number.__hash__()

    @property
    def Symbol(self):
        """Element symbol."""
        return Atom.NumberToSymbol[self.Number]

    @property
    def Number(self):
        """Atomic number."""
        return self._Number

    @property
    def Mass(self):
        """Atomic mass (g/mol)."""
        return Atom.NumberToMass[self.Number]

    # === REPORTING METHODS
    def __repr__(self):
        """Atom representation."""
        return "Atom({})".format(self.Symbol)
