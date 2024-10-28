"""
Tools for working with chemical formulas

Objects
-------

- Formula

Exceptions
----------

- InvalidFormula

"""

import numpy as np
import os.path
from string import digits
from typing import Dict, Final, Tuple, Union
import utils


EM: Final[float] = 0.00054858  # electron mass


class Isotope:
    """
    Representation of an Isotope.

    Attributes
    ----------
    z: int
        Atomic number
    n: int
        Neutron number
    a: int
        Mass number
    m: float
        Exact mass.
    defect: float
        Difference between the exact mass and mass number.
    abundance: float
        Relative abundance of the isotope.

    """

    __slots__ = ("z", "n", "a", "m", "defect", "abundance")

    def __init__(self, z: int, a: int, m: float, abundance: float):
        self.z = z
        self.n = a - z
        self.a = a
        self.m = m
        self.defect = m - a
        self.abundance = abundance

    def __str__(self):
        return "{}{}".format(self.a, self.get_symbol())

    def __repr__(self):
        return "Isotope({})".format(str(self))

    def get_element(self) -> "Element":
        return PeriodicTable().get_element(self.z)

    def get_symbol(self) -> str:
        return self.get_element().symbol


class Element(object):
    """
    Representation of a chemical element.

    Attributes
    ----------
    name : str
        Element name.
    symbol : str
        Element symbol
    isotopes : Dict[int, Isotope]
        Mapping from mass number to an isotope
    z : int
        Atomic number.
    nominal_mass : int
        Mass number of the most abundant isotope

    """

    def __init__(self, symbol: str, name: str, isotopes: Dict[int, Isotope]):
        self.name = name
        self.symbol = symbol
        self.isotopes = isotopes
        monoisotope = self.get_monoisotope()
        self.z = monoisotope.z
        self.nominal_mass = monoisotope.a
        self.monoisotopic_mass = monoisotope.m
        self.mass_defect = self.monoisotopic_mass - self.nominal_mass

    def __repr__(self):
        return "Element({})".format(self.symbol)

    def __str__(self):  # pragma: no cover
        return self.symbol

    def get_abundances(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Returns the Mass number, exact mass and abundance of each Isotope.

        Returns
        -------
        m: array[int]
            Mass number of each isotope.
        M: array[float]
            Exact mass of each isotope.
        p: array[float]
            Abundance of each isotope.

        """
        isotopes = list(self.isotopes.values())
        m = np.array([x.a for x in isotopes], dtype=int)
        M = np.array([x.m for x in isotopes])
        p = np.array([x.abundance for x in isotopes])
        return m, M, p

    def get_mmi(self) -> Isotope:
        """
        Returns the isotope with the lowest atomic mass.

        """
        return min(self.isotopes.values(), key=lambda x: x.a)

    def get_monoisotope(self) -> Isotope:
        """
        Returns the most abundant isotope.

        """
        return max(self.isotopes.values(), key=lambda x: x.abundance)


def PeriodicTable():
    """
    Reference the PeriodicTable object.

    Examples
    --------
    >>> import tidyms as ms
    >>> ptable = ms.chem.PeriodicTable()

    """
    if _PeriodicTable.instance is None:
        _PeriodicTable.instance = _PeriodicTable()
    return _PeriodicTable.instance


class _PeriodicTable:
    """
    Periodic Table representation. Contains element and isotope information.

    Methods
    -------
    get_element
    get_isotope

    """

    instance = None

    def __init__(self):
        self._symbol_to_element = _make_periodic_table()
        self._z_to_element = {v.z: v for v in self._symbol_to_element.values()}
        self._za_to_isotope = dict()
        self._str_to_isotope = dict()
        for el_str in self._symbol_to_element:
            el = self._symbol_to_element[el_str]
            for isotope in el.isotopes.values():
                self._za_to_isotope[(isotope.z, isotope.a)] = isotope
                self._str_to_isotope[str(isotope.a) + el_str] = isotope

    def get_element(self, element: Union[str, int]) -> Element:
        """
        Returns an Element object using its symbol or atomic number.

        Parameters
        ----------
        element : str or int
            element symbol or atomic number.

        Returns
        -------
        Element

        Examples
        --------
        >>> import tidyms as ms
        >>> ptable = ms.chem.PeriodicTable()
        >>> h = ptable.get_element("H")
        >>> c = ptable.get_element(6)

        """
        if isinstance(element, int):
            element = self._z_to_element[element]
        else:
            element = self._symbol_to_element[element]
        return element

    def __iter__(self):
        for el in self._symbol_to_element.values():
            yield el

    def get_isotope(self, x: str, copy: bool = False) -> Isotope:
        """
        Returns an isotope object from a string representation.

        Parameters
        ----------
        x : str
            A string representation of an isotope. If only the symbol is
            provided in the string, the monoisotope is returned.
        copy : bool
            If True creates a new Isotope object.

        Returns
        -------
        Isotope

        Examples
        --------
        >>> import tidyms as ms
        >>> ptable = ms.chem.PeriodicTable()
        >>> d = ptable.get_isotope("2H")
        >>> cl35 = ptable.get_isotope("Cl")

        """
        try:
            if x[0] in digits:
                isotope = self._str_to_isotope[x]
            else:
                isotope = self.get_element(x).get_monoisotope()
            if copy:
                isotope = Isotope(isotope.z, isotope.a, isotope.m, isotope.abundance)
            return isotope
        except KeyError:
            msg = "{} is not a valid input.".format(x)
            raise InvalidIsotope(msg)


def _make_periodic_table() -> Dict[str, Element]:
    this_dir, _ = os.path.split(__file__)

    element_data = {'Xx': 'Dummy', 'H': 'Hydrogen', 'He': 'Helium', 'Li': 'Lithium', 'Be': 'Beryllium', 'B': 'Boron', 'C': 'Carbon', 'N': 'Nitrogen', 'O': 'Oxygen', 'F': 'Fluorine', 'Ne': 'Neon', 'Na': 'Sodium', 'Mg': 'Magnesium', 'Al': 'Aluminium', 'Si': 'Silicon', 'P': 'Phosphorus', 'S': 'Sulfur', 'Cl': 'Chlorine', 'Ar': 'Argon', 'K': 'Potassium', 'Ca': 'Calcium', 'Sc': 'Scandium', 'Ti': 'Titanium', 'V': 'Vanadium', 'Cr': 'Chromium', 'Mn': 'Manganese', 'Fe': 'Iron', 'Co': 'Cobalt', 'Ni': 'Nickel', 'Cu': 'Copper', 'Zn': 'Zinc', 'Ga': 'Gallium', 'Ge': 'Germanium', 'As': 'Arsenic', 'Se': 'Selenium', 'Br': 'Bromine', 'Kr': 'Krypton', 'Rb': 'Rubidium', 'Sr': 'Strontium', 'Y': 'Yttrium', 'Zr': 'Zirconium', 'Nb': 'Niobium', 'Mo': 'Molybdenum', 'Tc': 'Technetium', 'Ru': 'Ruthenium', 'Rh': 'Rhodium', 'Pd': 'Palladium', 'Ag': 'Silver', 'Cd': 'Cadmium', 'In': 'Indium', 'Sn': 'Tin', 'Sb': 'Antimony', 'Te': 'Tellurium', 'I': 'Iodine', 'Xe': 'Xenon', 'Cs': 'Caesium', 'Ba': 'Barium', 'La': 'Lanthanum', 'Ce': 'Cerium', 'Pr': 'Praseodymium', 'Nd': 'Neodymium', 'Pm': 'Promethium', 'Sm': 'Samarium', 'Eu': 'Europium', 'Gd': 'Gadolinium', 'Tb': 'Terbium', 'Dy': 'Dysprosium', 'Ho': 'Holmium', 'Er': 'Erbium', 'Tm': 'Thulium', 'Yb': 'Ytterbium', 'Lu': 'Lutetium', 'Hf': 'Hafnium', 'Ta': 'Tantalum', 'W': 'Tungsten', 'Re': 'Rhenium', 'Os': 'Osmium', 'Ir': 'Iridium', 'Pt': 'Platinum', 'Au': 'Gold', 'Hg': 'Mercury', 'Tl': 'Thallium', 'Pb': 'Lead', 'Bi': 'Bismuth', 'Po': 'Polonium', 'At': 'Astatine', 'Rn': 'Radon', 'Fr': 'Francium', 'Ra': 'Radium', 'Ac': 'Actinium', 'Th': 'Thorium', 'Pa': 'Protactinium', 'U': 'Uranium', 'Np': 'Neptunium', 'Pu': 'Plutonium', 'Am': 'Americium', 'Cm': 'Curium', 'Bk': 'Berkelium', 'Cf': 'Californium', 'Es': 'Einsteinium', 'Fm': 'Fermium', 'Md': 'Mendelevium', 'No': 'Nobelium', 'Lr': 'Lawrencium', 'Rf': 'Rutherfordium', 'Db': 'Dubnium', 'Sg': 'Seaborgium', 'Bh': 'Bohrium', 'Hs': 'Hassium', 'Mt': 'Meitnerium', 'Ds': 'Darmstadtium', 'Rg': 'Roentgenium', 'Cn': 'Copernicium', 'Uut': 'Ununtrium', 'Fl': 'Flerovium', 'Uup': 'Ununpentium', 'Lv': 'Livermorium', 'Uus': 'Ununseptium', 'Uuo': 'Ununoctium'}


    isotope_data = {'H': [{'a': 1, 'abundance': 0.999885, 'm': 1.007825032, 'z': 1}, {'a': 2, 'abundance': 0.000115, 'm': 2.014101778, 'z': 1}], 'He': [{'a': 3, 'abundance': 1.37e-06, 'm': 3.016029319, 'z': 2}, {'a': 4, 'abundance': 0.99999863, 'm': 4.002603254, 'z': 2}], 'Li': [{'a': 6, 'abundance': 0.0759, 'm': 6.015122795, 'z': 3}, {'a': 7, 'abundance': 0.9240999999999999, 'm': 7.01600455, 'z': 3}], 'Be': [{'a': 9, 'abundance': 1.0, 'm': 9.0121822, 'z': 4}], 'B': [{'a': 10, 'abundance': 0.19899999999999998, 'm': 10.012937, 'z': 5}, {'a': 11, 'abundance': 0.8009999999999999, 'm': 11.0093054, 'z': 5}], 'C': [{'a': 12, 'abundance': 0.9893000000000001, 'm': 12.0, 'z': 6}, {'a': 13, 'abundance': 0.010700000000000001, 'm': 13.00335484, 'z': 6}], 'N': [{'a': 14, 'abundance': 0.9963200000000001, 'm': 14.003074, 'z': 7}, {'a': 15, 'abundance': 0.00368, 'm': 15.0001089, 'z': 7}], 'O': [{'a': 16, 'abundance': 0.9975700000000001, 'm': 15.99491462, 'z': 8}, {'a': 17, 'abundance': 0.00037999999999999997, 'm': 16.9991317, 'z': 8}, {'a': 18, 'abundance': 0.0020499999999999997, 'm': 17.999161, 'z': 8}], 'F': [{'a': 19, 'abundance': 1.0, 'm': 18.99840322, 'z': 9}], 'Ne': [{'a': 20, 'abundance': 0.9048, 'm': 19.99244018, 'z': 10}, {'a': 21, 'abundance': 0.0027, 'm': 20.99384668, 'z': 10}, {'a': 22, 'abundance': 0.0925, 'm': 21.99138511, 'z': 10}], 'Na': [{'a': 23, 'abundance': 1.0, 'm': 22.98976928, 'z': 11}], 'Mg': [{'a': 24, 'abundance': 0.7898999999999999, 'm': 23.9850417, 'z': 12}, {'a': 25, 'abundance': 0.1, 'm': 24.98583692, 'z': 12}, {'a': 26, 'abundance': 0.1101, 'm': 25.98259293, 'z': 12}], 'Al': [{'a': 27, 'abundance': 1.0, 'm': 26.98153863, 'z': 13}], 'Si': [{'a': 28, 'abundance': 0.9222969999999999, 'm': 27.97692653, 'z': 14}, {'a': 29, 'abundance': 0.046832000000000006, 'm': 28.9764947, 'z': 14}, {'a': 30, 'abundance': 0.030872, 'm': 29.97377017, 'z': 14}], 'P': [{'a': 31, 'abundance': 1.0, 'm': 30.97376163, 'z': 15}], 'S': [{'a': 32, 'abundance': 0.9493, 'm': 31.972071, 'z': 16}, {'a': 33, 'abundance': 0.0076, 'm': 32.97145876, 'z': 16}, {'a': 34, 'abundance': 0.0429, 'm': 33.9678669, 'z': 16}, {'a': 36, 'abundance': 0.0002, 'm': 35.96708076, 'z': 16}], 'Cl': [{'a': 35, 'abundance': 0.7578, 'm': 34.96885268, 'z': 17}, {'a': 37, 'abundance': 0.2422, 'm': 36.96590259, 'z': 17}], 'Ar': [{'a': 36, 'abundance': 0.0033650000000000004, 'm': 35.96754511, 'z': 18}, {'a': 38, 'abundance': 0.0006320000000000001, 'm': 37.9627324, 'z': 18}, {'a': 40, 'abundance': 0.9960030000000001, 'm': 39.96238312, 'z': 18}], 'K': [{'a': 39, 'abundance': 0.932581, 'm': 38.96370668, 'z': 19}, {'a': 40, 'abundance': 0.000117, 'm': 39.96399848, 'z': 19}, {'a': 41, 'abundance': 0.067302, 'm': 40.96182576, 'z': 19}], 'Ca': [{'a': 40, 'abundance': 0.96941, 'm': 39.96259098, 'z': 20}, {'a': 42, 'abundance': 0.00647, 'm': 41.95861801, 'z': 20}, {'a': 43, 'abundance': 0.00135, 'm': 42.9587666, 'z': 20}, {'a': 44, 'abundance': 0.02086, 'm': 43.9554818, 'z': 20}, {'a': 46, 'abundance': 4e-05, 'm': 45.9536926, 'z': 20}, {'a': 48, 'abundance': 0.00187, 'm': 47.952534, 'z': 20}], 'Sc': [{'a': 45, 'abundance': 1.0, 'm': 44.9559119, 'z': 21}], 'Ti': [{'a': 46, 'abundance': 0.0825, 'm': 45.9526316, 'z': 22}, {'a': 47, 'abundance': 0.07440000000000001, 'm': 46.9517631, 'z': 22}, {'a': 48, 'abundance': 0.7372, 'm': 47.9479463, 'z': 22}, {'a': 49, 'abundance': 0.0541, 'm': 48.94787, 'z': 22}, {'a': 50, 'abundance': 0.0518, 'm': 49.9447912, 'z': 22}], 'V': [{'a': 50, 'abundance': 0.0025, 'm': 49.9471585, 'z': 23}, {'a': 51, 'abundance': 0.9975, 'm': 50.9439595, 'z': 23}], 'Cr': [{'a': 50, 'abundance': 0.043449999999999996, 'm': 49.9460442, 'z': 24}, {'a': 52, 'abundance': 0.83789, 'm': 51.9405075, 'z': 24}, {'a': 53, 'abundance': 0.09501, 'm': 52.9406494, 'z': 24}, {'a': 54, 'abundance': 0.02365, 'm': 53.9388804, 'z': 24}], 'Mn': [{'a': 55, 'abundance': 1.0, 'm': 54.9380451, 'z': 25}], 'Fe': [{'a': 54, 'abundance': 0.058449999999999995, 'm': 53.9396105, 'z': 26}, {'a': 56, 'abundance': 0.91754, 'm': 55.9349375, 'z': 26}, {'a': 57, 'abundance': 0.02119, 'm': 56.935394, 'z': 26}, {'a': 58, 'abundance': 0.0028199999999999996, 'm': 57.9332756, 'z': 26}], 'Co': [{'a': 59, 'abundance': 1.0, 'm': 58.933195, 'z': 27}], 'Ni': [{'a': 58, 'abundance': 0.680769, 'm': 57.9353429, 'z': 28}, {'a': 60, 'abundance': 0.262231, 'm': 59.9307864, 'z': 28}, {'a': 61, 'abundance': 0.011399, 'm': 60.931056, 'z': 28}, {'a': 62, 'abundance': 0.036345, 'm': 61.9283451, 'z': 28}, {'a': 64, 'abundance': 0.009256, 'm': 63.927966, 'z': 28}], 'Cu': [{'a': 63, 'abundance': 0.6917, 'm': 62.9295975, 'z': 29}, {'a': 65, 'abundance': 0.30829999999999996, 'm': 64.9277895, 'z': 29}], 'Zn': [{'a': 64, 'abundance': 0.4863, 'm': 63.9291422, 'z': 30}, {'a': 66, 'abundance': 0.27899999999999997, 'm': 65.9260334, 'z': 30}, {'a': 67, 'abundance': 0.040999999999999995, 'm': 66.9271273, 'z': 30}, {'a': 68, 'abundance': 0.1875, 'm': 67.9248442, 'z': 30}, {'a': 70, 'abundance': 0.0062, 'm': 69.9253193, 'z': 30}], 'Ga': [{'a': 69, 'abundance': 0.60108, 'm': 68.9255736, 'z': 31}, {'a': 71, 'abundance': 0.39892000000000005, 'm': 70.9247013, 'z': 31}], 'Ge': [{'a': 70, 'abundance': 0.2084, 'm': 69.9242474, 'z': 32}, {'a': 72, 'abundance': 0.2754, 'm': 71.9220758, 'z': 32}, {'a': 73, 'abundance': 0.07730000000000001, 'm': 72.9234589, 'z': 32}, {'a': 74, 'abundance': 0.3628, 'm': 73.9211778, 'z': 32}, {'a': 76, 'abundance': 0.0761, 'm': 75.9214026, 'z': 32}], 'As': [{'a': 75, 'abundance': 1.0, 'm': 74.9215965, 'z': 33}], 'Se': [{'a': 74, 'abundance': 0.0089, 'm': 73.9224764, 'z': 34}, {'a': 76, 'abundance': 0.09369999999999999, 'm': 75.9192136, 'z': 34}, {'a': 77, 'abundance': 0.07629999999999999, 'm': 76.919914, 'z': 34}, {'a': 78, 'abundance': 0.2377, 'm': 77.9173091, 'z': 34}, {'a': 80, 'abundance': 0.4961, 'm': 79.9165213, 'z': 34}, {'a': 82, 'abundance': 0.0873, 'm': 81.9166994, 'z': 34}], 'Br': [{'a': 79, 'abundance': 0.5069, 'm': 78.9183371, 'z': 35}, {'a': 81, 'abundance': 0.49310000000000004, 'm': 80.9162906, 'z': 35}], 'Kr': [{'a': 78, 'abundance': 0.0034999999999999996, 'm': 77.9203648, 'z': 36}, {'a': 80, 'abundance': 0.022799999999999997, 'm': 79.916379, 'z': 36}, {'a': 82, 'abundance': 0.1158, 'm': 81.9134836, 'z': 36}, {'a': 83, 'abundance': 0.1149, 'm': 82.914136, 'z': 36}, {'a': 84, 'abundance': 0.57, 'm': 83.911507, 'z': 36}, {'a': 86, 'abundance': 0.17300000000000001, 'm': 85.91061073, 'z': 36}], 'Rb': [{'a': 85, 'abundance': 0.7217, 'm': 84.91178974, 'z': 37}, {'a': 87, 'abundance': 0.2783, 'm': 86.90918053, 'z': 37}], 'Sr': [{'a': 84, 'abundance': 0.005600000000000001, 'm': 83.913425, 'z': 38}, {'a': 86, 'abundance': 0.0986, 'm': 85.9092602, 'z': 38}, {'a': 87, 'abundance': 0.07, 'm': 86.9088771, 'z': 38}, {'a': 88, 'abundance': 0.8258, 'm': 87.9056121, 'z': 38}], 'Y': [{'a': 89, 'abundance': 1.0, 'm': 88.9058483, 'z': 39}], 'Zr': [{'a': 90, 'abundance': 0.5145000000000001, 'm': 89.9047044, 'z': 40}, {'a': 91, 'abundance': 0.11220000000000001, 'm': 90.9056458, 'z': 40}, {'a': 92, 'abundance': 0.17149999999999999, 'm': 91.9050408, 'z': 40}, {'a': 94, 'abundance': 0.17379999999999998, 'm': 93.9063152, 'z': 40}, {'a': 96, 'abundance': 0.027999999999999997, 'm': 95.9082734, 'z': 40}], 'Nb': [{'a': 93, 'abundance': 1.0, 'm': 92.9063781, 'z': 41}], 'Mo': [{'a': 92, 'abundance': 0.1484, 'm': 91.906811, 'z': 42}, {'a': 94, 'abundance': 0.0925, 'm': 93.9050883, 'z': 42}, {'a': 95, 'abundance': 0.1592, 'm': 94.9058421, 'z': 42}, {'a': 96, 'abundance': 0.1668, 'm': 95.9046795, 'z': 42}, {'a': 97, 'abundance': 0.0955, 'm': 96.9060215, 'z': 42}, {'a': 98, 'abundance': 0.2413, 'm': 97.9054082, 'z': 42}, {'a': 100, 'abundance': 0.09630000000000001, 'm': 99.907477, 'z': 42}], 'Ru': [{'a': 96, 'abundance': 0.0554, 'm': 95.907598, 'z': 44}, {'a': 98, 'abundance': 0.0187, 'm': 97.905287, 'z': 44}, {'a': 99, 'abundance': 0.1276, 'm': 98.9059393, 'z': 44}, {'a': 100, 'abundance': 0.126, 'm': 99.9042195, 'z': 44}, {'a': 101, 'abundance': 0.17059999999999997, 'm': 100.9055821, 'z': 44}, {'a': 102, 'abundance': 0.3155, 'm': 101.9043493, 'z': 44}, {'a': 104, 'abundance': 0.1862, 'm': 103.905433, 'z': 44}], 'Rh': [{'a': 103, 'abundance': 1.0, 'm': 102.905504, 'z': 45}], 'Pd': [{'a': 102, 'abundance': 0.0102, 'm': 101.905609, 'z': 46}, {'a': 104, 'abundance': 0.1114, 'm': 103.904036, 'z': 46}, {'a': 105, 'abundance': 0.22329999999999997, 'm': 104.905085, 'z': 46}, {'a': 106, 'abundance': 0.2733, 'm': 105.903486, 'z': 46}, {'a': 108, 'abundance': 0.2646, 'm': 107.903892, 'z': 46}, {'a': 110, 'abundance': 0.11720000000000001, 'm': 109.905153, 'z': 46}], 'Ag': [{'a': 107, 'abundance': 0.51839, 'm': 106.905097, 'z': 47}, {'a': 109, 'abundance': 0.48161000000000004, 'm': 108.904752, 'z': 47}], 'Cd': [{'a': 106, 'abundance': 0.0125, 'm': 105.906459, 'z': 48}, {'a': 108, 'abundance': 0.0089, 'm': 107.904184, 'z': 48}, {'a': 110, 'abundance': 0.1249, 'm': 109.9030021, 'z': 48}, {'a': 111, 'abundance': 0.128, 'm': 110.9041781, 'z': 48}, {'a': 112, 'abundance': 0.2413, 'm': 111.9027578, 'z': 48}, {'a': 113, 'abundance': 0.1222, 'm': 112.9044017, 'z': 48}, {'a': 114, 'abundance': 0.2873, 'm': 113.9033585, 'z': 48}, {'a': 116, 'abundance': 0.07490000000000001, 'm': 115.904756, 'z': 48}], 'In': [{'a': 113, 'abundance': 0.0429, 'm': 112.904058, 'z': 49}, {'a': 115, 'abundance': 0.9571, 'm': 114.903878, 'z': 49}], 'Sn': [{'a': 112, 'abundance': 0.0097, 'm': 111.904818, 'z': 50}, {'a': 114, 'abundance': 0.0066, 'm': 113.902779, 'z': 50}, {'a': 115, 'abundance': 0.0034000000000000002, 'm': 114.903342, 'z': 50}, {'a': 116, 'abundance': 0.1454, 'm': 115.901741, 'z': 50}, {'a': 117, 'abundance': 0.0768, 'm': 116.902952, 'z': 50}, {'a': 118, 'abundance': 0.2422, 'm': 117.901603, 'z': 50}, {'a': 119, 'abundance': 0.0859, 'm': 118.903308, 'z': 50}, {'a': 120, 'abundance': 0.3258, 'm': 119.9021947, 'z': 50}, {'a': 122, 'abundance': 0.0463, 'm': 121.903439, 'z': 50}, {'a': 124, 'abundance': 0.0579, 'm': 123.9052739, 'z': 50}], 'Sb': [{'a': 121, 'abundance': 0.5721, 'm': 120.9038157, 'z': 51}, {'a': 123, 'abundance': 0.4279, 'm': 122.904214, 'z': 51}], 'Te': [{'a': 120, 'abundance': 0.0009, 'm': 119.90402, 'z': 52}, {'a': 122, 'abundance': 0.0255, 'm': 121.9030439, 'z': 52}, {'a': 123, 'abundance': 0.0089, 'm': 122.90427, 'z': 52}, {'a': 124, 'abundance': 0.047400000000000005, 'm': 123.9028179, 'z': 52}, {'a': 125, 'abundance': 0.0707, 'm': 124.9044307, 'z': 52}, {'a': 126, 'abundance': 0.1884, 'm': 125.9033117, 'z': 52}, {'a': 128, 'abundance': 0.31739999999999996, 'm': 127.9044631, 'z': 52}, {'a': 130, 'abundance': 0.3408, 'm': 129.9062244, 'z': 52}], 'I': [{'a': 127, 'abundance': 1.0, 'm': 126.904473, 'z': 53}], 'Xe': [{'a': 124, 'abundance': 0.0009, 'm': 123.905893, 'z': 54}, {'a': 126, 'abundance': 0.0009, 'm': 125.904274, 'z': 54}, {'a': 128, 'abundance': 0.0192, 'm': 127.9035313, 'z': 54}, {'a': 129, 'abundance': 0.2644, 'm': 128.9047794, 'z': 54}, {'a': 130, 'abundance': 0.0408, 'm': 129.903508, 'z': 54}, {'a': 131, 'abundance': 0.2118, 'm': 130.9050824, 'z': 54}, {'a': 132, 'abundance': 0.26890000000000003, 'm': 131.9041535, 'z': 54}, {'a': 134, 'abundance': 0.10439999999999999, 'm': 133.9053945, 'z': 54}, {'a': 136, 'abundance': 0.08869999999999999, 'm': 135.907219, 'z': 54}], 'Cs': [{'a': 133, 'abundance': 1.0, 'm': 132.9054519, 'z': 55}], 'Ba': [{'a': 130, 'abundance': 0.00106, 'm': 129.9063208, 'z': 56}, {'a': 132, 'abundance': 0.00101, 'm': 131.9050613, 'z': 56}, {'a': 134, 'abundance': 0.024169999999999997, 'm': 133.9045084, 'z': 56}, {'a': 135, 'abundance': 0.06591999999999999, 'm': 134.9056886, 'z': 56}, {'a': 136, 'abundance': 0.07854, 'm': 135.9045759, 'z': 56}, {'a': 137, 'abundance': 0.11231999999999999, 'm': 136.9058274, 'z': 56}, {'a': 138, 'abundance': 0.71698, 'm': 137.9052472, 'z': 56}], 'La': [{'a': 138, 'abundance': 0.0009, 'm': 137.907112, 'z': 57}, {'a': 139, 'abundance': 0.9991, 'm': 138.9063533, 'z': 57}], 'Ce': [{'a': 136, 'abundance': 0.00185, 'm': 135.907172, 'z': 58}, {'a': 138, 'abundance': 0.00251, 'm': 137.905991, 'z': 58}, {'a': 140, 'abundance': 0.8845000000000001, 'm': 139.9054387, 'z': 58}, {'a': 142, 'abundance': 0.11114, 'm': 141.909244, 'z': 58}], 'Pr': [{'a': 141, 'abundance': 1.0, 'm': 140.9076528, 'z': 59}], 'Nd': [{'a': 142, 'abundance': 0.272, 'm': 141.9077233, 'z': 60}, {'a': 143, 'abundance': 0.122, 'm': 142.9098143, 'z': 60}, {'a': 144, 'abundance': 0.23800000000000002, 'm': 143.9100873, 'z': 60}, {'a': 145, 'abundance': 0.083, 'm': 144.9125736, 'z': 60}, {'a': 146, 'abundance': 0.172, 'm': 145.9131169, 'z': 60}, {'a': 148, 'abundance': 0.057, 'm': 147.916893, 'z': 60}, {'a': 150, 'abundance': 0.055999999999999994, 'm': 149.920891, 'z': 60}], 'Sm': [{'a': 144, 'abundance': 0.030699999999999998, 'm': 143.911999, 'z': 62}, {'a': 147, 'abundance': 0.1499, 'm': 146.9148979, 'z': 62}, {'a': 148, 'abundance': 0.1124, 'm': 147.9148227, 'z': 62}, {'a': 149, 'abundance': 0.1382, 'm': 148.9171847, 'z': 62}, {'a': 150, 'abundance': 0.0738, 'm': 149.9172755, 'z': 62}, {'a': 152, 'abundance': 0.2675, 'm': 151.9197324, 'z': 62}, {'a': 154, 'abundance': 0.2275, 'm': 153.9222093, 'z': 62}], 'Eu': [{'a': 151, 'abundance': 0.4781, 'm': 150.9198502, 'z': 63}, {'a': 153, 'abundance': 0.5219, 'm': 152.9212303, 'z': 63}], 'Gd': [{'a': 152, 'abundance': 0.002, 'm': 151.919791, 'z': 64}, {'a': 154, 'abundance': 0.0218, 'm': 153.9208656, 'z': 64}, {'a': 155, 'abundance': 0.14800000000000002, 'm': 154.922622, 'z': 64}, {'a': 156, 'abundance': 0.2047, 'm': 155.9221227, 'z': 64}, {'a': 157, 'abundance': 0.1565, 'm': 156.9239601, 'z': 64}, {'a': 158, 'abundance': 0.2484, 'm': 157.9241039, 'z': 64}, {'a': 160, 'abundance': 0.2186, 'm': 159.9270541, 'z': 64}], 'Tb': [{'a': 159, 'abundance': 1.0, 'm': 158.9253468, 'z': 65}], 'Dy': [{'a': 156, 'abundance': 0.0006, 'm': 155.924283, 'z': 66}, {'a': 158, 'abundance': 0.001, 'm': 157.924409, 'z': 66}, {'a': 160, 'abundance': 0.023399999999999997, 'm': 159.9251975, 'z': 66}, {'a': 161, 'abundance': 0.1891, 'm': 160.9269334, 'z': 66}, {'a': 162, 'abundance': 0.2551, 'm': 161.9267984, 'z': 66}, {'a': 163, 'abundance': 0.249, 'm': 162.9287312, 'z': 66}, {'a': 164, 'abundance': 0.2818, 'm': 163.9291748, 'z': 66}], 'Ho': [{'a': 165, 'abundance': 1.0, 'm': 164.9303221, 'z': 67}], 'Er': [{'a': 162, 'abundance': 0.0014000000000000002, 'm': 161.928778, 'z': 68}, {'a': 164, 'abundance': 0.0161, 'm': 163.9292, 'z': 68}, {'a': 166, 'abundance': 0.3361, 'm': 165.9302931, 'z': 68}, {'a': 167, 'abundance': 0.2293, 'm': 166.9320482, 'z': 68}, {'a': 168, 'abundance': 0.26780000000000004, 'm': 167.9323702, 'z': 68}, {'a': 170, 'abundance': 0.1493, 'm': 169.9354643, 'z': 68}], 'Tm': [{'a': 169, 'abundance': 1.0, 'm': 168.9342133, 'z': 69}], 'Yb': [{'a': 168, 'abundance': 0.0013, 'm': 167.933897, 'z': 70}, {'a': 170, 'abundance': 0.0304, 'm': 169.9347618, 'z': 70}, {'a': 171, 'abundance': 0.14279999999999998, 'm': 170.9363258, 'z': 70}, {'a': 172, 'abundance': 0.2183, 'm': 171.9363815, 'z': 70}, {'a': 173, 'abundance': 0.1613, 'm': 172.9382108, 'z': 70}, {'a': 174, 'abundance': 0.31829999999999997, 'm': 173.9388621, 'z': 70}, {'a': 176, 'abundance': 0.1276, 'm': 175.9425717, 'z': 70}], 'Lu': [{'a': 175, 'abundance': 0.9741, 'm': 174.9407718, 'z': 71}, {'a': 176, 'abundance': 0.0259, 'm': 175.9426863, 'z': 71}], 'Hf': [{'a': 174, 'abundance': 0.0016, 'm': 173.940046, 'z': 72}, {'a': 176, 'abundance': 0.0526, 'm': 175.9414086, 'z': 72}, {'a': 177, 'abundance': 0.18600000000000003, 'm': 176.9432207, 'z': 72}, {'a': 178, 'abundance': 0.2728, 'm': 177.9436988, 'z': 72}, {'a': 179, 'abundance': 0.1362, 'm': 178.9458161, 'z': 72}, {'a': 180, 'abundance': 0.3508, 'm': 179.94655, 'z': 72}], 'Ta': [{'a': 180, 'abundance': 0.00012, 'm': 179.9474648, 'z': 73}, {'a': 181, 'abundance': 0.99988, 'm': 180.9479958, 'z': 73}], 'W': [{'a': 180, 'abundance': 0.0012, 'm': 179.946704, 'z': 74}, {'a': 182, 'abundance': 0.265, 'm': 181.9482042, 'z': 74}, {'a': 183, 'abundance': 0.1431, 'm': 182.950223, 'z': 74}, {'a': 184, 'abundance': 0.3064, 'm': 183.9509312, 'z': 74}, {'a': 186, 'abundance': 0.2843, 'm': 185.9543641, 'z': 74}], 'Re': [{'a': 185, 'abundance': 0.374, 'm': 184.952955, 'z': 75}, {'a': 187, 'abundance': 0.626, 'm': 186.9557531, 'z': 75}], 'Os': [{'a': 184, 'abundance': 0.0002, 'm': 183.9524891, 'z': 76}, {'a': 186, 'abundance': 0.0159, 'm': 185.9538382, 'z': 76}, {'a': 187, 'abundance': 0.0196, 'm': 186.9557505, 'z': 76}, {'a': 188, 'abundance': 0.1324, 'm': 187.9558382, 'z': 76}, {'a': 189, 'abundance': 0.16149999999999998, 'm': 188.9581475, 'z': 76}, {'a': 190, 'abundance': 0.2626, 'm': 189.958447, 'z': 76}, {'a': 192, 'abundance': 0.4078, 'm': 191.9614807, 'z': 76}], 'Ir': [{'a': 191, 'abundance': 0.373, 'm': 190.960594, 'z': 77}, {'a': 193, 'abundance': 0.627, 'm': 192.9629264, 'z': 77}], 'Pt': [{'a': 190, 'abundance': 0.00014000000000000001, 'm': 189.959932, 'z': 78}, {'a': 192, 'abundance': 0.00782, 'm': 191.961038, 'z': 78}, {'a': 194, 'abundance': 0.32966999999999996, 'm': 193.9626803, 'z': 78}, {'a': 195, 'abundance': 0.33832, 'm': 194.9647911, 'z': 78}, {'a': 196, 'abundance': 0.25242000000000003, 'm': 195.9649515, 'z': 78}, {'a': 198, 'abundance': 0.07163, 'm': 197.967893, 'z': 78}], 'Au': [{'a': 197, 'abundance': 1.0, 'm': 196.9665687, 'z': 79}], 'Hg': [{'a': 196, 'abundance': 0.0015, 'm': 195.965833, 'z': 80}, {'a': 198, 'abundance': 0.09970000000000001, 'm': 197.966769, 'z': 80}, {'a': 199, 'abundance': 0.16870000000000002, 'm': 198.9682799, 'z': 80}, {'a': 200, 'abundance': 0.231, 'm': 199.968326, 'z': 80}, {'a': 201, 'abundance': 0.1318, 'm': 200.9703023, 'z': 80}, {'a': 202, 'abundance': 0.2986, 'm': 201.970643, 'z': 80}, {'a': 204, 'abundance': 0.0687, 'm': 203.9734939, 'z': 80}], 'Tl': [{'a': 203, 'abundance': 0.29524, 'm': 202.9723442, 'z': 81}, {'a': 205, 'abundance': 0.7047599999999999, 'm': 204.9744275, 'z': 81}], 'Pb': [{'a': 204, 'abundance': 0.013999999999999999, 'm': 203.9730436, 'z': 82}, {'a': 206, 'abundance': 0.24100000000000002, 'm': 205.9744653, 'z': 82}, {'a': 207, 'abundance': 0.221, 'm': 206.9758969, 'z': 82}, {'a': 208, 'abundance': 0.524, 'm': 207.9766521, 'z': 82}], 'Bi': [{'a': 209, 'abundance': 1.0, 'm': 208.9803987, 'z': 83}], 'Th': [{'a': 232, 'abundance': 1.0, 'm': 232.0380553, 'z': 90}], 'Pa': [{'a': 231, 'abundance': 1.0, 'm': 231.035884, 'z': 91}], 'U': [{'a': 234, 'abundance': 5.4999999999999995e-05, 'm': 234.0409521, 'z': 92}, {'a': 235, 'abundance': 0.0072, 'm': 235.0439299, 'z': 92}, {'a': 238, 'abundance': 0.992745, 'm': 238.0507882, 'z': 92}]}

    periodic_table = dict()
    for element in isotope_data:
        element_isotopes = isotope_data[element]
        isotopes = {x["a"]: Isotope(**x) for x in element_isotopes}
        name = element_data[element]
        periodic_table[element] = Element(element, name, isotopes)
    return periodic_table


class InvalidIsotope(ValueError):
    pass


def cartesian_product(*args):
    res = None
    for x in args:
        if res is None:
            # initialize cartesian product array
            res = np.array(x)
            res = res.reshape((res.size, 1))
        else:
            x = np.array(x)
            row, col = res.shape
            new_res_shape = (row * x.size, col + 1)
            new_res = np.zeros(shape=new_res_shape, dtype=res.dtype)
            ind = np.repeat(np.arange(row), x.size)
            new_col = np.tile(x, row)
            new_res[:, :col] = res[ind]
            new_res[:, -1] = new_col
            res = new_res
    return res



import numpy as np
from functools import lru_cache
from scipy.stats import multinomial
from typing import Dict, Optional, Tuple



def make_envelope_arrays(
    isotope: Isotope, n_min: int, n_max: int, max_length: int, p=None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Creates an array of exact mass and abundance for homonuclear formulas.

    Parameters
    ----------
    isotope : Isotope
    n_min : int
        Minimum formula coefficient
    n_max : int
        Maximum formula coefficient
    max_length : int
        Length of the envelope
    p : array or None, default=None
        Element abundance. If None, the natural abundance is used.

    Returns
    -------
    M : (n_max - n_min + 1, max_length) array
        Coefficients exact mass.
    p : (n_max - n_min + 1, max_length) array
        Coefficients abundance.


    """
    rows = n_max - n_min + 1
    M_arr = np.zeros((rows, max_length))
    p_arr = np.zeros((rows, max_length))
    for k in range(n_min, n_max + 1):
        Mk, pk = _get_n_atoms_envelope(isotope, k, max_length, p=p)
        M_arr[k - n_min] = Mk
        p_arr[k - n_min] = pk
    return M_arr, p_arr


def find_formula_envelope(
    composition: Dict[Isotope, int],
    max_length: int,
    p: Optional[Dict[str, np.ndarray]] = None,
    min_p: float = 1e-10,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes the isotopic envelope for a formula.

    """
    if p is None:
        p = dict()

    # initialize an empty envelope for the formula
    Mf = np.zeros((1, max_length), dtype=float)
    pf = np.zeros((1, max_length), dtype=float)
    pf[0, 0] = 1

    for i, coeff in composition.items():
        i_p = p.get(i.get_symbol())
        Mi, pi = _get_n_atoms_envelope(i, coeff, max_length, p=i_p)
        Mi = Mi.reshape((1, Mi.size))
        pi = pi.reshape((1, pi.size))
        Mf, pf = combine_envelopes(Mf, pf, Mi, pi)
    valid_p_mask = pf >= min_p
    pf = pf[valid_p_mask].flatten()
    Mf = Mf[valid_p_mask].flatten()
    return Mf, pf


def combine_envelopes(
    M1: np.ndarray,
    p1: np.ndarray,
    M2: np.ndarray,
    p2: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Combines exact mass and abundance of two envelopes.

    All arrays must be 2-dimensional and have the same shape.

    """
    shape = M1.shape
    M = np.zeros(shape, dtype=float)
    p = np.zeros(shape, dtype=float)
    # Ignore zero division errors when normalizing by pk
    with np.errstate(divide='ignore', invalid='ignore'):
        for k in range(shape[1]):
            pk = (p1[:, : k + 1] * p2[:, k::-1]).sum(axis=1)
            k1 = k + 1
            k2 = k
            Mk = (p1[:, :k1] * M1[:, :k1] * p2[:, k2::-1]) + (
                p1[:, :k1] * M2[:, k2::-1] * p2[:, k2::-1]
            )
            M[:, k] = Mk.sum(axis=1) / pk
            p[:, k] = pk
    np.nan_to_num(M, copy=False)
    return M, p


def _get_n_atoms_envelope(
    isotope: Isotope, n: int, max_length: int, p: Optional[np.ndarray] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes the nominal mass, exact mass and abundance of n atoms.

    If the isotope is the monoisotope and p is ``None``, the natural abundances
    for the element are used.

    If the isotope is the monoisotope and custom abundance `p` is provided, the
    envelope is computed using this value instead of the natural abundances.

    If the isotopes is not the monoisotope, it is assumed that only this
    isotope contributes to the envelope.

    """
    symbol = isotope.get_symbol()
    element = PeriodicTable().get_element(symbol)
    is_monoisotope = isotope.a == element.get_monoisotope().a
    n_isotopes = len(element.isotopes)
    if is_monoisotope and (n_isotopes > 1):
        if n == 0:
            M, p = _get_n_isotopes_envelope(isotope, n, max_length)
        elif p is None:
            M, p = _get_n_atoms_natural_abundance(symbol, n, max_length)
        else:
            m, M, _ = element.get_abundances()
            _validate_abundance(p, m, symbol)
            M, p = _get_n_atoms_envelope_aux(m, M, p, n, max_length)
    else:
        M, p = _get_n_isotopes_envelope(isotope, n, max_length)
    return M, p


@lru_cache
def _get_n_atoms_natural_abundance(symbol: str, n: int, max_length: int):
    """
    Computes the envelope of n atoms using the natural abundance.

    aux function to _get_n_atoms_envelope

    """
    m, M, p = PeriodicTable().get_element(symbol).get_abundances()
    return _get_n_atoms_envelope_aux(m, M, p, n, max_length)


def _get_n_atoms_envelope_aux(
    m: np.ndarray, M: np.ndarray, p: np.ndarray, n: int, max_length: int
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes the envelope of n atoms.

    aux function to _get_n_atoms_envelope.

    """
    n_isotopes = p.size
    # find combinations of isotopes that sum n
    combinations = _find_n_isotope_combination(n_isotopes, n)

    # find m, M and p for each combination of isotopes
    multinomial_dist = multinomial(n, p)
    m = np.matmul(combinations, m)
    M = np.matmul(combinations, M)
    p = multinomial_dist.pmf(combinations)

    # sort by exact mass
    sorted_index = np.argsort(M)
    m, M, p = m[sorted_index], M[sorted_index], p[sorted_index]

    # merge values with the same nominal mass
    _, first_occurrence = np.unique(m, return_index=True)
    m_unique = np.zeros(max_length, dtype=m.dtype)
    M_unique = np.zeros(max_length, dtype=M.dtype)
    p_unique = np.zeros(max_length, dtype=p.dtype)
    # add the length of m_unique to include all nominal mass values
    n_unique = first_occurrence.size
    first_occurrence = list(first_occurrence)
    first_occurrence.append(m.size)
    m0 = m[0]
    for k in range(max_length):
        if k < n_unique:
            start = first_occurrence[k]
            end = first_occurrence[k + 1]
            mk = m[start]
            i = mk - m0
            if i < max_length:
                m_unique[i] = mk
                pk = np.sum(p[start:end])
                p_unique[i] = pk
                M_unique[i] = np.sum(M[start:end] * p[start:end]) / pk
    p_unique = p_unique / np.sum(p_unique)
    return M_unique, p_unique


def _fill_missing_nominal(
    m: np.ndarray, M: np.ndarray, p: np.ndarray, max_length: int
) -> Tuple[np.ndarray, np.ndarray]:
    rel_m = m - m[0]
    dm = np.arange(max_length)
    M_filled = np.zeros(max_length, dtype=M.dtype)
    p_filled = np.zeros(max_length, dtype=p.dtype)
    if not np.array_equal(rel_m, dm):
        for k, rel_m_k in enumerate(rel_m):
            if 0 <= rel_m_k < max_length:
                M_filled[rel_m_k] = M[k]
                p_filled[rel_m_k] = p[k]
            else:
                break
        M, p = M_filled, p_filled
    return M, p


def _find_n_isotope_combination(n_isotopes, n):
    """
    Finds combinations of isotopes such that the sum is n.

    aux function to _find_n_atoms_abundances.

    """
    n_ranges = [range(x) for x in ([n + 1] * n_isotopes)]
    combinations = utils.cartesian_product(*n_ranges).astype(int)
    valid_combinations = combinations.sum(axis=1) == n
    combinations = combinations[valid_combinations, :]
    return combinations


def _validate_abundance(p: np.ndarray, m: np.ndarray, symbol: str):
    """
    Checks that user-created abundances are non-negative, normalized to 1 and
    has the same length as the number of stable isotopes.

    aux function to _get_n_atoms_envelope.

    """
    is_all_non_negative = (p >= 0.0).all()
    is_normalized = np.isclose(p.sum(), 1.0)
    is_same_size = p.size == m.size
    if not is_same_size:
        msg = "{} has {} stable isotopes. `p` must have the same size."
        raise ValueError(msg.format(symbol, m.size))
    elif not (is_normalized and is_all_non_negative):
        msg = "`p` elements must be non-negative and their sum normalized to 1."
        raise ValueError(msg)


def _get_n_isotopes_envelope(
    isotope: Isotope, n: int, max_length: int
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Creates the isotopic envelope for n isotopes.

    aux function to _get_n_atoms_envelope.

    """
    M = np.zeros(max_length, dtype=float)
    p = np.zeros(max_length, dtype=float)
    M[0] = isotope.m * n
    p[0] = 1.0
    return M, p



































import numpy as np
import string
from collections import Counter
from copy import copy
from typing import List, Tuple, Dict, Optional





class Formula:
    def __init__(self, *args):
        if len(args) == 1:
            formula_str = args[0]
            formula_str, charge = _parse_charge(formula_str)
            composition = _parse_formula(formula_str)
        else:
            f_composition, charge = args
            ptable = PeriodicTable()
            composition = Counter()
            for k, v in f_composition.items():
                if isinstance(k, str):
                    isotope = ptable.get_isotope(k)
                elif isinstance(k, Isotope):
                    isotope = k
                else:
                    msg = "Composition keys must be a string representation of an isotope or an isotope object"
                    raise InvalidIsotope(msg)

                if not isinstance(v, int) or (v < 1):
                    msg = "Formula coefficients must be positive integers"
                    raise ValueError(msg)
                composition[isotope] = v

            if not isinstance(charge, int):
                msg = "Charge must be an integer"
                raise ValueError(msg)
        self.charge = charge
        self.composition = composition

    def __add__(self, other: "Formula") -> "Formula":
        if not isinstance(other, Formula):
            msg = "sum operation is defined only for Formula objects"
            raise ValueError(msg)
        else:
            # copy Formula object and composition
            sum_composition = copy(self.composition)
            sum_composition.update(other.composition)
            sum_charge = self.charge + other.charge
            sum_f = Formula(sum_composition, sum_charge)
            return sum_f

    def __sub__(self, other: "Formula") -> "Formula":
        if not isinstance(other, Formula):
            msg = "subtraction operation is defined only for Formula objects"
            raise ValueError(msg)
        else:
            comp = copy(self.composition)
            comp.subtract(other.composition)
            charge = self.charge - other.charge
            min_coeff = min(comp.values())
            if min_coeff < 0:
                msg = "subtraction cannot generate negative coefficients"
                raise ValueError(msg)
            comp = Counter({k: v for k, v in comp.items() if v > 0})
            diff_f = Formula(comp, charge)
            return diff_f

    def __eq__(self, other: "Formula"):
        return (self.charge == other.charge) and (self.composition == other.composition)

    def get_exact_mass(self) -> float:
        exact_mass = sum(x.m * k for x, k in self.composition.items())
        exact_mass -= EM * self.charge
        return exact_mass

    def get_isotopic_envelope(
        self,
        n: int = 10,
        p: Optional[Dict[str, np.ndarray]] = None,
        min_p: float = 1e-10,
    ) -> Tuple[np.ndarray, np.ndarray]:
        M, p = find_formula_envelope(self.composition, n, p=p, min_p=min_p)
        M -= EM * self.charge
        return M, p

    def __repr__(self):
        return "Formula({})".format(str(self))

    def __str__(self):
        return _get_formula_str(self.composition, self.charge)


_matching_parenthesis = {"(": ")", "[": "]"}


def _multiply_formula_coefficients(composition: Counter, multiplier: int):
    if multiplier != 1:
        for k in composition:
            composition[k] *= multiplier


def _parse_charge(formula: str) -> Tuple[str, int]:
    """
    compute the charge state of a formula and remove the charge from the formula
    string

    Parameters
    ----------
    formula: str
        molecular formula
    Returns
    -------
    formula_without_charge, charge: str, int

    """

    # check if there's a charge in the formula
    if formula[-1] == "+":
        sign = 1
    elif formula[-1] == "-":
        sign = -1
    else:
        sign = 0

    # for charge with an absolute value greater than 1 enforces the use
    # of parenthesis to prevent ambiguity with a formula coefficient
    if sign and (formula[-2] in string.digits):
        try:
            matching = _matching_parenthesis[formula[0]]
            start = 1
            end = formula.rfind(matching)
            charge = sign * int(formula[end + 1 : -1])
            return formula[start:end], charge
        except (KeyError, ValueError):
            raise InvalidFormula
    elif sign:
        return formula[:-1], sign
    else:
        return formula, 0


def _get_token_type(formula: str, ind: int) -> int:
    """
    assigns 0 to elements, 1 to isotopes and 2 to expressions.
    Return the token type and a matching parenthesis if necessary
    """
    c = formula[ind]
    if c in string.ascii_uppercase:
        token_type = 0
    elif c in "(":
        if formula[ind + 1] in string.digits:
            token_type = 1
        else:
            token_type = 2
    elif c == "[":
        token_type = 2
    else:
        raise InvalidFormula
    return token_type


def _find_matching_parenthesis(formula: str, ind: int):
    parenthesis_open = formula[ind]
    parenthesis_close = _matching_parenthesis[parenthesis_open]
    match_ind = ind + 1
    level = 1
    try:
        while level > 0:
            c = formula[match_ind]
            if c == parenthesis_open:
                level += 1
            elif c == parenthesis_close:
                level -= 1
            match_ind += 1
        return match_ind - 1
    except IndexError:
        msg = "Formula string has non-matching parenthesis"
        raise InvalidFormula(msg)


def _get_coefficient(formula: str, ind: int):
    """
    traverses a formula string to compute a coefficient. ind is a position
    after an element or expression.

    Returns
    -------
    coefficient : int
    new_ind : int, new index to continue parsing the formula
    """
    length = len(formula)
    if (ind >= length) or (formula[ind] not in string.digits):
        coefficient = 1
        new_ind = ind
    else:
        end = ind + 1
        while (end < length) and (formula[end] in string.digits):
            end += 1
        coefficient = int(formula[ind:end])
        new_ind = end
    return coefficient, new_ind


def _tokenize_element(formula: str, ind: int):
    length = len(formula)
    if (ind < length - 1) and (formula[ind + 1] in string.ascii_lowercase):
        end = ind + 2
    else:
        end = ind + 1
    symbol = formula[ind:end]
    isotope = PeriodicTable().get_isotope(symbol)
    coefficient, end = _get_coefficient(formula, end)
    token = {isotope: coefficient}
    return token, end


def _tokenize_isotope(formula: str, ind: int):
    """
    Convert an isotope substring starting at `ind` index into a token.

    Returns
    -------
    token, new_ind

    """
    end = _find_matching_parenthesis(formula, ind)
    isotope = PeriodicTable().get_isotope(formula[ind + 1 : end])
    coefficient, end = _get_coefficient(formula, end + 1)
    token = {isotope: coefficient}
    return token, end


def _parse_formula(formula: str):
    """
    Parse a formula string into a Counter that maps isotopes to formula
    coefficients.
    """
    ind = 0
    n = len(formula)
    composition = Counter()
    while ind < n:
        token_type = _get_token_type(formula, ind)
        if token_type == 0:
            token, ind = _tokenize_element(formula, ind)
        elif token_type == 1:
            token, ind = _tokenize_isotope(formula, ind)
        else:
            # expression type evaluated recursively
            exp_end = _find_matching_parenthesis(formula, ind)
            token = _parse_formula(formula[ind + 1 : exp_end])
            exp_coefficient, ind = _get_coefficient(formula, exp_end + 1)
            _multiply_formula_coefficients(token, exp_coefficient)
        composition.update(token)
    return composition


# functions to get a formula string from a Formula


def _arg_sort_elements(symbol_list: List[str], mass_number_list: List[int]):
    """
    Return the sorted index for a list of elements symbols and mass numbers. If
    there are repeated elements, they are sorted by mass number.
    """
    zipped = list(zip(symbol_list, mass_number_list))
    return sorted(range(len(symbol_list)), key=lambda x: zipped[x])


class InvalidFormula(ValueError):
    pass


def _get_formula_str(composition, charge):

    # get C str
    c12 = PeriodicTable().get_isotope("12C")
    c13 = PeriodicTable().get_isotope("13C")
    h1 = PeriodicTable().get_isotope("1H")
    h2 = PeriodicTable().get_isotope("2H")
    ch = [c12, c13, h1, h2]

    # ADD C and H to formula str
    f_str = ""
    for i in ch:
        if i in composition:
            f_str += isotope_coeff_to_f_str(i, composition[i])

    # add other elements, sorted alphabetically
    isotopes = set(composition)
    isotopes = isotopes.difference(ch)
    for i in sorted(isotopes, key=lambda x: x.get_symbol() + str(x.a)):
        f_str += isotope_coeff_to_f_str(i, composition[i])

    if charge:
        charge_str = _get_charge_str(charge)
        f_str = "[{}]{}".format(f_str, charge_str)

    return f_str


def isotope_coeff_to_f_str(isotope: Isotope, coeff: int) -> str:
    coeff_str = str(coeff) if coeff > 1 else ""
    element = isotope.get_element()
    if isotope.a == element.nominal_mass:
        isotope_str = element.symbol
    else:
        isotope_str = "({})".format(isotope)
    return "{}{}".format(isotope_str, coeff_str)


def _get_charge_str(q: int) -> str:
    qa = abs(q)
    q_sign = "+" if q > 0 else "-"
    q_str = str(qa) if qa > 1 else ""
    q_str = q_str + q_sign
    return q_str
