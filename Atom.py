from Helper import Vector3


__author__ = "Ilia Kichev"
__credits__ = ["Ilia Kichev", "Lyuben Borislavov", "Alia Tadjer"]
__version__ = "1.1.0"
__maintainer__ = "Ilia Kichev"
__email__ = "ikichev@uni-sofia.bg"


class Atom:
    def __init__(self, coord: Vector3 = Vector3(), atomic_number: int = 0,
                 atomic_symbol: str = "", atomic_mass: float = 0, index: int = 0):
        self.__coord: Vector3 = coord
        self.__atomicMass: float = atomic_mass
        self.__atomicNumber: int = atomic_number
        self.__symbol: str = atomic_symbol
        self.__index: int = index
    
    def __str__(self):
        return "{}\t{}".format(self.symbol, self.coord)
        
    # Getters and setters
    @property
    def coord(self) -> Vector3:
        return self.__coord
    
    @coord.setter
    def coord(self, coordinates: Vector3):
        if not isinstance(coordinates, Vector3):
            raise TypeError("Parameter coordinates is not Vector3")
        self.__coord = coordinates

    @property
    def symbol(self) -> str:
        return self.__symbol
    
    @symbol.setter
    def symbol(self, s: str):
        if not isinstance(s, str):
            raise TypeError("Parameter s is not string")
        self.__symbol = s

    @property
    def atomic_num(self) -> int:
        return self.__atomicNumber
    
    @atomic_num.setter
    def atomic_num(self, z: int):
        if not isinstance(z, int):
            raise TypeError("Parameter z is not int")
        self.__atomicNumber = z
    
    @property
    def atomic_mass(self) -> float:
        return self.__atomicMass
    
    @atomic_mass.setter
    def atomic_mass(self, am: float):
        if not isinstance(am, float):
            raise TypeError
        self.__atomicMass = am
    
    @property
    def id(self) -> int:
        return self.__index
    
    @id.setter
    def id(self, ind: int):
        if not isinstance(ind, int):
            raise TypeError("The parameter ind is not int")
        self.__index = ind

    def __eq__(self, other: 'Atom'):
        return self.__coord == other.coord and self.atomic_num == other.atomic_num and \
               self.symbol == other.symbol and self.id == other.id
    