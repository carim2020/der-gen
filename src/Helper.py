from __future__ import annotations
import numpy as np
import math


__author__ = "Ilia Kichev"
__credits__ = ["Ilia Kichev", "Lyuben Borislavov", "Alia Tadjer"]
__version__ = "1.1.0"
__maintainer__ = "Ilia Kichev"
__email__ = "ikichev@uni-sofia.bg"


class Vector3:
    MIN_DIST = 0.0000001
    
    def __init__(self, x: float = 0., y: float = 0., z: float = 0.):
        """
                Constructor of the Vector3 class
                @param x: The x coordinate
                @param y: The y coordinate
                @param z: The z coordinate"""
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    @classmethod
    def from_numpy(cls, nparray: np.array):
        return cls(float(nparray[0]), float(nparray[1]), float(nparray[2]))

    def __str__(self):
        return "{}\t{}\t{}".format(self.x, self.y, self.z)
    
    def length_2(self):
        return (self.x ** 2) + (self.y ** 2) + (self.z ** 2)
    
    def length(self):
        return math.sqrt(self.length_2())

    def to_numpy(self) -> np.array:
        return np.array([self.x, self.y, self.z])

    # Operator override
    def __add__(self, other):
        x = self.x + other.x
        y = self.y + other.y
        z = self.z + other.z
        return Vector3(x, y, z)
    
    def __sub__(self, other):
        x = self.x - other.x
        y = self.y - other.y
        z = self.z - other.z
        return Vector3(x, y, z)
    
    def __mul__(self, num: float):
        x = self.x * num
        y = self.y * num
        z = self.z * num
        return Vector3(x, y, z)
    
    def __truediv__(self, num: float):
        x = self.x / num
        y = self.y / num
        z = self.z / num
        return Vector3(x, y, z)
    
    def __iadd__(self, other):
        self.x += other.x
        self.y += other.y
        self.z += other.z
        return self
    
    def __isub__(self, other):
        self.x -= other.x
        self.y -= other.y
        self.z -= other.z
        return self
    
    def __imul__(self, num: float):
        self.x *= num
        self.y *= num
        self.z *= num
        return self
    
    def __itruediv__(self, num: float):
        self.x /= num
        self.y /= num
        self.z /= num
        return self

    def __eq__(self, other):
        temporary = self - other
        if temporary.length_2() < Vector3.MIN_DIST:
            return True
        else:
            return False

    def rotate_around_axis(self, axis: 'Vector3', angle: float) -> 'Vector3':
        t = 1 - math.cos(angle)
        C = math.cos(angle)
        S = math.sin(angle)
        axis_norm = axis / axis.length()

        R = np.array([[(t * axis_norm.x**2 + C), (t * axis_norm.x * axis_norm.y - S * axis_norm.z), (t * axis_norm.x * axis_norm.z + S * axis_norm.y)],
                      [(t * axis_norm.x * axis_norm.y + S * axis_norm.z), (t * axis_norm.y**2 + C), (t * axis_norm.y * axis_norm.z - S * axis_norm.x)],
                      [(t * axis_norm.x * axis_norm.z - S * axis_norm.y), (t * axis_norm.y * axis_norm.z + S * axis_norm.x), (t * axis_norm.z**2 + C)]])
        return Vector3.from_numpy(np.transpose(R @ np.transpose(self.to_numpy())))


def cross_product(lhs: Vector3, rhs: Vector3) -> Vector3:
    return Vector3.from_numpy(np.cross(lhs.to_numpy(), rhs.to_numpy()))


def dot_product(lhs: Vector3, rhs: Vector3) -> float:
    return float(np.dot(lhs.to_numpy(), rhs.to_numpy()))


