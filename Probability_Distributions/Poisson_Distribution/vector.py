from math import sqrt

class Vector(object):
    def __init__(self, x, y):
        self.x, self.y = float(x), float(y)
    def __repr__(self):
        return 'Vector({}, {})'.format(self.x, self.y)
    
    def normal_vector(self):
        v = Vector(self.y, - self.x)
        return v * (1 / v.norm())
    
    def norm(self):
        return sqrt(self.x ** 2 + self.y **2)
    
    def normalize(self):
        return self * (1 / self.norm())
    
    def __add__(self, operand):
        return Vector(x=self.x+operand.x, y=self.y+operand.y)
    
    def __sub__(self, operand):
        return Vector(self.x - operand.x, self.y - operand.y)
    
    def __mul__(self, operand):
        if isinstance(operand, Vector): # Dot product
            return (self.x * operand.x + self.y * operand.y)
        elif isinstance(operand, int) or isinstance(operand, float): #scalar product
            return Vector(self.x * operand, self.y * operand)

    def __rmul__(operand, self):
        return self.__mul__(operand)
    
    def T(self):
        return Vector(self.y, self.x)
    
    
    def get_xy(self):
        return (self.x, self.y)
