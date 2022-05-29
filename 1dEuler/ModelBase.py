from __future__ import annotations


class PhysicalPropertyBase():
    def __init__(self, length:int):
        if type(length) is not int:
            raise TypeError
        
        self.value = [0. for _ in range(length)]
        self.leftBoundaryValue = [0. for _ in range(length+1)]
        self.rightBoundaryValue = [0. for _ in range(length+1)]
        
    def F(self) -> float:
        raise NotImplementedError
            