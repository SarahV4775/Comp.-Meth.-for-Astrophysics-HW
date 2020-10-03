# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 19:28:24 2020

@author: Sarah
"""

class Matrix:
    
    def __init__(self, dim, fill_value):
        #initializing the matrix
        self.height = dim[0]
        self.width = dim[1]
        self.rows = [[fill_value] * self.width for i in range(self.height)]
        
        
    def __str__(self) -> str:
        #will make it print the matrix in a readable format
        return "\n".join(map(str, self.rows))
    
    
    def __add__(self, other):
        #Adding two matrix or a matrix and an int or float
        #Create a new matrix
        C = Matrix( dim = (self.height, self.width), fill_value = 0)
        
        #Check if the other object is of type Matrix
        if isinstance (other, Matrix):
			#Add the corresponding element of 1 matrices to another
            for i in range(self.height):
                for j in range(self.width):
                    C.rows[i][j] = self.rows[i][j] + other.rows[i][j]
        
        #if the other object is an int or float
        elif isinstance (other, (int, float)):
            for i in range(self.height):
                for j in range(self.width):
                    #add the other object to ever instance of self
                    C.rows[i][j] = self.rows[i][j] + other
        return C
    
    
    def __mult__(self, other):
        #multiplying two matrix or a matrix and an int or float
        #Create a new matrix
        C = Matrix( dim = (self.height, self.width), fill_value = 0)
        
        #matrix*matrix
        #Check if the other object is of type Matrix
        if isinstance (other, Matrix):
			#multiply the corresponding element of 1 matrices to another
            #row of one * col of the two
            for i in range(self.height):
                for j in range(self.width):
                    new = 0
                    for k in range(self.height):
                        new += self.rows[i][k] * other.rows[k][j]
                    C.rows[i][j] = new
        
        #if the other object is an int or float
        elif isinstance (other, (int, float)):
            for i in range(self.height):
                for j in range(self.width):
                    #multiply the other object to ever instance of self
                    C.rows[i][j] = self.rows[i][j] * other
        return C
    
    def __tran__(self):
        #return the transpose of the matrix
        #Create a new matrix
        C = Matrix( dim = (self.height, self.width), fill_value = 0)
        
        for i in range(self.height):
            for j in range(self.width):
                C.rows[i][j] = self.rows[j][i]
        
        return C
    
    def __invert__(self):
        #return the inverse of the matrix
        import numpy as np
        I = np.linalg.inv(self)
        return I
    
    def __trace__(self):
        #return the inverse of the matrix
        import numpy as np
        T = np.trace(self)
        return T
    
    def determinant(self):
        #base case for 2x2 matrix
        if len(self) == 2:
            return self[0][0]*self[1][1]-self[0][1]*self[1][0]
    
        det = 0
        for c in range(len(self)):
            det += ((-1)**c)*self[0][c]*Matrix.determinant(self[1:])
        return det

    def lu_decomposition(self):
        import numpy as np
        n = len(self.height)
        
        #create new and empty matrices for L and U
        L = Matrix( dim = (self.height, self.width), fill_value = 0)
        U = Matrix( dim = (self.height, self.width), fill_value = 0)
        
        P = np.identity(n)
        PA = Matrix.__mult__(P, self)
        
        for j in range(n):
            L[j][j] = 1.0
            
            for i in range(j+1):
                sum1 = sum(U[k][j] * L[i][k] for k in range(i))
                U[i][j] = PA[i][j] - sum1
            for i in range(j, n):
                sum2 = sum(U[k][j] * L[i][k] for k in range(j))
                L[i][j] = (PA[i][j] - sum2)/U[j][j]
                
        return (L, U)
