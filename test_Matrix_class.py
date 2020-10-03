# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 13:38:43 2020

@author: Sarah
"""

import unittest
import numpy as np
from Matrix_class import Matrix as M

#two 3x3 matrices for testing the matrix class methods
m1 = np.array([[1, 1, 1], [2, 2, 2], [3, 3, 3]])
m2 = np.array([[4, 4, 4], [5, 5, 5], [6, 6, 6]])


class TestMatrix(unittest.TestCase):
    """ Unit tester for Matric_class.py """

    def test_add(self):
        """ Test add Method """
        self.assertAlmostEqual(M.__add__(m1,m2), np.add(m1,m2))
    
    def test_Mult(self):
        """ Test Mult Method """
        self.assertAlmostEqual(M.__mult__(m1,m2), np.dot(m1,m2))

    def test_transpose(self):
        """ Test Tran Method """
        self.assertAlmostEqual(M.__tran__(m1), m1.T)
    
    def test_invert(self):
        """ Test invert Method """
        self.assertAlmostEqual(M.__invert__(m1), np.linalg.inv(m1))
        
    def test_trace(self):
        """ Test trace Method """
        self.assertAlmostEqual(M.__trace__(m1), np.trace(m1))
        
    def test_determinant(self):
        """ Test determinant Method """
        self.assertAlmostEqual(M.determinant(m1), np.linalg.det(m1))
        
    def test_lu_decomp(self):
        """ Test lu_decomposition Method """
        self.assertAlmostEqual(M.lu_decomposition(m1), np.linalg.lu(m1))
        

if __name__ == '__main__':
    unittest.main()      