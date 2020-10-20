# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 14:07:22 2020

@author: Sarah
"""
import matplotlib.pyplot as plt
import numpy as np

# Octtree Node Class
class Node:

    def __init__(self, size, cen, points, ids, masses, leaves=[]):
        '''
        Initialize the the structure to reproduce and store the galaxies in 
        each node
        ----------------
        size: the size of the global node encompaceing all the data
        cen: coordinates of the center node where we start
        points: array of the coordinates of the points or galaxies
        ids: identifies of the nodes
        masses: array of the mass of each galaxy (for this project all masses
                are the same)
        leaves: list of all the leaves initalized to be empty to keep track of
                whenever a node reaches a leaf
        '''
        self.size = size
        self.cen = cen                  
        #initalize to no children
        self.children = []
        
        #while there is only one point in the node, it is a leaf and we will
        #store stuff in the node
        if len(points) == 1:
            #update the leaves list
            leaves.append(self)
            #calculate the Center of Mass (COM) and the mass (mass)
            self.COM = points[0]
            self.mass = masses[0]
            #keep track of the id
            self.id = ids[0]
            self.g = np.zeros(3)
            
        #more then one point in the node
        else:
            #spawn children
            self.get_children(points, masses, ids, leaves)
            
            com_total = np.zeros(3) # running total for COM
            m_total = 0.            # running total for masses
            for c in self.children:
                m = c.mass
                com = c.COM
                #update the tatal mass and total COM
                m_total += m
                com_total += com * m
            self.mass = m_total
            self.COM = com_total / self.mass  
 
    def get_children(self, points, masses, ids, leaves):
        """
        Generate the Children of the node
        ----------
        points: array of the coordinates of the points or galaxies
        ids: identifies of the nodes
        masses: array of the mass of each galaxy (for this project all masses
                are the same)
        leaves: list of all the leaves to keep track of whenever a node reaches a leaf
        """
        #index of the children
        index = (points > self.cen)
        
        #loop through 8 children
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    
                    #check if there is a galaxy in the child node
                    in_node = np.all(index == np.bool_([i,j,k]), axis=0)
                    
                    #if theres no point in the node
                    if not np.any(in_node): 
                        continue
                    # offset between parent and child box centers
                    offset = 0.5*self.size*(np.array([i,j,k])-0.5)   
                    self.children.append(Node(self.cen + offset,
                                              self.size/2,
                                              masses[in_node],
                                              points[in_node],
                                              ids[in_node],
                                              leaves))        

thetamax=0.7
G=1.0
#force softening threshold
fs = 0.001
def tree_contribution(node, other):
    """
    contribution of a node on another node
    ---------------
    node: the node that will effect the other object
    other: the node whose information is effected by the first node
    """
    # difference between the two nodes COM
    dx = node.COM - other.COM 
    # distance
    r = np.sqrt(np.sum(dx**2))
    if r>0:
        #if the node only has one particle or theta is small enough,
        #add the field contribution to value stored in node.g
        if (len(node.children)==0):
            #add force softening
            #if its within the threshold force sofenting is applied
            if r <= fs:
                other.g += G * node.mass * dx/(r**3 + (r*fs**2))
            elif r > fs:
                other.g += G * node.mass * dx/r**3
 
        elif (node.size/r < thetamax):
            if r <= fs:
                other.g += G * node.mass * dx/(r**3 + (r*fs**2))
            elif r > fs:
                other.g += G * node.mass * dx/r**3
    
        else:
            # otherwise split up the node and repeat
            for c in node.children: tree_contribution(c, other)

def Accel(points, masses):
    '''
    Calculates the Gravitational acceleration of each point based on the 
    other points
    ---------
    points: array of the coordinates of the points or galaxies
    masses: array of the mass of each galaxy (for this project all masses
            are the same)
    
    returns a list of accelerations for each point
    '''
    # Center of the box the the points are bound in
    center = (np.max(points, axis=0) + np.min(points, axis=0))
    # Size of box
    bound_size = np.max(np.max(points, axis=0) - np.min(points, axis=0))
    #initalize the empty list of leave to keep track of them
    leaves = []
    # Node that will build the tree up from
    Start_node = Node(bound_size, center, points, np.arange(len(masses)), masses, leaves)
    
    #initialize the place to store the Acceleration
    accel = np.empty_like(points)
    for i,leaf in enumerate(leaves):
        #sum of all the acceleration on all the galaxies in the array
        tree_contribution(Start_node, leaf)
        accel[leaf.id] = leaf.g
    
    return accel
        