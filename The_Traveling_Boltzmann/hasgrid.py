from collections import defaultdict
import numpy as np
import networkx as nx

class HashGrid(object):
    def __init__(self, G):
        self.xy = np.asarray([G.nodes[n]['pos'] for n in range(len(G.nodes))])
        x_min, x_max = np.min(self.xy[:,0]), np.max(self.xy[:,0])
        y_min, y_max = np.min(self.xy[:,1]), np.max(self.xy[:,1])
        self.grid_size = np.sqrt((x_max - x_min) * (y_max - y_min) / len(G.nodes))
    
        self.i_max = int(x_max // self.grid_size)
        self.i_min = int(x_min // self.grid_size)
        self.j_max = int(y_max // self.grid_size)
        self.j_min = int(y_min // self.grid_size)
        
        self.dict = defaultdict(set)
        self.fill(self.xy)
            
    def xy_to_ij(self, x, y):
        return (x // self.grid_size,
                y // self.grid_size)
    
    def ij_to_key(self, i, j):
        return self.i_max * j + i
    
    def xy_to_key(self, x, y):
        return self.ij_to_key(*self.xy_to_ij(x,y))
    
    def key_to_ij(self, key):
        return (key % self.i_max,
                key // self.i_max)
    
    def fill(self, xy):
        node = 0
        for coord in  xy: # don't know how to enumerate np.arrays
            x, y = coord
            key = self.xy_to_key(x, y)
            self.dict[key].add(node)
            node += 1
        
    
    
    
        
    def get_nearby_cities(self, node, reach=1):
        x, y = self.xy[node]
        key = self.xy_to_key(x, y)
        i, j = self.key_to_ij(key)
        r = int(reach // self.grid_size)
        
        nbrs = self.dict[key]
        for ii in range(int(max(self.i_min, i - r)),
                        int(min(self.i_max, i + r) + 1)):
            for jj in range(int(max(self.j_min, j - r)),
                            int(min(self.j_max, j + r) + 1)):
                nbrs = set.union(nbrs,self.dict[self.ij_to_key(ii,jj)])
                
        return nbrs - {node}
