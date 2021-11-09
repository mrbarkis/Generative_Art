from vector import Vector
class Edge(object):
    def __init__(self, v1, v2):
        self.v1, self.v2 = v1, v2
        self.t = self.v2 - self.v1 # tangent vector
        
    def edge_point(self, s):
        return self.v1 + (self.v2 - self.v1) * s
        
        
class Polygon(object):
    def __init__(self, vs, edge_color=[0,0,0], fill_color=None):
        self.vertices = vs
        self.edges = [Edge(vs[i], vs[(i + 1) % len(vs)]) for i in range(len(vs))]
        self.edge_color, self.fill_color = edge_color, fill_color
        
        
class Box(Polygon):
    def __init__(self, x, y, w, h, edge_color=[0,0,0], fill_color=None):
        super(Box, self).__init__(
            vs = [Vector(x, y), Vector(x+w, y), Vector(x+w, y+h), Vector(x, y+h)],
            edge_color=edge_color,
            fill_color=fill_color
            )

class Bar(Box):
    def __init__(self, x, y, w, h, h_target, speed, edge_color=[0,0,0], fill_color=[0,0,0]):
        super(Bar, self).__init__(x, y, w, h, edge_color, fill_color)
        self.x, self.y, self.w, self.h = x, y, w, h
        self.h_target = h_target
        self.speed = speed
        
    def get_target_bar(self):
        if self.h == self.h_target:
            return None
        target_step = self.h_target - self.h
        target_step_length = abs(target_step)
        possible_step_size = min(target_step_length, self.speed)
        step = target_step / target_step_length * possible_step_size
        return Bar(self.x, self.y - step, self.w, self.h + step, self.h_target, self.speed, self.edge_color, self.fill_color)

class BarChart(object):
    def __init__(self,pmf, pos, height_, width_, speed, n_bins, n_samples, edge_color=[0,0,0], fill_color=[0,0,0], active_color = 3*[0]):
          w_bar = width_ / n_bins
          
          self.H = height_
          self.N_b = n_bins
          self.N_s = n_samples
          self.fill_color = fill_color
          self.edge_color = edge_color
          self.active_color = active_color
  
          self.bars = []
          
          for bin in range(self.N_b):
              h_i = pmf(bin) * self.H
              self.bars.append(Bar(
                              x=pos[0] + bin * w_bar,
                              y=pos[1] + height_ - h_i, 
                              w=w_bar,
                              h=h_i,
                              h_target=h_i,
                              speed = speed,
                              edge_color = edge_color,
                              fill_color = fill_color))
              
    def update_pmf(self, pmf):
        for bin in range(self.N_b):
            self.bars[bin].h_target = pmf(bin) * self.H
        
              
    def add_to(self, bin):
        if bin >= self.N_b: # Larger values are truncated to the largest bin.
            return
        
        if bin == 0: # Lover all the bars to make room for a new one
            
            for bar in self.bars:
                bar.edge_color = self.edge_color
                bar.fill_color = self.fill_color
                bar.h_target *= float(self.N_s) / (self.N_s + 1)
                
        else: # Lover only the previous bar
            self.bars[bin - 1].edge_color = self.edge_color
            self.bars[bin - 1].fill_color = self.fill_color
            self.bars[bin - 1].h_target -= 1.0 / (self.N_s + 1) * self.H
        
        #self.bars[bin].edge_color = self.active_color
        self.bars[bin].fill_color = self.active_color  
        self.bars[bin].h_target += 1.0 / (self.N_s + 1) * self.H
        
                              
          
    
