from abc import ABCMeta, abstractmethod

class Particle:
    __metaclass__ = ABCMeta #processing supports only python 2.7.
    
    @abstractmethod
    def move(self, dframe=1):
        pass
    
    @abstractmethod
    def draw_myself(self):
        pass
    
    @abstractmethod
    def on_colliding(self):
        ''' Customize appearance or trigger counters when a CollisionEvent is handled.
        Note that the collision event takes care of the collision physics.
        '''
        pass


class Ball(Particle):
    def __init__(self, p, v, m, r, c, painter):
        self.p, self.v, self.m, self.r = p, v, m, r
        self.painter = painter
        self.c = c # RGB colors
        
    def move(self, dframe=1):
        self.p += self.v * dframe
        
    def draw_myself(self):
        self.painter(self)
    
    def on_colliding(self):
        pass
        
class Pollen(Ball):
    def __init__(self, p, v, m, r, c, painter, interval, n_samples, when_colliding, when_resetting):
        super(Pollen, self).__init__(p, v, m, r, c[0], painter)
        self.interval = interval
        self.frame = 1
        self.interval_count = 0
        self.c_t = c[2]
        self.c_arc = c[1]
        self.n_collisions = 0
        self.rate_ave = 0
        self.n_samples = n_samples
        self.when_colliding = when_colliding
        self.when_resetting = when_resetting
        # self.clicker = clicker
        # self.counter = counter
        
    def move(self, dframe=1):
        self.p += self.v * dframe
        self.frame += dframe
        if self.frame // self.interval > self.interval_count:
            self.interval_count += 1
            self.on_reset()
    
    def on_reset(self):
        self.rate_ave += float(self.n_collisions - self.rate_ave) / self.n_samples
        self.n_collisions = 0
        self.c, self.c_arc = self.c_arc, self.c # Flip colors
        for f in self.when_resetting:
            f(self)
        # self.counter(0)
        
    
    def on_colliding(self):
        self.n_collisions += 1
        for f in self.when_colliding:
            f(self)
            
        # self.counter(self.n_collisions)
        # self.clicker()
        
        #print("n_collisions = {}".format(self.n_collisions))
        
        
    
        
        
