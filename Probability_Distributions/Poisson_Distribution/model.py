from heapq import merge

def truncate(f):
    return float(int(f * 100)) / 100

class Model(object):
    def __init__(self, balls, pollens, bounds):
        self.balls = balls + pollens
        self.pollens = pollens
        self.bounds = bounds
        self.collisions = []
        self.next_collision = None
        
    def advance(self, dframe):
        dframe = truncate(dframe) # Ensures that we never overshoot due to floating point errors
        for ball in self.balls:
                ball.move(dframe=dframe)
        for c in self.collisions:
                c.advance(dframe=dframe)
                
        #self.next_collision.advance(dframe=dframe)
                
    def remove_collisions(self, collisions):
        self.collisions = filter(lambda c: c not in collisions, self.collisions)
    
    def merge_collisions(self, collisions):
        self.collisions = list(merge(self.collisions, collisions))    
