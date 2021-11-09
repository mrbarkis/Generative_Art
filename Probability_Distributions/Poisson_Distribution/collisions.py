from abc import ABCMeta, abstractmethod
from particles import Ball

class CollisionEvent:
    __metaclass__ = ABCMeta #processing supports only python 2.7.
    
    def __init__(self, eta, c1, c2):
        self.eta, self.c1, self.c2 = eta, c1, c2
        
    def advance(self, dframe):
        self.eta -= dframe
    
    def __lt__(self, operand):
        if operand == None:
            return True
        return self.eta < operand.eta
    
    # def __eq__(self, operand):
    #     return 
    # def get_common_balls(self, collision):
    #     ''' Return balls that also belong to the other collision. '''
    #     return self.get_balls{}.intersection(collision.get_balls{})
    
    def involves(self, other):
        return self.c1 == other or self.c2 == other
        
    def get_balls(self):
        return {b for b in [self.c1, self.c2] if isinstance(b, Ball)}
    
    @abstractmethod
    def handle(self):
        pass

class BallBallCollision(CollisionEvent):
    def handle(self):
        ball1, ball2 = self.c1, self.c2
        c = (ball2.p - ball1.p);
        c = c * (1 / c.norm())
        delta_v = ball2.v - ball1.v
        m = ball1.m + ball2.m
        delta_projected =  c * (delta_v * c) * 2 * (1./ m)
        #print('c: {}, delta_v: {}, division: {}'.format(c, delta_v, (m)))
        self.c1.v, self.c2.v = (ball1.v + delta_projected * ball2.m,
                            ball2.v - delta_projected * ball1.m)
        self.c1.on_colliding()
        self.c2.on_colliding()
        
        
        
    def __repr__(self):
        return "BallBallCollision: eta {:.2f}".format(self.eta)
    


class BallEdgeCollision(CollisionEvent):
    def handle(self):
        '''Reverse the velocity component perpendicular to the edge.'''
        ball, edge = self.c1, self.c2
        n = edge.t.normal_vector()
        ball.v -=  n * (2 * (ball.v * n))
        
    def __repr__(self):
        return "BallEdgeCollision: eta {:.2f}".format(self.eta)
    
class BallVertexCollision(CollisionEvent):
    def handle(self):
        '''Reverse the velocity component perpendicular to the tangent at the hit point'''
        ball, v = self.c1, self.c2
        n = ball.p + ball.v * self.eta - v
        n = n * (1 / n.norm())
        ball.v -=  n * (2 * (ball.v * n))
        
    def __repr__(self):
        return "BallVertexCollision: eta {:.2f}".format(self.eta)
