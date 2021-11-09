from collisions import BallEdgeCollision, BallVertexCollision, BallBallCollision
from custom_exceptions import BoundException
from model import Model

### Help Functions ###
def ball_to_edge_distance(ball, edge):
    p, a = ball.p, edge.v1, 
    n = edge.t * (1 / edge.t.norm())
    return ((p - a) - n * ((p - a) * n)).norm()

def intersection(ball, edge):
    ''' Solve when and where the center of mass of the ball meets the edge.
    That is, return t, and s such that:
        ball.p + ball.v * t == edge.v1 + edge.t * s
    Return None if there is no solution.
    '''
    denom = ball.v.x * edge.t.y - ball.v.y * edge.t.x
    try:        
        return (((edge.v1.x - ball.p.x) * edge.t.y - (edge.v1.y - ball.p.y) * edge.t.x) / denom,
                ((ball.p.y - edge.v1.y) * ball.v.x - (ball.p.x - edge.v1.x) * ball.v.y) / denom)
    except ZeroDivisionError:
        return None, None
    

def predict_ball_hit(ball1, ball2):
    """Solve when the balls collide. Returns None if they do not."""
    if ball1 == ball2: return None
    r_1 = ball1.r
    r_2 = ball2.r
    dv = ball1.v - ball2.v
    dq = ball1.p - ball2.p
    #a = (vx_1 - vx_2)**2 + (vy_1 - vy_2)**2
    a = dv * dv
    #b = 2*((vx_1- vx_2)*(x_1 - x_2) + (vy_1- vy_2)*(y_1 -y_2))
    b = (dv * dq) * 2
    #c = -(r_1 + r_2)**2 + (x_1 - x_2)**2 + (y_1 - y_2)**2
    c = dq * dq - (r_1 + r_2)**2 
    D = -4*a*c + b**2 # discriminant of the quadratic equation
    if D < 0 or a == 0:
        return None
    
    t1 = -(b + sqrt(D))/(2*a) # Time of the first collision
    t2 = -(b - sqrt(D))/(2*a) # Time of the second collision if the balls could overlap    
    
    if t1 < 0 and t2 < 0: # "Collision in the past"
        return None
    
    if t1 > 0: # Collision in the future
        return BallBallCollision(eta=t1, c1=ball1, c2=ball2)
    elif t2 >0:# and b < 0: # Balls overlap and still getting closer
        raise BoundException([ball1, ball2], message=' overlapped by {} pixels'.format(dq.norm() - ball1.r - ball2.r))
    
    
def predict_earliest_ball_hit(ball, model):
    earliest_c = None
    for b in model.balls:
        c = predict_ball_hit(ball, b)
        if c and c < earliest_c:
            earliest_c = c
    return earliest_c
        
    
    
def predict_vertex_hit(ball, v):
    n = ball.v.normal_vector()
    dr = (v - ball.p)
    d_now = dr.norm() # distance between cm of the ball and the vertex.
    if d_now < ball.r:
        raise BoundException([ball, v], message=' overlapped by {} pixels'.format(d_now-ball.r))
    
    d_min = abs(n * dr) # distance between trajectory of the ball and the vertex.

    if d_min >= ball.r:
        return None
    
    t = (ball.v * dr) / ball.v.norm()**2 - sqrt(ball.r**2 - d_min**2) / ball.v.norm()
    if t < 0:
        return None
    
    return BallVertexCollision(eta = t, c1 = ball, c2 = v)


def predict_edge_hit(ball, edge):
    """ Predict if the ball hits the edge. Return None if it does not;
    Raise a BoundException if the ball overlaps with the edge. Vertex collisions
    are ignored, as they are handled separately."""
    
   # Special Case 1: Overlap -> BoundException
    if overlaps(ball, edge):
        raise BoundException([ball, edge], message='Ball overlaps an edge.')
    

    # Special Case 2: No collision in the future -> None
    t, s = intersection(ball, edge) # time when cm trajectory meets edge:
        
    if t == None or t < 0: return None

    
    # Backtrack t and s, to account for the radius of the ball
    abs_t, abs_v = edge.t.norm(), ball.v.norm()
    cosa = (edge.t * ball.v / abs_t / abs_v) # Cosine of the intersection angle
    sina = sqrt(1 - cosa ** 2) # Sine of the intersection angle
    
    try: # Fails if the ball and the edge are too parallel
        t -= ball.r / abs_v / sina 
    except ZeroDivisionError:
        return None # Collision at infinity
        
    try: # Fails if the ball hits the edge headon
        tana = sina / cosa
        s -= ball.r / abs_t / tana
    except ZeroDivisionError:
        pass # No need to backtrack s during headon collisions.
    
    # Check if the ball hit the edge
    if 0 < s < 1 and t > 0: # Ball does not pass the edge.
        return BallEdgeCollision(eta = t, c1 = ball, c2 = edge)
    else:
        return None

def overlaps(ball, edge):
    q, t = edge.v1 - ball.p, edge.t
    a = t * t
    if a == 0: # Edge is actually a point and we can ignore it.
        return False
    
    b =  (q * t) * 2
    c = q * q - ball.r ** 2
    
    d = b ** 2 - 4 * a * c
    if d < 0:
        return False
    
    s1 = (-b - sqrt(d)) / (2 * a)
    s2 = (-b + sqrt(d)) / (2 * a)
    if 01 < s1 < 1 or 0 < s2 < 1:
        return True


def predict_earliest_bound_hit(ball, bounds):
    if isinstance(bounds, Model):
        bounds = bounds.bounds
    
    earliest_c = None # Alternatively: BallEdgeCollision(eta=float('inf'), c1=None, c2=None)
    for p in bounds:
        for e in p.edges:
            c = predict_edge_hit(ball, e)
            if c and c < earliest_c:
                earliest_c = c
                
    for p in bounds:
        for v in p.vertices:
            c = predict_vertex_hit(ball, v)
            if c and c < earliest_c:
                earliest_c = c        
                            
    return earliest_c
