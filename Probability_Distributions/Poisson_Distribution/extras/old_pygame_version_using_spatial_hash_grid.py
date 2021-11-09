import pygame
import numpy as np
from itertools import combinations
import time
import random

# Options
height = 360
width = 640
particle_count = 100
radius = 8
velocity = 5 #  initial velocity = (np.random.rand(2) -0.5) * velocity
acceleration = 0

# Initialize pygame and global variables
pygame.init()
clock = pygame.time.Clock()
particles = []
old_rects = []

# Initialize grid used to speed collision detection.
CELL_SPAN = radius / np.sqrt(2) # Only one ball will fit a cell.
MAX_COL, MAX_ROW = (np.ceil(width / CELL_SPAN).astype(int),
                    np.ceil(height / CELL_SPAN).astype(int))
grid = [None for i in (range(MAX_ROW * MAX_COL))] # Will contain balls

print(MAX_ROW * MAX_COL)
# How many nearby cells should be checked for potential collisions.
REACH =  np.ceil(2* velocity / CELL_SPAN).astype(int)

# Set logo, caption, and surface for the window.
logo = pygame.image.load("logo.png")
pygame.display.set_icon(logo)
pygame.display.set_caption("idealgas.py")
surface = pygame.display.set_mode(size=(width, height))
CHARCOAL = (21,27,31)
GRAY = (40, 44, 52)
BLUE = (97, 175, 239)
PINK = (198, 117, 168)
ORANGE = (203, 123, 70)
RED = (224, 108, 117)
GREEN = (145, 190, 97)
COLORS = [PINK, ORANGE, BLUE, RED, GREEN]
COLORS = [RED]
BACKGROUND = GRAY
BALL_COLOR = BLUE
UNIQUE_BALL = RED

class Ball(object):

    def __init__(self, surface, color, center, radius, mass, velocity,
                 acceleration=0):

        self.surface = surface
        self.color = color # RGB[A] tuple, where A is optional
        self.center = center
        self.radius = radius
        self.mass = mass
        self.velocity = velocity

        self.acceleration = acceleration

    def draw(self):
        return pygame.draw.circle(self.surface, self.color,
                                  self.center.astype(int), self.radius)

    def move(self, dtime=1):
        grid[self.grid_ind()] = None

        self.center += self.velocity * dtime
        self.velocity[1] += self.acceleration * dtime
        grid[self.grid_ind()] = self

    def grid_ind(self):
        ind_x, ind_y = self.xy_ind()
        return MAX_COL * ind_y + ind_x

    def xy_ind(self):
        return (int(self.center[0]/CELL_SPAN),
                        int(self.center[1]/CELL_SPAN))

    def distance(self, ball):
        return (np.linalg.norm(self.center - ball.center)
                - self.radius - ball.radius)

    def is_clear(self, include_center=True):
        for ball in nearest_nbrs(*self.xy_ind(), include_center=True):
            if ball != self and self.distance(ball) <= 0:
                return False

        return True


def exit_time(pos, velocity, radius, width, height):
    """Calculate the time it takes to move out of a box"""
    x, y = pos
    v_x, v_y = velocity

    if v_x > 0:
        dt_x = (width - x - radius) / v_x
    elif v_x < 0:
        dt_x = (radius - x) / v_x
    else:
        dt_x = np.inf

    if v_y > 0:
        dt_y = (height - y - radius) / v_y
    elif v_y < 0:
        dt_y = (radius- y) / v_y
    else:
        dt_y = np.inf

    return dt_x, dt_y



def nearest_nbrs(col, row, reach=3, include_center=False):
    """ List the content of grid withing reach of position (col, row)"""
    col_start, col_end = max(0, col-reach), min(col+reach, MAX_COL-1)
    row_start, row_end = max(0, row-reach), min(row+reach, MAX_ROW-1)

    if include_center == True:
        return [grid[c + MAX_COL * r]
                for c in range(col_start, col_end+1)
                for r in range(row_start, row_end+1)
                if grid[c + MAX_COL * r] != None]
    else:
        return [grid[c + MAX_COL * r]
                for c in range(col_start, col_end+1)
                for r in range(row_start, row_end+1)
                if grid[c + MAX_COL * r] != None and
                   ( c != col or r != row)]


def next_ball_collision(balls, brute_force=False):
    """ Search which balls will collide first."""
    min_t = np.inf
    closest = None
    if brute_force is True:
        for b1, b2 in combinations(balls,2):
            t = ball_collision_time(b1, b2)
            if t <= min_t:
                min_t = t
                closest = (b1, b2)
    else:
        tested = set()
        for b1 in balls:
            for b2 in nearest_nbrs(*b1.xy_ind(), REACH):
                pair = frozenset([b1, b2])
                if  pair not in tested:
                    tested.add(pair)
                    t = ball_collision_time(b1, b2)
                    if t <= min_t:
                        min_t = t
                        closest = (b1, b2)

    return min_t, closest


def next_wall_collision(balls):
    """Search which ball will collide with a wall first"""
    min_t = np.inf
    closest = None
    for b in balls:
        t_x, t_y = exit_time(b.center, b.velocity, b.radius, width, height)
        t = min(t_x,t_y)
        if t <= min_t:
            min_t = t
            if t_x < t_y:
                str = 'x-wall'
            elif t_x > t_y:
                str = 'y-wall'
            else:
                str = 'both'

            closest = (b, str)

    return min_t, closest

def ball_collision_time(ball1, ball2):
    """When will the balls collide. Returns np.inf if they do not."""

    r_1 = ball1.radius
    r_2 = ball2.radius
    dv = ball1.velocity - ball2.velocity
    dq = ball1.center - ball2.center
    #a = (vx_1 - vx_2)**2 + (vy_1 - vy_2)**2
    a = np.dot(dv, dv)
    #b = 2*((vx_1- vx_2)*(x_1 - x_2) + (vy_1- vy_2)*(y_1 -y_2))
    b = 2 * np.dot(dv, dq)
    #c = -(r_1 + r_2)**2 + (x_1 - x_2)**2 + (y_1 - y_2)**2
    c = np.dot(dq, dq) - (r_1 + r_2)**2
    D = -4*a*c + b**2

    if D >= 0 and a != 0: # Collision point exists
        t1 = -(b + np.sqrt(D))/(2*a) # Time of the first collision
        t2 = -(b - np.sqrt(D))/(2*a)
        if t1 >= 0: # Collision in the future
            return t1
        elif t2 >=0 and b < 0:
            return 0 # Balls overlap and still getting closer
        else:   # collision in the past
            return np.inf
    else:
        return np.inf


def spawn_test(surface):
    """ Spawn two balls for testing the engine """

    ball1 = Ball(surface=surface,
                      color=UNIQUE_BALL,
                      center=np.array([width/10+2*radius, height/2]),
                      radius=radius,
                      mass=100000,
                      velocity=np.array([-CELL_SPAN, 0])/10)
    particles.append(ball1)
    grid[ball1.grid_ind()] = ball1
    ball2 = Ball(surface=surface,
                      color=BALL_COLOR,
                      center=np.array([width/10-2*radius, height/2]),
                      radius=radius,
                      mass=1,
                      velocity=np.array([-CELL_SPAN, 0])/10)
    particles.append(ball2)
    grid[ball2.grid_ind()] = ball2


def spawn_particles(surface, color, mass, velocity, acceleration, max_attepts):

    spawns = [] # list of balls that act as swawning locations
    first_ball = Ball(surface=surface,
                      color=UNIQUE_BALL,
                      center=np.array([width/2, height/2]),
                      radius=radius,
                      mass=mass,
                      velocity=(np.random.rand(2) - 0.5) * velocity,
                      acceleration=acceleration)
    particles.append(first_ball)
    spawns.append(first_ball)
    grid[first_ball.grid_ind()] = first_ball

    while len(particles) < particle_count and len(spawns) != 0:

        spawn = np.random.choice(spawns)

        attempts = 0
        while attempts < max_attepts:
            # pick a random position near the spawn
            angle = np.random.rand()*2*np.pi
            r = (np.random.rand()+1) * 2*radius
            rnd_r = (spawn.center
                    + r * np.array([np.cos(angle), np.sin(angle)]))
            rnd_r = np.minimum(rnd_r, np.array([width,height])-radius - 1 )
            #rnd_r = np.maximum(rnd_r, radius + 1)

            new_ball = Ball(surface=surface,
                            color=random.choice(COLORS), #color,
                            center=rnd_r,
                            radius=radius,
                            mass=mass,
                            velocity=np.array(
                                     [np.cos(angle), np.sin(angle)]) * velocity,
                            acceleration=acceleration)
                            #velocity= velocity.copy())

            if new_ball.is_clear():
                spawns.append(new_ball)
                particles.append(new_ball)
                grid[new_ball.grid_ind()] = new_ball
                break

            elif attempts == max_attepts-1:
                spawns.remove(spawn)
                attempts += 1
            else:
                attempts += 1

def collide(b1, b2):

    c = (b2.center - b1.center);
    c = c / np.linalg.norm(c)
    delta_v = b2.velocity - b1.velocity
    M = b1.mass + b2.mass
    delta_projected = 2 / M * np.dot(delta_v, c) * c
    b1.velocity, b2.velocity = (b1.velocity + b2.mass * delta_projected,
                                b2.velocity - b1.mass * delta_projected)

def collide_with_balls(balls):
    for b1, b2 in balls:
        collide(b1,b2)

def collide_with_walls(balls):
    b, wall = balls
    if wall == 'x-wall':
        b.velocity[0] *= -1
    elif wall == 'y-wall':
        b.velocity[1] *= -1
    else:
        b.velocity *= -1


def draw_all(particles):

    global old_rects

    if old_rects == []:
        surface.fill(BACKGROUND)
        pygame.display.update()
    else:
        for rect in old_rects:
            surface.fill(BACKGROUND,rect)

    new_rects = [particle.draw() for particle in particles]
    pygame.display.update(old_rects + new_rects)
    old_rects = new_rects


def move_all(paricles, dtime=1):

    for particle in particles:
        particle.move(dtime)


def handle_events():

    for event in pygame.event.get():
        #print(event)
        if event.type == pygame.QUIT:
            return False
        elif event.type == pygame.KEYDOWN and event.key == pygame.K_q:
            return False

    return True


def main():

    spawn_particles(surface,
                      color=BALL_COLOR,
                      mass=1,
                      velocity=velocity,
                      acceleration=acceleration,
                      max_attepts=500)
    #spawn_test(surface)

    mean_time = 0

    running = True
    while running:
        # start_time = time.time()
        clock.tick(40) # Sets the maximum framerate
        # total_time = time.time()-start_time
        # mean_time = (mean_time * 10 + total_time) / 11
        # print(f"Time delta: {mean_time}")

        running = handle_events()

        dtime = 1
        while dtime > 0: # 20 ms / 5.7 ms  = 3.5 times @ 200 particles 720 p

            tw, closest = next_wall_collision(particles) # 0.67 ms
            tc, pairs = next_ball_collision(particles, brute_force=False) # 4.7 ms
            dt = min ([dtime, tw, tc])
            dtime -= dt

            move_all(particles, dt) # rest 1.3 ms
            if dt == tc:
                collide(*pairs)
            elif dt == tw:
                collide_with_walls(closest)

        draw_all(particles) # 18 ms at 200 particles


if __name__ == "__main__":
    main()
