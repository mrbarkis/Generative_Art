import random
import os
from math import sqrt, exp, pi, atan, tan, acos

####    PICK A SEED    ####
# My favorite seeds: 406, 302, 248, 373, 347#287#
SEED = 1649#None#302#302# None #248 #302  # Set to None, if you want a random one
if SEED == None:
    SEED = random.randint(42, 4200)
random.seed(SEED)
print("This tree of gauss grew from a seed {}".format(SEED))

### BRANCH PARAMETERS ###
IS_UNIFORM = False
CONSTANT_FORK = True
TREE_OFFSET = -20 # from the bottom of the window
TRUNK_HEIGHT = 100
THICKNESS_PER_LEVEL = 10

BRANCH_SHRINKAGE = 0.82#0.85  # Rate at which the braches get shorter
MAX_LEVEL = 6 # Levels of branches
TRUNK_THICKNESS = THICKNESS_PER_LEVEL * MAX_LEVEL
N_BRANCHES = 4
N_SAMPLES = N_BRANCHES ** MAX_LEVEL / 5
VAR = 5000  # A common variance for the branches

### LEAF PARAMETERS ###
LEAF_WIDTH = 6
LEAF_LENGTH =8
BIN_WIDTH =  40#60 #
LEAF_TERMINAL_VELOCITY_RANGE = (3, 5)
LEAF_SPIN_RATE_RANGE = (2 * pi / 150, 2 * pi / 250)
LEAF_SPIN_RADIUS_RANGE = (20, 30)


### HILL ###
HILL_OFFSET = -80
HILL_HEIGHT = 35
PILE_HEIGHT = 35
MAX_ZOOM = 5

### NOISE PARAMETERS ###
GUST_SPEED = 1.0/120
BRANCH_RESPONSE = 5
LEAF_RESPONSE = 2 * pi / 15
SHEEN_RESPONSE = 0.8
NOISE_BIAS = 0.3 # Wind noise seems to have range from 0 to 0.7, use this to center it
noiseDetail(8, 0.2)

###    COLORS    ####
RED, GREEN, BLUE, = [114, 0, 16], [73, 114, 0], [0, 73, 114]
OLIVE, BROWN, WHITE = [114, 98, 0], [114, 41, 0], [255, 255, 255]
AUTUM1, AUTUM2, AUTUM3 = [156, 39, 6], [212, 91, 18], [243, 188, 46]

def mix_colors(c1, c2, ratio):
    return [(1 - ratio) * a + ratio * b for a, b in zip(c1, c2)]
BG_COLOR = mix_colors(BLUE, WHITE, 1)
LEAF_COLOR = BLUE
LEAF_SAMPLE_COLOR = BLUE
LEAF_AUTUM_COLORS = [RED]#, AUTUM3]
TRUNK_COLOR  = mix_colors(BROWN, WHITE, 0.4)
HILL_COLOR = mix_colors(LEAF_SAMPLE_COLOR, WHITE, 0.7)
HILL_COLOR_2 = mix_colors(LEAF_AUTUM_COLORS[0], WHITE, 0.7)
GROUND_COLOR = mix_colors(BROWN, WHITE, 1)
print(BG_COLOR)
print(BLUE)

### ANIMATION ###
RECORD = False
FOLDER = os.path.join('.','saved_frames\\')
FRAMES_PER_LEVEL = 3 * 30 # how fast the tree grows
FRAMES_BEFORE_AUTUM = 3 * 30 # frames at full hight before leaves start changing color

MAX_DELAY_TO_AUTUM = 5 * 30 # for an individual leaf
FRAMES_TO_CHANGE_COLOR = 3 * 30 # for an individual leaf
MAX_DELAY_TO_FALL = 30 * 30 # how long can an autum leaf hang on
DELAY_AUTUM_NOT_SAMPLED = 0 * 30
DELAY_FALLING_NOT_SAMPLED = 3 * 30
STOP_DELAY = 20 * 30 
#DELAY_NOT_SAMPLED = 0 * 30

### NO MORE OPTIONS BELOW, JUST INITIALIZATION STUFF ###
N_LEAVES = N_BRANCHES ** MAX_LEVEL
STD = sqrt(VAR * MAX_LEVEL)  # Standard deviation of the hill, i.e. expected std of the leaf pile
print("Standard deviation of a branch is {:.2f}".format(sqrt(VAR)))
print("Standard deviation of leaves after {} branhings is  sqrt({}) times that, i.e. {:.2f}".format(MAX_LEVEL, MAX_LEVEL, STD))
PILE_SCALE = sqrt(2 * pi ) * STD / BIN_WIDTH * float(PILE_HEIGHT) / N_LEAVES # Scale the heigth of the leaf pile




### PRE-RANDOMIZE EVERYTHING  ###

print(N_LEAVES)
autum_starts = [(MAX_LEVEL + 1) * FRAMES_PER_LEVEL + FRAMES_BEFORE_AUTUM
                + random.randint(0, MAX_DELAY_TO_AUTUM) for _ in range(N_LEAVES)]
falling_starts = [kf_autum + FRAMES_TO_CHANGE_COLOR 
                  + random.randint(0, MAX_DELAY_TO_FALL) for kf_autum in autum_starts]
v_drop = [random.randint(*LEAF_TERMINAL_VELOCITY_RANGE) for _ in range(N_LEAVES)]
w_drop = [random.uniform(*LEAF_SPIN_RATE_RANGE) for _ in range(N_LEAVES)]
a_drop = [random.randint(*LEAF_SPIN_RADIUS_RANGE) for _ in range(N_LEAVES)]
leaf_colors = [random.choice(LEAF_AUTUM_COLORS) for _ in range(N_LEAVES)]



#### GLOBAL COUNTERS, CACHES, AND CONTAINERS ####
frame = 0 # Current frame

# As the branches will be drawn recursively (a branch at a time instead of level at a time),
# the nth branch drawn changes when ever we go deeper in the recursion. This causes the nth call to a
# random generator to happen at different branching points. We can use a key
branch_id = [0] * (MAX_LEVEL+1) # to store unique tuple: (current level, nt_branch on the level)
# and a dictionary
branch_cache = {} # to store the randomization.

leaf_xys = [(0,0)] * N_LEAVES # will contain coordinate tuples for the leaves.
leaf_angles = [0] * N_LEAVES # will contain the angles of the leaves.
leaf_stops = [10000 * 30] * N_LEAVES # will contain the frames at which leaves hit the ground.
leaf_detach_xys = [(0,0)] * N_LEAVES
leaf_p = [(0)] * N_LEAVES
leaf_p_max = [0] * (MAX_LEVEL + 1)
# leaf_bins that contain the leaf pile has to be declared in the setup() function

####    HELP FUNCTIONS    ####
def sample_leafs(probabilities, N=N_SAMPLES):
    global falling_starts, autum_starts, leaf_colors
    assert len(probabilities) >= N
    
    # Iterators that give random autum and detach keyframes in ascending order
    kf_autum = ( a for _, a in sorted((zip(falling_starts, autum_starts))))
    kf_detach = iter(sorted(falling_starts))
    
    # Divide leaves into sampled and not sampled sets
    sample_inds = sample(list(cdf(probabilities)), N=N)
    not_sampled_inds = set(range(len(probabilities))) - sample_inds
    
    
    # The sampled leaves drop first and have a different color
    for i in sample_inds:#sorted(sample_inds, key=lambda ind: probabilities[ind], reverse=True):
        leaf_colors[i] = LEAF_SAMPLE_COLOR
        autum_starts[i] = next(kf_autum)
        falling_starts[i] = next(kf_detach)

    for i in not_sampled_inds:#sorted(not_sampled_inds, key=lambda ind: probabilities[ind], reverse=True):
        autum_starts[i] = next(kf_autum) + DELAY_AUTUM_NOT_SAMPLED
        falling_starts[i] = next(kf_detach) + DELAY_AUTUM_NOT_SAMPLED + DELAY_FALLING_NOT_SAMPLED
    
    
def select_element(xs, ind):
    return xs[ind]
    

def cdf(probabilities):
    '''Cumulative distribution'''
    sum = 0
    for p in probabilities:
        sum += p
        yield sum

def find_supremum_index(y, bound):
    ''' Finds the smallest index i, s.t.  bound < y[i]
    assumes that y is monotonically increasing list, such as a cumulative distirbution'''
    
    if  bound <= y[0]:
        return 0
       
    i_inf, i_sup = 0, len(y) - 1 
    i = (i_sup - i_inf) // 2 
    while (i_sup - i_inf) // 2 > 0:
        if y[i] >= bound:
            i_sup = i
            i = i_sup - (i_sup - i_inf) // 2 
        else:
            i_inf = i
            i = i_inf + (i_sup - i_inf) // 2 
        
    return i_sup

def sample(cdf, N=1):
    assert len(cdf) >= N
    s = set()
    while len(s) != N:
        s.add(find_supremum_index(cdf, bound=random.random()))
        
    return s
    
    
def drawCurve(xs, ys, col=[114, 41, 0], filled=True):
    if filled:
        fill(*col)
    else:
        noFill()
    strokeWeight(2)
    stroke(*col)
    beginShape()
    curveVertex(xs[0], ys[0])  # the first control point
    for x, y in zip(xs, ys):
        curveVertex(x, y)
    curveVertex(xs[-1], ys[-1])  # the last control point
    endShape()


def generatePMF(N, var=100, is_uniform=IS_UNIFORM):
    """Greates a random probability mass function (PMF) with
    specific variance and zero mean."""
    
    # Random weighted shifts ( shift x times porbability p)
    xps = [random.random() for _ in range(N)]
    mean = sum(xps) / N
    xps = [xp - mean for xp in xps]
    
    if is_uniform:
        ps = [1.0/N] * N
    else:# Random porabilities
        ps = [random.random() for _ in range(N)]
        total = sum(ps)
        ps = [p / total for p in ps]  # Normalize probabilites

    # Solve the shifts
    xs = [ xp / p for xp, p in zip(xps, ps)]

    # Scale the variance
    current_variance = sum([p * x ** 2 for x, p in zip(xs, ps)])
    scale_factor = sqrt(var / current_variance)
    xs = [x * scale_factor for x in xs]
    
    # Sort to make left most branch start from left
    ps = [ p for _, p in sorted(zip(xs, ps))]
    xs = sorted(xs)
    
    return xs, ps


def gaussian(x, mu, sigma):
    return exp(-0.5 * (x - mu) ** 2 / sigma ** 2)  # / (sigma * sqrt(2 * pi))


def drawLeaf(col=LEAF_COLOR, leaf_length=LEAF_LENGTH, leaf_width=LEAF_WIDTH, sheen=0, leaf_alpha=255):
    l, w = leaf_length, leaf_width
    noStroke()
    # One half
    c = mix_colors(col, WHITE, sheen)
    fill(c[0], c[1], c[2], leaf_alpha)
    beginShape()
    vertex(0, 0)
    bezierVertex(-w, 0.1 * l, -w * 0.5, -l, 0, -l)
    endShape()
    # The other half
    c = mix_colors(col, WHITE, sheen * 0.8)
    fill(c[0], c[1], c[2], leaf_alpha)
    beginShape()
    vertex(0, 0)
    bezierVertex(w, 0.1 * l, w * 0.5, -l, 0, -l)
    endShape()


def branch(branch_length, branch_progress, x, p, level, leaf_coord, thickness, previous_point=(0, 0), depth_id=0, draw_leaves=True):
    global branch_id, leaf_id, leaf_xys, leaf_angle, frame, leaf_p_max
    wind_noise = noise(leaf_coord[0] + GUST_SPEED * frame, leaf_coord[1]) - NOISE_BIAS
    if level == 1: # The growing branch has not reached its full extent
        #branch_progres = max(0.10, branch_progress)

        x = [branch_progress * value for value in x]
        branch_length = branch_progress * branch_length
        
        
    if level > 0:
        shift, p_previous = 0, -1 # used to align new branches side by side.
        for xi, pi in zip(x, p):
            with pushMatrix():
                branch_id[depth_id] += 1
                fork_id = (depth_id, branch_id[depth_id])
                
                # rnd forks
                if fork_id in branch_cache:
                    # Load y-shift
                    _ = random.uniform(0.1, 2.5)
                    if CONSTANT_FORK: # load only y shifts
                        rnd = branch_cache[fork_id]
                    else: # load x and y shifts
                        _ , _ = generatePMF(N=N_BRANCHES, var=VAR)
                        rnd, x, p = branch_cache[fork_id]        
                else:
                    # Save y-shift
                    rnd = random.uniform(0.1, 2.5)
                    if CONSTANT_FORK:
                        branch_cache[fork_id] = rnd
                    else:
                        x, p = generatePMF(N=N_BRANCHES, var=VAR)
                        branch_cache[fork_id] = (rnd, x, p)
                
                
                #print("id {} {} {}".format(depth_id, branch_id[depth_id], rnd))
                dy =  rnd * branch_length

                # Thickness equals the probability of the branch:
                strokeWeight(thickness * pi)

                # Shift starting position next to previus branch:
                shift, p_previous = shift + thickness * (pi + p_previous) / 2, pi
                
                xi = xi + sqrt(depth_id) * BRANCH_RESPONSE * wind_noise
                if xi * wind_noise > 0: # Branch bent the same direction as wind
                    dy = dy + sqrt(depth_id) * BRANCH_RESPONSE * wind_noise # push down
                else:
                    dy = dy - sqrt(depth_id) * BRANCH_RESPONSE * wind_noise # push up
                    
                
                 
                # Draw the branch:
                stroke(*TRUNK_COLOR)
                noFill()
                beginShape()
                curveVertex(*previous_point)  # To make smooth bends.
                curveVertex(shift, 0)
                curveVertex(xi, -dy)
                curveVertex(xi, -dy)
                endShape()
                
                
                # Continue the branch:
                translate(xi, -dy)
                branch(
                    branch_length * BRANCH_SHRINKAGE,  # Shorthen the next branch
                    branch_progress,
                    x,
                    p, 
                    level=level - 1,  # End condition for recursion
                    leaf_coord=(leaf_coord[0] + xi, leaf_coord[1] - dy),
                    thickness=thickness * pi,  # Slim down the next branch
                    previous_point=(-xi, +dy),
                    depth_id=depth_id + 1,
                    draw_leaves=draw_leaves
                )
    else:
        total_leaf_probability = thickness / TRUNK_THICKNESS
        if total_leaf_probability > leaf_p_max[depth_id]:
            leaf_p_max[depth_id] = total_leaf_probability
        
        if draw_leaves:
            rotate(-atan(previous_point[1] / (0.1 + previous_point[0])) + LEAF_RESPONSE * wind_noise)
            drawLeaf(sheen=SHEEN_RESPONSE*(wind_noise + NOISE_BIAS))
                     #leaf_alpha = 255 * total_leaf_probability / leaf_p_max[depth_id] )
            
        # Store leaf positions into ends dictionary
        leaf_xys[leaf_id] = leaf_coord
        leaf_p[leaf_id] = total_leaf_probability
        leaf_angles[leaf_id] = -atan(previous_point[1] / (0.1 + previous_point[0]))
        leaf_id += 1
        
    

            
            
def draw_hill(y_offset, h, col=HILL_COLOR, x_offset=0, std=STD):
    mean = width / 2 + x_offset # Hill at the center of the screen.
    xs = range(width) # Hill extents across the entire window
    ys = [height + y_offset - h * gaussian(x, mean, std) for x in xs]
    drawCurve(xs, ys, col=col, filled=True)


def draw_tree(y_offset, branch_length, branch_progress, thickness, col, level, draw_leaves=True):
    # Generate rnd x-shifts for the branches
    branch_x, branch_p = generatePMF(N=N_BRANCHES, var=VAR)
    # Start the three from the bottom of the screen
    translate(width / 2, height + y_offset)
    # Draw the "root" stem
    strokeWeight(thickness)
    stroke(*col)
    strokeCap(SQUARE) # PROJECT, SQUARE, ROUND
    if level == 0: # The truck is growing
        branch_length = branch_length * branch_progress
    
    line(0, 0, 0, -branch_length)
    translate(0, -branch_length)
    branch(
        branch_length,
        branch_progress,
        branch_x,
        branch_p,
        level=level,
        leaf_coord=(width / 2, height + y_offset),
        thickness=thickness,
        previous_point=(0, branch_length),
        draw_leaves=draw_leaves
    )            
            
            

def setup():
    #size(400, 600)
    #size(1280, 720)
    size(1920, 1080)
    fullScreen()
    frameRate(30)
    
    
    
    global leaf_bins, leaves_x_offset, leaves_std, STOP_FRAME
    leaf_bins = [0] * (width // BIN_WIDTH)
    STOP_FRAME = max(falling_starts) + height / LEAF_TERMINAL_VELOCITY_RANGE[0] + STOP_DELAY

    # Calculate the position and variance of the unsampled leaves
    if CONSTANT_FORK: # then generate it and calculate the metrics
        random.seed(SEED)
        noiseSeed(SEED)
        branch_x, branch_p = generatePMF(N=N_BRANCHES, var=VAR)
        x_ave, x2_ave = (sum(branch_x)/N_BRANCHES,
                         sum([x**2 for x in branch_x])/N_BRANCHES)
        leaves_x_offset, leaves_std = (x_ave * MAX_LEVEL,
                                       sqrt((x2_ave - x_ave**2) ))
    else: # everything will average to zero
        leaves_x_offset, leaves_std = 0, STD / sqrt( MAX_LEVEL)
            
    
def draw():
    global frame, branch_id, leaf_id, leaf_xys, leaf_stops, RECORD
    frame += 1
    if frame == STOP_FRAME:
        print('stopped')
        RECORD = False
        exit()
    
    
    # FRESH CANVAS FOR EACH FRAME
    random.seed(SEED) # A fresh start every frame
    noiseSeed(SEED)
    branch_id = [0] * (MAX_LEVEL + 1) # trunk is the 0th level
    leaf_p_max = [0] * (MAX_LEVEL + 1) 
    leaf_id = 0

    
    # DRAW BACKGROUND
    
    zoom = min(sqrt(MAX_LEVEL+1), sqrt(float(frame) / FRAMES_PER_LEVEL)) # 0 to max level
    background(*BG_COLOR)
    fill(*GROUND_COLOR)
    rect(0, height + HILL_OFFSET* (1 + zoom), width, -HILL_OFFSET* (1 + zoom))
    
    draw_hill(y_offset=HILL_OFFSET * (1 + zoom),
              h=HILL_HEIGHT* zoom,
              x_offset=leaves_x_offset,
              std=leaves_std * zoom,
              col=HILL_COLOR_2
              )
    draw_hill(y_offset=HILL_OFFSET * (1 + zoom),
              h=HILL_HEIGHT * float(N_SAMPLES) / N_LEAVES* zoom,
              std=float(STD) / sqrt(MAX_LEVEL) * zoom)
    

    if frame < FRAMES_PER_LEVEL * (1 + MAX_LEVEL): # Growing the tree
        draw_tree(
            y_offset=TREE_OFFSET * (MAX_ZOOM - zoom),
            branch_length=TRUNK_HEIGHT,
            branch_progress=float(frame % FRAMES_PER_LEVEL) / (FRAMES_PER_LEVEL) ,
            #thickness=min(TRUNK_THICKNESS, floor(TRUNK_THICKNESS * growth_progress)),
            thickness=TRUNK_THICKNESS * float(frame) / (FRAMES_PER_LEVEL * (1 + MAX_LEVEL)),
            col=TRUNK_COLOR,
            level=frame // FRAMES_PER_LEVEL,
            draw_leaves = False
            )
        
        for ((x, y), angle) in zip(leaf_xys, leaf_angles):
            with pushMatrix():
                translate(x - width / 2, y - height - TREE_OFFSET* (MAX_ZOOM - zoom))
                rotate(angle + LEAF_RESPONSE * noise(x + GUST_SPEED *frame, y))
                drawLeaf(sheen= SHEEN_RESPONSE * noise(x + GUST_SPEED *frame,y))
            
        
                                                                
    elif frame == FRAMES_PER_LEVEL * (1 + MAX_LEVEL):
        draw_tree(
            y_offset=TREE_OFFSET *  (MAX_ZOOM - zoom),
            branch_length=TRUNK_HEIGHT,
            branch_progress=1 ,
            thickness=TRUNK_THICKNESS,
            col=TRUNK_COLOR,
            level=MAX_LEVEL,
            draw_leaves=False
            )
        for ((x, y), angle) in zip(leaf_xys, leaf_angles):
            with pushMatrix():
                translate(x - width / 2, y - height - TREE_OFFSET* (MAX_ZOOM - zoom))
                rotate(angle + LEAF_RESPONSE * noise(x + GUST_SPEED *frame, y))
                drawLeaf(sheen= SHEEN_RESPONSE * noise(x + GUST_SPEED *frame,y))
                
        sample_leafs(leaf_p, N=N_SAMPLES)
 
        
        
    else:
        # First the fullgrown tree without its leaves
        draw_tree(
            y_offset=TREE_OFFSET*  (MAX_ZOOM - zoom),
            branch_length=TRUNK_HEIGHT,
            branch_progress=1 ,
            thickness=TRUNK_THICKNESS,
            col=TRUNK_COLOR,
            level=MAX_LEVEL,
            draw_leaves=False
            )
        

        
        for i, ((x, y), angle, col, kf_autum, kf_detach, v, w, a) in enumerate(zip(
                                                                leaf_xys,
                                                                leaf_angles,
                                                                leaf_colors,
                                                                autum_starts,
                                                                falling_starts,
                                                                v_drop, w_drop, a_drop)):
            
            #print(sum(leaf_p))
            #print(sum([p * c[0] for c, p in zip(leaf_xys,leaf_p)]) - width / 2)
            with pushMatrix():
                
                width_factor = 1
                
                # Determine color for each leaf at this frame
                if frame < kf_autum:
                    current_color = LEAF_COLOR
                elif frame < kf_autum + FRAMES_TO_CHANGE_COLOR:
                    current_color = mix_colors(LEAF_COLOR, col, float(frame - kf_autum) / FRAMES_TO_CHANGE_COLOR)
                else:
                    current_color = col
                
                   
                # Determine y and x postions of each leaf at this frame
                if frame == kf_detach: # store where the leaf detaches
                    leaf_detach_xys[i] = (x, y)
                    
                if frame > leaf_stops[i]: # Has hit the pile
                    x, y = leaf_detach_xys[i]
                    dy = v * (leaf_stops[i] - kf_detach) # total shift from tree to pile
                    x += -a * 0.5 * (sin(w * dy + acos(tan(2*pi/360*angle))) - sin( acos(tan(2*pi/360*angle))))
                    y += dy
                elif frame > kf_detach: # Is falling
                    x, y = leaf_detach_xys[i]
                    dy = v*(frame - kf_detach)
                    width_factor = 0.5 * (sin(w * dy + acos(tan(2*pi/360*angle))) + sin( acos(tan(2*pi/360*angle))))
                    
                    angle =  atan(a * w * cos(w *dy + acos(tan(2*pi/360*angle))))
                    y += dy
                    bin_nbr = int(min(max(0, x // BIN_WIDTH), len(leaf_bins)-1)) # truncate leafs that exit the window
                    x += -a * 0.5 * (sin(w * dy + acos(tan(2*pi/360*angle))) - sin( acos(tan(2*pi/360*angle))))
                    if y - TRUNK_HEIGHT - TREE_OFFSET* (MAX_ZOOM - zoom) -float(TRUNK_THICKNESS) / 2 > height - leaf_bins[bin_nbr]: # Record when it hits the pile of leaves
                        bin_on_left, bin_on_right = (int(max(0, bin_nbr - 1)),
                                                     int(min(len(leaf_bins)-1, bin_nbr + 1)))
                        stack_increment = min(float(LEAF_WIDTH)/2, PILE_SCALE * zoom)
                        leaf_bins[bin_nbr] +=  stack_increment / 2
                        leaf_bins[bin_on_left] += stack_increment / 4
                        leaf_bins[bin_on_right] += stack_increment / 4
                        
                        leaf_stops[i] = frame
                        
                   

                    
                translate(x - width / 2, y - height - TREE_OFFSET* (MAX_ZOOM - zoom))
                rotate(angle + LEAF_RESPONSE * noise(x + GUST_SPEED *frame, y))
                drawLeaf(current_color,
                         leaf_width=LEAF_WIDTH *  width_factor,
                         leaf_length=LEAF_LENGTH,
                         sheen= SHEEN_RESPONSE * noise(x + GUST_SPEED *frame,y),
                         #leaf_alpha= 255 * leaf_p[i] / max(leaf_p)
                         )
                
                
    if RECORD:
        saveFrame("{}####.png".format(FOLDER))
        

                
                
