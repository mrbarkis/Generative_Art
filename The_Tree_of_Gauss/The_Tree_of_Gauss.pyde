add_library("pdf")  # modify setup() to print a pdf
import random
from math import sqrt, exp, pi, atan


####    PICK A SEED    ####
# My favorite seeds: 201, 403, 223, 67, 160
seed = 160  # Set to None, if you want a random one
if seed == None:
    seed = random.randint(42, 420)
random.seed(seed)
print("This tree of gauss grew from a seed {}".format(seed))

### PICK BRANCH PARAMETERS
l = 63  # Height of the first branch, i.e. of the trunk
l_rate = 0.93  # Rate at which the braches get shorter
thickness = 40  # Thickness of the first branch, i.e. of the trunk
level = 6  # Levels of branches
N_branches = 3
var = 1500  # Variance (in x direction) of each branching)
std = sqrt(var * level)  # of the hill
print("Standard deviation of a branch is {}".format(sqrt(var)))
print("Standard deviation of the tree {}".format(sqrt(var * level)))


###    PICK COLORS    ####
def mix_colors(c1, c2, ratio):
    return [(1 - ratio) * a + ratio * b for a, b in zip(c1, c2)]


blue_c = [0, 73, 114]
green_c = [73, 114, 0]
olive_c = [114, 98, 0]
brown_c = [114, 41, 0]
red_c = [114, 0, 16]
leaf_c = blue_c  # green_c #blue_c #olive_c
autum1_c = olive_c
autum2_c = red_c
white = [200, 200, 200]
trunc_c = mix_colors(brown_c, white, 0.4)
hill_c = mix_colors(brown_c, white, 0.95)

####    HELP FUNCTIONS    ####
ends = {}  # aux variable to which to store the leaf positions


def drawCurve(xs, ys, col=[114, 41, 0], filled=True):
    if filled:
        fill(*col)
    else:
        noFill()
    stroke(*col)

    beginShape()
    curveVertex(xs[0], ys[0])  # the first control point
    for x, y in zip(xs, ys):
        curveVertex(x, y)
    curveVertex(xs[-1], ys[-1])  # the last control point
    endShape()


def generatePMF(N, var=100):
    """Greates a random probability mass function (PMF) with
    spesific variance and zero mean."""
    x = random.sample(range(101), N)  # N discrete values.
    mean = sum(x) / len(x)
    x = [xi - mean for xi in x]  # Shift mean to zero
    x = sort(x)
    # Generate a random probabilities for the N_branches:
    p = [random.random() for _ in x]
    total = sum(p)
    p = [pi / total for pi in p]  # Normalize probabilites
    # Scale the variance
    current_variance = sum([pi * xi ** 2 for xi, pi in zip(x, p)])
    ratio = sqrt(var / current_variance)
    x = [xi * ratio for xi in x]
    return x, p


def gaussian(x, mu, sigma):
    return exp(-0.5 * (x - mu) ** 2 / sigma ** 2)  # / (sigma * sqrt(2 * pi))


def drawLeaf(col, l, w):
    noStroke()
    rnd = random.random()
    fill(*mix_colors(col, white, rnd))
    beginShape()
    vertex(0, 0)
    bezierVertex(-w, 0.1 * l, -w * 0.5, -l, 0, -l)
    endShape()
    fill(*mix_colors(col, white, rnd * 0.8))
    beginShape()
    vertex(0, 0)
    bezierVertex(w, 0.1 * l, w * 0.5, -l, 0, -l)
    endShape()


def branch(b_length, x, p, level, end_position, thickness, previous_point=(0, 0)):
    if level > 0:  # exit condition for recursion
        shift, p_previous = 0, -1
        for xi, pi in zip(x, p):
            with pushMatrix():
                # Decide y shift for artistic purposes (so that branches won't overlap)
                dy = random.uniform(0.3, 2.5) * b_length

                # Thickness equals the probability of the branch:
                strokeWeight(thickness * pi)  # Reduce the thickess of branches

                # Shift starting position next to previus branch:
                shift, p_previous = shift + thickness * (pi + p_previous) / 2, pi

                # Draw the branch:
                stroke(*trunc_c)
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
                    b_length * l_rate,  # Shorthen the next branch
                    x,
                    p,  # Use the same brach probabilities or generate new ones
                    level=level - 1,  # End condition for recursion
                    end_position=end_position + xi,  # Track branch position
                    thickness=thickness * pi,  # Slim down the next branch
                    previous_point=(-xi, +dy),
                )
    else:
        # Draw a leaf
        rotate(-atan(previous_point[1] / (0.1 + previous_point[0])))
        drawLeaf(leaf_c, l=15, w=10)
        # Store leaf positions into ends dictionary
        end_position = round(end_position - width / 2)
        if end_position in ends:
            ends[end_position] += 1
        else:
            ends[end_position] = 1


def setup():

    size(595, 842)  # A1 poster in mm
    #size(595, 842, PDF, "tree.pdf")
    noLoop()
    # fullScreen()


def draw():
    background(255)

    ### Draw a bell curve ###
    stroke(*blue_c)
    mean = width / 2  # Hill at the center of the screen.
    xs = range(width)
    ys = [height - 90 - 50 * gaussian(x, mean, std) for x in xs]
    drawCurve(xs, ys, col=hill_c, filled=True)

    ### Draw the tree ###
    branch_x, branch_p = generatePMF(N=N_branches, var=var)
    print("Horizontal branch shifts: {}".format(branch_x))
    print("Branch probabilites: {}".format(branch_p))

    stroke(*trunc_c)
    # Start the three from the bottom of the screen
    translate(width / 2, height - 30)
    # Draw the "root" stem
    strokeWeight(thickness)
    strokeCap(SQUARE) # affects how the branches are connected
    # PROJECT, SQUARE, ROUND
    line(0, 0, 0, -l)
    translate(0, -l)
    branch(
        l,
        branch_x,
        branch_p,
        level=level,
        end_position=width / 2,
        thickness=thickness,
        previous_point=(0, l),
    )

    drawLeafHist(ends)


def drawLeafHist(dict):
    bin_max_height = 50  # Height of the highest leaf bin.
    bin_width_p = 0.2  # bin_width = bin_height * bin_width_p
    translate(0, l + 30)
    max_count = max(dict.values())
    stroke(200, 0, 0)
    strokeWeight(10)
    for x, bin_count in dict.items():
        bin_height = bin_max_height * bin_count / max_count
        bin_width = bin_height * bin_width_p 
        for _ in range(bin_count):  # Fill a bin
            with pushMatrix():
                translate(
                    x + (random.random() - 0.5) * bin_width,
                    -random.random() * bin_height,
                )
                rotate(random.random() * 2 * pi)
                autum_c = red_c  # mix_colors(autum1_c, autum2_c, random.random())
                drawLeaf(autum_c, l=15, w=10)
