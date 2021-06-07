add_library('pdf')
import pickle
from math import sqrt, cos
import os

###    PICK COLORS    ####
def mix_colors(c1, c2, ratio):
    return [(1 - ratio) * a + ratio * b for a, b in zip(c1, c2)]

gray = [200, 200, 200]
blue_c = [0, 73, 114]
green_c = [73, 114, 0]
olive_c = [114, 98, 0]
brown_c = [114, 41, 0]
red_c = [114, 0, 16]
gray_c = mix_colors(brown_c, gray, 0.95) 
node_c = [255, 255, 255, 255] 
path_c = blue_c + [100] #mix_colors(blue_c, gray, 0.01) + [100]
fill_c = [255, 255, 255, 0]


def drawPoints(xs, ys, radius=10, col=[200, 100, 100], filled=True):
    if filled:
        fill(*col)
    else:
        noFill()
        
    noStroke()
    for x, y in zip(xs, ys):
        circle(x, y, radius)

def drawCurve(xs, ys, col=[100, 100, 200], thickness=3, filled=True):
    if filled:
        fill(*fill_c)
    else:
        noFill()
    
    stroke(*col)
    strokeWeight(thickness)

    beginShape()
    curveVertex(xs[0], ys[0])  # the first control point
    for x, y in zip(xs, ys):
        curveVertex(x, y)
    curveVertex(xs[0], ys[0])  
    curveVertex(xs[1], ys[1])  # the last control point
    endShape()



def loadData():
    
    path = pickle.load(open(os.path.join('..', 'data', 'best_route_8586206.p'), "rb" ) )
    #path = pickle.load(open(os.path.join('..', 'data', 'best_route_8545015.p'), "rb" ) )

    xs, ys =[], []
    for node in path:
        xs.append(node[0])
        ys.append(node[1])

    return xs, ys

def average(xs):
    return sum(xs) / len(xs)

def variance(xs):
    mu = average(xs)
    x_2 = sum([x ** 2 for x in xs])
    return x_2 / len(xs) - mu ** 2
    
def normalize(xs):
    mu = average(xs)
    std = sqrt(variance(xs))
    return [(x - mu) / std for x in xs]
    
def scale_and_move(xs, k, b):
    return [ k * x + b for x in xs]

def draw_background():
    img_path = os.path.join('..', 'media', 'map.png')
    image(loadImage(img_path),0, 0, 800, 1200)


def setup():

    size(800, 1200)
    #size(800, 1200, PDF, "trip_4by6.pdf")
    noLoop()
    # fullScreen()

def draw():
    background(255)
    draw_background()
    xs, ys = loadData()
    xs, ys = normalize(xs), normalize(ys)
    xs, ys = (scale_and_move(xs, 117, 375),
              scale_and_move(ys, -206.5, 894))
    
    drawCurve(xs, ys,thickness=5, col=path_c, filled=True)
    drawPoints(xs, ys,radius=4, col=node_c)
