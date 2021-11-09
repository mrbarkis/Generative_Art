add_library('sound')

from vector import Vector
from particles import Ball, Pollen
from boundaries import Box, Polygon, BarChart
from subroutines import mix_colors, poisson, save_list
from oracles import PolygonOracle, BallOracle
from view import ProcessingView
from model import Model
from presenter import Presenter
import time
from timeit import default_timer as timer
from oracle_functions import predict_earliest_ball_hit, predict_earliest_bound_hit
import random
import math
import os
import pickle

WINDOW_SIZE =(1920, 1080) # (640, 360)# 

RECORD = True
FOLDER_FRAMES = os.path.join('.','saved_frames\\')
FILE_COLLISION_TIMES = os.path.join('.', 'script_to_generate_audiotrack', 'collision_times.pickle')
print(FILE_COLLISION_TIMES)
MARGINS = (-20, -20, -20, 0) # left, top, right, bottom margins
#N = 3 # number of balls
SPAWN_INTERVAL = 10 # Spawn balls every x frames.
SPAWN_LOCATIONS = [(-10, WINDOW_SIZE[1] / 2), (WINDOW_SIZE[0] + 10, WINDOW_SIZE[1] / 2)]
SPAWN_START_FRAME = 10 * 60
# PHYSICS
AVERAGE_ENERGY = 100 # Pixels per frame
RADIUS_RANGE = (3, 6)
DENSITY_RANGE = (1, 3)
###    COLORS    ####

TINT_MIN = .5
RED, BLUE, BROWN, WHITE, CHAR= [114, 0, 16], [0, 73, 114], [114, 41, 0], [255, 255, 255], [54, 69, 79]
BG_COLOR = WHITE #mix_colors(BROWN, WHITE, 0.9)
TEXT_COLOR = mix_colors(BROWN, WHITE, 0.5)
BOX_EDGE_COLOR = mix_colors(BROWN, WHITE, 0.4)
BOX_FILL_COLOR = mix_colors(BROWN, WHITE, 0.95)
CLOCK_COLOR =  mix_colors(BLUE, WHITE, 0.5)
PREDICTION_COLOR = mix_colors(BROWN, WHITE, 0.6) #BROWN + [100]
CHART_COLOR = BLUE + [150]
CHART_ACTIVE_COLOR = BLUE + [200]
POLLEN_COLORS = (mix_colors(BROWN, WHITE, 0.4), mix_colors(BROWN, WHITE, 0.5), CHART_ACTIVE_COLOR)
BALL_COLORS = [RED, BLUE]



def setup():
    global presenter, chart, chart_prediction, paused, img, spawn_counter
    paused = False
    spawn_counter = -SPAWN_START_FRAME
    
    size(*WINDOW_SIZE)
    sound_player = SoundFile(this, "click.wav").play
    #print(sf)
    #sf.play()
    lambda_img = loadImage("lambda.png")
    view = ProcessingView(bg_color=BG_COLOR,lambda_color=PREDICTION_COLOR, lambda_img= lambda_img, sound_player=sound_player)

    container = Box(
                x=MARGINS[0],
                y=MARGINS[1],
                w=WINDOW_SIZE[0] - MARGINS[0]-MARGINS[2],
                h=WINDOW_SIZE[1] - MARGINS[1] - MARGINS[3],
                edge_color=BOX_EDGE_COLOR,
                fill_color=None
                )
    
    obstacle = Box(
                x=100,
                y=100,
                w=200,
                h=100,
                edge_color=BOX_EDGE_COLOR,
                fill_color=None
                )

    chart = BarChart(
                    pmf=lambda k: (k==0) , # delta pmf at k=0
                    pos=(0, height/2),
                    height_= height/2,
                    width_=width,
                    speed=5,
                    n_bins=20,
                    n_samples=20,
                    edge_color=WHITE,
                    fill_color=CHART_COLOR,
                    active_color=CHART_ACTIVE_COLOR
                    )
    chart_prediction = BarChart(
                    pmf=lambda k: poisson(k, 0) , # delta pmf at k=0
                    pos=(0, height/2),
                    height_= height/2,
                    width_=width,
                    speed=1,
                    n_bins=20,
                    n_samples=20,
                    edge_color=WHITE,
                    fill_color=PREDICTION_COLOR,
                    )
    global recording
    recording = []
    when_colliding = [lambda pollen: chart.add_to(pollen.n_collisions),
                      lambda pollen: view.play_click(),
                      lambda pollen: recording.append(pollen.frame)]
    pollen =  Pollen(
                p=Vector(WINDOW_SIZE[0] / 2,
                y=WINDOW_SIZE[1] / 2),
                v=Vector(3.1, 0.0001),
                m = .5 * 50 ** 2, r=50, c=POLLEN_COLORS,
                painter = view.draw_pollen,
                interval=60,
                n_samples = 40,
                when_colliding = when_colliding,
                when_resetting = [lambda pollen: chart.add_to(0),
                                  lambda pollen: chart_prediction.update_pmf(lambda k: poisson(k, pollen.rate_ave)),
                                  ],
                )
    
    
    oracle_funcs = {
        (Ball, Ball):predict_earliest_ball_hit,
        (Ball, Polygon): predict_earliest_bound_hit
    }
    
    model = Model(balls=[], pollens= [pollen], bounds=[container])
    presenter = Presenter(model, view, oracle_funcs)
    presenter.set_chart(chart)
    presenter.set_prediction(chart_prediction)
    

    fullScreen()
    
def draw():
    global paused, spawn_counter

    presenter.show()
    if mousePressed:
        presenter.add_rnd_ball(x=mouseX,
                     y=mouseY,
                     r_range=RADIUS_RANGE,
                     rho_range=DENSITY_RANGE,
                     e_ave=0.001,
                     tint_min=TINT_MIN,
                     colors=BALL_COLORS)
        
    if keyPressed:
        time.sleep(.1) # To avoid multiple presses
        if key in ('p','P'):
            paused = not paused
        if key in ('q', 'Q'):
            save_list(list_=recording, path=FILE_COLLISION_TIMES)
            exit()
  
    if paused:
        return
              
    if spawn_counter == SPAWN_INTERVAL:
        spawn_counter = 0
        x, y = random.choice(SPAWN_LOCATIONS)
        presenter.add_rnd_ball(x=x,
                     y=y,
                     r_range=RADIUS_RANGE,
                     rho_range=DENSITY_RANGE,
                     e_ave=AVERAGE_ENERGY,
                     tint_min=TINT_MIN,
                     colors=BALL_COLORS)
    else:
        spawn_counter += 1
        


   
    start = timer()
    presenter.advance(duration=1)
    end = timer()
    presenter.update_chart()
    presenter.update_prediction()
    if RECORD:
        saveFrame("{}####.png".format(FOLDER_FRAMES))
        if len(presenter.model.balls) % 100 == 0:
            save_list(list_=recording, path=FILE_COLLISION_TIMES)
        
        
    #print(end - start)

    
