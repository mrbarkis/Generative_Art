import random
import functools
from vector import Vector
from particles import Ball
from boundaries import Polygon
from oracles import BoundException
from custom_exceptions import UnboundException
from collections import defaultdict
import math
from subroutines import mix_colors

def get_rnd_unit_vector():
    return Vector(x=random.uniform(-1, 1), y=random.uniform(-1, 1)).normalize()
    

class Presenter(object):
    
    def __init__(self, model, view, oracles):
        self.model = model
        self.view = view
        self.oracles = oracles
        self.predict_collisions()
        self.chart= None
        self.prediction = None

    def add_rnd_ball(self, x, y, r_range, rho_range, e_ave, tint_min, colors):
        r_rnd = random.randrange(*r_range)
        rho_rnd = random.randrange(*rho_range)
        m_rnd = rho_rnd * math.pi * r_rnd ** 2
        v_rnd = math.sqrt(2 * e_ave / m_rnd)
        tint_rnd = (1 - tint_min) / (rho_range[1] - rho_range[0]) * (rho_range[0] - rho_rnd) + tint_min
        color_rnd = mix_colors(random.choice(colors), [255,255,255], tint_rnd)
        self.add_ball(x=x,
                           y=y,
                           v_abs=v_rnd,
                           m=m_rnd,
                           r=r_rnd,
                           c=color_rnd)
        
        
    def show(self):
        '''Draw the current model.'''
        self.view.draw_model(self.model, self.prediction)
                
  
    def predict_collisions(self):
        ''' Repredict the possible collisions.'''
        collisions = [self.predict_next_collision(b) for b in self.model.balls]         
        self.model.collisions = sorted(collisions)
    
    
    def find_affected(self, collision):
        ''' Find collisions that share colliders with the collision.'''
        affected_collisions = []
        affected_balls = collision.get_balls()
        for b in collision.get_balls():
            for c in self.model.collisions:
                if c.involves(b):
                    affected_collisions.append(c)
                    affected_balls = affected_balls.union(c.get_balls())
                
        return affected_balls, affected_collisions 


    def update_collisions(self, collision):
        ''' Repredict only the affected collisions.'''
        affected_balls, affected_collisions = self.find_affected(collision)
        new_collisions = [self.predict_next_collision(b) for b in affected_balls]
        self.model.remove_collisions(affected_collisions)
        self.model.merge_collisions(sorted(new_collisions))
        #self.model.next_collision = self.model.collisions.pop(0)


    def predict_next_collision(self, ball):
        next_collision = None
        for o in self.oracles.values():
            c = o(ball, self.model)
            if c and c < next_collision:
                next_collision = c
        if not next_collision:
            raise UnboundException(ball)
            
        return next_collision
    
    def find_involved_with_bound(self, bound):
        '''Find collisions involved with the bound'''
        involved_collisions = []
        
        for collision in self.model.collisions:
            for edge in bound.edges:
                if collision.involves(edge):
                    involved_collisions.append(collision)
            for vertex in bound.vertices:
                if collision.involves(vertex):
                    involved_collisions.append(collision)
        return involved_collisions
    

    def update_chart(self):
        '''Attempt to move the chart to its target value.'''
        if self.chart == None:
            return
        for n, bar in enumerate(self.chart.bars):
            target_bar = bar.get_target_bar()
            if target_bar == None:
                continue
            update_succeeded = self.add_bound(target_bar)
            if update_succeeded:
                self.chart.bars[n] = target_bar
                self.model.bounds.remove(bar)
        
        
    def update_prediction(self):
        ''' Move the chart to its target value.'''
        if self.chart == None:
            return
        for n, bar in enumerate(self.prediction.bars):
            target_bar = bar.get_target_bar()
            if target_bar == None:
                continue
            self.prediction.bars[n] = target_bar
                
    def set_chart(self, chart):
        self.chart = chart
        for bar in self.chart.bars:
            self.add_bound(bar)
    
    def set_prediction(self, chart):
        self.prediction = chart
           
    def add_bound(self, bound):
        new_collisions = []
        predict_bound_collision = self.oracles[(Ball, Polygon)]
        try:
            for ball in self.model.balls:
                c = predict_bound_collision(ball, [bound])
                if c:
                    new_collisions.append(c)
                    
        except BoundException:
            return False
        
        old_collisions = self.find_involved_with_bound(bound)
        self.model.remove_collisions(old_collisions)
        self.model.merge_collisions(sorted(new_collisions))
        self.model.bounds.append(bound)
        return True
                
    
    def add_ball(self, x, y, v_abs, m, r, c=[0,0,0]):
        new_ball = Ball(p=Vector(x, y), v=get_rnd_unit_vector() * v_abs, m=m, r=r, c=c, painter=self.view.draw_ball)
        try:
            new_collision = self.predict_next_collision(new_ball)
        except BoundException:
            return
        
        self.model.balls.append(new_ball)
        self.model.merge_collisions([new_collision])    
        
        
        
    def advance(self, duration):
            
        while duration > 0:
            if self.model.collisions[0].eta > duration: # Nothing happends during this timeframe
                self.model.advance(dframe=duration)
                duration = 0
            else:
                duration -= self.model.collisions[0].eta 
                self.model.advance(dframe=self.model.collisions[0].eta)
                self.model.collisions[0].handle()
                self.update_collisions(self.model.collisions[0])

                
        
