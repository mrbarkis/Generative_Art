
from abc import ABCMeta, abstractmethod
from math import pi


class View:
    __metaclass__ = ABCMeta #processing supports only python 2.7.
    

      
    # @abstractmethod
    # def draw_ball(self, ball):
    #     pass
    # @abstractmethod
    # def draw_polygon(self, ball):
    #     pass
    
    @abstractmethod
    def play_click(self, ball):
        pass
    
    @abstractmethod
    def draw_ball(self, ball):
        pass

    @abstractmethod
    def draw_pollen(self, ball):
        pass
        
    @abstractmethod
    def draw_model(self, ball, chart):
        pass
        
    
class ProcessingView(View):

    def __init__(self, bg_color=[255, 255, 255], lambda_color=3*[0], lambda_img=None, sound_player=None):
        self.bg_color = bg_color
        self.lambda_img = lambda_img
        self.sound_player = sound_player
        self.lambda_color = lambda_color
        
    def play_click(self):
        self.sound_player()
        
    def draw_ball(self, ball):
        noStroke();
        fill(*ball.c)
        circle(ball.p.x, ball.p.y, 2 * ball.r)
        
    def draw_pollen(self, pollen):
        self.draw_ball(pollen)
        fill(*pollen.c_arc)
        arc(pollen.p.x, pollen.p.y, 2*pollen.r, 2*pollen.r, -pi / 2 , -pi / 2 + 2 * pi * float(pollen.frame % pollen.interval) / pollen.interval )
        fill(*pollen.c_t);
        textSize(pollen.r);
        textAlign(CENTER, CENTER)
        text(pollen.n_collisions, pollen.p.x, pollen.p.y - 0.1 * pollen.r);
        
    
    def draw_polygon(self, polygon):
        if polygon.fill_color == None:
            noFill()
        else:
            fill(*polygon.fill_color)
        
        if polygon.edge_color == None:
            noStroke()
        else:
            stroke(*polygon.edge_color)

        beginShape();
        for v in polygon.vertices:
            vertex(v.x, v.y);
        endShape(CLOSE);
        

        
    def draw_model(self, model, chart=None):
        ''' Draw the background, bounds, and the balls.'''
        background(*self.bg_color)
        if model.pollens[0]:
            font_size= model.pollens[0].r
        else:
            font_size = 20
        
        scale_img = float(font_size) / self.lambda_img.height
        textSize(font_size);
        textAlign(LEFT, TOP)
        fill(*self.lambda_color)
        text("{0:.2f}".format(model.pollens[0].rate_ave), 0 + scale_img * self.lambda_img.width, 0);
        
        image(self.lambda_img, 0, 0, scale_img * self.lambda_img.width, font_size)
        
        for bar in chart.bars:
            self.draw_polygon(bar)
            
        for bound in model.bounds:
            self.draw_polygon(bound)
        
        
        
        for ball in model.balls:
            ball.draw_myself()
            
        
