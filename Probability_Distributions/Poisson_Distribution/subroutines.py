import math
import pickle

def mix_colors(c1, c2, ratio):
    return [(1 - ratio) * a + ratio * b for a, b in zip(c1, c2)]

def poisson(k, rate):
    return rate ** k * math.exp(-rate) / math.factorial(k) 
            
def save_list(list_, path):
    with open(path, 'wb') as file:
        pickle.dump(list_, file)
    
