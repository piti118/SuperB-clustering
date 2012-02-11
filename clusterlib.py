import csv
import numpy as np
from collections import defaultdict
from matplotlib.colors import *
from pylab import *
from matplotlib.widgets import Slider, Button, RadioButtons
import matplotlib.pyplot as plt
from math import pow
def readfile():
    f = open('qa-recluster-x2-crystals.txt')
    v = [x for x in csv.reader(f,delimiter=' ')]
    v = [ [int(x[0]), int(x[3]), int(x[6]), float(x[9])] for x in v]
    return v
v=readfile()
event_index = 0
theta_index = 1
phi_index = 2
e_index = 3
theta_min = min(x[theta_index] for x in v)
assert(theta_min==0)
theta_max = max(x[theta_index] for x in v)
num_theta = theta_max+1
phi_min = min(x[phi_index] for x in v)
assert(phi_min==0)
phi_max = max(x[phi_index] for x in v)
num_phi = phi_max+1
event_min = min(x[event_index] for x in v)
assert(event_min==0)
event_max = max(x[event_index] for x in v)
num_event = event_max+1
print num_event
barrel_theta_min=0
barrel_theta_max=130
barrel_phi_min=15
barrel_phi_max=70
num_gui_fig = 4


class Cluster:
    vmin = 1e-2
    vmax = 1e-1
    cutoff = 1e-5
    img=None
    def __init__(self):
        #hits is structred as theta,phi
        self.hits = np.zeros((num_theta,num_phi))
        pass
    def acc(self,x): #x are array(eventno thetaindex, phiindex, energy) from v
        self.hits[x[theta_index],x[phi_index]]+=x[e_index]
    def draw(self,ax=None):
        if ax is None: ax = gca()
        toshow=np.copy(self.hits)
        #print np.amax(self.hits)
        toshow[toshow<self.cutoff]=None
        self.img = ax.imshow(toshow,
               #norm=LogNorm(vmin=1e-8,vmax=0.1),
               norm=Normalize(vmin=self.vmin,vmax=self.vmin),
               vmin=self.vmin,vmax=self.vmax,
               aspect='auto',interpolation='nearest')
        ax.set_xlim((barrel_theta_min,barrel_theta_max))
        ax.set_ylim((barrel_phi_min,barrel_phi_max))
        #colorbar()
        #print dir(self.img)
        ax.grid(True)
    def update(self):
        toshow=np.copy(self.hits)
        #print np.amax(self.hits)
        print self.cutoff
        toshow[toshow<self.cutoff]=None
        self.img.set_data(toshow)
    def min():
        return np.amin(self.hits,norm=Normalize(vmin=self.vmin,vmax=self.vmin),
           vmin=self.vmin,vmax=self.vmax,
           aspect='auto',interpolation='nearest')

def main():
    clusters = [Cluster() for i in range(num_event)]
    for x in v: clusters[x[event_index]].acc(x)
    fig = plt.figure()
    for i in range(4):
        fig.add_subplot(2,2,i+1)
        clusters[i].draw()
    fig.subplots_adjust(left=0.25, bottom=0.25)
    axcolor = 'lightgoldenrodyellow'
        
    axcut = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
    cutoff = Slider(axcut, 'Cutoff', -5, -1, valinit=-3)

    def update(v):
        for i in range(4):
            clusters[i].cutoff=pow(10,cutoff.val)
            clusters[i].update()
        
    cutoff.on_changed(update)

    plt.show()

if __name__ == '__main__':
    main()