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
#print num_event
barrel_theta_min=0#inclusive
barrel_theta_max=120#exclusive
barrel_phi_min=20#inclusive
barrel_phi_max=68#exclusive
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
    def min(self):
        return np.amin(self.hits,norm=Normalize(vmin=self.vmin,vmax=self.vmin),
           vmin=self.vmin,vmax=self.vmax,
           aspect='auto',interpolation='nearest')

#these operations return array the same size as a
#note that user is responsible for selecting only the region that makes sense
#usually it's a[1:-1,1:-1]
#or your could filter out the deges by 
# a[0:,:]=0;a[:-1,:]=0;a[:,0:]=0;a[:,:-1]=0 etc.
#0,0 is defined as top left first index is row and second index is column
#what this does is the following
#  |ul | u | ur |
#  | l | c | r |
#  | bl| b | br|
# forevery element center at c new value = c_arg*c + ul*ul_arg .... br*br_arg
def nine_op(a,ul=0,u=0,ur=0,l=0,c=0,r=0,bl=0,b=0,br=0):

   ret = np.copy(a)
   ret*=c

   ret[1:,1:] += ul*a[:-1,:-1]
   ret[1:,:] += u*a[:-1,:]
   ret[1:,:-1] += ur*a[:-1,1:]

   ret[:,1:] += l*a[:,:-1]
   ret[:,:-1] += r*a[:,1:]

   ret[:-1,1:] += bl*a[1:,:-1]
   ret[:-1,:] += b*a[1:,:]
   ret[:-1,:-1] += br*a[1:,1:]
   return ret

def nine_opa(a,c):
   return nine_op(a,c[0,0],c[0,1],c[0,2],c[1,0],c[1,1],c[1,2],c[2,0],c[2,1],c[2,2])

#perform cross operation this is just to save some time
def cross_op(a,u=0,l=0,c=0,r=0,b=0):
   ret = np.copy(a)
   ret*=c
   ret[1:,:] += u*a[:-1,:]
   ret[:,1:] += l*a[:,:-1]
   ret[:,:-1] += r*a[:,1:]
   ret[:-1,:] += b*a[1:,:]
   return ret

#compute laplacian; poor man edge detection
def laplacian(a):
   return cross_op(a,c=4,r=-1,l=-1,u=-1,b=-1)

#return magnitude of gradient
def grad2(a):
    g_x = cross_op(a,l=1,r=-1)
    g_y = cross_op(a,u=1,b=-1)
    return np.sqrt(g_x*g_x+g_y*g_y)

#3x3 normalized gaussian blur
def gaussian_blur(a,sigma=1):
   x = np.array([[-1,0,1],[-1,0,1],[-1,0,1]])
   y = np.array([[-1,-1,-1],[0,0,0],[1,1,1]])
   unnorm_gau = np.exp(-(x**x+y**y)/(2*sigma));
   norm_gau = unnorm_gau/np.sum(unnorm_gau)
   return nine_opa(a,norm_gau)

def test_op():
   a = np.zeros((5,5))
   a[:,:]=0.5
   print a
   print laplacian(a)
   print gaussian_blur(a)
#test_op()

def plot_barrel(a,ax,theta_min=barrel_theta_min, theta_max=barrel_theta_max, phi_min=barrel_phi_min,phi_max=barrel_phi_max):
    if ax is None: ax=gca()
    self.img = ax.imshow(toshow,
           aspect='auto',interpolation='nearest')  
    ax.set_xlim((theta_min,theta_max))
    ax.set_ylim((phi_min,phi_max))
    return ax

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