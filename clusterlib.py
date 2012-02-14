import csv
import numpy as np
import sys
from collections import defaultdict
from matplotlib.colors import *
from pylab import *
from matplotlib.widgets import Slider, Button, RadioButtons
import matplotlib.pyplot as plt
from collections import defaultdict,OrderedDict
from math import pow

barrel_phi_min=0#inclusive
barrel_phi_max=120#exclusive
barrel_theta_min=20#inclusive
barrel_theta_max=68#exclusive

class HitMapFile:
    v=None
    event_index = 0
    theta_index = 1
    phi_index = 2
    e_index = 3
    
    theta_min=0 #inclusive
    theta_max=77#exclusive
    num_theta=77
    
    phi_min=0 #inclusive
    phi_max=265 #exclusive
    num_phi=265 
    def __init__(self,fname=None):
        if fname is not None:
            self.readfile(fname)
        
    def readfile(self,fname):
        f = open(fname)
        tmp = [x for x in csv.reader(f,delimiter=' ')]
        for i in f:
            print self.tmp
        tmp = [ [int(x[0]), int(x[3]), int(x[6]), float(x[9])] for x in tmp] 
        if self.v is None: 
            self.v = tmp
        else:
            self.v+=tmp
        f.close()
    
    def hitmaps(self):
        hitmaps = [HitMap(self.num_phi,self.num_theta) for i in range(self.numEvent())]
        for x in self.v: hitmaps[x[self.event_index]].acc(x[self.phi_index],x[self.theta_index],x[self.e_index])
        for c in hitmaps: c.compute_laplacian()
        return hitmaps
            
    def numEvent(self):
        self.event_max = max(x[self.event_index] for x in self.v)
        self.num_event = self.event_max+1
        return self.num_event
         
    def process(self): #find all the bounds
        # self.theta_min = min(x[self.theta_index] for x in self.v)
        # assert(self.theta_min==0)
        # self.theta_max = max(x[self.theta_index] for x in self.v)
        # self.num_theta = self.theta_max+1
        # self.phi_min = min(x[self.phi_index] for x in self.v)
        # assert(self.phi_min==0)
        # self.phi_max = max(x[self.phi_index] for x in self.v)
        # self.num_phi = self.phi_max+1
        # self.event_min = min(x[self.event_index] for x in self.v)
        # assert(self.event_min==0)
        self.event_max = max(x[self.event_index] for x in self.v)
        self.num_event = self.event_max+1 

class HitMap:
    
    def __init__(self,num_phi,num_theta):
        #hits is structred as theta,phi
        self.hits = np.zeros((num_phi,num_theta))
        
    def acc(self,phi,theta,energy): #x are array(eventno thetaindex, phiindex, energy) from v
        self.hits[phi,theta]+=energy #imshow maps first index to y axis and second to x axis
        
    def compute_laplacian(self):
        self.lpc = laplacian(self.hits)
        
    def sumE(self,cluster):
        return self.sumE_from_hits(self.hits,cluster)
    
    @classmethod
    def sumE_from_hits(self,hits,cluster):
        xlist,ylist = zip(*cluster)
        return np.sum(hits[xlist,ylist])

class Clustering:
    seed_cutoff = 0.025 #25MeV default
    expand_cutoff = 0.005 #5MeV
    directions = [[-1,0],[0,-1],[1,0],[0,1]] #how cluster look around uldr
    def __init__(self):
        pass
    #return ordereddict of seed
    def find_seed(self,hitmap):
        hm=hitmap
        seedlist = zip(*np.where(hm>self.seed_cutoff))
        seedlist.sort(key=lambda x: hm[x])        
        od = OrderedDict()
        for p in seedlist: 
            od[p] = hm[p]
        return od
    #return list of set of tuple of corrdinates
    #each set is one cluster
    #also note that seedod is passed by reference and will have it value change(empty upon return)
    def find_clusters(self,hitmap,seedod=None):
        hm = hitmap.hits
        #making seedlist put in ordered dict
        if seedod is None: seedod = self.find_seed(hm)
        clusters = []
        while len(seedod)!=0:
            seed,E = seedod.popitem()
            #print seed,E
            cluster_so_far = set()
            cluster_so_far.update([seed])
            this_cluster = self.expand_cluster(seed,hm,cluster_so_far,seed)
            #remove seed if seed in this cluster
            for hit_pos in this_cluster:
                if hit_pos in seedod: del seedod[hit_pos]
            clusters.append(this_cluster)
        return clusters

    #recursively expand cluster from given seed
    #note that cluster_so_far is passed by reference and it acts as accumulator
    #last argument is the original seed from this cluster(for calculating upper cutoff)
    def expand_cluster(self,seed,hits,cluster_so_far,org_seed):
        #look in 4 direction
        for direction in self.directions:
            neighbor = self.add_direction(seed,direction)
            if neighbor not in cluster_so_far and self.passcut(neighbor,hits,org_seed):
                cluster_so_far.update([neighbor])
                self.expand_cluster(neighbor,hits,cluster_so_far,org_seed)
        return cluster_so_far

    def add_direction(self,org,direc):
        return tuple([org[0]+direc[0],org[1]+direc[1]])

    def passcut(self,pos,hits,org_seed):
        return hits[pos]>self.expand_cutoff

#assume seed list is sorted
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

class Visualizer:
    @classmethod
    def show_hits(self,hits, ax=None,cutoff=0.0005):
        if ax is None: ax = gca()
        toshow = np.copy(hits)
        if(cutoff is not None):
            toshow[toshow<cutoff] = None
        img = ax.imshow(toshow,aspect='auto',interpolation='nearest',origin='lower')  
        ax.set_xlim((barrel_theta_min,barrel_theta_max))
        ax.set_ylim((barrel_phi_min,barrel_phi_max))
        ax.grid(True,which='both')
        return ax,img
    
    @classmethod
    def show_cluster(self,clusters,ax=None,hits=None):
        if ax is None: ax = gca()
        p=[]
        
        for cluster in clusters:
            
            cx,cy = zip(*cluster)
            #print cx,cy
            #yep y then x it's the way imshow works
            q = ax.plot(cy,cx,'o')
            p.append(q)
            if hits is not None:
                Ecl = HitMap.sumE_from_hits(hits,cluster)
                #print Ecl,cluster
                ax.annotate('%5.4f'%Ecl,xy=(cy[0],cx[0]))
        ax.set_xlim((barrel_theta_min,barrel_theta_max))
        ax.set_ylim((barrel_phi_min,barrel_phi_max))
        ax.grid(True)
        return ax,p
    
    @classmethod
    def show_seeds(self,seeds,hits, ax=None,cutoff=None):
        if ax is None: ax = gca()
        p=[]
        for seed in seeds:
            #print seed
            q = ax.plot(seed[1],seed[0],'x')
            p.append(q)
        ax.set_xlim((barrel_theta_min,barrel_theta_max))
        ax.set_ylim((barrel_phi_min,barrel_phi_max))
        ax.grid(True,which='both')
        return ax,p
