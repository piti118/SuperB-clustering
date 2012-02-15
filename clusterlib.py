import csv
import numpy as np
import sys
from collections import defaultdict
from matplotlib.colors import *
from pylab import *
from matplotlib.widgets import Slider, Button, RadioButtons
import matplotlib.pyplot as plt
from collections import defaultdict,OrderedDict
from math import pow,pi,exp,sqrt
from itertools import product, combinations
barrel_phi_min=0#inclusive
barrel_phi_max=120#exclusive
barrel_theta_min=20#inclusive
barrel_theta_max=68#exclusive

def gauss_peak_norm_functor(mean,sigma,peak):
    my_mean = mean
    my_sigma = sigma
    my_peak = peak
    my_norm = peak/exp(-1*(0.0-mean)**2/(2.0*sigma**2))
    def ret(x):
        return exp(-1*(x-my_mean)**2/(2.0*my_sigma**2))*my_norm
    return ret

def distance(x,y):
    return sqrt(x**2+y**2)
def dis(xy1,xy2):
    return sqrt((xy1[0]-xy2[0])**2+(xy1[1]-xy2[1])**2)
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
        tmp = [filter(lambda z: z!='', x) for x in tmp]#remove all the empty string
        tmp = [ [int(x[0]), int(x[1]), int(x[2]), float(x[3])] for x in tmp] 
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
    @classmethod
    def in_barrel(self,phi,theta):
        return barrel_phi_min <= phi < barrel_phi_max and barrel_theta_min <= theta < barrel_theta_max
    
"""
    Cluster class contains cluster as set of tuple (each one represent the point on the hitmap)
    and the seed is the original seed used in making this cluster
"""
class Cluster:
    def __init__(self,seed,cluster):
        self.cluster = cluster
        self.seed = seed

class Clustering:
    seed_cutoff = 0.025 #25MeV default
    expand_cutoff = 0.005 #1MeV
    directions = [[-1,0],[0,-1],[1,0],[0,1]] #how cluster look around uldr
    moire_r = 3.6/5.0 #in the unit of the crystal face length
    
    #the allow region is enveloped by two gaussian
    #upper is a very wide one with peak normalized to the seed energy
    #lower one is a narrow one with peak normalized to half the seed energy
    #there is also an absolute cutoff at expand_cutoff at 1 MeV
    
    upper_sigma_factor_cutoff = 10.0
    upper_norm_cutoff =  2.0
    
    lower_low_cutoff = 0.001
    lower_high_cutoff = 0.005
    lower_rising_factor = 50.0
    #lower_sigma_factor_cutoff = 0.5
    #lower_norm_cutoff = 0.5
    
    def __init__(self):
        pass
        
    def use9x9dir(self):
        self.directions = [[-1,0],[0,-1],[1,0],[0,1],[-1,-1],[1,1],[-1,1],[1,-1]]
    
    def use25x25dir(self):
        x = range(-2,3)
        self.directions = [ (y,z) for y in x for z in x]
    
    def useDiamondDir(self):
        #9x9 plus the 2 straight in each direction
        self.directions = [[-1,0],[0,-1],[1,0],[0,1],[-1,-1],[1,1],[-1,1],[1,-1],[2,0],[-2,0],[0,2],[0,-2]]
        
    #return ordereddict of seed
    def find_seed(self,hitmap):
        hm=hitmap
        seedlist = zip(*np.where(hm>self.seed_cutoff))
        seedlist.sort(key=lambda x: hm[x])        
        od = OrderedDict()
        for p in seedlist: 
            od[p] = hm[p]
        return od
    #return list of cluster object
    #also note that seedod is passed by reference and will have it value change(empty upon return)
    def find_clusters(self,hitmap,seedod=None):
        hm = hitmap.hits
        #making seedlist put in ordered dict
        if seedod is None: seedod = self.find_seed(hm)
        clusters = []
        seeds_ret = [] #list of seed used for each seed used by cluster
        while len(seedod)!=0:
            seed,E = seedod.popitem()
            if not HitMap.in_barrel(*seed): continue; #skip endcaps
            #print seed,E
            cluster_so_far = set()
            cluster_so_far.update([seed])
            this_cluster = self.expand_cluster(seed,hm,cluster_so_far,seed)
            #remove seed if seed in this cluster
            for hit_pos in this_cluster:
                if hit_pos in seedod: del seedod[hit_pos]
            clusters.append(Cluster(seed,this_cluster))
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
        #the modulo is for wrapping
        return tuple([(org[0]+direc[0])%barrel_phi_max,org[1]+direc[1]])

    def calculate_lower_cutoff(self,seedE,dis):
        #____       _____
        #    \_____/
        vpart = abs(dis/self.moire_r/self.lower_rising_factor)*seedE
        low_cutoff = min(self.lower_high_cutoff,max(vpart,self.lower_low_cutoff))
        return low_cutoff
        
    def calculate_upper_cutoff(self,seedE,dis):
        high_cutoff_g = gauss_peak_norm_functor(0,self.moire_r*self.upper_sigma_factor_cutoff,self.upper_norm_cutoff*seedE)(dis)
        high_cutoff = min(high_cutoff_g,seedE)
        return high_cutoff

    def passcut(self,pos,hits,org_seed):
        dis = distance((pos[0]-org_seed[0])%barrel_phi_max,pos[1]-org_seed[1])
        
        #low_cutoff_g =  gauss_peak_norm_functor(0,self.moire_r*self.lower_sigma_factor_cutoff,self.lower_norm_cutoff*hits[org_seed])(dis)
        #low_cutoff_g = dis*abs(hits[org_seed])/50
        #low_cutoff_g = 0
        #low_cutoff = min(low_cutoff_g,self.expand_cutoff)
        low_cutoff = self.calculate_lower_cutoff(dis,hits[org_seed])
        #high_cutoff_g = gauss_peak_norm_functor(0,self.moire_r*self.upper_sigma_factor_cutoff,self.upper_norm_cutoff*hits[org_seed])(dis)
        #high_cutoff = min(high_cutoff_g,hits[org_seed])
        #high_cutoff = 100
        #print dis,low_cutoff,low_cutoff_g,high_cutoff,high_cutoff_g
        #high_cutoff = hits[org_seed]
        high_cutoff = self.calculate_upper_cutoff(dis,hits[org_seed])
        return high_cutoff > hits[pos] > low_cutoff
    
    #note that this is pass by reference clusters will be changed
    def reduce_clusters(self,clusters,hits):
        #taking care of overlapping clusters
        #it first find all the crystals that is in two or more clusters
        #then for each one of them calculate the figure of merit based on seed energy and distance for the seed
        #the one with the highest figure of merits gets the whole crystal
        
        #first put them all in a map by seed
        cl_map = {c.seed : c.cluster for c in clusters}
        
        #dupe_map is map from crystal_pos to list of seed that has this crystal in it
        dupe_map = defaultdict(set)
        #build dupe_map
        for lhs_cl, rhs_cl in combinations(clusters, 2):
            #find intersection for all the intersection
            intersection = lhs_cl.cluster.intersection(rhs_cl.cluster)
            for overlap in intersection:
                dupe_map[overlap].update([lhs_cl.seed,rhs_cl.seed])
        #for each overlapped crystal compute the expected E and remove it from the low expected E clusters
        for crystal_pos, seed_list in dupe_map.items():
            expected_E = [(seed,gauss_peak_norm_functor(0,self.moire_r,hits[seed])(dis(crystal_pos,seed))) for seed in seed_list]
            best_seed = max(expected_E,key=lambda x:x[1])
            #now that we get the seed remove it from other list
            for seed in seed_list:
                if(seed!=best_seed[0]): 
                    cl_map[seed]-=set([crystal_pos])
        #now clusters should have no overlap
        #check if there is a seed with empty cluster
        clusters = [cl for cl in clusters if len(cl.cluster)!=0]
        return clusters
    
    def draw_cutoff(self,E=0.100,ax=None):
        if ax is None: ax = gca()
        #E=0.100
        #high_cutoff_g = gauss_peak_norm_functor(0,self.moire_r*self.upper_sigma_factor_cutoff,self.upper_norm_cutoff*E)
        #vhg = np.vectorize(high_cutoff_g)
        #low_cutoff_g =  gauss_peak_norm_functor(0,self.moire_r*self.lower_sigma_factor_cutoff,self.lower_norm_cutoff*E)
        #lhg = np.vectorize(low_cutoff_g)
        high_f = lambda x: self.calculate_upper_cutoff(E,x)
        low_f = lambda x: self.calculate_lower_cutoff(E,x)
        lv = np.vectorize(low_f)
        hv = np.vectorize(high_f)
        x = np.linspace(-20,20,1000)
        ax.plot(x,lv(x),label='low cutoff')
        ax.plot(x,hv(x),label='high cutoff')
        #ax.set_ylim(ymin=0,ymax=0.01)
        #ax.plot(x,[self.expand_cutoff]*len(x),label='absolute low')
        #ax.plot(x,[E]*len(x),label='absolute high') 
        ax.minorticks_on()
        ax.grid(True,which='both')
        ax.set_xlabel('distance (#of crystal)')
        ax.set_ylabel('E(GeV)')
        ax.legend()
        return ax       
        
        

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
            
            cx,cy = zip(*cluster.cluster)
            #print cx,cy
            #yep y then x it's the way imshow works
            q = ax.plot(cy,cx,'s',ms=2.5,alpha=0.7)
            p.append(q)
            if hits is not None:
                Ecl = HitMap.sumE_from_hits(hits,cluster.cluster)
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
