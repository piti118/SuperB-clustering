import numpy as np
import matplotlib.pyplot as plt
from clusterlib import *
def main():
    #load data
    hitmaps = [HitMap() for i in range(num_event)]
    for x in v: hitmaps[x[event_index]].acc(x)
    for c in hitmaps: c.compute_laplacian()
    
    num_row=2
    num_col=3
    raw_im = [None]*num_row
    seed_im = [None]*num_row
    cluster_im = [None]*num_row
    clustering = Clustering()
    
    fig, axs = plt.subplots(num_row, 3, sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.05,wspace=0.05)
    fig.subplots_adjust(bottom=0.25) #reserve space for controls
    
    for ax in axs[-1,:]: ax.set_xlabel('theta')
    for ax in axs[:,0]: ax.set_ylabel('phi')
    for i in range(num_row):
        hits = hitmaps[i].hits
        _,raw_im[i] = Visualizer.show_hits(hits,axs[i,0],cutoff=0.0005)
        seeds = clustering.find_seed(hits)
        Visualizer.show_seeds(seeds,hits,axs[i,1])
        cluster = clustering.find_clusters(hitmaps[i],seeds)
        Visualizer.show_cluster(cluster,axs[i,2],hits)
        fig.canvas.draw()
    axcolor = 'lightgoldenrodyellow'
    ax_seed_cutoff = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
    ax_min_expand  = fig.add_axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
    
    s_seed_cutoff = Slider(ax_seed_cutoff, 'seed_E_cut_off(GeV)', 0, 0.1, valinit=clustering.seed_cutoff)
    s_min_expand = Slider(ax_min_expand, 'expand_min_E(GeV)', 0, 0.020, valinit=clustering.expand_cutoff)
    def update(v):
        clustering.seed_cutoff=s_seed_cutoff.val
        clustering.expand_cutoff=s_min_expand.val
        print clustering.seed_cutoff, clustering.expand_cutoff
        for i in range(num_row):
            axs[i,1].clear()
            axs[i,2].clear()
            seeds = clustering.find_seed(hitmaps[i].hits)
            Visualizer.show_seeds(seeds,hitmaps[i].hits,axs[i,1])
            cluster = clustering.find_clusters(hitmaps[i],seeds)
            Visualizer.show_cluster(cluster,axs[i,2],hitmaps[i].hits)
            
        
    s_seed_cutoff.on_changed(update)
    s_min_expand.on_changed(update)
    plt.show()
if __name__ == '__main__':
    main()
