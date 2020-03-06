import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def draw_conformations(conf, structs_per_row = 5):
    n_beads = conf.x.shape[0]
    n_confs = conf.x.shape[1]
    for i in range(n_confs):
        
        plt.subplot(int(n_confs / structs_per_row)+1, structs_per_row, i+1)

        plt.plot(conf.x[:,i], c='k')

        for j_ind, j in enumerate(conf.d[:,i]):
            if j == 1:
                plt.scatter([j_ind],[conf.x[j_ind,i]],marker=10,c='red')
            elif j == -1:
                plt.scatter([j_ind],[conf.x[j_ind,i]],marker=11,c='orange')

        for j_ind, j in enumerate(conf.p[:,i]):
            if j != -1:
                if j_ind > j:
                    plt.plot([j, j_ind],[conf.x[j_ind,i],conf.x[j_ind,i]],c='blue', linewidth=8)

        plt.text(0,n_beads,"%.2f" % conf.energies[i])
        
        plt.xlim([-1,n_beads+1])
        plt.ylim([-n_beads,n_beads])
        
        plt.axhline(0,c='grey',linewidth=0.5)
        plt.axis('off')

def draw_traces(conf, boltz_weight=False):
    n_beads = conf.x.shape[0]
    n_confs = conf.x.shape[1]
    
    min_E = np.min(conf.energies)
    for i in range(n_confs):
        if boltz_weight:
            wt = np.exp(-0.01*(conf.energies[i]+min_E))
        else:
            wt=1
        plt.plot(conf.x[:,i]+np.random.normal(scale=0.1,size=n_beads), linewidth=wt, c='k', alpha=0.5)
        for j_ind, j in enumerate(conf.p[:,i]):
            if j != -1:
                if j_ind > j:
                    plt.plot([j, j_ind],[conf.x[j_ind,i],conf.x[j_ind,i]], c='blue', linewidth=wt, alpha=0.5)

    plt.axis('off')

def plot_seqlogo(sequences):
    '''sequences: list of sequences'''

    U_c, A_c, C_c, G_c = sns.color_palette('pastel',4)
    N = len(sequences[0])
    dct = {k:{'A':0,'C':0,'G':0,'U':0} for k in range(N)}
    
    for j, seq in enumerate(sequences):
        for i in range(N):
            dct[i][seq[i]] += 1
    plt.figure(figsize=(1*N,1))
    for i in range(N):
        plt.bar(i,1,color=A_c, label='A')
        plt.bar(i,1-dct[i]['A']/len(sequences),color=U_c, label='U')
        plt.bar(i,(dct[i]['C']+dct[i]['G'])/len(sequences),color=G_c, label='G')
        plt.bar(i,(dct[i]['C'])/len(sequences),color=C_c, label='C')

    plt.xlim([-1,N])
    #plt.legend(bbox_to_anchor=(1,1))

    #todo: get legend working
    