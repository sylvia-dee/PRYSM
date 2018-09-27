import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def to_plot(datain_a, name, unit, monthly_ave = False):
    matplotlib.style.use('ggplot')
    if monthly_ave == True:
        datain_b = datain_a.reshape(12, datain_a.shape[0]//12)
        mon_ave = datain_b.mean(1)
        data_toplot = mon_ave
    else:
        datain_b = datain_a.reshape(datain_a.shape[0]//12, 12)
        yr_ave = datain_b.mean(1)
        data_toplot = yr_ave
    fig1 = plt.figure(figsize = (8, 8))
    ax1 = fig1.add_subplot(1,1,1)
    plt.plot(data_toplot,  linewidth = 2)
    if monthly_ave == True:
        plt.xlabel('Month', fontsize = 14)
        plt.xlim(-1,13)
        plt.xticks(np.arange(0,12), np.arange(0,12) +1,fontsize = 15)
        name = name + ' monthly average'
    else:
        plt.xlabel('Year', fontsize = 14)
        name = name + ' yearly average'
    plt.ylabel(unit, fontsize = 18)
    plt.yticks(fontsize = 15)
    plt.title(name, fontsize = 17)
    plt.grid(True)
    plt.show()

leafwax_proxy = np.load('leafwax_dDp_sample.npy')

# Russell et al., 2018
leafwax_proxy = np.load('leafwax_dDp_sample.npy')
to_plot(leafwax_proxy, 'leafwax', 'unit', monthly_ave = False)
to_plot(leafwax_proxy, 'leafwax', 'unit', monthly_ave = True)