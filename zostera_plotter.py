#o-*- coding: utf-8 -*-

#this is ploting function for the output of the  simulation function of full_class.py
#this function makes a series of plots in such a way that it represents the steps of a simulation
#each frame a plot is made of each individual in the meadow over a world
#and secondary  plots that show other information about the meadow  are also made.
#to save the animations this function needs 'ffmpeg', whose complete route may need to be specified

from __future__ import division
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec
from matplotlib import ticker
import sys
import os
import scipy.stats
import numpy as np
import matplotlib.ticker as tkr #to change the tick names

def plot_meadow(output, irradiance, world_stuff, save, show):
    (grid_x, grid_y, world_map) = world_stuff
    #output has the elements points, ages, lengths, demogra, zos_num
    #extract the information
    points = output[0]
    ages = output[1]
    lengths = output[2]
    demogra = output[3]
    zos_num = output[4]
    date_str = demogra[3]
    #the plot could be of the whole world, but sometimes when it is too large, the individual wont be appreciated,
    #the next lines set the limit of the area to be plotted as the area where there are individuals.
    lvl1 = [item for sublist in points for item in sublist]#flatten  dimensions of points matrix 
    lvl2 = [item for sublist in lvl1 for item in sublist]#
    max_x = max([item[0] for item in lvl2])
    min_x = min([item[0] for item in lvl2])
    max_y = max([item[1] for item in lvl2])
    min_y = min([item[1] for item in lvl2])

    fig = plt.figure()#make the plot canvas
    fig.set_size_inches(7,5)#set size

    def animate(t):#the animation will be made using the fuctionanimation
        #subplots
        fig.clear()#
        ax3 = plt.subplot2grid((2, 3), (0, 2))#make 3 subplots, one for the plot of the individuals, other for a 
        #traceplot of the population size, and other for an histogram of the age structure of the population
        ax2 = plt.subplot2grid((2, 3), (1, 2))
        ax1 = plt.subplot2grid((2, 3), (0, 0), colspan=2, rowspan=2)
        plt.subplots_adjust(hspace=0.8,wspace=0.8)#
        
        #the first subplot is an histogram of the age of each phytomer at the time 
        ax3.set_title('')
        datos_p_histogram = ages[t]
        #to determine the appropriate number of bins use the e freedman-diaconis rule
        iqr = scipy.stats.iqr(datos_p_histogram) # inter quantile range  the difference between the  75 and 25 quantile
        bin_width = max(1,2*(iqr/(len(datos_p_histogram)**(1./3))))#if it is less than one make it 1
        if sys.version_info[0] > 2:  #  normalized is for python 2, density for python 3
            ax3.hist(datos_p_histogram, density=True, facecolor='navy', 
                bins=np.arange(min(datos_p_histogram), max(datos_p_histogram) + bin_width, bin_width))
        else:
            ax3.hist(datos_p_histogram, normed=True, facecolor='navy', 
                bins=np.arange(min(datos_p_histogram), max(datos_p_histogram) + bin_width, bin_width))
        ax3.set_ylabel('Frequency')
        from matplotlib.ticker import FormatStrFormatter
        ax3.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))# set tick name decimal significant decimal points
        ax3.set_xlabel('Age (plastochrones)')
        ax3.grid(True)
        ax3.set_title('Age of internodes')
        #########################################################
        #scatter plot of the meadow size
        ti = range(0, t + 1)
        individuos = demogra[0][0:t + 1]
        #for all the points to be visible 
        #the min and max are because the range could be 0
        rango = max(max(demogra[0]) - min(demogra[0]),max(demogra[0]))
        ax2.plot(ti, individuos, 'b-', linewidth=3, label='No. of Individuals')
        ax2.set_ylim([max(0,rango*0.125), rango*1.25])
        ax2.grid(True)
        ax2.set_ylabel('Number')
        ax2.set_xlabel('Fortnight')
        ax2.set_title('Population Size')
        #########################################################################
        #draw the meadow over the map 

        ax1.set_title(date_str[t])        
        def numformat(x, pos): # function to divide the tick value  by 1000 to convert from mm to m
            s = '{}'.format(x / 1000)
            return s

        yfmt = tkr.FuncFormatter(numformat)#using the function to change the ticks
        ax1.set_xlim([min_x,max_x])
        ax1.set_ylim([min_y,max_y])#set the minimum and maximum of map, if want the whole map of the world comment this
        #and the previous line
        ax1.yaxis.set_major_formatter(yfmt)#adjust the values of the ticks
        ax1.xaxis.set_major_formatter(yfmt)
        ax1.set_xlabel('m')
        ax1.set_ylabel('m')
        maximo_irrad =100 #the maximun value to be used  in the pcolor plot
        im = ax1.pcolor(grid_y, grid_x,irradiance[:,:,t], vmin=0, vmax=maximo_irrad, cmap="viridis", edgecolors='k',)#make
        #a pseucolor plot of the world map 
        cont = ax1.contour(grid_y, grid_x,world_map/1000, 15, colors='k')# make a contour plot of the depth of the  world
        #to plot the bathymetry of the world
        plt.clabel(cont, fontsize=7, inline=1)
        #to set the properties of the color bar
        bar = fig.colorbar(im, ax=ax1, orientation='horizontal')#create the bar at the position wanted
        loca = ticker.MaxNLocator(nbins=7)
        bar.locator = loca
        bar.update_ticks()
        bar.set_label(r'Irradiance ($KW\ fortnight \ m^{-2}$)')
        #there are two option to plot the individuals, one is to plot the all green and make the age of phytomer a shade
        #of green darker the older it is
        #the other is to plot each individual as a unique color, the age of each phytomer will be a shade of the color
        #of the color the older it is, this is useful in testing or debugging 
        #if you dont want to plot each individual as a unique color comment the next two lines
        #colorz = list(set(zos_num[t]))
        #numberz = range(0, len(colorz))
        for z in points[t]:  # plot each point
            index = points[t].index(z)
            #to plot all the points in shades of green use the next line
            ax1.plot([z[0][0], z[1][0]], [z[0][1], z[1][1]], marker='.',
                linestyle='-', color=[0,1 - (ages[t][index] / 100.0), 0], linewidth=2)
            #to plot the points in different colors for each individual use the next line (and the previous 4 and 5 five lines
            #must be uncommented)
            #col = (numberz[colorz.index(zos_num[t][index])]) / (1.0*len(numberz))
 #          ax1.plot([z[0][0], z[1][0]], [z[0][1], z[1][1]], marker='.',
 #               linestyle='-', color=[col,
 #                    1 - (ages[t][index] / 100.0), 1-col],
 #               linewidth=3)
        ax1.set_aspect('equal')
        #plt.tight_layout()

    ani = animation.FuncAnimation(fig, animate, frames=len(ages))#make animation
    if show is True:#if the animation is to be showed
        print(("Preparing animation..."))
        plt.show()#show animation
    if save is True:#if the animation is to be saved
        fps=len(output[2])/20#set the number of frames per second
        print("Saving animation...")
        try:
            Writer = animation.writers['ffmpeg']
        except RuntimeError:#ffmpeg program must be installed. The path to the program can  be specified in the next line
            plt.rcParams['animation.ffmpeg_path'] = 'C:/bla/bla/ffmpeg' # write the path to ffmpeg
            Writer = animation.writers['ffmpeg']
        writer = Writer(fps=fps, bitrate=1800)
        if not (os.path.exists(os.getcwd() + os.sep + 'outputs')):
            os.mkdir('outputs')

        ani.save(os.getcwd() + os.sep + 'outputs' + os.sep + 'animation_meadow.mp4', writer=writer)

