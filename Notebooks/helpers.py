#!/usr/bin/env python
# coding: utf-8

import sys, re
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib import colors, cm
from matplotlib.offsetbox import AnchoredText
plt.rcParams["font.family"] = "Arial"

alphadict = {0.05:"*", 0.01:"**", 0.001:"***", "NS":"(NS)"}

def tpm_plot(countsDF, column, target, ax, title, report="median", sleuthDF=None, titlefontsize=18, fontsize=12, xticks='full', twolinetitle=False):
    
    ax.set_xticks(range(0, 54, 2))
    if xticks=='reduced':
        ax.set_xticklabels(["","","","","10","","", "16", "","","22","","","4",""], size=fontsize)
    else:
        ax.set_xticklabels([i%24 for i in range(0, 54, 2)], size=fontsize)
    

    subset = countsDF[countsDF[column]==target]      
    subset = subset.groupby(["individual", "time_point"]).sum()["tpm"].reset_index().groupby("time_point").describe()

    if report=="median":
        ax.plot([10, 16, 22, 28], subset["tpm"]["50%"][:4], "o", color="black", linewidth=2.0)
        ax.errorbar(x=[10, 16, 22, 28], y=subset["tpm"]["50%"][:4], 
                    yerr=(subset["tpm"]["50%"][:4]-subset["tpm"]["25%"][:4], 
                            subset["tpm"]["75%"][:4]-subset["tpm"]["50%"][:4]), 
                    color="black", ls='-', capsize=5.0, linewidth=2.0)
        ax.plot([10, 16, 22, 28], subset["tpm"]["50%"][4:], "o", color="xkcd:reddish", linewidth=2.0)
        ax.errorbar(x=[10, 16, 22, 28], y=subset["tpm"]["50%"][4:], 
                    yerr=(subset["tpm"]["50%"][4:]-subset["tpm"]["25%"][4:], 
                            subset["tpm"]["75%"][4:]-subset["tpm"]["50%"][4:]), 
                    color="xkcd:reddish", ls='-', capsize=5.0, linewidth=2.0)
    elif report=="mean":
        ax.plot([10, 16, 22, 28], subset["tpm"]["mean"][:4], "o", color="black", linewidth=2.0)
        ax.errorbar(x=[10, 16, 22, 28], y=subset["tpm"]["mean"][:4], 
                    yerr=(subset["tpm"]["std"][:4]), color="black", ls='-', capsize=5.0, linewidth=2.0)
        ax.plot([10, 16, 22, 28], subset["tpm"]["mean"][4:], "o", color="xkcd:reddish")
        ax.errorbar(x=[10, 16, 22, 28], y=subset["tpm"]["mean"][4:], 
                    yerr=(subset["tpm"]["std"][4:]), color="xkcd:reddish", ls='-', capsize=5.0, linewidth=2.0)
    else:
        print("report must be either 'median' or 'mean'")
        return
    ax.set_xlim(8, 30)
    
    patch_height = ax.get_ylim()[1]
    rect1 = mpatches.Rectangle(xy=(0, 0), width=6, height=patch_height, color='grey', alpha=0.5)
    rect2 = mpatches.Rectangle(xy=(20, 0), width=10, height=patch_height, color='grey', alpha=0.5)
    ax.add_patch(rect1)
    ax.add_patch(rect2)
    
    ymax = ax.get_ylim()[1]
    ax.set_ylim(0,ymax)
    oldyticks = ax.get_yticks()
    newyticks = np.append(oldyticks, oldyticks[-1]*2-oldyticks[-2])
    ax.set_yticks(newyticks)
    if max(newyticks) < 10:
        yticklabels = ["{:1.1f}".format(s) for s in newyticks]
    else:
        yticklabels = ["{}".format(int(s)) for s in newyticks]
    ax.set_yticklabels(yticklabels, size=fontsize)
    ax.set_ylim(0,ymax)
    
    if type(sleuthDF)==pd.core.frame.DataFrame:
        q = sleuthDF.loc[sleuthDF["target_id"]==target, "qval"].values[0]
        if q<0.001: alpha = 0.001
        elif q<0.01: alpha = 0.01
        elif q<0.05: alpha = 0.05
        else: 
            alpha="NS"

    if twolinetitle:
        ax.set_title("{}\n{}".format(target, title), size=titlefontsize)
    else:
        ax.set_title("{} | {}".format(target, title), size=titlefontsize)
    anchored_text = AnchoredText(alphadict[alpha], loc=1, frameon=False, 
                                 prop={"size":fontsize, "weight":"bold"})
    ax.add_artist(anchored_text)

def nearest_square(n):
    return int(np.sqrt(round(np.sqrt(n))**2))

def multiTPMplot(targets, countsDF, countCol, sleuthDF, shape="maxcols", fontsize=10, titlefontsize=12, save=False, path=None, dpi=100, transparent=False, subplot_w=4, subplot_h=4, xticks='full', twolinetitle=False):
    """Creates a set of subplots of the genes using tpm_plot"""
    
    if type(sleuthDF)==pd.core.frame.DataFrame: targets = targets.values
    targetDict = dict(zip(targets[:,0], targets[:,1]))
    
    expressed = sleuthDF[sleuthDF["target_id"].isin(targetDict.keys())]["target_id"].values
    
    if shape == "square": ncols, nrows = (nearest_square(len(expressed)), nearest_square(len(expressed)))
    elif shape == "maxcols": 
        if len(expressed) > 4:
            ncols = 4
        else:
            ncols = len(expressed)
        nrows = int(np.ceil(len(expressed)/4))
    else: ncols, nrows = shape
    
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols*subplot_w, nrows*subplot_h), sharey=False, sharex=False)
    if len(expressed) == 0:
        print("Warning: no targets were expressed")
        return
    elif len(expressed) == 1:
        tpm_plot(countsDF=countsDF, column=countCol, target=expressed[0], ax=axes, title=targetDict[expressed[0]], sleuthDF=sleuthDF, fontsize=fontsize, titlefontsize=titlefontsize, xticks=xticks, twolinetitle=twolinetitle)
    else:
        i = 1
        for i,ax in enumerate(axes.flat):
            if i < len(expressed):
                tpm_plot(countsDF=countsDF, column=countCol, target=expressed[i], ax=ax, title=targetDict[expressed[i]], fontsize=fontsize, titlefontsize=titlefontsize, sleuthDF=sleuthDF, xticks=xticks, twolinetitle=twolinetitle)
            else:
                ax.axis("off")
            i+=1
            
    plt.tight_layout()
    if save:
        plt.savefig(path, dpi=dpi, transparent=transparent, bbox_inches="tight")

def lighten_color(color, amount=0.5):
    """
    Lightens (or darkerns) the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.
    
    This function was written by Ian Hincks:
    https://gist.github.com/ihincks/6a420b599f43fcd7dbd79d56798c4e5a
    https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
    
    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

flatten = lambda t: [item for sublist in t for item in sublist]

def remove_singletons(graph):
    keep = []
    for k,v in graph.degree(graph.nodes):
        if v>0: keep.append(k)
    return graph.subgraph(keep)

def get_plot_height(ax, expansion=1.0):
    plot_height = 0
    for line in ax.lines:
        if max(line.get_ydata()) > plot_height: plot_height = max(line.get_ydata())
        
    return plot_height*expansion

class Vividict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)() # retain local pointer to value
        return value                     # faster to return than dict lookup