"""
mplrc.py

matplotlib rc params and axes rectangles to generate figures of appropriate
size for different types of publication.

see: http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples
"""
from numpy import mod
#import pyplot as plt

dashes = ['--', #    : dashed line
          '-', #     : solid line
          '-.', #   : dash-dot line
          ':', #    : dotted line
          '-']

colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546"]


# Generate a linestyle based on a list
def fdashes(i):
    return dashes[mod(i,len(dashes))]

# Generate a color based on a list
def fcolors(i):
    return colors[mod(i,len(colors))]



from math import sqrt
fig_width_pt = 246.0                    # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (sqrt(5.0)-1.0)/2.0       # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
#aps = {'params': {'axes.labelsize': 10,
#                  'text.fontsize': 10,
#                  'legend.fontsize': 10,
#                  'xtick.labelsize': 8,
#                  'ytick.labelsize': 8,
#                  'font.family': 'serif',
#                  'font.serif': 'Computer Modern Roman',
#                  'test.usetex': True,
#                  'figure.figsize': fig_size,
#                  'xtick.major.size': 4,
#                  'xtick.minor.size': 2,
#                  'xtick.major.pad': 4,
#                  'xtick.minor.pad': 4,
#                  'ytick.major.size': 4,
#                  'ytick.minor.size': 2,
#                  'ytick.major.pad': 4,
#                  'ytick.minor.pad': 4,
#                  'axes': [0.13,0.2,0.95-0.13,0.95-0.2]}}
#plt.rcParams.update(aps['params'])

aps = {'params': {'axes.labelsize': 10,
                  'text.fontsize': 10,
                  'legend.fontsize': 10,
                  'xtick.labelsize': 8,
                  'ytick.labelsize': 8,
                  'font.family': 'serif',
                  'font.size': 10,
                  'font.serif': 'Computer Modern Roman',
                  'test.usetex': True,
                  'figure.figsize': fig_size,
                  #'xtick.major.size': 10,
                  #'xtick.minor.size': 5,
                  #'xtick.major.pad': 10,
                  #'xtick.minor.pad': 10,
                  #'ytick.major.size': 10,
                  #'ytick.minor.size': 5,
                  #'ytick.major.pad': 10,
                  #'ytick.minor.pad': 10,
                  'axes': [0.13,0.2,0.95-0.13,0.95-0.2]}}
#plt.rcParams.update(aps['params'])
