import matplotlib.pyplot as plt
import seaborn as sns


def set_style_paper():
    sns.set_context('paper')
    sns.set_style(style='ticks')
    fontsize = 7
    plt.rcParams['figure.figsize'] = [2.67, 2]
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['legend.fontsize'] = fontsize
    plt.rcParams['legend.title_fontsize'] = fontsize
    plt.rcParams['figure.titlesize'] = fontsize
    plt.rcParams['xtick.labelsize'] = fontsize
    plt.rcParams['ytick.labelsize'] = fontsize
    plt.rcParams['axes.titlesize'] = fontsize
    plt.rcParams['axes.labelsize'] = fontsize
    plt.rcParams['svg.fonttype'] = 'none'


def set_style_talk():
    sns.set_context('talk')
    sns.set_style(style='ticks')
    plt.rcParams['svg.fonttype'] = 'none'


# palettes
saturation = 0.75
colors = 8
palette = sns.color_palette('viridis_r', colors)
palette_qual = sns.color_palette(['#b2df8a', '#1f78b4', '#a6cee3'])
