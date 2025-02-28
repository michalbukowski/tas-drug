# https://github.com/michalbukowski/tas-drug
# (c) 2023-2025 Michał Bukowski (michal.bukowski@tuta.io, m.bukowski@uj.edu.pl)
# Department of Analytical Biochemistry, Faculty of Biochemistry, Biophysics and Biotechnology,
# Jagiellonian University, Krakow, Poland
# Distributed under GPL-3.0 license
# This notebook is a part of the TAS-Drug repository and is associated with the following study:
# Bukowski M, Banasik M, Chlebicka K, Bednarczyk K, Bonar E, Sokołowska D, Żądło T,
# Dubin G, Władyka B. (2025) Analysis of co-occurrence of type II toxin–antitoxin systems
# and antibiotic resistance determinants in Staphylococcus aureus. mSystems 0:e00957-24.
# https://doi.org/10.1128/msystems.00957-24

import numpy as np
import pandas as pd
from collections import defaultdict
from scipy.stats import hypergeom
from matplotlib import pyplot as plt

def fetch_seqids(fpath: str) -> list[str]:
    '''Scans a FASTA format file for IDs of sequence records present in that file.
       Arguments:
       fpath: str        -- a path to a FASTA format file with sequences
       Returns:
       seqids: list[str] -- a list with IDs for all sequences found in the file
    '''
    
    seqids = []
    f = open(fpath)
    for line in f:
        if line[0] == '>':
            seqid = line[1:line.find(' ')]
            seqid = seqid.rsplit('|', 1)[-1]
            seqids.append(seqid)
    f.close()
    return seqids


def sf_10(sample_size: int, total_cases: int, total_size: int) -> float:
    '''Calculates a probability value for hypergeometric distribution survival function (cumulative
       probability) for 10 or more occurrences of a given class cases in a sample of `sample_size`,
       assuming the general population of `total_size` and the number of `total_cases`
       of the class cases in that population.
       Arguments:
       sample_size: int -- the sample size.
       total_cases: int -- the total count of a given class cases in the general population.
       total_size:  int -- the size of the general population.
       Returns:
       prob: float -- hypergeometric distribution survival function value.
    '''
    
    prob = hypergeom.sf(10-1, total_size, total_cases, sample_size)
    
    return prob


def hyper(sample_cases: int, sample_size: int, total_cases: int, total_size: int) -> float:
    '''Calculates a probability value for hypergeometric cumulative distribution function (CDF)
       or survival function (SF), whichever lower, for `sample_cases` or, respectively,
       less or more occurrences of a given class cases in a sample of `sample_size`,
       assuming the general population of `total_size` and the number of `total_cases`
       of the class cases in that population. The returned value is negative for CDF
       (under-representation) and positive for SF (over-representation).
       Arguments:
       sample_cases: int -- the total count of a given class cases in the sample.
       sample_size:  int -- the sample size.
       total_cases:  int -- the total count of a given class cases in the general population.
       total_size:   int -- the size of the general population.
       Returns:
       prob: float -- the value of a corresponding cumulative probability function.
    '''
    
    prob_1 = hypergeom.cdf(sample_cases,    total_size, total_cases, sample_size)
    prob_2 = hypergeom.sf(sample_cases - 1, total_size, total_cases, sample_size) 
    
    if prob_1 <= prob_2: 
        return -1*prob_1
    else: 
        return prob_2


def get_color(ratio: float) -> np.ndarray[ [float, float, float] ]:
    '''Maps a normalised ratio value (in a range -1.0 to 1.0) to a blue-red color scale.
       Arguments:
       ratio: float -- a normalised ration value from a -1.0 to 1.0 range.
       Returns:
       color: np.ndarray[ [float, float, float] ] -- an Numpy Array with RGB color values
                                                     mapped to 0.0 to 1.0 range.
    '''
    
    if ratio == 0.0:
        color = np.array([1.0, 1.0, 1.0])
    elif ratio > 0.0:
        color = np.array([1.0, 1.0-ratio, 1.0-ratio])
    else:
        color = np.array([1.0+ratio, 1.0+ratio, 1.0])
        
    return color

    
def calc_ratios(fin_df: pd.DataFrame) -> pd.DataFrame:
    '''Given a Pandas DataFrame adds new `ratio` column with values of ratio of ratios.
       The first ratio Pandas Series is calculated for the `corr` and `ta_count` columns,
       the second for the `drug_count` and `tot` columns. The final Pandas Series
       of ratio of ratio values is calculated by dividing the the first ratio by the other
       (positive fold value). However, if the other is grater than the first,
       the values are swapped and expressed as negative (negative fold value).
       Importantly, the other ratio must never be 0 since the general `drug_count` of a given
       drug resistance determinant and the value of population of size `tot`
       are never equal to 0. The general `ta_count` is also never equal to 0. It is a consequence
       of the fact that only existing TA systems and drug resistance determinants
       in the general non-empty population are analysed, thus all those 3 values must never be 0.
       However the value of `corr` that expresses how many times a given TA system co-occurs with
       a given drug resistance determinant might be 0. In such a case `ratio` and `log2_ratio`
       values will remain `-np.inf` (negative infinity) and will be dealt with during a subsequent
       stage of the analysis.
       Another added column `log2_ratio` contains log2 values of such ratios of ratios.
       Arguments:
       fin_df: pd.DataFrame -- the input Pandas DataFrame with `corr`, `ta_count`, `drug_count`
                               and `tot` columns containing integer values.
       Returns:
       fin_df: pd.DataFrame -- a copy of the input Pandas DataFrame with the new `ratio`
                               and `log2_ratio` columns.
    '''
    
    fin_df = fin_df.copy()
    # fraction of a given TA-carriers with a given resistance determinant
    # to the general count of TA-carriers
    first_ratio = ( fin_df['corr'] /  fin_df['ta_count'] )
    # the general count of carriers of a given drug resistance determinant
    # to the size of the general population
    other_ratio = ( fin_df['drug_count'] / fin_df['tot'] )
    
    positive = first_ratio >= other_ratio
    fin_df.loc[ positive, 'ratio'] =   first_ratio / other_ratio
    fin_df.loc[~positive, 'ratio'] = - other_ratio / first_ratio
    
    fin_df.loc[fin_df['ratio'] > 0, 'log2_ratio'] =   np.log2(  fin_df.loc[fin_df['ratio'] > 0, 'ratio'])
    fin_df.loc[fin_df['ratio'] < 0, 'log2_ratio'] = - np.log2(- fin_df.loc[fin_df['ratio'] < 0, 'ratio'])
    
    return fin_df


def unstack_ratios(
    fin_df: pd.DataFrame, ta_count: defaultdict[str, int], drug_count: defaultdict[str, int],
    pval10_th: float, ratio_th: float, ratio_max: float, adjust: bool
) -> tuple[pd.DataFrame, pd.DataFrame, float]:
    '''Given a Pandas DataFrame with the `log2_ratio` and `pval_10` columns, creates
       two new DataFrames, the rows of which refere to TA systems, and the columns,
       to drug resistance determinants. The values of the new DataFrames denotes, respectively,
       the `log2_ratio` column values or those values mapped to a color scale. Before creation,
       the input DataFrame is filtered in respect to `pval_10` column values and `-np.inf` values
       are replaced either by `ratio_max` value of its adjusted value. The `log2_ratio` values
       that do not meet the `pval_10` threshold (`pval10_th`) are set to 0 (not analysed further).
       If the adjustment is chosen (the `adjust` argument set to `True`), the `-np.inf` values
       will be replaced by `log2(abs_max)`. The `abs_max` is the absolute maximal value
       of the `log2_ratio` column without taking into the consideration the `-np.inf` values.
       Such na adjustment might be applied to avoid any overestimation of a negative co-occurrence
       ratio of ratios for cases in which a given TA systems does not occur at all with
       a given drug resistant determinant (the `corr` value equals to 0). Such situations
       are usually a consequence of low count of one genetic element or the other (or both),
       and the the absence of the co-occurrence might be a caused by non-ranodomness of genetic data
       instead of being a results of any underlying biological mechanism. It is recommended
       to decide whether to apply the adjustment and to select the `ratio_max` value after inspecting
       the results obtained without adjustment and the returned `abs_max` value.
       The finally selected `log2_values` that exceeds ratio_max or -ratio_max values are reduced to
       `ratio_max` value. This might be useful for visualisation purposes. The rows and columns
       of the final DataFrames are ordered by the count values of TA systems and
       drug resistance determinants in the general population in a descending manner.
       Arguments:
       fin_df: pd.DataFrame              -- the input DataFrame that contains `log2_ratio` column.
       ta_count: defaultdict[str, int]   -- a dafaultdict mapping TA system name to its total count.
       drug_count: defaultdict[str, int] -- a dafaultdict mapping drug resistance determinant name
                                            to its total count.
       pval10_th: float                  -- a minimal threshold for the `pval_10` column values.
       ratio_th: float                   -- log_2 ratio minimal threshold values.
       ratio_max: float                  -- the limit of for the absolute of log2 ratio value.
                                            Values exceeding the threshold will be reduced to this
                                            threshold value.
       adjust: bool                      -- use (True) or not (False) the log2 adjustment for
                                            the `abs_max` before replacing with it `-np.inf` values
                                            of `log2_ratio`.
       Returns:
       ratio_df: pd.DataFrame -- the final DataFrame with adjusted `log2_ratio` values.
       color_df: pd.DataFrame -- the final DataFrame with color values (`color`) that
                                 are mapped to `log2_ratio` values.
       abs_max: float         -- the absolute maximal `log2_ratio` values (not accounting for
                                 -np.inf values) before any adjustment is done.
    '''
    
    fin_df = fin_df.copy()
    fin_df.loc[
        (fin_df['pval_10']  < pval10_th) | (fin_df['log2_ratio'].abs() < ratio_th), 'log2_ratio'
    ] = 0.0
    
    abs_max = fin_df.loc[fin_df['log2_ratio'] != -np.inf, 'log2_ratio'].abs().max()
    
    adj_max = np.log2(abs_max) if adjust else ratio_max
    fin_df.loc[fin_df['log2_ratio'] == -np.inf, 'log2_ratio'] = -adj_max
    
    fin_df.loc[fin_df['log2_ratio'] >  ratio_max,  'log2_ratio'] =  ratio_max
    fin_df.loc[fin_df['log2_ratio'] < -ratio_max,  'log2_ratio'] = -ratio_max
    
    fin_df['color'] = (fin_df['log2_ratio']/ratio_max).apply(get_color)
    ratio_df, color_df = (
        fin_df.set_index('ta drug'.split())[col].unstack()
        for col in 'log2_ratio color'.split()
    )
    ratio_df.sort_index(key=lambda ind: [ta_count[key] for key in ind  ],
                        ascending=False, inplace=True)
    ratio_df.sort_index(key=lambda ind: [drug_count[key] for key in ind],
                        ascending=False, axis=1, inplace=True)
    
    color_df = color_df.loc[ratio_df.index, ratio_df.columns]
    
    return ratio_df, color_df, abs_max


def final_filter(
    ratio_df: pd.DataFrame, remove: list[str], keep: list[str], filter_blank: bool
) -> pd.DataFrame:
    '''Performs final Pandas DataFrame filtering by removing rows and columns of given indices
       as well as those the sum values of which are equal to zero (empty rows and columns).
       Arguments:
       ratio_df: pd.DataFrame -- a Pandas DataFrame to be filtered.
       remove: list[str]      -- a list of indices names (labels) of rows and columns
                                 that are to be removed
       keep: list[str]        -- a list of indices names (labels) of rows and columns.
                                 that are to be kept even if other criteria are not met.
       filter_blank: bool     -- if true, remove empty rows and columns (the values of
                                 which sums up to zero).
       Returns:
       filtered_df: pd.DataFrame -- the final filtered DataFrame.
    '''
    
    filtered_df = ratio_df.loc[
        ~ratio_df.index.isin(remove),
        ~ratio_df.columns.isin(remove)
    ].copy()
    if filter_blank:
        filtered_df = filtered_df.loc[
            (filtered_df.sum(axis=1) != 0) | filtered_df.index.isin(keep),
            (filtered_df.sum()       != 0) | filtered_df.columns.isin(keep)
        ]
    return filtered_df


def plot_ratios(
    ratio_df: pd.DataFrame, color_df: pd.DataFrame,
    ta_count: defaultdict[str, int], drug_count: defaultdict[str, int],
    alt_names: dict[str, str], fig_width: float, ax_left: float, ax_right: float,
    ax_top: float, ax_bottom: float, bar_count: float, bar_scale: float,
    bar_left: float, bar_bottom: float, ratio_max: float
) -> tuple[plt.Figure, plt.Axes, plt.Axes]:
    '''Given Pandas DataFrames with ratio of ratio values and colors mapped to them,
       plots the final visualisation, a heatmap and a heatmap scale bar.
       Arguments:
       ratio_df: pd.DataFrame -- a Pandas DataFrame with ratio of ratios values.
       color_df: pd.DataFrame -- a Pandas DataFrame with color values mapped
                                 to the ratio of ratios values.
       ta_count: defaultdict[str, int]   -- a dafaultdict mapping TA systems names
                                            to their total count.
       drug_count: defaultdict[str, int] -- a dafaultdict mapping drug resistance determinants
                                            names to their total count.
       alt_names: dict[str, str]         -- a dict mapping rows and columns names (labels)
                                            to replacement ones that are used only
                                            for visualisation purposes
       fig_width: float  -- the width of the final figure in inches. The height is calculated
                            automatically to keep the heat map basic units (cells) square.
       ax_left: float    -- left spacing of the heatmap axis from the figure edge expressed
                            as a fraction of the total figure width.
       ax_right: float   -- right spacing of the heatmap axis from the figure edge expressed
                            as a fraction of the total figure width.
       ax_top: float     -- top spacing of the heatmap axis from the figure edge expressed
                            as a fraction of the total figure height.
       ax_bottom: float  -- bottom spacing of the heatmap axis from the figure edge expressed
                            as a fraction of the total figure height.
       bar_count: float  -- the count of basic units (cells) of the heatmap scale bar.
       bar_scale: float  -- general scaling factor for the heatmap scale bar.
       bar_left: float   -- left spacing of the heatmap scale bar axis from the figure edge
                            expressed as a fraction of the total figure width.
       bar_bottom: float -- bottom spacing of the heatmap scale bar axis from the figure edge
                            expressed as a fraction of the total figure height.
       ratio_max: float  -- maximal absolute value to be used for the scale bar. Should be the same
                            as one used for scaling the data in the input color_df DataFrame.
       Returns:
       fig: plt.Figure  -- the reference to the final Figure.
       ax: plt.Axes     -- the reference to the heatmap Axes.
       bar_ax: plt.Axes -- the reference to the heatmap scale bar Axes.
    '''
    
    data = np.stack(
        color_df.to_numpy().reshape(color_df.shape[0]*color_df.shape[1])
    ).reshape( (color_df.shape[0], color_df.shape[1], 3) )

    fig_height = ax_bottom + ax_top + (fig_width-ax_left-ax_right) \
                 * data.shape[0] / data.shape[1]

    ax_pwidth  = 1.0-(ax_left+ax_right)/fig_width
    ax_pheight = 1.0-(ax_top+ax_bottom)/fig_height
    ax_pleft   = ax_left/fig_width
    ax_pbottom = ax_bottom/fig_height

    fig = plt.figure(figsize=(fig_width, fig_height), facecolor='white', dpi=150)
    ax = fig.add_axes([ax_pleft, ax_pbottom, ax_pwidth, ax_pheight])

    for y in range(data.shape[0]+1):
        ax.hlines(-0.5+y, -0.5, data.shape[1], linewidth=0.1, color='black')
    for x in range(data.shape[1]+1):
        ax.vlines(-0.5+x, -0.5, data.shape[0], linewidth=0.1, color='black')

    ax.imshow(data, aspect='auto')

    ratios = ratio_df.abs().to_numpy()
    
    indices = np.nonzero(np.abs(ratios) != 0)
    for row, col in zip(*indices):
        color = 'black' if ratios[row, col] <= ratio_max/2 + 1 else 'white'
        ax.text(col, row, f'{ratios[row, col]:0.1f}', ha='center', va='center', color=color, fontsize=11.0)
        
    indices = np.nonzero(np.abs(ratios) == 0)
    for row, col in zip(*indices):
        ax.scatter(col, row, s=0.75, linewidth=0.0, color='lightgray')
    
    get_name = lambda drug_name: alt_names.get(drug_name, drug_name).replace('_', '\_')
    ta_labels   = [ f'{alt_names.get(ta_name, ta_name)} ({ta_count[ta_name]})'
                    for ta_name   in color_df.index   ]
    drug_labels = [ f'$\it {get_name(drug_name)}$ ({drug_count[drug_name]})'
                    for drug_name in color_df.columns ]

    ax.set_xticks(range(data.shape[1]), drug_labels, fontsize=13.0, rotation=90)
    ax.set_yticks(range(data.shape[0]), ta_labels, fontsize=13.0)
    ax.spines['top bottom right left'.split()].set_visible(False)
    ax.tick_params(axis='both', which='both', length=0)
    ax.set_xlim((0-0.5, data.shape[1]-0.5))
    ax.set_ylim((data.shape[0]-0.5, 0-0.5))
    ax.set_title('Occurrence rate ratio of drug determinants and TA systems', fontsize=16)

    bar_pwidth  = ax_pwidth  / data.shape[1] * bar_count * bar_scale
    bar_pheight = ax_pheight / data.shape[0] * bar_scale

    bar_ax = fig.add_axes([bar_left, bar_bottom, bar_pwidth, bar_pheight])

    for y in range(2):
        bar_ax.hlines(-0.5+y, -0.5, data.shape[1], linewidth=0.1, color='black')
    for x in range(bar_count+1):
        bar_ax.vlines(-0.5+x, -0.5, data.shape[0], linewidth=0.1, color='black')

    bar_data = np.array(
        [ get_color(ratio) for ratio in np.linspace(-1, 1, bar_count) ]
    ).reshape(3*bar_count).reshape((1, bar_count, 3))

    bar_ax.spines['top bottom right left'.split()].set_visible(False)
    bar_ax.tick_params(axis='both', which='both', length=0)
    bar_ax.set_yticks([])
    bar_labels = [
        f'{ratio:d}' for ratio in np.linspace(
            -ratio_max, ratio_max, int( (bar_count-1)/2 ) + 1
        ).round(0).astype(int)
    ]
    bar_ticks  = np.linspace(0, bar_count-1, len(bar_labels)).round(0).astype(int)
    bar_ax.set_xticks(bar_ticks, bar_labels, fontsize=12.0, ha='center', va='top')
    bar_ax.imshow(bar_data, aspect='auto')
    bar_ax.set_title(f'Ratio (log2)', fontsize=14.0, pad=5.0)
    
    return fig, ax, bar_ax
