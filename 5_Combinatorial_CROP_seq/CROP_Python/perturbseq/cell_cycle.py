# Perturbseq library for loading and manipulating single-cell experiments
# Copyright (C) 2019  Thomas Norman

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

import pandas as pd
import numpy as np
from collections import OrderedDict
import seaborn as sns

def group_corr(population, gene_list):
    """Measures the correlation of all genes within a list to the average expression of all genes within that
    list (used for cell cycle position calling)
    
    Args:
        population: CellPopulation to pull expression from
        gene_list: list of gene names
        
    Returns:
        Correlation coefficient of each gene with the mean expression of all
    """
    # returns list of correlations of each gene within a list of genes with the total expression of the group
    expression_matrix = population.where(genes='gene_name in @gene_list', gene_names=True, gene_list=gene_list)
    expression_matrix.loc[:,'total'] = expression_matrix.mean(axis=1)
    return expression_matrix.corr()['total'][:-1] # delete last element which is total

def refine_gene_list(population, gene_list, threshold, return_corrs=False):
    """Refines a list of genes by removing those that don't correlate well with the average expression of 
    those genes
    
    Args:
        population: CellPopulation to pull expression from
        gene_list: list of gene names
        threshold: threshold on correlation coefficient used to discard genes (expression of each gene is
            compared to the bulk expression of the group and any gene with a correlation coefficient less
            than this is discarded)
        return_corrs: whether to return the correlations along with the gene names (default: False)
        
    Returns:
        Refined list of genes that are well correlated with the average expression trend
    """
    corrs = group_corr(population, gene_list)    
    if (return_corrs):
        return corrs[corrs >= threshold].reset_index()
    else:
        return pd.Series(corrs[corrs >= threshold].index.values)
    
def group_score(population, gene_list):
    """Scores cells within population for expression of a set of genes. Raw expression data are first 
    log transformed, then the values are summed, and then scores are Z-normalized across all cells.
    
    Args:
        population: CellPopulation to pull expression from
        gene_list: list of gene names    
        
    Returns:
        Z-scored expression data
    """
    expression_matrix = population.where(genes='gene_name in @gene_list', gene_names=True, gene_list=gene_list)
    scores = expression_matrix.apply(lambda x: np.log2(x+1)).sum(axis=1)
    scores = (scores - scores.mean())/scores.std()
    return scores

def batch_group_score(population, gene_lists):
    """Scores cells within population for expression of sets of genes. Raw expression data are first 
    log transformed, then the values are summed, and then scores are Z-normalized across all cells. 
    Returns an OrderedDict of each score.
    
    Args:
        population: CellPopulation to pull expression from
        gene_lists: list of lists of gene names    
    """
    batch_scores = OrderedDict()
    for gene_list in gene_lists:
        batch_scores[gene_list] = group_score(population, gene_lists[gene_list])
    return batch_scores

def get_cell_phase_genes(pop, refine=False, threshold=0):
    """Returns a list of cell-cycle-regulated marker genes, filtered for coherence
    
    Args:
        refine: whether to refine the gene lists based on how consistent the expression is among 
            the groups
        threshold: threshold on correlation coefficient used to discard genes (expression of each 
            gene is compared to the bulk expression of the group and any gene with a correlation 
            coefficient less than this is discarded)
            
    Returns:
        a list of cell-cycle-regulated marker genes that show strong co-expression
    """
    cell_phase_genes = OrderedDict()
    cell_phase_genes['G1-S'] = pd.Series(['ARGLU1', 'BRD7', 'CDC6', 'CLSPN', 'ESD', 'GINS2', 
                                          'GMNN', 'LUC7L3', 'MCM5', 'MCM6', 'NASP', 'PCNA', 
                                          'PNN', 'SLBP', 'SRSF7', 'SSR3', 'ZRANB2'])
    cell_phase_genes['S'] = pd.Series(['ASF1B', 'CALM2', 'CDC45', 'CDCA5', 'CENPM', 'DHFR',
                                        'EZH2', 'FEN1', 'HIST1H2AC', 'HIST1H4C', 'NEAT1',
                                        'PKMYT1', 'PRIM1', 'RFC2', 'RPA2', 'RRM2', 'RSRC2',
                                        'SRSF5', 'SVIP', 'TOP2A', 'TYMS', 'UBE2T', 'ZWINT'])
    cell_phase_genes['G2-M'] = pd.Series(['AURKB', 'BUB3', 'CCNA2', 'CCNF', 'CDCA2', 'CDCA3',
                                           'CDCA8', 'CDK1', 'CKAP2', 'DCAF7', 'HMGB2', 'HN1',
                                           'KIF5B', 'KIF20B', 'KIF22', 'KIF23', 'KIFC1', 'KPNA2',
                                           'LBR', 'MAD2L1', 'MALAT1', 'MND1', 'NDC80', 'NUCKS1',
                                           'NUSAP1', 'PIF1', 'PSMD11', 'PSRC1', 'SMC4', 'TIMP1',
                                           'TMEM99', 'TOP2A', 'TUBB', 'TUBB4B', 'VPS25'])
    cell_phase_genes['M'] = pd.Series(['ANP32B', 'ANP32E', 'ARL6IP1', 'AURKA', 'BIRC5', 'BUB1',
                                        'CCNA2', 'CCNB2', 'CDC20', 'CDC27', 'CDC42EP1', 'CDCA3',
                                        'CENPA', 'CENPE', 'CENPF', 'CKAP2', 'CKAP5', 'CKS1B',
                                        'CKS2', 'DEPDC1', 'DLGAP5', 'DNAJA1', 'DNAJB1', 'GRK6',
                                        'GTSE1', 'HMG20B', 'HMGB3', 'HMMR', 'HN1', 'HSPA8',
                                        'KIF2C', 'KIF5B', 'KIF20B', 'LBR', 'MKI67', 'MZT1',
                                        'NUF2', 'NUSAP1', 'PBK', 'PLK1', 'PRR11', 'PSMG3', 'PWP1',
                                        'RAD51C', 'RBM8A', 'RNF126', 'RNPS1', 'RRP1', 'SFPQ',
                                        'SGOL2', 'SMARCB1', 'SRSF3', 'TACC3', 'THRAP3', 'TPX2',
                                        'TUBB4B', 'UBE2D3', 'USP16', 'WIBG', 'YWHAH', 'ZNF207'])
    cell_phase_genes['M-G1'] = pd.Series(['AMD1', 'ANP32E', 'CBX3', 'CDC42', 'CNIH4', 'CWC15',
                                           'DKC1', 'DNAJB6', 'DYNLL1', 'EIF4E', 'FXR1', 'GRPEL1',
                                           'GSPT1', 'HMG20B', 'HSPA8', 'ILF2', 'KIF5B', 'KPNB1',
                                           'LARP1', 'LYAR', 'MORF4L2', 'MRPL19', 'MRPS2', 'MRPS18B',
                                           'NUCKS1', 'PRC1', 'PTMS', 'PTTG1', 'RAN', 'RHEB', 'RPL13A',
                                           'SRSF3', 'SYNCRIP', 'TAF9', 'TMEM138', 'TOP1', 'TROAP',
                                           'UBE2D3', 'ZNF593'])

    if (refine):
        for phase in cell_phase_genes:
            cell_phase_genes[phase] = refine_gene_list(pop, cell_phase_genes[phase], threshold)

    return cell_phase_genes

def get_cell_phase(pop, gene_list=None, refine=False, threshold=0):
    """Compute cell cycle phase scores for cells in the population
    
    Args:
        gene_list: OrderedDict of marker genes to use for cell cycle phases. If None, the default 
            list will be used.
        refine: whether to refine the gene lists based on how consistent the expression is among 
            the groups
        threshold: threshold on correlation coefficient used to discard genes (expression of each 
            gene is compared to the bulk expression of the group and any gene with a correlation 
            coefficient less than this is discarded)

    Returns:
        Cell cycle scores indicating the likelihood a given cell is in a given cell cycle phase 
    """
    # get list of genes if one is not provided
    if gene_list is None: 
        cell_phase_genes = get_cell_phase_genes(pop, refine=True, threshold=0.3)
    else:
        cell_phase_genes = gene_list
    
    # score each cell cycle phase and Z-normalize
    phase_scores = pd.DataFrame(batch_group_score(pop, cell_phase_genes))
    normalized_phase_scores = phase_scores.sub(phase_scores.mean(axis=1), axis=0).div(phase_scores.std(axis=1), axis=0)

    normalized_phase_scores_corr = normalized_phase_scores.transpose()
    normalized_phase_scores_corr['G1-S'] = [1, 0, 0, 0, 0]
    normalized_phase_scores_corr['S'] = [0, 1, 0, 0, 0]
    normalized_phase_scores_corr['G2-M'] = [0, 0, 1, 0, 0]
    normalized_phase_scores_corr['M'] = [0, 0, 0, 1, 0]
    normalized_phase_scores_corr['M-G1'] = [0, 0, 0, 0, 1]

    phase_list = ['G1-S', 'S', 'G2-M', 'M', 'M-G1']

    # final scores for each phaase are correlation of expression profile with vectors defined above
    cell_cycle_scores = normalized_phase_scores_corr.corr()[-len(phase_list):].transpose()[:-len(phase_list)]
    
    # pick maximal score as the phase for that cell
    cell_cycle_scores['cell_cycle_phase'] = cell_cycle_scores.idxmax(axis=1)
    cell_cycle_scores['cell_cycle_phase'] = cell_cycle_scores['cell_cycle_phase'].astype('category')
    cell_cycle_scores['cell_cycle_phase'].cat.set_categories(phase_list, inplace=True)

    def progress_ratio(x, phase_list):
        ind = phase_list.index(x['cell_cycle_phase'])
        return x[phase_list[(ind - 1) % len(phase_list)]] - x[phase_list[(ind + 1) % len(phase_list)]]

    # interpolate position within given cell cycle phase
    cell_cycle_scores['cell_cycle_progress'] = cell_cycle_scores.apply(lambda x: progress_ratio(x, list(phase_list)), axis=1)
    cell_cycle_scores.sort_values(['cell_cycle_phase', 'cell_cycle_progress'],
                                  ascending=[True, False],
                                  inplace=True)

    # order of cell within cell cycle phase
    cell_cycle_scores['cell_cycle_order'] = cell_cycle_scores.groupby('cell_cycle_phase').cumcount()
    cell_cycle_scores['cell_cycle_order'] = cell_cycle_scores.groupby('cell_cycle_phase')['cell_cycle_order'].apply(lambda x: x/(len(x) - 1))

    return cell_cycle_scores

def add_cell_cycle_scores(pop, gene_list=None):
    """Call cell cycle positions for cells within the population. If more direct control is desired,
    use get_cell_phase.
    
    Args:
        gene_list: OrderedDict of marker genes to use for cell cycle phases. If None, the default 
            list will be used.
    """
    phase_list = ['G1-S', 'S', 'G2-M', 'M', 'M-G1']
    
    cell_cycle_scores = get_cell_phase(pop, refine=True, gene_list=gene_list, threshold=0.2)
    pop.add_property(cells=cell_cycle_scores)
    pop.cells['cell_cycle_phase'] = cell_cycle_scores['cell_cycle_phase'].astype('category')
    pop.cells['cell_cycle_phase'].cat.set_categories(phase_list, inplace=True)
    
def cell_cycle_position_heatmap(pop, cells=None, **kwargs):
    """Plot a heatmap of cells ordered by cell cycle position
    
    Args:
        pop: CellPopulation instance
        cells: query string for cell properties (i.e. executes pop.cells.query(cells=cells))
        **kwargs: all other keyword arguments are passed to pop.where
    """
    if cells is None:     
        cell_cycle_scores = pop.cells[['G1-S', 'S', 'G2-M', 'M', 'M-G1','cell_cycle_phase', 'cell_cycle_progress']].dropna()
    else:
        celllist = pop.where(cells=cells, **kwargs).index
        cell_cycle_scores = pop.cells.loc[celllist][['G1-S', 'S', 'G2-M', 'M', 'M-G1','cell_cycle_phase', 'cell_cycle_progress']].dropna()

    cell_cycle_scores.sort_values(['cell_cycle_phase', 'cell_cycle_progress'],
                                  ascending=[True, False],
                                  inplace=True)
    ax = sns.heatmap(cell_cycle_scores[cell_cycle_scores.columns[:-2]].transpose(), annot=False, xticklabels=False, linewidths=0)

