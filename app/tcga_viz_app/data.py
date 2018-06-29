import os
import numpy as np
import pandas as pd

from matplotlib.colors import rgb2hex, hex2color
from colorspacious import cspace_convert


def desaturate(hex_color, n):
    if n == 1:
        return hex_color
    cspace = 'CIELab'
    dim = 0
    jch_color = cspace_convert([hex2color(hex_color)], 'sRGB1', cspace)
    desats = np.repeat(jch_color, n, axis=0)
    orig_sat = jch_color[0, dim]
    desats[:, dim] = orig_sat * np.linspace(0.25, 1, n)
    return np.clip(cspace_convert(desats, cspace, 'sRGB1'), 0, 1)


#dropdown_headers = ['Tumor type', 'Cancer type', 'Cancer subtype']
#dropdown_columns = ['tumor.type.matched', 'cancertype', 'subtype']
dropdown_headers = ['Tumor type', 'Organ system', 'Cancer subtype']
#dropdown_columns = ['tumor.type.matched', 'organ_system', 'subtype']
dropdown_columns = ['tumor.type', 'organ_system', 'subtype_full']

def read_color_data():
    data_file = os.path.join(os.path.dirname(__file__),
                             'data/cancer_data.csv')
    meta_data = pd.read_csv(data_file, header=0, index_col=0)

    return meta_data

def read_data():
    data_file = os.path.join(os.path.dirname(__file__),
                             'data/sample_dr.csv')
    sample_data = pd.read_csv(data_file, header=0, index_col=0)
    sample_data.fillna('None', inplace=True)

    sample_types = sample_data['sample.type'].astype(int).values
    tumor_mask = sample_types < 11
    sample_data.loc[tumor_mask, 'tumor.type.matched'] =\
        sample_data.loc[tumor_mask, 'tumor.type'] + ' T'
    sample_data.loc[~tumor_mask, 'tumor.type.matched'] =\
        sample_data.loc[~tumor_mask, 'tumor.type'] + ' N'
    sample_data.loc[tumor_mask, 'match'] = 'tumor'
    sample_data.loc[~tumor_mask, 'match'] = 'normal'

    sample_data['subtype_full'] = sample_data[['tumor.type', 'subtype']].apply(lambda x: ' - '.join(x), axis=1)
    none_mask = sample_data['subtype_full'].str.contains('None')
    sample_data.loc[none_mask, 'subtype_full'] = 'None'

    color_data = read_color_data()
    gray = '#e8e8e8'
    sample_data['hue'] = gray
    sample_data['shade'] = gray
    sample_data['subtype_hue'] = gray
    for abbrev, row in color_data.iterrows():
        inds = sample_data['tumor.type'] == abbrev
        sample_data.loc[inds, 'organ_system'] = row.organ_system.replace('_', ' ').capitalize()
        grouped = sample_data.loc[inds].groupby('subtype')
        n_subtypes = grouped.ngroups
        if n_subtypes > 1:
            subtype_hues = map(rgb2hex, desaturate(row.hue, n_subtypes))
            for (subtype, st_inds), h in zip(grouped, subtype_hues):
                sample_data.loc[st_inds.index, 'subtype_hue'] = h
        sample_data.loc[inds, 'shade'] = row.shade
        sample_data.loc[inds, 'hue'] = row.hue
    #sample_data.loc[~tumor_mask, 'shade'] = gray
    sample_data.loc[~tumor_mask, 'subtype_hue'] = gray

    return sample_data

#def get_tumortype_checkboxes(dataframe):
#    # matched normals are anything above 11
#    subdata = dataframe[['tumor.type', 'sample.type']]
#    checkboxes = {

def read_heatmap_data(fname):
    color_data = read_color_data()

    data_file = os.path.join(os.path.dirname(__file__), fname)
    table = pd.read_csv(data_file, header=0, index_col=0)
    table = sort_rectangular_table(clean_table(table), color_data.index)

    colors_x, colors_y = get_heatmap_colors(table, color_data, 'hue', other_col=None)
    colors_x_bar, colors_y_bar = get_heatmap_colors(table, color_data, 'shade', other_col=None)

    return (table, colors_x, colors_y)

def read_primary_ext_heatmap_data():
    return read_heatmap_data('data/primary_ext_contingency.csv')

def read_primary_heatmap_data():
    return read_heatmap_data('data/primary_cv_contingency.csv')

def read_subtype_heatmap_data(subtype):
    color_data = read_color_data()

    fname = f'data/{subtype}_contingency_table.csv'
    data_file = os.path.join(os.path.dirname(__file__), fname)
    table = pd.read_csv(data_file, header=0, index_col=0)
    table = clean_table(table)

    subtype_colors = desaturate(color_data.loc[subtype, 'hue'], len(table))

    return (table, subtype_colors, subtype_colors)


def get_available_checkboxes(dataframe):

    available_checkboxes = dict((col, sorted(dataframe[col].unique()))
                                for col in dropdown_columns)

    for col, checkboxes in available_checkboxes.items():
        if 'None' not in checkboxes: continue
        new_checkboxes = ['None'] + sorted(set(checkboxes) - set(['None']))
        available_checkboxes[col] = new_checkboxes

    for col, checkboxes in available_checkboxes.items():
        #lengths = [sum(dataframe[col] == cb) for cb in checkboxes]
        #new_checkboxes = [{'label': f'{cb}  ({length})', 'value': cb}
        #                  for cb, length in zip(checkboxes, lengths)]
        #available_checkboxes[col] = new_checkboxes
        new_checkboxes = [{'label': f'{cb}', 'value': cb} for cb in checkboxes]
        available_checkboxes[col] = new_checkboxes

    return available_checkboxes

def sort_rectangular_table(table, reference_sort):
    reference_sort = pd.Index(reference_sort)
    # sort by reference_sort
    cols = table.columns
    inds = table.index
    table = table.loc[reference_sort[reference_sort.isin(inds)],
                      reference_sort[reference_sort.isin(cols)]]

    cols = table.columns
    inds_in_cols = table.index[table.index.isin(cols)]
    inds_not_in_cols = table.index[~table.index.isin(cols)]

    inds = inds_in_cols.tolist() + inds_not_in_cols.tolist()
    return table.loc[inds, cols]

def clean_table(table, blacklist=[]):
    table = table.fillna(0)
    empty_rows = table.sum(axis=1) < 1e-2
    empty_cols = table.sum(axis=0) < 1e-2

    blacklist = table.columns.isin(blacklist)
    empty_cols = empty_cols | blacklist

    return table.loc[~empty_rows, ~empty_cols]

def get_heatmap_colors(table, reference, hue_or_shade='hue', other_col='organ_system'):
    if other_col:
        color_reference = reference.drop_duplicates('organ_system').set_index('organ_system')
    else:
        color_reference = reference

    x_colors = color_reference.loc[table.columns, hue_or_shade]
    y_colors = color_reference.loc[table.index, hue_or_shade]
    return x_colors, y_colors

def group_metrics(contingency, groupby, metrics=['pos pred value', 'sensitivity', 'specificity']):
    names = groupby.groups.keys()
    groups = list(groupby.groups.values())

    L = len(groups)
    cont_reduced = np.zeros((L, L), dtype=int)
    for i, j in product(range(L), repeat=2):
        a = groups[i]
        b = groups[j]
        ai = contingency.index.isin(a)
        bi = contingency.columns.isin(b)

        if (not ai.any()) or (not bi.any()): continue
        cont_reduced[i, j] = contingency.loc[ai, bi].values.sum()
    cont_reduced = clean_table(pd.DataFrame(cont_reduced, index=names, columns=names))

    ms = pd.DataFrame(0, index=cont_reduced.columns, columns=metrics)
    for k in range(cont_reduced.shape[1]):
        ak = cont_reduced.columns == cont_reduced.columns[k]
        bk = cont_reduced.index == cont_reduced.columns[k]
        # columns are reference, so [~ak, a] is FN and [ak, ~ak] is FP
        tp, fn, fp, tn = (cont_reduced.loc[bk, ak].values.sum(),
                          cont_reduced.loc[~bk, ak].values.sum(),
                          cont_reduced.loc[bk, ~ak].values.sum(),
                          cont_reduced.loc[~bk, ~ak].values.sum())
        #print(cont_reduced.columns[k], tp, fn, fp, tn)
        ms.iloc[k, :] = [metric_functions[metric](tp, fn, fp, tn)
                         for metric in metrics]
        #ms.iloc[k, :] = [tp / (tp + fp), tp / (tp + fn), tn / (tn + fp)]

    #print(cont_reduced.shape, ms.shape)
    return cont_reduced, clean_table(ms)
