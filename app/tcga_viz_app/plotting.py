import numpy as np
from itertools import product
import plotly.graph_objs as go
import plotly.figure_factory as ff
from plotly import tools

umap_scene = {
    'domain': {'x': [0.0, 1.0], 'y': [0.0, 1.0]},
    'xaxis': {'range': [-15, 15], 'nticks': 7},
    'yaxis': {'range': [-15, 15], 'nticks': 7},
    'zaxis': {'range': [-15, 15], 'nticks': 7},
    'aspectratio': dict(x=1, y=1, z=1),
   }
tsne_scene = {
    'domain': {'x': [0.0, 1.0], 'y': [0.0, 1.0]},
    'xaxis': {'range': [-60, 60], 'nticks': 7},
    'yaxis': {'range': [-60, 60], 'nticks': 7},
    'zaxis': {'range': [-60, 60], 'nticks': 7},
    'aspectratio': dict(x=1, y=1, z=1),
   }

default_layout = {
    'scene': umap_scene,
    #'title': 'TCGA Cancer (sub)types',
    'hovermode': 'closest',
    'height': 680,
    'margin': dict(l=0, r=0, b=0, t=0)
}

color_cols = {'tumor.type': 'hue',
              'organ_system': 'shade',
              'subtype_full': 'subtype_hue'}

def generate_dr_figure(subset_data, color_col, old_layout=None, dr_type='umap'):
    data = []
    if dr_type == 'umap':
        scene = umap_scene
    else:
        scene = tsne_scene

    # data will come in grouped by:
    # match, tumor.type, organ_system, subtype
    color_col = color_cols[color_col]
    for (match, tt, organ, subtype), subset in subset_data:
        color = subset[color_col]

        if dr_type == 'umap':
            data_vals = dict(x=subset['umap_0'],
                             y=subset['umap_1'],
                             z=subset['umap_2'])
        else:
            data_vals = dict(x=subset['tsne_0'],
                             y=subset['tsne_1'],
                             z=subset['tsne_2'])

        marker = 'circle' if match == 'tumor' else 'cross'#star-square'

        text_match = match.capitalize()
        hovertext = ""
        if subtype != 'None':
            hovertext += f"Tumor type: {subtype}"
        else:
            hovertext += f"Tumor type: {tt}"
        hovertext += f"<br>Organ system: {organ}<br>{text_match} sample"

        trace = go.Scatter3d(
                    visible=True,
                    mode='markers', name="",
                    marker=dict(color=color, size=3, symbol=marker,
                                line=dict(width=0, color=color),),
                    hoverinfo='text', hovertext=[hovertext]*len(subset),
                    showlegend=False,
                    **data_vals)

        data.append(trace)

    if old_layout == None:
        layout = default_layout.copy()
        layout.update(scene=scene)
    else:
        layout = old_layout
        layout['scene'].update(**scene)

    return go.Figure(data=data, layout=layout)


def marginal_heatmap(data_m, #data_x, data_y,
                     colors_x, colors_y, #colors_x_bar, colors_y_bar,
                     marg_xlabel='Samples per class',
                     marg_ylabel='Pos Pred Value',
                     title=None, margins=dict(l=50, r=50, b=100, t=100),
                     size=(1000, 700)):

    normed = data_m / data_m.sum(axis=0)
    normed.fillna(0, inplace=True)

    colorscale = [[0, '#ffffff'], [0.00001, '#fcfcfc'], [1, '#333333']]
    font_colors = ['#222222', '#222222', '#eeeeee']

    normed = normed.iloc[::-1]
    data_m = data_m.iloc[::-1]
    hover_text = [[f'Reference: {y}<br>Predicted: {x}<br>Count: {data_m.values[i,j]}'
                   for j, y in enumerate(normed.columns)]
                  for i, x in enumerate(normed.index)]

    heatmap_params = dict(z=normed.values,
                          y=list(normed.index), x=list(normed.columns),
                          annotation_text=data_m.values,
                          colorscale=colorscale, font_colors=font_colors,
                          y0=normed.columns[0],
                          text=hover_text, hoverinfo='text',
                          )
    figure = ff.create_annotated_heatmap(**heatmap_params)

    # add colorscale
    for i in range(len(figure['data'])):
        figure['data'][i]['xaxis'] = 'x1'
        figure['data'][i]['yaxis'] = 'y1'

    if size is not None:
        figure['layout'].update(width=size[0], height=size[1], autosize=False)
    if title is not None:
        figure['layout']['annotations'] = []
        figure['layout']['annotations'].append({'text': title,
                                                'showarrow': False,
                                                'ax':0, 'ay':0,
                                                'x':0, 'y':1,
                                                'font': {'size': 16}})
        figure['layout']['title'] = title
    if margins:
        figure['layout']['margin'] = go.Margin(**margins)
    figure['layout']['xaxis1'].update(**{'domain': [0.025, 1]})
    figure['layout']['xaxis2'] = {'domain': [0, 0.025]}
    figure['layout']['yaxis2'] = {'anchor': 'x2'}

    return figure
    #ax_x.bar(np.arange(data_x.shape[0]) + 0.5, data_x,
    #         width=0.9, color=colors_x)
    #ax_y.barh(np.arange(data_y.shape[0]) + 0.5, data_y,
    #          height=0.9, color=colors_y)

    #for k, color in enumerate(colors_x_bar):
    #    ax_h1.barh(0, 1, left=k, color=color, height=1)
    #    ax_h2.barh(0, 1, left=k, color=color, height=1)
    #for k, color in enumerate(colors_y_bar):
    #    ax_v1.bar(0, 1, bottom=k, color=color, width=1)
    #    ax_v2.bar(0, 1, bottom=k, color=color, width=1)
    #for (i, ci), (j, cj) in product(enumerate(colors_x_bar), enumerate(colors_y_bar)):
    #    if ci != cj: continue
    #    patch = mpatches.Rectangle([i, j], 1, 1, color=ci,
    #                               alpha=0.2, lw=0)
    #    ax_m.add_artist(patch)

    #transfer_ticks(ax_m, ax_h1, which='x', rotation=90)
    #transfer_ticks(ax_m, ax_v1, which='y')

    #fix_spines(ax_h1, [], keep_ticks=True)
    #fix_spines(ax_v1, [], keep_ticks=True)
    #fix_spines(ax_h2, [], keep_ticks=False)
    #fix_spines(ax_v2, [], keep_ticks=False)
    #fix_spines(ax_x, ['left'], keep_ticks=False)
    #fix_spines(ax_y, ['bottom'], keep_ticks=False)

    #ax_h1.set_xlabel('Reference')
    #ax_v1.set_ylabel('Prediction')

    #ax_x.set_ylabel(marg_xlabel)
    #ax_y.set_xlabel(marg_ylabel)
