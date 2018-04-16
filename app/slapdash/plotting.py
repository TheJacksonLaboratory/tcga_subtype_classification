import numpy as np
from matplotlib import cm
from matplotlib.colors import rgb2hex
import plotly.graph_objs as go
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
    'height': 700,
    'margin': dict(l=0, r=0, b=0, t=0)
}

# retinue of 12 + 15 = 27 colors
# Gray saved for Nones
COLORS = sum(zip(*zip(cm.tab20c.colors[:16:4],  cm.tab20c.colors[1:16:4],
                      cm.tab20c.colors[2:16:4], cm.tab20b.colors[::4],
                      cm.tab20b.colors[1::4],   cm.tab20b.colors[2::4])), ())
COLORS = list(map(rgb2hex, COLORS))
COLORS = COLORS + COLORS + COLORS

def generate_figure(sample_data, options, values, col):
    usable_colors = COLORS[:len(options)]

    data = []
    for k, item in enumerate(options):
        value = item.get('value')
        if value not in values: continue

        subdata = sample_data[sample_data[col] == value]
        color = usable_colors[k]
        if value == 'None': color = '#e8e8e8'

        trace_umap = go.Scatter3d(
                    x=subdata['umap_0'],
                    y=subdata['umap_1'],
                    z=subdata['umap_2'],
                    visible=True,
                    mode='markers', name=f'{value}', text=[],
                    marker=dict(color=color, size=2, symbol='circle',
                                line=dict(width=0, color=color),),
                    legendgroup=f'{value}', showlegend=False)
        trace_tsne = go.Scatter3d(
                    x=subdata['tsne_0'],
                    y=subdata['tsne_1'],
                    z=subdata['tsne_2'],
                    visible=False,
                    mode='markers', name=f'{value}', text=[],
                    marker=dict(color=color, size=2, symbol='circle',
                                line=dict(width=0, color=color),),
                    legendgroup=f'{value}', showlegend=False)

        data.append(trace_umap)
        data.append(trace_tsne)

    umap_inds = list(np.tile([1,0], len(data)//2).astype(bool))
    tsne_inds = umap_inds[::-1]


    updatemenus = [dict(
        buttons=[
            dict(
                args=[dict(visible=umap_inds),
                      umap_scene,
                    ],
                label='UMAP', method='update'),
            dict(
                args=[dict(visible=tsne_inds),
                      tsne_scene,
                    ],
                label='t-SNE', method='update')],
            type='buttons',
            showactive=True
            )]
    default_layout.update(updatemenus=updatemenus)
    print(default_layout.items())

    data = go.Data(data, visible=umap_inds)
    return go.Figure(data=data, layout=default_layout)
