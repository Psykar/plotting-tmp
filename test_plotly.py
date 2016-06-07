import plotly.offline as py
import plotly.graph_objs as go

import numpy as np

# Creating the data
x = np.linspace(-5, 5, 50)
y = np.linspace(-5, 5, 50)
xGrid, yGrid = np.meshgrid(y, x)
R = np.sqrt(xGrid ** 2 + yGrid ** 2)
z = np.sin(R)

# Creating the plot
lines = []
line_marker = dict(color='#0066FF', width=2)
for i, j, k in zip(xGrid, yGrid, z):
    lines.append(go.Scatter3d(x=i, y=j, z=k, mode='lines', line=line_marker))

layout = go.Layout(
    title='Wireframe Plot',
    scene=dict(
        xaxis=dict(
            gridcolor='rgb(255, 255, 255)',
            zerolinecolor='rgb(255, 255, 255)',
            showbackground=True,
            backgroundcolor='rgb(230, 230,230)'
        ),
        yaxis=dict(
            gridcolor='rgb(255, 255, 255)',
            zerolinecolor='rgb(255, 255, 255)',
            showbackground=True,
            backgroundcolor='rgb(230, 230,230)'
        ),
        zaxis=dict(
            gridcolor='rgb(255, 255, 255)',
            zerolinecolor='rgb(255, 255, 255)',
            showbackground=True,
            backgroundcolor='rgb(230, 230,230)'
        )
    ),
    showlegend=False,
)
fig = go.Figure(data=lines, layout=layout)
py.offline.plot(fig, filename='elevations-3d-surface')
