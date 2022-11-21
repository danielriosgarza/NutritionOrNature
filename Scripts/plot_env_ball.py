# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 06:46:39 2022

@author: u0139894
"""

import plotly.io as pio
pio.renderers.default='browser'
import plotly.graph_objects as go

from Env_ball_class import Env_ball

eb=Env_ball(1000)

mc = eb.get_main_components()


def makePlot(mc, figPath):
    x = mc.embedding_.T[0]
    y = mc.embedding_.T[1]
    z = mc.embedding_.T[2]
    
    data = []

    data.append(go.Scatter3d(x=x, y=y, z=z,mode='markers', name='environment compositions', marker=dict(size=3.0, color='#39FF14')))
    
    fig = go.Figure(data=data)

    fig.update_layout(title_text='Environment ball')
    fig.update_layout(legend= {'itemsizing': 'constant'})
    fig.write_image(figPath)
    fig.show()
