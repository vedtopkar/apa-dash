import os
from dash import Dash, html, dcc
from dash.dependencies import Input, Output
import pandas as pd
import numpy as np
import plotly.express as px
import json

text_style = {
    'color': "#506784",
    'font-family': 'Open Sans'
}

app = Dash(__name__)

def description():
    return 'An interactive in-browser track viewer.'

def header_colors():
    return {
        'bg_color': '#0F5BA7',
        'font_color': 'white',
    }

# Differentially expressed genes (identified in R, see assets/data/rna/README.md)

df = pd.read_csv('20220201_counted_pas_for_deseq.csv')
volcano = px.scatter(x=df['log2FoldChange'], y=-np.log(df['padj']), hover_name=df['pas_name'])
volcano.update_layout(title_text='log2FC pAdj')
volcano.update_layout(clickmode='event')

active_pas = 'ENSRNOG00000016516-1'
active_df = df[df['pas_name'] == active_pas]
projection_tpm = float(active_df['Projection_Mean_TPM'])
soma_tpm = float(active_df['Soma_Mean_TPM'])

bars = px.bar(x=[0, 1], y=[soma_tpm, projection_tpm])

app.layout = html.Div([
    html.Div([
        dcc.Graph(
            id='volcano',
            figure=volcano
        )
    ], style={'width': '49%', 'display': 'inline-block', 'padding': '0 20'}),
    html.Div([
        dcc.Graph(
            id='bars',
            figure=bars
        )
    ], style={'display': 'inline-block', 'width': '49%'}),
    html.Div([
        dcc.Markdown("""
            **Click Data**

            Click on points in the graph.
        """),
        html.Pre(id='click-data')])
])



@app.callback(
    Output('click-data', 'children'),
    Input('volcano', 'clickData'))
def display_click_data(clickData):
    return df.iloc[clickData['points'][0]['pointNumber']]['pas_name']

@app.callback(
    Output('bars', 'figure'),
    Input('volcano', 'clickData')
)
def update_expression_bars(clickData):
    active_df = df.iloc[clickData['points'][0]['pointNumber']]
    projection_tpm = float(active_df['Projection_Mean_TPM'])
    soma_tpm = float(active_df['Soma_Mean_TPM'])
    bars = px.bar(x=[0, 1], y=[soma_tpm, projection_tpm])
    return bars


if __name__ == '__main__':
    app.run_server(debug=True)