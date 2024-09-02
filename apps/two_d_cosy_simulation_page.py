import pandas as pd
import numpy as np
import argparse
import json
import pathlib
import base64
import io
import re
import itertools
import plotly.graph_objects as go
import statsmodels.api as sm
from scipy.stats import linregress
import copy
import tracemalloc

import colorsys

import dash
from dash import dcc
from dash import html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State, MATCH, ALL, ClientsideFunction
from dash import dash_table
# import dash_uploader as du

from simulate_2D.construct_hmdb_avg_cons import get_hmdb_normal_avg_cons, get_hmdb_abnormal_avg_cons
from simulate_2D.sample_concentrations import simulate_concentrations, simulate_continuous_concentrations

from simulate_2D.read_parameters import read_param
from simulate_2D.read_2d_spectra import read_2d_cosy
from simulate_2D.match_names import db_match_cons, input_match_db, format_input_mixture, input_cons_match_db
from simulate_2D.match_names import input_corr_match_db, db_names_match_hmdb_names
from simulate_2D.preprocess_2d_spectra import remove_water_calibration, filter_noise, smooth_data, normalize_data
from simulate_2D.plot_2d_spectra import plot_2d_cosy_spectra, plot_cosy_repli_spectra, plot_cosy_repli_spectra_with_ph
from simulate_2D.calculate_2d_without_peak_shift import simulate_cosy_mixture_for_all_repli, \
    simulate_continuous_cosy_mixture_for_all_repli
from simulate_2D.calculate_peak_shift_for_cosy import get_mixture_data_for_all_replicates, \
    conti_get_mixture_data_for_all_replicates


from app import app
tracemalloc.start()

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--parameter", dest="parameter_file", required=True, help="Input Parameter Files")
args = parser.parse_args()
param_dict = read_param(args.parameter_file)
file_path_2d = param_dict['file_path_cosy']

base_path = pathlib.Path(__file__).resolve().parents[1]

# cons_df_1 = pd.read_csv(base_path.joinpath("Input/cons_df_1.csv"), index_col=0)
# cons_df_2 = pd.read_csv(base_path.joinpath("Input/cons_df_2.csv"), index_col=0)
protons_df = pd.read_csv(base_path.joinpath("Input/hmdb_protons.csv"), index_col=0)

data_dict, x_scale, y_scale = read_2d_cosy(file_path_2d)
# print(np.max(data_dict["Citric Acid"]))

with open(base_path.joinpath("Input/hmdb_id_names.json")) as json_file_1:
    hmdb_dict = json.load(json_file_1)

# match_data_dict = db_match_cons(data_dict, cons_df_1, hmdb_dict)
match_data_dict = db_names_match_hmdb_names(data_dict, hmdb_dict)
# print(np.max(match_data_dict["citric acid"]))

with open(base_path.joinpath("Input/hmdb_normal_concentrations.json")) as json_file_2:
    hmdb_norm_cons_dict = json.load(json_file_2)
with open(base_path.joinpath("Input/hmdb_abnormal_concentrations.json")) as json_file_3:
    hmdb_abnorm_cons_dict = json.load(json_file_3)
with open(base_path.joinpath("Input/hmdb_id_pka.json")) as json_file_4:
    hmdb_id_pka_dict = json.load(json_file_4)


page_3_tabs_styles = {
    'height': '50px'
}
page_3_tab_style = {
    'fontWeight': 'bold'
}
page_3_tab_selected_style = {
    'fontWeight': 'bold',
    'color': 'steelblue',
}

#### tab 1
page_3_tab_1 = dcc.Tab(
    id="page-3-tab-1",
    label="STEP 1.   Select Metabolites",
    value="page-3-tab-1",
    style=page_3_tab_style, selected_style=page_3_tab_selected_style,
    children=[
        html.Div(
            [
                html.Br(),
                html.H5("Select Metabolites for Mixture", style={"font-weight": "bold", 'color': 'steelblue'}),
                html.Br(),
                dbc.Card(
                    [dbc.CardBody(
                        [
                            html.Div([
                                html.P("Upload Mixture File", style={"font-weight": "bold"}),
                                dcc.Upload(
                                    id='page-3-upload-mixture',
                                    children=html.Div(
                                        ["Upload A Mixture TXT File."]
                                    ),
                                    style={
                                        'width': '80%',
                                        'height': '50px',
                                        'lineHeight': '40px',
                                        'borderWidth': '1px',
                                        'borderStyle': 'dashed',
                                        'borderRadius': '5px',
                                        'textAlign': 'center',
                                        'margin': '10px',
                                    },
                                    multiple=False,
                                ),
                                dcc.Store(
                                    id='page-3-store-upload-file'
                                ),
                            ], className="input__container",),
                            html.Br(),
                            html.Div([
                                html.P("Or Select Metabolites", style={"font-weight": "bold"}),
                                dcc.Dropdown(
                                    id='page-3-metabolite-select',
                                    options=[{"label": i, "value": i} for i in list(match_data_dict.keys())],
                                    value=['citric acid', 'd-arabitol', 'd-lactate', 'l-alanine', 'l-serine'],
                                    multi=True,),
                            ]),
                            html.Br(),
                        ]
                    )], className="w-75 mb-3",
                ),

                html.Br(),
                html.Div([
                    dbc.Button(id='page-3-submit-meta-btn', children='Submit', n_clicks=0, style={"font-weight": "bold"},
                               size="lg", color='secondary'),
                ]),
                html.Br(),
                # show selected metabolites in the db or not in the db
                html.Br(),
                html.Div(id="page-3-hide-select-metabolites",
                         style={"display": "none"}),
                # show metabolite names that correspond to more than one HMDB IDs
                html.Br(),
                html.Div(id="page-3-hide-multi-hmdb-names",
                         style={"display": "none"}),
                # show and confirm the selections for multiple names
                html.Br(),
                html.Div(id="page-3-hide-confirm-select",
                         style={"display": "none"}),
                # show the name-hmdb_id-hmdb_name table
                html.Br(),
                html.Div(id="page-3-hide-mixture-table",
                         style={"display": "none"}),
                # show the conform button
                html.Br(),
                html.Div(
                    id="page-3-hide-confirm-btn",
                    children=[
                        dbc.Button(id='page-3-confirm-meta-btn', children='Confirm', n_clicks=0,
                                   style={"font-weight": "bold"}, size='lg', color='primary')
                    ],
                    style={"display": "none"},
                ),
                html.Br(),
                # store metabolite names and corresponding HMDB IDs
                dcc.Store(id='page-3-db-names-hmdb-ids-dict'),
                html.Br(),
            ],
            className="container__1",
        )
    ]
)

# -------------------------------------- tab 2 -----------------------------------------
page_3_select_upload_cons_div = html.Div(
    id="page-3-select-upload-cons-div",
    children=[
        dbc.Row(
            [
                dbc.Col(html.H6("Please decide the source of concentration data:",
                                style={"font-weight": "bold", 'font-style': 'italic', 'color': 'darkred'})),
                dbc.Col(dbc.RadioItems(
                    id="page-3-select-upload-cons",
                    options=[
                        {"label": "Use the concentration simulator", "value": "not upload"},
                        {"label": "Upload your own concentration data", "value": "upload"}
                    ],
                    inline=True,
                    labelCheckedClassName="font-weight-bold text-danger",
                    inputCheckedClassName="border border-primary bg-danger",),
                )
            ],
            justify="start",
            className="g-0",
        ),
        html.Br(),
    ],
    style={"display": "none"},
)

page_3_template_df_1 = pd.DataFrame(
    {
        ("Group 1", "Replicate_1"): {
            "name 1": "...",
            "name 2": "..."
        },
        ("Group 1", "..."): {
            "name 1": "...",
            "name 2": "..."
        },
        ("Group 1", "Replicate_n"): {
            "name 1": "...",
            "name 2": "..."
        },
        ("Group 2", "Replicate_1"): {
            "name 1": "...",
            "name 2": "..."
        },
        ("Group 2", "..."): {
            "name 1": "...",
            "name 2": "..."
        },
        ("Group 2", "Replicate_n"): {
            "name 1": "...",
            "name 2": "..."
        },
    }
)
page_3_template_df_1.index.set_names("Name", inplace=True)

page_3_template_table_1 = dbc.Table.from_dataframe(
    page_3_template_df_1, striped=True, bordered=True, hover=True, index=True
)

page_3_upload_discrete_cons_div = html.Div(
    id="page-3-upload-discrete-cons-data-div",
    children=[
        dbc.Card(
            dbc.CardBody(
                [
                    html.P("Please upload a .csv file following the structures:", style={"font-weight": "bold"}),
                    html.Li("The index is the metabolite names;"),
                    html.Li("The column is the replicates in each group"),
                    page_3_template_table_1,
                    html.Br(),
                    dcc.Upload(
                        id='page-3-upload-final-discrete-cons-file',
                        children=dbc.Button("Please Upload A .CSV File", color="primary", outline=True, size="lg"),
                        multiple=False,
                    ),
                    html.Br(),
                ]
            )
        ),
        html.Br(),
        dcc.Loading(html.Div(
            id="page-3-upload-discrete-cons-table-div",
            children=[
                dash_table.DataTable(id="page-3-upload-discrete-cons-table",
                                     merge_duplicate_headers=True, editable=True, virtualization=True,
                                     style_cell={
                                         'font_family': 'arial',
                                         'font_size': '16px',
                                         'text_align': 'center'
                                     },
                                     style_data_conditional=[
                                         {
                                             "if": {"column_id": "meta_name"},
                                             "fontWeight": "bold",
                                         },
                                     ],
                                     fixed_rows={'headers': True},
                                     style_table={'height': 250},
                                     css=[
                                         {"selector": ".dash-spreadsheet tr th", "rule": "height: 35px;"},
                                         {"selector": ".dash-spreadsheet tr td", "rule": "height: 35px;"},
                                     ]
                                    ),
            ]
        ),),
        html.Br(),
        html.Div(
            id="page-3-confirm-upload-discrete-cons-btn-div",
            children=[
                dbc.Button(id='page-3-confirm-upload-discrete-cons-btn', children='Confirm', n_clicks=0, style={"font-weight": "bold"},
                           size='lg', color='primary'),
            ],
            style={"display": "none"}
        ),
    ],
    style={"display": "none"},
)

page_3_template_df_2 = pd.DataFrame(
    {"Name": ["Y", "name_1", "name_2"],
     "Replicate_1": ["...", "...", "..."],
     "...": ["...", "...", "..."],
     "Replicate_n": ["...", "...", "..."],}
)

page_3_template_table_2 = dbc.Table.from_dataframe(page_3_template_df_2, striped=True, bordered=True, hover=True)

page_3_upload_conti_cons_div = html.Div(
    id="page-3-upload-conti-cons-data-div",
    children=[
        dbc.Card(
            dbc.CardBody(
                [
                    html.P("Please upload a .csv file following the structures:", style={"font-weight": "bold"}),
                    html.Li("The index is the metabolite names;"),
                    html.Li("The column is all the replicates"),
                    page_3_template_table_2,
                    html.Br(),
                    dcc.Upload(
                        id='page-3-upload-final-conti-cons-file',
                        children=dbc.Button("Please Upload A .CSV File", color="primary", outline=True, size="lg"),
                        multiple=False,
                    ),
                    html.Br(),
                ]
            )
        ),
        html.Br(),
        html.Div(
            id="page-3-upload-conti-cons-table-div",
            children=[
                dash_table.DataTable(id="page-3-upload-conti-cons-table",
                                     editable=True, virtualization=True,
                                     style_cell={
                                         'font_family': 'arial',
                                         'font_size': '16px',
                                         'text_align': 'center'
                                     },
                                     style_data_conditional=[
                                         {
                                             "if": {"column_id": "meta_name"},
                                             "fontWeight": "bold",
                                         },
                                     ],
                                     fixed_rows={'headers': True},
                                     style_table={'height': 250},
                                     css=[
                                         {"selector": ".dash-spreadsheet tr th", "rule": "height: 35px;"},
                                         {"selector": ".dash-spreadsheet tr td", "rule": "height: 35px;"},
                                     ]
                ),
            ]
        ),
        html.Br(),
        html.Div(
            id="page-3-confirm-upload-conti-cons-btn-div",
            children=[
                dbc.Button(id='page-3-confirm-upload-conti-cons-btn', children='Confirm', n_clicks=0, style={"font-weight": "bold"},
                           size='lg', color='primary'),
            ],
            style={"display": "none"}
        ),
    ],
    style={"display": "none"},
)

page_3_upload_final_cons_div = html.Div(
    id="page-3-upload-final-cons-div",
    children=[
        dbc.Row(
            [dbc.Col(html.P("Please select the type of concentration data:",
                            style={"font-weight": "bold", 'font-style': 'italic', 'color': 'green'}), ),
             dbc.Col(dbc.RadioItems(
                 id="page-3-select-upload-cons-type",
                 options=[
                     {"label": "Discrete outcomes (two groups)", "value": "discrete"},
                     {"label": "Continuous outcomes (one group)", "value": "continuous"}
                 ],
                 # value="group",
                 inline=True,
                 labelCheckedClassName="font-weight-bold text-success",
                 inputCheckedClassName="border border-primary bg-success", ),
                 # width=4,
             )],
            justify="start",
            className="g-0",
        ),
        html.Br(),
        page_3_upload_discrete_cons_div,
        page_3_upload_conti_cons_div,
    ],
    style={"display": "none"},
)

page_3_select_type_cons_div = html.Div(
    id="page-3-select-cons-type-div",
    children=[
        dbc.Row(
            [dbc.Col(html.P("Please select the type of spectra simulation:",
                     style={"font-weight": "bold", 'font-style': 'italic', 'color': 'darkred'}), width=3),
            dbc.Col(dbc.RadioItems(
                id="page-3-select-simulation-type",
                options=[
                    {"label": "Simulate with two groups", "value": "group"},
                    {"label": "Simulate continuous outcomes", "value": "continuous"}
                ],
                inline=True,
                labelCheckedClassName="font-weight-bold text-primary",
                inputCheckedClassName="border border-primary bg-primary",),
                width=4,
            )],
            justify="start",
            className="g-0",
        ),
        html.Br(),
    ],
    style={"display": "none"},
)

# ------ tab 2: group results div --------
page_3_part_1_group_1_card = dbc.Card(
    id="page-3-part-1-group-1-card",
    children=[
        dbc.CardHeader(html.H6("Group 1 (Normal Group)",
                               style={"font-weight": "bold", "color": "steelblue"})),
        dbc.CardBody(
            [
                html.Div(
                    [
                        dbc.RadioItems(
                            id="page-3-part-1-group-1-hmdb-data",
                            options=[
                                {"label": "Use concentration data in HMDB", "value": "hmdb"},
                                {"label": "Upload your own concentration data", "value": "not hmdb"},
                            ],
                            value="hmdb",
                            inline=True,
                        ),
                    ]
                ),
                html.Br(),
                html.Div(
                    id="page-3-group-1-hmdb-cons-data-select",
                    children=[
                        html.P("Please select the type of biospecimen", style={"font-weight": "bold"}),
                        dbc.RadioItems(
                            id="page-3-group-1-bio-type",
                            options=[
                                {"label": "Blood", "value": "Blood"},
                                {"label": "Urine", "value": "Urine"},
                                {"label": "Cerebrospinal Fluid (CSF)", "value": "Cerebrospinal Fluid (CSF)"}
                            ],
                            value="Blood",
                            inline=True,
                        )
                    ],
                    style={"display": "none"},
                ),
                html.Div(
                    id="page-3-group-1-not-hmdb-cons-data-upload",
                    children=[
                        html.P("Please upload a concentration file containing mean & std", style={"font-weight": "bold"}),
                        dcc.Upload(
                            id='page-3-group-1-upload-mean-std',
                            children=dbc.Button("upload a concentration file", color="primary", outline=True),
                            multiple=False,
                        ),
                    ],
                    style={"display": "none"},
                ),
                html.Br(),
            ]
        )
    ]
)

page_3_part_1_group_2_card = dbc.Card(
    id="page-3-part-1-group-2-card",
    children=[
        dbc.CardHeader(html.H6("Group 2 (Abnormal Group)",
                               style={"font-weight": "bold", "color": "steelblue"})),
        dbc.CardBody(
            [
                html.Div(
                    [
                        dbc.RadioItems(
                            id="page-3-part-1-group-2-hmdb-data",
                            options=[
                                {"label": "Use concentration data in HMDB", "value": "hmdb"},
                                {"label": "Upload your own concentration data", "value": "not hmdb"},
                            ],
                            value="hmdb",
                            inline=True,
                        ),
                    ]
                ),
                html.Br(),
                html.Div(
                    id="page-3-group-2-hmdb-cons-data-select",
                    children=[
                        html.P("Please select the type of biospecimen", style={"font-weight": "bold"}),
                        dbc.RadioItems(
                            id="page-3-group-2-bio-type",
                            options=[
                                {"label": "Blood (Heart Transplant)", "value": "Blood"},
                                {"label": "Urine", "value": "Urine"},
                                {"label": "Cerebrospinal Fluid (CSF)", "value": "Cerebrospinal Fluid (CSF)"}
                            ],
                            value="Blood",
                            inline=True,
                        )
                    ],
                    style={"display": "none"},
                ),
                html.Div(
                    id="page-3-group-2-not-hmdb-cons-data-upload",
                    children=[
                        html.P("Please upload a concentration file containing mean & std", style={"font-weight": "bold"}),
                        dcc.Upload(
                            id='page-3-group-2-upload-mean-std',
                            children=dbc.Button("upload a concentration file", color="primary", outline=True),
                            multiple=False,
                        ),
                    ],
                    style={"display": "none"},
                ),
                html.Br(),
            ]
        )
    ]
)

page_3_part_1_div = html.Div(
    id="page-3-group-part-1-div",
    children=[
        html.H6("Part I: determine the mean & std of concentrations", style={"font-weight": "bold"}),
        dbc.Row(
            [
                dbc.Col(page_3_part_1_group_1_card),
                dbc.Col(page_3_part_1_group_2_card)
            ]
        ),
        html.Br(),
        html.Div([
                    dbc.Button(id='page-3-submit-part-1-btn', children='OK', n_clicks=0, style={"font-weight": "bold"},
                               size='lg', color='secondary'),
        ]),
        html.Br(),
        html.Br(),
    ],
    style={"display": "none"},
)

page_3_hide_cons_mean_std_table_div = html.Div(
    id="page-3-hide-cons-mean-std-table-div",
    children=[
        dash_table.DataTable(id="page-3-cons-mean-std-table", export_format="csv", export_headers="display",
                             merge_duplicate_headers=True, editable=True, virtualization=True,
                             style_cell={
                               'font_family': 'arial',
                               'font_size': '16px',
                               'text_align': 'center'
                             },
                             style_data_conditional=[
                               {
                                   "if": {"column_id": "meta_name"},
                                   "fontWeight": "bold",
                               },
                             ],
                             fixed_rows={'headers': True},
                             style_table={'height': 250},
                             css=[
                               {"selector": ".dash-spreadsheet tr th", "rule": "height: 35px;"},
                               {"selector": ".dash-spreadsheet tr td", "rule": "height: 35px;"},
                               {"selector": ".export", "rule": "position: absolute; right:0%; top:-40px; "
                                                               "background-color: lightsteelblue; "
                                                               # "border-color: rgb(107,174,214); "
                                                               "font-weight: 300; font-size: 1.5rem;"},
                             ]),
        html.Br(),
        html.Br(),
    ],
    style={"display": "none"}
)

page_3_part_2_group_1_card = dbc.Card(
    id="page-3-part-2-group-1-card",
    children=[
        dbc.CardHeader(html.H6("Group 1 (Normal Group)", style={"font-weight": "bold", "color": "steelblue"})),
        dbc.CardBody(
            [
                html.Div(
                    [
                        html.P("Number of Replicates", style={"font-weight": "bold"}),
                        dcc.Input(
                            id='page-3-group-1-num-repli',
                            type='number',
                            debounce=True,
                            value=5,
                            min=1,
                            step=1,
                            size="lg",
                        ),
                    ],
                ),
                html.Br(),
                html.Div(
                    [
                        html.P("Inter-Metabolites Correlation", style={"font-weight": "bold"}),
                        dbc.RadioItems(
                            id='page-3-group-1-corr-flag',
                            options=[
                                {'label': 'False', 'value': False},
                                {'label': 'True', 'value': True},
                            ],
                            value=False,
                            inline=True,
                        ),
                    ],
                ),
                html.Br(),
                html.Div(
                    id="page-3-group-1-upload-corr-file-div",
                    children=[
                        dcc.Upload(
                            id='page-3-group-1-corr-file',
                            children=dbc.Button("If True, select a txt file to upload.", color="primary", outline=True),
                            multiple=False,
                        ),
                    ],
                    style={"display": "none"},
                ),
                html.Br(),
            ]
        )
    ]
)

page_3_part_2_group_2_card = dbc.Card(
    id="page-3-part-2-group-2-card",
    children=[
        dbc.CardHeader(html.H6("Group 2 (Abnormal Group)", style={"font-weight": "bold", "color": "steelblue"})),
        dbc.CardBody(
            [
                html.Div(
                    [
                        html.P("Number of Replicates", style={"font-weight": "bold"}),
                        dcc.Input(
                            id='page-3-group-2-num-repli',
                            type='number',
                            debounce=True,
                            value=5,
                            min=1,
                            step=1,
                            size="lg",
                        ),
                    ],
                ),
                html.Br(),
                html.Div(
                    [
                        html.P("Inter-Metabolites Correlation", style={"font-weight": "bold"}),
                        dbc.RadioItems(
                            id='page-3-group-2-corr-flag',
                            options=[
                                {'label': 'False', 'value': False},
                                {'label': 'True', 'value': True},
                            ],
                            value=False,
                            inline=True,
                        ),
                    ],
                ),
                html.Br(),
                html.Div(
                    id="page-3-group-2-upload-corr-file-div",
                    children=[
                        dcc.Upload(
                            id='page-3-group-2-corr-file',
                            children=dbc.Button("If True, select a txt file to upload.", color="primary", outline=True),
                            multiple=False,
                        ),
                    ],
                    style={"display": "none"},
                ),
                html.Br(),
            ]
        )
    ]
)

page_3_part_2_div = html.Div(
    id="page-3-group-part-2-div",
    children=[
        html.H6("Part II: set parameters for simulating concentrations", style={"font-weight": "bold"}),
        dbc.Row(
            [
                dbc.Col(page_3_part_2_group_1_card),
                dbc.Col(page_3_part_2_group_2_card)
            ]
        ),
        html.Br(),
        html.Div([
                    dbc.Button(id='page-3-submit-part-2-btn', children='OK', n_clicks=0, style={"font-weight": "bold"},
                               size='lg', color='secondary'),
        ]),
        html.Br(),
        html.Br()
    ],
    style={"display": "none"}
)

page_3_hide_simulated_cons_table_div = html.Div(
    id="page-3-hide-simulated-cons-table-div",
    children=[
        dash_table.DataTable(id="page-3-simulated-cons-table", export_format="csv", export_headers="display", editable=True,
                             merge_duplicate_headers=True, virtualization=True,
                             style_cell={
                               'font_family': 'arial',
                               'font_size': '16px',
                               'text_align': 'center'
                             },
                             style_data_conditional=[
                               {
                                   "if": {"column_id": "meta_name"},
                                   "fontWeight": "bold",
                               },
                             ],
                             fixed_rows={'headers': True},
                             style_table={'height': 250},
                             css=[
                               {"selector": ".dash-spreadsheet tr th", "rule": "height: 35px;"},
                               {"selector": ".dash-spreadsheet tr td", "rule": "height: 35px;"},
                               {"selector": ".export", "rule": "position: absolute; right:0%; top:-1px; "
                                                               "background-color: lightsteelblue; "
                                                               "font-weight: 300; font-size: 1.5rem;"},
                             ]),
        html.Br(),
        html.Br(),
    ],
    style={"display": "none"}
)

page_3_hide_box_plot_div = html.Div(
    id="page-3-hide-box-plot-div",
    children=[
        dcc.Graph(id="page-3-groups-cons-box-plot"),
        html.Br(),
        html.Br(),
    ],
    style={"display": "none"}
)

page_3_group_type_div = html.Div(
    id="page-3-group-sim-div",
    children=[
        page_3_part_1_div,
        page_3_hide_cons_mean_std_table_div,
        page_3_part_2_div,
        page_3_hide_simulated_cons_table_div,
        page_3_hide_box_plot_div,
        html.Div(
            id="page-3-confirm-group-cons-btn-div",
            children=[
                dbc.Button(id='page-3-confirm-group-cons-btn', children='Confirm', n_clicks=0, style={"font-weight": "bold"},
                           size='lg', color='primary'),
            ],
            style={"display": "none"}
        ),
    ],
    style={"display": "none"},
)

# ------ tab 2: continuous results div --------
page_3_conti_select_hmdb_cons_data_div = html.Div(
    id="page-3-conti-hmdb-cons-data-select-div",
    children=[
        html.P("Please select the type of biospecimen", style={"font-style": 'italic'}),
        dbc.RadioItems(
            id="page-3-conti-bio-type",
            options=[
                {"label": "Blood", "value": "Blood"},
                {"label": "Urine", "value": "Urine"},
                {"label": "Cerebrospinal Fluid (CSF)", "value": "Cerebrospinal Fluid (CSF)"}
            ],
            value="Blood",
            inline=True,
        ),
        html.Br(),
    ],
    style={"display": "none"},
)

page_3_conti_upload_cons_data_div = html.Div(
    id="page-3-conti-upload-cons-data-div",
    children=[
        html.P("Please upload a concentration file containing 5 columns (name, mean, std, a, b)",
               style={"font-style": 'italic'}),
        dcc.Upload(
            id='page-3-conti-upload-mean-std-file',
            children=dbc.Button("upload a concentration file (.csv / .txt)", color="primary", outline=True),
            multiple=False,
        ),
        html.Br(),
    ],
    style={"display": "none"}
)

page_3_conti_part_1_div = html.Div(
    id="page-3-conti-part-1-div",
    children=[
        html.P("The linear relationship between y and x can be written as: y = ax + b; where x: concentration, "
               "y: the outcome (e.g., age or BMI)", style={"font-weight": "bold", "color": "olivedrab"}),
        # html.Br(),
        html.H6("Part I: determine the source of concentration data", style={"font-weight": "bold"}),
        dbc.Card(
            [dbc.CardBody(
                [
                    dbc.RadioItems(
                        id="page-3-conti-part-1-hmdb-data",
                        options=[
                            {"label": "Use concentration data in HMDB", "value": "hmdb"},
                            {"label": "Upload your own concentration data", "value": "not hmdb"},
                        ],
                        value="hmdb",
                        inline=True,
                    ),
                    html.Br(),
                    page_3_conti_select_hmdb_cons_data_div,
                    page_3_conti_upload_cons_data_div,
                    # html.Br(),
                ]
            )], className="w-75 mb-3",
        ),
        html.Br(),
        html.Div([
            dbc.Button(id='page-3-conti-submit-part-1-btn', children='OK', n_clicks=0, style={"font-weight": "bold"},
                       size='lg', color='secondary'),
        ]),
        html.Br(),
        html.Br(),
    ],
    style={"display": "none"}
)

page_3_conti_hide_cons_mean_std_table_div = html.Div(
    id="page-3-conti-hide-cons-mean-std-table-div",
    children=[
        dash_table.DataTable(id="page-3-conti-mean-std-ab-table", export_format="csv", export_headers="display",
                             editable=True, virtualization=True,
                             style_cell={
                               'font_family': 'arial',
                               'font_size': '16px',
                               'text_align': 'center',
                               'minWidth': '180px', 'width': '180px', 'maxWidth': '180px',
                             },
                             style_header={
                                'backgroundColor': 'whitesmoke',
                                'fontWeight': 'bold'
                             },
                             style_data_conditional=[
                               {
                                   "if": {"column_id": "meta_name"},
                                   "fontWeight": "bold",
                               },
                             ],
                             fixed_rows={'headers': True},
                             style_table={'height': 250},
                             css=[
                               {"selector": ".dash-spreadsheet tr th", "rule": "height: 35px;"},
                               {"selector": ".dash-spreadsheet tr td", "rule": "height: 35px;"},
                               {"selector": ".export", "rule": "position: absolute; right:0%; top:-40px; "
                                                               "background-color: lightsteelblue; "
                                                               # "border-color: rgb(107,174,214); "
                                                               "font-weight: 300; font-size: 1.5rem;"},
                             ]),
        html.Br(),
        # html.Br(),
    ],
    style={"display": "none"}
)

page_3_conti_normal_distribution_div = html.Div(
    id="page-3-conti-normal-div",
    children=[
        html.P("Please specify the mean and std of y", style={"font-style": 'italic'}),
        dbc.Row(
            [
                dbc.Col(dbc.InputGroup([
                    dbc.InputGroupText("mean=", style={"font-weight": "bold"}),
                    dbc.Input(
                        id="page-3-y-mean",
                        type="number",
                        value=50,
                    )
                ], size="lg")),
                dbc.Col(dbc.InputGroup([
                    dbc.InputGroupText("std=", style={"font-weight": "bold"}),
                    dbc.Input(
                        id="page-3-y-std",
                        type="number",
                        value=20
                    )
                ], size="lg")),
            ]
        ),
    ],
    style={"display": "none"}
)

page_3_conti_uniform_distribution_div = html.Div(
    id="page-3-conti-uniform-div",
    children=[
        html.P("Please specify the min and max of y", style={"font-style": 'italic'}),
        dbc.Row(
            [
                dbc.Col(dbc.InputGroup([
                    dbc.InputGroupText("min=", style={"font-weight": "bold"}),
                    dbc.Input(
                        id="page-3-y-min",
                        type="number",
                        value=18,
                    )
                ], size="lg")),
                dbc.Col(dbc.InputGroup([
                    dbc.InputGroupText("max=", style={"font-weight": "bold"}),
                    dbc.Input(
                        id="page-3-y-max",
                        type="number",
                        value=100
                    )
                ], size="lg")),
            ]
        ),
    ],
    style={"display": "none"}
)

page_3_conti_part_2_div = html.Div(
    id="page-3-conti-part-2-div",
    children=[
        html.H6("Part II: set the parameters for simulating concentrations", style={"font-weight": "bold"}),
        dbc.Card(
            [dbc.CardBody(
                [
                    html.P("Define the number of replicates of y: ", style={"font-weight": "bold"}),
                    dcc.Input(
                        id='page-3-y-num-repli',
                        type='number',
                        debounce=True,
                        value=5,
                        min=1,
                        step=1,
                        size="lg",
                    ),
                    html.Br(),
                    html.Br(),
                    html.P("Select the distribution of y: ", style={"font-weight": "bold"}),
                    dbc.RadioItems(
                        id="page-3-conti-part-2-select-distribution",
                        options=[
                            {"label": "Normal distribution", "value": "normal"},
                            {"label": "Uniform distribution", "value": "uniform"},
                        ],
                        value="normal",
                        inline=True,
                    ),
                    html.Br(),
                    page_3_conti_normal_distribution_div,
                    page_3_conti_uniform_distribution_div,
                    html.Br(),
                ]
            )], className="w-75 mb-3",
        ),
        html.Br(),
        html.Div([
            dbc.Button(id='page-3-conti-submit-part-2-btn', children='OK', n_clicks=0, style={"font-weight": "bold"},
                       size='lg', color='secondary'),
        ]),
        html.Br(),
        html.Br(),
    ],
    style={"display": "none"}
)

page_3_conti_hide_simulated_cons_table_div = html.Div(
    id="page-3-conti-hide-simulated-cons-table-div",
    children=[
        dash_table.DataTable(id="page-3-conti-simulated-cons-table", export_format="csv", export_headers="display",
                             editable=True, virtualization=True,
                             style_cell={
                               'font_family': 'arial',
                               'font_size': '16px',
                               'text_align': 'center',
                               'minWidth': '180px', 'width': '180px', 'maxWidth': '180px',
                             },
                             style_header={
                                'backgroundColor': 'whitesmoke',
                                'fontWeight': 'bold'
                             },
                             style_data_conditional=[
                               {
                                   "if": {"column_id": "meta_name"},
                                   "fontWeight": "bold",
                               },
                             ],
                             fixed_rows={'headers': True},
                             style_table={'height': 250},
                             css=[
                               {"selector": ".dash-spreadsheet tr th", "rule": "height: 35px;"},
                               {"selector": ".dash-spreadsheet tr td", "rule": "height: 35px;"},
                               {"selector": ".export", "rule": "position: absolute; right:0%; top:-40px; "
                                                               "background-color: lightsteelblue; "
                                                               # "border-color: rgb(107,174,214); "
                                                               "font-weight: 300; font-size: 1.5rem;"},
                             ]),
        html.Br(),
    ],
    style={"display": "none"}
)

page_3_conti_hide_box_plot_div = html.Div(
    id="page-3-conti-hide-plot-div",
    children=[
        dbc.Row(
            [
                dbc.Col(dcc.Graph(id="page-3-conti-line-plot")),
                dbc.Col(dcc.Graph(id="page-3-conti-box-plot")),
            ]
        ),
        html.Br(),
        html.P("Table of R-squared between Y and metabolites:", style={"font-weight": "bold"}),
        dbc.Table(id="page-3-conti-r-squared-table", bordered=True, color="light"),
        html.Br(),
    ],
    style={"display": "none"}
)

page_3_conti_type_div = html.Div(
    id="page-3-conti-sim-div",
    children=[
        page_3_conti_part_1_div,
        page_3_conti_hide_cons_mean_std_table_div,
        page_3_conti_part_2_div,
        page_3_conti_hide_simulated_cons_table_div,
        page_3_conti_hide_box_plot_div,
        html.Div(
            id="page-3-confirm-conti-cons-btn-div",
            children=[
                dbc.Button(id='page-3-confirm-conti-cons-btn', children='Confirm', n_clicks=0, style={"font-weight": "bold"},
                           size='lg', color='primary'),
            ],
            style={"display": "none"}),
        html.Br(),
        html.Br(),
    ],
    style={"display": "none"},
)

# tab 2
page_3_tab_2 = dcc.Tab(
    id="page-3-tab-2",
    label="STEP 2.   Simulate Concentrations",
    value="page-3-tab-2",
    style=page_3_tab_style, selected_style=page_3_tab_selected_style,
    children=[
        html.Div(
            [
                html.Br(),
                html.H5("Simulate Concentrations", style={"font-weight": "bold", 'color': 'steelblue'}),
                html.Br(),
                page_3_select_upload_cons_div,
                page_3_upload_final_cons_div,
                page_3_select_type_cons_div,
                page_3_group_type_div,
                page_3_conti_type_div,
                html.Br(),
            ],
            className="container__1",
        )
    ],
)

# -------------------------------------- tab 3 -----------------------------------------
page_3_select_peak_shift_div = html.Div(
    id="page-3-select-peak-shift-div",
    children=[
        dbc.Row(
            [dbc.Col(html.P("Please select the type of spectra simulation:",
                     style={"font-weight": "bold", 'font-style': 'italic', 'color': 'darkred'}), width=3),
            dbc.Col(dbc.RadioItems(
                id="page-3-select-spectra-simulation-type",
                options=[
                    {"label": "Simulate without peak shift", "value": "no_peak_shift"},
                    {"label": "Simulate with peak shift", "value": "peak_shift"}
                ],
                inline=True,
                labelCheckedClassName="font-weight-bold text-primary",
                inputCheckedClassName="border border-primary bg-primary",),
                width=4,
            )],
            justify="start",
            className="g-0",
        ),
        html.Br(),
    ],
    style={"display": "none"},
)

# ---------- tab 3: group no peak shift results div ----------
page_3_advance_parameters_setting_card = dbc.Card(
    id="page-3-ad-param-setting-card",
    children=[dbc.CardBody([
        html.H6("Advanced Parameters Setting", style={"font-weight": "bold", "color": "steelblue"}),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(dbc.Card(
                    [
                        dbc.CardHeader(html.P("Remove exclusion regions", style={"font-weight": "bold"})),
                        dbc.CardBody(
                            html.Div(
                                [
                                    html.P("Range of calibration peaks (ppm)", style={"font-style": 'italic'}),
                                    dcc.RangeSlider(
                                        id="page-3-calibration-range",
                                        min=-0.3, max=0.3, step=0.1,
                                        value=[-0.3, 0.3],
                                    ),
                                    html.Br(),
                                    html.P("Range of water suppression (ppm)", style={"font-style": 'italic'}),
                                    dcc.RangeSlider(
                                        id="page-3-water-range",
                                        min=4.5, max=5.5, step=0.1,
                                        value=[4.5, 5.0]
                                    )
                                ]
                            ),
                        )
                    ], style={'height': '25vh'}, color="primary", outline=True
                )),
                dbc.Col(dbc.Card(
                    [
                        dbc.CardHeader(html.P("Filter noise", style={"font-weight": "bold"})),
                        dbc.CardBody(
                            html.Div(
                                [
                                    html.P("Parameters for filtering noise", style={"font-style": 'italic'}),
                                    dbc.InputGroup([
                                        dbc.InputGroupText("threshold = ", style={"font-weight": "bold"}),
                                        dbc.Select(
                                            id="page-3-filter-thres",
                                            options=[
                                                {"label": "0.01", "value": 0.01},
                                                {"label": "0.05", "value": 0.05},
                                            ],
                                            value=0.05,
                                        )
                                    ], size="lg"),
                                ]
                            ),
                        )
                    ], style={'height': '25vh'}, color="success", outline=True
                )),
                dbc.Col(dbc.Card(
                    [
                        dbc.CardHeader(html.P("Smooth spectrum", style={"font-weight": "bold"})),
                        dbc.CardBody(
                            html.Div(
                                [
                                    html.P("Thresholds for smoothing", style={"font-style": 'italic'}),
                                    dbc.InputGroup([
                                        dbc.InputGroupText("m = ", style={"font-weight": "bold"}),
                                        dbc.Input(
                                            id="page-3-smooth-thres-m",
                                            type="number",
                                            value=1,
                                        )
                                    ], size="lg"),
                                    html.Br(),
                                    dbc.InputGroup([
                                        dbc.InputGroupText("n = ", style={"font-weight": "bold"}),
                                        dbc.Input(
                                            id="page-3-smooth-thres-n",
                                            type="number",
                                            value=1,
                                        )
                                    ], size="lg"),
                                ]
                            ),
                        )
                    ], style={'height': '25vh'}, color="danger", outline=True
                )),
                dbc.Col(dbc.Card(
                    [
                        dbc.CardHeader(html.P("Add noise", style={"font-weight": "bold"})),
                        dbc.CardBody(
                            html.Div(
                                [
                                    html.P("SNR for adding noise", style={"font-style": 'italic'}),
                                    dbc.RadioItems(
                                        id="page-3-snr-noise",
                                        options=[
                                            {"label": "1000", "value": 1000},
                                            {"label": "5000", "value": 5000},
                                            {"label": "10000", "value": 10000},
                                        ],
                                        value=1000,
                                        inline=True,
                                    ),
                                ]
                            ),
                        )
                    ], style={'height': '25vh'}, color="warning", outline=True
                )),
            ]
        ),
        html.Br(),
        html.Div(
            [
                dbc.Button(id='page-3-confirm-param-btn', children='Simulate spectra', n_clicks=0,
                           style={"font-weight": "bold"}, size='lg', color='primary')
            ],
            className="d-grid gap-2 d-md-flex justify-content-md-center",
        ),
        dcc.Store(id='page-3-processed-spectra-data-dict'),
        html.Br(),
    ])]
)

page_3_all_metabolites_spectra_fig_div = dbc.Card(
    id='page-3-all-spectra-div',
    children=[
        dbc.CardBody(
            [
                html.H6("Metabolite spectrum in the selected mixture", style={"font-weight": "bold", "color": "steelblue"}),
                html.Br(),
                dcc.Loading(children=[
                    dcc.Dropdown(id='page-3-select-mix', multi=False,),
                    dcc.Graph(id='page-3-all-meta-fig')
                ], type="circle"),
            ]
        )
    ],
)

page_3_group_1_results_div = dbc.Card(
    id="page-3-group-1-spectra-div",
    children=[
        dbc.CardBody(
            [
                html.H6("Results of Normal Group (Group 1)", style={"font-weight": "bold", "color": "steelblue"}),
                dcc.Loading(children=[
                    dcc.Dropdown(id='page-3-select-repli-1', multi=False,),
                    dcc.Graph(id="page-3-group-1-spectra-fig"),
                    dcc.Store(id='page-3-group-1-final-dict'),
                    # html.P("Contour Level Slider:"),
                    # dcc.RangeSlider(id='page-3-contour-level-1', allowCross=False,
                    #                 min=0, max=100),
                    # html.Div(id='page-3-update-level-container-1', style={'margin-top': 10}),
                    html.Br(),
                    html.Button("Download simulated spectra as csv", id="page-3-group-1-btn-csv"),
                    dcc.Download(id="page-3-group-1-download-spectra-csv"),
                ]),
            ]
        )
    ]
)

page_3_group_2_results_div = dbc.Card(
    id="page-3-group-2-spectra-div",
    children=[
        dbc.CardBody(
            [
                html.H6("Results of Abnormal Group (Group 2)", style={"font-weight": "bold", "color": "steelblue"}),
                # html.Br(),
                dcc.Loading(children=[
                    dcc.Dropdown(id='page-3-select-repli-2', multi=False,),
                    dcc.Graph(id="page-3-group-2-spectra-fig"),
                    dcc.Store(id='page-3-group-2-final-dict'),
                    # html.P("Contour Level Slider:"),
                    # dcc.RangeSlider(id='page-3-contour-level-2', allowCross=False,
                    #                 min=0, max=100),
                    # html.Div(id='page-3-update-level-container-2', style={'margin-top': 10}),
                    html.Br(),
                    html.Button("Download simulated spectra as csv", id="page-3-group-2-btn-csv"),
                    dcc.Download(id="page-3-group-2-download-spectra-csv"),
                ]),
            ]
        )
    ]
)

page_3_group_no_peak_shift_results_div = html.Div(
    id="page-3-group-no-peak-shift-results-div",
    children=[
        page_3_advance_parameters_setting_card,
        html.Br(),
        page_3_all_metabolites_spectra_fig_div,
        html.Br(),
        dbc.Row(
            [
                dbc.Col(page_3_group_1_results_div),
                dbc.Col(page_3_group_2_results_div),
            ]
        ),
    ],
    style={"display": "none"},
)

# ---------- tab 3: group with peak shift results div ----------
page_3_group_1_pH_same_div = html.Div(
    id="page-3-group-1-ph-same-div",
    children=[
        html.P("Please input the pH:", style={"font-style": 'italic'}),
        dbc.InputGroup([
            dbc.InputGroupText("pH = ", style={"font-weight": "bold"}),
            dbc.Input(
                id="page-3-group-1-same-ph",
                type="number",
                value=6.5,
            )
        ], size="lg")
    ],
    style={"display": "none"}
)

page_3_group_1_pH_not_same_div = html.Div(
    id="page-3-group-1-ph-not-same-div",
    children=[
        html.P("Please input the mean and std for pH:", style={"font-style": 'italic'}),
        dbc.Row(
            [
                dbc.Col(dbc.InputGroup([
                    dbc.InputGroupText("mean = ", style={"font-weight": "bold"}),
                    dbc.Input(
                        id="page-3-group-1-ph-mean",
                        type="number",
                        value=7.4,
                    )
                ], size="lg")),
                dbc.Col(dbc.InputGroup([
                    dbc.InputGroupText("std = ", style={"font-weight": "bold"}),
                    dbc.Input(
                        id="page-3-group-1-ph-std",
                        type="number",
                        value=1.5
                    )
                ], size="lg")),
            ]
        ),
    ],
    style={"display": "none"}
)

page_3_group_1_set_pH_card = dbc.Card(
    [
        dbc.CardHeader(html.H6("Group 1 (Normal Group)",
                               style={"font-weight": "bold", "color": "steelblue"})),
        dbc.CardBody(
            [
                html.P("Set the same pH for all replicates: ", style={"font-weight": "bold"}),
                dbc.RadioItems(
                    id="page-3-group-1-ph-same-flag",
                    options=[
                        {"label": "True", "value": "true"},
                        {"label": "False", "value": "false"},
                    ],
                    value="false",
                    inline=True,
                ),
                html.Br(),
                page_3_group_1_pH_same_div,
                page_3_group_1_pH_not_same_div,
                html.Br(),
            ]
        ),
    ]
)

page_3_group_2_pH_same_div = html.Div(
    id="page-3-group-2-ph-same-div",
    children=[
        html.P("Please input the pH:", style={"font-style": 'italic'}),
        dbc.InputGroup([
            dbc.InputGroupText("pH = ", style={"font-weight": "bold"}),
            dbc.Input(
                id="page-3-group-2-same-ph",
                type="number",
                value=6.5,
            )
        ], size="lg")
    ],
    style={"display": "none"}
)

page_3_group_2_pH_not_same_div = html.Div(
    id="page-3-group-2-ph-not-same-div",
    children=[
        html.P("Please input the mean and std for pH:", style={"font-style": 'italic'}),
        dbc.Row(
            [
                dbc.Col(dbc.InputGroup([
                    dbc.InputGroupText("mean = ", style={"font-weight": "bold"}),
                    dbc.Input(
                        id="page-3-group-2-ph-mean",
                        type="number",
                        value=7.4,
                    )
                ], size="lg")),
                dbc.Col(dbc.InputGroup([
                    dbc.InputGroupText("std = ", style={"font-weight": "bold"}),
                    dbc.Input(
                        id="page-3-group-2-ph-std",
                        type="number",
                        value=1.5
                    )
                ], size="lg")),
            ]
        ),
    ],
    style={"display": "none"}
)

page_3_group_2_set_pH_card = dbc.Card(
    [
        dbc.CardHeader(html.H6("Group 2 (Abnormal Group)",
                               style={"font-weight": "bold", "color": "steelblue"})),
        dbc.CardBody(
            [
                html.P("Set the same pH for all replicates: ", style={"font-weight": "bold"}),
                dbc.RadioItems(
                    id="page-3-group-2-ph-same-flag",
                    options=[
                        {"label": "True", "value": "true"},
                        {"label": "False", "value": "false"},
                    ],
                    value="false",
                    inline=True,
                ),
                html.Br(),
                page_3_group_2_pH_same_div,
                page_3_group_2_pH_not_same_div,
                html.Br(),
            ]
        ),
    ]
)

page_3_group_cons_and_ph_table_div = html.Div(
    id="page-3-group-cons-ph-table-div",
    children=[
        dash_table.DataTable(id="page-3-group-cons-ph-table", export_format="csv", export_headers="display", editable=True,
                             virtualization=True, merge_duplicate_headers=True,
                             style_cell={
                               'font_family': 'arial',
                               'font_size': '16px',
                               'text_align': 'center'
                             },
                             style_data_conditional=[
                               {
                                   "if": {"column_id": "meta_name"},
                                   "fontWeight": "bold",
                               },
                             ],
                             fixed_rows={'headers': True},
                             style_table={'height': 250},
                             css=[
                               {"selector": ".dash-spreadsheet tr th", "rule": "height: 35px;"},
                               {"selector": ".dash-spreadsheet tr td", "rule": "height: 35px;"},
                               {"selector": ".export", "rule": "position: absolute; right:0%; top:-5px; "
                                                               "background-color: lightsteelblue; "
                                                               "font-weight: 300; font-size: 1.5rem;"},
                             ],
                             ),
        html.Br(),
        html.Br(),
    ],
    style={"display": "none"}
)

page_3_group_peak_shift_part_2_div = html.Div(
    id="page-3-group-peak-shift-part-2-div",
    children=[
        html.H6("Part II: simulate spectra with peak shift", style={"font-weight": "bold"}),
        dbc.Card(
            id="page-3-group-peak-shift-ad-param-setting-card",
            children=[dbc.CardBody([
                html.H6("Advanced Parameters Setting", style={"font-weight": "bold", "color": "steelblue"}),
                html.Br(),
                dbc.Row(
                    [
                        dbc.Col(dbc.Card(
                            [
                                dbc.CardHeader(html.P("Remove exclusion regions", style={"font-weight": "bold"})),
                                dbc.CardBody(
                                    html.Div(
                                        [
                                            html.P("Range of calibration peaks (ppm)", style={"font-style": 'italic'}),
                                            dcc.RangeSlider(
                                                id="page-3-group-peak-shift-calibration-range",
                                                min=-0.3, max=0.3, step=0.1,
                                                value=[-0.3, 0.3],
                                            ),
                                            html.Br(),
                                            html.P("Range of water suppression (ppm)", style={"font-style": 'italic'}),
                                            dcc.RangeSlider(
                                                id="page-3-group-peak-shift-water-range",
                                                min=4.5, max=5.5, step=0.1,
                                                value=[4.5, 5.0]
                                            )
                                        ]
                                    ),
                                )
                            ], style={'height': '25vh'}, color="primary", outline=True
                        )),
                        dbc.Col(dbc.Card(
                            [
                                dbc.CardHeader(html.P("Filter noise", style={"font-weight": "bold"})),
                                dbc.CardBody(
                                    html.Div(
                                        [
                                            html.P("Parameters for filtering noise", style={"font-style": 'italic'}),
                                            dbc.InputGroup([
                                                dbc.InputGroupText("threshold = ", style={"font-weight": "bold"}),
                                                dbc.Select(
                                                    id="page-3-group-peak-shift-filter-thres",
                                                    options=[
                                                        {"label": "0.01", "value": 0.01},
                                                        {"label": "0.05", "value": 0.05},
                                                    ],
                                                    value=0.05,
                                                )
                                            ], size="lg"),
                                        ]
                                    ),
                                )
                            ], style={'height': '25vh'}, color="success", outline=True
                        )),
                        dbc.Col(dbc.Card(
                            [
                                dbc.CardHeader(html.P("Smooth spectrum", style={"font-weight": "bold"})),
                                dbc.CardBody(
                                    html.Div(
                                        [
                                            html.P("Thresholds for smoothing", style={"font-style": 'italic'}),
                                            dbc.InputGroup([
                                                dbc.InputGroupText("m = ", style={"font-weight": "bold"}),
                                                dbc.Input(
                                                    id="page-3-group-peak-shift-smooth-thres-m",
                                                    type="number",
                                                    value=1,
                                                )
                                            ], size="lg"),
                                            html.Br(),
                                            dbc.InputGroup([
                                                dbc.InputGroupText("n = ", style={"font-weight": "bold"}),
                                                dbc.Input(
                                                    id="page-3-group-peak-shift-smooth-thres-n",
                                                    type="number",
                                                    value=1,
                                                )
                                            ], size="lg"),
                                        ]
                                    ),
                                )
                            ], style={'height': '25vh'}, color="danger", outline=True
                        )),
                        dbc.Col(dbc.Card(
                            [
                                dbc.CardHeader(html.P("Add noise", style={"font-weight": "bold"})),
                                dbc.CardBody(
                                    html.Div(
                                        [
                                            html.P("SNR for adding noise", style={"font-style": 'italic'}),
                                            dbc.RadioItems(
                                                id="page-3-group-peak-shift-snr-noise",
                                                options=[
                                                    {"label": "1000", "value": 1000},
                                                    {"label": "5000", "value": 5000},
                                                    {"label": "10000", "value": 10000},
                                                ],
                                                value=1000,
                                                inline=True,
                                            ),
                                        ]
                                    ),
                                )
                            ], style={'height': '25vh'}, color="warning", outline=True
                        )),
                    ]
                ),
                html.Br(),
                html.Div(
                    [
                        dbc.Button(id='page-3-group-peak-shift-confirm-param-btn', children='Simulate spectra', n_clicks=0,
                                   style={"font-weight": "bold"}, size='lg', color='primary')
                    ],
                    className="d-grid gap-2 d-md-flex justify-content-md-center",
                ),
                dcc.Store(id='page-3-group-peak-shift-processed-spectra-data-dict'),
                html.Br(),
            ])]
        )
    ],
    style={"display": "none"},
)

page_3_group_peak_shift_group_1_results_div = dbc.Card(
    id="page-3-group-peak-shift-group-1-spectra-div",
    children=[
        dbc.CardBody(
            [
                html.H6("Results of Normal Group (Group 1)", style={"font-weight": "bold", "color": "steelblue"}),
                # html.Br(),
                dcc.Loading(children=[
                    dcc.Dropdown(id='page-3-group-peak-shift-select-repli-1', multi=False,),
                    dcc.Graph(id="page-3-group-peak-shift-group-1-spectra-fig"),
                    dcc.Store(id='page-3-group-peak-shift-group-1-final-dict'),
                    # html.P("Contour Level Slider:"),
                    # dcc.RangeSlider(id='page-3-group-peak-shift-contour-level-1', allowCross=False,
                    #                 min=0, max=100),
                    # html.Div(id='page-3-group-peak-shift-update-level-container-1', style={'margin-top': 10}),
                    html.Br(),
                    # dcc.Graph(id="page-3-group-peak-shift-group-1-stacked-fig"),
                    # html.P("Vertical Spacing Slider:"),
                    # dcc.Slider(
                    #     id='page-3-group-peak-shift-group-1-vs-slider',
                    #     min=0,
                    #     max=3,
                    #     step=0.01,
                    #     marks={i: '{}'.format(10 ** i) for i in range(4)},
                    #     value=1,
                    # ),
                    # html.Div(id='page-3-group-peak-shift-update-vs-container-1', style={'margin-top': 10}),
                    html.Br(),
                    html.Button("Download simulated spectra as csv", id="page-3-group-peak-shift-group-1-btn-csv"),
                    dcc.Download(id="page-3-group-peak-shift-group-1-download-spectra-csv"),
                ]),
            ]
        )
    ]
)

page_3_group_peak_shift_group_2_results_div = dbc.Card(
    id="page-3-group-peak-shift-group-2-spectra-div",
    children=[
        dbc.CardBody(
            [
                html.H6("Results of Abnormal Group (Group 2)", style={"font-weight": "bold", "color": "steelblue"}),
                # html.Br(),
                dcc.Loading(children=[
                    dcc.Dropdown(id='page-3-group-peak-shift-select-repli-2', multi=False,),
                    dcc.Graph(id="page-3-group-peak-shift-group-2-spectra-fig"),
                    dcc.Store(id='page-3-group-peak-shift-group-2-final-dict'),
                    # html.P("Contour Level Slider:"),
                    # dcc.RangeSlider(id='page-3-group-peak-shift-contour-level-2', allowCross=False,
                    #                 min=0, max=100),
                    # html.Div(id='page-3-group-peak-shift-update-level-container-2', style={'margin-top': 10}),
                    html.Br(),
                    # dcc.Graph(id="page-3-group-peak-shift-group-2-stacked-fig"),
                    # html.P("Vertical Spacing Slider:"),
                    # dcc.Slider(
                    #     id='page-3-group-peak-shift-group-2-vs-slider',
                    #     min=0,
                    #     max=3,
                    #     step=0.01,
                    #     marks={i: '{}'.format(10 ** i) for i in range(4)},
                    #     value=1,
                    # ),
                    # html.Div(id='page-3-group-peak-shift-update-vs-container-2', style={'margin-top': 10}),
                    html.Br(),
                    html.Button("Download simulated spectra as csv", id="page-3-group-peak-shift-group-2-btn-csv"),
                    dcc.Download(id="page-3-group-peak-shift-group-2-download-spectra-csv"),
                ]),
            ]
        )
    ]
)

page_3_group_peak_shift_spectra_div = html.Div(
    id="page-3-group-peak-shift-spectra-div",
    children=[
        dbc.Card(
            id="page-3-group-peak-shift-all-spectra-card",
            children=[
                dbc.CardBody(
                    [
                        html.H6("Metabolite spectrum in the selected mixture",
                                style={"font-weight": "bold", "color": "steelblue"}),
                        html.Br(),
                        dcc.Loading(children=[
                            dcc.Dropdown(id='page-3-group-peak-shift-select-mix', multi=False,),
                            dcc.Graph(id='page-3-group-peak-shift-all-meta-fig')
                        ], type="circle"),
                    ]
                )
            ]
        ),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(page_3_group_peak_shift_group_1_results_div),
                dbc.Col(page_3_group_peak_shift_group_2_results_div),
            ]
        ),
    ],
    style={"display": "none"},
)

page_3_group_peak_shift_results_div = html.Div(
    id="page-3-group-peak-shift-results-div",
    children=[
        html.H6("Part I: define the pH for each group", style={"font-weight": "bold"}),
        dbc.Row(
            [
                dbc.Col(page_3_group_1_set_pH_card),
                dbc.Col(page_3_group_2_set_pH_card),
            ]
        ),
        html.Br(),
        html.Div([
                    dbc.Button(id='page-3-tab-3-group-submit-part-1-btn', children='OK', n_clicks=0, style={"font-weight": "bold"},
                               size='lg', color='secondary'),
        ]),
        html.Br(),
        html.Br(),
        page_3_group_cons_and_ph_table_div,
        page_3_group_peak_shift_part_2_div,
        html.Br(),
        page_3_group_peak_shift_spectra_div,
        html.Br(),
    ],
    style={"display": "none"},
)

# ---------- tab 3: continuous no peak shift results div ----------
page_3_conti_advance_parameters_setting_card = dbc.Card(
    id="page-3-conti-ad-param-setting-card",
    children=[
        dbc.CardBody([
        html.H6("Advanced Parameters Setting", style={"font-weight": "bold", "color": "steelblue"}),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(dbc.Card(
                    [
                        dbc.CardHeader(html.P("Remove exclusion regions", style={"font-weight": "bold"})),
                        dbc.CardBody(
                            html.Div(
                                [
                                    html.P("Range of calibration peaks (ppm)", style={"font-style": 'italic'}),
                                    dcc.RangeSlider(
                                        id="page-3-conti-calibration-range",
                                        min=-0.3, max=0.3, step=0.1,
                                        value=[-0.3, 0.3],
                                    ),
                                    html.Br(),
                                    html.P("Range of water suppression (ppm)", style={"font-style": 'italic'}),
                                    dcc.RangeSlider(
                                        id="page-3-conti-water-range",
                                        min=4.5, max=5.5, step=0.1,
                                        value=[4.5, 5.0]
                                    )
                                ]
                            ),
                        )
                    ], style={'height': '25vh'}, color="primary", outline=True
                )),
                dbc.Col(dbc.Card(
                    [
                        dbc.CardHeader(html.P("Filter noise", style={"font-weight": "bold"})),
                        dbc.CardBody(
                            html.Div(
                                [
                                    html.P("Parameters for filtering noise", style={"font-style": 'italic'}),
                                    dbc.InputGroup([
                                        dbc.InputGroupText("threshold = ", style={"font-weight": "bold"}),
                                        dbc.Select(
                                            id="page-3-conti-filter-thres",
                                            options=[
                                                {"label": "0.01", "value": 0.01},
                                                {"label": "0.05", "value": 0.05},
                                            ],
                                            value=0.05,
                                        )
                                    ], size="lg"),
                                ]
                            ),
                        )
                    ], style={'height': '25vh'}, color="success", outline=True
                )),
                dbc.Col(dbc.Card(
                    [
                        dbc.CardHeader(html.P("Smooth spectrum", style={"font-weight": "bold"})),
                        dbc.CardBody(
                            html.Div(
                                [
                                    html.P("Thresholds for smoothing", style={"font-style": 'italic'}),
                                    dbc.InputGroup([
                                        dbc.InputGroupText("m = ", style={"font-weight": "bold"}),
                                        dbc.Input(
                                            id="page-3-conti-smooth-thres-m",
                                            type="number",
                                            value=1,
                                        )
                                    ], size="lg"),
                                    html.Br(),
                                    dbc.InputGroup([
                                        dbc.InputGroupText("n = ", style={"font-weight": "bold"}),
                                        dbc.Input(
                                            id="page-3-conti-smooth-thres-n",
                                            type="number",
                                            value=1,
                                        )
                                    ], size="lg"),
                                ]
                            ),
                        )
                    ], style={'height': '25vh'}, color="danger", outline=True
                )),
                dbc.Col(dbc.Card(
                    [
                        dbc.CardHeader(html.P("Add noise", style={"font-weight": "bold"})),
                        dbc.CardBody(
                            html.Div(
                                [
                                    html.P("SNR for adding noise", style={"font-style": 'italic'}),
                                    dbc.RadioItems(
                                        id="page-3-conti-snr-noise",
                                        options=[
                                            {"label": "1000", "value": 1000},
                                            {"label": "5000", "value": 5000},
                                            {"label": "10000", "value": 10000},
                                        ],
                                        value=1000,
                                        inline=True,
                                    ),
                                ]
                            ),
                        )
                    ], style={'height': '25vh'}, color="warning", outline=True
                )),
            ]
        ),
        html.Br(),
        html.Div(
            [
                dbc.Button(id='page-3-conti-confirm-param-btn', children='Simulate spectra', n_clicks=0,
                           style={"font-weight": "bold"}, size='lg', color='primary')
            ],
            className="d-grid gap-2 d-md-flex justify-content-md-center",
        ),
        dcc.Store(id='page-3-conti-processed-spectra-data-dict'),
        html.Br(),
    ])]
)

page_3_conti_all_metabolites_spectra_fig_div = dbc.Card(
    id='page-3-conti-all-spectra-div',
    children=[
        dbc.CardBody(
            [
                html.H6("Metabolite spectrum in the selected mixture", style={"font-weight": "bold", "color": "steelblue"}),
                html.Br(),
                dcc.Loading(children=[
                    dcc.Dropdown(id='page-3-conti-select-mix', multi=False,),
                    dcc.Graph(id='page-3-conti-all-meta-fig')
                ]),
            ]
        )
    ],
)

page_3_conti_meta_y_plot_div = dbc.Card(
    id="page-3-conti-meta-y-plot-div",
    children=[
        dcc.Loading(dcc.Graph(id='page-3-conti-meta-y-fig')),
        html.Br(),
    ],
    style={'height': '72vh'}
)

page_3_conti_mix_spectra_div = dbc.Card(
    id="page-3-conti-mix-spectra-div",
    children=[
        dbc.CardBody(
        [
            dcc.Loading(children=[
                dcc.Store(id="page-3-conti-final-replicate-dict"),
                dcc.Dropdown(
                    id="page-3-conti-select-replicate",
                    multi=False
                ),
                dcc.Graph(id="page-3-conti-mix-spectra-fig"),
                # html.P("Contour Level Slider:"),
                # dcc.RangeSlider(id='page-3-conti-contour-level-1', allowCross=False,
                #                 min=0, max=100),
                # html.Div(id='page-3-conti-update-level-container-1', style={'margin-top': 10}),
                html.Br(),
            ]),
        ]
        )
    ],
    style={'height': '72vh'}
)

page_3_conti_no_peak_shift_results_div = html.Div(
    id="page-3-conti-no-peak-shift-results-div",
    children=[
        page_3_conti_advance_parameters_setting_card,
        html.Br(),
        page_3_conti_all_metabolites_spectra_fig_div,
        html.Br(),
        dbc.Row(
            [
                dbc.Col(page_3_conti_meta_y_plot_div),
                dbc.Col(page_3_conti_mix_spectra_div),
            ]
        ),
        html.Br(),
        html.Div(children=[html.Button("Download simulated spectra as csv", id="page-3-conti-btn-csv")],
                 className="d-grid gap-2 d-md-flex justify-content-md-left"),
        dcc.Download(id="page-3-conti-download-spectra-csv"),
    ],
    style={"display": "none"},
)

# ---------- tab 3: continuous with peak shift results div ----------
page_3_conti_pH_same_div = html.Div(
    id="page-3-conti-ph-same-div",
    children=[
        html.P("Please input the pH:", style={"font-style": 'italic'}),
        dbc.InputGroup([
            dbc.InputGroupText("pH = ", style={"font-weight": "bold"}),
            dbc.Input(
                id="page-3-conti-same-ph",
                type="number",
                value=6.5,
            )
        ], size="lg")
    ],
    style={"display": "none"}
)

page_3_conti_pH_not_same_div = html.Div(
    id="page-3-conti-ph-not-same-div",
    children=[
        html.P("Please input the mean and std for pH:", style={"font-style": 'italic'}),
        dbc.Row(
            [
                dbc.Col(dbc.InputGroup([
                    dbc.InputGroupText("mean = ", style={"font-weight": "bold"}),
                    dbc.Input(
                        id="page-3-conti-ph-mean",
                        type="number",
                        value=7.4,
                    )
                ], size="lg")),
                dbc.Col(dbc.InputGroup([
                    dbc.InputGroupText("std = ", style={"font-weight": "bold"}),
                    dbc.Input(
                        id="page-3-conti-ph-std",
                        type="number",
                        value=1.5
                    )
                ], size="lg")),
            ]
        ),
    ],
    style={"display": "none"}
)

page_3_conti_set_pH_card = dbc.Card(
    [
    dbc.CardBody(
        [
            html.P("Set the same pH for all replicates: ", style={"font-weight": "bold"}),
            dbc.RadioItems(
                id="page-3-conti-ph-same-flag",
                options=[
                    {"label": "True", "value": "true"},
                    {"label": "False", "value": "false"},
                ],
                value="false",
                inline=True,
            ),
            html.Br(),
            page_3_conti_pH_same_div,
            page_3_conti_pH_not_same_div,
            html.Br(),
        ]
    )], className="w-75 mb-3",
)

page_3_conti_cons_and_ph_table_div = html.Div(
    id="page-3-conti-cons-ph-table-div",
    children=[
        dash_table.DataTable(id="page-3-conti-cons-ph-table", export_format="csv", export_headers="display",
                             editable=True, virtualization=True,
                             style_cell={
                               'font_family': 'arial',
                               'font_size': '16px',
                               'text_align': 'center',
                               'minWidth': '180px', 'width': '180px', 'maxWidth': '180px',
                             },
                             style_header={
                                'backgroundColor': 'whitesmoke',
                                'fontWeight': 'bold'
                             },
                             style_data_conditional=[
                               {
                                   "if": {"column_id": "meta_name"},
                                   "fontWeight": "bold",
                               },
                             ],
                             fixed_rows={'headers': True},
                             style_table={'height': 250},
                             css=[
                               {"selector": ".dash-spreadsheet tr th", "rule": "height: 35px;"},
                               {"selector": ".dash-spreadsheet tr td", "rule": "height: 35px;"},
                               {"selector": ".export", "rule": "position: absolute; right:0%; top:-40px; "
                                                               "background-color: lightsteelblue; "
                                                               "font-weight: 300; font-size: 1.5rem;"},
                             ]),
        html.Br(),
        html.Br(),
    ],
    style={"display": "none"}
)

page_3_conti_peak_shift_part_2_div = html.Div(
    id="page-3-conti-peak-shift-part-2-div",
    children=[
        html.H6("Part II: simulate spectra with peak shift", style={"font-weight": "bold"}),
        dbc.Card(
            id="page-3-conti-peak-shift-ad-param-setting-card",
            children=[
                dbc.CardBody([
                html.H6("Advanced Parameters Setting", style={"font-weight": "bold", "color": "steelblue"}),
                html.Br(),
                dbc.Row(
                    [
                        dbc.Col(dbc.Card(
                            [
                                dbc.CardHeader(html.P("Remove exclusion regions", style={"font-weight": "bold"})),
                                dbc.CardBody(
                                    html.Div(
                                        [
                                            html.P("Range of calibration peaks (ppm)", style={"font-style": 'italic'}),
                                            dcc.RangeSlider(
                                                id="page-3-conti-peak-shift-calibration-range",
                                                min=-0.3, max=0.3, step=0.1,
                                                value=[-0.3, 0.3],
                                            ),
                                            html.Br(),
                                            html.P("Range of water suppression (ppm)", style={"font-style": 'italic'}),
                                            dcc.RangeSlider(
                                                id="page-3-conti-peak-shift-water-range",
                                                min=4.5, max=5.5, step=0.1,
                                                value=[4.5, 5.0]
                                            )
                                        ]
                                    ),
                                )
                            ], style={'height': '25vh'}, color="primary", outline=True
                        )),
                        dbc.Col(dbc.Card(
                            [
                                dbc.CardHeader(html.P("Filter noise", style={"font-weight": "bold"})),
                                dbc.CardBody(
                                    html.Div(
                                        [
                                            html.P("Parameters for filtering noise", style={"font-style": 'italic'}),
                                            dbc.InputGroup([
                                                dbc.InputGroupText("threshold = ", style={"font-weight": "bold"}),
                                                dbc.Select(
                                                    id="page-3-conti-peak-shift-filter-thres",
                                                    options=[
                                                        {"label": "0.01", "value": 0.01},
                                                        {"label": "0.05", "value": 0.05},
                                                    ],
                                                    value=0.05,
                                                )
                                            ], size="lg"),
                                        ]
                                    ),
                                )
                            ], style={'height': '25vh'}, color="success", outline=True
                        )),
                        dbc.Col(dbc.Card(
                            [
                                dbc.CardHeader(html.P("Smooth spectrum", style={"font-weight": "bold"})),
                                dbc.CardBody(
                                    html.Div(
                                        [
                                            html.P("Thresholds for smoothing", style={"font-style": 'italic'}),
                                            dbc.InputGroup([
                                                dbc.InputGroupText("m = ", style={"font-weight": "bold"}),
                                                dbc.Input(
                                                    id="page-3-conti-peak-shift-smooth-thres-m",
                                                    type="number",
                                                    value=1,
                                                )
                                            ], size="lg"),
                                            html.Br(),
                                            dbc.InputGroup([
                                                dbc.InputGroupText("n = ", style={"font-weight": "bold"}),
                                                dbc.Input(
                                                    id="page-3-conti-peak-shift-smooth-thres-n",
                                                    type="number",
                                                    value=1,
                                                )
                                            ], size="lg"),
                                        ]
                                    ),
                                )
                            ], style={'height': '25vh'}, color="danger", outline=True
                        )),
                        dbc.Col(dbc.Card(
                            [
                                dbc.CardHeader(html.P("Add noise", style={"font-weight": "bold"})),
                                dbc.CardBody(
                                    html.Div(
                                        [
                                            html.P("SNR for adding noise", style={"font-style": 'italic'}),
                                            dbc.RadioItems(
                                                id="page-3-conti-peak-shift-snr-noise",
                                                options=[
                                                    {"label": "1000", "value": 1000},
                                                    {"label": "5000", "value": 5000},
                                                    {"label": "10000", "value": 10000},
                                                ],
                                                value=1000,
                                                inline=True,
                                            ),
                                        ]
                                    ),
                                )
                            ], style={'height': '25vh'}, color="warning", outline=True
                        )),
                    ]
                ),
                html.Br(),
                html.Div(
                    [
                        dbc.Button(id='page-3-conti-peak-shift-confirm-param-btn', children='Simulate spectra', n_clicks=0,
                                   style={"font-weight": "bold"}, size='lg', color='primary')
                    ],
                    className="d-grid gap-2 d-md-flex justify-content-md-center",
                ),
                dcc.Store(id='page-3-conti-peak-shift-processed-spectra-data-dict'),
                html.Br(),
            ])]
        )
    ]
)

page_3_conti_peak_shift_spectra_div = html.Div(
    id="page-3-conti-peak-shift-spectra-div",
    children=[
        dbc.Card(
            id='page-3-conti-peak-shift-all-spectra-div',
            children=[
                dbc.CardBody(
                    [
                        html.H6("Metabolite spectrum in the selected mixture", style={"font-weight": "bold", "color": "steelblue"}),
                        html.Br(),
                        dcc.Loading(children=[
                            dcc.Dropdown(id='page-3-conti-peak-shift-select-mix', multi=False,),
                            dcc.Graph(id='page-3-conti-peak-shift-all-meta-fig')
                        ]),
                    ]
                )
            ],
        ),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    dbc.Card(
                        id="page-3-conti-peak-shift-meta-y-plot-div",
                        children=[
                            dcc.Loading(dcc.Graph(id='page-3-conti-peak-shift-meta-y-fig')),
                            html.Br(),
                        ],
                        style={'height': '80vh'}
                    ),
                ),
                dbc.Col(
                    dbc.Card(
                        id="page-3-conti-peak-shift-mix-spectra-div",
                        children=[
                            dbc.CardBody(
                                [
                                    dcc.Loading(children=[
                                        dcc.Store(id="page-3-conti-peak-shift-final-replicate-dict"),
                                        dcc.Dropdown(
                                            id="page-3-conti-peak-shift-select-replicate",
                                            multi=False
                                        ),
                                        dcc.Graph(id="page-3-conti-peak-shift-mix-spectra-fig"),
                                        # html.P("Contour Level Slider:"),
                                        # dcc.RangeSlider(id='page-3-conti-peak-shift-contour-level-1', allowCross=False,
                                        #                 min=0, max=100),
                                        # html.Div(id='page-3-conti-peak-shift-update-level-container-1', style={'margin-top': 10}),
                                        html.Br(),
                                    ]),
                                ]
                            )
                        ],
                        style={'height': '80vh'}
                    )
                )
            ]
        ),
        html.Br(),
        html.Button("Download simulated spectra as csv", id="page-3-conti-peak-shift-btn-csv"),
        dcc.Download(id="page-3-conti-peak-shift-download-spectra-csv"),
    ],
    style={"display": "none"}
)

page_3_conti_peak_shift_results_div = html.Div(
    id="page-3-conti-peak-shift-results-div",
    children=[
        html.H6("Part I: define the pH for all replicates", style={"font-weight": "bold"}),
        page_3_conti_set_pH_card,
        html.Br(),
        html.Div([
                    dbc.Button(id='page-3-tab-3-conti-submit-part-1-btn', children='OK', n_clicks=0, style={"font-weight": "bold"},
                               size='lg', color='secondary'),
        ]),
        html.Br(),
        html.Br(),
        page_3_conti_cons_and_ph_table_div,
        page_3_conti_peak_shift_part_2_div,
        html.Br(),
        page_3_conti_peak_shift_spectra_div,
        html.Br(),
    ],
    style={"display": "none"},
)

# ---------- tab 3: upload group no peak shift results div ----------
page_3_upload_group_advance_parameters_setting_card = dbc.Card(
    id="page-3-upload-group-no-peak-shift-ad-param-setting-card",
    children=[dbc.CardBody([
        html.H6("Advanced Parameters Setting", style={"font-weight": "bold", "color": "steelblue"}),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(dbc.Card(
                    [
                        dbc.CardHeader(html.P("Remove exclusion regions", style={"font-weight": "bold"})),
                        dbc.CardBody(
                            html.Div(
                                [
                                    html.P("Range of calibration peaks (ppm)", style={"font-style": 'italic'}),
                                    dcc.RangeSlider(
                                        id="page-3-upload-group-no-peak-shift-calibration-range",
                                        min=-0.3, max=0.3, step=0.1,
                                        value=[-0.3, 0.3],
                                    ),
                                    html.Br(),
                                    html.P("Range of water suppression (ppm)", style={"font-style": 'italic'}),
                                    dcc.RangeSlider(
                                        id="page-3-upload-group-no-peak-shift-water-range",
                                        min=4.5, max=5.5, step=0.1,
                                        value=[4.5, 5.0]
                                    )
                                ]
                            ),
                        )
                    ], style={'height': '25vh'}, color="primary", outline=True
                )),
                dbc.Col(dbc.Card(
                    [
                        dbc.CardHeader(html.P("Filter noise", style={"font-weight": "bold"})),
                        dbc.CardBody(
                            html.Div(
                                [
                                    html.P("Parameters for filtering noise", style={"font-style": 'italic'}),
                                    dbc.InputGroup([
                                        dbc.InputGroupText("threshold = ", style={"font-weight": "bold"}),
                                        dbc.Select(
                                            id="page-3-upload-group-no-peak-shift-filter-thres",
                                            options=[
                                                {"label": "0.01", "value": 0.01},
                                                {"label": "0.05", "value": 0.05},
                                            ],
                                            value=0.05,
                                        )
                                    ], size="lg"),
                                ]
                            ),
                        )
                    ], style={'height': '25vh'}, color="success", outline=True
                )),
                dbc.Col(dbc.Card(
                    [
                        dbc.CardHeader(html.P("Smooth spectrum", style={"font-weight": "bold"})),
                        dbc.CardBody(
                            html.Div(
                                [
                                    html.P("Thresholds for smoothing", style={"font-style": 'italic'}),
                                    dbc.InputGroup([
                                        dbc.InputGroupText("m = ", style={"font-weight": "bold"}),
                                        dbc.Input(
                                            id="page-3-upload-group-no-peak-shift-smooth-thres-m",
                                            type="number",
                                            value=1,
                                        )
                                    ], size="lg"),
                                    html.Br(),
                                    dbc.InputGroup([
                                        dbc.InputGroupText("n = ", style={"font-weight": "bold"}),
                                        dbc.Input(
                                            id="page-3-upload-group-no-peak-shift-smooth-thres-n",
                                            type="number",
                                            value=1,
                                        )
                                    ], size="lg"),
                                ]
                            ),
                        )
                    ], style={'height': '25vh'}, color="danger", outline=True
                )),
                dbc.Col(dbc.Card(
                    [
                        dbc.CardHeader(html.P("Add noise", style={"font-weight": "bold"})),
                        dbc.CardBody(
                            html.Div(
                                [
                                    html.P("SNR for adding noise", style={"font-style": 'italic'}),
                                    dbc.RadioItems(
                                        id="page-3-upload-group-no-peak-shift-snr-noise",
                                        options=[
                                            {"label": "1000", "value": 1000},
                                            {"label": "5000", "value": 5000},
                                            {"label": "10000", "value": 10000},
                                        ],
                                        value=1000,
                                        inline=True,
                                    ),
                                ]
                            ),
                        )
                    ], style={'height': '25vh'}, color="warning", outline=True
                )),
            ]
        ),
        html.Br(),
        html.Div(
            [
                dbc.Button(id='page-3-upload-group-no-peak-shift-confirm-param-btn', children='Simulate spectra', n_clicks=0,
                           style={"font-weight": "bold"}, size='lg', color='primary')
            ],
            className="d-grid gap-2 d-md-flex justify-content-md-center",
        ),
        dcc.Store(id='page-3-upload-group-no-peak-shift-processed-spectra-data-dict'),
        html.Br(),
    ])]
)

page_3_upload_group_no_peak_shift_spectra_div = html.Div(
    id="page-3-upload-group-no-peak-shift-spectra-div",
    children=[
        dbc.Card(
            id="page-3-upload-group-no-peak-shift-all-spectra-card",
            children=[
                dbc.CardBody(
                    [
                        html.H6("Metabolite spectrum in the selected mixture",
                                style={"font-weight": "bold", "color": "steelblue"}),
                        html.Br(),
                        dcc.Loading(children=[
                            dcc.Dropdown(id='page-3-upload-group-no-peak-shift-select-mix', multi=False,),
                            dcc.Graph(id="page-3-upload-group-no-peak-shift-all-meta-fig")
                        ])
                    ]
                )
            ]
        ),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    dbc.Card(
                        id="page-3-upload-group-no-peak-shift-group-1-spectra-div",
                        children=[
                            dbc.CardBody(
                                [
                                    html.H6("Results of Normal Group (Group 1)", style={"font-weight": "bold", "color": "steelblue"}),
                                    # html.Br(),
                                    dcc.Loading(children=[
                                        dcc.Dropdown(id='page-3-upload-group-no-peak-shift-select-repli-1', multi=False,),
                                        dcc.Graph(id="page-3-upload-group-no-peak-shift-group-1-spectra-fig"),
                                        dcc.Store(id='page-3-upload-group-no-peak-shift-group-1-final-dict'),
                                        # html.P("Contour Level Slider:"),
                                        # dcc.RangeSlider(
                                        #     id='page-3-upload-group-no-peak-shift-contour-level-1',
                                        #     allowCross=False,
                                        #     min=0, max=100
                                        # ),
                                        # html.Div(id='page-3-upload-group-no-peak-shift-update-level-container-1', style={'margin-top': 10}),
                                        html.Br(),
                                        html.Button("Download simulated spectra as csv", id="page-3-upload-group-no-peak-shift-group-1-btn-csv"),
                                        dcc.Download(id="page-3-upload-group-no-peak-shift-group-1-download-spectra-csv"),]),
                                ]
                            )
                        ]
                    )
                ),
                dbc.Col(
                    dbc.Card(
                        id="page-3-upload-group-no-peak-shift-group-2-spectra-div",
                        children=[
                            dbc.CardBody(
                                [
                                    html.H6("Results of Abnormal Group (Group 2)", style={"font-weight": "bold", "color": "steelblue"}),
                                    # html.Br(),
                                    dcc.Loading(children=[
                                        dcc.Dropdown(id='page-3-upload-group-no-peak-shift-select-repli-2', multi=False,),
                                        dcc.Graph(id="page-3-upload-group-no-peak-shift-group-2-spectra-fig"),
                                        dcc.Store(id='page-3-upload-group-no-peak-shift-group-2-final-dict'),
                                        # html.P("Contour Level Slider:"),
                                        # dcc.RangeSlider(
                                        #     id='page-3-upload-group-no-peak-shift-contour-level-2',
                                        #     allowCross=False,
                                        #     min=0, max=100
                                        # ),
                                        # html.Div(id='page-3-upload-group-no-peak-shift-update-level-container-2', style={'margin-top': 10}),
                                        html.Br(),
                                        html.Button("Download simulated spectra as csv", id="page-3-upload-group-no-peak-shift-group-2-btn-csv"),
                                        dcc.Download(id="page-3-upload-group-no-peak-shift-group-2-download-spectra-csv"),]),
                                ]
                            )
                        ]
                    )
                ),
            ]
        )
    ],
    # style={"display": "none"},
)

page_3_upload_group_no_peak_shift_results_div = html.Div(
    id="page-3-upload-group-no-peak-shift-results-div",
    children=[
        page_3_upload_group_advance_parameters_setting_card,
        html.Br(),
        page_3_upload_group_no_peak_shift_spectra_div,
        html.Br(),
    ],
    style={"display": "none"},
)

# ---------- tab 3: upload group with peak shift results div ----------
page_3_upload_group_1_pH_same_div = html.Div(
    id="page-3-upload-group-1-ph-same-div",
    children=[
        html.P("Please input the pH:", style={"font-style": 'italic'}),
        dbc.InputGroup([
            dbc.InputGroupText("pH = ", style={"font-weight": "bold"}),
            dbc.Input(
                id="page-3-upload-group-1-same-ph",
                type="number",
                value=6.5,
            )
        ], size="lg")
    ],
    style={"display": "none"}
)

page_3_upload_group_1_pH_not_same_div = html.Div(
    id="page-3-upload-group-1-ph-not-same-div",
    children=[
        html.P("Please input the mean and std for pH:", style={"font-style": 'italic'}),
        dbc.Row(
            [
                dbc.Col(dbc.InputGroup([
                    dbc.InputGroupText("mean = ", style={"font-weight": "bold"}),
                    dbc.Input(
                        id="page-3-upload-group-1-ph-mean",
                        type="number",
                        value=7.4,
                    )
                ], size="lg")),
                dbc.Col(dbc.InputGroup([
                    dbc.InputGroupText("std = ", style={"font-weight": "bold"}),
                    dbc.Input(
                        id="page-3-upload-group-1-ph-std",
                        type="number",
                        value=1.5
                    )
                ], size="lg")),
            ]
        ),
    ],
    style={"display": "none"}
)

page_3_upload_group_1_set_pH_card = dbc.Card(
    [
        dbc.CardHeader(html.H6("Group 1 (Normal Group)",
                               style={"font-weight": "bold", "color": "steelblue"})),
        dbc.CardBody(
            [
                html.P("Set the same pH for all replicates: ", style={"font-weight": "bold"}),
                dbc.RadioItems(
                    id="page-3-upload-group-1-ph-same-flag",
                    options=[
                        {"label": "True", "value": "true"},
                        {"label": "False", "value": "false"},
                    ],
                    value="false",
                    inline=True,
                ),
                html.Br(),
                page_3_upload_group_1_pH_same_div,
                page_3_upload_group_1_pH_not_same_div,
                html.Br(),
            ]
        ),
    ]
)

page_3_upload_group_2_pH_same_div = html.Div(
    id="page-3-upload-group-2-ph-same-div",
    children=[
        html.P("Please input the pH:", style={"font-style": 'italic'}),
        dbc.InputGroup([
            dbc.InputGroupText("pH = ", style={"font-weight": "bold"}),
            dbc.Input(
                id="page-3-upload-group-2-same-ph",
                type="number",
                value=6.5,
            )
        ], size="lg")
    ],
    style={"display": "none"}
)

page_3_upload_group_2_pH_not_same_div = html.Div(
    id="page-3-upload-group-2-ph-not-same-div",
    children=[
        html.P("Please input the mean and std for pH:", style={"font-style": 'italic'}),
        dbc.Row(
            [
                dbc.Col(dbc.InputGroup([
                    dbc.InputGroupText("mean = ", style={"font-weight": "bold"}),
                    dbc.Input(
                        id="page-3-upload-group-2-ph-mean",
                        type="number",
                        value=7.4,
                    )
                ], size="lg")),
                dbc.Col(dbc.InputGroup([
                    dbc.InputGroupText("std = ", style={"font-weight": "bold"}),
                    dbc.Input(
                        id="page-3-upload-group-2-ph-std",
                        type="number",
                        value=1.5
                    )
                ], size="lg")),
            ]
        ),
    ],
    style={"display": "none"}
)

page_3_upload_group_2_set_pH_card = dbc.Card(
    [
        dbc.CardHeader(html.H6("Group 2 (Abnormal Group)",
                               style={"font-weight": "bold", "color": "steelblue"})),
        dbc.CardBody(
            [
                html.P("Set the same pH for all replicates: ", style={"font-weight": "bold"}),
                dbc.RadioItems(
                    id="page-3-upload-group-2-ph-same-flag",
                    options=[
                        {"label": "True", "value": "true"},
                        {"label": "False", "value": "false"},
                    ],
                    value="false",
                    inline=True,
                ),
                html.Br(),
                page_3_upload_group_2_pH_same_div,
                page_3_upload_group_2_pH_not_same_div,
                html.Br(),
            ]
        ),
    ]
)

page_3_upload_group_cons_and_ph_table_div = html.Div(
    id="page-3-upload-group-cons-ph-table-div",
    children=[
        dash_table.DataTable(id="page-3-upload-group-cons-ph-table", export_format="csv", export_headers="display", editable=True,
                             virtualization=True, merge_duplicate_headers=True,
                             style_cell={
                               'font_family': 'arial',
                               'font_size': '16px',
                               'text_align': 'center'
                             },
                             style_data_conditional=[
                               {
                                   "if": {"column_id": "meta_name"},
                                   "fontWeight": "bold",
                               },
                             ],
                             fixed_rows={'headers': True},
                             style_table={'height': 250},
                             css=[
                               {"selector": ".dash-spreadsheet tr th", "rule": "height: 35px;"},
                               {"selector": ".dash-spreadsheet tr td", "rule": "height: 35px;"},
                               {"selector": ".export", "rule": "position: absolute; right:0%; top:-40px; "
                                                               "background-color: lightsteelblue; "
                                                               "font-weight: 300; font-size: 1.5rem;"},
                             ],
                             ),
        html.Br(),
        html.Br(),
    ],
    style={"display": "none"}
)

page_3_upload_group_peak_shift_part_2_div = html.Div(
    id="page-3-upload-group-peak-shift-part-2-div",
    children=[
        html.H6("Part II: simulate spectra with peak shift", style={"font-weight": "bold"}),
        dbc.Card(
            id="page-3-upload-group-peak-shift-ad-param-setting-card",
            children=[dbc.CardBody([
                html.H6("Advanced Parameters Setting", style={"font-weight": "bold", "color": "steelblue"}),
                html.Br(),
                dbc.Row(
                    [
                        dbc.Col(dbc.Card(
                            [
                                dbc.CardHeader(html.P("Remove exclusion regions", style={"font-weight": "bold"})),
                                dbc.CardBody(
                                    html.Div(
                                        [
                                            html.P("Range of calibration peaks (ppm)", style={"font-style": 'italic'}),
                                            dcc.RangeSlider(
                                                id="page-3-upload-group-peak-shift-calibration-range",
                                                min=-0.3, max=0.3, step=0.1,
                                                value=[-0.3, 0.3],
                                            ),
                                            html.Br(),
                                            html.P("Range of water suppression (ppm)", style={"font-style": 'italic'}),
                                            dcc.RangeSlider(
                                                id="page-3-upload-group-peak-shift-water-range",
                                                min=4.5, max=5.5, step=0.1,
                                                value=[4.5, 5.0]
                                            )
                                        ]
                                    ),
                                )
                            ], style={'height': '25vh'}, color="primary", outline=True
                        )),
                        dbc.Col(dbc.Card(
                            [
                                dbc.CardHeader(html.P("Filter noise", style={"font-weight": "bold"})),
                                dbc.CardBody(
                                    html.Div(
                                        [
                                            html.P("Parameters for filtering noise", style={"font-style": 'italic'}),
                                            dbc.InputGroup([
                                                dbc.InputGroupText("threshold = ", style={"font-weight": "bold"}),
                                                dbc.Select(
                                                    id="page-3-upload-group-peak-shift-filter-thres",
                                                    options=[
                                                        {"label": "0.01", "value": 0.01},
                                                        {"label": "0.05", "value": 0.05},
                                                    ],
                                                    value=0.05,
                                                )
                                            ], size="lg"),
                                        ]
                                    ),
                                )
                            ], style={'height': '25vh'}, color="success", outline=True
                        )),
                        dbc.Col(dbc.Card(
                            [
                                dbc.CardHeader(html.P("Smooth spectrum", style={"font-weight": "bold"})),
                                dbc.CardBody(
                                    html.Div(
                                        [
                                            html.P("Thresholds for smoothing", style={"font-style": 'italic'}),
                                            dbc.InputGroup([
                                                dbc.InputGroupText("m = ", style={"font-weight": "bold"}),
                                                dbc.Input(
                                                    id="page-3-upload-group-peak-shift-smooth-thres-m",
                                                    type="number",
                                                    value=1,
                                                )
                                            ], size="lg"),
                                            html.Br(),
                                            dbc.InputGroup([
                                                dbc.InputGroupText("n = ", style={"font-weight": "bold"}),
                                                dbc.Input(
                                                    id="page-3-upload-group-peak-shift-smooth-thres-n",
                                                    type="number",
                                                    value=1,
                                                )
                                            ], size="lg"),
                                        ]
                                    ),
                                )
                            ], style={'height': '25vh'}, color="danger", outline=True
                        )),
                        dbc.Col(dbc.Card(
                            [
                                dbc.CardHeader(html.P("Add noise", style={"font-weight": "bold"})),
                                dbc.CardBody(
                                    html.Div(
                                        [
                                            html.P("SNR for adding noise", style={"font-style": 'italic'}),
                                            dbc.RadioItems(
                                                id="page-3-upload-group-peak-shift-snr-noise",
                                                options=[
                                                    {"label": "1000", "value": 1000},
                                                    {"label": "5000", "value": 5000},
                                                    {"label": "10000", "value": 10000},
                                                ],
                                                value=1000,
                                                inline=True,
                                            ),
                                        ]
                                    ),
                                )
                            ], style={'height': '25vh'}, color="warning", outline=True
                        )),
                    ]
                ),
                html.Br(),
                html.Div(
                    [
                        dbc.Button(id='page-3-upload-group-peak-shift-confirm-param-btn', children='Simulate spectra', n_clicks=0,
                                   style={"font-weight": "bold"}, size='lg', color='primary')
                    ],
                    className="d-grid gap-2 d-md-flex justify-content-md-center",
                ),
                dcc.Store(id='page-3-upload-group-peak-shift-processed-spectra-data-dict'),
                html.Br(),
            ])]
        )
    ],
    style={"display": "none"},
)

page_3_upload_group_peak_shift_spectra_div = html.Div(
    id="page-3-upload-group-peak-shift-spectra-div",
    children=[
        dbc.Card(
            id="page-3-upload-group-peak-shift-all-spectra-card",
            children=[
                dbc.CardBody(
                    [
                        html.H6("Metabolite spectrum in the selected mixture",
                                style={"font-weight": "bold", "color": "steelblue"}),
                        html.Br(),
                        # dcc.Graph(id='group-peak-shift-all-meta-fig'),
                        dcc.Loading(children=[
                            dcc.Dropdown(id='page-3-upload-group-peak-shift-select-mix', multi=False,),
                            dcc.Graph(id='page-3-upload-group-peak-shift-all-meta-fig')],
                                    type="default"),
                    ]
                )
            ]
        ),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    dbc.Card(
                        id="page-3-upload-group-peak-shift-group-1-spectra-card",
                        children=[
                            dbc.CardBody(
                                [
                                    html.H6("Results of Normal Group (Group 1)", style={"font-weight": "bold", "color": "steelblue"}),
                                    # html.Br(),
                                    dcc.Loading(children=[
                                        dcc.Dropdown(id='page-3-upload-group-peak-shift-select-repli-1', multi=False,),
                                        dcc.Graph(id="page-3-upload-group-peak-shift-group-1-spectra-fig"),
                                        dcc.Store(id='page-3-upload-group-peak-shift-group-1-final-dict'),
                                        # html.P("Contour Level Slider:"),
                                        # dcc.RangeSlider(id='page-3-upload-group-peak-shift-contour-level-1', allowCross=False,
                                        #                 min=0, max=100),
                                        # html.Div(id='page-3-upload-group-peak-shift-update-level-container-1',
                                        #          style={'margin-top': 10}),
                                        html.Br(),
                                        # dcc.Graph(id="page-3-upload-group-peak-shift-group-1-stacked-fig"),
                                        # html.P("Vertical Spacing Slider:"),
                                        # dcc.Slider(
                                        #     id='page-3-upload-group-peak-shift-group-1-vs-slider',
                                        #     min=0,
                                        #     max=3,
                                        #     step=0.01,
                                        #     marks={i: '{}'.format(10 ** i) for i in range(4)},
                                        #     value=1,
                                        # ),
                                        # html.Div(id='page-3-upload-group-peak-shift-update-vs-container-1', style={'margin-top': 10}),
                                        html.Br(),
                                        html.Button("Download simulated spectra as csv", id="page-3-upload-group-peak-shift-group-1-btn-csv"),
                                        dcc.Download(id="page-3-upload-group-peak-shift-group-1-download-spectra-csv")],
                                                type="default"),
                                ]
                            )
                        ]
                    )
                ),
                dbc.Col(
                    dbc.Card(
                        id="page-3-upload-group-peak-shift-group-2-spectra-card",
                        children=[
                            dbc.CardBody(
                                [
                                    html.H6("Results of Abnormal Group (Group 2)", style={"font-weight": "bold", "color": "steelblue"}),
                                    # html.Br(),
                                    dcc.Loading(children=[
                                        dcc.Dropdown(id='page-3-upload-group-peak-shift-select-repli-2', multi=False, ),
                                        dcc.Graph(id="page-3-upload-group-peak-shift-group-2-spectra-fig"),
                                        dcc.Store(id='page-3-upload-group-peak-shift-group-2-final-dict'),
                                        # html.P("Contour Level Slider:"),
                                        # dcc.RangeSlider(id='page-3-upload-group-peak-shift-contour-level-2',
                                        #                 allowCross=False,
                                        #                 min=0, max=100),
                                        # html.Div(id='page-3-upload-group-peak-shift-update-level-container-2',
                                        #          style={'margin-top': 10}),
                                        html.Br(),
                                        # dcc.Graph(id="page-3-upload-group-peak-shift-group-2-stacked-fig"),
                                        # html.P("Vertical Spacing Slider:"),
                                        # dcc.Slider(
                                        #     id='page-3-upload-group-peak-shift-group-2-vs-slider',
                                        #     min=0,
                                        #     max=3,
                                        #     step=0.01,
                                        #     marks={i: '{}'.format(10 ** i) for i in range(4)},
                                        #     value=1,
                                        # ),
                                        # html.Div(id='page-3-upload-group-peak-shift-update-vs-container-2', style={'margin-top': 10}),
                                        html.Br(),
                                        html.Button("Download simulated spectra as csv", id="page-3-upload-group-peak-shift-group-2-btn-csv"),
                                        dcc.Download(id="page-3-upload-group-peak-shift-group-2-download-spectra-csv")
                                        ],
                                                type="default"),
                                ]
                            )
                        ]
                    )
                )
            ]
        )
    ],
    style={"display": "none"},
)

page_3_upload_group_peak_shift_results_div = html.Div(
    id="page-3-upload-group-peak-shift-results-div",
    children=[
        html.H6("Part I: define the pH for each group", style={"font-weight": "bold"}),
        dbc.Row(
            [
                dbc.Col(page_3_upload_group_1_set_pH_card),
                dbc.Col(page_3_upload_group_2_set_pH_card),
            ]
        ),
        html.Br(),
        html.Div([
                    dbc.Button(id='page-3-tab-3-upload-group-submit-part-1-btn', children='OK', n_clicks=0, style={"font-weight": "bold"},
                               size='lg', color='secondary'),
        ]),
        html.Br(),
        html.Br(),
        page_3_upload_group_cons_and_ph_table_div,
        page_3_upload_group_peak_shift_part_2_div,
        html.Br(),
        page_3_upload_group_peak_shift_spectra_div,
        html.Br(),
    ],
    style={"display": "none"},
)

# ---------- tab 3: upload continuous no peak shift results div ----------
page_3_upload_conti_advance_parameters_setting_card = dbc.Card(
    id="page-3-upload-conti-ad-param-setting-card",
    children=[dbc.CardBody([
        html.H6("Advanced Parameters Setting", style={"font-weight": "bold", "color": "steelblue"}),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(dbc.Card(
                    [
                        dbc.CardHeader(html.P("Remove exclusion regions", style={"font-weight": "bold"})),
                        dbc.CardBody(
                            html.Div(
                                [
                                    html.P("Range of calibration peaks (ppm)", style={"font-style": 'italic'}),
                                    dcc.RangeSlider(
                                        id="page-3-upload-conti-calibration-range",
                                        min=-0.3, max=0.3, step=0.1,
                                        value=[-0.3, 0.3],
                                    ),
                                    html.Br(),
                                    html.P("Range of water suppression (ppm)", style={"font-style": 'italic'}),
                                    dcc.RangeSlider(
                                        id="page-3-upload-conti-water-range",
                                        min=4.5, max=5.5, step=0.1,
                                        value=[4.5, 5.0]
                                    )
                                ]
                            ),
                        )
                    ], style={'height': '25vh'}, color="primary", outline=True
                )),
                dbc.Col(dbc.Card(
                    [
                        dbc.CardHeader(html.P("Filter noise", style={"font-weight": "bold"})),
                        dbc.CardBody(
                            html.Div(
                                [
                                    html.P("Parameters for filtering noise", style={"font-style": 'italic'}),
                                    dbc.InputGroup([
                                        dbc.InputGroupText("threshold = ", style={"font-weight": "bold"}),
                                        dbc.Select(
                                            id="page-3-upload-conti-filter-thres",
                                            options=[
                                                {"label": "0.01", "value": 0.01},
                                                {"label": "0.05", "value": 0.05},
                                            ],
                                            value=0.05,
                                        )
                                    ], size="lg"),
                                ]
                            ),
                        )
                    ], style={'height': '25vh'}, color="success", outline=True
                )),
                dbc.Col(dbc.Card(
                    [
                        dbc.CardHeader(html.P("Smooth spectrum", style={"font-weight": "bold"})),
                        dbc.CardBody(
                            html.Div(
                                [
                                    html.P("Thresholds for smoothing", style={"font-style": 'italic'}),
                                    dbc.InputGroup([
                                        dbc.InputGroupText("m = ", style={"font-weight": "bold"}),
                                        dbc.Input(
                                            id="page-3-upload-conti-smooth-thres-m",
                                            type="number",
                                            value=1,
                                        )
                                    ], size="lg"),
                                    html.Br(),
                                    dbc.InputGroup([
                                        dbc.InputGroupText("n = ", style={"font-weight": "bold"}),
                                        dbc.Input(
                                            id="page-3-upload-conti-smooth-thres-n",
                                            type="number",
                                            value=1,
                                        )
                                    ], size="lg"),
                                ]
                            ),
                        )
                    ], style={'height': '25vh'}, color="danger", outline=True
                )),
                dbc.Col(dbc.Card(
                    [
                        dbc.CardHeader(html.P("Add noise", style={"font-weight": "bold"})),
                        dbc.CardBody(
                            html.Div(
                                [
                                    html.P("SNR for adding noise", style={"font-style": 'italic'}),
                                    dbc.RadioItems(
                                        id="page-3-upload-conti-snr-noise",
                                        options=[
                                            {"label": "1000", "value": 1000},
                                            {"label": "5000", "value": 5000},
                                            {"label": "10000", "value": 10000},
                                        ],
                                        value=1000,
                                        inline=True,
                                    ),
                                ]
                            ),
                        )
                    ], style={'height': '25vh'}, color="warning", outline=True
                )),
            ]
        ),
        html.Br(),
        html.Div(
            [
                dbc.Button(id='page-3-upload-conti-confirm-param-btn', children='Simulate spectra', n_clicks=0,
                           style={"font-weight": "bold"}, size='lg', color='primary')
            ],
            className="d-grid gap-2 d-md-flex justify-content-md-center",
        ),
        dcc.Store(id='page-3-upload-conti-processed-spectra-data-dict'),
        html.Br(),
    ])]
)

page_3_upload_conti_no_peak_shift_spectra_div = html.Div(
    id="page-3-upload-conti-no-peak-shift-spectra-div",
    children=[
        dbc.Card(
            id='page-3-upload-conti-no-peak-shift-all-spectra-div',
            children=[
                dbc.CardBody(
                    [
                        html.H6("Metabolite spectrum in the selected mixture",
                                style={"font-weight": "bold", "color": "steelblue"}),
                        html.Br(),
                        dcc.Loading(children=[
                            dcc.Dropdown(id='page-3-upload-conti-no-peak-shift-select-mix', multi=False,),
                            dcc.Graph(id='page-3-upload-conti-no-peak-shift-all-meta-fig')
                        ]),
                    ]
                )
            ],
        ),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    dbc.Card(
                        id="page-3-upload-conti-no-peak-shift-meta-y-plot-div",
                        children=[
                            dcc.Loading(dcc.Graph(id='page-3-upload-conti-no-peak-shift-meta-y-fig')),
                            html.Br(),
                        ],
                        style={'height': '72vh'}
                    )
                ),
                dbc.Col(
                    dbc.Card(
                        id="page-3-upload-conti-no-peak-shift-mix-spectra-div",
                        children=[
                            dcc.Loading(children=[dcc.Store(id="page-3-upload-conti-no-peak-shift-final-replicate-dict"),
                                                dcc.Dropdown(
                                                    id="page-3-upload-conti-no-peak-shift-select-replicate",
                                                    multi=False
                                                ),
                                                dcc.Graph(id="page-3-upload-conti-no-peak-shift-mix-spectra-fig"),
                                                # html.P("Contour Level Slider:"),
                                                # dcc.RangeSlider(id='page-3-upload-conti-no-peak-shift-contour-level-1', allowCross=False,
                                                #                 min=0, max=100),
                                                # html.Div(id='page-3-upload-conti-no-peak-shift-update-level-container-1', style={'margin-top': 10}),
                                                html.Br(),]),
                        ],
                        style={'height': '72vh'}
                    )
                ),
            ]
        ),
        html.Br(),
        html.Div(children=[html.Button("Download simulated spectra as csv", id="page-3-upload-conti-no-peak-shift-btn-csv")],
                 className="d-grid gap-2 d-md-flex justify-content-md-left"),
        dcc.Download(id="page-3-upload-conti-no-peak-shift-download-spectra-csv"),
    ]
)

page_3_upload_conti_no_peak_shift_results_div = html.Div(
    id="page-3-upload-conti-no-peak-shift-results-div",
    children=[
        page_3_upload_conti_advance_parameters_setting_card,
        html.Br(),
        page_3_upload_conti_no_peak_shift_spectra_div,
        html.Br(),
    ],
    style={"display": "none"},
)

# ---------- tab 3: upload continuous with peak shift results div ----------
page_3_upload_conti_pH_same_div = html.Div(
    id="page-3-upload-conti-ph-same-div",
    children=[
        html.P("Please input the pH:", style={"font-style": 'italic'}),
        dbc.InputGroup([
            dbc.InputGroupText("pH = ", style={"font-weight": "bold"}),
            dbc.Input(
                id="page-3-upload-conti-same-ph",
                type="number",
                value=6.5,
            )
        ], size="lg")
    ],
    style={"display": "none"}
)

page_3_upload_conti_pH_not_same_div = html.Div(
    id="page-3-upload-conti-ph-not-same-div",
    children=[
        html.P("Please input the mean and std for pH:", style={"font-style": 'italic'}),
        dbc.Row(
            [
                dbc.Col(dbc.InputGroup([
                    dbc.InputGroupText("mean = ", style={"font-weight": "bold"}),
                    dbc.Input(
                        id="page-3-upload-conti-ph-mean",
                        type="number",
                        value=7.4,
                    )
                ], size="lg")),
                dbc.Col(dbc.InputGroup([
                    dbc.InputGroupText("std = ", style={"font-weight": "bold"}),
                    dbc.Input(
                        id="page-3-upload-conti-ph-std",
                        type="number",
                        value=1.5
                    )
                ], size="lg")),
            ]
        ),
    ],
    style={"display": "none"}
)

page_3_upload_conti_set_pH_card = dbc.Card(
    [
    dbc.CardBody(
        [
            html.P("Set the same pH for all replicates: ", style={"font-weight": "bold"}),
            dbc.RadioItems(
                id="page-3-upload-conti-ph-same-flag",
                options=[
                    {"label": "True", "value": "true"},
                    {"label": "False", "value": "false"},
                ],
                value="false",
                inline=True,
            ),
            html.Br(),
            page_3_upload_conti_pH_same_div,
            page_3_upload_conti_pH_not_same_div,
            html.Br(),
        ]
    )], className="w-75 mb-3",
)

page_3_upload_conti_cons_and_ph_table_div = html.Div(
    id="page-3-upload-conti-cons-ph-table-div",
    children=[
        dash_table.DataTable(id="page-3-upload-conti-cons-ph-table", export_format="csv", export_headers="display",
                             editable=True, virtualization=True,
                             style_cell={
                               'font_family': 'arial',
                               'font_size': '16px',
                               'text_align': 'center',
                               'minWidth': '180px', 'width': '180px', 'maxWidth': '180px',
                             },
                             style_header={
                                'backgroundColor': 'whitesmoke',
                                'fontWeight': 'bold'
                             },
                             style_data_conditional=[
                               {
                                   "if": {"column_id": "meta_name"},
                                   "fontWeight": "bold",
                               },
                             ],
                             fixed_rows={'headers': True},
                             style_table={'height': 250},
                             css=[
                               {"selector": ".dash-spreadsheet tr th", "rule": "height: 35px;"},
                               {"selector": ".dash-spreadsheet tr td", "rule": "height: 35px;"},
                               {"selector": ".export", "rule": "position: absolute; right:0%; top:-40px; "
                                                               "background-color: lightsteelblue; "
                                                               "font-weight: 300; font-size: 1.5rem;"},
                             ]),
        html.Br(),
        html.Br(),
    ],
    style={"display": "none"}
)

page_3_upload_conti_peak_shift_part_2_div = html.Div(
    id="page-3-upload-conti-peak-shift-part-2-div",
    children=[
        html.H6("Part II: simulate spectra with peak shift", style={"font-weight": "bold"}),
        dbc.Card(
            id="page-3-upload-conti-peak-shift-ad-param-setting-card",
            children=[dbc.CardBody([
                html.H6("Advanced Parameters Setting", style={"font-weight": "bold", "color": "steelblue"}),
                html.Br(),
                dbc.Row(
                    [
                        dbc.Col(dbc.Card(
                            [
                                dbc.CardHeader(html.P("Remove exclusion regions", style={"font-weight": "bold"})),
                                dbc.CardBody(
                                    html.Div(
                                        [
                                            html.P("Range of calibration peaks (ppm)", style={"font-style": 'italic'}),
                                            dcc.RangeSlider(
                                                id="page-3-upload-conti-peak-shift-calibration-range",
                                                min=-0.3, max=0.3, step=0.1,
                                                value=[-0.3, 0.3],
                                            ),
                                            html.Br(),
                                            html.P("Range of water suppression (ppm)", style={"font-style": 'italic'}),
                                            dcc.RangeSlider(
                                                id="page-3-upload-conti-peak-shift-water-range",
                                                min=4.5, max=5.5, step=0.1,
                                                value=[4.5, 5.0]
                                            )
                                        ]
                                    ),
                                )
                            ], style={'height': '25vh'}, color="primary", outline=True
                        )),
                        dbc.Col(dbc.Card(
                            [
                                dbc.CardHeader(html.P("Filter noise", style={"font-weight": "bold"})),
                                dbc.CardBody(
                                    html.Div(
                                        [
                                            html.P("Parameters for filtering noise", style={"font-style": 'italic'}),
                                            dbc.InputGroup([
                                                dbc.InputGroupText("threshold = ", style={"font-weight": "bold"}),
                                                dbc.Select(
                                                    id="page-3-upload-conti-peak-shift-filter-thres",
                                                    options=[
                                                        {"label": "0.01", "value": 0.01},
                                                        {"label": "0.05", "value": 0.05},
                                                    ],
                                                    value=0.05,
                                                )
                                            ], size="lg"),
                                        ]
                                    ),
                                )
                            ], style={'height': '25vh'}, color="success", outline=True
                        )),
                        dbc.Col(dbc.Card(
                            [
                                dbc.CardHeader(html.P("Smooth spectrum", style={"font-weight": "bold"})),
                                dbc.CardBody(
                                    html.Div(
                                        [
                                            html.P("Thresholds for smoothing", style={"font-style": 'italic'}),
                                            dbc.InputGroup([
                                                dbc.InputGroupText("m = ", style={"font-weight": "bold"}),
                                                dbc.Input(
                                                    id="page-3-upload-conti-peak-shift-smooth-thres-m",
                                                    type="number",
                                                    value=1,
                                                )
                                            ], size="lg"),
                                            html.Br(),
                                            dbc.InputGroup([
                                                dbc.InputGroupText("n = ", style={"font-weight": "bold"}),
                                                dbc.Input(
                                                    id="page-3-upload-conti-peak-shift-smooth-thres-n",
                                                    type="number",
                                                    value=1,
                                                )
                                            ], size="lg"),
                                        ]
                                    ),
                                )
                            ], style={'height': '25vh'}, color="danger", outline=True
                        )),
                        dbc.Col(dbc.Card(
                            [
                                dbc.CardHeader(html.P("Add noise", style={"font-weight": "bold"})),
                                dbc.CardBody(
                                    html.Div(
                                        [
                                            html.P("SNR for adding noise", style={"font-style": 'italic'}),
                                            dbc.RadioItems(
                                                id="page-3-upload-conti-peak-shift-snr-noise",
                                                options=[
                                                    {"label": "1000", "value": 1000},
                                                    {"label": "5000", "value": 5000},
                                                    {"label": "10000", "value": 10000},
                                                ],
                                                value=1000,
                                                inline=True,
                                            ),
                                        ]
                                    ),
                                )
                            ], style={'height': '25vh'}, color="warning", outline=True
                        )),
                    ]
                ),
                html.Br(),
                html.Div(
                    [
                        dbc.Button(id='page-3-upload-conti-peak-shift-confirm-param-btn', children='Simulate spectra', n_clicks=0,
                                   style={"font-weight": "bold"}, size='lg', color='primary')
                    ],
                    className="d-grid gap-2 d-md-flex justify-content-md-center",
                ),
                dcc.Store(id='page-3-upload-conti-peak-shift-processed-spectra-data-dict'),
                html.Br(),
            ])]
        )
    ],
    style={"display": "none"}
)

page_3_upload_conti_peak_shift_spectra_div = html.Div(
    id="page-3-upload-conti-peak-shift-spectra-div",
    children=[
        dbc.Card(
            id="page-3-upload-conti-peak-shift-all-spectra-card",
            children=[
                dbc.CardBody(
                    [
                        html.H6("Metabolite spectrum in the selected mixture",
                                style={"font-weight": "bold", "color": "steelblue"}),
                        html.Br(),
                        dcc.Loading(children=[
                            dcc.Dropdown(id='page-3-upload-conti-peak-shift-select-mix', multi=False,),
                            dcc.Graph(id='page-3-upload-conti-peak-shift-all-meta-fig')
                        ], type="default"),
                    ]
                )
            ]
        ),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(
                    # dbc.Card(
                    #     id="page-3-upload-conti-peak-shift-stacked-fig-card",
                    #     children=[
                    #         dbc.CardBody(
                    #             [
                    #                 dcc.Loading(children=[dcc.Graph(id='page-3-upload-conti-peak-shift-stacked-fig'),
                    #                                     html.P("Vertical Spacing Slider:"),
                    #                                     dcc.Slider(
                    #                                         id='page-3-upload-conti-peak-shift-vs-slider',
                    #                                         min=0,
                    #                                         max=3,
                    #                                         step=0.01,
                    #                                         marks={i: '{}'.format(10 ** i) for i in range(4)},
                    #                                         value=1,
                    #                                     ),
                    #                                     html.Div(id='page-3-upload-conti-peak-shift-update-vs-container', style={'margin-top': 10})],
                    #                             type="default"),
                    #             ]
                    #         )
                    #     ],
                    #     style={'height': '72vh'}
                    # ),
                ),
                dbc.Col(
                    dbc.Card(
                        id="page-3-upload-conti-peak-shift-mix-spectra-div",
                        children=[
                            dcc.Loading(children=[dcc.Store(id="page-3-upload-conti-peak-shift-final-replicate-dict"),
                                                dcc.Dropdown(
                                                    id="page-3-upload-conti-peak-shift-select-replicate",
                                                    multi=False
                                                ),
                                                dcc.Graph(id="page-3-upload-conti-peak-shift-mix-spectra-fig"),
                                                # html.P("Contour Level Slider:"),
                                                # dcc.RangeSlider(id='page-3-upload-conti-peak-shift-contour-level-1',
                                                #               allowCross=False,
                                                #               min=0, max=100),
                                                # html.Div(id='page-3-upload-conti-peak-shift-update-level-container-1',
                                                #             style={'margin-top': 10}),
                                                html.Br(),
                                                ],
                                        type="default"),
                        ],
                        style={'height': '72vh'}
                    )
                ),
            ]
        ),
        html.Br(),
        dbc.Card(
            id="page-3-upload-conti-peak-shift-meta-y-plot-div",
            children=[
                dcc.Loading(dcc.Graph(id='page-3-upload-conti-peak-shift-meta-y-fig'), type="default"),
                # dcc.Graph(id='conti-peak-shift-meta-y-fig'),
            ],
            # style={'height': '55vh'}
        ),
        html.Br(),
        html.Button("Download simulated spectra as csv", id="page-3-upload-conti-peak-shift-btn-csv"),
        dcc.Download(id="page-3-upload-conti-peak-shift-download-spectra-csv"),
    ],
    style={"display": "none"}
)

page_3_upload_conti_peak_shift_results_div = html.Div(
    id="page-3-upload-conti-peak-shift-results-div",
    children=[
        html.H6("Part I: define the pH for all replicates", style={"font-weight": "bold"}),
        page_3_upload_conti_set_pH_card,
        html.Br(),
        html.Div([
                    dbc.Button(id='page-3-tab-3-upload-conti-submit-part-1-btn', children='OK', n_clicks=0, style={"font-weight": "bold"},
                               size='lg', color='secondary'),
        ]),
        html.Br(),
        html.Br(),
        page_3_upload_conti_cons_and_ph_table_div,
        page_3_upload_conti_peak_shift_part_2_div,
        html.Br(),
        page_3_upload_conti_peak_shift_spectra_div,
        html.Br(),
    ],
    style={"display": "none"},
)

# tab 3
page_3_tab_3 = dcc.Tab(
    id="page-3-tab-3",
    label="STEP 3.   Simulation Results",
    value="page-3-tab-3",
    style=page_3_tab_style, selected_style=page_3_tab_selected_style,
    children=[
        html.Div(
            [
                html.Br(),
                html.H5("Simulate Spectra", style={"font-weight": "bold", 'color': 'steelblue'}),
                html.Br(),
                page_3_select_peak_shift_div,
                page_3_group_no_peak_shift_results_div,
                page_3_group_peak_shift_results_div,
                page_3_conti_no_peak_shift_results_div,
                page_3_conti_peak_shift_results_div,
                page_3_upload_group_no_peak_shift_results_div,
                page_3_upload_group_peak_shift_results_div,
                page_3_upload_conti_no_peak_shift_results_div,
                page_3_upload_conti_peak_shift_results_div,
                html.Br(),
            ],
            className="container__1",
        )
    ],
)


layout = html.Div([
    html.H3("2D COSY Simulator", style={"font-style": 'italic', "font-weight": "bold", "color": "dimgrey"}),
    dcc.Tabs(
        id="page-3-whole-tabs",
        children=[page_3_tab_1, page_3_tab_2, page_3_tab_3, ],
        value="page-3-tab-1",
    )
])

# --------------------------------------------- start callbacks ---------------------------------------------
# ----------------------------- tab 1 callbacks -----------------------------
@app.callback(
    [
        Output("page-3-metabolite-select", "value"),
        Output("page-3-store-upload-file", "data"),
        # Output("store-upload-file", "clear_data")
    ],
    [
        Input("page-3-upload-mixture", "filename"),
        Input("page-3-upload-mixture", "contents"),
    ]
)
def page_3_upload_input_mixture_test(input_filename, input_contents):
    if input_contents is not None and input_filename is not None:
        content_type, content_string = input_contents.split(',')
        decoded = base64.b64decode(content_string)
        input_mixtures = pd.read_csv(io.StringIO(decoded.decode("utf-8")), sep="\n", header=None).iloc[:, 0].values.tolist()
        format_mixture_list = format_input_mixture(input_mixtures)
        in_db_mixture_list = input_match_db(input_mixtures, match_data_dict, hmdb_dict)
        return in_db_mixture_list, format_mixture_list
    else:
        return ['citric acid', 'd-arabitol', 'd-lactate', 'l-alanine', 'l-serine'], None


@app.callback(
    [
        Output("page-3-hide-select-metabolites", "children"),
        Output("page-3-hide-select-metabolites", "style"),
        Output("page-3-store-upload-file", "clear_data")
    ],
    [Input("page-3-submit-meta-btn", 'n_clicks')],
    [
        State("page-3-metabolite-select", "value"),
        State("page-3-store-upload-file", "data")
    ]
)
def page_3_show_selected_metabolites_test(n_clicks, in_db_mixture_list, all_mixture_list):
    if n_clicks == 0:
        default_mixture_list = ['citric acid', 'd-arabitol', 'd-lactate', 'l-alanine', 'l-serine']
        return "no selected metabolites", {"display": "none"}, True
    else:
        if all_mixture_list:
            not_in_db_list = [x for x in all_mixture_list if x not in in_db_mixture_list]
            if not_in_db_list:
                exist_string = ', '.join(in_db_mixture_list)
                not_exist_string = ', '.join(not_in_db_list)
                print_text = [html.P([html.Span(f'These metabolites are in the database: '),
                                     html.Span(f'{exist_string}', style={"font-weight": "bold", 'color': 'green '})]),
                              html.P([html.Span(f'These metabolites are not in the database: '),
                                      html.Span(f'{not_exist_string}', style={"font-weight": "bold", 'color': 'red '})])]
                return print_text, {"display": "block"}, True
            else:
                exist_string = ', '.join(in_db_mixture_list)
                print_text = html.P([html.Span(f'These metabolites are all in database: '),
                                    html.Span(f'{exist_string}', style={"font-weight": "bold", 'color': 'green '})])
                return print_text, {"display": "block"}, True

        else:
            text_string = ', '.join(in_db_mixture_list)
            print_text = html.P([html.Span(f'Selected metabolites are: '),
                                 html.Span(f'{text_string}', style={"font-weight": "bold", 'color': 'green '})])
            return print_text, {"display": "block"}, True


@app.callback(
    [
        Output("page-3-hide-multi-hmdb-names", "children"),
        Output("page-3-hide-multi-hmdb-names", "style"),
        Output("page-3-hide-confirm-select", "children"),
        Output("page-3-hide-confirm-select", "style"),
    ],
    [Input("page-3-submit-meta-btn", 'n_clicks')],
    [
        State("page-3-metabolite-select", "value"),
    ]
)
def page_3_update_hmdb_selections(n_clicks, mixture_list):
    if n_clicks == 0:
        return [], {"display": "none"}, [], {"display": "none"}
    else:
        db_name_hmdb_id_dict = page_3_get_name_hmdb_id_dict(mixture_list)
        multi_hmdb_id_select = []
        confirm_div = []
        for db_name, idx_dict_list in db_name_hmdb_id_dict.items():
            if len(idx_dict_list) > 1:
                multi_hmdb_id_select.append(db_name)
                options = []
                for idx_name_dict in idx_dict_list:
                    hmdb_id = list(idx_name_dict.keys())[0]
                    hmdb_name = list(idx_name_dict.values())[0]
                    options.append({"label": str(hmdb_id)+', '+str(hmdb_name),
                                    "value": str(db_name)+', '+str(hmdb_id)+', '+str(hmdb_name)})
                confirmitems = html.Div(
                    [
                        html.P([html.Span("Please confirm the HMDB ID for "),
                                html.Span(f"{db_name}", style={"font-weight": "bold", 'color': 'steelblue '})]),
                        dbc.RadioItems(id={"type": "radioitems", "index": db_name},
                                       options=options,
                                       value=db_name+', '+list(idx_dict_list[0].keys())[0]+", "+list(idx_dict_list[0].values())[0]),
                        html.Br(),
                    ]
                )
                confirm_div.append(confirmitems)

        if multi_hmdb_id_select:
            multi_text_string = ', '.join(multi_hmdb_id_select)
            print_text = html.P([html.Span(f'These metabolite names correspond to more than one HMDB IDs: '),
                                     html.Span(f'{multi_text_string}', style={"font-weight": "bold", 'color': 'steelblue '})])
            return print_text, {"display": "block"}, confirm_div, {"display": "block"}
        else:
            return [], {"display": "none"}, [], {"display": "none"}


def page_3_get_name_hmdb_id_dict(in_db_mixture_list):
    db_name_hmdb_id_dict = dict()
    for name in in_db_mixture_list:
        idx_list = []
        for idx, name_list in hmdb_dict.items():
            if name in name_list:
                idx_list.append({idx: hmdb_dict[idx][0]})
        db_name_hmdb_id_dict[name] = idx_list
    return db_name_hmdb_id_dict


@app.callback(
    [
        Output("page-3-hide-mixture-table", "children"),
        Output("page-3-hide-mixture-table", "style"),
        Output("page-3-hide-confirm-btn", "className"),
        Output("page-3-hide-confirm-btn", "style"),
        Output("page-3-db-names-hmdb-ids-dict", "data")
    ],
    [
        Input("page-3-submit-meta-btn", 'n_clicks'),
        Input({"type": "radioitems", "index": ALL}, "value"),
    ],
    [
        State("page-3-metabolite-select", "value"),
    ]
)
def page_3_update_names_hmdb_table(n_clicks, selected_list, mixture_list):
    if n_clicks == 0:
        return [], {"display": "none"}, None, {"display": "none"}, None
    else:
        db_name_hmdb_id_dict = page_3_get_name_hmdb_id_dict(mixture_list)
        table_header = [html.Thead(html.Tr([html.Th("Name"), html.Th("HMDB_ID"), html.Th("HMDB_Name")]))]
        rows_list = []
        final_name_hmdb_id_dict = dict()
        for db_name, idx_dict_list in db_name_hmdb_id_dict.items():
            if len(idx_dict_list) == 1:
                hmdb_id = list(idx_dict_list[0].keys())[0]
                hmdb_name = list(idx_dict_list[0].values())[0]
                row = html.Tr([html.Td(db_name), html.Td(hmdb_id), html.Td(hmdb_name)])
                rows_list.append(row)
                final_name_hmdb_id_dict[db_name] = hmdb_id

        for selected_str in selected_list:
            db_name, hmdb_id, hmdb_name = re.findall(r'(.+), (HMDB\d+), (.+)', selected_str)[0]
            row = html.Tr([html.Td(db_name), html.Td(hmdb_id), html.Td(hmdb_name)],
                          style={"font-weight": "bold", "color": "darkred"})
            rows_list.append(row)
            final_name_hmdb_id_dict[db_name] = hmdb_id

        # print(final_name_hmdb_id_dict)
        table_body = [html.Tbody(rows_list)]
        table = dbc.Table(table_header + table_body, bordered=True, color='primary')
        return table, {"display": "block"}, \
            "d-grid gap-2 d-md-flex justify-content-md-center", {"display": "block"}, \
            final_name_hmdb_id_dict


# ----------------------------- tab 2 callbacks -----------------------------
@app.callback(
    Output("page-3-select-upload-cons-div", "style"),
    Input("page-3-confirm-meta-btn", "n_clicks"),
)
def page_3_update_tab_2_div(n_clicks):
    if n_clicks == 0:
        return {"display": "none"}
    else:
        return {"display": "block"}


@app.callback(
    [
        Output("page-3-select-cons-type-div", "style"),
        Output("page-3-upload-final-cons-div", "style"),
    ],
    Input("page-3-select-upload-cons-div", "style"),
    Input("page-3-select-upload-cons", "value")
)
def page_3_select_upload_own_cons(upload_div, upload_or_not):
    if upload_div == {"display": "none"}:
        return {"display": "none"}, {"display": "none"}
    else:
        if upload_or_not == "not upload":
            return {"display": "block"}, {"display": "none"}
        elif upload_or_not == "upload":
            return {"display": "none"}, {"display": "block"}
        else:
            return {"display": "none"}, {"display": "none"}


@app.callback(
    [
        Output("page-3-upload-discrete-cons-data-div", "style"),
        Output("page-3-upload-conti-cons-data-div", "style")
    ],
    Input("page-3-upload-final-cons-div", "style"),
    Input("page-3-select-upload-cons-type", "value")
)
def page_3_select_the_type_of_upload_cons(upload_div, cons_type):
    if upload_div == {"display": "none"}:
        return {"display": "none"}, {"display": "none"}
    else:
        if cons_type == "discrete":
            return {"display": "block"}, {"display": "none"}
        elif cons_type == "continuous":
            return {"display": "none"}, {"display": "block"}
        else:
            return {"display": "none"}, {"display": "none"}


@ app.callback(
    [
        Output("page-3-group-sim-div", "style"),
        Output("page-3-conti-sim-div", "style"),
    ],
    Input("page-3-select-cons-type-div", "style"),
    Input("page-3-select-simulation-type", "value"),
)
def page_3_select_simulation_type(select_div, sim_type):
    if select_div == {"display": "none"}:
        return {"display": "none"}, {"display": "none"}
    else:
        if sim_type == "group":
            return {"display": "block"}, {"display": "none"}
        elif sim_type == "continuous":
            return {"display": "none"}, {"display": "block"}
        else:
            return {"display": "none"}, {"display": "none"}


# -------------- tab 2 callbacks: upload group cons table --------------
@app.callback(
    [
        Output("page-3-upload-discrete-cons-table", "data"),
        Output("page-3-upload-discrete-cons-table", "columns"),
        Output("page-3-confirm-upload-discrete-cons-btn-div", "style"),
        Output("page-3-confirm-upload-discrete-cons-btn-div", "className")
    ],
    Input("page-3-upload-final-discrete-cons-file", "filename"),
    Input("page-3-upload-final-discrete-cons-file", "contents")
)
def page_3_show_upload_discrete_cons_table(input_filename, input_contents):
    if input_contents is not None and input_filename is not None:
        content_type, content_string = input_contents.split(',')
        decoded = base64.b64decode(content_string)
        df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), header=[0, 1])

        table_columns = [{"name": ["", "Name"], "id": "meta_name"}] + \
            [{"name": list(col), "id": re.findall(r'Group (\d+)+', col[0])[0]+"_"+col[1].lower()} for col in df.columns[1:]]

        id_list = [dic['id'] for dic in table_columns]
        table_data = []
        for d in df.to_dict('records'):
            table_data.append(dict(zip(id_list, d.values())))
        return table_data, table_columns, {"display": "block"}, "d-grid gap-2 d-md-flex justify-content-md-center"
    else:
        return [{}], [], {"display": "none"}, None


# -------------- tab 2 callbacks: upload continuous cons table --------------
@app.callback(
    [
        Output("page-3-upload-conti-cons-table", "data"),
        Output("page-3-upload-conti-cons-table", "columns"),
        Output("page-3-confirm-upload-conti-cons-btn-div", "style"),
        Output("page-3-confirm-upload-conti-cons-btn-div", "className")
    ],
    Input("page-3-upload-final-conti-cons-file", "filename"),
    Input("page-3-upload-final-conti-cons-file", "contents")
)
def page_3_show_upload_continuous_cons_table(input_filename, input_contents):
    if input_contents is not None and input_filename is not None:
        content_type, content_string = input_contents.split(',')
        decoded = base64.b64decode(content_string)
        df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))

        table_columns = [{"name": "Name", "id": "meta_name"}] + \
                        [{"name": col, "id": col.lower()} for col in df.columns[1:]]

        id_list = [dic['id'] for dic in table_columns]
        table_data = []
        for d in df.to_dict('records'):
            table_data.append(dict(zip(id_list, d.values())))
        return table_data, table_columns, {"display": "block"}, "d-grid gap-2 d-md-flex justify-content-md-center"
    else:
        return [{}], [], {"display": "none"}, None


# -------------- tab 2 callbacks: simulate group cons table --------------
@app.callback(
    Output("page-3-group-part-1-div", "style"),
    Input("page-3-group-sim-div", "style")
)
def page_3_show_group_part_1_div(sim_type_style):
    if sim_type_style == {"display": "block"}:
        return {"display": "block"}
    else:
        return {"display": "none"}


@app.callback(
    [
        Output("page-3-group-1-hmdb-cons-data-select", "style"),
        Output("page-3-group-1-not-hmdb-cons-data-upload", "style"),
    ],
    Input("page-3-part-1-group-1-hmdb-data", "value")
)
def page_3_update_hmdb_data_select_1(select_hmdb_data):
    if select_hmdb_data == "hmdb":
        return {"display": "block"}, {"display": "none"}
    elif select_hmdb_data == "not hmdb":
        return {"display": "none"}, {"display": "block"}
    else:
        return {"display": "none"}, {"display": "none"}


@app.callback(
    [
        Output("page-3-group-2-hmdb-cons-data-select", "style"),
        Output("page-3-group-2-not-hmdb-cons-data-upload", "style"),
    ],
    Input("page-3-part-1-group-2-hmdb-data", "value")
)
def page_3_update_hmdb_data_select_2(select_hmdb_data):
    if select_hmdb_data == "hmdb":
        return {"display": "block"}, {"display": "none"}
    elif select_hmdb_data == "not hmdb":
        return {"display": "none"}, {"display": "block"}
    else:
        return {"display": "none"}, {"display": "none"}


@app.callback(
    [
        Output("page-3-hide-cons-mean-std-table-div", "style"),
        Output("page-3-cons-mean-std-table", "data"),
        Output("page-3-cons-mean-std-table", "columns"),
        Output("page-3-group-part-2-div", "style")
    ],
    Input("page-3-submit-part-1-btn", "n_clicks"),
    [
        State("page-3-db-names-hmdb-ids-dict", "data"),
        State("page-3-part-1-group-1-hmdb-data", "value"),
        State("page-3-group-1-bio-type", "value"),
        State("page-3-group-1-upload-mean-std", "filename"),
        State("page-3-group-1-upload-mean-std", "contents"),
        State("page-3-part-1-group-2-hmdb-data", "value"),
        State("page-3-group-2-bio-type", "value"),
        State("page-3-group-2-upload-mean-std", "filename"),
        State("page-3-group-2-upload-mean-std", "contents"),
    ]
)
def page_3_update_cons_mean_std_table(n_clicks, final_name_hmdb_id_dict,
                               select_hmdb_data_1, bio_type_1, upload_cons_filename_1, upload_cons_content_1,
                               select_hmdb_data_2, bio_type_2, upload_cons_filename_2, upload_cons_content_2):
    if n_clicks == 0:
        return {"display": "none"}, [{}], [], {"display": "none"}
    else:
        table_columns = [
            {"name": ["", "Name"], "id": "meta_name"},
            {"name": ["", "HMDB ID"], "id": "hmdb_id"},
            {"name": ["Group 1 (Normal)", "Mean (uM)"], "id": "mean_1"},
            {"name": ["Group 1 (Normal)", "Std (uM)"], "id": "std_1"},
            {"name": ["Group 2 (Abnormal)", "Mean (uM)"], "id": "mean_2"},
            {"name": ["Group 2 (Abnormal)", "Std (uM)"], "id": "std_2"},
        ]
        if select_hmdb_data_1 == "hmdb" and select_hmdb_data_2 == "hmdb":
            table_data = []
            for meta_name, hmdb_id in final_name_hmdb_id_dict.items():
                avg_mean_1, avg_std_1, unit_1 = get_hmdb_normal_avg_cons(hmdb_norm_cons_dict, hmdb_id, bio_type_1)
                avg_mean_2, avg_std_2, unit_2 = get_hmdb_abnormal_avg_cons(hmdb_abnorm_cons_dict, hmdb_id, bio_type_2)
                temp_dict = {"meta_name": meta_name, "hmdb_id": hmdb_id, "mean_1": avg_mean_1, "std_1": avg_std_1,
                             "mean_2": avg_mean_2, "std_2": avg_std_2}
                table_data.append(temp_dict)
            return {"display": "block"}, table_data, table_columns, {"display": "block"}

        elif select_hmdb_data_1 == "hmdb" and select_hmdb_data_2 == "not hmdb":
            table_data = []
            temp_cons_df = page_3_parse_cons_contents(upload_cons_content_2, upload_cons_filename_2)
            temp_cons_df.columns = ["meta_name", "mean_2", "std_2"]
            temp_cons_dict = temp_cons_df.to_dict('records')
            for meta_name, hmdb_id in final_name_hmdb_id_dict.items():
                avg_mean_1, avg_std_1, unit_1 = get_hmdb_normal_avg_cons(hmdb_norm_cons_dict, hmdb_id, bio_type_1)
                temp_dict = list(filter(lambda d: d['meta_name'] == meta_name, temp_cons_dict))[0]
                temp_dict["hmdb_id"] = hmdb_id
                temp_dict["mean_1"] = avg_mean_1
                temp_dict["std_1"] = avg_std_1
                table_data.append(temp_dict)
            return {"display": "block"}, table_data, table_columns, {"display": "block"}

        elif select_hmdb_data_1 == "not hmdb" and select_hmdb_data_2 == "hmdb":
            table_data = []
            temp_cons_df = page_3_parse_cons_contents(upload_cons_content_1, upload_cons_filename_1)
            temp_cons_df.columns = ["meta_name", "mean_1", "std_1"]
            temp_cons_dict = temp_cons_df.to_dict('records')
            for meta_name, hmdb_id in final_name_hmdb_id_dict.items():
                avg_mean_2, avg_std_2, unit_2 = get_hmdb_abnormal_avg_cons(hmdb_abnorm_cons_dict, hmdb_id, bio_type_2)
                temp_dict = list(filter(lambda d: d['meta_name'] == meta_name, temp_cons_dict))[0]
                temp_dict["hmdb_id"] = hmdb_id
                temp_dict["mean_2"] = avg_mean_2
                temp_dict["std_2"] = avg_std_2
                table_data.append(temp_dict)
            return {"display": "block"}, table_data, table_columns, {"display": "block"}

        elif select_hmdb_data_1 == "not hmdb" and select_hmdb_data_2 == "not hmdb":
            temp_cons_df_1 = page_3_parse_cons_contents(upload_cons_content_1, upload_cons_filename_1)
            temp_cons_df_2 = page_3_parse_cons_contents(upload_cons_content_2, upload_cons_filename_2)
            temp_cons_df_1.columns = ["meta_name", "mean_1", "std_1"]
            temp_cons_df_2.columns = ["meta_name", "mean_2", "std_2"]
            table_data = pd.merge(temp_cons_df_1, temp_cons_df_2, on="meta_name").to_dict("records")
            final_table_data = []
            for meta_name, hmdb_id in final_name_hmdb_id_dict.items():
                temp_dict = list(filter(lambda d: d['meta_name'] == meta_name, table_data))[0]
                temp_dict["hmdb_id"] = hmdb_id
                final_table_data.append(temp_dict)
            return {"display": "block"}, final_table_data, table_columns, {"display": "block"}


def page_3_parse_cons_contents(contents, filename):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    if 'csv' in filename:
        # Assume that the user uploaded a CSV file
        df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
        format_df = input_cons_match_db(df, match_data_dict, hmdb_dict)
        return format_df
    elif 'xls' in filename:
        # Assume that the user uploaded an excel file
        df = pd.read_excel(io.BytesIO(decoded))
        format_df = input_cons_match_db(df, match_data_dict, hmdb_dict)
        return format_df
    elif 'txt' in filename:
        df = pd.read_csv(io.StringIO(decoded.decode("utf-8")), sep=",", header=0)
        format_df = input_cons_match_db(df, match_data_dict, hmdb_dict)
        return format_df


@app.callback(
    Output("page-3-group-1-upload-corr-file-div", "style"),
    Input("page-3-group-1-corr-flag", "value")
)
def page_3_update_upload_corr_file_div_1(corr_flag):
    if corr_flag:
        return {"display": "block"}
    else:
        return {"display": "none"}


@app.callback(
    Output("page-3-group-2-upload-corr-file-div", "style"),
    Input("page-3-group-2-corr-flag", "value")
)
def page_3_update_upload_corr_file_div_2(corr_flag):
    if corr_flag:
        return {"display": "block"}
    else:
        return {"display": "none"}


@app.callback(
    [
        Output("page-3-hide-simulated-cons-table-div", "style"),
        Output("page-3-simulated-cons-table", "data"),
        Output("page-3-simulated-cons-table", "columns"),
        Output("page-3-confirm-group-cons-btn-div", "style"),
        Output("page-3-confirm-group-cons-btn-div", "className")
    ],
    Input("page-3-submit-part-2-btn", "n_clicks"),
    [
        State("page-3-cons-mean-std-table", "data"),
        State("page-3-db-names-hmdb-ids-dict", "data"),
        State("page-3-group-1-num-repli", "value"),
        State("page-3-group-1-corr-flag", "value"),
        State("page-3-group-1-corr-file", "filename"),
        State("page-3-group-1-corr-file", "contents"),
        State("page-3-group-2-num-repli", "value"),
        State("page-3-group-2-corr-flag", "value"),
        State("page-3-group-2-corr-file", "filename"),
        State("page-3-group-2-corr-file", "contents")
    ]
)
def page_3_update_simulated_cons_table(n_clicks, table_rows, final_name_hmdb_id_dict,
                                num_repli_1, corr_flag_1, corr_filename_1, corr_contents_1,
                                num_repli_2, corr_flag_2, corr_filename_2, corr_contents_2,):
    if n_clicks == 0:
        return {"display": "none"}, [{}], [], {"display": "none"}, None
    else:
        corr_df_1 = page_3_parse_corr_df(corr_filename_1, corr_contents_1)
        corr_df_2 = page_3_parse_corr_df(corr_filename_2, corr_contents_2)
        # print("table_rows: ", table_rows)
        simulated_cons_df_1 = simulate_concentrations(table_rows, num_repli_1, corr_flag_1, corr_df_1, "1")
        simulated_cons_df_2 = simulate_concentrations(table_rows, num_repli_2, corr_flag_2, corr_df_2, "2")
        table_columns = [{"name": ["", "Name"], "id": "meta_name"}, {"name": ["", "HMDB ID"], "id": "hmdb_id"}] + \
                        [{"name": ["Group 1 (Normal)", "Replicate_{}".format(i+1)], "id": "1_replicate_{}".format(i+1),
                          "hideable": True} for i in range(num_repli_1)] + \
                        [{"name": ["Group 2 (Abnormal)", "Replicate_{}".format(i+1)],
                          "id": "2_replicate_{}".format(i+1), "hideable": True} for i in range(num_repli_2)]
        simulated_cons_df = pd.merge(simulated_cons_df_1, simulated_cons_df_2, on="meta_name")
        table_data = simulated_cons_df.to_dict("records")
        final_table_data = []
        for meta_name, hmdb_id in final_name_hmdb_id_dict.items():
            temp_dict = list(filter(lambda d: d['meta_name'] == meta_name, table_data))[0]
            temp_dict["hmdb_id"] = hmdb_id
            final_table_data.append(temp_dict)
        return {"display": "block"}, table_data, table_columns, {"display": "block"}, \
            "d-grid gap-2 d-md-flex justify-content-md-center"


def page_3_parse_corr_df(input_filename, input_contents):
    if input_filename is not None and input_contents is not None:
        content_type, content_string = input_contents.split(',')
        decoded = base64.b64decode(content_string)
        corr_df = pd.read_csv(io.StringIO(decoded.decode("utf-8")), sep="\t", header=0)
        corr_df.set_index(corr_df.columns, inplace=True)
        format_corr_df = input_corr_match_db(corr_df, match_data_dict, hmdb_dict)
        return format_corr_df
    else:
        return None


@app.callback(
    [
        Output("page-3-hide-box-plot-div", "style"),
        Output("page-3-groups-cons-box-plot", "figure"),
    ],
    [
        Input("page-3-simulated-cons-table", "data"),
        Input("page-3-group-1-num-repli", "value"),
        Input("page-3-group-2-num-repli", "value"),
    ]
)
def page_3_update_cons_box_plots(table_rows, num_repli_1, num_repli_2):
    if table_rows != [{}]:
        # print(table_rows)
        fig = page_3_plot_cons_box_plot(table_rows, num_repli_1, num_repli_2)
        return {"display": "block"}, fig
    else:
        return {"display": "none"}, {"data": [], "layout": {}, "frames": []}


def page_3_plot_cons_box_plot(table_rows, num_repli_1, num_repli_2):
    fig = go.Figure()
    mixture_list = list(map(lambda d: d["meta_name"], table_rows))

    x_1 = [meta_name for meta_name in mixture_list for _ in range(num_repli_1)]
    x_2 = [meta_name for meta_name in mixture_list for _ in range(num_repli_2)]

    # plot group 1
    fig.add_trace(go.Box(
        y=[d[key] for d in table_rows for key, item in d.items() if key.startswith("1_replicate")],
        x=x_1,
        name="Group 1",
        marker_color='rgba(93, 164, 214, 0.5)',
    ))
    # plot group 2
    fig.add_trace(go.Box(
        y=[d[key] for d in table_rows for key, item in d.items() if key.startswith("2_replicate")],
        x=x_2,
        name="Group 2",
        marker_color='lightcoral',
    ))

    fig.update_layout(
        title_text="Simulated concentrations for two groups",
        template="plotly_white",
        yaxis_title='simulated concentrations (uM)',
        boxmode="group",
        # height=500,
        legend=dict(
            itemclick="toggleothers",
            itemdoubleclick="toggle"
        )
    )
    fig.update_xaxes(tickfont=dict(size=16), tickwidth=3)
    return fig


# -------------- tab 2 callbacks: simulate continuous cons table --------------
@app.callback(
    Output("page-3-conti-part-1-div", "style"),
    Input("page-3-conti-sim-div", "style")
)
def page_3_show_conti_part_1_div(sim_type_style):
    if sim_type_style == {"display": "block"}:
        return {"display": "block"}
    else:
        return {"display": "none"}


@app.callback(
    [
        Output("page-3-conti-hmdb-cons-data-select-div", "style"),
        Output("page-3-conti-upload-cons-data-div", "style"),
    ],
    Input("page-3-conti-part-1-hmdb-data", "value")
)
def page_3_show_conti_upload_mean_std_table_div(data_source):
    if data_source == "hmdb":
        return {"display": "block"}, {"display": "none"}
    elif data_source == "not hmdb":
        return {"display": "none"}, {"display": "block"}
    else:
        return {"display": "none"}, {"display": "none"}


@app.callback(
    [
        Output("page-3-conti-hide-cons-mean-std-table-div", "style"),
        Output("page-3-conti-mean-std-ab-table", "data"),
        Output("page-3-conti-mean-std-ab-table", "columns"),
        Output("page-3-conti-part-2-div", "style")
    ],
    Input("page-3-conti-submit-part-1-btn", "n_clicks"),
    [
        State("page-3-db-names-hmdb-ids-dict", "data"),
        State("page-3-conti-part-1-hmdb-data", "value"),
        State("page-3-conti-bio-type", "value"),
        State("page-3-conti-upload-mean-std-file", "filename"),
        State("page-3-conti-upload-mean-std-file", "contents"),
    ]
)
def page_3_update_continuous_cons_mean_std_table(n_clicks, final_name_hmdb_id_dict, select_hmdb_data, conti_bio_type,
                                          conti_upload_cons_filename, conti_upload_cons_content):
    if n_clicks == 0:
        return {"display": "none"}, [{}], [], {"display": "none"}
    else:
        table_columns = [
            {"name": "Name", "id": "meta_name"},
            {"name": "HMDB ID", "id": "hmdb_id"},
            {"name": "Mean (uM)", "id": "mean"},
            {"name": "Std (uM)", "id": "std"},
            {"name": "a", "id": "a"},
            {"name": "Std of random error", "id": "std_error"},
        ]
        if select_hmdb_data == "hmdb":
            table_data = []
            for meta_name, hmdb_id in final_name_hmdb_id_dict.items():
                avg_mean, avg_std, unit = get_hmdb_normal_avg_cons(hmdb_norm_cons_dict, hmdb_id, conti_bio_type)
                temp_dict = {"meta_name": meta_name, "hmdb_id": hmdb_id, "mean": avg_mean, "std": avg_std, "a": 0, "std_error": 0}
                table_data.append(temp_dict)
            return {"display": "block"}, table_data, table_columns, {"display": "block"}
        elif select_hmdb_data == "not hmdb":
            temp_cons_df = page_3_parse_cons_contents(conti_upload_cons_content, conti_upload_cons_filename)
            table_data = temp_cons_df.to_dict("records")
            final_table_data = []
            for meta_name, hmdb_id in final_name_hmdb_id_dict.items():
                temp_dict = list(filter(lambda d: d['meta_name'] == meta_name, table_data))[0]
                temp_dict["hmdb_id"] = hmdb_id
                final_table_data.append(temp_dict)
            return {"display": "block"}, final_table_data, table_columns, {"display": "block"}


@app.callback(
    [
        Output("page-3-conti-normal-div", "style"),
        Output("page-3-conti-uniform-div", "style"),
    ],
    Input("page-3-conti-part-2-select-distribution", "value")
)
def page_3_show_conti_distribution_div(dist_type):
    if dist_type == "normal":
        return {"display": "block"}, {"display": "none"}
    elif dist_type == "uniform":
        return {"display": "none"}, {"display": "block"}
    else:
        return {"display": "none"}, {"display": "none"}


@app.callback(
    [
        Output("page-3-conti-hide-simulated-cons-table-div", "style"),
        Output("page-3-conti-simulated-cons-table", "data"),
        Output("page-3-conti-simulated-cons-table", "columns"),
        Output("page-3-conti-simulated-cons-table", "style_data_conditional"),
        Output("page-3-conti-hide-plot-div", "style"),
        Output("page-3-conti-line-plot", "figure"),
        Output("page-3-conti-box-plot", "figure"),
        Output("page-3-conti-r-squared-table", "children"),
        Output("page-3-confirm-conti-cons-btn-div", "style"),
        Output("page-3-confirm-conti-cons-btn-div", "className")
    ],
    Input("page-3-conti-submit-part-2-btn", "n_clicks"),
    [
        State("page-3-conti-mean-std-ab-table", "data"),
        State("page-3-y-num-repli", "value"),
        State("page-3-conti-part-2-select-distribution", "value"),
        State("page-3-y-mean", "value"),
        State("page-3-y-std", "value"),
        State("page-3-y-min", "value"),
        State("page-3-y-max", "value"),
    ]
)
def page_3_update_conti_simulated_cons_table(n_clicks, table_rows, y_num_repli, dist_type, y_mean, y_std, y_min, y_max):
    if n_clicks == 0:
        return {"display": "none"}, [{}], [], [{}], {"display": "none"}, {"data": [], "layout": {}, "frames": []}, \
               {"data": [], "layout": {}, "frames": []}, [], {"display": "none"}, None
    else:
        table_columns = [{"name": "Name", "id": "meta_name"}, {"name": "HMDB ID", "id": "hmdb_id"}] + \
                        [{"name": "Replicate_{}".format(i+1), "id": "replicate_{}".format(i+1)}
                         for i in range(y_num_repli)]
        if dist_type == "normal":
            sample_y = np.round(np.random.normal(y_mean, y_std, y_num_repli), 2)
        elif dist_type == "uniform":
            sample_y = np.round(np.random.uniform(y_min, y_max, y_num_repli), 2)
        else:
            sample_y = []
        y_temp_dict = dict({"meta_name": "Y", "hmdb_id": "/"},
                           **dict(zip(["replicate_{}".format(i + 1) for i in range(y_num_repli)], sample_y)))
        x_cons_dict_list, meta_with_y_list, meta_not_y_list = \
            simulate_continuous_concentrations(table_rows, y_num_repli, sample_y)
        table_data = [y_temp_dict] + x_cons_dict_list
        style_data_conditional = [{
            "if": {"row_index": 0},
            "backgroundColor": "lightsteelblue"
        }] + [{"if": {"row_index": i + 1}, "backgroundColor": "lightblue"} for i in range(len(meta_with_y_list))] + \
                                 [{"if": {"row_index": j + 1 + len(meta_with_y_list)}, "backgroundColor": "powderblue"}
                                  for j in range(len(meta_not_y_list))]
        # if len(meta_with_y_list) != 0:
        fig1 = page_3_plot_conti_line_plot(table_data, meta_with_y_list)
        fig2 = page_3_plot_conti_box_plot(x_cons_dict_list)

        table_header = [html.Thead(html.Tr([html.Th(meta_name) for meta_name in list(map(lambda d: d["meta_name"], x_cons_dict_list))]))]
        row_list = []
        for temp_dict in x_cons_dict_list:
            sample_x = [v for k, v in temp_dict.items() if k != "meta_name" and k != "hmdb_id"]
            slope, intercept, r_value, p_value, std_err = linregress(sample_x, sample_y)
            r_squared = round(r_value ** 2, 2)
            row_list.append(r_squared)
        table_body = [html.Tbody([html.Tr([html.Td(r_s) for r_s in row_list])])]

        return {"display": "block"}, table_data, table_columns, style_data_conditional, {"display": "block"}, \
               fig1, fig2, table_header+table_body, {"display": "block"}, "d-grid gap-2 d-md-flex justify-content-md-center"


def page_3_get_N_HexCol(N):
    HSV_tuples = [(x * 1.0 / N, 0.5, 0.5) for x in range(N)]
    RGB_tuples = list(map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples))
    rgb_out = list(map(lambda x: "".join(("rgb", str(x))), RGB_tuples))
    return rgb_out


def page_3_plot_conti_line_plot(table_data, meta_with_y_list):
    fig = go.Figure()
    y_dict = list(filter(lambda t: t['meta_name'] == "Y", table_data))[0]
    y_values = [value for key, value in y_dict.items() if key != "meta_name" and key != "hmdb_id"]
    c = ['#636EFA','#EF553B','#00CC96','#AB63FA','#FFA15A','#19D3F3','#FF6692','#B6E880','#FF97FF', '#FECB52']
    N = len(table_data)
    rgb_out = page_3_get_N_HexCol(N)
    i = 0
    for meta_dict in meta_with_y_list:
        temp_dict = list(filter(lambda t: t['meta_name'] == meta_dict["meta_name"], table_data))[0]
        x_values = [value for key, value in temp_dict.items() if key != "meta_name" and key != "hmdb_id"]
        y_best_fit = sm.OLS(y_values, sm.add_constant(x_values)).fit().fittedvalues
        fig.add_trace(go.Scatter(
            x=x_values, y=y_values, name=temp_dict["meta_name"],
            mode='markers',
            marker_color=rgb_out[i]
            # c[i]
        ))
        fig.add_trace(go.Scatter(
            x=x_values, y=y_best_fit, name=temp_dict["meta_name"],
            mode='lines',
            line=dict(color=rgb_out[i])
        ))
        i = i + 1
        # print(x_values, y_values)
    fig.update_layout(
        title_text="Line chart between metabolites (X) and outcome (Y)",
        template="plotly_white",
        xaxis_title='X concentrations (uM)',
        yaxis_title='Y',
        boxmode="group",
        # height=500,
        legend=dict(
            itemclick="toggleothers",
            itemdoubleclick="toggle"
        )
    )
    return fig


def page_3_plot_conti_box_plot(table_data):
    fig = go.Figure()
    mixture_list = list(map(lambda d: d["meta_name"], table_data))
    for meta_name in mixture_list:
        temp_dict = list(filter(lambda t: t['meta_name'] == meta_name, table_data))[0]
        fig.add_trace(go.Box(
            y=[value for key, value in temp_dict.items() if key != "meta_name" and key != "hmdb_id"],
            name=meta_name,
        ))
    fig.update_layout(
        title_text="Simulated concentrations for metabolites",
        template="plotly_white",
        yaxis_title='simulated concentrations (uM)',
        boxmode="group",
        # height=500,
        legend=dict(
            itemclick="toggleothers",
            itemdoubleclick="toggle"
        )
    )
    return fig


# ----------------------------- tab 3 callbacks -----------------------------
@app.callback(
    Output("page-3-select-peak-shift-div", "style"),
    Input("page-3-confirm-group-cons-btn", "n_clicks"),
    Input("page-3-confirm-conti-cons-btn", "n_clicks"),
    Input("page-3-confirm-upload-discrete-cons-btn", "n_clicks"),
    Input("page-3-confirm-upload-conti-cons-btn", "n_clicks"),
)
def page_3_update_tab_3_div(group_n_clicks, conti_n_clicks, upload_group_n_clicks, upload_conti_n_clicks):
    if group_n_clicks == 0 and conti_n_clicks == 0 and upload_group_n_clicks == 0 and upload_conti_n_clicks == 0:
        return {"display": "none"}
    else:
        # print("here", upload_group_n_clicks)
        return {"display": "block"}


@app.callback(
    [
        Output("page-3-group-no-peak-shift-results-div", "style"),
        Output("page-3-conti-no-peak-shift-results-div", "style"),
        Output("page-3-group-peak-shift-results-div", "style"),
        Output("page-3-conti-peak-shift-results-div", "style"),
        Output("page-3-upload-group-no-peak-shift-results-div", "style"),
        Output("page-3-upload-group-peak-shift-results-div", "style"),
        Output("page-3-upload-conti-no-peak-shift-results-div", "style"),
        Output("page-3-upload-conti-peak-shift-results-div", "style"),
    ],
    Input("page-3-confirm-group-cons-btn", "n_clicks"),
    Input("page-3-confirm-conti-cons-btn-div", "n_clicks"),
    Input("page-3-confirm-upload-discrete-cons-btn", "n_clicks"),
    Input("page-3-confirm-upload-conti-cons-btn", "n_clicks"),
    Input("page-3-select-spectra-simulation-type", "value"),
    State("page-3-select-simulation-type", "value"),
    State("page-3-select-upload-cons-type", "value")
)
def page_3_update_simulation_results_div(group_n_clicks, conti_n_clicks, upload_group_n_clicks, upload_conti_n_clicks,
                                         sim_spetcra_type, sim_cons_type, upload_cons_type):
    if group_n_clicks == 0 and conti_n_clicks == 0 and upload_group_n_clicks == 0 and upload_conti_n_clicks == 0:
        return {"display": "none"}, {"display": "none"}, {"display": "none"}, {"display": "none"}, \
               {"display": "none"}, {"display": "none"}, {"display": "none"}, {"display": "none"}
    elif group_n_clicks != 0 and sim_cons_type == "group" and sim_spetcra_type == "no_peak_shift":
        return {"display": "block"}, {"display": "none"}, {"display": "none"}, {"display": "none"}, \
               {"display": "none"}, {"display": "none"}, {"display": "none"}, {"display": "none"}
    elif group_n_clicks != 0 and sim_cons_type == "group" and sim_spetcra_type == "peak_shift":
        return {"display": "none"}, {"display": "none"}, {"display": "block"}, {"display": "none"}, \
               {"display": "none"}, {"display": "none"}, {"display": "none"}, {"display": "none"}
    elif conti_n_clicks != 0 and sim_cons_type == "continuous" and sim_spetcra_type == "no_peak_shift":
        return {"display": "none"}, {"display": "block"}, {"display": "none"}, {"display": "none"}, \
               {"display": "none"}, {"display": "none"}, {"display": "none"}, {"display": "none"}
    elif conti_n_clicks != 0 and sim_cons_type == "continuous" and sim_spetcra_type == "peak_shift":
        return {"display": "none"}, {"display": "none"}, {"display": "none"}, {"display": "block"}, \
               {"display": "none"}, {"display": "none"}, {"display": "none"}, {"display": "none"}
    elif upload_group_n_clicks != 0 and upload_cons_type == "discrete" and sim_spetcra_type == "no_peak_shift":
        return {"display": "none"}, {"display": "none"}, {"display": "none"}, {"display": "none"}, \
               {"display": "block"}, {"display": "none"}, {"display": "none"}, {"display": "none"}
    elif upload_group_n_clicks != 0 and upload_cons_type == "discrete" and sim_spetcra_type == "peak_shift":
        return {"display": "none"}, {"display": "none"}, {"display": "none"}, {"display": "none"}, \
               {"display": "none"}, {"display": "block"}, {"display": "none"}, {"display": "none"}
    elif upload_conti_n_clicks != 0 and upload_cons_type == "continuous" and sim_spetcra_type == "no_peak_shift":
        return {"display": "none"}, {"display": "none"}, {"display": "none"}, {"display": "none"}, \
               {"display": "none"}, {"display": "none"}, {"display": "block"}, {"display": "none"}
    elif upload_conti_n_clicks != 0 and upload_cons_type == "continuous" and sim_spetcra_type == "peak_shift":
        return {"display": "none"}, {"display": "none"}, {"display": "none"}, {"display": "none"}, \
               {"display": "none"}, {"display": "none"}, {"display": "none"}, {"display": "block"}
    else:
        return {"display": "none"}, {"display": "none"}, {"display": "none"}, {"display": "none"}, \
               {"display": "none"}, {"display": "none"}, {"display": "none"}, {"display": "none"}


# -------------- tab 3 callbacks: simulate group spectra no peak shift --------------
@app.callback(
    [
        Output("page-3-select-mix", "options"),
        Output("page-3-select-mix", "value"),
    ],
    Input("page-3-confirm-param-btn", "n_clicks"),
    State("page-3-db-names-hmdb-ids-dict", "data"),
)
def page_3_update_all_fig_dropdown(n_clicks, final_name_hmdb_id_dict):
    if n_clicks == 0:
        return [], None
    else:
        options = list(final_name_hmdb_id_dict.keys())
        value = options[0]
    return options, value


@app.callback(
    Output("page-3-processed-spectra-data-dict", "data"),
    Input("page-3-confirm-param-btn", "n_clicks"),
    [
        State("page-3-calibration-range", "value"),
        State("page-3-water-range", "value"),
        State("page-3-filter-thres", "value"),
        State("page-3-smooth-thres-m", "value"),
        State("page-3-smooth-thres-n", "value"),
        State("page-3-simulated-cons-table", "data"),
    ]
)
def page_3_update_group_process_data(n_clicks, cal_range, water_range, filter_thres, smooth_thres_m, smooth_thres_n, table_rows):
    if n_clicks == 0:
        return None
    else:
        # mixture_list = list(map(lambda d: d["meta_name"], table_rows))
        mixture_dict = {item['meta_name']: item['hmdb_id'] for item in table_rows}
        mixture_list = list(mixture_dict.keys())
        temp_raw_data_dict = {name: match_data_dict[name] for name in mixture_list}

        # preprocess 2d data
        removed_data_dict = remove_water_calibration(temp_raw_data_dict, x_scale, y_scale, water_range, cal_range)
        filtered_data_dict = filter_noise(removed_data_dict, filter_thres)
        smooth_data_dict = smooth_data(filtered_data_dict, smooth_thres_m, smooth_thres_n)
        norm_data_dict = normalize_data(smooth_data_dict)

        return norm_data_dict


@app.callback(
    Output("page-3-all-meta-fig", "figure"),
    [
        Input("page-3-confirm-param-btn", "n_clicks"),
        Input("page-3-processed-spectra-data-dict", "data"),
        Input("page-3-select-mix", "value"),
    ]
)
def page_3_update_group_results_all_fig(n_clicks, norm_data_dict, select_meta):
    if n_clicks == 0 or norm_data_dict is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        fig = plot_2d_cosy_spectra(select_meta, norm_data_dict, x_scale, y_scale)
        return fig


@app.callback(
    [
        Output("page-3-group-1-final-dict", "data"),
        Output("page-3-group-2-final-dict", "data"),
    ],
    [
        Input("page-3-confirm-param-btn", "n_clicks"),
        Input("page-3-processed-spectra-data-dict", "data"),
    ],
    [
        State("page-3-snr-noise", "value"),
        State("page-3-group-1-num-repli", "value"),
        State("page-3-group-2-num-repli", "value"),
        State("page-3-simulated-cons-table", "data"),
    ]
)
def page_3_update_groups_mix_fig(n_clicks, norm_data_dict, snr, num_repli_1, num_repli_2, table_rows):
    if n_clicks == 0 or norm_data_dict is None:
        return None, None
    else:
        mixture_dict = {item['meta_name']: item['hmdb_id'] for item in table_rows}
        final_data_dict_1 = simulate_cosy_mixture_for_all_repli(num_repli_1, mixture_dict, norm_data_dict,
                                                           table_rows, protons_df, snr, "1")
        final_data_dict_2 = simulate_cosy_mixture_for_all_repli(num_repli_2, mixture_dict, norm_data_dict,
                                                           table_rows, protons_df, snr, "2")
        # print(final_data_dict_1)
        # print(final_data_dict_2)
        return final_data_dict_1, final_data_dict_2


@app.callback(
    [
        Output("page-3-select-repli-1", "options"),
        Output("page-3-select-repli-1", "value"),
        Output("page-3-select-repli-2", "options"),
        Output("page-3-select-repli-2", "value"),
    ],
    [
        Input("page-3-confirm-param-btn", "n_clicks"),
        Input("page-3-group-1-final-dict", "data"),
        Input("page-3-group-2-final-dict", "data")
    ]
)
def page_3_update_group_replicate_selection(n_clicks, final_data_dict_1, final_data_dict_2):
    if n_clicks == 0 or final_data_dict_1 is None or final_data_dict_2 is None:
        # print("not arrive here!!!!!!!")
        return [], None, [], None
    else:
        # print("arrive here!!!!!!!!!!")
        options_1 = [{'label': i, 'value': i} for i in final_data_dict_1.keys()]
            # list(final_data_dict_1.keys())
        # [{'label': i, 'value': i} for i in final_data_dict_1.keys()]
        value_1 = "replicate_1"
        options_2 = [{'label': i, 'value': i} for i in final_data_dict_2.keys()]
            # list(final_data_dict_2.keys())
        # [{'label': i, 'value': i} for i in final_data_dict_2.keys()]
        value_2 = "replicate_1"
        # print(options_1, options_2)
        return options_1, value_1, options_2, value_2


@app.callback(
    Output("page-3-group-1-spectra-fig", "figure"),
    [
        Input("page-3-confirm-param-btn", "n_clicks"),
        Input("page-3-group-1-final-dict", "data"),
        Input("page-3-select-repli-1", "value"),
        # Input("page-3-contour-level-1", "value")
    ]
)
def page_3_update_group_1_figure(n_clicks, final_data_dict_1, select_repli_1):
    if n_clicks == 0 or final_data_dict_1 is None or select_repli_1 is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        return plot_cosy_repli_spectra(final_data_dict_1, select_repli_1, x_scale, y_scale)


@app.callback(
    Output("page-3-group-2-spectra-fig", "figure"),
    [
        Input("page-3-confirm-param-btn", "n_clicks"),
        Input("page-3-group-2-final-dict", "data"),
        Input("page-3-select-repli-2", "value"),
        # Input("page-3-contour-level-2", "value")
    ]
)
def page_3_update_group_2_figure(n_clicks, final_data_dict_2, select_repli_2):
    if n_clicks == 0 or final_data_dict_2 is None or select_repli_2 is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        return plot_cosy_repli_spectra(final_data_dict_2, select_repli_2, x_scale, y_scale)


@app.callback(
    Output("page-3-group-1-download-spectra-csv", "data"),
    Input("page-3-group-1-btn-csv", "n_clicks"),
    State("page-3-group-1-final-dict", "data"),
    prevent_initial_call=True,
)
def page_3_download_group_1_simulated_spectra(n_clicks, final_replicate_dict):

    final_df = page_3_construct_2d_df_without_peak_shift(final_replicate_dict)
    return dcc.send_data_frame(final_df.to_csv, "group_1_simulated_jres_spectra_without_peak_shift.csv")


def page_3_construct_2d_df_without_peak_shift(final_replicate_dict):
    repli_list = [[key] * len(y_scale) for key in final_replicate_dict.keys()]
    index_1 = list(itertools.chain(*repli_list))
    index_2 = list(y_scale) * len(final_replicate_dict)
    arrays = [index_1, index_2]

    df_list = []
    for key, value in final_replicate_dict.items():
        df_list.append(pd.DataFrame(np.array(value)))

    final_df = pd.concat(df_list, axis=0)
    final_df.columns = x_scale
    final_df.index = arrays
    return final_df


@app.callback(
    Output("page-3-group-2-download-spectra-csv", "data"),
    Input("page-3-group-2-btn-csv", "n_clicks"),
    State("page-3-group-2-final-dict", "data"),
    prevent_initial_call=True,
)
def page_3_download_group_2_simulated_spectra(n_clicks, final_replicate_dict):
    final_df = page_3_construct_2d_df_without_peak_shift(final_replicate_dict)
    return dcc.send_data_frame(final_df.to_csv, "group_2_simulated_jres_spectra_without_peak_shift.csv")


# -------------- tab 3 callbacks: simulate group spectra with peak shift --------------
@app.callback(
    [
        Output("page-3-group-1-ph-same-div", "style"),
        Output("page-3-group-1-ph-not-same-div", "style"),
    ],
    Input("page-3-group-1-ph-same-flag", "value")
)
def page_3_show_group_1_ph_select_div(ph_same_flag):
    if ph_same_flag == "true":
        return {"display": "block"}, {"display": "none"}
    elif ph_same_flag == "false":
        return {"display": "none"}, {"display": "block"}
    else:
        return {"display": "none"}, {"display": "none"}


@app.callback(
    [
        Output("page-3-group-2-ph-same-div", "style"),
        Output("page-3-group-2-ph-not-same-div", "style"),
    ],
    Input("page-3-group-2-ph-same-flag", "value")
)
def page_3_show_group_2_ph_select_div(ph_same_flag):
    if ph_same_flag == "true":
        return {"display": "block"}, {"display": "none"}
    elif ph_same_flag == "false":
        return {"display": "none"}, {"display": "block"}
    else:
        return {"display": "none"}, {"display": "none"}


@app.callback(
    [
        Output("page-3-group-cons-ph-table-div", "style"),
        Output("page-3-group-cons-ph-table", "data"),
        Output("page-3-group-cons-ph-table", "columns"),
        Output("page-3-group-cons-ph-table", "style_data_conditional"),
        Output("page-3-group-peak-shift-part-2-div", "style"),
        Output("page-3-group-peak-shift-spectra-div", "style")
    ],
    Input("page-3-tab-3-group-submit-part-1-btn", "n_clicks"),
    [
        State("page-3-simulated-cons-table", "data"),
        State("page-3-simulated-cons-table", "columns"),
        State("page-3-group-1-num-repli", "value"),
        State("page-3-group-1-ph-same-flag", "value"),
        State("page-3-group-1-same-ph", "value"),
        State("page-3-group-1-ph-mean", "value"),
        State("page-3-group-1-ph-std", "value"),
        State("page-3-group-2-num-repli", "value"),
        State("page-3-group-2-ph-same-flag", "value"),
        State("page-3-group-2-same-ph", "value"),
        State("page-3-group-2-ph-mean", "value"),
        State("page-3-group-2-ph-std", "value")
    ]
)
def page_3_update_group_cons_ph_table(n_clicks, cons_table_data, cons_table_columns, num_repli_1, ph_same_flag_1, same_ph_1,
                               ph_mean_1, ph_std_1, num_repli_2, ph_same_flag_2, same_ph_2, ph_mean_2, ph_std_2):
    if n_clicks == 0:
        return {"display": "none"}, [{}], [], [{}], {"display": "none"}, {"display": "none"}
    else:
        ph_list_1 = page_3_sample_ph_list(num_repli_1, ph_same_flag_1, same_ph_1, ph_mean_1, ph_std_1)
        ph_list_2 = page_3_sample_ph_list(num_repli_2, ph_same_flag_2, same_ph_2, ph_mean_2, ph_std_2)
        all_ph_list = np.concatenate((["pH", "/"], ph_list_1, ph_list_2))
        name_list = ["meta_name", "hmdb_id"] + ["1_replicate_" + str(i + 1) for i in range(num_repli_1)] + \
                    ["2_replicate_" + str(i + 1) for i in range(num_repli_2)]
        temp_dict = dict(zip(name_list, all_ph_list))
        table_data = [temp_dict] + cons_table_data
        style_data_conditional = [{"if": {"row_index": 0}, "color": "green", "fontWeight": "bold"}]
        return {"display": "block"}, table_data, cons_table_columns, style_data_conditional, {"display": "block"}, {"display": "block"}


def page_3_sample_ph_list(num_replicate, ph_same_flag, same_ph, ph_mean, ph_std):
    if ph_same_flag == "true":
        ph_list = np.round([same_ph] * num_replicate, 2)
    elif ph_same_flag == "false":
        ph_list = np.round(page_3_get_positive_ph(ph_mean, ph_std, num_replicate), 2)
    else:
        ph_list = [7.4] * num_replicate
    return ph_list


def page_3_get_positive_ph(mean, std, num_replicates):
    x = np.random.normal(mean,std, num_replicates)
    if np.all(x > 0):
        return x
    else:
        return page_3_get_positive_ph(mean, std, num_replicates)


@app.callback(
    Output("page-3-group-peak-shift-processed-spectra-data-dict", "data"),
    Input("page-3-group-peak-shift-confirm-param-btn", "n_clicks"),
    [
        State("page-3-group-peak-shift-calibration-range", "value"),
        State("page-3-group-peak-shift-water-range", "value"),
        State("page-3-group-peak-shift-filter-thres", "value"),
        State("page-3-group-peak-shift-smooth-thres-m", "value"),
        State("page-3-group-peak-shift-smooth-thres-n", "value"),
        State("page-3-group-cons-ph-table", "data"),
    ]
)
def page_3_group_peak_shift_update_group_process_data(n_clicks, cal_range, water_range, filter_thres,
                                                      smooth_thres_m, smooth_thres_n, table_rows):
    if n_clicks == 0:
        return None
    else:
        mixture_list = list(map(lambda d: d["meta_name"], table_rows))
        mixture_list.remove('pH')
        temp_raw_data_dict = {name: match_data_dict[name] for name in mixture_list}

        # preprocess 2d data
        removed_data_dict = remove_water_calibration(temp_raw_data_dict, x_scale, y_scale, water_range, cal_range)
        filtered_data_dict = filter_noise(removed_data_dict, filter_thres)
        smooth_data_dict = smooth_data(filtered_data_dict, smooth_thres_m, smooth_thres_n)
        norm_data_dict = normalize_data(smooth_data_dict)
        return norm_data_dict


@app.callback(
    [
        Output("page-3-group-peak-shift-select-mix", "options"),
        Output("page-3-group-peak-shift-select-mix", "value"),
    ],
    Input("page-3-group-peak-shift-confirm-param-btn", "n_clicks"),
    State("page-3-db-names-hmdb-ids-dict", "data"),
)
def page_3_group_peak_shift_update_all_fig_dropdown(n_clicks, final_name_hmdb_id_dict):
    if n_clicks == 0:
        return [], None
    else:
        options = list(final_name_hmdb_id_dict.keys())
        value = options[0]
    return options, value


@app.callback(
    Output("page-3-group-peak-shift-all-meta-fig", "figure"),
    [
        Input("page-3-group-peak-shift-confirm-param-btn", "n_clicks"),
        Input("page-3-group-peak-shift-processed-spectra-data-dict", "data"),
        Input("page-3-group-peak-shift-select-mix", "value"),
    ]
)
def page_3_group_peak_shift_update_group_results_all_fig(n_clicks, norm_data_dict, select_meta):
    if n_clicks == 0 or norm_data_dict is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        fig = plot_2d_cosy_spectra(select_meta, norm_data_dict, x_scale, y_scale)
        return fig


@app.callback(
    [
        Output("page-3-group-peak-shift-group-1-final-dict", "data"),
        Output("page-3-group-peak-shift-group-2-final-dict", "data"),
    ],
    [
        Input("page-3-group-peak-shift-confirm-param-btn", "n_clicks"),
        Input("page-3-group-peak-shift-processed-spectra-data-dict", "data"),
    ],
    [
        State("page-3-db-names-hmdb-ids-dict", "data"),
        State("page-3-group-peak-shift-snr-noise", "value"),
        State("page-3-group-cons-ph-table", "data"),
    ]
)
def page_3_group_peak_shift_update_groups_repli_data(n_clicks, norm_data_dict, final_names_hmdb_id_dict, snr, cons_ph_table_data):
    if n_clicks == 0 or norm_data_dict is None:
        return None, None
    else:
        # mixture_list = list(map(lambda d: d["meta_name"], cons_ph_table_data))
        # mixture_list.remove('pH')
        mixture_dict = {item['meta_name']: item['hmdb_id'] for item in cons_ph_table_data}
        del mixture_dict["pH"]
        mixture_list = list(mixture_dict.keys())

        mixture_pka_dict = page_3_get_mixture_pka_dict(mixture_list, final_names_hmdb_id_dict)
        group_repli_ph_dict = list(filter(lambda d: d["meta_name"] == 'pH', cons_ph_table_data))[0]

        final_mix_data_dict_1 = get_mixture_data_for_all_replicates("1", group_repli_ph_dict, mixture_dict,
                                                                    norm_data_dict, mixture_pka_dict, x_scale,
                                                                    cons_ph_table_data, protons_df, snr)
        final_mix_data_dict_2 = get_mixture_data_for_all_replicates("2", group_repli_ph_dict, mixture_dict,
                                                                    norm_data_dict, mixture_pka_dict, x_scale,
                                                                    cons_ph_table_data, protons_df, snr)
        # print(final_mix_data_dict_1["replicate_1"])
        return final_mix_data_dict_1, final_mix_data_dict_2


def page_3_get_mixture_pka_dict(mixture_list, final_names_hmdb_id_dict):
    mixture_pka_dict = dict()
    for meta_name in mixture_list:
        temp_pka_dict = hmdb_id_pka_dict[final_names_hmdb_id_dict[meta_name]]
        try:
            temp_pka_1 = float(temp_pka_dict["pka_strongest_basic"])
        except:
            temp_pka_1 = None
        try:
            temp_pka_2 = float(temp_pka_dict["pka_strongest_acidic"])
        except:
            temp_pka_2 = None
        if temp_pka_1 is not None and temp_pka_2 is not None:
            temp_pka = (temp_pka_1 + temp_pka_2) / 2
        elif temp_pka_1 is None and temp_pka_2 is not None:
            temp_pka = temp_pka_2
        elif temp_pka_1 is not None and temp_pka_2 is None:
            temp_pka = temp_pka_1
        else:
            temp_pka = np.random.normal(6.013, 2.972, 1)[0]
        mixture_pka_dict[meta_name] = temp_pka
    return mixture_pka_dict


@app.callback(
    [
        Output("page-3-group-peak-shift-select-repli-1", "options"),
        Output("page-3-group-peak-shift-select-repli-1", "value"),
        Output("page-3-group-peak-shift-select-repli-2", "options"),
        Output("page-3-group-peak-shift-select-repli-2", "value"),
    ],
    Input("page-3-group-peak-shift-confirm-param-btn", "n_clicks"),
    [
        State("page-3-group-1-num-repli", "value"),
        State("page-3-group-2-num-repli", "value"),
    ]
)
def page_3_group_peak_shift_update_group_replicate_selection(n_clicks, num_repli_1, num_repli_2):
    if n_clicks == 0:
        # print("not arrive here!!!!!!!")
        return [], None, [], None
    else:
        # print("arrive here!!!!!!!!!!")
        options_1 = [{'label': f"replicate_{i+1}", 'value': f"replicate_{i+1}"} for i in range(num_repli_1)]
        value_1 = "replicate_1"
        options_2 = [{'label': f"replicate_{i+1}", 'value': f"replicate_{i+1}"} for i in range(num_repli_2)]
        value_2 = "replicate_1"
        # print(options_1, options_2)
        return options_1, value_1, options_2, value_2


@app.callback(

    Output("page-3-group-peak-shift-group-1-spectra-fig", "figure"),
    [
        Input("page-3-group-peak-shift-confirm-param-btn", "n_clicks"),
        Input("page-3-group-peak-shift-group-1-final-dict", "data"),
        Input("page-3-group-peak-shift-select-repli-1", "value"),
    ]
)
def page_3_group_peak_shift_update_group_1_figure(n_clicks, final_data_dict_1, select_repli_1):
    if n_clicks == 0 or final_data_dict_1 is None or select_repli_1 is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        return plot_cosy_repli_spectra_with_ph(final_data_dict_1, select_repli_1, x_scale, y_scale)


@app.callback(

    Output("page-3-group-peak-shift-group-2-spectra-fig", "figure"),
    [
        Input("page-3-group-peak-shift-confirm-param-btn", "n_clicks"),
        Input("page-3-group-peak-shift-group-2-final-dict", "data"),
        Input("page-3-group-peak-shift-select-repli-2", "value"),
    ]
)
def page_3_group_peak_shift_update_group_2_figure(n_clicks, final_data_dict_2, select_repli_2):
    if n_clicks == 0 or final_data_dict_2 is None or select_repli_2 is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        return plot_cosy_repli_spectra_with_ph(final_data_dict_2, select_repli_2, x_scale, y_scale)


@app.callback(
    Output("page-3-group-peak-shift-group-1-download-spectra-csv", "data"),
    Input("page-3-group-peak-shift-group-1-btn-csv", "n_clicks"),
    State("page-3-group-peak-shift-group-1-final-dict", "data"),
    prevent_initial_call=True,
)
def page_3_group_peak_shift_download_group_1_simulated_spectra(n_clicks, final_replicate_dict):
    final_df = page_3_construct_2d_df_with_peak_shift(final_replicate_dict)
    return dcc.send_data_frame(final_df.to_csv, "group_1_simulated_cosy_spectra_with_peak_shift.csv")


def page_3_construct_2d_df_with_peak_shift(final_replicate_dict):
    repli_list = [[key] * len(y_scale) for key in final_replicate_dict.keys()]
    index_1 = list(itertools.chain(*repli_list))
    index_2 = list(y_scale) * len(final_replicate_dict)
    arrays = [index_1, index_2]

    df_list = []
    for key, value in final_replicate_dict.items():
        df_list.append(pd.DataFrame(np.array(value[1])))

    final_df = pd.concat(df_list, axis=0)
    final_df.columns = x_scale
    final_df.index = arrays
    return final_df


@app.callback(
    Output("page-3-group-peak-shift-group-2-download-spectra-csv", "data"),
    Input("page-3-group-peak-shift-group-2-btn-csv", "n_clicks"),
    State("page-3-group-peak-shift-group-2-final-dict", "data"),
    prevent_initial_call=True,
)
def page_3_group_peak_shift_download_group_2_simulated_spectra(n_clicks, final_replicate_dict):
    final_df = page_3_construct_2d_df_with_peak_shift(final_replicate_dict)
    return dcc.send_data_frame(final_df.to_csv, "group_2_simulated_cosy_spectra_with_peak_shift.csv")


# -------------- tab 3 callbacks: simulate continuous spectra no peak shift --------------
@app.callback(
    Output("page-3-conti-processed-spectra-data-dict", "data"),
    Input("page-3-conti-confirm-param-btn", "n_clicks"),
    [
        State("page-3-conti-calibration-range", "value"),
        State("page-3-conti-water-range", "value"),
        State("page-3-conti-filter-thres", "value"),
        State("page-3-conti-smooth-thres-m", "value"),
        State("page-3-conti-smooth-thres-n", "value"),
        State("page-3-conti-simulated-cons-table", "data"),
    ]
)
def page_3_update_continuous_process_data(n_clicks, cal_range, water_range, filter_thres, smooth_thres_m, smooth_thres_n, table_rows):
    if n_clicks == 0:
        return None
    else:
        mixture_list = list(map(lambda d: d["meta_name"], table_rows))
        mixture_list.remove("Y")
        temp_raw_data_dict = {name: match_data_dict[name] for name in mixture_list}

        # preprocess 2d data
        removed_data_dict = remove_water_calibration(temp_raw_data_dict, x_scale, y_scale, water_range, cal_range)
        filtered_data_dict = filter_noise(removed_data_dict, filter_thres)
        smooth_data_dict = smooth_data(filtered_data_dict, smooth_thres_m, smooth_thres_n)
        norm_data_dict = normalize_data(smooth_data_dict)
        return norm_data_dict


@app.callback(
    [
        Output("page-3-conti-select-mix", "options"),
        Output("page-3-conti-select-mix", "value"),
    ],
    Input("page-3-conti-confirm-param-btn", "n_clicks"),
    State("page-3-db-names-hmdb-ids-dict", "data"),
)
def page_3_update_continuous_all_fig_dropdown(n_clicks, final_name_hmdb_id_dict):
    if n_clicks == 0:
        return [], None
    else:
        options = list(final_name_hmdb_id_dict.keys())
        value = options[0]
    return options, value


@app.callback(
    Output("page-3-conti-all-meta-fig", "figure"),
    [
        Input("page-3-conti-confirm-param-btn", "n_clicks"),
        Input("page-3-conti-processed-spectra-data-dict", "data"),
        Input("page-3-conti-select-mix", "value"),
    ]
)
def page_3_update_continuous_results_all_fig(n_clicks, norm_data_dict, select_meta):
    if n_clicks == 0 or norm_data_dict is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        fig = plot_2d_cosy_spectra(select_meta, norm_data_dict, x_scale, y_scale)
        return fig


@app.callback(
    [
        Output("page-3-conti-final-replicate-dict", "data"),
        Output("page-3-conti-select-replicate", "options"),
        Output("page-3-conti-select-replicate", "value")
    ],
    [
        Input("page-3-conti-confirm-param-btn", "n_clicks"),
        Input("page-3-conti-processed-spectra-data-dict", "data")
    ],
    [
        State("page-3-conti-snr-noise", "value"),
        State("page-3-y-num-repli", "value"),
        State("page-3-conti-simulated-cons-table", "data"),
    ]
)
def page_3_update_continuous_final_replicate_dict(n_clicks, norm_data_dict, snr, y_num_repli, table_rows):
    if n_clicks == 0 or norm_data_dict is None:
        return None, [], None
    else:
        # mixture_list = list(map(lambda d: d["meta_name"], table_rows))
        # mixture_list.remove("Y")
        mixture_dict = {item['meta_name']: item['hmdb_id'] for item in table_rows}
        del mixture_dict["Y"]
        final_mix_data_dict = simulate_continuous_cosy_mixture_for_all_repli(y_num_repli, mixture_dict,
                                                                             norm_data_dict, table_rows,
                                                                             protons_df, snr)
        options = [{"label": i, "value": i} for i in list(final_mix_data_dict.keys())]
        value = list(final_mix_data_dict.keys())[0]
        return final_mix_data_dict, options, value


@app.callback(
    Output("page-3-conti-mix-spectra-fig", "figure"),
    Input("page-3-conti-select-replicate", "value"),
    [
        State("page-3-conti-confirm-param-btn", "n_clicks"),
        State("page-3-conti-final-replicate-dict", "data"),
    ]
)
def page_3_update_continuous_mix_fig(select_replicate, n_clicks, final_data_dict):
    if n_clicks == 0 or select_replicate is None or final_data_dict is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        fig = plot_cosy_repli_spectra(final_data_dict, select_replicate, x_scale, y_scale)
        return fig


@app.callback(
    Output("page-3-conti-meta-y-fig", "figure"),
    Input("page-3-conti-select-replicate", "value"),
    [
        State("page-3-conti-confirm-param-btn", "n_clicks"),
        State("page-3-conti-simulated-cons-table", "data"),
    ]
)
def page_3_update_continuous_meta_y_line_plot(select_replicate, n_clicks, cons_table_rows):
    if n_clicks == 0 or select_replicate is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        fig = page_3_plot_conti_results_line_plot(cons_table_rows, select_replicate)
        return fig


def page_3_plot_conti_results_line_plot(table_data, select_replicate):
    fig = go.Figure()
    y_dict = list(filter(lambda t: t['meta_name'] == "Y", table_data))[0]
    y_values = [value for key, value in y_dict.items() if key != "meta_name" and key != "hmdb_id"]
    # c = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
    N = len(table_data)
    rgb_out = page_3_get_N_HexCol(N)
    # print(rgb_out)
    mixture_list = list(map(lambda d: d["meta_name"], table_data))
    mixture_list.remove("Y")
    i = 0
    for meta_name in mixture_list:
        temp_dict = list(filter(lambda t: t['meta_name'] == meta_name, table_data))[0]
        x_values = [value for key, value in temp_dict.items() if key != "meta_name" and key != "hmdb_id"]
        y_best_fit = sm.OLS(y_values, sm.add_constant(x_values)).fit().fittedvalues
        fig.add_trace(go.Scatter(
            x=x_values, y=y_values, name=temp_dict["meta_name"],
            mode='markers',
            marker_color=rgb_out[i]
        ))
        fig.add_trace(go.Scatter(
            x=x_values, y=y_best_fit, name=temp_dict["meta_name"],
            mode='lines',
            line=dict(color=rgb_out[i])
        ))
        i = i + 1
    highlight_y = float(y_dict[select_replicate])
    fig.add_hline(
        y=highlight_y,
        line_dash="dash",
        line_color="darkred"
    )
    fig.update_layout(
        title_text="Line chart between metabolites (X) and outcome (Y)",
        template="plotly_white",
        xaxis_title='X concentrations (uM)',
        yaxis_title='Y',
        legend=dict(
            itemclick="toggleothers",
            itemdoubleclick="toggle"
        ),
        height=700,
    )
    return fig


@app.callback(
    Output("page-3-conti-download-spectra-csv", "data"),
    Input("page-3-conti-btn-csv", "n_clicks"),
    State("page-3-conti-final-replicate-dict", "data"),
    prevent_initial_call=True,
)
def page_3_download_continuous_simulated_spectra(n_clicks, final_replicate_dict):
    final_df = page_3_construct_2d_df_without_peak_shift(final_replicate_dict)
    return dcc.send_data_frame(final_df.to_csv, "continuous_simulated_cosy_spectra_without_peak_shift.csv")


# -------------- tab 3 callbacks: simulate continuous spectra with peak shift --------------
@app.callback(
    [
        Output("page-3-conti-ph-same-div", "style"),
        Output("page-3-conti-ph-not-same-div", "style"),
    ],
    Input("page-3-conti-ph-same-flag", "value")
)
def page_3_show_conti_ph_select_div(ph_same_flag):
    if ph_same_flag == "true":
        return {"display": "block"}, {"display": "none"}
    elif ph_same_flag == "false":
        return {"display": "none"}, {"display": "block"}
    else:
        return {"display": "none"}, {"display": "none"}


@app.callback(
    [
        Output("page-3-conti-cons-ph-table-div", "style"),
        Output("page-3-conti-cons-ph-table", "data"),
        Output("page-3-conti-cons-ph-table", "columns"),
        Output("page-3-conti-cons-ph-table", "style_data_conditional"),
        Output("page-3-conti-peak-shift-part-2-div", "style"),
        Output("page-3-conti-peak-shift-spectra-div", "style")
    ],
    Input("page-3-tab-3-conti-submit-part-1-btn", "n_clicks"),
    [
        State("page-3-conti-simulated-cons-table", "data"),
        State("page-3-conti-simulated-cons-table", "columns"),
        State("page-3-y-num-repli", "value"),
        State("page-3-conti-ph-same-flag", "value"),
        State("page-3-conti-same-ph", "value"),
        State("page-3-conti-ph-mean", "value"),
        State("page-3-conti-ph-std", "value"),
    ]
)
def page_3_update_continuous_cons_ph_table(n_clicks, cons_table_data, cons_table_columns, y_num_repli, ph_same_flag, same_ph,
                                    ph_mean, ph_std):
    if n_clicks == 0:
        return {"display": "none"}, [{}], [], [{}], {"display": "none"}, {"display": "none"}
    else:
        ph_list = page_3_sample_ph_list(y_num_repli, ph_same_flag, same_ph, ph_mean, ph_std)
        all_ph_list = np.concatenate((["pH", "/"], ph_list))
        name_list = ["meta_name", "hmdb_id"] + ["replicate_" + str(i + 1) for i in range(y_num_repli)]
        temp_dict = dict(zip(name_list, all_ph_list))
        table_data = [temp_dict] + cons_table_data
        style_data_conditional = [{"if": {"row_index": 0}, "color": "green", "fontWeight": "bold"}]
        return {"display": "block"}, table_data, cons_table_columns, style_data_conditional, {"display": "block"}, {
            "display": "block"}


@app.callback(
    Output("page-3-conti-peak-shift-processed-spectra-data-dict", "data"),
    Input("page-3-conti-peak-shift-confirm-param-btn", "n_clicks"),
    [
        State("page-3-conti-peak-shift-calibration-range", "value"),
        State("page-3-conti-peak-shift-water-range", "value"),
        State("page-3-conti-peak-shift-filter-thres", "value"),
        State("page-3-conti-peak-shift-smooth-thres-m", "value"),
        State("page-3-conti-peak-shift-smooth-thres-n", "value"),
        State("page-3-conti-cons-ph-table", "data"),
    ]
)
def page_3_continuous_peak_shift_update_processed_data(n_clicks, cal_range, water_range, filter_thres,
                                                      smooth_thres_m, smooth_thres_n, table_rows):
    if n_clicks == 0:
        return None
    else:
        mixture_list = list(map(lambda d: d["meta_name"], table_rows))
        mixture_list.remove("Y")
        mixture_list.remove('pH')
        temp_raw_data_dict = {name: match_data_dict[name] for name in mixture_list}
        # preprocess 2d data
        removed_data_dict = remove_water_calibration(temp_raw_data_dict, x_scale, y_scale, water_range, cal_range)
        filtered_data_dict = filter_noise(removed_data_dict, filter_thres)
        smooth_data_dict = smooth_data(filtered_data_dict, smooth_thres_m, smooth_thres_n)
        norm_data_dict = normalize_data(smooth_data_dict)
        return norm_data_dict


@app.callback(
    [
        Output("page-3-conti-peak-shift-select-mix", "options"),
        Output("page-3-conti-peak-shift-select-mix", "value"),
    ],
    Input("page-3-conti-peak-shift-confirm-param-btn", "n_clicks"),
    State("page-3-db-names-hmdb-ids-dict", "data"),
)
def page_3_continuous_peak_shift_update_all_fig_dropdown(n_clicks, final_name_hmdb_id_dict):
    if n_clicks == 0:
        return [], None
    else:
        options = list(final_name_hmdb_id_dict.keys())
        value = options[0]
    return options, value


@app.callback(
    Output("page-3-conti-peak-shift-all-meta-fig", "figure"),
    [
        Input("page-3-conti-peak-shift-confirm-param-btn", "n_clicks"),
        Input("page-3-conti-peak-shift-processed-spectra-data-dict", "data"),
        Input("page-3-conti-peak-shift-select-mix", "value"),
    ]
)
def page_3_continuous_peak_shift_update_group_results_all_fig(n_clicks, norm_data_dict, select_meta):
    if n_clicks == 0 or norm_data_dict is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        fig = plot_2d_cosy_spectra(select_meta, norm_data_dict, x_scale, y_scale)
        return fig


@app.callback(
    Output("page-3-conti-peak-shift-final-replicate-dict", "data"),
    [
        Input("page-3-conti-peak-shift-confirm-param-btn", "n_clicks"),
        Input("page-3-conti-peak-shift-processed-spectra-data-dict", "data"),
    ],
    [
        State("page-3-db-names-hmdb-ids-dict", "data"),
        State("page-3-conti-peak-shift-snr-noise", "value"),
        State("page-3-conti-cons-ph-table", "data"),
    ]
)
def page_3_continous_peak_shift_update_groups_repli_data(n_clicks, norm_data_dict, final_names_hmdb_id_dict, snr, cons_ph_table_data):
    if n_clicks == 0 or norm_data_dict is None:
        return None, None
    else:
        mixture_dict = {item['meta_name']: item['hmdb_id'] for item in cons_ph_table_data}
        del mixture_dict["pH"]
        del mixture_dict["Y"]
        mixture_list = list(mixture_dict.keys())

        mixture_pka_dict = page_3_get_mixture_pka_dict(mixture_list, final_names_hmdb_id_dict)
        conti_repli_ph_dict = list(filter(lambda d: d["meta_name"] == 'pH', cons_ph_table_data))[0]

        final_mix_data_dict = conti_get_mixture_data_for_all_replicates(conti_repli_ph_dict, mixture_dict,
                                                                        norm_data_dict, mixture_pka_dict, x_scale,
                                                                        cons_ph_table_data, protons_df, snr)
        return final_mix_data_dict


@app.callback(
    [
        Output("page-3-conti-peak-shift-select-replicate", "options"),
        Output("page-3-conti-peak-shift-select-replicate", "value"),
    ],
    Input("page-3-conti-peak-shift-confirm-param-btn", "n_clicks"),
    State("page-3-y-num-repli", "value"),
)
def page_3_group_peak_shift_update_group_replicate_selection(n_clicks, y_num_repli):
    if n_clicks == 0:
        return [], None
    else:
        options = [{'label': f"replicate_{i + 1}", 'value': f"replicate_{i + 1}"} for i in range(y_num_repli)]
        value = "replicate_1"
        return options, value


@app.callback(
    Output("page-3-conti-peak-shift-mix-spectra-fig", "figure"),
    [
        Input("page-3-conti-peak-shift-confirm-param-btn", "n_clicks"),
        Input("page-3-conti-peak-shift-final-replicate-dict", "data"),
        Input("page-3-conti-peak-shift-select-replicate", "value"),
    ]
)
def page_3_continuous_peak_shift_update_mix_figure(n_clicks, final_data_dict, select_repli):
    if n_clicks == 0 or final_data_dict is None or select_repli is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        return plot_cosy_repli_spectra_with_ph(final_data_dict, select_repli, x_scale, y_scale)


@app.callback(
    Output("page-3-conti-peak-shift-meta-y-fig", "figure"),
    Input("page-3-conti-peak-shift-select-replicate", "value"),
    [
        State("page-3-conti-peak-shift-confirm-param-btn", "n_clicks"),
        State("page-3-conti-simulated-cons-table", "data"),
    ]
)
def page_3_update_continuous_peak_shift_meta_y_line_plot(select_replicate, n_clicks, cons_table_rows):
    if n_clicks == 0 or select_replicate is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        fig = page_3_plot_conti_results_line_plot(cons_table_rows, select_replicate)
        return fig


@app.callback(
    Output("page-3-conti-peak-shift-download-spectra-csv", "data"),
    Input("page-3-conti-peak-shift-btn-csv", "n_clicks"),
    State("page-3-conti-peak-shift-final-replicate-dict", "data"),
    prevent_initial_call=True,
)
def page_3_download_conti_peak_shift_simulated_spectra(n_clicks, final_replicate_dict):
    final_df = page_3_construct_2d_df_with_peak_shift(final_replicate_dict)
    return dcc.send_data_frame(final_df.to_csv, "continuous_simulated_cosy_spectra_with_peak_shift.csv")


# -------------- tab 3 callbacks: simulate upload group spectra no peak shift --------------
@app.callback(
    [
        Output("page-3-upload-group-no-peak-shift-select-mix", "options"),
        Output("page-3-upload-group-no-peak-shift-select-mix", "value"),
    ],
    Input("page-3-upload-group-no-peak-shift-confirm-param-btn", "n_clicks"),
    State("page-3-db-names-hmdb-ids-dict", "data"),
)
def page_3_update_upload_all_fig_dropdown(n_clicks, final_name_hmdb_id_dict):
    if n_clicks == 0:
        return [], None
    else:
        options = list(final_name_hmdb_id_dict.keys())
        value = options[0]
    return options, value,


@app.callback(

    Output("page-3-upload-group-no-peak-shift-processed-spectra-data-dict", "data"),
    Input("page-3-upload-group-no-peak-shift-confirm-param-btn", "n_clicks"),
    [
        State("page-3-upload-group-no-peak-shift-calibration-range", "value"),
        State("page-3-upload-group-no-peak-shift-water-range", "value"),
        State("page-3-upload-group-no-peak-shift-filter-thres", "value"),
        State("page-3-upload-group-no-peak-shift-smooth-thres-m", "value"),
        State("page-3-upload-group-no-peak-shift-smooth-thres-n", "value"),
        State("page-3-upload-discrete-cons-table", "data"),
    ]
)
def page_3_update_upload_group_process_data(n_clicks, cal_range, water_range, filter_thres, smooth_thres_m, smooth_thres_n, table_rows):
    if n_clicks == 0:
        return None
    else:
        mixture_list = list(map(lambda d: d["meta_name"], table_rows))
        temp_raw_data_dict = {name: match_data_dict[name] for name in mixture_list}

        # preprocess 2d data
        removed_data_dict = remove_water_calibration(temp_raw_data_dict, x_scale, y_scale, water_range, cal_range)
        filtered_data_dict = filter_noise(removed_data_dict, filter_thres)
        smooth_data_dict = smooth_data(filtered_data_dict, smooth_thres_m, smooth_thres_n)
        norm_data_dict = normalize_data(smooth_data_dict)
        return norm_data_dict


@app.callback(
    Output("page-3-upload-group-no-peak-shift-all-meta-fig", "figure"),
    [
        Input("page-3-upload-group-no-peak-shift-confirm-param-btn", "n_clicks"),
        Input("page-3-upload-group-no-peak-shift-processed-spectra-data-dict", "data"),
        Input("page-3-upload-group-no-peak-shift-select-mix", "value"),
    ]
)
def page_3_update_upload_group_results_all_fig(n_clicks, norm_data_dict, select_meta):
    if n_clicks == 0 or norm_data_dict is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        fig = plot_2d_cosy_spectra(select_meta, norm_data_dict, x_scale, y_scale)
        return fig


@app.callback(
    [
        Output("page-3-upload-group-no-peak-shift-group-1-final-dict", "data"),
        Output("page-3-upload-group-no-peak-shift-group-2-final-dict", "data"),
    ],
    [
        Input("page-3-upload-group-no-peak-shift-confirm-param-btn", "n_clicks"),
        Input("page-3-upload-group-no-peak-shift-processed-spectra-data-dict", "data"),
    ],
    [
        State("page-3-upload-group-no-peak-shift-snr-noise", "value"),
        State("page-3-upload-discrete-cons-table", "data"),
    ]
)
def page_3_update_upload_groups_mix_fig(n_clicks, norm_data_dict, snr, table_rows):
    num_repli_1 = sum(1 for i in table_rows[0].keys() if i[0] == "1")
    num_repli_2 = sum(1 for i in table_rows[0].keys() if i[0] == "2")
    if n_clicks == 0 or norm_data_dict is None:
        return None, None
    else:
        mixture_dict = {item['meta_name']: item['hmdb_id'] for item in table_rows}
        final_data_dict_1 = simulate_cosy_mixture_for_all_repli(num_repli_1, mixture_dict, norm_data_dict,
                                                                table_rows, protons_df, snr, "1")
        final_data_dict_2 = simulate_cosy_mixture_for_all_repli(num_repli_2, mixture_dict, norm_data_dict,
                                                                table_rows, protons_df, snr, "2")
        return final_data_dict_1, final_data_dict_2


@app.callback(
    [
        Output("page-3-upload-group-no-peak-shift-select-repli-1", "options"),
        Output("page-3-upload-group-no-peak-shift-select-repli-1", "value"),
        Output("page-3-upload-group-no-peak-shift-select-repli-2", "options"),
        Output("page-3-upload-group-no-peak-shift-select-repli-2", "value"),
    ],
    [
        Input("page-3-upload-group-no-peak-shift-confirm-param-btn", "n_clicks"),
        Input("page-3-upload-group-no-peak-shift-group-1-final-dict", "data"),
        Input("page-3-upload-group-no-peak-shift-group-2-final-dict", "data")
    ]
)
def page_3_update_upload_group_replicate_selection(n_clicks, final_data_dict_1, final_data_dict_2):
    if n_clicks == 0 or final_data_dict_1 is None or final_data_dict_2 is None:
        # print("not arrive here!!!!!!!")
        return [], None, [], None
    else:
        # print("arrive here!!!!!!!!!!")
        options_1 = [{'label': i, 'value': i} for i in final_data_dict_1.keys()]
        value_1 = "replicate_1"
        options_2 = [{'label': i, 'value': i} for i in final_data_dict_2.keys()]
        value_2 = "replicate_1"
        return options_1, value_1, options_2, value_2


@app.callback(
    Output("page-3-upload-group-no-peak-shift-group-1-spectra-fig", "figure"),
    [
        Input("page-3-upload-group-no-peak-shift-confirm-param-btn", "n_clicks"),
        Input("page-3-upload-group-no-peak-shift-group-1-final-dict", "data"),
        Input("page-3-upload-group-no-peak-shift-select-repli-1", "value"),
    ]
)
def page_3_update_upload_group_1_figure(n_clicks, final_data_dict_1, select_repli_1):
    if n_clicks == 0 or final_data_dict_1 is None or select_repli_1 is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        return plot_cosy_repli_spectra(final_data_dict_1, select_repli_1, x_scale, y_scale)


@app.callback(
    Output("page-3-upload-group-no-peak-shift-group-2-spectra-fig", "figure"),
    [
        Input("page-3-upload-group-no-peak-shift-confirm-param-btn", "n_clicks"),
        Input("page-3-upload-group-no-peak-shift-group-2-final-dict", "data"),
        Input("page-3-upload-group-no-peak-shift-select-repli-2", "value"),
    ]
)
def page_3_update_upload_group_2_figure(n_clicks, final_data_dict_2, select_repli_2):
    if n_clicks == 0 or final_data_dict_2 is None or select_repli_2 is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        return plot_cosy_repli_spectra(final_data_dict_2, select_repli_2, x_scale, y_scale)


@app.callback(
    Output("page-3-upload-group-no-peak-shift-group-1-download-spectra-csv", "data"),
    Input("page-3-upload-group-no-peak-shift-group-1-btn-csv", "n_clicks"),
    State("page-3-upload-group-no-peak-shift-group-1-final-dict", "data"),
    prevent_initial_call=True,
)
def page_3_download_upload_group_1_simulated_spectra(n_clicks, final_replicate_dict):
    final_df = page_3_construct_2d_df_without_peak_shift(final_replicate_dict)
    return dcc.send_data_frame(final_df.to_csv, "uploaded_group_1_simulated_cosy_spectra_without_peak_shift.csv")


@app.callback(
    Output("page-3-upload-group-no-peak-shift-group-2-download-spectra-csv", "data"),
    Input("page-3-upload-group-no-peak-shift-group-2-btn-csv", "n_clicks"),
    State("page-3-upload-group-no-peak-shift-group-2-final-dict", "data"),
    prevent_initial_call=True,
)
def page_3_download_upload_group_2_simulated_spectra(n_clicks, final_replicate_dict):
    final_df = page_3_construct_2d_df_without_peak_shift(final_replicate_dict)
    return dcc.send_data_frame(final_df.to_csv, "uploaded_group_2_simulated_cosy_spectra_without_peak_shift.csv")


# -------------- tab 3 callbacks: simulate upload group spectra with peak shift --------------
@app.callback(
    [
        Output("page-3-upload-group-1-ph-same-div", "style"),
        Output("page-3-upload-group-1-ph-not-same-div", "style"),
    ],
    Input("page-3-upload-group-1-ph-same-flag", "value")
)
def page_3_show_upload_group_1_ph_select_div(ph_same_flag):
    if ph_same_flag == "true":
        return {"display": "block"}, {"display": "none"}
    elif ph_same_flag == "false":
        return {"display": "none"}, {"display": "block"}
    else:
        return {"display": "none"}, {"display": "none"}


@app.callback(
    [
        Output("page-3-upload-group-2-ph-same-div", "style"),
        Output("page-3-upload-group-2-ph-not-same-div", "style"),
    ],
    Input("page-3-upload-group-2-ph-same-flag", "value")
)
def page_3_show_upload_group_2_ph_select_div(ph_same_flag):
    if ph_same_flag == "true":
        return {"display": "block"}, {"display": "none"}
    elif ph_same_flag == "false":
        return {"display": "none"}, {"display": "block"}
    else:
        return {"display": "none"}, {"display": "none"}


@app.callback(
    [
        Output("page-3-upload-group-cons-ph-table-div", "style"),
        Output("page-3-upload-group-cons-ph-table", "data"),
        Output("page-3-upload-group-cons-ph-table", "columns"),
        Output("page-3-upload-group-cons-ph-table", "style_data_conditional"),
        Output("page-3-upload-group-peak-shift-part-2-div", "style"),
        Output("page-3-upload-group-peak-shift-spectra-div", "style")
    ],
    Input("page-3-tab-3-upload-group-submit-part-1-btn", "n_clicks"),
    [
        State("page-3-upload-discrete-cons-table", "data"),
        State("page-3-upload-discrete-cons-table", "columns"),
        State("page-3-upload-group-1-ph-same-flag", "value"),
        State("page-3-upload-group-1-same-ph", "value"),
        State("page-3-upload-group-1-ph-mean", "value"),
        State("page-3-upload-group-1-ph-std", "value"),
        State("page-3-upload-group-2-ph-same-flag", "value"),
        State("page-3-upload-group-2-same-ph", "value"),
        State("page-3-upload-group-2-ph-mean", "value"),
        State("page-3-upload-group-2-ph-std", "value")
    ]
)
def page_3_update_upload_group_cons_ph_table(n_clicks, cons_table_data, cons_table_columns, ph_same_flag_1, same_ph_1,
                               ph_mean_1, ph_std_1, ph_same_flag_2, same_ph_2, ph_mean_2, ph_std_2):
    num_repli_1 = sum(1 for i in cons_table_data[0].keys() if i[0] == "1")
    num_repli_2 = sum(1 for i in cons_table_data[0].keys() if i[0] == "2")
    if n_clicks == 0:
        return {"display": "none"}, [{}], [], [{}], {"display": "none"}, {"display": "none"}
    else:
        ph_list_1 = page_3_sample_ph_list(num_repli_1, ph_same_flag_1, same_ph_1, ph_mean_1, ph_std_1)
        ph_list_2 = page_3_sample_ph_list(num_repli_2, ph_same_flag_2, same_ph_2, ph_mean_2, ph_std_2)
        all_ph_list = np.concatenate((["pH", "/"], ph_list_1, ph_list_2))
        name_list = ["meta_name", "hmdb_id"] + ["1_replicate_" + str(i + 1) for i in range(num_repli_1)] + \
                    ["2_replicate_" + str(i + 1) for i in range(num_repli_2)]
        temp_dict = dict(zip(name_list, all_ph_list))
        table_data = [temp_dict] + cons_table_data
        style_data_conditional = [{"if": {"row_index": 0}, "color": "green", "fontWeight": "bold"}]
        return {"display": "block"}, table_data, cons_table_columns, style_data_conditional, {"display": "block"}, {"display": "block"}


@app.callback(
    Output("page-3-upload-group-peak-shift-processed-spectra-data-dict", "data"),
    Input("page-3-upload-group-peak-shift-confirm-param-btn", "n_clicks"),
    [
        State("page-3-upload-group-peak-shift-calibration-range", "value"),
        State("page-3-upload-group-peak-shift-water-range", "value"),
        State("page-3-upload-group-peak-shift-filter-thres", "value"),
        State("page-3-upload-group-peak-shift-smooth-thres-m", "value"),
        State("page-3-upload-group-peak-shift-smooth-thres-n", "value"),
        State("page-3-upload-group-cons-ph-table", "data"),
    ]
)
def page_3_upload_group_peak_shift_update_group_process_data(n_clicks, cal_range, water_range, filter_thres,
                                                      smooth_thres_m, smooth_thres_n, table_rows):
    if n_clicks == 0:
        return None
    else:
        mixture_list = list(map(lambda d: d["meta_name"], table_rows))
        mixture_list.remove('pH')
        temp_raw_data_dict = {name: match_data_dict[name] for name in mixture_list}

        # preprocess 2d data
        removed_data_dict = remove_water_calibration(temp_raw_data_dict, x_scale, y_scale, water_range, cal_range)
        filtered_data_dict = filter_noise(removed_data_dict, filter_thres)
        smooth_data_dict = smooth_data(filtered_data_dict, smooth_thres_m, smooth_thres_n)
        norm_data_dict = normalize_data(smooth_data_dict)
        return norm_data_dict


@app.callback(
    [
        Output("page-3-upload-group-peak-shift-select-mix", "options"),
        Output("page-3-upload-group-peak-shift-select-mix", "value"),
    ],
    Input("page-3-upload-group-peak-shift-confirm-param-btn", "n_clicks"),
    State("page-3-db-names-hmdb-ids-dict", "data"),
)
def page_3_upload_group_peak_shift_update_all_fig_dropdown(n_clicks, final_name_hmdb_id_dict):
    if n_clicks == 0:
        return [], None
    else:
        options = list(final_name_hmdb_id_dict.keys())
        value = options[0]
    return options, value,


@app.callback(
    Output("page-3-upload-group-peak-shift-all-meta-fig", "figure"),
    [
        Input("page-3-upload-group-peak-shift-confirm-param-btn", "n_clicks"),
        Input("page-3-upload-group-peak-shift-processed-spectra-data-dict", "data"),
        Input("page-3-upload-group-peak-shift-select-mix", "value"),
    ]
)
def page_3_upload_group_peak_shift_update_group_results_all_fig(n_clicks, norm_data_dict, select_meta):
    if n_clicks == 0 or norm_data_dict is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        fig = plot_2d_cosy_spectra(select_meta, norm_data_dict, x_scale, y_scale)
        return fig


@app.callback(
    [
        Output("page-3-upload-group-peak-shift-group-1-final-dict", "data"),
        Output("page-3-upload-group-peak-shift-group-2-final-dict", "data"),
    ],
    [
        Input("page-3-upload-group-peak-shift-confirm-param-btn", "n_clicks"),
        Input("page-3-upload-group-peak-shift-processed-spectra-data-dict", "data"),
    ],
    [
        State("page-3-db-names-hmdb-ids-dict", "data"),
        State("page-3-upload-group-peak-shift-snr-noise", "value"),
        State("page-3-upload-group-cons-ph-table", "data"),
    ]
)
def page_3_upload_group_peak_shift_update_groups_repli_data(n_clicks, norm_data_dict, final_names_hmdb_id_dict, snr, cons_ph_table_data):
    if n_clicks == 0 or norm_data_dict is None:
        return None, None
    else:
        mixture_dict = {item['meta_name']: item['hmdb_id'] for item in cons_ph_table_data}
        del mixture_dict["pH"]
        mixture_list = list(mixture_dict.keys())

        mixture_pka_dict = page_3_get_mixture_pka_dict(mixture_list, final_names_hmdb_id_dict)
        group_repli_ph_dict = list(filter(lambda d: d["meta_name"] == 'pH', cons_ph_table_data))[0]

        final_mix_data_dict_1 = get_mixture_data_for_all_replicates("1", group_repli_ph_dict, mixture_dict,
                                                                    norm_data_dict, mixture_pka_dict, x_scale,
                                                                    cons_ph_table_data, protons_df, snr)
        final_mix_data_dict_2 = get_mixture_data_for_all_replicates("2", group_repli_ph_dict, mixture_dict,
                                                                    norm_data_dict, mixture_pka_dict, x_scale,
                                                                    cons_ph_table_data, protons_df, snr)
        return final_mix_data_dict_1, final_mix_data_dict_2



@app.callback(
    [
        Output("page-3-upload-group-peak-shift-select-repli-1", "options"),
        Output("page-3-upload-group-peak-shift-select-repli-1", "value"),
        Output("page-3-upload-group-peak-shift-select-repli-2", "options"),
        Output("page-3-upload-group-peak-shift-select-repli-2", "value"),
    ],
    [
        Input("page-3-upload-group-peak-shift-confirm-param-btn", "n_clicks"),
        Input("page-3-upload-group-peak-shift-group-1-final-dict", "data"),
        Input("page-3-upload-group-peak-shift-group-2-final-dict", "data")
    ]
)
def page_3_upload_group_peak_shift_update_group_replicate_selection(n_clicks, final_data_dict_1, final_data_dict_2):
    if n_clicks == 0 or final_data_dict_1 is None or final_data_dict_2 is None:
        return [], None, [], None
    else:
        options_1 = [{'label': i, 'value': i} for i in final_data_dict_1.keys()]
        value_1 = "replicate_1"
        options_2 = [{'label': i, 'value': i} for i in final_data_dict_2.keys()]
        value_2 = "replicate_1"
        return options_1, value_1, options_2, value_2


@app.callback(
    Output("page-3-upload-group-peak-shift-group-1-spectra-fig", "figure"),
    [
        Input("page-3-upload-group-peak-shift-confirm-param-btn", "n_clicks"),
        Input("page-3-upload-group-peak-shift-group-1-final-dict", "data"),
        Input("page-3-upload-group-peak-shift-select-repli-1", "value"),
    ]
)
def page_3_upload_group_peak_shift_update_group_1_figure(n_clicks, final_data_dict_1, select_repli_1):
    if n_clicks == 0 or final_data_dict_1 is None or select_repli_1 is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        return plot_cosy_repli_spectra_with_ph(final_data_dict_1, select_repli_1, x_scale, y_scale)


@app.callback(
    Output("page-3-upload-group-peak-shift-group-2-spectra-fig", "figure"),
    [
        Input("page-3-upload-group-peak-shift-confirm-param-btn", "n_clicks"),
        Input("page-3-upload-group-peak-shift-group-2-final-dict", "data"),
        Input("page-3-upload-group-peak-shift-select-repli-2", "value"),
    ]
)
def page_3_upload_group_peak_shift_update_group_2_figure(n_clicks, final_data_dict_2, select_repli_2):
    if n_clicks == 0 or final_data_dict_2 is None or select_repli_2 is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        return plot_cosy_repli_spectra_with_ph(final_data_dict_2, select_repli_2, x_scale, y_scale)


@app.callback(
    Output("page-3-upload-group-peak-shift-group-1-download-spectra-csv", "data"),
    Input("page-3-upload-group-peak-shift-group-1-btn-csv", "n_clicks"),
    State("page-3-upload-group-peak-shift-group-1-final-dict", "data"),
    prevent_initial_call=True,
)
def page_3_upload_group_peak_shift_download_group_1_simulated_spectra(n_clicks, final_replicate_dict):
    final_df = page_3_construct_2d_df_with_peak_shift(final_replicate_dict)
    return dcc.send_data_frame(final_df.to_csv, "uploaded_group_1_simulated_cosy_spectra_with_peak_shift.csv")


@app.callback(
    Output("page-3-upload-group-peak-shift-group-2-download-spectra-csv", "data"),
    Input("page-3-upload-group-peak-shift-group-2-btn-csv", "n_clicks"),
    State("page-3-upload-group-peak-shift-group-2-final-dict", "data"),
    prevent_initial_call=True,
)
def page_3_upload_group_peak_shift_download_group_2_simulated_spectra(n_clicks, final_replicate_dict):
    final_df = page_3_construct_2d_df_with_peak_shift(final_replicate_dict)
    return dcc.send_data_frame(final_df.to_csv, "uploaded_group_2_simulated_cosy_spectra_with_peak_shift.csv")


# -------------- tab 3 callbacks: simulate upload continuous spectra no peak shift --------------
@app.callback(
    Output("page-3-upload-conti-processed-spectra-data-dict", "data"),
    Input("page-3-upload-conti-confirm-param-btn", "n_clicks"),
    [
        State("page-3-upload-conti-calibration-range", "value"),
        State("page-3-upload-conti-water-range", "value"),
        State("page-3-upload-conti-filter-thres", "value"),
        State("page-3-upload-conti-smooth-thres-m", "value"),
        State("page-3-upload-conti-smooth-thres-n", "value"),
        State("page-3-upload-conti-cons-table", "data"),
    ]
)
def page_3_update_upload_continuous_process_data(n_clicks, cal_range, water_range, filter_thres, smooth_thres_m, smooth_thres_n, table_rows):
    if n_clicks == 0:
        return None
    else:
        mixture_list = list(map(lambda d: d["meta_name"], table_rows))
        mixture_list.remove("Y")
        temp_raw_data_dict = {name: match_data_dict[name] for name in mixture_list}

        # preprocess 2d data
        removed_data_dict = remove_water_calibration(temp_raw_data_dict, x_scale, y_scale, water_range, cal_range)
        filtered_data_dict = filter_noise(removed_data_dict, filter_thres)
        smooth_data_dict = smooth_data(filtered_data_dict, smooth_thres_m, smooth_thres_n)
        norm_data_dict = normalize_data(smooth_data_dict)
        return norm_data_dict


@app.callback(
    [
        Output("page-3-upload-conti-no-peak-shift-select-mix", "options"),
        Output("page-3-upload-conti-no-peak-shift-select-mix", "value"),
    ],
    Input("page-3-upload-conti-confirm-param-btn", "n_clicks"),
    State("page-3-db-names-hmdb-ids-dict", "data"),
)
def page_3_update_upload_continuous_all_fig_dropdown(n_clicks, final_name_hmdb_id_dict):
    if n_clicks == 0:
        return [], None
    else:
        options = list(final_name_hmdb_id_dict.keys())
        value = options[0]
    return options, value


@app.callback(
    Output("page-3-upload-conti-no-peak-shift-all-meta-fig", "figure"),
    [
        Input("page-3-upload-conti-confirm-param-btn", "n_clicks"),
        Input("page-3-upload-conti-processed-spectra-data-dict", "data"),
        Input("page-3-upload-conti-no-peak-shift-select-mix", "value"),
    ]
)
def page_3_update_upload_continuous_results_all_fig(n_clicks, norm_data_dict, select_meta):
    if n_clicks == 0 or norm_data_dict is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        fig = plot_2d_cosy_spectra(select_meta, norm_data_dict, x_scale, y_scale)
        return fig


@app.callback(
    [
        Output("page-3-upload-conti-no-peak-shift-final-replicate-dict", "data"),
        Output("page-3-upload-conti-no-peak-shift-select-replicate", "options"),
        Output("page-3-upload-conti-no-peak-shift-select-replicate", "value")
    ],
    [
        Input("page-3-upload-conti-confirm-param-btn", "n_clicks"),
        Input("page-3-upload-conti-processed-spectra-data-dict", "data")
    ],
    [
        State("page-3-upload-conti-snr-noise", "value"),
        State("page-3-upload-conti-cons-table", "data"),
    ]
)
def page_3_update_upload_continuous_final_replicate_dict(n_clicks, norm_data_dict, snr, table_rows):
    y_num_repli = len(table_rows[0]) - 1
    if n_clicks == 0 or norm_data_dict is None:
        return None, [], None
    else:
        mixture_dict = {item['meta_name']: item['hmdb_id'] for item in table_rows}
        del mixture_dict["Y"]
        final_mix_data_dict = simulate_continuous_cosy_mixture_for_all_repli(y_num_repli, mixture_dict,
                                                                             norm_data_dict, table_rows,
                                                                             protons_df, snr)
        options = [{"label": i, "value": i} for i in list(final_mix_data_dict.keys())]
        value = list(final_mix_data_dict.keys())[0]
        return final_mix_data_dict, options, value


@app.callback(
    Output("page-3-upload-conti-no-peak-shift-mix-spectra-fig", "figure"),
    Input("page-3-upload-conti-no-peak-shift-select-replicate", "value"),
    [
        State("page-3-upload-conti-confirm-param-btn", "n_clicks"),
        State("page-3-upload-conti-no-peak-shift-final-replicate-dict", "data"),
    ]
)
def page_3_update_upload_continuous_mix_fig(select_replicate, n_clicks, final_data_dict):
    if n_clicks == 0 or select_replicate is None or final_data_dict is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        fig = plot_cosy_repli_spectra(final_data_dict, select_replicate, x_scale, y_scale)
        return fig


@app.callback(
    Output("page-3-upload-conti-no-peak-shift-meta-y-fig", "figure"),
    Input("page-3-upload-conti-no-peak-shift-select-replicate", "value"),
    [
        State("page-3-upload-conti-confirm-param-btn", "n_clicks"),
        State("page-3-upload-conti-cons-table", "data"),
    ]
)
def page_3_update_upload_continuous_meta_y_line_plot(select_replicate, n_clicks, cons_table_rows):
    if n_clicks == 0 or select_replicate is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        fig = page_3_plot_conti_results_line_plot(cons_table_rows, select_replicate)
        return fig


@app.callback(
    Output("page-3-upload-conti-no-peak-shift-download-spectra-csv", "data"),
    Input("page-3-upload-conti-no-peak-shift-btn-csv", "n_clicks"),
    State("page-3-upload-conti-no-peak-shift-final-replicate-dict", "data"),
    prevent_initial_call=True,
)
def page_3_download_upload_continuous_simulated_spectra(n_clicks, final_replicate_dict):
    final_df = page_3_construct_2d_df_without_peak_shift(final_replicate_dict)
    return dcc.send_data_frame(final_df.to_csv, "uploaded_continuous_simulated_cosy_spectra_without_peak_shift.csv")


# -------------- tab 3 callbacks: simulate upload continuous spectra with peak shift --------------
@app.callback(
    [
        Output("page-3-upload-conti-ph-same-div", "style"),
        Output("page-3-upload-conti-ph-not-same-div", "style"),
    ],
    Input("page-3-upload-conti-ph-same-flag", "value")
)
def page_3_show_upload_conti_ph_select_div(ph_same_flag):
    if ph_same_flag == "true":
        return {"display": "block"}, {"display": "none"}
    elif ph_same_flag == "false":
        return {"display": "none"}, {"display": "block"}
    else:
        return {"display": "none"}, {"display": "none"}


@app.callback(
    [
        Output("page-3-upload-conti-cons-ph-table-div", "style"),
        Output("page-3-upload-conti-cons-ph-table", "data"),
        Output("page-3-upload-conti-cons-ph-table", "columns"),
        Output("page-3-upload-conti-cons-ph-table", "style_data_conditional"),
        Output("page-3-upload-conti-peak-shift-part-2-div", "style"),
        Output("page-3-upload-conti-peak-shift-spectra-div", "style")
    ],
    Input("page-3-tab-3-upload-conti-submit-part-1-btn", "n_clicks"),
    [
        State("page-3-upload-conti-cons-table", "data"),
        State("page-3-upload-conti-cons-table", "columns"),
        State("page-3-upload-conti-ph-same-flag", "value"),
        State("page-3-upload-conti-same-ph", "value"),
        State("page-3-upload-conti-ph-mean", "value"),
        State("page-3-upload-conti-ph-std", "value"),
    ]
)
def page_3_update_upload_continuous_cons_ph_table(n_clicks, cons_table_data, cons_table_columns, ph_same_flag, same_ph,
                                    ph_mean, ph_std):
    y_num_repli = len(cons_table_data[0]) - 1
    if n_clicks == 0:
        return {"display": "none"}, [{}], [], [{}], {"display": "none"}, {"display": "none"}
    else:
        ph_list = page_3_sample_ph_list(y_num_repli, ph_same_flag, same_ph, ph_mean, ph_std)
        all_ph_list = np.concatenate((["pH", "/"], ph_list))
        name_list = ["meta_name", "hmdb_id"] + ["replicate_" + str(i + 1) for i in range(y_num_repli)]
        temp_dict = dict(zip(name_list, all_ph_list))
        table_data = [temp_dict] + cons_table_data
        style_data_conditional = [{"if": {"row_index": 0}, "color": "green", "fontWeight": "bold"}]
        return {"display": "block"}, table_data, cons_table_columns, style_data_conditional, {"display": "block"}, {
            "display": "block"}


@app.callback(
    Output("page-3-upload-conti-peak-shift-processed-spectra-data-dict", "data"),
    Input("page-3-upload-conti-peak-shift-confirm-param-btn", "n_clicks"),
    [
        State("page-3-upload-conti-peak-shift-calibration-range", "value"),
        State("page-3-upload-conti-peak-shift-water-range", "value"),
        State("page-3-upload-conti-peak-shift-filter-thres", "value"),
        State("page-3-upload-conti-peak-shift-smooth-thres-m", "value"),
        State("page-3-upload-conti-peak-shift-smooth-thres-n", "value"),
        State("page-3-upload-conti-cons-ph-table", "data"),
    ]
)
def page_3_upload_continuous_peak_shift_update_processed_data(n_clicks, cal_range, water_range, filter_thres,
                                                      smooth_thres_m, smooth_thres_n, table_rows):
    if n_clicks == 0:
        return None
    else:
        mixture_list = list(map(lambda d: d["meta_name"], table_rows))
        mixture_list.remove("Y")
        mixture_list.remove('pH')
        temp_raw_data_dict = {name: match_data_dict[name] for name in mixture_list}
        # preprocess 2d data
        removed_data_dict = remove_water_calibration(temp_raw_data_dict, x_scale, y_scale, water_range, cal_range)
        filtered_data_dict = filter_noise(removed_data_dict, filter_thres)
        smooth_data_dict = smooth_data(filtered_data_dict, smooth_thres_m, smooth_thres_n)
        norm_data_dict = normalize_data(smooth_data_dict)
        return norm_data_dict


@app.callback(
    [
        Output("page-3-upload-conti-peak-shift-select-mix", "options"),
        Output("page-3-upload-conti-peak-shift-select-mix", "value"),
    ],
    Input("page-3-upload-conti-peak-shift-confirm-param-btn", "n_clicks"),
    State("page-3-db-names-hmdb-ids-dict", "data"),
)
def page_3_upload_continuous_peak_shift_update_all_fig_dropdown(n_clicks, final_name_hmdb_id_dict):
    if n_clicks == 0:
        return [], None
    else:
        options = list(final_name_hmdb_id_dict.keys())
        value = options[0]
    return options, value


@app.callback(
    Output("page-3-upload-conti-peak-shift-all-meta-fig", "figure"),
    [
        Input("page-3-upload-conti-peak-shift-confirm-param-btn", "n_clicks"),
        Input("page-3-upload-conti-peak-shift-processed-spectra-data-dict", "data"),
        Input("page-3-upload-conti-peak-shift-select-mix", "value"),
    ]
)
def page_3_upload_continuous_peak_shift_update_group_results_all_fig(n_clicks, norm_data_dict, select_meta):
    if n_clicks == 0 or norm_data_dict is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        fig = plot_2d_cosy_spectra(select_meta, norm_data_dict, x_scale, y_scale)
        return fig


@app.callback(
    Output("page-3-upload-conti-peak-shift-final-replicate-dict", "data"),
    [
        Input("page-3-upload-conti-peak-shift-confirm-param-btn", "n_clicks"),
        Input("page-3-upload-conti-peak-shift-processed-spectra-data-dict", "data"),
    ],
    [
        State("page-3-db-names-hmdb-ids-dict", "data"),
        State("page-3-upload-conti-peak-shift-snr-noise", "value"),
        State("page-3-upload-conti-cons-ph-table", "data"),
    ]
)
def page_3_upload_continous_peak_shift_update_groups_repli_data(n_clicks, norm_data_dict, final_names_hmdb_id_dict, snr, cons_ph_table_data):
    if n_clicks == 0 or norm_data_dict is None:
        return None, None
    else:
        mixture_dict = {item['meta_name']: item['hmdb_id'] for item in cons_ph_table_data}
        del mixture_dict["pH"]
        del mixture_dict["Y"]
        mixture_list = list(mixture_dict.keys())

        mixture_pka_dict = page_3_get_mixture_pka_dict(mixture_list, final_names_hmdb_id_dict)
        conti_repli_ph_dict = list(filter(lambda d: d["meta_name"] == 'pH', cons_ph_table_data))[0]

        final_mix_data_dict = conti_get_mixture_data_for_all_replicates(conti_repli_ph_dict, mixture_dict,
                                                                        norm_data_dict, mixture_pka_dict, x_scale,
                                                                        cons_ph_table_data, protons_df, snr)
        return final_mix_data_dict


@app.callback(
    [
        Output("page-3-upload-conti-peak-shift-select-replicate", "options"),
        Output("page-3-upload-conti-peak-shift-select-replicate", "value"),
    ],
    [
        Input("page-3-upload-conti-peak-shift-confirm-param-btn", "n_clicks"),
        Input("page-3-upload-conti-peak-shift-final-replicate-dict", "data"),
    ]
)
def page_3_upload_group_peak_shift_update_group_replicate_selection(n_clicks, final_data_dict):
    if n_clicks == 0 or final_data_dict is None:
        return [], None
    else:
        options = [{'label': i, 'value': i} for i in final_data_dict.keys()]
        value = "replicate_1"
        return options, value


@app.callback(
    Output("page-3-upload-conti-peak-shift-mix-spectra-fig", "figure"),
    [
        Input("page-3-upload-conti-peak-shift-confirm-param-btn", "n_clicks"),
        Input("page-3-upload-conti-peak-shift-final-replicate-dict", "data"),
        Input("page-3-upload-conti-peak-shift-select-replicate", "value"),
    ]
)
def page_3_upload_continuous_peak_shift_update_mix_figure(n_clicks, final_data_dict, select_repli):
    if n_clicks == 0 or final_data_dict is None or select_repli is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        return plot_cosy_repli_spectra_with_ph(final_data_dict, select_repli, x_scale, y_scale)


@app.callback(
    Output("page-3-upload-conti-peak-shift-meta-y-fig", "figure"),
    Input("page-3-upload-conti-peak-shift-select-replicate", "value"),
    [
        State("page-3-upload-conti-peak-shift-confirm-param-btn", "n_clicks"),
        State("page-3-upload-conti-cons-table", "data"),
    ]
)
def page_3_update_upload_continuous_peak_shift_meta_y_line_plot(select_replicate, n_clicks, cons_table_rows):
    if n_clicks == 0 or select_replicate is None:
        return {"data": [], "layout": {}, "frames": []}
    else:
        fig = page_3_plot_conti_results_line_plot(cons_table_rows, select_replicate)
        return fig


@app.callback(
    Output("page-3-upload-conti-peak-shift-download-spectra-csv", "data"),
    Input("page-3-upload-conti-peak-shift-btn-csv", "n_clicks"),
    State("page-3-upload-conti-peak-shift-final-replicate-dict", "data"),
    prevent_initial_call=True,
)
def page_3_download_upload_conti_peak_shift_simulated_spectra(n_clicks, final_replicate_dict):
    final_df = page_3_construct_2d_df_with_peak_shift(final_replicate_dict)
    return dcc.send_data_frame(final_df.to_csv, "uploaded_continuous_simulated_cosy_spectra_with_peak_shift.csv")


# ----------------------------------- switch from tab to tab -----------------------------------
@app.callback(
    Output("page-3-whole-tabs", "value"),
    [
        Input("page-3-confirm-meta-btn", "n_clicks"),
        Input("page-3-confirm-group-cons-btn", "n_clicks"),
        Input("page-3-confirm-conti-cons-btn", "n_clicks"),
        Input("page-3-confirm-upload-discrete-cons-btn", "n_clicks"),
        Input("page-3-confirm-upload-conti-cons-btn", "n_clicks"),
    ]
)
def page_3_switch_from_tab_to_tab(n_clicks_1, n_clicks_2, n_clicks_3, n_clicks_4, n_clicks_5):
    ctx = dash.callback_context
    if not ctx.triggered:
        # print("not triggered")
        return "page-3-tab-1"
    else:
        changed_id = [p['prop_id'] for p in ctx.triggered][0]
        # print(changed_id)
        if "page-3-confirm-meta-btn" in changed_id:
            return "page-3-tab-2"
        elif "page-3-confirm-group-cons-btn" in changed_id:
            return "page-3-tab-3"
        elif "page-3-confirm-conti-cons-btn" in changed_id:
            return "page-3-tab-3"
        elif "page-3-confirm-upload-discrete-cons-btn" in changed_id:
            return "page-3-tab-3"
        elif "page-3-confirm-upload-conti-cons-btn" in changed_id:
            return "page-3-tab-3"


# ----------------------------------- reset the buttons -----------------------------------
@app.callback(
    [
        Output("page-3-select-simulation-type", "value"),
        Output("page-3-select-spectra-simulation-type", "value"),
    ],
    Input("page-3-confirm-meta-btn", "n_clicks"),
    prevent_initial_call=True,
)
def page_3_update_select_simulation_btn(tab_1_n_clicks):
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
    if "confirm-meta-btn" in changed_id:
        return None, None


@app.callback(
    [
        Output("page-3-submit-part-1-btn", "n_clicks"),
        Output("page-3-submit-part-2-btn", "n_clicks"),
        Output("page-3-conti-submit-part-1-btn", "n_clicks"),
        Output("page-3-conti-submit-part-2-btn", "n_clicks"),
    ],
    [
        Input("page-3-confirm-meta-btn", "n_clicks"),
        Input("page-3-select-simulation-type", "value"),
    ],
    prevent_initial_call=True,
)
def page_3_update_tab_2_buttons(tab_1_n_clicks, cons_type):
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
    if "confirm-meta-btn" in changed_id:
        return 0, 0, 0, 0
    else:
        if cons_type == "group":
            return 0, 0, 0, 0
        elif cons_type == "continuous":
            return 0, 0, 0, 0


@app.callback(
    [
        Output("page-3-confirm-param-btn", "n_clicks"),
        Output("page-3-conti-confirm-param-btn", "n_clicks"),
        Output("page-3-tab-3-group-submit-part-1-btn", "n_clicks"),
        Output("page-3-group-peak-shift-confirm-param-btn", "n_clicks"),
        Output("page-3-tab-3-conti-submit-part-1-btn", "n_clicks"),
        Output("page-3-conti-peak-shift-confirm-param-btn", "n_clicks"),
    ],
    [
        Input("page-3-confirm-meta-btn", "n_clicks"),
        Input("page-3-select-spectra-simulation-type", "value")
    ],
    prevent_initial_call=True,
)
def page_3_update_tab_3_buttons(tab_1_n_clicks, spectra_type):
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
    if "confirm-meta-btn" in changed_id:
        return 0, 0, 0, 0, 0, 0
    else:
        if spectra_type == "no_peak_shift":
            return 0, 0, 0, 0, 0, 0
        elif spectra_type == "peak_shift":
            return 0, 0, 0, 0, 0, 0


current, peak = tracemalloc.get_traced_memory()
print(f"2D COSY page: Current memory usage is {current / 10**6} MB; Peak was {peak / 10**6} MB")
tracemalloc.stop()
