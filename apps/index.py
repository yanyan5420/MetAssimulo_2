import dash
from dash import dcc
from dash import html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State, MATCH, ALL, ClientsideFunction
from dash import dash_table

from apps import home_page, one_d_simulation_page, two_d_simulation_page, two_d_cosy_simulation_page
from app import app


# import sys
# sys.path.insert(1, "/Users/yanyan/Desktop/Version_2")


# the style arguments for the sidebar. We use position:fixed and a fixed width
# SIDEBAR_STYLE = {
#     "position": "fixed",
#     "top": 0,
#     "left": 0,
#     "bottom": 0,
#     "width": "16rem",
#     "padding": "2rem 1rem",
#     "background-color": "#f8f9fa",
# }

# the styles for the main content position it to the right of the sidebar and
# add some padding.
CONTENT_STYLE = {
    "margin-left": "5rem",
    "margin-right": "5rem",
    "padding": "2rem 1rem",
}

# sidebar = html.Div(
#     [
#         html.H5("MetAssimulo 2.0", style={"font-weight": "bold", 'color': 'steelblue'}),
#         html.Hr(),
#         # html.P(
#         #     "A simple sidebar layout with navigation links", className="lead"
#         # ),
#         dbc.Nav(
#             [
#                 # html.I(className="fa fa-shield"),
#                 dbc.NavLink("Home", href="/", active="exact"),
#                 dbc.NavLink("1H Simulation", href="/page-1", active="exact"),
#                 dbc.NavLink("JRes Simulation", href="/page-2", active="exact"),
#             ],
#             vertical=True,
#             pills=True,
#         ),
#     ],
#     style=SIDEBAR_STYLE,
# )

navbar = dbc.Navbar(
    dbc.Container(
        children=[
            html.A(
                dbc.NavbarBrand("MetaAssimulo 2.0", className="ms-2", style={'fontSize': 'large',
                                                                             'fontWeight': 'bold',
                                                                             "padding": "0rem 4rem 0rem"}),
                href="/",
            ),
            dbc.NavbarToggler(id="navbar-toggler2"),
            dbc.Collapse(
                dbc.Nav(
                    [
                        # html.I(className="fa fa-shield"),
                        dbc.NavLink("Home", href="/", ),
                        dbc.NavLink("1H Simulation", href="/page-1", ),
                        dbc.NavLink("JRes Simulation", href="/page-2", ),
                        dbc.NavLink("COSY Simulation", href="/page-3", ),
                        # dbc.NavLink("Combined Simulation", href="/page-4", )
                    ],
# active="exact"
                    # className="ml-auto",
                    # className='navbar-nav',
                    className="g-50 ms-auto flex-nowrap mt-3 mt-md-0",
                    navbar=True,
                    # pills=True,
                    style={'fontSize': 'large', 'fontWeight': 'bold', "padding": "0rem 4rem 0rem"}
                ),
                id="navbar-collapse2",
                navbar=True,
            ),
        ],
        className="ml-auto",
        fluid=True,
    ),
    color="primary",
    dark=True,
    # className="ml-auto",
    # className="mb-5 mb-5",
    style={"height": "5.5rem"}
)


# def toggle_navbar_collapse(n, is_open):
#     if n:
#         return not is_open
#     return is_open
#
#
# for i in [2]:
#     app.callback(
#         Output(f"navbar-collapse{i}", "is_open"),
#         [Input(f"navbar-toggler{i}", "n_clicks")],
#         [State(f"navbar-collapse{i}", "is_open")],
#     )(toggle_navbar_collapse)

content = html.Div(id="page-content", style=CONTENT_STYLE)
app.layout = html.Div([dcc.Location(id="url"), navbar, content])


@app.callback(
    Output("page-content", "children"),
    [Input("url", "pathname")]
)
def render_page_content(pathname):
    if pathname == "/":
        return home_page.layout
    elif pathname == "/page-1":
        return one_d_simulation_page.layout
    elif pathname == "/page-2":
        return two_d_simulation_page.layout
            # "JRes simulation page!!!!!!!"
    elif pathname == "/page-3":
        return two_d_cosy_simulation_page.layout
    # elif pathname == "/page-4":
    #     return combined_simulation_page.layout
    else:
        return html.H1("404: Not found", className="text-danger")
    # If the user tries to reach a different page, return a 404 message
    # return dbc.Jumbotron(
    #     [
    #         html.H1("404: Not found", className="text-danger"),
    #         html.Hr(),
    #         html.P(f"The pathname {pathname} was not recognised..."),
    #     ]
    # )


if __name__ == "__main__":
    app.run_server(debug=True, port="8060")

