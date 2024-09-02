from dash import dcc
import dash_bootstrap_components as dbc
from dash import html
import tracemalloc
import pathlib
import base64


# tracemalloc.start()
# base_path = pathlib.Path(__file__).resolve().parents[1]
#
# image_filename = 'my-image.png' # replace with your own image
# encoded_image_1 = base64.b64encode(open(base_path.joinpath("Input/MetAssimulo_2_workflow.png"), 'rb').read())
# encoded_image_2 = base64.b64encode(open(base_path.joinpath("Input/peak_shift_leucine.png"), 'rb').read())
# encoded_image_3 = base64.b64encode(open(base_path.joinpath("Input/continuous_line_chart.png"), 'rb').read())

layout = html.Div(
    [
        html.H3("Welcome to MetAssimulo 2.0", style={"font-weight": "bold", 'color': 'steelblue'}),
        html.Br(),
        dbc.Carousel(
            items=[
                {
                    "key": "1",
                    "src": "/assets/MetAssimulo_2_workflow.jpg",
                    "header": "Workflow of MetAssimulo 2.0",
                    "img_style": {"height": "500px", "width": "500px"}
                    # "caption": "and caption",
                },
                {
                    "key": "2",
                    "src": "/assets/peak_shift_leucine.jpg",
                    "header": "Example of peak shift",
                    "caption": "leucine",
                    "img_style": {"height": "500px", "width": "500px"}
                },
                {
                    "key": "3",
                    "src": "/assets/continuous_line_chart.jpg",
                    "header": "Example of simulating continuous outcome",
                    # "caption": "This slide has a caption only",
                    "img_style": {"height": "500px", "width": "500px"}
                },
            ],
            variant="dark",
            controls=True,
            indicators=True,
            interval=2000,
        )
    ]
)



# current, peak = tracemalloc.get_traced_memory()
# print(f"Current memory usage is {current / 10**6} MB; Peak was {peak / 10**6} MB")
# tracemalloc.stop()
