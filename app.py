import io
import numpy as np
import pandas as pd
from scipy import stats

import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from flask import Flask, render_template, request, Response

from qPCR import qPCR_plot

# Create app
app = Flask(__name__)

# Home route
@app.route("/")
def index():
    ''' Return the homepage '''
    return render_template("index.html")

# Route for data processing and bar graph plot
@app.route("/qPCR", methods=["POST"])
def data_proc_plot():

    # Retrieve ".csv" file(s)
    file_list = request.files.getlist("file-sel")
    # Retrieve names of ref gene and control group
    ref_list = request.form.getlist("file-ref-gene")
    ctrl_list = request.form.getlist("file-ctrl-group")

    # Generate list to store colors
    # https://stackoverflow.com/questions/22408237/named-colors-in-matplotlib
    color_list = list(colors.CSS4_COLORS.values())
    # Dict for the default color for bar plot
    color_dict = {}
    for i in range(len(color_list)):
        color_dict[f'Bar Color{i}'] = color_list[i]
    
    # Default values of variables for bar plot
    params_dict = {
        "params-sort": "alphabet_asc",
        "params-sample-groups": color_dict,
        "params-total-bars": 30,
        "params-graph-title": "qPCR",
        "params-axis-label": "Avg. Rel. Tx/Ctrl",
        "params-break-thold": 10,
        "params-alpha-transparency": 0.5,
        "params-capsize": 4,
        "params-legend": 0
    }

    # Update values of variables for bar plot
    for param in list(params_dict.keys()):
        # For those values have been updated in index.html ...
        if request.form.get(param):
            # If no. of sample groups has been updated ...
            if param == "params-sample-groups":
                # List to store updated bar colors
                bar_colors = request.form.getlist("params-bar-colors")
                # Update bar colors in "params_dict"
                for i in range(len(bar_colors)):
                    params_dict[param][f'Bar Color{i}'] = bar_colors[i]
            # For those other than bar colors ...
            else:
                try:
                    # Update info. in "params_dict"
                    params_dict[param] = float(request.form.get(param))
                except ValueError:
                    params_dict[param] = request.form.get(param)

    # Call "qPCR_plot()" function to generate bar graph and save it as "fig"
    fig = qPCR_plot(
        file_list=file_list, ref_list=ref_list, ctrl_list=ctrl_list, bar_color=params_dict["params-sample-groups"], 
        sort_by=params_dict["params-sort"], thold_hbar_ct=params_dict["params-total-bars"], 
        title=params_dict["params-graph-title"], value_label=params_dict["params-axis-label"], 
        break_thold=params_dict["params-break-thold"], alpha=params_dict["params-alpha-transparency"], 
        capsize=params_dict["params-capsize"], legend_loc=params_dict["params-legend"])

    # Pass figure from Flask directly to HTML through in-memory binary streams
    # https://stackoverflow.com/questions/50728328/python-how-to-show-matplotlib-in-flask
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype="image/png")

# Run app
if __name__ == "__main__":
    app.run(debug=True)