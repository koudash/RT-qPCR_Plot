import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import glob

def qPCR_plot(sort_by="target_asc", break_thold=10, s_multi_t = 30, title="qPCR", value_label="Rel. Avg. Tx/Ctrl", \
        alpha=0.5, ecolor="black", capsize=4, legend_loc="best"):

    # ***************************************** #
    #               DATA READ-IN                #
    # ***************************************** # 
    # Identify csv files in "data" folder
    filelist = [filename for filename in glob.glob("./data/*.csv")]

    # DataFrame to hold columns for bar plot
    data = pd.DataFrame(columns=["Sample Name", "Target Name", "Rel. Avg. Tx/Ctrl", "Stdev"])

    # Loop through "filelist"
    for file in filelist:

        # Name of ref gene and control group
        ref = input(f'Enter name of ref gene in {file}: ')
        ctrl = input(f'Enter name of control group in {file}: ')

        # Read data
        raw_data = pd.read_csv(file)

        # Create DataFrame to hold columns for data processing
        df = pd.DataFrame(columns=["Sample Name", "Target Name", "Concat Name", "CT", "Target/Ref", "Rel. Tx/Ctrl", "P"])
        # Append values to columns of "df"
        for col in ["Sample Name", "Target Name", "CT"]:
            df[col] = raw_data[col]
        # Change the name of control group to "Ctrl"
        df.loc[df["Sample Name"] == ctrl, "Sample"] = "Ctrl"
        # Determine "Concat Name" of "df"
        df["Concat Name"] = raw_data["Sample Name"] + "_" + raw_data["Target Name"]
        # Set value of "Undetermined" ct as 50
        df.loc[df["CT"] == "Undetermined", "CT"] = 50
        # Typecast value in "CT" column as float
        df["CT"] = df["CT"].astype(float)

        # List to hold index of outliers in "CT" column
        outlier_index = []

        # Loop through the first repeated ct value of "df"
        # Note that the default parallel repeat for each target name is 3
        for i in range(0, len(df), 3):

            # Save 3 repeated ct values in ascending order in "ct_list"
            ct_list = sorted([df["CT"][i], df["CT"][i + 1], df["CT"][i + 2]])

            # Calculate the distance of extreme values towards their adjacent, respectively
            delta_upper = ct_list[2] - ct_list[1]
            delta_lower = ct_list[1] - ct_list[0]

            # Since there are only 3 parallel ct values, consider the one that are far from the other two as the potential outlier
            if delta_upper > delta_lower:
                # delta ct = 1 (2 folds) is set as the threshold for determining outlier ct value
                if delta_upper > 1:
                    # Append index of the maximum value among 3 ct values to "outlier_index"
                    [outlier_index.append(i + j) for j in range(3) if df["CT"][i + j] == ct_list[2]]
            else:
                if delta_lower > 1:
                # Append index of the minimum value among 3 ct values to "outlier_index"
                    [outlier_index.append(i + j) for j in range(3) if df["CT"][i + j] == ct_list[0]]

        # Drop rows with outlier ct values
        df = df.drop(outlier_index)
        # Reset index of "df"
        df.reset_index(drop=True, inplace=True)

        # Group "df" by "Concat Name" for mean ct values
        ct_mean = df.groupby(["Concat Name"]).mean()
        # Reset index for "ct_mean"
        ct_mean.reset_index(inplace=True)

        # Loop through "Sample Name"
        for s_name in df["Sample Name"].unique():

            # Calculate mean ct value of ref gene for iterated "Sample Name"
            s_name_ct_mean = ct_mean.loc[ct_mean["Concat Name"] == f'{s_name}_{ref}', "CT"].values[0]

            # Calculate ratio of Target/Ref for each repeat of all target genes
            df.loc[df["Sample Name"] == s_name, "Target/Ref"] = pow(2, (s_name_ct_mean - df["CT"]))

        # Typecast values in "Target/Ref" column from non-null object to float
        df["Target/Ref"] = df["Target/Ref"].astype(float)

        # Remove ref gene from "df"
        df = df.loc[df["Target Name"] != ref, :]

        # DataFrame to hold mean Target/Ref ratio of Target Names in control group
        ctrl_mean_df = df.loc[df["Sample Name"] == "Ctrl", :].groupby(["Concat Name"])["Target/Ref"].mean()

        # Loop through "Target Name"
        for t_name in df["Target Name"].unique():

            # Loop through "Sample Name"
            for s_name in df["Sample Name"].unique():

                # Calculate relative ratio of Target/Ref (as compared to that of avg. value from Ctrl) for each repeat of all target genes
                df.loc[df["Concat Name"] == f'{s_name}_{t_name}', "Rel. Tx/Ctrl"] = df["Target/Ref"] / ctrl_mean_df[f'Ctrl_{t_name}']

        # Typecast values in "Rel. Tx/Ctrl" column from non-null object to float
        df["Rel. Tx/Ctrl"] = df["Rel. Tx/Ctrl"].astype(float)

        # DataFrame to hold columns of iterated csv file for bar plot
        file_data = pd.DataFrame(columns=["Sample Name", "Target Name", "Concat Name", "Rel. Avg. Tx/Ctrl", "Stdev", "P"])

        # Append values to columns of "file_data"
        file_data["Concat Name"] = df.groupby(["Concat Name"])["Rel. Tx/Ctrl"].mean().index
        file_data["Rel. Avg. Tx/Ctrl"] = df.groupby(["Concat Name"])["Rel. Tx/Ctrl"].mean().values
        file_data["Stdev"] = df.groupby(["Concat Name"])["Rel. Tx/Ctrl"].std().values

        # Lists to store "Sample Name" and "Target Name"
        s_name_list = []
        t_name_list = []

        # Loop through "file_data" to append sample and target names for each row to their designated lists
        for i in range(len(file_data)):
            s_name_list.append(file_data["Concat Name"][i].split("_")[0])
            t_name_list.append(file_data["Concat Name"][i].split("_")[1])

        # Add values to columns of "Sample Name" and "Target Name" from "file_data"
        file_data["Sample Name"] = s_name_list
        file_data["Target Name"] = t_name_list

        # Loop through target names of "df"
        for t_name in df["Target Name"].unique():

            # DataFrame to hold rows of "t_name"
            t_name_df = df.loc[df["Target Name"] == t_name, :]
            # Reset index for "t_name_df"
            t_name_df.reset_index(drop=True, inplace=True)

            # Dict to hold possibility of ttest for "Rel. Tx/Ctrl"
            p_df = {}            

            # Loop through "t_name_df"
            for i in range(len(t_name_df)):

                # Variable for "Sample Name"
                s = t_name_df["Sample Name"][i]

                # Check if "p_df[s]" exists
                if not s in list(p_df.keys()):
                    p_df[s] = []
                # Append "Rel. Tx/Ctrl" value to "p_df[s]"
                p_df[s].append(t_name_df["Rel. Tx/Ctrl"][i])

            # Loop through keys of "p_df"
            for s in list(p_df.keys()):

                # Perform t-test on group(s) other than "Ctrl"
                if s != "Ctrl":

                    # Consider default variances for values between two series are equal
                    equal_var = True

                    # Calculate p_cdf for f-test (cumulative distribution function)
                    p_cdf = stats.f.cdf(np.var(p_df[s]) / np.var(p_df["Ctrl"]), len(p_df[s]) - 1, len(p_df["Ctrl"]) - 1)

                    # Make sure variance1 is greater than variance2 to determine "p_cdf"
                    if np.var(p_df["Ctrl"]) > np.var(p_df[s]):
                        p_cdf = 1 - p_cdf
                    
                    # Check whether p-cdf is less than alpha which is by default 0.05 and if true ...
                    if p_cdf < 0.05:
                        # Reject the null hypothesis and accept that the variances between two series are not equal
                        equal_var = False                  

                    # Calculate p value for t-test
                    p = stats.ttest_ind(p_df[s], p_df["Ctrl"], equal_var=equal_var)[1]

                    # Variable for asterisk to be shown on bar plot
                    asterisk = ""
                    # Determine "asterisk" based on "p"
                    if p < 0.001:
                        asterisk = "***"
                    elif p < 0.01:
                        asterisk = "**"
                    elif p < 0.05:
                        asterisk = "*"
                    # Append "asterisk" to designated "P" column of "file_data"
                    file_data.loc[file_data["Concat Name"] == f'{s}_{t_name}', "P"] = asterisk
                
                else:
                    # Append "" to "P" column of "Ctrl" in "file_data"
                    file_data.loc[file_data["Concat Name"] == f'{s}_{t_name}', "P"] = ""

        # Drop "Concat Name" column from "file_data"
        file_data = file_data.drop(columns=["Concat Name"])
        
        # Add "file_data" to "data"
        data = pd.concat([data, file_data], sort=True, ignore_index=True)

    # ************************************** #
    #               BAR PLOT                 #
    # ************************************** # 
    # DataFrame is sorted by "Target Name" in descending order by default after all csv files have been processed
    data = data.sort_values(by=["Target Name"])

    # Dict to store DataFrame for different sample names
    s_df = {}

    # List to store legend for plots of different sample names
    s_plots_legend = []

    # Dict for upper limit on axis of value
    upper_limit = {}

    # Boolean for bar break
    bar_break = False

    # Dicts for break values on axis of value
    upper_break = {}
    lower_break = {}

    # List to hold unique sample names (control + treatment groups) of "data"
    s_names = list(data["Sample Name"].unique())

    # Make sure "Ctrl" is the first element in "s_names"
    s_names.remove("Ctrl")
    s_names.insert(0, "Ctrl")

    # Divide "data" into different DataFrames according to unique elements in "s_names" and reset index for each "s_name"
    for s_name in s_names:
        s_df[s_name] = data.loc[data["Sample Name"] == s_name, :].reset_index(drop=True)

    # Variable for total bars
    total_bars = len(s_names) * len(s_df["Ctrl"]["Target Name"])

    # Basic position on axis of target names for bar graph
    basic_pos = range(len(s_df["Ctrl"]))
    # Width of the bars
    width = 1 / (len(s_names) + 1)
    # List to hold coords of target names for bar graph
    value_ticks = []
    # Offside towards the basic position (one bar width plus half of that for total bars)
    [value_ticks.append(ele + (1 + len(s_names) / 2) * width) for ele in basic_pos]

    # Check the order of bars for the plot
    if len(s_names) == 2 and sort_by.split("_")[0] == "folds":

        # Boolean for ascending parameter
        sort_asc = False
        if sort_by.split("_")[1] == "asc":
            sort_asc = True

        # Sort "s_df" of tx group by sort_by variable and reset index
        s_df[s_names[1]] = s_df[s_names[1]].sort_values(by=["Rel. Avg. Tx/Ctrl"], ascending=sort_asc).reset_index(drop=True)
        # Sort "s_df" of "Ctrl" group accordingly
        # https://stackoverflow.com/questions/45576800/how-to-sort-dataframe-based-on-a-column-in-another-dataframe-in-pandas
        s_df["Ctrl"] = s_df["Ctrl"].set_index("Target Name")
        s_df["Ctrl"] = s_df["Ctrl"].reindex(index=s_df[s_names[1]]["Target Name"])
        s_df["Ctrl"] = s_df["Ctrl"].reset_index()

    # <---------- DETERMINATION OF UPPER LIMIT, UPPER BREAK, AND LOWER BREAK FOR EACH SAMPLE NAME DATAFREAME ----------> #
    # Loop through "s_names"
    for i in range(len(s_names)):

        # Sort iterated "s_df" by "Rel. Avg. Tx/Ctrl" in ascending order, reset index, and temporarily store in "s_df_iter" variable 
        s_df_iter = s_df[s_names[i]].sort_values(by=["Rel. Avg. Tx/Ctrl"], ascending=True).reset_index(drop=True)

        # Determine upper limit on axis of value for iterated "s_df"
        upper_limit[s_names[i]] = round(s_df_iter.iloc[-1]["Rel. Avg. Tx/Ctrl"] * 1.5)

        # Loop through "Target Name" of "s_df_iter"
        for j in range(len(s_df_iter) - 1):

            # Check if break will be applied on axis of value
            if s_df_iter["Rel. Avg. Tx/Ctrl"][j + 1] / s_df_iter["Rel. Avg. Tx/Ctrl"][j] > break_thold:

                # Change "bar_break" to True
                bar_break = True

                # Determine lower break value on axis of value for bar of iterated "Target Name"
                lower_break[s_names[i]] = round(s_df_iter["Rel. Avg. Tx/Ctrl"][j] * 1.5)

                # Loop through the rest of "Target Name", including j
                for k in range(j, len(s_df_iter) - 1):

                    # Look for upper break value on axis of value for bar of iterated "Target Name"
                    if s_df_iter["Rel. Avg. Tx/Ctrl"][k + 1] / s_df_iter["Rel. Avg. Tx/Ctrl"][k] > break_thold:

                        # Determine upper break value on axis of value for bar of iterated "Target Name"
                        upper_break[s_names[i]] = round(s_df_iter["Rel. Avg. Tx/Ctrl"][k + 1] * 0.5)

                        # Jump out "for loop"
                        break

    # Determine upper limit on axis of value
    ul = max(upper_limit.values())

    # Check if bar graph will be plotted horizontally
    if total_bars > s_multi_t:

        # Variable to store orientation of bar graph
        orien = "H"
    
        # Check if plot bar graph with break
        if bar_break:

            # Bar plots
            fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=[total_bars / 7.5, total_bars / 3], dpi=100)

            # Counter for offside width for "pos"
            counter = 0

            # Loop through "s_names" in reverse order
            for i in range(len(s_names))[::-1]:

                # List for "self" parameter of bar plot
                pos = []
                [pos.append(ele + width * (counter + 1)) for ele in basic_pos]

                # Add 1 to "counter"
                counter += 1
            
                # Determine bar plots for different sample names
                bp = ax1.barh(pos, s_df[s_names[i]]["Rel. Avg. Tx/Ctrl"], width, xerr=s_df[s_names[i]]["Stdev"], \
                    align="edge", alpha=alpha, ecolor=ecolor, capsize=capsize)
                ax2.barh(pos, s_df[s_names[i]]["Rel. Avg. Tx/Ctrl"], width, xerr=s_df[s_names[i]]["Stdev"], \
                    align="edge", alpha=alpha, ecolor=ecolor, capsize=capsize)

                # Add asterisk(s) for bars of treatment group(s) if applicable
                if s_names[i] != "Ctrl":

                    # Loop through "s_df" of iterated sample name
                    for j in range(len(s_df[s_names[i]])):

                        # Asterisk(s) for "s" parameter of "Axes.text()"
                        text_s = s_df[s_names[i]]["P"][j]

                        # x-coord for upper cap of iterated error bar
                        text_x = s_df[s_names[i]]["Rel. Avg. Tx/Ctrl"][j] + s_df[s_names[i]]["Stdev"][j]

                        # Check if iterated bar has break and if not ...
                        if text_x < lower_break[s_names[i]]:
                            # Add text label to bar in "ax1"
                            # Set space between upper cap of error bar and text as 1/20 of the width in "ax1"
                            ax1.text(text_x + lower_break[s_names[i]] / 20, pos[j] + width / 8, text_s)
                        # If break exists ...
                        else:
                            # Add text label to bar in "ax2"
                            # Set space between upper cap of error bar and text as 1/20 of the width in "ax2"
                            ax2.text(text_x + (upper_limit[s_names[i]] - upper_break[s_names[i]]) /20, pos[j] + width / 8, text_s)

                # Determine legend for plots of different sample names
                s_plots_legend.append(bp)

            # Determine upper and lower break limits on axis of value
            ub = max(upper_break.values())
            lb = min(lower_break.values())

            # Set limit on axis of value
            ax1.set_xlim(0, lb)
            ax2.set_xlim(ub, ul)
            # Set limit on axis of target names
            ax1.set_ylim(0, len(s_df["Ctrl"]) + width) 

            # Set ticks, ticklables, and lables on axis of value
            ax1.tick_params(axis="x", labeltop="on", top=True)
            ax2.tick_params(axis="x", labeltop="on", top=True)
            ax1.set_xlabel(value_label)

            # Set ticks and ticklabels on axis of target names
            ax2.tick_params(axis="y", left=False)
            ax1.set_yticks(value_ticks)
            ax1.set_yticklabels(list(s_df["Ctrl"]["Target Name"]))

            # Set title of bar graph
            ax2.set_title(title)

            # Set legend of bar graph
            ax2.legend(s_plots_legend, s_names[::-1], loc=legend_loc) 

            # Set spines of the bar graph for break area as invisible
            ax1.spines["right"].set_visible(False)
            ax2.spines["left"].set_visible(False)

            # How big to make the diagonal lines in axes coordinates
            d= 0.005
            # Plot "/" on right limit of axis of target name
            kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
            ax1.plot((1-d, 1+d), (-d, +d), **kwargs)      # Bottom-left diagonal
            ax1.plot((1-d, 1+d), (1-d, 1+d), **kwargs)    # Top-left diagonal
            # Plot "/" on left limit of axis of target name
            kwargs.update(transform=ax2.transAxes)  
            ax2.plot((-d, +d), (1-d, 1+d), **kwargs)      # Top-right diagonal
            ax2.plot((-d, +d), (-d, +d), **kwargs)        # Bottom-right diagonal
            
            # Show the bar graph
            plt.show()

            # Save bar graph
            fig.savefig(f'./figures/{title}_sortby-{sort_by.split("_")[0]}.{sort_by.split("_")[1]}_orien-{orien}.png', \
                bbox_inches="tight")

        # No break in horizontal bar graph
        else:

            # Bar plots
            fig, ax = plt.subplots(figsize=[total_bars / 7.5, total_bars / 3], dpi=100)

            # Counter for offside width for "pos"
            counter = 0

            # Loop through "s_names" in reverse order
            for i in range(len(s_names))[::-1]:

                # List for "self" parameter of bar plot
                pos = []
                [pos.append(ele + width * (counter + 1)) for ele in basic_pos]

                # Add 1 to "counter"
                counter += 1
            
                # Determine bar plots for different sample names
                bp = ax.barh(pos, s_df[s_names[i]]["Rel. Avg. Tx/Ctrl"], width, xerr=s_df[s_names[i]]["Stdev"], \
                    align="edge", alpha=alpha, ecolor=ecolor, capsize=capsize)

                # Add asterisk(s) for bars of treatment group(s) if applicable
                if s_names[i] != "Ctrl":

                    # Loop through "s_df" of iterated sample name
                    for j in range(len(s_df[s_names[i]])):

                        # Asterisk(s) for "s" parameter of "Axes.text()"
                        text_s = s_df[s_names[i]]["P"][j]

                        # x-coord for upper cap of iterated error bar
                        text_x = s_df[s_names[i]]["Rel. Avg. Tx/Ctrl"][j] + s_df[s_names[i]]["Stdev"][j]

                        # Add text label to bar in "ax1"
                        # Set space between upper cap of error bar and text as 1/20 of the width in "ax1"
                        ax.text(text_x + upper_limit[s_names[i]] / 40, pos[j] + width / 8, text_s)

                # Determine legend for plots of different sample names
                s_plots_legend.append(bp)

            # Set limit on axis of value
            ax.set_xlim(0, ul)
            # Set limit on axis of target names
            ax.set_ylim(0, len(s_df["Ctrl"]) + width) 

            # Set ticks, ticklables, and lables on axis of value
            ax.tick_params(axis="x", labeltop="on", top=True)
            ax.set_xlabel(value_label)

            # Set ticks and ticklabels on axis of target names
            ax.set_yticks(value_ticks)
            ax.set_yticklabels(list(s_df["Ctrl"]["Target Name"]))

            # Set title of bar graph
            ax.set_title(title)

            # Set legend of bar graph
            ax.legend(s_plots_legend, s_names[::-1], loc=legend_loc) 
          
            # Show the bar graph
            plt.show()

            # Save bar graph
            fig.savefig(f'./figures/{title}_sortby-{sort_by.split("_")[0]}.{sort_by.split("_")[1]}_orien-{orien}.png', \
                bbox_inches="tight")

    # Bar graph will be plotted vertically
    else:
        # Variable to store orientation of bar graph
        orien = "V"

        # Check if plot bar graph with break
        if bar_break:

            # Bar plots
            fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=[total_bars / 3, total_bars / 7.5], dpi=100)

            # Counter for offside width for "pos"
            counter = 0

            # Loop through "s_names" in reverse order
            for i in range(len(s_names)):

                # List for "self" parameter of bar plot
                pos = []
                [pos.append(ele + width * (counter + 1)) for ele in basic_pos]

                # Add 1 to "counter"
                counter += 1
            
                # Determine bar plots for different sample names
                bp = ax2.bar(pos, s_df[s_names[i]]["Rel. Avg. Tx/Ctrl"], width, yerr=s_df[s_names[i]]["Stdev"], \
                    align="edge", alpha=alpha, ecolor=ecolor, capsize=capsize)
                ax1.bar(pos, s_df[s_names[i]]["Rel. Avg. Tx/Ctrl"], width, yerr=s_df[s_names[i]]["Stdev"], \
                    align="edge", alpha=alpha, ecolor=ecolor, capsize=capsize)

                # Add asterisk(s) for bars of treatment group(s) if applicable
                if s_names[i] != "Ctrl":

                    # Loop through "s_df" of iterated sample name
                    for j in range(len(s_df[s_names[i]])):

                        # Asterisk(s) for "s" parameter of "Axes.text()"
                        text_s = s_df[s_names[i]]["P"][j]

                        # x-coord for upper cap of iterated error bar
                        text_y = s_df[s_names[i]]["Rel. Avg. Tx/Ctrl"][j] + s_df[s_names[i]]["Stdev"][j]

                        # Check if iterated bar has break and if not ...
                        if text_y < lower_break[s_names[i]]:
                            # Add text label to bar in "ax2"
                            # Set space between upper cap of error bar and text as 1/20 of the width in "ax2"
                            ax2.text(pos[j] + width / 2, text_y + lower_break[s_names[i]] / 20, text_s, horizontalalignment='center')
                        # If break exists ...
                        else:
                            # Add text label to bar in "ax1"
                            # Set space between upper cap of error bar and text as 1/20 of the width in "ax1"
                            ax1.text(pos[j] + width / 2, text_y + (upper_limit[s_names[i]] - upper_break[s_names[i]]) /20, text_s, \
                                horizontalalignment='center')

                # Determine legend for plots of different sample names
                s_plots_legend.append(bp)

            # Determine upper and lower break limits on axis of value
            ub = max(upper_break.values())
            lb = min(lower_break.values())

            # Set limit on axis of value
            ax2.set_ylim(0, lb)
            ax1.set_ylim(ub, ul)
            # Set limit on axis of target names
            ax2.set_xlim(0, len(s_df["Ctrl"]) + width) 

            # Set lable on axis of value
            ax1.set_ylabel(value_label)

            # Set ticks and ticklabels on axis of target names
            ax1.tick_params(axis="x", bottom=False)
            ax2.tick_params(axis="x", labelrotation=30)
            ax2.set_xticks(value_ticks)
            ax2.set_xticklabels(list(s_df["Ctrl"]["Target Name"]))

            # Set title of bar graph
            ax1.set_title(title)

            # Set legend of bar graph
            ax1.legend(s_plots_legend, s_names, loc=legend_loc) 

            # Set spines of the bar graph for break area as invisible
            ax2.spines["top"].set_visible(False)
            ax1.spines["bottom"].set_visible(False)

            # How big to make the diagonal lines in axes coordinates
            d= 0.005
            # Plot "/" on upper limit of axis of target name
            kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
            ax1.plot((-d, +d), (-d, +d), **kwargs)      # Top-left diagonal
            ax1.plot((1-d, 1+d), (-d, +d), **kwargs)    # Top-right diagonal
            # Plot "/" on lower limit of axis of target name
            kwargs.update(transform=ax2.transAxes)  
            ax2.plot((-d, +d), (1-d, 1+d), **kwargs)      # Bottom-left diagonal
            ax2.plot((1-d, 1+d), (1-d, 1+d), **kwargs)        # Bottom-right diagonal
            
            # Show the bar graph
            plt.show()

            # Save bar graph
            fig.savefig(f'./figures/{title}_sortby-{sort_by.split("_")[0]}.{sort_by.split("_")[1]}_orien-{orien}.png', \
                bbox_inches="tight")

        # No break in bar graph
        else:
            # Bar plots
            fig, ax = plt.subplots(figsize=[total_bars / 3, total_bars / 7.5], dpi=100)

            # Counter for offside width for "pos"
            counter = 0

            # Loop through "s_names" in reverse order
            for i in range(len(s_names)):

                # List for "self" parameter of bar plot
                pos = []
                [pos.append(ele + width * (counter + 1)) for ele in basic_pos]

                # Add 1 to "counter"
                counter += 1
            
                # Determine bar plots for different sample names
                bp = ax.bar(pos, s_df[s_names[i]]["Rel. Avg. Tx/Ctrl"], width, yerr=s_df[s_names[i]]["Stdev"], \
                    align="edge", alpha=alpha, ecolor=ecolor, capsize=capsize)

                # Add asterisk(s) for bars of treatment group(s) if applicable
                if s_names[i] != "Ctrl":

                    # Loop through "s_df" of iterated sample name
                    for j in range(len(s_df[s_names[i]])):

                        # Asterisk(s) for "s" parameter of "Axes.text()"
                        text_s = s_df[s_names[i]]["P"][j]

                        # x-coord for upper cap of iterated error bar
                        text_y = s_df[s_names[i]]["Rel. Avg. Tx/Ctrl"][j] + s_df[s_names[i]]["Stdev"][j]

                        # Add text label to bar in "ax"
                        # Set space between upper cap of error bar and text as 1/40 of the width in "ax"
                        ax.text(pos[j] + width / 2, text_y + upper_limit[s_names[i]] / 40, text_s, horizontalalignment='center')

                # Determine legend for plots of different sample names
                s_plots_legend.append(bp)

            # Set limit on axis of value
            ax.set_ylim(0, ul)
            # Set limit on axis of target names
            ax.set_xlim(0, len(s_df["Ctrl"]) + width) 

            # Set lables on axis of value
            ax.set_ylabel(value_label)

            # Set ticks and ticklabels on axis of target names
            ax.tick_params(axis="x", labelrotation=30)
            ax.set_xticks(value_ticks)
            ax.set_xticklabels(list(s_df["Ctrl"]["Target Name"]))

            # Set title of bar graph
            ax.set_title(title)

            # Set legend of bar graph
            ax.legend(s_plots_legend, s_names, loc=legend_loc) 

            # Show the bar graph
            plt.show()

            # Save bar graph
            fig.savefig(f'./figures/{title}_sortby-{sort_by.split("_")[0]}.{sort_by.split("_")[1]}_orien-{orien}.png', \
                bbox_inches="tight")