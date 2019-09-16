import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

'''
Encode function named "qPCR_plot()" that takes in info. of qPCR design (reference gene and control groups for 
results from each expt. saved as ".csv" file) as well as parameters for the display of bar graph. Generated figure 
is returned as output of the function.
'''

# Function to read ".csv" file(s) and generate bar graph accordingly
def qPCR_plot(
    file_list, ref_list, ctrl_list, bar_color, 
    ct_repeat, sort_by="alphabet_asc", thold_hbar_ct=30, 
    title="qPCR", value_label="Avg. Rel. Tx/Ctrl", break_thold=10, 
    alpha=0.5, capsize=4, legend_loc="0"):

    # ***************************************** #
    #               DATA READ-IN                #
    # ***************************************** # 
    # DataFrame for bar plot
    data = pd.DataFrame(columns=["Sample Name", "Target Name", "Avg. Rel. Tx/Ctrl", "Stdev"])

    # Iterations to retrieve column info. for "data"
    for i in range(len(file_list)):

        # >>>>>>>>>> Ori. data read in <<<<<<<<<< #
        # Names of ref gene and control group(s) in each ".csv" file
        ref = ref_list[i]
        ctrl = ctrl_list[i]
        # Read data in iterated ".csv" file
        raw_data = pd.read_csv(file_list[i])

        # Create DataFrame to hold columns for data processing
        df = pd.DataFrame(columns=["Sample Name", "Target Name", "Concat Name", "CT", "Target/Ref", "Rel. Tx/Ctrl"])

        # Append values to columns of "df"
        for col in ["Sample Name", "Target Name", "CT"]:
            df[col] = raw_data[col]
        # Unify the name of control group as "Ctrl" in case more than one ctrl group names are encountered
        df.loc[df["Sample Name"] == ctrl, "Sample Name"] = "Ctrl"
        # Determine "Concat Name" of "df"
        df["Concat Name"] = df["Sample Name"] + "_" + df["Target Name"]
        # Set value of "Undetermined" ct as 50
        df.loc[df["CT"] == "Undetermined", "CT"] = 50
        # Typecast value in "CT" column as float
        df["CT"] = df["CT"].astype(float)

        # Check if outlier removal is needed
        if ct_repeat == 3:

            # >>>>>>>>>> Removal of ct outliers <<<<<<<<<< #
            # List to hold index of outliers in "CT" column
            outlier_index = []

            # Note that the default parallel repeat for each target name is 3
            # Loop through the first repeated ct value of "df"
            for i in range(0, len(df), 3):
                # Save 3 repeated ct values in ascending order in "ct_list"
                ct_list = sorted([df["CT"][i], df["CT"][i + 1], df["CT"][i + 2]])
                # Calculate the distance in between each value in "ct_list"
                delta_upper = ct_list[2] - ct_list[1]
                delta_lower = ct_list[1] - ct_list[0]
                # Since there are only 3 parallel ct values, consider the most distant one as potential outlier
                if delta_upper > delta_lower:
                    # delta ct = 1 (2 folds) is set as the threshold for determining outlier ct value
                    if delta_upper > 1:
                        # Append index of the maximum ct value to "outlier_index"
                        [outlier_index.append(i + j) for j in range(3) if df["CT"][i + j] == ct_list[2]]
                else:
                    if delta_lower > 1:
                        # Append index of the minimum ct value to "outlier_index"
                        [outlier_index.append(i + j) for j in range(3) if df["CT"][i + j] == ct_list[0]]

            # Drop rows with outlier ct values
            df = df.drop(outlier_index)
            # Reset index of "df"
            df.reset_index(drop=True, inplace=True)

        # >>>>>>>>>> Column calc in "df" <<<<<<<<<< #
        # ..... 1. "Target/Ref" column ..... #
        # Calculate mean ct values for experimental repeats
        ct_mean = df.groupby(["Concat Name"]).mean()
        # Reset index for "ct_mean"
        ct_mean.reset_index(inplace=True)

        # Loop through sample groups including both ctrl and treatments
        for s_name in df["Sample Name"].unique():
            # Calculate mean ct value of ref gene for iterated "Sample Name"
            ct_mean_ref = ct_mean.loc[ct_mean["Concat Name"] == f'{s_name}_{ref}', "CT"].values[0]
            # Calculate ratio of Target/Ref for each record with iterated "Sample Name" 
            df.loc[df["Sample Name"] == s_name, "Target/Ref"] = pow(2, (ct_mean_ref - df["CT"]))

        # Typecast values in "Target/Ref" column from non-null object to float
        df["Target/Ref"] = df["Target/Ref"].astype(float)
        # Remove rows of ref gene from "df"
        df = df.loc[df["Target Name"] != ref, :]

        # ..... 2. "Rel. Tx/Ctrl" column ..... #
        # DataFrame to hold mean Target/Ref ratio for each "Target Name" in control group
        t_r_mean_ctrl = df.loc[df["Sample Name"] == "Ctrl", :].groupby(["Concat Name"])["Target/Ref"].mean()

        # Loop through "Target Name"
        for t_name in df["Target Name"].unique():
            # Loop through "Sample Name"
            for s_name in df["Sample Name"].unique():
                # Calculate relative ratio of Target/Ref (as compared to that of avg. value from Ctrl) for each repeat of target genes
                df.loc[df["Concat Name"] == f'{s_name}_{t_name}', "Rel. Tx/Ctrl"] = df["Target/Ref"] / t_r_mean_ctrl[f'Ctrl_{t_name}']

        # Typecast values in "Rel. Tx/Ctrl" column from non-null object to float
        df["Rel. Tx/Ctrl"] = df["Rel. Tx/Ctrl"].astype(float)

        # >>>>>>>>>> Column calc in "file_data" for bar plot <<<<<<<<<< #
        # DataFrame to hold columns of iterated ".csv" file for bar plot
        file_data = pd.DataFrame(columns=["Sample Name", "Target Name", "Concat Name", "Avg. Rel. Tx/Ctrl", "Stdev", "P"])

        # Calculate values in columns of "file_data"
        file_data["Concat Name"] = df.groupby(["Concat Name"])["Rel. Tx/Ctrl"].mean().index
        file_data["Avg. Rel. Tx/Ctrl"] = df.groupby(["Concat Name"])["Rel. Tx/Ctrl"].mean().values
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

        # Loop through "Target Name" of "df"
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
                # Create key of "s" in "p_df" if not exists
                if not s in list(p_df.keys()):
                    p_df[s] = []
                # Append "Rel. Tx/Ctrl" value to "p_df[s]"
                p_df[s].append(t_name_df["Rel. Tx/Ctrl"][i])

            # Loop through keys of "p_df"
            for s in list(p_df.keys()):

                # Perform t-test between ctrl and treatment(s) groups
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
    # Processed data is sorted by "Target Name" in descending order by default after all ".csv" files have been processed
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

    # Re-order "s_df" of ctrl and treatment groups if "sort by value" option is available and selected
    if len(s_names) == 2 and sort_by.split("_")[0] == "value":
        # Boolean for ascending parameter
        sort_asc = False
        if sort_by.split("_")[1] == "asc":
            sort_asc = True
        # Sort "s_df" of treatment group by "sort_by" and reset index
        # Note that "Ctrl" is the first value in "s_names"
        s_df[s_names[1]] = s_df[s_names[1]].sort_values(by=["Avg. Rel. Tx/Ctrl"], ascending=sort_asc).reset_index(drop=True)
        # Sort "s_df" of ctrl group accordingly
        # https://stackoverflow.com/questions/45576800/how-to-sort-dataframe-based-on-a-column-in-another-dataframe-in-pandas
        s_df["Ctrl"] = s_df["Ctrl"].set_index("Target Name").reindex(index=s_df[s_names[1]]["Target Name"]).reset_index()

    # Basic position on axis of target names for bar graph
    basic_pos = range(len(s_df["Ctrl"]))
    # Bar width (extra one on denominator represents gap between bunches of bars from different sample groups)
    width = 1 / (len(s_names) + 1)
    # List to hold coords on axis of target names for bar graph
    value_ticks = []
    # Coords for ticks on axis of target names (basic position + offsides)
    # Note that offsides towards basic position equals to gap (1 * width) plus half of that for all bars (middle position of bars)
    [value_ticks.append(num + (1 + len(s_names) / 2) * width) for num in basic_pos]

    # Determine break limits (if applicable) and upper limit on axis of value
    for i in range(len(s_names)):

        # Sort iterated "s_df" by "Avg. Rel. Tx/Ctrl" in ascending order, reset index, and temporarily store in "s_df_iter" variable 
        s_df_iter = s_df[s_names[i]].sort_values(by=["Avg. Rel. Tx/Ctrl"], ascending=True).reset_index(drop=True)

        # Determine upper limit on axis of value for iterated "s_df" (1.2 folds of the maximum "Avg. Rel. Tx/Ctrl" value)
        upper_limit[s_names[i]] = round((s_df_iter.iloc[-1]["Avg. Rel. Tx/Ctrl"] + s_df_iter.iloc[-1]["Stdev"])* 1.2)

        # Determine if break is applicable for each target name
        for j in range(len(s_df_iter) - 1):
            # Check if break will be applied on axis of value
            if s_df_iter["Avg. Rel. Tx/Ctrl"][j + 1] / s_df_iter["Avg. Rel. Tx/Ctrl"][j] > break_thold:
                # Only make breaks for values greater than 1 (ctrl group)
                if s_df_iter["Avg. Rel. Tx/Ctrl"][j] > 1:
                    # Change "bar_break" to True
                    bar_break = True
                    # Determine coords of lower and upper break on axis of value for iterated target name
                    lower_break[s_names[i]] = round((s_df_iter["Avg. Rel. Tx/Ctrl"][j] + s_df_iter["Stdev"][j]) * 1.05)
                    upper_break[s_names[i]] = round((s_df_iter["Avg. Rel. Tx/Ctrl"][j + 1] - s_df_iter["Stdev"][j + 1]) * 0.95)
                    break

    # Determine upper limit for all data on axis of value
    ul = max(upper_limit.values())

    # Variable for total bars
    total_bars = len(s_names) * len(basic_pos)

    # Check if bar graph will be plotted horizontally
    if total_bars > thold_hbar_ct:
    
        # Check if plot bar graph with break
        if bar_break:

            # Bar plots
            if total_bars > 20:
                fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=[total_bars / 7.5, total_bars / 3], dpi=100)
            else:
                fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, dpi=100)

            # Counter for offside width for "pos"
            counter = 0

            # Loop through "s_names" in reverse order
            for i in range(len(s_names))[::-1]:

                # List for coords on axis of target names in bar plot
                pos = []
                [pos.append(num + (counter + 1.5) * width) for num in basic_pos]

                # Add 1 to "counter" after coord appending to "pos"
                counter += 1

                # Plot bar graph with right error bars only            
                # https://stackoverflow.com/questions/45752981/removing-the-bottom-error-caps-only-on-matplotlib
                # Determine bar plots of "ax1" for different sample names
                bp = ax1.barh(pos, s_df[s_names[i]]["Avg. Rel. Tx/Ctrl"], width, color=bar_color[f'Bar Color{i}'], alpha=alpha)
                # Set up error bar in "ax1"
                plotline1, caplines1, barlinecols1 = ax1.errorbar(s_df[s_names[i]]["Avg. Rel. Tx/Ctrl"], pos, \
                    xerr=s_df[s_names[i]]["Stdev"], xlolims=True, ls="None", color="k")
                # Set up shape and size of error bar cap in "ax1"
                caplines1[0].set_marker("|")
                caplines1[0].set_markersize(capsize)
                # Determine bar plots of "ax2" for different sample names
                ax2.barh(pos, s_df[s_names[i]]["Avg. Rel. Tx/Ctrl"], width, color=bar_color[f'Bar Color{i}'], alpha=alpha)
                # Set up error bar in "ax2"
                plotline2, caplines2, barlinecols2 = ax2.errorbar(s_df[s_names[i]]["Avg. Rel. Tx/Ctrl"], pos, \
                    xerr=s_df[s_names[i]]["Stdev"], xlolims=True, ls="None", color="k")
                # Set up shape and size of error bar cap in "ax2"                
                caplines2[0].set_marker("|")
                caplines2[0].set_markersize(capsize)

                # Add asterisk(s) for bars of treatment group(s) if applicable
                if s_names[i] != "Ctrl":

                    # Loop through "s_df" of iterated sample name
                    for j in range(len(s_df[s_names[i]])):
                        # Asterisk(s) for "s" parameter of "Axes.text()"
                        text_s = s_df[s_names[i]]["P"][j]
                        # x-coord for upper cap of iterated error bar
                        text_x = s_df[s_names[i]]["Avg. Rel. Tx/Ctrl"][j] + s_df[s_names[i]]["Stdev"][j]
                        # Check if iterated bar has break and if not ...
                        if text_x < lower_break[s_names[i]]:
                            # Add text label to bar in "ax1"
                            # Set space between upper cap of error bar and text as 1/20 of the width in "ax1"
                            ax1.text(text_x + lower_break[s_names[i]] / 20, pos[j] - 3 * width / 8, text_s)
                        # If break exists ...
                        else:
                            # Add text label to bar in "ax2"
                            # Set space between upper cap of error bar and text as 1/20 of the width in "ax2"
                            ax2.text(text_x + (upper_limit[s_names[i]] - upper_break[s_names[i]]) /20, pos[j] - 3 * width / 8, text_s)

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

        # No break in horizontal bar graph
        else:

            # Bar plots
            if total_bars > 20:
                fig, ax = plt.subplots(figsize=[total_bars / 7.5, total_bars / 3], dpi=100)
            else:
                fig, ax = plt.subplots(dpi=100)

            # Counter for offside width for "pos"
            counter = 0

            # Loop through "s_names" in reverse order
            for i in range(len(s_names))[::-1]:

                # List for coords on axis of target names in bar plot
                pos = []                
                [pos.append(ele + (counter + 1.5) * width) for ele in basic_pos]

                # Add 1 to "counter"
                counter += 1

                # Plot bar graph with right error bars only            
                # Determine bar plots for different sample names
                bp = ax.barh(pos, s_df[s_names[i]]["Avg. Rel. Tx/Ctrl"], width, color=bar_color[f'Bar Color{i}'], alpha=alpha)
                # Set up error bar
                plotline, caplines, barlinecols = ax.errorbar(s_df[s_names[i]]["Avg. Rel. Tx/Ctrl"], pos, \
                    xerr=s_df[s_names[i]]["Stdev"], xlolims=True, ls="None", color="k")
                # Set up shape and size of error bar cap
                caplines[0].set_marker("|")
                caplines[0].set_markersize(capsize) 

                # Add asterisk(s) for bars of treatment group(s) if applicable
                if s_names[i] != "Ctrl":

                    # Loop through "s_df" of iterated sample name
                    for j in range(len(s_df[s_names[i]])):
                        # Asterisk(s) for "s" parameter of "Axes.text()"
                        text_s = s_df[s_names[i]]["P"][j]
                        # x-coord for upper cap of iterated error bar
                        text_x = s_df[s_names[i]]["Avg. Rel. Tx/Ctrl"][j] + s_df[s_names[i]]["Stdev"][j]
                        # Add text label to bar in "ax1"
                        # Set space between upper cap of error bar and text as 1/20 of the width in "ax1"
                        ax.text(text_x + upper_limit[s_names[i]] / 40, pos[j] - 3 * width / 8, text_s)

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

    # Bar graph will be plotted vertically
    else:

        # Check if plot bar graph with break
        if bar_break:

            # Bar plots
            if total_bars > 20:
                fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=[total_bars / 3, total_bars / 7.5], dpi=100)
            else:
                fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, dpi=100)

            # Counter for offside width for "pos"
            counter = 0

            # Loop through "s_names" in reverse order
            for i in range(len(s_names)):

                # List for coords on axis of target names in bar plot
                pos = []
                [pos.append(ele + (counter + 1.5) * width) for ele in basic_pos]

                # Add 1 to "counter"
                counter += 1
            
                # Plot bar graph with top error bars only            
                # Determine bar plots of "ax1" for different sample names
                bp = ax1.bar(pos, s_df[s_names[i]]["Avg. Rel. Tx/Ctrl"], width, color=bar_color[f'Bar Color{i}'], alpha=alpha)
                # Set up error bar in "ax1"
                plotline1, caplines1, barlinecols1 = ax1.errorbar(pos, s_df[s_names[i]]["Avg. Rel. Tx/Ctrl"], \
                    yerr=s_df[s_names[i]]["Stdev"], lolims=True, ls="None", color="k")
                # Set up shape and size of error bar cap in "ax1"
                caplines1[0].set_marker("_")
                caplines1[0].set_markersize(capsize)
                # Determine bar plots of "ax2" for different sample names
                ax2.bar(pos, s_df[s_names[i]]["Avg. Rel. Tx/Ctrl"], width, color=bar_color[f'Bar Color{i}'], alpha=alpha)
                # Set up error bar in "ax2"
                plotline2, caplines2, barlinecols2 = ax2.errorbar(pos, s_df[s_names[i]]["Avg. Rel. Tx/Ctrl"], \
                    yerr=s_df[s_names[i]]["Stdev"], lolims=True, ls="None", color="k")
                # Set up shape and size of error bar cap in "ax2"                
                caplines2[0].set_marker("_")
                caplines2[0].set_markersize(capsize)

                # Add asterisk(s) for bars of treatment group(s) if applicable
                if s_names[i] != "Ctrl":

                    # Loop through "s_df" of iterated sample name
                    for j in range(len(s_df[s_names[i]])):

                        # Asterisk(s) for "s" parameter of "Axes.text()"
                        text_s = s_df[s_names[i]]["P"][j]
                        # x-coord for upper cap of iterated error bar
                        text_y = s_df[s_names[i]]["Avg. Rel. Tx/Ctrl"][j] + s_df[s_names[i]]["Stdev"][j]
                        # Check if iterated bar has break and if not ...
                        if text_y < lower_break[s_names[i]]:
                            # Add text label to bar in "ax2"
                            # Set space between upper cap of error bar and text as 1/20 of the width in "ax2"
                            ax2.text(pos[j], text_y + lower_break[s_names[i]] / 20, text_s, horizontalalignment='center')
                        # If break exists ...
                        else:
                            # Add text label to bar in "ax1"
                            # Set space between upper cap of error bar and text as 1/20 of the width in "ax1"
                            ax1.text(pos[j], text_y + (upper_limit[s_names[i]] - upper_break[s_names[i]]) /20, text_s, \
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
            if total_bars > 15:
                ax.tick_params(axis="x", labelrotation=30)
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

        # No break in bar graph
        else:
            # Bar plots
            if total_bars > 20:
                fig, ax = plt.subplots(figsize=[total_bars / 3, total_bars / 7.5], dpi=100)
            else:
                fig, ax = plt.subplots(dpi=100)

            # Counter for offside width for "pos"
            counter = 0

            # Loop through "s_names" in reverse order
            for i in range(len(s_names)):

                # List for coords on axis of target names in bar plot
                pos = []
                [pos.append(ele + (counter + 1.5) * width) for ele in basic_pos]

                # Add 1 to "counter"
                counter += 1

                # Plot bar graph with top error bars only            
                # Determine bar plots for different sample names
                bp = ax.bar(pos, s_df[s_names[i]]["Avg. Rel. Tx/Ctrl"], width, color=bar_color[f'Bar Color{i}'], alpha=alpha)
                # Set up error bar in "ax1"
                plotline, caplines, barlinecols = ax.errorbar(pos, s_df[s_names[i]]["Avg. Rel. Tx/Ctrl"], \
                    yerr=s_df[s_names[i]]["Stdev"], lolims=True, ls="None", color="k")
                # Set up shape and size of error bar cap
                caplines[0].set_marker("_")
                caplines[0].set_markersize(capsize)                

                # Add asterisk(s) for bars of treatment group(s) if applicable
                if s_names[i] != "Ctrl":

                    # Loop through "s_df" of iterated sample name
                    for j in range(len(s_df[s_names[i]])):

                        # Asterisk(s) for "s" parameter of "Axes.text()"
                        text_s = s_df[s_names[i]]["P"][j]
                        # x-coord for upper cap of iterated error bar
                        text_y = s_df[s_names[i]]["Avg. Rel. Tx/Ctrl"][j] + s_df[s_names[i]]["Stdev"][j]
                        # Add text label to bar in "ax"
                        # Set space between upper cap of error bar and text as 1/40 of the width in "ax"
                        ax.text(pos[j], text_y + upper_limit[s_names[i]] / 40, text_s, horizontalalignment='center')

                # Determine legend for plots of different sample names
                s_plots_legend.append(bp)

            # Set limit on axis of value
            ax.set_ylim(0, ul)
            # Set limit on axis of target names
            ax.set_xlim(0, len(s_df["Ctrl"]) + width) 

            # Set lables on axis of value
            ax.set_ylabel(value_label)

            # Set ticks and ticklabels on axis of target names
            if total_bars > 15:
                ax.tick_params(axis="x", labelrotation=30)

            ax.set_xticks(value_ticks)
            ax.set_xticklabels(list(s_df["Ctrl"]["Target Name"]))

            # Set title of bar graph
            ax.set_title(title)

            # Set legend of bar graph
            ax.legend(s_plots_legend, s_names, loc=legend_loc) 

    # Make sure rotated xlables are included in "fig"
    fig.tight_layout()    

    return fig