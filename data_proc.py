import numpy as np
import pandas as pd
from scipy import stats

'''
Encode function named "data_proc()" that takes in info. of qPCR design (reference gene, control groups, and experimental repeat for 
results from each expt. saved as ".csv" file). Returns processed data in pandas DataFrame.
'''

# Function to read ".csv" file(s) and generate bar graph accordingly
def data_proc(file_list, ref_list, ctrl_list, ct_repeat, ctrl_name="WT"):

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
        # Unify the name of control group as ctrl_name in case more than one ctrl group names are encountered
        df.loc[df["Sample Name"] == ctrl, "Sample Name"] = ctrl_name
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
            for j in range(0, len(df), 3):
                # Save 3 repeated ct values in ascending order in "ct_list"
                ct_list = sorted([df["CT"][j], df["CT"][j + 1], df["CT"][j + 2]])
                # Calculate the distance in between each value in "ct_list"
                delta_upper = ct_list[2] - ct_list[1]
                delta_lower = ct_list[1] - ct_list[0]
                # Since there are only 3 parallel ct values, consider the most distant one as potential outlier
                if delta_upper > delta_lower:
                    # delta ct = 1 (2 folds) is set as the threshold for determining outlier ct value
                    if delta_upper > 1:
                        # Append index of the maximum ct value to "outlier_index"
                        [outlier_index.append(j + k) for k in range(3) if df["CT"][j + k] == ct_list[2]]
                else:
                    if delta_lower > 1:
                        # Append index of the minimum ct value to "outlier_index"
                        [outlier_index.append(j + k) for k in range(3) if df["CT"][j + k] == ct_list[0]]

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
        t_r_mean_ctrl = df.loc[df["Sample Name"] == ctrl_name, :].groupby(["Concat Name"])["Target/Ref"].mean()

        # Loop through "Target Name"
        for t_name in df["Target Name"].unique():
            # Loop through "Sample Name"
            for s_name in df["Sample Name"].unique():
                # Calculate relative ratio of Target/Ref (as compared to that of avg. value from Ctrl) for each repeat of target genes
                df.loc[df["Concat Name"] == f'{s_name}_{t_name}', "Rel. Tx/Ctrl"] = \
                    df["Target/Ref"] / t_r_mean_ctrl[f'{ctrl_name}_{t_name}']

        # Typecast values in "Rel. Tx/Ctrl" column from non-null object to float
        df["Rel. Tx/Ctrl"] = df["Rel. Tx/Ctrl"].astype(float)

        # >>>>>>>>>> Column calc in "file_data" for bar plot <<<<<<<<<< #
        # DataFrame to hold columns of iterated ".csv" file for bar plot
        file_data = pd.DataFrame(columns=["Sample Name", "Target Name", "Concat Name", "Avg. Rel. Tx/Ctrl", "Stdev", "P", "P_asterisk"])

        # Calculate values in columns of "file_data"
        file_data["Concat Name"] = df.groupby(["Concat Name"])["Rel. Tx/Ctrl"].mean().index
        file_data["Avg. Rel. Tx/Ctrl"] = df.groupby(["Concat Name"])["Rel. Tx/Ctrl"].mean().values
        file_data["Stdev"] = df.groupby(["Concat Name"])["Rel. Tx/Ctrl"].std().values

        # Lists to store "Sample Name" and "Target Name"
        s_name_list = []
        t_name_list = []

        # Loop through "file_data" to append sample and target names for each row to their designated lists
        for j in range(len(file_data)):
            s_name_list.append(file_data["Concat Name"][j].split("_")[0])
            t_name_list.append(file_data["Concat Name"][j].split("_")[1])
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
            for j in range(len(t_name_df)):
                # Variable for "Sample Name"
                s = t_name_df["Sample Name"][j]
                # Create key of "s" in "p_df" if not exists
                if not s in list(p_df.keys()):
                    p_df[s] = []
                # Append "Rel. Tx/Ctrl" value to "p_df[s]"
                p_df[s].append(t_name_df["Rel. Tx/Ctrl"][j])

            # Loop through keys of "p_df"
            for s in list(p_df.keys()):

                # Perform t-test between ctrl and treatment(s) groups
                if s != ctrl_name:
                    # Consider default variances for values between two series are equal
                    equal_var = True
                    # Calculate p_cdf for f-test (cumulative distribution function)
                    p_cdf = stats.f.cdf(np.var(p_df[s]) / np.var(p_df[ctrl_name]), len(p_df[s]) - 1, len(p_df[ctrl_name]) - 1)
                    # Make sure variance1 is greater than variance2 to determine "p_cdf"
                    if np.var(p_df[ctrl_name]) > np.var(p_df[s]):
                        p_cdf = 1 - p_cdf                    
                    # Check whether p-cdf is less than alpha which is by default 0.05 and if true ...
                    if p_cdf < 0.05:
                        # Reject the null hypothesis and accept that the variances between two series are not equal
                        equal_var = False
                    # Calculate p value for t-test
                    p = stats.ttest_ind(p_df[s], p_df[ctrl_name], equal_var=equal_var)[1]
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
                    file_data.loc[file_data["Concat Name"] == f'{s}_{t_name}', "P"] = p
                    file_data.loc[file_data["Concat Name"] == f'{s}_{t_name}', "P_asterisk"] = asterisk                
                else:
                    # Append "" to "P" column of ctrl_name in "file_data"
                    file_data.loc[file_data["Concat Name"] == f'{s}_{t_name}', "P"] = ""
                    file_data.loc[file_data["Concat Name"] == f'{s}_{t_name}', "P_asterisk"] = ""

        # Drop "Concat Name" column from "file_data"
        file_data = file_data.drop(columns=["Concat Name"])
        # Add "file_data" to "data"
        data = pd.concat([data, file_data], sort=True, ignore_index=True)

    # Return processed data
    return data