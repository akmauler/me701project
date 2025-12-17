import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import io
import math
import os

# Define the files and their respective multiplication factors (Dose normalization)
FILES_AND_MULTIPLIERS = {
    "Result_ATP-AM_VARIAN.out": 72,
    "Result_ATP-AM_VARIAN-field2.out": 71,
    "Result_ATP-AM_VARIAN-field3.out": 69,
    "Result_ATP-AM_VARIAN-field4.out": 95
}

def parse_phits_file(filename, multiplier=1.0):
    """
    Parses a single PHITS output file to extract the data table (Dose and R.Err).
    Applies the multiplier to the 'all' (Dose) column.
    """
    
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return None

    data_lines = []
    in_data_block = False
    header_start = "#  num    reg     volume       all       r.err"
    data_end = "#   sum over"

    for line in lines:
        if line.strip() == header_start:
            in_data_block = True
            continue
        if line.strip().startswith(data_end):
            in_data_block = False
            break

        if in_data_block and line.strip():
            # get rid of headers that may be in line
            if not line.strip().startswith('#'):
                data_lines.append(line)

    if not data_lines:
        print(f"No region data block found in {filename}.")
        return None

    try:
        data_string = "".join(data_lines)
        df = pd.read_csv(
            io.StringIO(data_string),
            sep='\s+',
            header=None,
            # 'all' is the Dose value
            names=['num', 'reg', 'volume', 'all', 'r.err']
        )
        
        df['reg'] = df['reg'].astype(int)
        df.set_index('reg', inplace=True)
        # Apply multiplication factor 
        df['all'] = df['all'] * multiplier
        return df
        
    except Exception as e:
        print(f"Error parsing region data block in {filename}: {e}")
        return None

def extract_summary_data(filename, multiplier):
    """
    Extracts and scales the 'sum over' line (Total Volume, Dose, and R.Err)
    from a PHITS output file for whole body dose.
    """
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        print(f"Error accessing summary data in {filename}: {e}")
        return None 

    data_end = "#   sum over"
    
    for line in reversed(lines):
        if line.strip().startswith(data_end):
            parts = line.split()
            if len(parts) >= 6:
                try:
                    volume = float(parts[3])   
                    dose_raw = float(parts[4])  
                    r_err_raw = float(parts[5])
                    
                    # multiplier
                    dose_scaled = dose_raw * multiplier
                    # scale to be absolute error
                    abs_err_scaled = dose_scaled * r_err_raw
                    
                    return {
                        'filename': filename,
                        'volume': volume,
                        'dose': dose_scaled,
                        'r_err_raw': r_err_raw,
                        'abs_err': abs_err_scaled
                    }
                except ValueError as ve:
                    print(f"Error parsing numbers in summary line of {filename}: {ve}")
                    return None
            return None 
    return None 

def calculate_combined_data(file_multiplier_dict):
    """
    Reads all files, combines doses per region, and calculates the propagated error.
    Returns the combined data with total dose and total relative error.
    """
    all_dfs = []
    
    for f, mult in file_multiplier_dict.items():
        df = parse_phits_file(f, mult)
        if df is not None:
            all_dfs.append(df)

    if not all_dfs:
        print("Error: Could not process any files :(")
        return None

    # Step 1: Calculate Absolute Error for each file/region
    # Abs. Error = Dose * Relative Error
    for df in all_dfs:
        df['abs_err'] = df['all'] * df['r.err']

    # Step 2: Combine data by region (sum doses and square of absolute errors)
    # Concat will join together the strings
    combined_data = pd.concat(all_dfs)
    # Group the data by region so the regions all match up when summing.
    grouped_data = combined_data.groupby(level='reg')
    
    # Sum the dose data for all the files
    total_dose = grouped_data['all'].sum()
    # Sum the squares of absolute errors
    sum_sq_abs_err = grouped_data['abs_err'].apply(lambda x: (x**2).sum())
    
    # Step 3: Calculate the combined absolute error
    total_abs_err = np.sqrt(sum_sq_abs_err)
    
    combined_df = pd.DataFrame({
        'total_dose': total_dose,
        'total_abs_err': total_abs_err
    })
    
    # Step 4: Calculate the new combined relative error
    # R.Err = Abs.Err / Dose_Total
    combined_df['total_r_err'] = (
        combined_df['total_abs_err'] / combined_df['total_dose']
    ).fillna(0).replace(np.inf, 0) # Handle division by zero/inf cases
    
    return combined_df

def process_and_print_summary_data(file_multiplier_dict):
    """
    Extracts, prints individual file summaries, and calculates the combined total summary.
    """
    print("\nSummary Data")
    
    summary_data_list = []
    
    for f, mult in file_multiplier_dict.items():
        data = extract_summary_data(f, mult)
        if data:
            summary_data_list.append(data)

    if not summary_data_list:
        print("No summary data could be processed.")
        return
    
    # Combined summary data
    
    # Volume is expected to be the same across all files for PHITS T-Deposit
    combined_volume = summary_data_list[0]['volume'] if summary_data_list else 0

    total_combined_dose = sum(data['dose'] for data in summary_data_list)
    sum_sq_abs_err = sum(data['abs_err']**2 for data in summary_data_list)
    total_abs_err = math.sqrt(sum_sq_abs_err)
    
    if total_combined_dose != 0:
        total_r_err = total_abs_err / total_combined_dose
    else:
        total_r_err = 0.0
        
    print("\nTotal combined summary results")
    print(f"Total Volume (sum over): {combined_volume:.4E}")
    print(f"Total Combined Dose:     {total_combined_dose:.4E} Gy")
    print(f"Total Absolute Error:    {total_abs_err:.4E} Gy")
    print(f"Total Relative Error:    {total_r_err:.4f}")


def plot_combined_dose(combined_df):
    """
    Generates a single plot of Region Number vs. Total Dose (Gy) 
    using a linear y-axis scale. It automatically adapts to all non-zero regions.
    """
    if combined_df is None or combined_df.empty:
        print("No combined data available for plotting.")
        return
        
    # Use all regions in the combined data that have a total dose > 0
    plot_df = combined_df[combined_df['total_dose'] > 0].copy()
    
    if plot_df.empty:
        print("No regions with non-zero dose found for plotting.")
        return

    # Sort the data by region number 
    plot_df.sort_index(inplace=True)
    
    np_regions = plot_df.index.values
    np_total_doses = plot_df['total_dose'].values
    np_total_abs_errors = plot_df['total_abs_err'].values

    # Plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    # Error bars
    ax.errorbar(
        np_regions, 
        np_total_doses, 
        yerr=np_total_abs_errors,
        marker='o', 
        linestyle='', 
        color="#1b590f",
        capsize=5,
        label='Total Combined Dose (Gy)',
        zorder=3
    )
    
    x_min = np_regions.min()
    x_max = np_regions.max()
    ax.set_xlim(x_min - 0.5, x_max + 0.5)
    ax.set_xticks(np_regions)
    ax.set_xticklabels([str(r) for r in np_regions], rotation=45, ha='right', fontsize=9)
    y_max = np_total_doses.max()
    y_min = np_total_doses[np_total_doses > 0].min() if np_total_doses[np_total_doses > 0].size > 0 else 0
    y_range = y_max - y_min
    y_buffer = y_range * 0.1 if y_range > 0 else 0.1
    
    ax.set_ylim(max(0, y_min - y_buffer), y_max + y_buffer)
    
    ax.set_title('Combined Simulated Dose vs. Region Number', fontsize=14)
    ax.set_ylabel('Total Combined Dose (Gy)', fontsize=12)
    ax.set_xlabel('Region Number', fontsize=12)
    ax.grid(True, which="major", ls=":", alpha=0.7, zorder=0)
    ax.legend(loc='upper right', fontsize=10)
    
    plt.tight_layout()
    plt.show()

def print_regional_data(combined_df):
    """
    Prints the combined dose and error for all regions with non-zero dose.
    """
    print("\nCombined region based dose results")
    if combined_df is None or combined_df.empty:
        print("No combined data available.")
        return
        
    # Filter for regions with non-zero total dose
    df_to_print = combined_df[combined_df['total_dose'] > 0].copy()
    
    if df_to_print.empty:
        print("All combined regional doses are zero. Nothing to report.")
        return
        
    df_to_print.index.name = 'Region'
    df_to_print['Total Dose (Gy)'] = df_to_print['total_dose'].map(lambda x: f"{x:.4E}")
    df_to_print['Abs. Error (Gy)'] = df_to_print['total_abs_err'].map(lambda x: f"{x:.4E}")
    df_to_print['Rel. Error (R.Err)'] = df_to_print['total_r_err'].map(lambda x: f"{x:.4f}")
    
    print(df_to_print[['Total Dose (Gy)', 'Abs. Error (Gy)', 'Rel. Error (R.Err)']].to_string())


if __name__ == "__main__":
        
    # 1. Process and Print Summary Data
    process_and_print_summary_data(FILES_AND_MULTIPLIERS)

    # 2. Calculate Combined Region Data 
    combined_df = calculate_combined_data(FILES_AND_MULTIPLIERS)

    # 3. Print All Regional Results 
    print_regional_data(combined_df)

    # 4. Plot
    plot_combined_dose(combined_df)
