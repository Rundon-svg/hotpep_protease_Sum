import os
import pandas as pd
import argparse

# Define the order of Merops family categories by their prefix letters
MEROPS_CATEGORIES = ["A", "C", "G", "M", "N", "P", "S", "T", "U"]

def parse_family(family):
    """Parse the Merops family to extract category, family, and subfamily for sorting purposes."""
    category = family[0]
    if len(family) == 1:
        return (category, "", "")
    if len(family) == 2:
        return (category, family[:2], "")
    return (category, family[:2], family)

def summarize_peptidase_results(input_path):
    # Create a DataFrame to combine results from all species
    combined_df = pd.DataFrame()
    
    # Iterate through each species folder
    for species_folder in os.listdir(input_path):
        species_path = os.path.join(input_path, species_folder)
        if os.path.isdir(species_path):
            summary_file = os.path.join(species_path, "peptidases", "summary.txt")
            if os.path.exists(summary_file):
                df = pd.read_csv(summary_file, sep="\t", usecols=["Merops family", "proteins"])
                df.rename(columns={"proteins": species_folder}, inplace=True)
                df.set_index("Merops family", inplace=True)
                combined_df = pd.concat([combined_df, df], axis=1)
    
    # Fill missing values with 0
    combined_df.fillna(0, inplace=True)
    
    # Parse Merops family and sort by category, family, and subfamily
    combined_df["Category"], combined_df["Family"], combined_df["SubFamily"] = zip(*combined_df.index.map(parse_family))
    combined_df["Category"] = pd.Categorical(combined_df["Category"], categories=MEROPS_CATEGORIES, ordered=True)
    combined_df.sort_values(["Category", "Family", "SubFamily"], inplace=True)
    
    # Remove temporary columns used for sorting
    combined_df.drop(columns=["Category", "Family", "SubFamily"], inplace=True)
    
    # Save the combined summary table
    combined_df.to_csv("combined_summary.csv")
    
    # Initialize summary statistics DataFrame
    summary_stats = pd.DataFrame()
    
    # Calculate total gene count
    total_row = pd.DataFrame(combined_df.sum(axis=0)).T
    total_row.index = ["Total"]
    summary_stats = pd.concat([summary_stats, total_row], axis=0)
    
    # Calculate gene count by Catalytic Type (Category), ensure all categories are represented even if zero
    for cat in MEROPS_CATEGORIES:
        cat_df = combined_df[combined_df.index.str.startswith(cat)]
        if not cat_df.empty:
            cat_row = pd.DataFrame(cat_df.sum(axis=0)).T
        else:
            # If the category is missing, create a row of zeros
            cat_row = pd.DataFrame([0] * combined_df.shape[1], index=combined_df.columns).T
        cat_row.index = [f"{cat}"]
        summary_stats = pd.concat([summary_stats, cat_row], axis=0)
    
    # Calculate gene count for each Family and SubFamily without renaming family
    for family in combined_df.index.unique():
        family_df = combined_df.loc[family]
        if isinstance(family_df, pd.Series):
            family_df = pd.DataFrame(family_df).T
        family_row = pd.DataFrame(family_df.sum(axis=0)).T
        family_row.index = [family]  # Use original family name
        summary_stats = pd.concat([summary_stats, family_row], axis=0)
    
    # Save the summary statistics table
    summary_stats.to_csv("summary_statistics.csv")
    print("Summary completed! 'combined_summary.csv' and 'summary_statistics.csv' have been generated.")

def merge_subfamilies(summary_stats_path):
    """Modify summary_statistics.csv to merge subfamilies into their family level if the -family argument is used."""
    df = pd.read_csv(summary_stats_path, index_col=0)
    
    # Initialize merged_df with the same columns as df
    merged_df = pd.DataFrame(columns=df.columns)
    
    # Process each family and subfamily
    for index in df.index:
        family = index[:3] if len(index) > 3 else index  # Extract the family name (first 3 characters)
        
        # If family is already in merged_df, sum the values, otherwise add the new family
        if family in merged_df.index:
            merged_df.loc[family] += df.loc[index]
        else:
            merged_df.loc[family] = df.loc[index]
    
    # Save the modified summary_statistics.csv
    merged_df.to_csv("summary_statistics.csv")
    print("Subfamily rows merged successfully! Updated 'summary_statistics.csv' has been saved.")

# Main entry point
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize peptidase prediction results across multiple species.")
    parser.add_argument("-in", "--input_path", required=True, help="Path to the folder containing all species results.")
    parser.add_argument("-family", action="store_true", help="Merge subfamilies into family level in summary statistics.")
    args = parser.parse_args()
    
    summarize_peptidase_results(args.input_path)
    
    # If -family argument is used, merge subfamilies in summary_statistics.csv
    if args.family:
        merge_subfamilies("summary_statistics.csv")
