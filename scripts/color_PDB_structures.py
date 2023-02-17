import subprocess
import sys

import pandas as pd 
import polyclonal

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

# snakemake variables
input_pdb_file = snakemake.input.input_pdb_file
env_chains = snakemake.params.env_chains
average_escape_model = snakemake.input.average_escape_model
output_pdb_file_name = snakemake.output.output_pdb_file_name
output_file_sub_name = snakemake.params.output_file_sub_name

# read pdb, as suggested here:
# https://gist.github.com/tbrittoborges/929ced78855945f3e296
colspecs = [(0, 6), (6, 11), (12, 16), (16, 17), (17, 20), (21, 22), (22, 26),
            (26, 27), (30, 38), (38, 46), (46, 54), (54, 60), (60, 66), (76, 78),
            (78, 80)]

names = ['ATOM', 'serial', 'name', 'altloc', 'resname', 'chainid', 'resseq',
         'icode', 'x', 'y', 'z', 'occupancy', 'tempfactor', 'element', 'charge']

pdb = pd.read_fwf(input_pdb_file, names=names, colspecs=colspecs)

escape_data = pd.read_csv(average_escape_model).query('times_seen>3')

chains = env_chains.split(' ')

for epitope, df in escape_data.groupby('epitope'):
    
    # make a dataframe with site escape averages for the eptiope for each chain
    dfs = []
    for chain in chains:
        chain_df = df.copy()
        chain_df['chain'] = chain
        dfs.append(chain_df)
    df = pd.concat(dfs)
    df = df.groupby(['site', 'chain']).mean().reset_index()
    
    # polyclonal.pdb_utils.reassign_b_factor can only currently take int inputs
    df = df[~df['site'].str.contains("a|b|c|d|e|f|g")]
    df['site'] = df['site'].astype(int)
    
    # save the new PDB file
    polyclonal.pdb_utils.reassign_b_factor(input_pdbfile=input_pdb_file, 
                                           output_pdbfile=f"{output_file_sub_name}_epitope_{epitope}.pdb", 
                                           df=df, 
                                           metric_col='escape_mean')