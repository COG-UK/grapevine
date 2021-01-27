
from collections import defaultdict
from collections import Counter
from collections import OrderedDict
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import os

import mutation_funcs as mfunk
import mapping as map_funks


parser = argparse.ArgumentParser(description='Run mutation report')

parser.add_argument('--date-data', help="When the report was run", dest="date_data")
parser.add_argument("--metadata-file", dest="metadata_file")
parser.add_argument("--snp-file", dest="snp_file")
parser.add_argument("--snp-list",help="list of snps desired in the report", dest="snp_list")
parser.add_argument("--figdir")
parser.add_argument("--snps-for-matrix", dest="snps_for_matrix", default=None)
parser.add_argument("--uk-map", dest="uk_map")
parser.add_argument("--ni-counties", dest="ni_counties")
parser.add_argument("--channel-islands", dest="channel_islands")
parser.add_argument("--outdir")

args = parser.parse_args()

date_data = args.date_data
metadata_file = args.metadata_file
snp_file = args.snp_file
snp_list = args.snp_list
figdir = args.figdir
snps_for_matrix = args.snps_for_matrix
outdir = args.outdir

uk = args.uk_map
ni = args.ni_counties
channels = args.channel_islands

if not os.path.exists(outdir):
    os.mkdir(outdir)

if outdir not in figdir:
    figdir_writing = f'{outdir}/{figdir}'
else:
    figdir_writing = figdir
    figdir = figdir_writing.lstrip(outdir)

if not os.path.exists(figdir_writing):
    os.mkdir(figdir_writing)

map_files = [uk,channels,ni]

all_uk = map_funks.generate_all_uk_dataframe(map_files)

if type(snp_list) == str:
    snp_list = snp_list.split(",")
if type(snps_for_matrix) == str:
    snps_for_matrix = snps_for_matrix.split(",")

if not snps_for_matrix:
    snps_for_matrix = snp_list


##Data procesing and manipulation
print("Parsing metadata")
taxon_dict = mfunk.parse_metadata(metadata_file)
query_to_snps, snp_to_queries = mfunk.parse_snp_data(snp_file, snp_list)

print("Making heatmap")
mfunk.make_heatmap(snps_for_matrix, query_to_snps, figdir_writing)

print("Making maps and sorting adm2 out")
adm2_perc_dict = defaultdict(dict)
adm2_count_dict = defaultdict(dict)

for snp in snp_list:
    relevant_taxa = {}
    relevant_taxa_names = snp_to_queries[snp]
    for tax in relevant_taxa_names:
        if tax in taxon_dict:
            relevant_taxa[tax] = taxon_dict[tax]
    output = map_funks.map_adm2(relevant_taxa, all_uk, figdir_writing, snp)
    if type(output) != bool:
        adm2_counter, adm2_percentages = output
        adm2_perc_dict[snp] = adm2_percentages
        adm2_count_dict[snp] = adm2_counter
    else:
        adm2_perc_dict[snp] = False

print("Making tables")
df, snp_to_dates = mfunk.make_snp_table(snp_to_queries, taxon_dict, adm2_count_dict)
df.to_csv(f"{outdir}/SNP_summary_table.csv")

print("Making line figure")
mfunk.make_overall_lines(taxon_dict,snp_to_dates,figdir_writing, None)

# print("Making heatmap")
# mfunk.make_heatmap(snps_for_matrix, query_to_snps, figdir_writing)

## Writing the report ##
print("Writing the report")
fw = open(f"{outdir}/mutation_report.md", 'w')

fw.write("# Mutation report\n")

fw.write(f'Date report run: {date_data} \n\n')

fw.write("Please note that adm2s can be joined together if there are sequences with ambiguous location data, so the adm2 number should only be used as an estimate.\n\n")

fw.write(df.to_markdown())
fw.write("\n\n")

fw.write("## SNP summaries\n")
fw.write("### COG-UK amino acid variants by date\n\n")

fw.write(f'![]({figdir}/all_snps_line.svg)')
fw.write("\n\n")

fw.write("### Co-occurence matrix\n\n")

fw.write(f'![]({figdir}/pairwise_cooccurance.svg)')
fw.write("\n\n")

for snp in snp_list:
    adm2_in_map = adm2_perc_dict[snp]
    adm2_counts = adm2_count_dict[snp]
    fw.write(f"### {snp}\n\n")

    small_snp_dict = defaultdict(list)
    small_snp_dict[snp] = snp_to_dates[snp]
    mfunk.make_overall_lines(taxon_dict,small_snp_dict,figdir_writing, snp)

    fw.write(f'![]({figdir}/{snp}_line.svg)')
    fw.write("\n\n")

    if adm2_in_map:
        sorted_adm2_in_map =  {k: v for k, v in sorted(adm2_in_map.items(), key=lambda item: item[1], reverse=True)}
        fw.write(f"There are sequences with {snp} from {str(len(adm2_in_map))} admin2 regions\n")
        if len(adm2_in_map) > 5:
            fw.write("The top five are:\n")
        else:
            fw.write("This is divided into:\n")
        
        count = 0
        for adm2,percentage in sorted_adm2_in_map.items():
            if count < 5:
                fw.write(f'- {percentage}% ({adm2_counts[adm2]}) in {adm2}\n')
                count += 1
        
        fw.write("\n")
        fw.write(f'![]({figdir}/{snp}_map.svg)')
        fw.write("\n\n")
    else:
        fw.write("There is no geographical data for this SNP\n\n")


fw.close()
