import datetime

date = datetime.date.today()

# rule combine_gisaid_and_cog:
#     input:
#         gisaid_fasta =
#         gisaid_metadata =
#         uk_fasta =
#         uk_metadata =
#     params:
#         prefix = config["output_path"] + "/cog_gisaid_%s_" %date
#     output:
#         fasta =
#         metadata =
#     shell:
#         """
#         fastafunk merge
#
#         num_seqs = $(cat {} | grep ">" | wc -l)
#         mv
#         """

