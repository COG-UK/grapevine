import csv
from collections import defaultdict
from collections import Counter
from collections import OrderedDict
import matplotlib.pyplot as plt
import datetime as dt
from epiweeks import Week,Year
import pandas as pd
import numpy as np
import seaborn as sns


class taxon():
    
    def __init__(self, record_id, uk_lineage, date, adm1_input, adm2, epiweek_input):
    
        self.id = record_id
        self.uk_lineage = uk_lineage
        self.adm1 = self.sort_adm1(adm1_input)

        if adm2 == "Needs_manual_curation":
            self.adm2 = ""
        else:
            self.adm2 = adm2

        if epiweek_input != "0" and epiweek_input != "":
            self.epiweek = Week(2020, int(float(epiweek_input)))
        elif epiweek_input == "0":
            self.epiweek = Week(2019, 52)
        elif epiweek_input == "":
            self.epiweek = "NA"
        
        if "/" in date:
            print("ERROR DATE FORMAT INCORRECT")
        
        if date == "None":
            self.date = "NA"
        self.date_dt = self.make_date_object(date)
       
    def sort_adm1(self, adm1_input):
        contract_dict = {"ENG":"England", "WLS": "Wales", "NIR": "Northern_Ireland", "SCT": "Scotland"}
        if "-" in adm1_input:
            country_prep = adm1_input.split("-")[1]
            adm1 = contract_dict[country_prep]
        else:
            adm1 = adm1_input
        return adm1
        

    def make_date_object(self,date):
        date_bits = date.split("-")
        if len(date_bits) == 3:
            date = dt.date(int(date_bits[0]), int(date_bits[1]), int(date_bits[2]))
        else:
            date = "NA"
        return date


def parse_metadata(metadata_file):

    taxon_dict = {}

    with open(metadata_file) as f:
        reader = csv.DictReader(f)
        in_data = [r for r in reader]
        for sequence in in_data:
            if sequence['country'] == "UK":
                seq_name = sequence['sequence_name']
                adm1 = sequence["adm1"]
                adm2 = sequence['adm2']
                date = sequence['sample_date']
                epiweek = sequence['epi_week']
                uk_lineage = sequence['uk_lineage']

                new = taxon(seq_name, uk_lineage, date, adm1, adm2, epiweek)

                taxon_dict[seq_name] = new

    return taxon_dict


def parse_snp_data(snp_file, snp_list):

    query_to_snps = defaultdict(list)
    snp_to_queries = defaultdict(list)

    with open(snp_file) as f:
        reader = csv.DictReader(f)
        data = [r for r in reader]
        for line in data:
            seq_name = line["query"]
            snps = line["variants"]

            for snp in snps.split("|"):
                actual_snp = snp.split(":")[1] #not sure about this, could just pull it including the gene name
                if actual_snp in snp_list:
                    query_to_snps[seq_name].append(actual_snp)
                    snp_to_queries[actual_snp].append(seq_name)

    for i in snp_list:
        if i not in snp_to_queries.keys():
            snp_to_queries[i] = []

            
    return query_to_snps, snp_to_queries

def make_snp_table(snp_to_queries, taxon_dict, adm2_count_dict):

    df_dict = defaultdict(list)
    snp_to_dates = defaultdict(list)

    for snp, query_list in snp_to_queries.items():
        adm2s = []
        dates = []
        num_seqs = len(query_list)
        
        if num_seqs > 0:

            for query in query_list:
                if query in taxon_dict:
                    dates.append(taxon_dict[query].date_dt)

            num_adm2s = len(adm2_count_dict[snp])

            if len(dates) > 0:
                earliest_date = min(dates)
                latest_date = max(dates)
                date_range = f'{earliest_date} to {latest_date}'
            else:
                date_range = "NA"

        else:
            num_adm2s = 0
            date_range = "NA"

        df_dict["SNP of interest"].append(snp)
        df_dict["Number of sequences in the UK"].append(num_seqs)
        df_dict["Date range"].append(date_range)
        df_dict["Number of adm2s present"].append(num_adm2s)

        snp_to_dates[snp] = dates

    df = pd.DataFrame(df_dict)
    df.set_index("SNP of interest", inplace=True)

    return df, snp_to_dates


def make_overall_lines(taxon_dict, snp_to_dates, figdir, focal_snp):

    ## Get total number of samples 
    all_dates = []
    for taxon in taxon_dict.values():
        if taxon.date_dt != "NA":
            all_dates.append(taxon.date_dt)


    total_df_dict = defaultdict(list)
    total_day_dict = {}

    all_date_counter = Counter(all_dates)
    all_date_counter_sorted = OrderedDict(sorted(all_date_counter.items()))

    for date, counts in all_date_counter_sorted.items():
        total_df_dict["Day"].append(date)
        total_df_dict["Counts"].append(counts)

    df_all = pd.DataFrame(total_df_dict)

    rolling_mean = df_all.Counts.rolling(window=7).mean()

    for mean, day in zip(rolling_mean, total_df_dict["Day"]):
        total_day_dict[day] = mean


    ## Get samples per snp per day ane seven day rolling average
    snp_dict = defaultdict(dict)

    for snp, dates in snp_to_dates.items():

        df_dict = defaultdict(list)
        day_dict = {}
        
        date_counter = Counter(dates)
        date_counter_sorted = OrderedDict(sorted(date_counter.items()))
        
        if len(date_counter_sorted) > 0:
            
            for date, counts in date_counter_sorted.items():
                df_dict["Day"].append(date)
                df_dict["Counts"].append(counts)

            df = pd.DataFrame(df_dict)
            rolling_mean = df.Counts.rolling(window=7).mean()
            for mean, day in zip(rolling_mean, df_dict["Day"]):
                day_dict[day] = mean

            snp_dict[snp] = day_dict
            
        else:
            for day in total_day_dict:
                day_dict[day] = 0
            
            snp_dict[snp] = day_dict


    ## Divide snp seqs by total number of seqs 
    new_snp_dict = defaultdict(dict)
    for snp, dictionary in snp_dict.items():
        new_dict = {}
        for day, mean in dictionary.items():
            freq = mean/total_day_dict[day]
            new_dict[day] = freq
            if np.isnan(mean):
                new_dict[day] = 0
            
        for days in total_day_dict.keys():
            if days not in new_dict:
                new_dict[days] = 0
                
        new_dict_sorted = OrderedDict(sorted(new_dict.items()))
        new_snp_dict[snp] = new_dict_sorted

    plt.rc('xtick',labelsize=20)
    plt.rc('ytick',labelsize=20)
        

    fig, ax = plt.subplots(1,1,figsize=(30,10))
    used_ys = []
    for snp, dictionary in new_snp_dict.items():
        
        x=list(dictionary.keys())
        y=list(dictionary.values())
        
        snp_label = snp

        a = ax.plot(x,y, linewidth=2, label=snp_label)
        color = a[-1].get_color()
        
        if y[-1] not in used_ys:
            ax.text(x[-1], y[-1], snp_label, fontsize=20, color=color)
            used_ys.append(y[-1])
        else:
            for i in np.arange(0,1,0.05):
                if y[-1]+i not in used_ys:
                    ax.text(x[-1], y[-1]+i, snp_label, fontsize=20, color=color)
                    used_ys.append(y[-1]+i)
                    break
        
        
    ax.legend(fontsize=20)
    ax.set_ylabel("Frequency of SNP (7 day rolling average)", fontsize=20)
    ax.set_xlabel("Date", fontsize=20)

    if not focal_snp:
        plt.savefig(f"{figdir}/all_snps_line.svg", format='svg')
    else:
        plt.savefig(f"{figdir}/{focal_snp}_line.svg", format='svg')


def make_heatmap(snps, query_to_snps, figdir):

    mut_dict = defaultdict(dict)
    for snp in snps:
        mut_dict[snp] = Counter()

    # add all pairwise hits to the mut dict
    for query, snp_list in query_to_snps.items():
        for i in snps:
            for j in snps:
                if i != j and i in snp_list and j in snp_list:
                    mut_dict[i][j]+=1
    
    # if there's no link then add it to the counter as 1, so that when it becomes logged it == 0
    for i in mut_dict:
        match = mut_dict[i]
        for mut in snps:
            if mut not in match:
                mut_dict[i][mut] = 1
    
    data = {}
    index = []
    
    for i in sorted(mut_dict):
        index = [j for j in sorted(mut_dict[i])] # sorted list of mut dict keys
        data[i] = [mut_dict[i][j] for j in sorted(mut_dict[i])] # putting the data into the dataframe, with col names the sorted list of mut_dict keys
    
    df = pd.DataFrame(data, columns = index, index=index)
    df_log = df.transform(lambda x : np.log10(x))

    muted_pal = sns.cubehelix_palette(as_cmap=True) # get nice colour palette
    sns.heatmap(df_log, cmap=muted_pal, annot=True, fmt='.2g',square=True, cbar_kws={'label': 'log(count)'})  

    # Heatmap args in order:
    # - the logged data frame
    # - the colour map
    # - annotate set to true to display number
    # - the fmt='.2g' bit rounds numbers to 2 decimal places
    # - square layout of the heatmap
    # - cbar_kws customises the heatmap on the side

    plt.savefig(f"{figdir}/pairwise_cooccurance.svg", format="svg")


