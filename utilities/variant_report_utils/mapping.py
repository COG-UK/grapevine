import pandas as pd 
from collections import defaultdict
import geopandas
import matplotlib.pyplot as plt
from collections import Counter

def generate_all_uk_dataframe(mapping_input):

    uk_file = mapping_input[0]
    channel_file = mapping_input[1]
    ni_file = mapping_input[2]

    UK = geopandas.read_file(uk_file)
    NI = geopandas.read_file(ni_file)
    channels = geopandas.read_file(channel_file)
    
    ###GET NI COUNTIES

    ni_name = []
    for i in range(len(NI["CountyName"])):
        ni_name.append("Northern Ireland C")

    NI["NAME_2"] = NI["CountyName"]
    NI["NAME_1"] = ni_name  

    all_uk = UK.append(channels).append(NI)

    return all_uk


def find_ambiguities(adm2s):

    ambiguous = []
    ambiguous_dict = defaultdict(set)
    clusters = []

    for adm2 in adm2s:
        if "|" in adm2:
            ambiguous.append(set(adm2.split("|")))

    for group in ambiguous:
        for group2 in ambiguous:
            if group & group2:
                group |= group2

        clusters.append(group)

    for cluster in clusters:
        for place in cluster:
            ambiguous_dict[place] = "|".join(sorted(cluster))
    
    return ambiguous_dict

def prep_mapping_data(all_uk, tax_dict):

    adm2s = []

    for tax in tax_dict.values():
        if tax.adm2 != "":
            adm2s.append(tax.adm2) #should already be upper and with underscores

    if len(adm2s) == 0:
        return False

    ambiguous_dict = find_ambiguities(adm2s)

    ###CAPITALISE GEOJSON DATA

    uppers = []
    for i in all_uk["NAME_2"]:
        uppers.append(i.upper().replace(" ","_").replace(",",""))

    all_uk["NAME_2"] = uppers

    original = all_uk.copy()

    ##DEAL WITH MERGED LOCATIONS EG WEST MIDLANDS


    for_merging = []

    for location in all_uk["NAME_2"]:
        if location in ambiguous_dict:
            for_merging.append(ambiguous_dict[location])
        else:
            for_merging.append(location)

    all_uk["Multi_loc"] = for_merging

    merged_locs = all_uk.dissolve(by="Multi_loc")

    mergeds = []
    for multi_loc in merged_locs.index:
        mergeds.append(multi_loc)

    merged_locs["NAME_2"] = mergeds    

    ###ADD MERGED AND NON-MERGED DATABASES TOGETHER

    result = pd.merge(merged_locs, original, how="outer")

    return all_uk, result, adm2s, ambiguous_dict

def make_centroids_get_counts(result, adm2s, ambiguous_dict):

    not_mappable = ["WALES", "OTHER", "UNKNOWN", "UNKNOWN_SOURCE", "NOT_FOUND", "GIBRALTAR", "FALKLAND_ISLANDS", "CITY_CENTRE"]

    centroid_df = defaultdict(list)
    centroid_dict = {}

    for name, geometry in zip(result["NAME_2"], result["geometry"]):
        centroid_dict[name.upper()] = geometry.centroid

    centroid_counts = Counter(adm2s)
    to_remove = set()

    for_iterating = centroid_counts.copy()

    for k,v in for_iterating.items():
        if k in ambiguous_dict:
            testing = ambiguous_dict[k]
            for location in centroid_counts.keys():
                if "|" in location:
                    if any([i for i in testing.split("|") if i in location.split("|")]): #in case it's in there in a different order to the value in the ambiguity dict
                        centroid_counts[location] += v
                        to_remove.add(k)
                        break
        elif "|" in k and k not in centroid_dict:
            for location in ambiguous_dict.values():
                if any([i for i in k.split("|") if i in location.split("|")]):
                    centroid_counts[location] += v
                    to_remove.add(k)
                    break

    for loc in to_remove:
        del centroid_counts[loc]


    for adm2, count in centroid_counts.items():
        centroid_df["Adm2"].append(adm2)
        centroid_df["geometry"].append(centroid_dict[adm2])
        centroid_df["seq_count"].append(count)

    centroid_geo = geopandas.GeoDataFrame(centroid_df)

    return centroid_geo, centroid_counts

def make_map(centroid_geo, all_uk, figdir, snp):

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(20, 15)

    all_uk = all_uk.to_crs("EPSG:3395")
    centroid_geo.crs = "EPSG:4326"
    centroids_final = centroid_geo.to_crs(all_uk.crs)

    base = all_uk.plot(ax=ax, color="steelblue")

    centroids_final.plot(ax=base, color="goldenrod", markersize=centroids_final["seq_count"]/10)

    ax.axis("off")

    plt.savefig(f'{figdir}/{snp}_map.svg', format="svg")



def map_adm2(tax_dict, mapping_json_files, figdir, snp): #So this takes adm2s and plots them onto the whole UK

    output = prep_mapping_data(mapping_json_files, tax_dict)

    if type(output) == bool:
        # print("None of the sequences provided have adequate adm2 data and so cannot be mapped")
        return False
    else:
        all_uk, result, adm2s, ambiguous_dict = output

    centroid_geo, adm2_counter = make_centroids_get_counts(result, adm2s, ambiguous_dict)

    make_map(centroid_geo, all_uk, figdir,snp)

    adm2_percentages = {}

    total = len(tax_dict)

    for adm2, count in adm2_counter.items():
        adm2_percentages[adm2] = round(((count/total)*100),2)

    return adm2_counter, adm2_percentages