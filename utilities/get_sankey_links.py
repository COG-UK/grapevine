import csv, sys, math
from collections import defaultdict

all_traits = sys.argv[1]
updated_traits = sys.argv[2]
out_file = sys.argv[3]

old_data = []
with open(all_traits) as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        old_data.append(row)

new_data = []
with open(updated_traits) as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        new_data.append(row)

threshold = math.floor(len(new_data) / 180)

taxa_to_lineage_old = {x['taxon']: x['uk_lineage'] for x in old_data if x['uk_lineage'][0:2] == 'UK'}
taxa_to_lineage_new = {x['taxon']: x['uk_lineage'] for x in new_data if x['uk_lineage'][0:2] == 'UK'}

################################################################################
target_dict = {}
for taxon in taxa_to_lineage_new.keys():
    new_lineage = taxa_to_lineage_new[taxon]
    if taxon in taxa_to_lineage_old:
        old_lineage = taxa_to_lineage_old[taxon]
    else:
        old_lineage = 'new'

    if new_lineage in target_dict:
        if old_lineage in target_dict[new_lineage]:
            target_dict[new_lineage][old_lineage] += 1
        else:
            target_dict[new_lineage][old_lineage] = 1
    else:
        target_dict[new_lineage] = {old_lineage: 1}

target_dict['lost'] = {}
for taxon in taxa_to_lineage_old.keys():
    old_lineage = taxa_to_lineage_old[taxon]
    if taxon not in taxa_to_lineage_new:
        if old_lineage in target_dict['lost']:
            target_dict['lost'][old_lineage] += 1
        else:
            target_dict['lost'][old_lineage] = 1


###############################################################################
links = []
# now aggregate over target_dict to populate the links data
for target, source_dictionary in target_dict.items():
    for source in source_dictionary:
        if source == 'new':
            continue
        value = source_dictionary[source]
        if value > threshold:
            links.append({'source': source, 'target': target, 'value': value})

# get a list of the impactful targets for later
targets_in_links = []
for link in links:
    targets_in_links.append(link['target'])
if 'lost' not in targets_in_links:
    targets_in_links.append('lost')

# get a list of the impactful sources for later
sources_in_links = []
for link in links:
    sources_in_links.append(link['source'])
if 'new' not in sources_in_links:
    sources_in_links.append('new')


# add in news now that we have populated the dict with impactful stuff
for target, source_dictionary in target_dict.items():
    for source in source_dictionary:
        value = source_dictionary[source]
        if source == 'new' and target in targets_in_links:
            links.append({'source': source, 'target': target, 'value': value})

# add in losts here too?
for target, source_dictionary in target_dict.items():
    for source in source_dictionary:
        value = source_dictionary[source]
        if value <= threshold:
            if target == 'lost' and source in sources_in_links:
                links.append({'source': source, 'target': target, 'value': value})



# populate not impactful source to not impactful target in links as other - other
other_other = []
counter_other_other = 0
for target, source_dictionary in target_dict.items():
    for source in source_dictionary:
        value = source_dictionary[source]

        if source not in sources_in_links:
            if target not in targets_in_links:
                counter_other_other+=value

if counter_other_other > 0:
    other_other.append({'source': 'other_source', 'target': 'other_target', 'value': counter_other_other})


# populate impactful source in links to other_target
source_otherTarget = []
source_otherTarget_dict = {}
for target, source_dictionary in target_dict.items():
    if target in targets_in_links:
        continue

    for source in source_dictionary:
        value = source_dictionary[source]

        if source in sources_in_links:
            if source in source_otherTarget_dict:
                source_otherTarget_dict[source] += value
            else:
                source_otherTarget_dict[source] = value

for source in source_otherTarget_dict:
    value = source_otherTarget_dict[source]
    if value > 0:
        source_otherTarget.append({'source': source, 'target': 'other_target', 'value': value})


# populate other_source in links to impactful target
otherSource_target = []
for target, source_dictionary in target_dict.items():
    if target not in targets_in_links:
        continue

    total_value=0
    for source in source_dictionary:
        if source == 'new':
            continue
        value = source_dictionary[source]
        if source in sources_in_links:
            if value <= threshold:
                total_value+=value
        else:
            total_value+=value

    if total_value > 0:
        otherSource_target.append({'source': 'other_source', 'target': target, 'value': total_value})


links = links + source_otherTarget + otherSource_target + other_other

links_sorted = sorted(links, key = lambda i: int(i['value']), reverse = True)


check_count = 0
for x in links_sorted:
    check_count += x['value']

print("The input file has " + str(len(taxa_to_lineage_new)) + " lines in it")
print("The sum of all sankey links is: " + str(check_count))
# print("The number of lost lineages is: " + str(sum([x for x in target_dict['lost']])))

with open(out_file, 'w') as f:
    f.write("source\ttarget\tvalue\n")
    for link in links_sorted:
        f.write(link['source'] + "\t" + link['target'] + "\t" + str(link['value']) + '\n')







#
