## CoG-UK Lineage Designation

### What are UK lineages?

UK lineages are groups of sequences which roughly approximate transmission chains after introductions from abroad. 

Due to the relatively low rate of evolution of SARS-CoV-2, it is very difficult to differentiate between groups of sequences that may or may not be from different introductions into the UK, as the earliest sequences in the lineage are identical to each other. The intention therefore is less about counting the precise number of introductions into the UK and more about monitoring the spread of the virus once it's there.

Some lineages are geographically localised within one of the four nation states of the UK, and some others are found across the UK. Most lineages are small, many with only one or a few samples. Others are much larger, with the largest lineage consisting of over 1000 sequences sampled from across the UK.   

### The Pipeline

#### Ancestral reconstruction:

For the whole phylogenetic tree we:

1. Collapse all branches with length < 5E-6. These result from differences in ambiguities among the sequences.

2. Create a binary trait from the location data to distinguish UK from non-UK samples

3. Reconstruct the ancestral states of this trait using the deltrans implementation of the Fitch algorithm.

	(NB: This would be valid on a bifurcating tree, but our trees are rife with polytomies.  For a valid parsimonious reconstruction, we would have to resolve these polytomies in some way, and would require estimating the number of UK introductions. In a way that’s sort of the whole point of the analysis. Everything that follows, results from us trying come up with reasonable heuristics for interpreting these reconstructions)
	
4. For deltrans ancestral states we traverse the tree and identify where the location changes from non-UK to UK. These are labeled as introductions and have the caveat that we label any descendent export and subsequent re-introduction as the ancestral introduction.

5. We then group introductions together by traversing the tree in a level-order traversal starting at the tips. Any introductions that stem from a polytomy are clustered together and given a `del_lineage` identifier.  This is allowed to happen once on every path to the root.  This is an arbitrary threshold, but it groups recent introductions together without combining all introductions into one lineage.

#### UK lineage (re-)naming

Taking as input:

1. A table of sequence names
2. The previous week's UK lineage designations 
3. This week's deltrans designations

A **`UK Lineage`** is defined as a relabeling of del_lineage that tries to keep labels as consistent as possible from week to week to aid with long term analyses.

Due to uncertainty in the tree, sequences are often reassigned to new deltrans designations, and lineages split and merge with each new dataset. Therefore at this point, lineage names are reassigned so that **one deltrans designation relates to one lineage**. 

The largest lineage name from the previous week becomes the new name for the lineage. 
If it is the largest in more than one new designation, the designation with the most sequences of that lineage claims the name. In this way the names stay relatively stable from week to week.

Example scenario: 
Group A contains 10 UK7 sequences and 5 UK2 sequences. Group B contains 8 UK7 sequences and 3 UK9 sequences. No other groups contain UK7 sequences. Group A claims the UK7 name, as it is the largest lineage in group A, and it also contains more UK7 sequences than group B. 
For group B, we then assess its suitability to be named UK9.  However, group C contains 100 UK9 sequences and 10 UK2 sequences. Therefore group C claims the name UK9. Group B has no more old lineages represented, so it acquires an entirely new name.

Sequences with curated lineage assignments and associated metadata are then passed into the report generator. This summarises the lineages to obtain relevant statistics based on size, timing and location, and produces figures for the whole of the UK, each constituent country, and sequencing centres.




#### Phylotyping




