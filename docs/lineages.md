## COG-UK Lineage Designation

### What are UK lineages?

UK lineages are groups of genetically similar sequences found in the UK. Groups of UK sequences that share a mutation with any samples found outside the UK will be classified as a distinct lineage. 

Due to the relatively low rate of evolution of SARS-CoV-2, it is very difficult to differentiate between groups of sequences that may or may not be from different introductions into the UK, as the earliest sequences in the lineage are identical to each other. It is possible (and likely common) for samples from distinct introductions to be grouped together. It is also possible for samples from the same introduction to be found in separate UK lineages. This happens when a UK transmission chain leads to an exportation where the exported sample shares a mutation with only some of the UK sequences. 

Some lineages are geographically localised within one of the four nation states of the UK, and some others are found across the UK. Most lineages are small, many with only one or a few samples. Others are much larger, with the largest lineage consisting of over 1000 sequences sampled from across the UK.   

### The Pipeline

#### Ancestral reconstruction:

For the whole phylogenetic tree we:

1. Collapse all branches with length < 5E-6. These branches result from differences in ambiguities among the sequences.

2. Create a binary trait from the location data to distinguish UK from non-UK samples

3. Reconstruct the ancestral states of this trait using the DELTRANS implementation of the Fitch algorithm.

	(NB: This would be valid on a bifurcating tree, but our trees are rife with polytomies.  For a valid parsimonious reconstruction, we would have to resolve these polytomies in some way, and would require estimating the number of UK introductions. Everything that follows, results from us trying come up with reasonable heuristics for interpreting these reconstructions)
	
4. For DELTRANS ancestral states we traverse the tree and identify where the location changes from non-UK to UK. These are labeled as introductions and have the caveat that we label any descendent export and subsequent re-introduction as the ancestral introduction.

5. We then group introductions together by traversing the tree in a level-order traversal starting at the tips. Any sibling introductions are clustered together and given a `del_lineage` identifier. This ensures that identical samples are given the same UK lineage.

6. A **`UK Lineage`** is a relabeling of `del_lineage` that tries to keep labels as consistent as possible from week to week to aid with long term analyses (see below).


#### UK lineage (re-)naming

Taking as input:

1. A table of sample names with:
2. Their previous UK lineage designations, and
3. Their new DELTRANS designations

Due to uncertainty in the tree and the addition of data from outside the UK, sequences are often reassigned to new deltrans designations, and lineages split and merge with each new dataset. Therefore at this point, lineage names are reassigned so that **one DELTRANS designation relates to one UK lineage**. 

The most-represented UK lineage name (by number of samples) from the previous week that maps to each new DELTRANS designation becomes the new name for that lineage. 

If one previous UK lineage is the most-represented lineage in more than one new designation, the designation with the most sequences of that lineage claims the name. In this way the names stay relatively stable from week to week.

**Example scenario:** 
Group A contains 10 UK7 sequences and 5 UK2 sequences. Group B contains 8 UK7 sequences and 3 UK9 sequences. No other groups contain UK7 sequences. Group A claims the UK7 name, as it is the largest lineage in group A, and it also contains more UK7 sequences than group B. 
For group B, we then assess its suitability to be named UK9.  However, group C contains 100 UK9 sequences and 10 UK2 sequences. Therefore group C claims the name UK9. Group B has no more old lineages represented, so it acquires an entirely new name.

Sequences with curated lineage assignments and associated metadata are then passed into the report generator. This summarises the lineages to obtain relevant statistics based on size, timing and location, and produces figures for the whole of the UK, each constituent country, and sequencing centres.


#### Phylotyping
Phylotypes codify the internal nodes on the tree. Samples that have the same phylotypes are allowed differ from each other only by unique mutations. Samples in different phylotypes share mutations with other sequences in the data set. An exception to this rule is that samples with different phylotypes may share a mutation if that mutation is a homoplasy (i.e. it arose multiple times on the tree). Phylotypes are defined for each UK lineage, and a UK lineage can include multiple phylotypes. The UK lineage definition respects phylotypes in that any samples in different UK lineages would necessarily be in different phylotypes.



