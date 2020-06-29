##CoG-UK Lineage Designation

### What are UK lineages?

UK lineages are groups of sequences which roughly approximate transmission chains after introductions from abroad. 

Due to the relatively low rate of evolution of SARS-CoV-2, it is very difficult to differentiate between groups of sequences that may or may not be from different introductions into the UK, as the earliest sequences in the lineage are identical to each other. The intention therefore is less about counting the precise number of introductions into the country and more about monitoring the spread of the virus internally.

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

a **`UK Lineage`** is defined as...


#### Phylotyping




