Scientific concepts
===================

Why KmerSutra exists
--------------------

The motivating Plasmodium benchmark showed that analytical sensitivity and
species-level interpretability are not the same outcome. Kraken2 and Metabuli
detected low-abundance targets sensitively, but their positive reports also
contained multiple non-expected Plasmodium labels. KmerSutra was developed to
test a more conservative approach in which every species-level conclusion can
be traced back to taxonomically contextualised marker evidence.

This does not mean that general-purpose classifiers perform poorly. They solve
a broader read-assignment problem and occupy a different operating point.

Outgroup-aware markers
----------------------

A candidate k-mer is assessed against the full configured reference context.
If it is shared beyond one species, it is downgraded to the most specific
supported taxonomic rank rather than being forced into a species call. Adding
near neighbours and broader outgroups therefore changes the interpretation of
markers, not merely the size of a lookup table.

Multiple k values
-----------------

The standard ladder is ``51, 77, 101, 151``. Shorter k values offer greater
tolerance of sequence divergence but may be less specific. Longer values offer
stronger specificity but are more vulnerable to sequencing error and
held-out-strain differences. KmerSutra retains the evidence contributed by each
k and can require cross-k support at the report layer.

Raw evidence versus reportable detection
----------------------------------------

KmerSutra deliberately separates:

1. observed marker evidence;
2. evidence passing numerical thresholds;
3. rank-aware lineage interpretation;
4. the final reportable species call.

A true organism can therefore have substantial raw evidence while remaining
``neighbour_lineage_evidence`` under a conservative report rule. That is not
equivalent to “no detection”; it is evidence that the current panel and
thresholds do not justify an unqualified species label.

Optional calibration
--------------------

The frozen calibration model operates on summarised evidence features. It is
optional, auditable and downstream of deterministic calls. It does not perform
per-read neural classification and must not be described as a replacement for
the rule-based framework.
