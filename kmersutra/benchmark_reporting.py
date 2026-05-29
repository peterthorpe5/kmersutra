"""Benchmark-aware reporting helpers for KmerSutra summaries.

These helpers sit above raw KmerSutra evidence. They are intended for
controlled spike-in benchmarks where the expected target species are known.
The functions keep raw taxa auditable while allowing summary reports to
separate expected targets, empirical-background candidates, expected-lineage
neighbours, and true reportable off-target species.
"""

from __future__ import annotations

import re
from collections.abc import Sequence


def normalise_taxon_name(*, value: object) -> str:
    """Normalise a taxon label for robust equality checks.

    Parameters
    ----------
    value : object
        Raw taxon label.

    Returns
    -------
    str
        Lower-case taxon label with underscores replaced by spaces and repeated
        whitespace collapsed.
    """
    if value is None:
        return ""
    text = str(value).replace("_", " ").strip().lower()
    text = re.sub(r"\s+", " ", text)
    return text


def extract_genus(*, taxon_name: object) -> str:
    """Extract a normalised genus token from a taxon label.

    Parameters
    ----------
    taxon_name : object
        Raw taxon label.

    Returns
    -------
    str
        First whitespace-delimited token, or an empty string when unavailable.
    """
    label = normalise_taxon_name(value=taxon_name)
    if not label:
        return ""
    return label.split(" ", maxsplit=1)[0]


def normalise_taxa(*, taxa: Sequence[object] | None) -> set[str]:
    """Normalise a sequence of taxon labels.

    Parameters
    ----------
    taxa : sequence of object or None
        Taxon labels.

    Returns
    -------
    set[str]
        Non-empty normalised labels.
    """
    return {
        normalise_taxon_name(value=taxon)
        for taxon in (taxa or [])
        if normalise_taxon_name(value=taxon)
    }


def expected_genera_from_targets(*, expected_targets: Sequence[object] | None) -> set[str]:
    """Return genera represented among the expected benchmark targets.

    Parameters
    ----------
    expected_targets : sequence of object or None
        Expected species labels for one benchmark sample.

    Returns
    -------
    set[str]
        Non-empty genus tokens.
    """
    return {
        extract_genus(taxon_name=target)
        for target in (expected_targets or [])
        if extract_genus(taxon_name=target)
    }


def is_expected_genus_neighbour(
    *,
    report_label: object,
    expected_targets: Sequence[object] | None,
    is_expected_target: bool,
    is_negative_sample: bool,
    is_background_candidate: bool,
    demote_expected_genus_neighbours: bool = True,
) -> bool:
    """Return whether a taxon should be treated as expected-lineage evidence.

    This function is benchmark-aware. It should only demote a non-expected
    same-genus species in positive spike-in samples where the expected species
    are known. Negative/no-spike samples are deliberately not demoted, because
    same-genus signal there is unexpected and should remain visible.

    Parameters
    ----------
    report_label : object
        Taxon label being evaluated.
    expected_targets : sequence of object or None
        Expected species labels for the sample.
    is_expected_target : bool
        Whether ``report_label`` is itself an expected target.
    is_negative_sample : bool
        Whether the benchmark sample is a no-spike/shuffled negative sample.
    is_background_candidate : bool
        Whether the taxon has already been classified as empirical-background
        candidate signal.
    demote_expected_genus_neighbours : bool, optional
        Whether benchmark-aware expected-genus demotion is enabled.

    Returns
    -------
    bool
        True when the taxon should be counted as expected-lineage neighbour
        evidence rather than as a strict off-target species.
    """
    if not demote_expected_genus_neighbours:
        return False
    if is_negative_sample or is_expected_target or is_background_candidate:
        return False
    expected_genera = expected_genera_from_targets(expected_targets=expected_targets)
    if not expected_genera:
        return False
    return extract_genus(taxon_name=report_label) in expected_genera


def reporting_layer_for_call(
    *,
    is_positive_call: bool,
    is_species_level: bool,
    is_expected_target: bool,
    is_background_candidate: bool,
    is_expected_genus_neighbour_call: bool,
) -> str:
    """Assign a benchmark-level report layer to one taxon call.

    Parameters
    ----------
    is_positive_call : bool
        Whether the row is a positive KmerSutra call.
    is_species_level : bool
        Whether the row represents species-level evidence.
    is_expected_target : bool
        Whether the row is one of the expected benchmark targets.
    is_background_candidate : bool
        Whether the row is empirical-background candidate signal.
    is_expected_genus_neighbour_call : bool
        Whether the row is expected-genus neighbouring-lineage evidence.

    Returns
    -------
    str
        One of ``expected_target``, ``background_candidate``,
        ``expected_genus_neighbour``, ``reportable_off_target``,
        ``non_species_context`` or ``not_positive``.
    """
    if not is_positive_call:
        return "not_positive"
    if not is_species_level:
        return "non_species_context"
    if is_expected_target:
        return "expected_target"
    if is_background_candidate:
        return "background_candidate"
    if is_expected_genus_neighbour_call:
        return "expected_genus_neighbour"
    return "reportable_off_target"
