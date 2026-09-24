"""
Module for antiSMASH JSON processing
"""

import sys
from typing import Optional, Any

from antismash.modules.nrps_pks.name_mappings import get_substrate_by_name
from antismash.modules.nrps_pks.pks_names import get_short_form

from poliwag.antismash_processing.polymer import PredictedPolymer
from poliwag.antismash_processing.location import Location
from poliwag.antismash_processing.domain import AdenylationDomain
import json



def convert_domain_name_to_type(domain_name: str) -> str:
    """Convert name of an AS domain to its domain type, e.g. 'nrpspksdomains_H4696_RS46465_PKS_AT.1' -> 'AT'

    :param domain_name: domain name of the form nrpspksdomains_H4696_RS46465_PKS_AT.1

    :returns: domain type of the form AT
    """
    domain_type = domain_name.split('_')[-1].split('.')[0]

    return domain_type


def domain_suffix(domain_id: str) -> str:
    """Return the trailing index of a domain id, e.g.
    'nrpspksdomains_H4696_RS46465_AMP-binding.2' -> '2'.

    This index is only unique *within a single gene* (every gene's first
    AMP-binding domain is '.1', etc.), so it must always be paired with a
    locus tag to identify a domain unambiguously. This is the same logic
    poliwag.amplicons.amp_binding_domains_to_fasta uses to build its
    AmpBindingDomain.domain_index -- kept here as the single source of truth
    so a polymer's domain provenance and an amplicon's FASTA header always
    agree on the same (locus_tag, domain_index) key.
    """
    return domain_id.rsplit(".", 1)[-1] if "." in domain_id else "1"


def modify_substrate(domains: list[str], module_types: list[str], base: str = "") -> str:  # pylint: disable=too-many-branches
    """ Builds a monomer including modifications from the given base.

        Arguments:
            domains: list of domains
            module_types: list of module types
            base: a string of the substrate (or an empty string in case of trans-AT)

        Returns:
            the modified substrate, or an empty string if no appropriate base was
            given
    """
    domains = list(map(convert_domain_name_to_type, domains))

    if "KS" in domains and "AT" not in domains:
        base = "mal"

    if not base:
        return ""

    if "pks" in module_types:
        if "KR" in domains:
            conversions = {"mal": "ohmal", "mmal": "ohmmal", "mxmal": "ohmxmal", "emal": "ohemal"}
            base = conversions.get(base, base)

        if {"DH", "DH2", "DHt"}.intersection(domains):
            conversions = {"ohmal": "ccmal", "ohmmal": "ccmmal", "ohmxmal": "ccmxmal", "ohemal": "ccemal"}
            base = conversions.get(base, base)

        if "ER" in domains:
            conversions = {"ccmal": "redmal", "ccmmal": "redmmal", "ccmxmal": "redmxmal", "ccemal": "redemal"}
            base = conversions.get(base, base)

    state: list[str] = []
    for domain in domains:

        if domain == "nMT":
            state.append("NMe")
        elif domain == "cMT":
            if "pks" in module_types:
                state.append("Me")
        elif domain == "oMT":
            state.append("OMe")

    if base.endswith("mmal"):
        state.append("Me")
        base = base.replace("mmal", "mal", 1)

    if "Epimerization" in domains:
        conversions = {"Ile": "aIle",
                       "aIle": "Ile",
                       "Thr": "aThr",
                       "aThr": "Thr",
                       "threonine": "allo-threonine",
                       "allo-threonine": "threonine",
                       "isoleucine": "allo-isoleucine",
                       "allo-isoleucine": "isoleucine"}
        base = conversions.get(base, base)

    state.append(base)

    if "Epimerization" in domains:
        state.insert(0, "D")
    return "-".join(state)


def parse_location(location: str) -> Location:
    """Return genomic location from a location string of the form [634592:635573](+) or join{[634592:635573](+), 635576:635910](+)}

    :param location: location string of the form [634592:635573](+)

    :returns Location object
    """
    if 'join' not in location:
        # Location is of the form [start:end](+)
        data = location.split(':')
        start = int(data[0][1:].strip('>').strip('<'))
        end, strand = data[1].split('](')
        end = int(end.strip('>').strip('<'))
        strand = strand[:-1]
        return Location(start, end, strand)
    else:
        # Location is of the form join{[start_1:end_1](strand), start_2:end_2](strand)}
        data = location.split('join{')[-1]
        data = data[:-1]
        locations = data.split(', ')
        overall_start = 500000000000000000000000000000
        overall_end = -10000
        overall_strand = None
        for loc in locations:
            data = loc.split(':')
            start = int(data[0][1:].strip('>').strip('<'))
            end, strand = data[1].split('](')
            end = int(end.strip('>').strip('<'))
            strand = strand[:-1]
            if overall_strand is not None:
                assert overall_strand == strand
            else:
                overall_strand = strand
            if end > overall_end:
                overall_end = end
            if start < overall_start:
                overall_start = start

        assert overall_strand is not None

        return Location(overall_start, overall_end, overall_strand)


def get_modules_from_cds_name(cds_name: str, record: dict, strand: str) -> list[dict[str, Any]]:
    """
    Obtain list of modules from a CDS name

    :param cds_name: name of the CDS
    :param record: antiSMASH record in JSON format
    :param strand: + or -

    :return: list of antiSMASH modules
    """
    module_features: list = []
    for feature in record["features"]:

        if feature["type"] == "aSModule":

            if cds_name in feature["qualifiers"]["locus_tags"]:
                module_features.append(feature)

    if strand == '+':
        module_features.sort(key = lambda x: parse_location(x["location"]).start)
    else:
        module_features.sort(key=lambda x: parse_location(x["location"]).start, reverse = True)

    return module_features


def get_cand_clusters_from_location(location: Location, record: dict[str, Any]) -> list[int]:
    """Get candidate cluster numbers from a location

    :param location: genomic location
    :param record: antiSMASH record in JSON format
    :return: list of antiSMASH candidate cluster numbers"""

    cand_clusters: list = []

    for feature in record["features"]:
        if feature["type"] == "cand_cluster":
            cand_location = parse_location(feature["location"])
            if location.is_contained_in(cand_location):
                cand_clusters.append(int(feature["qualifiers"]["candidate_cluster_number"][0]))

    return cand_clusters






def get_cds_from_name(cds_name: str, record: dict[str, Any]) -> dict[str, Any] | None:
    """
    Get CDS feature from a CDS name

    :param cds_name: name of the CDS
    :param record: antiSMASH record in JSON format

    :return: CDS feature
    """
    for feature in record["features"]:
        if feature["type"] == "CDS":
            if cds_name in feature["qualifiers"]["locus_tag"]:
                return feature
            if "protein_id" in feature["qualifiers"] and cds_name in feature["qualifiers"]["protein_id"]:
                return feature
            if "product" in feature["qualifiers"] and cds_name in feature["qualifiers"]["product"]:
                return feature

    return None


def get_sorted_modules(record: dict) -> dict[str, Any]:
    """Retrieve sorted modules from an antiSMASH record in JSON format.

    Each module is paired with the locus tag of the gene it belongs to, as
    ``(module, gene_locus_tag)``, since a module/domain's own index (e.g.
    'AMP-binding.1') is only unique within a single gene -- see
    ``modules_to_polymer`` and ``domain_suffix``.
    """

    region_to_cluster_to_modules = {}
    if "modules" in record:
        if "antismash.modules.nrps_pks" in record["modules"]:
            if "region_predictions" in record["modules"]["antismash.modules.nrps_pks"]:
                for region_nr, region_predictions in record["modules"]["antismash.modules.nrps_pks"]["region_predictions"].items():
                    region_to_cluster_to_modules[region_nr] = {}
                    for region_prediction in region_predictions:
                        cluster_number = region_prediction["sc_number"]
                        order = region_prediction["ordering"]
                        if order:

                            seen_modules = set()
                            sorted_modules: list[tuple[dict[str, Any], str]] = []

                            for gene in order:
                                cds = get_cds_from_name(gene, record)
                                gene_locus_tag = cds["qualifiers"]["locus_tag"][0]
                                location = parse_location(cds["location"])
                                modules = get_modules_from_cds_name(gene_locus_tag, record, location.strand)
                                for module in modules:
                                    module_string = tuple(module["qualifiers"]["domains"])

                                    if module_string not in seen_modules:
                                        sorted_modules.append((module, gene_locus_tag))

                                    seen_modules.add(module_string)

                            region_to_cluster_to_modules[region_nr][cluster_number] = sorted_modules

    return region_to_cluster_to_modules


def get_modules_per_gene(record: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    """
    Retrieve a mapping of genes to their modules from an antiSMASH record

    :param record: antiSMASH record in JSON format

    :return: mapping of genes to their modules
    """
    gene_to_modules: dict[str, list[dict[str, Any]]] = {}
    for feature in record["features"]:
        if feature["type"] == "CDS":
            modules = get_modules_from_cds_name(feature["qualifiers"]["locus_tag"][0], record, parse_location(feature["location"]).strand)
            gene_to_modules[feature["qualifiers"]["locus_tag"][0]] = modules
    return gene_to_modules

def get_gene_to_region(record: dict) -> dict:
    gene_to_region = {}
    for feature in record["features"]:
        if feature["type"] == "CDS":
            region = region_from_cds(feature, record)
            gene_to_region[feature["qualifiers"]["locus_tag"][0]] = region

    return gene_to_region


def modules_to_polymer(modules: list[dict[str, Any]],
                       record: dict,
                       region: int,
                       cluster: Optional[int] = None,
                       locus_tag: Optional[str] = None,
                       module_locus_tags: Optional[list[str]] = None) -> PredictedPolymer:
    """
    Obtain a building block polymer from a list of modules, in assembly-line order

    :param modules: list of antiSMASH modules
    :param record: antiSMASH record in JSON format
    :param region: region number
    :param cluster: candidate cluster number
    :param locus_tag: locus tag of the modules, when every module in ``modules``
        belongs to the same gene (the per-gene extraction path). Ignored for a
        given module if that module has an entry in ``module_locus_tags``.
    :param module_locus_tags: locus tag for each entry in ``modules``, when
        ``modules`` spans more than one gene (the per-region extraction path,
        as produced by ``get_sorted_modules``). Must be the same length as
        ``modules`` if given.
    """
    polymer = []
    domain_ids: list[Optional[tuple[str, str]]] = []
    for i, module in enumerate(modules):
        if "incomplete" in module["qualifiers"]:
            continue

        module_locus_tag = module_locus_tags[i] if module_locus_tags is not None else locus_tag

        domains = module["qualifiers"]["domains"]
        integrated_monomer = "X"
        domain_id: Optional[tuple[str, str]] = None
        for domain in domains:
            if "AMP-binding" in domain:
                if "paras" in record["modules"]["antismash.modules.nrps_pks"]["domain_predictions"][domain]:
                    substrate = get_substrate_by_name(
                        record["modules"]["antismash.modules.nrps_pks"]["domain_predictions"][domain]["paras"][
                            "predicted_substrate"])
                    integrated_monomer = modify_substrate(domains, module["qualifiers"]["type"], substrate.long)
                    if module_locus_tag is not None:
                        domain_id = (module_locus_tag, domain_suffix(domain))
            elif "PKS_AT" in domain:
                substrate = get_short_form(
                    record["modules"]["antismash.modules.nrps_pks"]["domain_predictions"][domain]["minowa_at"][
                        "predictions"][0][0])
                integrated_monomer = modify_substrate(domains, module["qualifiers"]["type"], substrate)
            elif "CAL" in domain:
                integrated_monomer = \
                record["modules"]["antismash.modules.nrps_pks"]["domain_predictions"][domain]["minowa_cal"][
                    "predictions"][0][0]

        polymer.append(integrated_monomer)
        domain_ids.append(domain_id)

    return PredictedPolymer(polymer, record["id"], int(region), cluster, locus_tag, domain_ids)

def get_regions(record: dict) -> list[dict[str, Any]]:
    """
    Obtain regions from antiSMASH record

    :param record: antiSMASH record in JSON format
    :return: list of regions
    """
    regions = []
    for feature in record["features"]:
        if feature["type"] == "region":
            regions.append(feature)
    return regions

def region_from_cds(cds: dict[str, Any], record: dict[str, Any]) -> int:
    """
    Obtain region number from CDS
    :param cds: CDS in antiSMASH JSON format
    :param record: antiSMASH record in JSON format
    :return: region number
    """
    cds_location = parse_location(cds["location"])

    for region in get_regions(record):
        region_location = parse_location(region["location"])
        if min(cds_location.start, cds_location.end) >= min(region_location.start, region_location.end) and max(cds_location.start, cds_location.end) <= max(region_location.start, region_location.end):
            return int(region["qualifiers"]["region_number"][0])

    else:
        return 0


def _as_to_polymers_per_gene(records: list[dict[str, Any]]) -> list[str]:
    """
    Return per-gene predicted NRPS/PKS polymers from a list of antiSMASH records

    :param records: list of antiSMASH records in JSON format
    :return: list of NRPS/PKS polymers, with each polymer representing the building blocks incorporated by a single gene
    """
    polymers = []
    for record in records:

        gene_to_region = get_gene_to_region(record)
        gene_to_modules = get_modules_per_gene(record)
        for gene, modules in gene_to_modules.items():
            if modules:
                region = gene_to_region[gene]
                if region == 0:
                    print(f'Warning: no regions found for CDS {gene}. Setting region number to 0')
                polymer = modules_to_polymer(modules, record, int(region), locus_tag = gene)
                if polymer.polymer:
                    polymers.append(polymer)
    return polymers

def _as_to_polymers_per_region(records: list[dict[str, Any]]):
    """
    Return per-region predicted NRPS/PKS polymers from a list of antiSMASH records

    :param records: list of antiSMASH records in JSON format
    :return: list of NRPS/PKS polymers, with each polymer representing the building blocks incorporated by the entire region
    """
    polymers = []
    for record in records:
        region_to_cluster_to_modules = get_sorted_modules(record)
        for region, cluster_to_modules in region_to_cluster_to_modules.items():
            for cluster, modules_with_loci in cluster_to_modules.items():
                modules = [module for module, _ in modules_with_loci]
                module_locus_tags = [gene_locus_tag for _, gene_locus_tag in modules_with_loci]
                polymer = modules_to_polymer(
                    modules, record, int(region), int(cluster), module_locus_tags=module_locus_tags
                )
                polymers.append(polymer)
    return polymers

def as_to_polymers(antismash_json_file: str, per_gene: bool) -> list[str]:
    """
    Return predicted NRPS/PKS polymers from antiSMASH JSON format

    :param antismash_json_file: antiSMASH JSON output file
    :param per_gene: whether to return predicted NRPS/PKS polymers per gene or per region
    :return: list of NRPS/PKS polymers
    """

    antismash_data = json.load(open(antismash_json_file))
    records = antismash_data['records']
    if per_gene:
        polymers = _as_to_polymers_per_gene(records)
    else:
        polymers = _as_to_polymers_per_region(records)

    return polymers

def as_to_a_domains(antismash_json_file: str) -> list[AdenylationDomain]:
    """
    Write antiSMASH A domains to file
    """
    antismash_data = json.load(open(antismash_json_file))
    records = antismash_data['records']
    a_domains: list[AdenylationDomain] = []

    for record in records:
        for feature in record["features"]:
            if feature["type"] == "aSDomain":
                domain_ids = feature["qualifiers"]["domain_id"]
                if domain_ids:
                    domain_id = domain_ids[0]
                    if 'AMP-binding' in domain_id:
                        location = parse_location(feature["location"])
                        translation = feature["qualifiers"]["translation"][0]
                        locus_tag = feature["qualifiers"]["locus_tag"][0]

                        protein_start = int(feature["qualifiers"]["protein_start"][0])
                        protein_end = int(feature["qualifiers"]["protein_end"][0])

                        cds = get_cds_from_name(locus_tag, record)

                        if cds is not None:
                            region = region_from_cds(cds, record)
                            cand_clusters = get_cand_clusters_from_location(location, record)



                            domain_predictions = record["modules"]["antismash.modules.nrps_pks"]["domain_predictions"]

                            paras_substrate = None
                            paras_confidence = None

                            if "paras" in domain_predictions.get(domain_id, {}):
                                paras_substrate = domain_predictions[domain_id]["paras"]["predicted_substrate"]
                                paras_confidence = domain_predictions[domain_id]["paras"]["confidence"]

                            a_domain = AdenylationDomain(
                                translation,
                                location,
                                region,
                                cand_clusters,
                                locus_tag,
                                protein_start,
                                protein_end,
                                paras_substrate,
                                paras_confidence
                            )
                            a_domains.append(a_domain)
    return a_domains


if __name__ == "__main__":
    print(as_to_a_domains(sys.argv[1]))
