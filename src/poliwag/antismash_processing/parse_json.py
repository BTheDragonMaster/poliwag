import sys
from typing import Optional

from PIL.Image import register_open
from antismash.common.serialiser import AntismashResults
from antismash.modules.nrps_pks.results import NRPS_PKS_Results
from antismash.modules.nrps_pks.name_mappings import get_substrate_by_name
from antismash.modules.nrps_pks.results import modify_substrate as modify_substrate_as
from antismash.modules.nrps_pks.pks_names import get_short_form

from poliwag.antismash_processing.polymer import PredictedPolymer
from dataclasses import dataclass
import json

@dataclass
class Location:
    start: int
    end: int
    strand: str

def convert_domain_name(domain_name: str) -> str:
    name = domain_name.split('_')[-1].split('.')[0]
    return name

def modify_substrate(domains: list[str], module_types: list[str], base: str = "") -> str:  # pylint: disable=too-many-branches
    """ Builds a monomer including modifications from the given base.

        Arguments:
            module: the Module holding the relevant domains
            base: a string of the substate (or an empty string in case of trans-AT)

        Returns:
            the modified substrate, or an empty string if no appropriate base was
            given
    """
    domains = list(map(convert_domain_name, domains))

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
            state.append("Me")
        elif domain == "oMT":
            state.append("OMe")

    if base.endswith("mmal"):
        state.append("Me")
        base = base.replace("mmal", "mal", 1)

    if "Epimerization" in domains:
        conversions = {"Ile": "aIle", "aIle": "Ile", "Thr": "aThr", "aThr": "Thr"}
        base = conversions.get(base, base)

    state.append(base)

    if "Epimerization" in domains:
        state.insert(0, "D")
    return "-".join(state)


def parse_location(location: str) -> Location:
    if 'join' not in location:
        data = location.split(':')
        start = int(data[0][1:].strip('>').strip('<'))
        end, strand = data[1].split('](')
        end = int(end.strip('>').strip('<'))
        strand = strand[:-1]
        return Location(start, end, strand)
    else:
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


def get_modules_from_cds_name(cds_name: str, record: dict, strand: str) -> list:
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


def get_cds_from_name(cds_name, record: dict):
    for feature in record["features"]:
        if feature["type"] == "CDS":
            if cds_name in feature["qualifiers"]["locus_tag"]:
                return feature
            if "protein_id" in feature["qualifiers"] and cds_name in feature["qualifiers"]["protein_id"]:
                return feature
            if "product" in feature["qualifiers"] and cds_name in feature["qualifiers"]["product"]:
                return feature


def get_sorted_modules(record: dict) -> dict:
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
                            sorted_modules = []

                            for gene in order:
                                cds = get_cds_from_name(gene, record)
                                location = parse_location(cds["location"])
                                modules = get_modules_from_cds_name(cds["qualifiers"]["locus_tag"][0], record, location.strand)
                                for module in modules:
                                    module_string = tuple(module["qualifiers"]["domains"])

                                    if module_string not in seen_modules:
                                        sorted_modules.append(module)

                                    seen_modules.add(module_string)

                            region_to_cluster_to_modules[region_nr][cluster_number] = sorted_modules

    return region_to_cluster_to_modules


def get_modules_per_gene(record: dict) -> dict:
    gene_to_modules = {}
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


def modules_to_polymer(modules: list, record: dict, region: int, cluster: Optional[int] = None,
                       locus_tag: Optional[str] = None) -> PredictedPolymer:
    polymer = []
    for i, module in enumerate(modules):
        if "incomplete" in module["qualifiers"]:
            continue

        domains = module["qualifiers"]["domains"]
        integrated_monomer = "X"
        for domain in domains:
            if "AMP-binding" in domain:
                if "paras" in record["modules"]["antismash.modules.nrps_pks"]["domain_predictions"][domain]:
                    substrate = get_substrate_by_name(
                        record["modules"]["antismash.modules.nrps_pks"]["domain_predictions"][domain]["paras"][
                            "predicted_substrate"])
                    integrated_monomer = modify_substrate(domains, module["qualifiers"]["type"], substrate.short)
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

    return PredictedPolymer(polymer, record["id"], int(region), cluster, locus_tag)

def get_regions(record: dict) -> list:
    regions = []
    for feature in record["features"]:
        if feature["type"] == "region":
            regions.append(feature)
    return regions

def region_from_cds(cds: dict, record) -> int:
    cds_location = parse_location(cds["location"])

    for region in get_regions(record):
        region_location = parse_location(region["location"])
        if min(cds_location.start, cds_location.end) >= min(region_location.start, region_location.end) and max(cds_location.start, cds_location.end) <= max(region_location.start, region_location.end):
            return int(region["qualifiers"]["region_number"][0])

    else:
        return 0



def as_to_polymer_blocks(antismash_json_file: str) -> list:
    antismash_data = json.load(open(antismash_json_file))
    records = antismash_data['records']
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

def parse_json(antismash_json_file: str):
    antismash_data = json.load(open(antismash_json_file))
    records = antismash_data['records']
    polymers = []
    for record in records:
        region_to_cluster_to_modules = get_sorted_modules(record)
        for region, cluster_to_modules in region_to_cluster_to_modules.items():
            for cluster, modules in cluster_to_modules.items():
                polymer = modules_to_polymer(modules, record, int(region), int(cluster))
                polymers.append(polymer)
    return polymers



def polymers_from_json(antismash_json_file: str) -> list[PredictedPolymer]:
    """
    Return a list of polymers based on PARAS predictions in antiSMASH JSON output
    """
    antismash_results = AntismashResults.from_file(antismash_json_file)
    polymers = []
    for i, record in enumerate(antismash_results.records):
        result = antismash_results.results[i]
        if 'antismash.modules.nrps_pks' in result:


            nrps_annotations = NRPS_PKS_Results.from_json(result['antismash.modules.nrps_pks'], record)
            nrps_annotations.add_to_record(record)

            for region in record.get_regions():
                region_number = str(region.get_region_number())
                if "NRPS" in region.product_categories:
                    if region_number in result["antismash.modules.nrps_pks"]["region_predictions"]:
                        region_predictions = result["antismash.modules.nrps_pks"]["region_predictions"][region_number]
                    else:
                        continue


                    for region_prediction in region_predictions:

                        cluster_number = region_prediction["sc_number"]
                        order = region_prediction["ordering"]
                        polymer = []
                        seen_modules = set()
                        for gene in order:
                            cds = record.get_cds_by_name(gene)
                            if cds.location.strand == 1:
                                sorted_modules = sorted(cds.modules, key=lambda x: x.location.start)
                            else:
                                sorted_modules = sorted(cds.modules, key=lambda x: x.location.start, reverse=True)

                            for module in sorted_modules:

                                if module not in seen_modules:

                                    for domain in module.domains:

                                        if domain.domain == 'AMP-binding':

                                            substrate = nrps_annotations.domain_predictions[domain.domain_id]['paras'].get_classification()[0]
                                            monomer = modify_substrate_as(module, substrate)
                                            polymer.append(monomer)
                                seen_modules.add(module)
                        predicted_polymer = PredictedPolymer(polymer, record.id,
                                                             region.get_region_number(),
                                                             cluster_number)
                        polymers.append(predicted_polymer)


    return polymers


if __name__ == "__main__":
    print(as_to_polymer_blocks(sys.argv[1]))