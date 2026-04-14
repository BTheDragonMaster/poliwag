from sys import argv

from antismash.common.serialiser import AntismashResults
from antismash.common.secmet.features import CDSFeature
from antismash.modules.nrps_pks.results import NRPS_PKS_Results, modify_substrate
from poliwag.utils import iterate_over_dir

class PredictedPolymer:
    def __init__(self, polymer: list[str],
                 accession: str,
                 region_number: int,
                 candidate_cluster: str):
        self.polymer = polymer
        self.accession = accession
        self.region = region_number
        self.candidate_cluster = candidate_cluster

    def __repr__(self):
        return ' | '.join(self.polymer)

    def __str__(self):
        return ' | '.join(self.polymer)


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
                    region_predictions = result["antismash.modules.nrps_pks"]["region_predictions"][region_number]

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
                                            monomer = modify_substrate(module, substrate)
                                            polymer.append(monomer)
                                seen_modules.add(module)
                        predicted_polymer = PredictedPolymer(polymer, record.id,
                                                             region.get_region_number(),
                                                             cluster_number)
                        polymers.append(predicted_polymer)


    return polymers

def write_polymers(antismash_output_folder: str, out_file: str):
    with open(out_file, 'w') as out:
        for _, folder_path in iterate_over_dir(antismash_output_folder, get_dirs=True):
            for genome_name, file_path in iterate_over_dir(folder_path, '.json'):
                try:
                    polymers = polymers_from_json(file_path)
                    for polymer in polymers:
                        out.write(f"{polymer.accession}\t{polymer.region}\t{polymer.candidate_cluster}\t{polymer}\n")
                except Exception as e:
                    print(f"Could not decode input file {genome_name}")
                    print(e)

def read_polymers(polymers_file: str) -> list[PredictedPolymer]:
    predicted_polymers = []
    with open(polymers_file, 'r') as polymers:
        for line in polymers:
            accession, region, candidate, polymer = line.split('\t')
            polymer = PredictedPolymer(polymer.strip().split(' | '), accession, int(region), candidate)
            predicted_polymers.append(polymer)

    return predicted_polymers


def get_unique_polymers(predicted_polymers: list[PredictedPolymer]) -> list[str]:
    unique_polymers = set()
    for polymer in predicted_polymers:
        unique_polymers.add(str(polymer))

    unique_polymers = list(unique_polymers)
    unique_polymers.sort()

    return unique_polymers

if __name__ == '__main__':
    # write_polymers(argv[1], argv[2])
    # print(polymers_from_json(argv[1]))
    for polymer in get_unique_polymers(read_polymers(argv[1])):
        print(polymer)
