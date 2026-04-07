from sys import argv


from antismash.common.secmet.record import Record
from antismash.modules.nrps_pks import run_on_record
from antismash.config import get_config


def get_gene_order(antismash_region_genbank: str, candidate_cluster_type: str = 'NRPS') -> list[str]:
    gene_orders = []
    records = Record.from_genbank(antismash_region_genbank)
    for record in records:

        nrps_pks_results = run_on_record(record, None, get_config())

        for domain_prediction in nrps_pks_results.domain_predictions:
            print(domain_prediction)

        for region, predictions in nrps_pks_results.region_predictions.items():

            for prediction in predictions:
                candidate_cluster = record.get_candidate_cluster(prediction.candidate_cluster_number)
                if candidate_cluster_type in candidate_cluster.product_categories and len(candidate_cluster.product_categories) == 1:
                    gene_orders.append(prediction.ordering)

    return gene_orders

def get_module_order(antismash_region_genbank: str, candidate_cluster_type: str = 'NRPS') ->

def get_modules_from_gbk():

    pass


if __name__ == "__main__":
    gene_orders = get_gene_order(argv[1])
    print(gene_orders)