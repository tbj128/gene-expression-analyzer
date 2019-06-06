from gene.core.constants import RAW_GENE_TABLE_LABELS_FILENAME
from gene.core.data_io import DataIO


class Taxonomy(object):
    """
    TODO: Use GO annotation table
    """

    OTU_COL = 0
    TAXONOMY_COL = 1

    def __init__(self, user_id, pid):
        self.user_id = user_id
        self.pid = pid
        self.taxonomy_map = {}
        self.__load_taxonomy()

    def __load_taxonomy(self):
        tax = DataIO.tsv_to_table(self.user_id, self.pid, RAW_GENE_TABLE_LABELS_FILENAME)
        headers = tax[0]
        self.taxonomy_map = self.__get_taxonomy_mapping_from_dict(headers)

    def get_taxonomy_map(self):
        return self.taxonomy_map

    def __get_taxonomy_mapping_from_dict(self, headers):
        taxonomyMap = {}
        for header in headers:
            taxonomyMap[header] = header
        return taxonomyMap
