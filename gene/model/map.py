import json
import os


class Map(object):

    def __init__(self, user_id, pid):
        self.user_id = user_id
        self.pid = pid
        self.project_name = ""
        self.orig_otu_table_name = ""
        self.orig_sample_metadata_name = ""
        self.taxonomy_type = ""
        self.matrix_type = "int"
        self.num_samples = 0
        self.num_otus = 0
        self.shared = "no"
        self.load()

    def load(self):
        if os.path.exists(self.__get_map_path()):
            with open(self.__get_map_path()) as f:
                map_from_json = json.load(f)
                self.project_name = map_from_json["project_name"]
                self.orig_otu_table_name = map_from_json["orig_otu_table_name"]
                self.orig_sample_metadata_name = map_from_json["orig_sample_metadata_name"]
                self.taxonomy_type = map_from_json["taxonomy_type"]
                self.matrix_type = map_from_json["matrix_type"] if "matrix_type" in map_from_json else "int"
                self.num_samples = map_from_json["num_samples"] if "num_samples" in map_from_json else 0
                self.num_otus = map_from_json["num_otus"] if "num_otus" in map_from_json else 0
                self.shared = map_from_json["shared"] if "shared" in map_from_json else "no"

    def save(self):
        map_from_json = {
            "project_name": self.project_name,
            "orig_otu_table_name": self.orig_otu_table_name,
            "orig_sample_metadata_name": self.orig_sample_metadata_name,
            "taxonomy_type": self.taxonomy_type,
            "matrix_type": self.matrix_type,
            "num_samples": self.num_samples,
            "num_otus": self.num_otus,
            "shared": self.shared
        }

        with open(self.__get_map_path(), 'w') as outfile:
            json.dump(map_from_json, outfile)

    def __get_map_path(self):
        project_dir = os.path.dirname(__file__)
        project_dir = os.path.abspath(os.path.join(project_dir, os.pardir))  # Gets the parent folder
        project_dir = os.path.join(project_dir, "data")
        project_dir = os.path.join(project_dir, self.user_id)
        project_dir = os.path.join(project_dir, self.pid)
        map_path = os.path.join(project_dir, "map.txt")
        return map_path