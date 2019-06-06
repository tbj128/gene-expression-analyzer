# ==============================================================================
#
# Utility functions used for data transformation or other common functionality
# @author: tbj128
#
# ==============================================================================

#
# Imports
#

import uuid

from io import StringIO
import numpy as np

from gene.core.data_io import DataIO
from gene.core.constants import RAW_GENE_TABLE_FILENAME, \
    SAMPLE_METADATA_FILENAME, \
    RAW_GENE_TABLE_LABELS_FILENAME
import csv
import os
import re
import logging
import shutil

from gene.model.map import Map
from gene.model.otu_table import OTUTable
from gene.model.user_request import UserRequest

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

OK = 0
GENERAL_ERROR = -1
BIOM_ERROR = -2
OTU_ERROR = -3
SAMPLE_METADATA_ERROR = -5
OTU_LABEL_NOT_IN_SAMPLE_METADATA_ERROR = -7
OTU_DATATYPE_ERROR = -8

class ProjectManager(object):
    BASE_DIRECTORY = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))  # Gets the parent folder
    DATA_DIRECTORY = os.path.join(BASE_DIRECTORY, "data")
    STAGING_DIRECTORY = os.path.join(DATA_DIRECTORY, "staging")

    def __init__(self, user_id):
        self.user_id = user_id

    def get_file_for_download(self, project_name, type):
        if type == "sample_metadata":
            return DataIO.tsv_to_table(self.user_id, project_name, SAMPLE_METADATA_FILENAME)
        elif type == "otu":
            table = DataIO.tsv_to_table(self.user_id, project_name, RAW_GENE_TABLE_FILENAME)
            labels = DataIO.tsv_to_table(self.user_id, project_name, RAW_GENE_TABLE_LABELS_FILENAME)
            new_headers = ["Sample Labels"]
            new_headers.extend(labels[0])
            full_table = [new_headers]
            i = 0
            while i < len(table):
                new_row = [labels[1][i] if i < len(labels[1]) else ""]
                new_row.extend(table[i])
                full_table.append(new_row)
                i += 1
            return full_table
        else:
            return []

    def get_filtering_info(self, pid, sampleFilter, sampleFilterRole, sampleFilterVals):
        """
        Returns information that will tell the user what samples will be removed and what the subsample value would be
        :param pid:
        :param sampleFilter:
        :param sampleFilterRole:
        :param sampleFilterVals:
        :return:
        """
        user_request = UserRequest(self.user_id, pid, "", "", "",
                                   "", [], sampleFilter, sampleFilterRole,
                                   sampleFilterVals, 0, "")
        map = Map(self.user_id, pid)

        table = OTUTable(self.user_id, pid, use_raw=True)
        orig_base = table.get_table()
        orig_headers = table.get_headers()
        orig_sample_labels = table.get_sample_labels()
        base, headers, sample_labels = table.filter_otu_table_by_metadata(orig_base, orig_headers, orig_sample_labels, user_request)
        initial_samples_removed = set(orig_sample_labels) - set(sample_labels)

        has_float = map.matrix_type == "float"

        orig_base = np.array(orig_base)
        orig_row_sums = orig_base.sum(axis=1)


        samples = {}
        i = 0
        while i < len(orig_base):
            row_sum = orig_row_sums[i]
            samples[orig_sample_labels[i]] = {
                "row_sum": row_sum,
                "removed": orig_sample_labels[i] in initial_samples_removed
            }
            i += 1

        return {
            "samples": samples,
            "has_float": has_float
        }


    def stage_project_from_tsv(self, project_name, otu_filename, sample_metadata_filename):
        # Creates a directory for this project
        pid = str(uuid.uuid4())
        project_dir = os.path.join(ProjectManager.DATA_DIRECTORY, self.user_id)
        project_dir = os.path.join(project_dir, pid)
        if not os.path.exists(project_dir):
            os.makedirs(project_dir)
        else:
            logger.exception("Cannot create project folder")
            raise Exception("Cannot create project folder as it already exists")

        # Renames the uploaded files to a standard file schema and moves to the project directory
        user_staging_dir = os.path.join(ProjectManager.STAGING_DIRECTORY, self.user_id)
        os.rename(os.path.join(user_staging_dir, sample_metadata_filename),
                  os.path.join(project_dir, SAMPLE_METADATA_FILENAME))

        sample_ids_from_sample_metadata = {}
        sample_metadata = DataIO.tsv_to_table(self.user_id, pid, SAMPLE_METADATA_FILENAME, accept_empty_headers=False)

        i = 0
        while i < len(sample_metadata):
            if i > 0:
                if len(sample_metadata[i]) > 0:
                    sample_ids_from_sample_metadata[sample_metadata[i][0]] = 1
            i += 1

        logger.info("Beginning to load the OTU table")
        base_arr = DataIO.tsv_to_table_from_path(os.path.join(user_staging_dir, otu_filename),
                                                 accept_empty_headers=False)

        # Processes the uploaded OTU file by removing unnecessary columns and extracting the headers and sample labels
        try:
            logger.info("Beginning to process the OTU table")
            base, headers, sample_labels, matrix_type = self.__process_base(self.user_id, pid, base_arr)
        except ValueError:
            logger.exception("OTU file contains non-integers")
            # Removes the project directory since the files in it are invalid
            shutil.rmtree(project_dir, ignore_errors=True)
            return OTU_DATATYPE_ERROR, ""
        except:
            logger.exception("Invalid OTU file format")
            # Removes the project directory since the files in it are invalid
            shutil.rmtree(project_dir, ignore_errors=True)
            return OTU_ERROR, ""

        # Creates map.txt file
        logger.info("Creating the map.txt file")
        map_file = Map(self.user_id, pid)
        map_file.project_name = project_name
        map_file.orig_otu_table_name = otu_filename
        map_file.orig_sample_metadata_name = sample_metadata_filename
        map_file.matrix_type = matrix_type
        map_file.num_samples = len(sample_labels)
        map_file.num_otus = len(headers)
        map_file.save()

        return OK, pid


    def __process_base(self, user_id, pid, base_arr, output_raw_otu_file_name=RAW_GENE_TABLE_FILENAME,
                       output_raw_otu_labels_file_name=RAW_GENE_TABLE_LABELS_FILENAME):
        '''
        Takes a TSV-separated base OTU file and extracts the header and the sample labels from the OTU file.
        Removes unnecessary columns if input is mothur derived file.
        Returns an OTU file (with only numeric values) and corresponding table header and sample labels
        :param user_id:
        :param pid:
        :return:
        '''

        project_dir = os.path.join(ProjectManager.DATA_DIRECTORY, user_id)
        project_dir = os.path.join(project_dir, pid)
        raw_table_path = os.path.join(project_dir, output_raw_otu_file_name)

        base = []

        matrix_type = "int"

        # Transpose the input base array to match the internal Mian engine format
        col_offset = 1
        row_offset = 1
        col = col_offset
        while col < len(base_arr[0]):
            new_row = []
            row = row_offset
            while row < len(base_arr):
                if base_arr[row][col] == "":
                    # Empty values will default to zero
                    new_row.append(0)
                else:
                    val = float(base_arr[row][col])
                    if val.is_integer():
                        new_row.append(int(float(base_arr[row][col])))
                    else:
                        new_row.append(float(base_arr[row][col]))
                        matrix_type = "float"
                row += 1
            base.append(new_row)
            col += 1

        headers = []
        for row in base_arr:
            headers.append(row[0])

        sample_labels = base_arr[0][col_offset:]

        labels = [headers, sample_labels]

        DataIO.table_to_tsv(base, user_id, pid, raw_table_path)
        DataIO.table_to_tsv(labels, user_id, pid, output_raw_otu_labels_file_name)

        return base, headers, sample_labels, matrix_type

    def create_project(self, pid, sampleFilter, sampleFilterRole, sampleFilterVals):

        map_file = Map(self.user_id, pid)

        project_dir = os.path.join(ProjectManager.DATA_DIRECTORY, self.user_id)
        project_dir = os.path.join(project_dir, pid)

        try:
            user_request = UserRequest(self.user_id, pid, "", "", "",
                                       "", [], sampleFilter, sampleFilterRole,
                                       sampleFilterVals, 0, "")

            table = OTUTable(self.user_id, pid, use_raw=True, use_np=False)
            orig_base = table.get_table()
            orig_headers = table.get_headers()
            orig_sample_labels = table.get_sample_labels()
            base, headers, sample_labels = table.filter_otu_table_by_metadata(orig_base, orig_headers, orig_sample_labels,
                                                                              user_request)
            initial_samples_removed = list(set(orig_sample_labels) - set(sample_labels))

            new_base = []
            i = 0
            while i < len(base):
                new_row = []
                j = 0
                while j < len(base[i]):
                    try:
                        if map_file.matrix_type == "float":
                            new_row.append(float(base[i][j]))
                        else:
                            new_row.append(int(base[i][j]))
                    except ValueError:
                        new_row.append(0)
                    j += 1
                new_base.append(new_row)
                i += 1

            logger.info("Subsampled file")


            # Updates map.txt file
            map_file.num_samples = len(sample_labels)
            map_file.num_otus = len(headers)
            map_file.save()

            return pid, ""
        except Exception as e:
            print(e)
            logger.exception("Error while processing the file format")
            # Removes the project directory since the files in it are invalid
            shutil.rmtree(project_dir, ignore_errors=True)
            return GENERAL_ERROR, ""

    def __validate_otu_table_sample_labels(self, sample_labels, sample_ids_from_metadata):
        missing = []
        for sample_label in sample_labels:
            if sample_label not in sample_ids_from_metadata:
                missing.append(sample_label)
        return missing
