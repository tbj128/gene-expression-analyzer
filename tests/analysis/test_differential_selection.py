import json

from gene.core.constants import SAMPLE_METADATA_FILENAME
from gene.analysis.differential_selection import DifferentialSelection
from tests.analysis.analysis_test_utils import AnalysisTestUtils
from gene.model.metadata import Metadata
import unittest


class TestDifferentialSelection(unittest.TestCase):

    def test_simple_differential_selection(self):

        user_request = AnalysisTestUtils.create_default_user_request()
        user_request.catvar = "Category"
        user_request.set_custom_attr("pvalthreshold", "0.01")
        user_request.set_custom_attr("pwVar1", "Control")
        user_request.set_custom_attr("pwVar2", "Disease")
        user_request.set_custom_attr("type", "ttest")

        otu_table = AnalysisTestUtils.get_test_input_as_table(AnalysisTestUtils.SIMPLE_TEST_CASE_ROOT)
        headers, sample_labels = AnalysisTestUtils.get_test_input_as_metadata(AnalysisTestUtils.SIMPLE_TEST_CASE_ROOT)
        metadata_table = AnalysisTestUtils.get_test_input_as_table(AnalysisTestUtils.SIMPLE_TEST_CASE_ROOT, SAMPLE_METADATA_FILENAME)
        metadata_col = AnalysisTestUtils.get_disease_metadata_values(AnalysisTestUtils.SIMPLE_TEST_CASE_ROOT)
        taxonomic_map = AnalysisTestUtils.get_test_taxonomy(AnalysisTestUtils.SIMPLE_TEST_CASE_ROOT)
        sample_ids_from_metadata = AnalysisTestUtils.get_sample_ids_from_metadata(
            AnalysisTestUtils.SIMPLE_TEST_CASE_ROOT)
        sample_id_to_metadata = {}
        i = 0
        while i < len(sample_ids_from_metadata):
            sample_id_to_metadata[sample_ids_from_metadata[i]] = metadata_col[i]
            i += 1

        metadata = Metadata("test", "test", False)
        metadata.set_table(metadata_table)

        plugin = DifferentialSelection()
        abundances = plugin.analyse(user_request, otu_table, headers, sample_labels, sample_id_to_metadata, taxonomic_map)
        print(json.dumps(abundances))
        expected_output = AnalysisTestUtils.get_expected_output(AnalysisTestUtils.SIMPLE_TEST_CASE_OUTPUT_ROOT,
                                                                "differential_selection_control_disease.json")
        comparison_output = AnalysisTestUtils.compare_two_objects(expected_output, abundances)
        if not comparison_output:
            print("Expected: ")
            print(expected_output)
            print("Actual: ")
            print(abundances)
        self.assertTrue(comparison_output)

    def test_simple_differential_selection_with_ancom(self):

        user_request = AnalysisTestUtils.create_default_user_request()
        user_request.catvar = "Category"
        user_request.set_custom_attr("pvalthreshold", "0.01")
        user_request.set_custom_attr("pwVar1", "Control")
        user_request.set_custom_attr("pwVar2", "Disease")

        otu_table = AnalysisTestUtils.get_test_input_as_table(AnalysisTestUtils.SIMPLE_TEST_CASE_ROOT)
        headers, sample_labels = AnalysisTestUtils.get_test_input_as_metadata(AnalysisTestUtils.SIMPLE_TEST_CASE_ROOT)
        metadata_table = AnalysisTestUtils.get_test_input_as_table(AnalysisTestUtils.SIMPLE_TEST_CASE_ROOT, SAMPLE_METADATA_FILENAME)
        metadata_col = AnalysisTestUtils.get_disease_metadata_values(AnalysisTestUtils.SIMPLE_TEST_CASE_ROOT)
        taxonomic_map = AnalysisTestUtils.get_test_taxonomy(AnalysisTestUtils.SIMPLE_TEST_CASE_ROOT)
        sample_ids_from_metadata = AnalysisTestUtils.get_sample_ids_from_metadata(
            AnalysisTestUtils.SIMPLE_TEST_CASE_ROOT)
        sample_id_to_metadata = {}
        i = 0
        while i < len(sample_ids_from_metadata):
            sample_id_to_metadata[sample_ids_from_metadata[i]] = metadata_col[i]
            i += 1

        metadata = Metadata("test", "test", False)
        metadata.set_table(metadata_table)

        plugin = DifferentialSelection()
        abundances = plugin.analyse_with_ancom(user_request, otu_table, headers, sample_labels, sample_id_to_metadata, taxonomic_map)
        print(json.dumps(abundances))
        expected_output = AnalysisTestUtils.get_expected_output(AnalysisTestUtils.SIMPLE_TEST_CASE_OUTPUT_ROOT,
                                                                "differential_selection_with_ancom_control_disease.json")
        comparison_output = AnalysisTestUtils.compare_two_objects(expected_output, abundances)
        if not comparison_output:
            print("Expected: ")
            print(expected_output)
            print("Actual: ")
            print(abundances)
        self.assertTrue(comparison_output)

if __name__ == '__main__':
    unittest.main()