# ===========================================
# 
# mian Analysis Alpha/Beta Diversity Library
# @author: tbj128
#
# ===========================================

#
# Imports
#

#
# ======== R specific setup =========
#
import logging

import rpy2.robjects as robjects
import rpy2.rlike.container as rlc
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

from scipy import stats, math
from skbio import TreeNode
from skbio.diversity import alpha_diversity
from io import StringIO

from gene.analysis.analysis_base import AnalysisBase
from gene.core.statistics import Statistics

from gene.model.otu_table import OTUTable
from gene.model.map import Map

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class AlphaDiversity(AnalysisBase):
    r = robjects.r

    #
    # ======== Main code begins =========
    #

    rcode = """
    
    alphaDiversity <- function(allOTUs, alphaType, alphaContext) {
        alphaDiv = diversity(allOTUs, index = alphaType)
        if (alphaContext == "evenness") {
            S <- specnumber(allOTUs)
            if (alphaType == "shannon") {
                J <- alphaDiv/log(S)
                return(J)
            } else if (alphaType == "simpson") {
                # simpson index = 1 - D
                J <- (1-alphaDiv)/S
                return(J)
            } else {
                # invsimpson index = 1/D
                J <- (1/alphaDiv)/S
                return(J)
            }
        } else if (alphaContext == "speciesnumber") {
            S <- specnumber(allOTUs)
            return(S)
        } else {
            return(alphaDiv)
        }
    }
    
    """

    veganR = SignatureTranslatedAnonymousPackage(rcode, "veganR")

    def run(self, user_request):
        table = OTUTable(user_request.user_id, user_request.pid)

        # No OTUs should be excluded for diversity analysis
        otu_table, headers, sample_labels = table.get_table_after_filtering_and_aggregation(user_request)

        metadata_values = table.get_sample_metadata().get_metadata_column_table_order(sample_labels, user_request.get_custom_attr("expvar"))
        if user_request.get_custom_attr("colorvar") != "None":
            color_metadata_values = table.get_sample_metadata().get_metadata_column_table_order(sample_labels, user_request.get_custom_attr("colorvar"))
        else:
            color_metadata_values = []

        if user_request.get_custom_attr("sizevar") != "None":
            size_metadata_values = table.get_sample_metadata().get_metadata_column_table_order(sample_labels, user_request.get_custom_attr("sizevar"))
        else:
            size_metadata_values = []

        sample_ids_to_metadata_map = table.get_sample_metadata().get_sample_id_to_metadata_map(user_request.get_custom_attr("expvar"))
        phylogenetic_tree = table.get_phylogenetic_tree()

        return self.analyse(user_request, otu_table, headers, sample_labels, metadata_values, color_metadata_values, size_metadata_values, sample_ids_to_metadata_map, phylogenetic_tree)

    def analyse(self, user_request, otu_table, headers, sample_labels, metadata_values, color_metadata_values, size_metadata_values, sample_ids_to_metadata_map, phylogenetic_tree):
        logger.info("Starting Alpha Diversity analysis")

        plotType = user_request.get_custom_attr("plotType")
        alphaType = user_request.get_custom_attr("alphaType")
        alphaContext = user_request.get_custom_attr("alphaContext")
        statisticalTest = user_request.get_custom_attr("statisticalTest")

        vals = []
        if alphaType == "faith_pd":
            if phylogenetic_tree == "":
                return {
                    "no_tree": True
                }

            project_map = Map(user_request.user_id, user_request.pid)
            if project_map.matrix_type == "float":
                return {
                    "has_float": True
                }

            otu_table = otu_table.astype(int)

            tree = TreeNode.read(StringIO(phylogenetic_tree))
            if len(tree.root().children) > 2:
                # Ensure that the tree is rooted if it is not already rooted
                tree = tree.root_at_midpoint()
            vals = alpha_diversity(alphaType, otu_table, ids=sample_labels, otu_ids=headers, tree=tree)

        else:
            # Creates an R-compatible dictionary of columns to vectors of column values WITHOUT headers
            allOTUs = []
            col = 0
            while col < len(otu_table[0]):
                colVals = []
                row = 0
                while row < len(otu_table):
                    sampleID = sample_labels[row]
                    if len(metadata_values) == 0 or sampleID in sample_ids_to_metadata_map:
                        colVals.append(otu_table[row][col])
                    row += 1
                allOTUs.append((headers[col], robjects.FloatVector(colVals)))
                col += 1

            logger.info("After creating an R-compatible dictionary")

            od = rlc.OrdDict(allOTUs)
            dataf = robjects.DataFrame(od)

            logger.info("Before vegan alpha diversity")

            vals = self.veganR.alphaDiversity(dataf, alphaType, alphaContext)

            logger.info("After vegan alpha diversity")


        if plotType == "boxplot":
            # Calculate the statistical p-value
            abundances = {}
            statsAbundances = {}
            i = 0
            while i < len(vals):
                obj = {}
                obj["s"] = str(sample_labels[i])
                if vals[i] == float('inf'):
                    raise ValueError("Cannot have infinite values")
                obj["a"] = round(vals[i], 6)
                meta = metadata_values[i] if len(metadata_values) > 0 else "All"

                # Group the abundance values under the corresponding metadata values
                if meta in statsAbundances:
                    statsAbundances[meta].append(vals[i])
                    abundances[meta].append(obj)
                else:
                    statsAbundances[meta] = [vals[i]]
                    abundances[meta] = [obj]

                i += 1
            statistics = Statistics.getTtest(statsAbundances, statisticalTest)

            logger.info("After T-test")

            abundancesObj = {}
            abundancesObj["abundances"] = abundances
            abundancesObj["stats"] = statistics
            return abundancesObj
        else:
            corrArr = []
            corrValArr1 = []
            corrValArr2 = []
            i = 0
            while i < len(vals):
                obj = {}
                obj["s"] = str(sample_labels[i])
                obj["c1"] = float(metadata_values[i])
                obj["c2"] = float(vals[i])
                obj["color"] = color_metadata_values[i] if len(color_metadata_values) == len(vals) else ""
                obj["size"] = float(size_metadata_values[i]) if len(size_metadata_values) == len(vals) else ""
                corrArr.append(obj)
                corrValArr1.append(float(metadata_values[i]))
                corrValArr2.append(float(vals[i]))
                i += 1

            coef, pval = stats.pearsonr(corrValArr1, corrValArr2)
            if math.isnan(coef):
                coef = 0
            if math.isnan(pval):
                pval = 1

            abundances_obj = {"corrArr": corrArr, "coef": coef, "pval": pval}
            return abundances_obj
