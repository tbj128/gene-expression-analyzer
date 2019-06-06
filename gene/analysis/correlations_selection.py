# ===========================================
#
# mian Analysis Data Mining/ML Library
# @author: tbj128
#
# ===========================================

#
# Imports
#

#
# ======== R specific setup =========
#

from scipy import stats, math
from gene.analysis.analysis_base import AnalysisBase
from gene.core.statistics import Statistics
from gene.model.otu_table import OTUTable


class CorrelationsSelection(AnalysisBase):

    def run(self, user_request):
        table = OTUTable(user_request.user_id, user_request.pid)
        otu_table, headers, sample_labels = table.get_table_after_filtering_and_aggregation(user_request)

        metadata = table.get_sample_metadata()
        taxonomy_map = table.get_otu_metadata().get_taxonomy_map()

        return self.analyse(user_request, otu_table, headers, sample_labels, metadata, taxonomy_map)

    def analyse(self, user_request, base, headers, sample_labels, metadata, taxonomy_map):
        otu_to_genus = {}
        if int(user_request.level) == -1:
            # We want to display a short hint for the OTU using the genus (column 5)
            for header in headers:
                if header in taxonomy_map and len(taxonomy_map[header]) > 5:
                    otu_to_genus[header] = taxonomy_map[header][5]
                else:
                    otu_to_genus[header] = ""
        expvar = user_request.get_custom_attr("expvar")
        metadata_values = metadata.get_metadata_column_table_order(sample_labels, expvar)
        metadata_values = list(map(float, metadata_values))

        correlations = []
        pvals = []

        c = 0
        while c < len(base[0]):
            otu_name = headers[c]
            otu_vals = []
            r = 0
            while r < len(base):
                otu_vals.append(float(base[r][c]))
                r += 1

            coef, pval = stats.pearsonr(metadata_values, otu_vals)
            if math.isnan(coef):
                coef = 0
            if math.isnan(pval):
                pval = 1

            if int(user_request.level) == -1 and otu_name in otu_to_genus:
                correlations.append({"otu": otu_name, "coef": coef, "pval": pval, "hint": otu_to_genus[otu_name]})
            else:
                correlations.append({"otu": otu_name, "coef": coef, "pval": pval})
            pvals.append(pval)

            c += 1

        qvals = Statistics.getFDRCorrection(pvals)

        i = 0
        while i < len(qvals):
            correlations[i]["qval"] = qvals[i]
            i += 1

        return {"correlations": correlations}
