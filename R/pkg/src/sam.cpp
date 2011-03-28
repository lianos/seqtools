#include "sam.h"

using namespace std;

RcppExport SEXP
collapse_sam_tag_list(SEXP taglist, SEXP the_sep) {
BEGIN_RCPP
    Rcpp::List tag_list(taglist);
    string sep = Rcpp::as<string>(the_sep);
    int ntags = tag_list.size();
    int i, j;

    vector<string> elt;
    vector<vector<string> > all_tags(ntags);
    for (i = 0; i < ntags; i++) {
        elt = Rcpp::as<vector<string> >(tag_list[i]);
        all_tags[i] = elt;
    }

    bool first;
    int nelts = all_tags[0].size();
    vector<string> ans(nelts);
    string current_tag, current_result;
    for (i = 0; i < nelts; i++) {
        first = true;
        current_result = "";
        for (j = 0; j < ntags; j++) {
            current_tag = all_tags[j][i];
            // cout << current_tag << " ";
            if (current_tag.length() > 0) {
                if (first) {
                    current_result = current_tag;
                    first = false;
                } else {
                    current_result += sep + current_tag;
                }
            }
        }
        ans[i] = current_result;
        // cout << "\n";
    }

    return Rcpp::wrap(ans);
END_RCPP
}
