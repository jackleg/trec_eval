/* 
   Copyright (c) 2008 - Chris Buckley. 

   Permission is granted for use and modification of this file for
   research, non-commercial purposes. 
*/

#include "common.h"
#include "sysfunc.h"
#include "trec_eval.h"
#include "functions.h"
#include "trec_format.h"

static int 
te_calc_mrr_cut (const EPI *epi, const REL_INFO *rel_info,
		    const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval);
static long long_cutoff_array[] = {1, 3, 5, 10};
static PARAMS default_mrr_cutoffs = {
    NULL, sizeof (long_cutoff_array) / sizeof (long_cutoff_array[0]),
    &long_cutoff_array[0]};

/* See trec_eval.h for definition of TREC_MEAS */
TREC_MEAS te_meas_mrr_cut =
    {"mrr_cut",
     "    Reciprocal Rank of the first relevant retrieved doc.\n\
    Measure is most useful for tasks in which there is only one relevant\n\
    doc, or the user only wants one relevant doc.\n",
     te_init_meas_a_float_cut_long,
     te_calc_mrr_cut,
     te_acc_meas_a_cut,
     te_calc_avg_meas_a_cut,
     te_print_single_meas_a_cut,
     te_print_final_meas_a_cut,
     (void *) &default_mrr_cutoffs, -1};

static int 
te_calc_mrr_cut (const EPI *epi, const REL_INFO *rel_info,
		    const RESULTS *results, const TREC_MEAS *tm,
		    TREC_EVAL *eval)
{
    RES_RELS res_rels;
    long i;
    long  *cutoffs = (long *) tm->meas_params->param_values;
    long cutoff_index = 0;
	char rel = 0;
	double res = 0.0;

    if (UNDEF == te_form_res_rels (epi, rel_info, results, &res_rels))
	return (UNDEF);

    for (i = 0; i < res_rels.num_ret; i++) {
		if (res_rels.results_rel_list[i] >= epi->relevance_level) {
			// 정답일 경우
			rel = i+1;
			break;
		}
	}

    for (i = 0; i < res_rels.num_ret; i++) {
		if (i+1 == cutoffs[cutoff_index]) {
			if (rel <= i+1 && rel) {
				eval->values[tm->eval_index + cutoff_index].value = (double) 1.0 / (double) rel;
			    //fprintf(stderr, "relevance_level=%ld %ld %d %d\n", epi->relevance_level, rel, cutoffs[cutoff_index], cutoff_index);
				res = (double) 1.0 / (double) rel;
			} else {
				eval->values[tm->eval_index + cutoff_index].value = (double) 0.0;
				res = (double)0.0;
			}
			if (++cutoff_index == tm->meas_params->num_params)
				break;
		}
    }

    while (cutoff_index < tm->meas_params->num_params) {
			eval->values[tm->eval_index + cutoff_index].value = res;
			cutoff_index++;
    }
    return (1);
}
