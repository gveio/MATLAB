/* ==========================================================
 * Precision-aware bitonic sorter MEX with per-step MSB tie
 * and tie-swap bias counters
 *
 * File:
 *   pa_step_tie_bias_mex.c
 *
 * MATLAB usage:
 *   [sorted_idx, tie_count, total_count, tie_swap_count] = ...
 *       pa_step_tie_bias_mex(y_soft)
 *
 * Purpose:
 *   tie_count[t]:
 *      number of MSB3 ties at step t
 *
 *   total_count[t]:
 *      number of valid comparisons at step t
 *
 *   tie_swap_count[t]:
 *      among MSB3 ties, number of cases where the full 5-bit
 *      magnitude comparator would have swapped
 *
 * Metrics:
 *   p_tie(t) = tie_count(t) / total_count(t)
 *
 *   p_swap_given_tie(t) =
 *      tie_swap_count(t) / tie_count(t)
 *
 * ========================================================== */

#include "mex.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static uint64_t next_pow2_u64(uint64_t x)
{
    uint64_t p = 1;
    while (p < x) p <<= 1;
    return p;
}

static uint64_t ilog2_u64(uint64_t x)
{
    uint64_t r = 0;
    while (x > 1) {
        x >>= 1;
        r++;
    }
    return r;
}

static uint64_t num_bitonic_steps_u64(uint64_t n_effective)
{
    uint64_t m = ilog2_u64(n_effective);
    return (m * (m + 1)) / 2;
}

/* ----------------------------------------------------------
 * MSB3 sorter with full-magnitude tie-swap observation
 * ---------------------------------------------------------- */
static void bitonic_sort_count_tie_bias_per_step(
    uint8_t  *arr_msb,
    uint8_t  *arr_full,
    uint32_t *ind_order,
    uint64_t  n,
    uint64_t  n_effective,
    uint64_t *tie_count,
    uint64_t *total_count,
    uint64_t *tie_swap_count
)
{
    uint64_t w, j, i, l;
    uint64_t step_idx = 0;

    for (w = 2; w <= n_effective; w <<= 1) {
        for (j = w >> 1; j > 0; j >>= 1, step_idx++) {

            for (i = 0; i < n_effective; i++) {
                l = i ^ j;

                if (l > i) {

                    int ascending = ((i & w) == 0);

                    /* Count only valid real-input comparisons */
                    if ((i < n) && (l < n)) {
                        total_count[step_idx]++;

                        if (arr_msb[i] == arr_msb[l]) {
                            tie_count[step_idx]++;

                            /*
                             * Ask:
                             * If this MSB3 tie were resolved by the
                             * full 5-bit magnitude comparator, would
                             * the CAE perform a swap?
                             */
                            int full_swap =
                                (ascending  && (arr_full[i] > arr_full[l])) ||
                                (!ascending && (arr_full[i] < arr_full[l]));

                            if (full_swap) {
                                tie_swap_count[step_idx]++;
                            }
                        }
                    }

                    /*
                     * Actual sorter path:
                     * still sort using MSB3 only.
                     */
                    if ((ascending  && (arr_msb[i] > arr_msb[l])) ||
                        (!ascending && (arr_msb[i] < arr_msb[l]))) {

                        uint8_t  temp_msb  = arr_msb[i];
                        uint8_t  temp_full = arr_full[i];
                        uint32_t temp_idx  = ind_order[i];

                        arr_msb[i]   = arr_msb[l];
                        arr_full[i]  = arr_full[l];
                        ind_order[i] = ind_order[l];

                        arr_msb[l]   = temp_msb;
                        arr_full[l]  = temp_full;
                        ind_order[l] = temp_idx;
                    }
                }
            }
        }
    }
}

/* ----------------------------------------------------------
 * Wrapper
 * ---------------------------------------------------------- */
static void pa_sorter_with_tie_bias_counters(
    const double *y_soft,
    uint64_t      n,
    uint32_t     *sorted_idx_out,
    uint64_t     *tie_count,
    uint64_t     *total_count,
    uint64_t     *tie_swap_count
)
{
    const double LLR_max = 31.0;
    const int    B       = 6;
    const int    B_mag   = B - 1;
    const int    MSB_NUM = 3;

    const uint64_t msb_shift = (uint64_t)(B_mag - MSB_NUM);
    const uint8_t  pad_msb   = (uint8_t)((1U << MSB_NUM) - 1U);
    const uint8_t  pad_full  = 31U;

    uint64_t n_effective = next_pow2_u64(n);

    uint8_t  *LLR_mag_msb  = (uint8_t *)calloc((size_t)n_effective, sizeof(uint8_t));
    uint8_t  *LLR_mag_full = (uint8_t *)calloc((size_t)n_effective, sizeof(uint8_t));
    uint32_t *sorted_list  = (uint32_t *)calloc((size_t)n_effective, sizeof(uint32_t));

    if (LLR_mag_msb == NULL || LLR_mag_full == NULL || sorted_list == NULL) {
        free(LLR_mag_msb);
        free(LLR_mag_full);
        free(sorted_list);
        mexErrMsgTxt("Memory allocation failed.");
    }

    for (uint64_t i = 0; i < n; i++) {
        double L = fmax(fmin(y_soft[i], LLR_max), -LLR_max);
        uint8_t mag_q = (uint8_t)llround(fabs(L));

        if (mag_q > 31U) mag_q = 31U;

        LLR_mag_full[i] = mag_q;
        LLR_mag_msb[i]  = (uint8_t)(mag_q >> msb_shift);
        sorted_list[i]  = (uint32_t)i;
    }

    for (uint64_t i = n; i < n_effective; i++) {
        LLR_mag_full[i] = pad_full;
        LLR_mag_msb[i]  = pad_msb;
        sorted_list[i]  = (uint32_t)i;
    }

    bitonic_sort_count_tie_bias_per_step(
        LLR_mag_msb,
        LLR_mag_full,
        sorted_list,
        n,
        n_effective,
        tie_count,
        total_count,
        tie_swap_count
    );

    for (uint64_t i = 0; i < n; i++) {
        sorted_idx_out[i] = sorted_list[i] + 1U;
    }

    free(LLR_mag_msb);
    free(LLR_mag_full);
    free(sorted_list);
}

/* ----------------------------------------------------------
 * MEX interface
 * ---------------------------------------------------------- */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 1) {
        mexErrMsgTxt("Usage: [sorted_idx, tie_count, total_count, tie_swap_count] = pa_step_tie_bias_mex(y_soft)");
    }

    if (nlhs != 4) {
        mexErrMsgTxt("Four outputs required: sorted_idx, tie_count, total_count, tie_swap_count.");
    }

    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgTxt("Input y_soft must be a real double vector.");
    }

    if (mxGetN(prhs[0]) != 1 && mxGetM(prhs[0]) != 1) {
        mexErrMsgTxt("Input y_soft must be a vector.");
    }

    uint64_t n = (uint64_t)mxGetNumberOfElements(prhs[0]);
    const double *y_soft = mxGetPr(prhs[0]);

    uint64_t n_effective = next_pow2_u64(n);
    uint64_t n_steps = num_bitonic_steps_u64(n_effective);

    plhs[0] = mxCreateNumericMatrix((mwSize)n, 1, mxUINT32_CLASS, mxREAL);
    uint32_t *sorted_idx = (uint32_t *)mxGetData(plhs[0]);

    plhs[1] = mxCreateNumericMatrix((mwSize)n_steps, 1, mxUINT64_CLASS, mxREAL);
    uint64_t *tie_count = (uint64_t *)mxGetData(plhs[1]);
    memset(tie_count, 0, (size_t)n_steps * sizeof(uint64_t));

    plhs[2] = mxCreateNumericMatrix((mwSize)n_steps, 1, mxUINT64_CLASS, mxREAL);
    uint64_t *total_count = (uint64_t *)mxGetData(plhs[2]);
    memset(total_count, 0, (size_t)n_steps * sizeof(uint64_t));

    plhs[3] = mxCreateNumericMatrix((mwSize)n_steps, 1, mxUINT64_CLASS, mxREAL);
    uint64_t *tie_swap_count = (uint64_t *)mxGetData(plhs[3]);
    memset(tie_swap_count, 0, (size_t)n_steps * sizeof(uint64_t));

    pa_sorter_with_tie_bias_counters(
        y_soft,
        n,
        sorted_idx,
        tie_count,
        total_count,
        tie_swap_count
    );
}