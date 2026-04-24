/* ==========================================================
 * Precision-aware bitonic sorter MEX with per-CAE MSB tie
 * and tie-swap bias counters
 *
 * File:
 *   pa_cae_tie_bias_mex.c
 *
 * MATLAB:
 *   [sorted_idx, tie_count, total_count, tie_swap_count] = ...
 *       pa_cae_tie_bias_mex(y_soft)
 *
 * Outputs:
 *   sorted_idx       : n x 1 uint32, MATLAB 1-based sorted indices
 *   tie_count        : n_steps x n_cae uint64
 *   total_count      : n_steps x n_cae uint64
 *   tie_swap_count   : n_steps x n_cae uint64
 *
 * Meaning:
 *   tie_count(step,cae):
 *      number of MSB3 ties at that CAE location
 *
 *   total_count(step,cae):
 *      number of valid comparisons at that CAE location
 *
 *   tie_swap_count(step,cae):
 *      among MSB3 ties, number of cases where the full 5-bit
 *      comparator would have swapped
 *
 * Metrics in MATLAB:
 *   p_tie = double(tie_count) ./ double(total_count);
 *
 *   p_swap_given_tie = double(tie_swap_count) ./ double(tie_count);
 *
 * ========================================================== */

#include "mex.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static uint64_t next_pow2_cae_bias_u64(uint64_t x)
{
    uint64_t p = 1;
    while (p < x) p <<= 1;
    return p;
}

static uint64_t ilog2_cae_bias_u64(uint64_t x)
{
    uint64_t r = 0;
    while (x > 1) {
        x >>= 1;
        r++;
    }
    return r;
}

static uint64_t num_bitonic_steps_cae_bias_u64(uint64_t n_effective)
{
    uint64_t m = ilog2_cae_bias_u64(n_effective);
    return (m * (m + 1)) / 2;
}

/* ----------------------------------------------------------
 * MSB3 sorter with per-CAE full-magnitude tie-swap observation
 * ---------------------------------------------------------- */
static void bitonic_sort_count_cae_tie_bias(
    uint8_t  *arr_msb,
    uint8_t  *arr_full,
    uint32_t *ind_order,
    uint64_t  n,
    uint64_t  n_effective,
    uint64_t  n_steps,
    uint64_t *tie_count,
    uint64_t *total_count,
    uint64_t *tie_swap_count
)
{
    uint64_t w, j, i, l;
    uint64_t step_idx = 0;

    for (w = 2; w <= n_effective; w <<= 1) {
        for (j = w >> 1; j > 0; j >>= 1, step_idx++) {

            uint64_t cae_idx = 0;

            for (i = 0; i < n_effective; i++) {
                l = i ^ j;

                if (l > i) {

                    /*
                     * MATLAB uses column-major layout.
                     * Matrix size is n_steps x n_cae.
                     * Element (step_idx, cae_idx) maps to:
                     *   idx = step_idx + cae_idx * n_steps
                     */
                    uint64_t idx = step_idx + cae_idx * n_steps;

                    int ascending = ((i & w) == 0);

                    if ((i < n) && (l < n)) {
                        total_count[idx]++;

                        if (arr_msb[i] == arr_msb[l]) {
                            tie_count[idx]++;

                            /*
                             * Diagnostic question:
                             * If MSB3 ties here, would the full 5-bit
                             * comparator have swapped?
                             */
                            if ((ascending  && (arr_full[i] > arr_full[l])) ||
                                (!ascending && (arr_full[i] < arr_full[l]))) {
                                tie_swap_count[idx]++;
                            }
                        }
                    }

                    /*
                     * Actual simulated sorter path:
                     * MSB3-only.
                     * If MSB3 ties, no swap is performed.
                     */
                    if ((ascending  && (arr_msb[i] > arr_msb[l])) ||
                        (!ascending && (arr_msb[i] < arr_msb[l]))) {

                        uint8_t  tmp_msb  = arr_msb[i];
                        uint8_t  tmp_full = arr_full[i];
                        uint32_t tmp_idx  = ind_order[i];

                        arr_msb[i]   = arr_msb[l];
                        arr_full[i]  = arr_full[l];
                        ind_order[i] = ind_order[l];

                        arr_msb[l]   = tmp_msb;
                        arr_full[l]  = tmp_full;
                        ind_order[l] = tmp_idx;
                    }

                    cae_idx++;
                }
            }
        }
    }
}

/* ----------------------------------------------------------
 * Wrapper
 * ---------------------------------------------------------- */
static void pa_sorter_with_cae_tie_bias_counters(
    const double *y_soft,
    uint64_t      n,
    uint32_t     *sorted_idx_out,
    uint64_t     *tie_count,
    uint64_t     *total_count,
    uint64_t     *tie_swap_count
)
{
    const double LLR_max = 31.0;

    const int B_mag   = 5;
    const int MSB_NUM = 3;

    const uint64_t msb_shift = (uint64_t)(B_mag - MSB_NUM);
    const uint8_t  pad_msb   = (uint8_t)((1U << MSB_NUM) - 1U);
    const uint8_t  pad_full  = 31U;

    uint64_t n_effective = next_pow2_cae_bias_u64(n);
    uint64_t n_steps     = num_bitonic_steps_cae_bias_u64(n_effective);

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

        if (mag_q > 31U) {
            mag_q = 31U;
        }

        LLR_mag_full[i] = mag_q;
        LLR_mag_msb[i]  = (uint8_t)(mag_q >> msb_shift);
        sorted_list[i]  = (uint32_t)i;
    }

    for (uint64_t i = n; i < n_effective; i++) {
        LLR_mag_full[i] = pad_full;
        LLR_mag_msb[i]  = pad_msb;
        sorted_list[i]  = (uint32_t)i;
    }

    bitonic_sort_count_cae_tie_bias(
        LLR_mag_msb,
        LLR_mag_full,
        sorted_list,
        n,
        n_effective,
        n_steps,
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
        mexErrMsgTxt("Usage: [sorted_idx, tie_count, total_count, tie_swap_count] = pa_cae_tie_bias_mex(y_soft)");
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

    uint64_t n_effective = next_pow2_cae_bias_u64(n);
    uint64_t n_steps     = num_bitonic_steps_cae_bias_u64(n_effective);
    uint64_t n_cae       = n_effective / 2;

    plhs[0] = mxCreateNumericMatrix((mwSize)n, 1, mxUINT32_CLASS, mxREAL);
    uint32_t *sorted_idx = (uint32_t *)mxGetData(plhs[0]);

    plhs[1] = mxCreateNumericMatrix((mwSize)n_steps, (mwSize)n_cae, mxUINT64_CLASS, mxREAL);
    uint64_t *tie_count = (uint64_t *)mxGetData(plhs[1]);
    memset(tie_count, 0, (size_t)n_steps * (size_t)n_cae * sizeof(uint64_t));

    plhs[2] = mxCreateNumericMatrix((mwSize)n_steps, (mwSize)n_cae, mxUINT64_CLASS, mxREAL);
    uint64_t *total_count = (uint64_t *)mxGetData(plhs[2]);
    memset(total_count, 0, (size_t)n_steps * (size_t)n_cae * sizeof(uint64_t));

    plhs[3] = mxCreateNumericMatrix((mwSize)n_steps, (mwSize)n_cae, mxUINT64_CLASS, mxREAL);
    uint64_t *tie_swap_count = (uint64_t *)mxGetData(plhs[3]);
    memset(tie_swap_count, 0, (size_t)n_steps * (size_t)n_cae * sizeof(uint64_t));

    pa_sorter_with_cae_tie_bias_counters(
        y_soft,
        n,
        sorted_idx,
        tie_count,
        total_count,
        tie_swap_count
    );
}