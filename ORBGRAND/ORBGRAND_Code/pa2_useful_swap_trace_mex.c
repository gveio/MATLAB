/* ==========================================================
 * PA2 useful-swap trace MEX
 *
 * MATLAB:
 * [idx_msb3, idx_pa2, tie_count, total_count, tie_swap_count, ...
 *  swap_step, swap_cae, swap_a, swap_b] = pa2_useful_swap_trace_mex(y_soft)
 *
 * swap_* traces only PA2 swaps caused by MSB3 ties in the LAST TWO stages.
 * Indices returned to MATLAB are 1-based.
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

/* ==========================================================
 * MSB3-only sorter
 * ========================================================== */
static void bitonic_sort_msb3(
    uint8_t  *arr_msb,
    uint32_t *ind_order,
    uint64_t  n_effective
)
{
    uint64_t w, j, i, l;

    for (w = 2; w <= n_effective; w <<= 1) {
        for (j = w >> 1; j > 0; j >>= 1) {
            for (i = 0; i < n_effective; i++) {

                l = i ^ j;

                if (l > i) {

                    int ascending = ((i & w) == 0);

                    int do_swap =
                        (ascending  && (arr_msb[i] > arr_msb[l])) ||
                        (!ascending && (arr_msb[i] < arr_msb[l]));

                    if (do_swap) {
                        uint8_t  tmp_msb = arr_msb[i];
                        uint32_t tmp_idx = ind_order[i];

                        arr_msb[i]   = arr_msb[l];
                        ind_order[i] = ind_order[l];

                        arr_msb[l]   = tmp_msb;
                        ind_order[l] = tmp_idx;
                    }
                }
            }
        }
    }
}

/* ==========================================================
 * PA2 sorter + MSB3-tie PA2-swap trace
 * ========================================================== */
static uint64_t bitonic_sort_pa2_trace(
    uint8_t  *arr_msb,
    uint8_t  *arr_full,
    uint32_t *ind_order,
    uint64_t  n,
    uint64_t  n_effective,
    uint64_t  n_steps,
    uint64_t *tie_count,
    uint64_t *total_count,
    uint64_t *tie_swap_count,
    uint32_t *swap_step,
    uint32_t *swap_cae,
    uint32_t *swap_a,
    uint32_t *swap_b
)
{
    uint64_t w, j, i, l;
    uint64_t step_idx = 0;
    uint64_t trace_len = 0;

    uint64_t m_bits = ilog2_u64(n_effective);
    uint64_t last2_start_step = ((m_bits - 2) * (m_bits - 1)) / 2;

    for (w = 2; w <= n_effective; w <<= 1) {
        for (j = w >> 1; j > 0; j >>= 1, step_idx++) {

            uint64_t cae_idx = 0;
            int in_last_two_stages = (step_idx >= last2_start_step);

            for (i = 0; i < n_effective; i++) {

                l = i ^ j;

                if (l > i) {

                    uint64_t mat_idx = step_idx + cae_idx * n_steps;
                    int ascending = ((i & w) == 0);

                    int valid_pair = ((i < n) && (l < n));

                    if (valid_pair) {
                        total_count[mat_idx]++;
                    }

                    int msb_tie = (arr_msb[i] == arr_msb[l]);

                    if (valid_pair && msb_tie) {
                        tie_count[mat_idx]++;
                    }

                    /*
                     * PA2 comparator:
                     *   if MSB differs -> compare MSB3
                     *   if MSB ties   -> compare full 5-bit magnitude
                     */
                    int do_swap = 0;

                    if (arr_msb[i] > arr_msb[l]) {
                        do_swap = ascending ? 1 : 0;
                    }
                    else if (arr_msb[i] < arr_msb[l]) {
                        do_swap = ascending ? 0 : 1;
                    }
                    else {
                        if (arr_full[i] > arr_full[l]) {
                            do_swap = ascending ? 1 : 0;
                        }
                        else if (arr_full[i] < arr_full[l]) {
                            do_swap = ascending ? 0 : 1;
                        }
                        else {
                            do_swap = 0;
                        }
                    }

                    /*
                     * Count/trace only PA2 swaps that happen because of MSB3 tie.
                     */
                    if (valid_pair && msb_tie && do_swap) {
                        tie_swap_count[mat_idx]++;

                        if (in_last_two_stages) {
                            swap_step[trace_len] = (uint32_t)(step_idx + 1);     /* MATLAB 1-based */
                            swap_cae[trace_len]  = (uint32_t)(cae_idx + 1);      /* MATLAB 1-based */
                            swap_a[trace_len]    = ind_order[i] + 1U;            /* original index, 1-based */
                            swap_b[trace_len]    = ind_order[l] + 1U;            /* original index, 1-based */
                            trace_len++;
                        }
                    }

                    if (do_swap) {
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

    return trace_len;
}

/* ==========================================================
 * Main wrapper
 * ========================================================== */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 1) {
        mexErrMsgTxt("Usage: [idx_msb3, idx_pa2, tie_count, total_count, tie_swap_count, swap_step, swap_cae, swap_a, swap_b] = pa2_useful_swap_trace_mex(y_soft)");
    }

    if (nlhs != 9) {
        mexErrMsgTxt("Nine outputs required.");
    }

    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgTxt("Input y_soft must be a real double vector.");
    }

    if (mxGetN(prhs[0]) != 1 && mxGetM(prhs[0]) != 1) {
        mexErrMsgTxt("Input y_soft must be a vector.");
    }

    const double *y_soft = mxGetPr(prhs[0]);
    uint64_t n = (uint64_t)mxGetNumberOfElements(prhs[0]);

    const double LLR_max = 31.0;
    const int B_mag = 5;
    const int MSB_NUM = 3;
    const uint64_t msb_shift = (uint64_t)(B_mag - MSB_NUM);

    uint64_t n_effective = next_pow2_u64(n);
    uint64_t n_steps = num_bitonic_steps_u64(n_effective);
    uint64_t n_cae = n_effective / 2;
    uint64_t max_trace = n_steps * n_cae;

    uint8_t *msb_msb3 = (uint8_t *)calloc((size_t)n_effective, sizeof(uint8_t));
    uint8_t *msb_pa2  = (uint8_t *)calloc((size_t)n_effective, sizeof(uint8_t));
    uint8_t *full_pa2 = (uint8_t *)calloc((size_t)n_effective, sizeof(uint8_t));

    uint32_t *idx_msb3 = (uint32_t *)calloc((size_t)n_effective, sizeof(uint32_t));
    uint32_t *idx_pa2  = (uint32_t *)calloc((size_t)n_effective, sizeof(uint32_t));

    if (!msb_msb3 || !msb_pa2 || !full_pa2 || !idx_msb3 || !idx_pa2) {
        mexErrMsgTxt("Memory allocation failed.");
    }

    for (uint64_t i = 0; i < n; i++) {
        double L = fmax(fmin(y_soft[i], LLR_max), -LLR_max);
        uint8_t mag_q = (uint8_t)llround(fabs(L));

        if (mag_q > 31U) {
            mag_q = 31U;
        }

        uint8_t msb = (uint8_t)(mag_q >> msb_shift);

        msb_msb3[i] = msb;
        msb_pa2[i]  = msb;
        full_pa2[i] = mag_q;

        idx_msb3[i] = (uint32_t)i;
        idx_pa2[i]  = (uint32_t)i;
    }

    for (uint64_t i = n; i < n_effective; i++) {
        msb_msb3[i] = 7U;
        msb_pa2[i]  = 7U;
        full_pa2[i] = 31U;

        idx_msb3[i] = (uint32_t)i;
        idx_pa2[i]  = (uint32_t)i;
    }

    /* Outputs */
    plhs[0] = mxCreateNumericMatrix((mwSize)n, 1, mxUINT32_CLASS, mxREAL);
    plhs[1] = mxCreateNumericMatrix((mwSize)n, 1, mxUINT32_CLASS, mxREAL);

    uint32_t *out_msb3 = (uint32_t *)mxGetData(plhs[0]);
    uint32_t *out_pa2  = (uint32_t *)mxGetData(plhs[1]);

    plhs[2] = mxCreateNumericMatrix((mwSize)n_steps, (mwSize)n_cae, mxUINT64_CLASS, mxREAL);
    plhs[3] = mxCreateNumericMatrix((mwSize)n_steps, (mwSize)n_cae, mxUINT64_CLASS, mxREAL);
    plhs[4] = mxCreateNumericMatrix((mwSize)n_steps, (mwSize)n_cae, mxUINT64_CLASS, mxREAL);

    uint64_t *tie_count      = (uint64_t *)mxGetData(plhs[2]);
    uint64_t *total_count    = (uint64_t *)mxGetData(plhs[3]);
    uint64_t *tie_swap_count = (uint64_t *)mxGetData(plhs[4]);

    memset(tie_count, 0, (size_t)n_steps * (size_t)n_cae * sizeof(uint64_t));
    memset(total_count, 0, (size_t)n_steps * (size_t)n_cae * sizeof(uint64_t));
    memset(tie_swap_count, 0, (size_t)n_steps * (size_t)n_cae * sizeof(uint64_t));

    uint32_t *tmp_swap_step = (uint32_t *)calloc((size_t)max_trace, sizeof(uint32_t));
    uint32_t *tmp_swap_cae  = (uint32_t *)calloc((size_t)max_trace, sizeof(uint32_t));
    uint32_t *tmp_swap_a    = (uint32_t *)calloc((size_t)max_trace, sizeof(uint32_t));
    uint32_t *tmp_swap_b    = (uint32_t *)calloc((size_t)max_trace, sizeof(uint32_t));

    if (!tmp_swap_step || !tmp_swap_cae || !tmp_swap_a || !tmp_swap_b) {
        mexErrMsgTxt("Trace allocation failed.");
    }

    /* Sort MSB3 */
    bitonic_sort_msb3(msb_msb3, idx_msb3, n_effective);

    /* Sort PA2 + trace PA2 swaps caused by MSB3 ties */
    uint64_t trace_len = bitonic_sort_pa2_trace(
        msb_pa2,
        full_pa2,
        idx_pa2,
        n,
        n_effective,
        n_steps,
        tie_count,
        total_count,
        tie_swap_count,
        tmp_swap_step,
        tmp_swap_cae,
        tmp_swap_a,
        tmp_swap_b
    );

    for (uint64_t i = 0; i < n; i++) {
        out_msb3[i] = idx_msb3[i] + 1U;
        out_pa2[i]  = idx_pa2[i] + 1U;
    }

    plhs[5] = mxCreateNumericMatrix((mwSize)trace_len, 1, mxUINT32_CLASS, mxREAL);
    plhs[6] = mxCreateNumericMatrix((mwSize)trace_len, 1, mxUINT32_CLASS, mxREAL);
    plhs[7] = mxCreateNumericMatrix((mwSize)trace_len, 1, mxUINT32_CLASS, mxREAL);
    plhs[8] = mxCreateNumericMatrix((mwSize)trace_len, 1, mxUINT32_CLASS, mxREAL);

    memcpy(mxGetData(plhs[5]), tmp_swap_step, (size_t)trace_len * sizeof(uint32_t));
    memcpy(mxGetData(plhs[6]), tmp_swap_cae,  (size_t)trace_len * sizeof(uint32_t));
    memcpy(mxGetData(plhs[7]), tmp_swap_a,    (size_t)trace_len * sizeof(uint32_t));
    memcpy(mxGetData(plhs[8]), tmp_swap_b,    (size_t)trace_len * sizeof(uint32_t));

    free(msb_msb3);
    free(msb_pa2);
    free(full_pa2);
    free(idx_msb3);
    free(idx_pa2);

    free(tmp_swap_step);
    free(tmp_swap_cae);
    free(tmp_swap_a);
    free(tmp_swap_b);
}