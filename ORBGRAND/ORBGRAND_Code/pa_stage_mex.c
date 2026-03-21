/* ==========================================================
 * Precision-aware bitonic sorter MEX with per-stage MSB tie counters
 *
 * Input:
 *   y_soft       : real double vector of channel LLRs
 *
 * Outputs:
 *   sorted_idx   : sorted reliability indices (1-based for MATLAB)
 *   tie_count    : number of MSB ties per bitonic stage
 *   total_count  : number of valid comparisons per bitonic stage
 *
 * Stage definition:
 *   Stage s corresponds to the outer bitonic loop w = 2,4,8,...,n_effective
 *   so the number of stages is log2(n_effective).
 *
 * Tie probability in MATLAB:
 *   p_stage = double(tie_count) ./ double(total_count);
 * ========================================================== */

#include "mex.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* ----------------------------------------------------------
 * Round n up to the next power of two
 * ---------------------------------------------------------- */
static uint64_t next_pow2_u64(uint64_t x)
{
    uint64_t p = 1;
    while (p < x) p <<= 1;
    return p;
}

/* ----------------------------------------------------------
 * Integer log2 for powers of two
 * ---------------------------------------------------------- */
static uint64_t ilog2_u64(uint64_t x)
{
    uint64_t r = 0;
    while (x > 1) {
        x >>= 1;
        r++;
    }
    return r;
}

/* ----------------------------------------------------------
 * Bitonic sorter with per-stage MSB tie counting
 *
 * arr_msb     : 3-MSB quantized values used by the sorter
 * ind_order   : associated indices
 * n           : number of valid inputs
 * n_effective : padded power-of-two length
 *
 * tie_count[s]   = number of MSB ties at stage s
 * total_count[s] = number of valid comparisons at stage s
 *
 * Note:
 * Comparisons involving padded entries are excluded from
 * the counters, so the statistics reflect only valid bits.
 * ---------------------------------------------------------- */
static void bitonic_sort_count_ties(
    uint8_t  *arr_msb,
    uint32_t *ind_order,
    uint64_t  n,
    uint64_t  n_effective,
    uint64_t *tie_count,
    uint64_t *total_count
)
{
    uint64_t w, j, i, l;
    uint64_t stage_idx = 0;

    /* Outer loop defines the bitonic stage: w = 2,4,8,... */
    for (w = 2; w <= n_effective; w <<= 1, stage_idx++) {

        /* Internal compare-exchange passes within each stage */
        for (j = w >> 1; j > 0; j >>= 1) {
            for (i = 0; i < n_effective; i++) {
                l = i ^ j;   /* XOR gives compare partner */

                if (l > i) {
                    /* Count only comparisons between valid entries */
                    if ((i < n) && (l < n)) {
                        total_count[stage_idx]++;

                        /* MSB tie event */
                        if (arr_msb[i] == arr_msb[l]) {
                            tie_count[stage_idx]++;
                        }
                    }

                    /* Standard bitonic compare-and-swap */
                    if ((((i & w) == 0) && (arr_msb[i] > arr_msb[l])) ||   /* ascending */
                        (((i & w) != 0) && (arr_msb[i] < arr_msb[l]))) {   /* descending */

                        uint8_t  temp_val = arr_msb[i];
                        uint32_t temp_idx = ind_order[i];

                        arr_msb[i]   = arr_msb[l];
                        arr_msb[l]   = temp_val;

                        ind_order[i] = ind_order[l];
                        ind_order[l] = temp_idx;
                    }
                }
            }
        }
    }
}

/* ----------------------------------------------------------
 * Precision-aware sorter wrapper
 *
 * Quantization:
 *   6-bit sign-magnitude input
 *   1 sign bit + 5 magnitude bits
 *   sort using only the 3 MSBs of the magnitude
 * ---------------------------------------------------------- */
static void pa_sorter_with_counters(
    const double *y_soft,
    uint64_t      n,
    uint32_t     *sorted_idx_out,
    uint64_t     *tie_count,
    uint64_t     *total_count
)
{
    /* Quantization parameters */
    const double LLR_max = 31.0;   /* clipped range for 6-bit sign-magnitude */
    const int    B       = 6;
    const int    B_mag   = B - 1;  /* 5-bit magnitude */
    const int    MSB_NUM = 3;      /* use 3 MSBs for reduced-precision sorting */

    const uint64_t msb_shift = (uint64_t)(B_mag - MSB_NUM); /* 5-3 = 2 */
    const uint8_t  pad_val   = (uint8_t)((1U << MSB_NUM) - 1U); /* max 3-bit value = 7 */

    /* Pad input size to power of two for bitonic sorting */
    uint64_t n_effective = next_pow2_u64(n);

    /* Arrays used by the sorter */
    uint8_t  *LLR_mag_q    = (uint8_t *)calloc((size_t)n_effective, sizeof(uint8_t));
    uint32_t *sorted_list_q = (uint32_t *)calloc((size_t)n_effective, sizeof(uint32_t));

    if (LLR_mag_q == NULL || sorted_list_q == NULL) {
        free(LLR_mag_q);
        free(sorted_list_q);
        mexErrMsgTxt("Memory allocation failed.");
    }

    /* ------------------------------------------------------
     * Quantize valid LLRs and initialize indices
     * ------------------------------------------------------ */
    for (uint64_t i = 0; i < n; i++) {

        /* Clip soft input to supported LLR range */
        double L = fmax(fmin(y_soft[i], LLR_max), -LLR_max);

        /* Quantize 5-bit magnitude: 0...31 */
        uint8_t mag_q = (uint8_t)llround(fabs(L));

        /* Keep only the 3 MSBs for precision-aware sorting */
        LLR_mag_q[i] = (uint8_t)(mag_q >> msb_shift);

        /* Internal indices remain 0-based */
        sorted_list_q[i] = (uint32_t)i;
    }

    /* ------------------------------------------------------
     * Pad unused entries with maximum value so they migrate
     * to the end of the sorted sequence
     * ------------------------------------------------------ */
    for (uint64_t i = n; i < n_effective; i++) {
        LLR_mag_q[i] = pad_val;
        sorted_list_q[i] = (uint32_t)i;
    }

    /* Run bitonic sorter and collect tie statistics */
    bitonic_sort_count_ties(
        LLR_mag_q,
        sorted_list_q,
        n,
        n_effective,
        tie_count,
        total_count
    );

    /* Return only valid sorted indices, converted to MATLAB 1-based indexing */
    for (uint64_t i = 0; i < n; i++) {
        sorted_idx_out[i] = sorted_list_q[i] + 1U;
    }

    free(LLR_mag_q);
    free(sorted_list_q);
}

/* ----------------------------------------------------------
 * MEX interface
 *
 * MATLAB usage:
 *   [sorted_idx, tie_count, total_count] = pa_sorter_ties_mex(y_soft)
 * ---------------------------------------------------------- */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Check number of inputs/outputs */
    if (nrhs != 1) {
        mexErrMsgTxt("Usage: [sorted_idx, tie_count, total_count] = pa_sorter_ties_mex(y_soft)");
    }

    if (nlhs != 3) {
        mexErrMsgTxt("Three outputs required: sorted_idx, tie_count, total_count.");
    }

    /* Input must be a real double vector */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgTxt("Input y_soft must be a real double vector.");
    }

    if (mxGetN(prhs[0]) != 1 && mxGetM(prhs[0]) != 1) {
        mexErrMsgTxt("Input y_soft must be a vector.");
    }

    /* Read input */
    uint64_t n = (uint64_t)mxGetNumberOfElements(prhs[0]);
    const double *y_soft = mxGetPr(prhs[0]);

    /* Determine number of bitonic stages */
    uint64_t n_effective = next_pow2_u64(n);
    uint64_t n_stages    = ilog2_u64(n_effective);

    /* ------------------------------------------------------
     * Output 1: sorted indices (n x 1), uint32
     * ------------------------------------------------------ */
    plhs[0] = mxCreateNumericMatrix((mwSize)n, 1, mxUINT32_CLASS, mxREAL);
    uint32_t *sorted_idx = (uint32_t *)mxGetData(plhs[0]);

    /* ------------------------------------------------------
     * Output 2: tie counts per stage (n_stages x 1), uint64
     * ------------------------------------------------------ */
    plhs[1] = mxCreateNumericMatrix((mwSize)n_stages, 1, mxUINT64_CLASS, mxREAL);
    uint64_t *tie_count = (uint64_t *)mxGetData(plhs[1]);
    memset(tie_count, 0, (size_t)n_stages * sizeof(uint64_t));

    /* ------------------------------------------------------
     * Output 3: total comparisons per stage (n_stages x 1), uint64
     * ------------------------------------------------------ */
    plhs[2] = mxCreateNumericMatrix((mwSize)n_stages, 1, mxUINT64_CLASS, mxREAL);
    uint64_t *total_count = (uint64_t *)mxGetData(plhs[2]);
    memset(total_count, 0, (size_t)n_stages * sizeof(uint64_t));

    /* Run sorter + counters */
    pa_sorter_with_counters(
        y_soft,
        n,
        sorted_idx,
        tie_count,
        total_count
    );
}