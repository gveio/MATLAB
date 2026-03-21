/* ==========================================================
 * Precision-aware bitonic sorter MEX with per-step MSB tie counters
 *
 * Input:
 *   y_soft       : real double vector of channel LLRs
 *
 * Outputs:
 *   sorted_idx   : sorted reliability indices (1-based for MATLAB)
 *   tie_count    : number of MSB ties per compare-exchange step
 *   total_count  : number of valid comparisons per compare-exchange step
 *
 * Step definition:
 *   One step corresponds to one internal bitonic pass, i.e. one value of j
 *   inside the nested loops:
 *       for (w = 2; w <= n_effective; w <<= 1)
 *           for (j = w>>1; j > 0; j >>= 1)
 *
 * Number of steps:
 *   m = log2(n_effective)
 *   n_steps = m(m+1)/2
 *
 * Tie probability in MATLAB:
 *   p_step = double(tie_count) ./ double(total_count);
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
 * Total number of bitonic compare-exchange steps
 * m = log2(n_effective), steps = m(m+1)/2
 * ---------------------------------------------------------- */
static uint64_t num_bitonic_steps_u64(uint64_t n_effective)
{
    uint64_t m = ilog2_u64(n_effective);
    return (m * (m + 1)) / 2;
}

/* ----------------------------------------------------------
 * Bitonic sorter with per-step MSB tie counting
 *
 * arr_msb      : 3-MSB quantized values used by the sorter
 * ind_order    : associated indices
 * n            : number of valid inputs
 * n_effective  : padded power-of-two length
 *
 * tie_count[t]   = number of MSB ties at CAE step t
 * total_count[t] = number of valid comparisons at CAE step t
 *
 * Note:
 * Comparisons involving padded entries are excluded from
 * the counters, so the statistics reflect only valid bits.
 * ---------------------------------------------------------- */
static void bitonic_sort_count_ties_per_step(
    uint8_t  *arr_msb,
    uint32_t *ind_order,
    uint64_t  n,
    uint64_t  n_effective,
    uint64_t *tie_count,
    uint64_t *total_count
)
{
    uint64_t w, j, i, l;
    uint64_t step_idx = 0;

    for (w = 2; w <= n_effective; w <<= 1) {
        for (j = w >> 1; j > 0; j >>= 1, step_idx++) {

            for (i = 0; i < n_effective; i++) {
                l = i ^ j;

                if (l > i) {
                    /* Count only comparisons between valid entries */
                    if ((i < n) && (l < n)) {
                        total_count[step_idx]++;

                        if (arr_msb[i] == arr_msb[l]) {
                            tie_count[step_idx]++;
                        }
                    }

                    /* Standard bitonic compare-and-swap */
                    if ((((i & w) == 0) && (arr_msb[i] > arr_msb[l])) ||
                        (((i & w) != 0) && (arr_msb[i] < arr_msb[l]))) {

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
    const double LLR_max = 31.0;   /* clipped range for 6-bit sign-magnitude */
    const int    B       = 6;
    const int    B_mag   = B - 1;  /* 5-bit magnitude */
    const int    MSB_NUM = 3;      /* use 3 MSBs for reduced-precision sorting */

    const uint64_t msb_shift = (uint64_t)(B_mag - MSB_NUM); /* 5-3 = 2 */
    const uint8_t  pad_val   = (uint8_t)((1U << MSB_NUM) - 1U); /* 7 */

    uint64_t n_effective = next_pow2_u64(n);

    uint8_t  *LLR_mag_q     = (uint8_t *)calloc((size_t)n_effective, sizeof(uint8_t));
    uint32_t *sorted_list_q = (uint32_t *)calloc((size_t)n_effective, sizeof(uint32_t));

    if (LLR_mag_q == NULL || sorted_list_q == NULL) {
        free(LLR_mag_q);
        free(sorted_list_q);
        mexErrMsgTxt("Memory allocation failed.");
    }

    /* Quantize valid LLRs and initialize indices */
    for (uint64_t i = 0; i < n; i++) {
        double L = fmax(fmin(y_soft[i], LLR_max), -LLR_max);
        uint8_t mag_q = (uint8_t)llround(fabs(L));   /* magnitude 0...31 */
        LLR_mag_q[i] = (uint8_t)(mag_q >> msb_shift); /* keep top 3 bits */
        sorted_list_q[i] = (uint32_t)i;
    }

    /* Pad unused entries with max value */
    for (uint64_t i = n; i < n_effective; i++) {
        LLR_mag_q[i] = pad_val;
        sorted_list_q[i] = (uint32_t)i;
    }

    bitonic_sort_count_ties_per_step(
        LLR_mag_q,
        sorted_list_q,
        n,
        n_effective,
        tie_count,
        total_count
    );

    /* Return only valid sorted indices, MATLAB 1-based */
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
    if (nrhs != 1) {
        mexErrMsgTxt("Usage: [sorted_idx, tie_count, total_count] = pa_sorter_ties_mex(y_soft)");
    }

    if (nlhs != 3) {
        mexErrMsgTxt("Three outputs required: sorted_idx, tie_count, total_count.");
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

    /* Output 1: sorted indices */
    plhs[0] = mxCreateNumericMatrix((mwSize)n, 1, mxUINT32_CLASS, mxREAL);
    uint32_t *sorted_idx = (uint32_t *)mxGetData(plhs[0]);

    /* Output 2: tie counts per step */
    plhs[1] = mxCreateNumericMatrix((mwSize)n_steps, 1, mxUINT64_CLASS, mxREAL);
    uint64_t *tie_count = (uint64_t *)mxGetData(plhs[1]);
    memset(tie_count, 0, (size_t)n_steps * sizeof(uint64_t));

    /* Output 3: total comparisons per step */
    plhs[2] = mxCreateNumericMatrix((mwSize)n_steps, 1, mxUINT64_CLASS, mxREAL);
    uint64_t *total_count = (uint64_t *)mxGetData(plhs[2]);
    memset(total_count, 0, (size_t)n_steps * sizeof(uint64_t));

    pa_sorter_with_counters(
        y_soft,
        n,
        sorted_idx,
        tie_count,
        total_count
    );
}