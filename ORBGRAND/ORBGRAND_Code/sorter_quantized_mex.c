/* ==========================================================
 * Bitonic sorter MEX for true-valued and quantized LLR sorting
 *
 * Input:
 *   y_soft : real double vector of channel LLRs
 *   mode   : scalar
 *            0 -> true-valued reference: sort abs(y_soft)
 *            1 -> baseline quantized (full 5-bit magnitude)
 *            2 -> MSB2
 *            3 -> MSB3
 *            4 -> MSB4
 *            5 -> PA1 = MSB3 + LSB1 tie-break in final 2 stages
 *            6 -> PA2 = MSB3 + LSB2 tie-break in final 2 stages
 *
 * Output:
 *   sorted_idx : sorted reliability indices (1-based for MATLAB)
 * ========================================================== */

#include "mex.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------- */
static uint64_t next_pow2_u64(uint64_t x)
{
    uint64_t p = 1;
    while (p < x) p <<= 1;
    return p;
}

/* ----------------------------------------------------------
 * True-valued bitonic sort
 * ---------------------------------------------------------- */
static void bitonic_sort_double(
    double   *arr,
    uint32_t *ind_order,
    uint64_t  n_effective
)
{
    for (uint64_t w = 2; w <= n_effective; w <<= 1) {
        for (uint64_t j = w >> 1; j > 0; j >>= 1) {
            for (uint64_t i = 0; i < n_effective; i++) {
                uint64_t l = i ^ j;
                if (l <= i) continue;

                if ((((i & w) == 0) && (arr[i] > arr[l])) ||
                    (((i & w) != 0) && (arr[i] < arr[l]))) {

                    double   tmp_val = arr[i];
                    uint32_t tmp_idx = ind_order[i];

                    arr[i] = arr[l];
                    arr[l] = tmp_val;

                    ind_order[i] = ind_order[l];
                    ind_order[l] = tmp_idx;
                }
            }
        }
    }
}

/* ----------------------------------------------------------
 * Build reduced-precision key for plain MSB modes
 * ---------------------------------------------------------- */
static uint8_t build_sort_key(uint8_t mag_q, uint32_t mode)
{
    switch (mode) {
        case 1: return mag_q;                  /* baseline full 5-bit */
        case 2: return (uint8_t)(mag_q >> 3);  /* MSB2 */
        case 3: return (uint8_t)(mag_q >> 2);  /* MSB3 */
        case 4: return (uint8_t)(mag_q >> 1);  /* MSB4 */
        default:
            mexErrMsgTxt("Invalid mode for build_sort_key.");
            return 0;
    }
}

/* ----------------------------------------------------------
 * Standard bitonic sort for quantized baseline/MSB2/MSB3/MSB4
 * ---------------------------------------------------------- */
static void bitonic_sort_indices_u8(
    uint8_t  *arr_key,
    uint32_t *ind_order,
    uint64_t  n_effective
)
{
    for (uint64_t w = 2; w <= n_effective; w <<= 1) {
        for (uint64_t j = w >> 1; j > 0; j >>= 1) {
            for (uint64_t i = 0; i < n_effective; i++) {
                uint64_t l = i ^ j;
                if (l <= i) continue;

                if ((((i & w) == 0) && (arr_key[i] > arr_key[l])) ||
                    (((i & w) != 0) && (arr_key[i] < arr_key[l]))) {

                    uint8_t  tmp_val = arr_key[i];
                    uint32_t tmp_idx = ind_order[i];

                    arr_key[i] = arr_key[l];
                    arr_key[l] = tmp_val;

                    ind_order[i] = ind_order[l];
                    ind_order[l] = tmp_idx;
                }
            }
        }
    }
}

/* ----------------------------------------------------------
 * PA helper functions
 * ---------------------------------------------------------- */
static inline uint32_t msb_key(uint32_t x, const int MSB_NUM, const int B_MAG)
{
    return x >> (B_MAG - MSB_NUM);
}

static inline uint32_t lsb_key(uint32_t x, const int LSB_NUM)
{
    return x & ((1u << LSB_NUM) - 1u);
}

static inline int key_gt(uint32_t a, uint32_t b, int use_tiebreak,
                         const int MSB_NUM, const int LSB_NUM, const int B_MAG)
{
    uint32_t ap = msb_key(a, MSB_NUM, B_MAG);
    uint32_t bp = msb_key(b, MSB_NUM, B_MAG);
    if (ap != bp) return (ap > bp);

    if (!use_tiebreak) return 0;
    uint32_t as = lsb_key(a, LSB_NUM);
    uint32_t bs = lsb_key(b, LSB_NUM);
    return (as > bs);
}

static inline int key_lt(uint32_t a, uint32_t b, int use_tiebreak,
                         const int MSB_NUM, const int LSB_NUM, const int B_MAG)
{
    uint32_t ap = msb_key(a, MSB_NUM, B_MAG);
    uint32_t bp = msb_key(b, MSB_NUM, B_MAG);
    if (ap != bp) return (ap < bp);

    if (!use_tiebreak) return 0;
    uint32_t as = lsb_key(a, LSB_NUM);
    uint32_t bs = lsb_key(b, LSB_NUM);
    return (as < bs);
}

/* ----------------------------------------------------------
 * Bitonic sort for PA1 / PA2
 * ---------------------------------------------------------- */
static void bitonic_sort_msb3_tiebreak_last2(
    uint32_t *mag_q,
    uint32_t *ind_order,
    uint64_t  n_effective,
    const int B_MAG,
    const int MSB_NUM,
    const int LSB_NUM
)
{
    for (uint64_t w = 2; w <= n_effective; w <<= 1) {
        int use_tiebreak = (w >= (n_effective >> 1)); /* last 2 stages */

        for (uint64_t j = w >> 1; j > 0; j >>= 1) {
            for (uint64_t i = 0; i < n_effective; i++) {
                uint64_t l = i ^ j;
                if (l <= i) continue;

                int asc = ((i & w) == 0);
                int do_swap = 0;

                if (asc) {
                    do_swap = key_gt(mag_q[i], mag_q[l], use_tiebreak,
                                     MSB_NUM, LSB_NUM, B_MAG);
                } else {
                    do_swap = key_lt(mag_q[i], mag_q[l], use_tiebreak,
                                     MSB_NUM, LSB_NUM, B_MAG);
                }

                if (do_swap) {
                    uint32_t tmp_mag = mag_q[i];
                    uint32_t tmp_idx = ind_order[i];

                    mag_q[i] = mag_q[l];
                    mag_q[l] = tmp_mag;

                    ind_order[i] = ind_order[l];
                    ind_order[l] = tmp_idx;
                }
            }
        }
    }
}

/* ----------------------------------------------------------
 * Main sorter core
 * ---------------------------------------------------------- */
static void sorter_mex_core(
    const double *y_soft,
    uint64_t      n,
    uint32_t      mode,
    uint32_t     *sorted_idx_out
)
{
    const double LLR_max = 31.0;
    const int    B_MAG   = 5;

    uint64_t n_effective = next_pow2_u64(n);

    uint32_t *ind_order = (uint32_t *)calloc((size_t)n_effective, sizeof(uint32_t));
    if (ind_order == NULL) {
        mexErrMsgTxt("Memory allocation failed.");
    }

    /* ---------------- true-valued reference ---------------- */
    if (mode == 0) {
        double *reliability = (double *)calloc((size_t)n_effective, sizeof(double));
        if (reliability == NULL) {
            free(ind_order);
            mexErrMsgTxt("Memory allocation failed.");
        }

        for (uint64_t i = 0; i < n; i++) {
            reliability[i] = fabs(y_soft[i]);
            ind_order[i]   = (uint32_t)i;
        }

        for (uint64_t i = n; i < n_effective; i++) {
            reliability[i] = 31.0; /* same padding style as your ORBGRAND code */
            ind_order[i]   = (uint32_t)i;
        }

        bitonic_sort_double(reliability, ind_order, n_effective);
        free(reliability);
    }

    /* ---------------- quantized family ---------------- */
    else {
        uint32_t *mag_q_full = (uint32_t *)calloc((size_t)n_effective, sizeof(uint32_t));
        if (mag_q_full == NULL) {
            free(ind_order);
            mexErrMsgTxt("Memory allocation failed.");
        }

        for (uint64_t i = 0; i < n; i++) {
            double L = fmax(fmin(y_soft[i], LLR_max), -LLR_max);
            uint32_t mag_q = (uint32_t)llround(fabs(L)); /* 0...31 */
            mag_q_full[i] = mag_q;
            ind_order[i]  = (uint32_t)i;
        }

        for (uint64_t i = n; i < n_effective; i++) {
            mag_q_full[i] = 31U;
            ind_order[i]  = (uint32_t)i;
        }

        if (mode >= 1 && mode <= 4) {
            uint8_t *sort_key = (uint8_t *)calloc((size_t)n_effective, sizeof(uint8_t));
            if (sort_key == NULL) {
                free(mag_q_full);
                free(ind_order);
                mexErrMsgTxt("Memory allocation failed.");
            }

            for (uint64_t i = 0; i < n_effective; i++) {
                sort_key[i] = build_sort_key((uint8_t)mag_q_full[i], mode);
            }

            bitonic_sort_indices_u8(sort_key, ind_order, n_effective);
            free(sort_key);
        }
        else if (mode == 5) {
            bitonic_sort_msb3_tiebreak_last2(mag_q_full, ind_order, n_effective,
                                             B_MAG, 3, 1); /* PA1 */
        }
        else if (mode == 6) {
            bitonic_sort_msb3_tiebreak_last2(mag_q_full, ind_order, n_effective,
                                             B_MAG, 3, 2); /* PA2 */
        }
        else {
            free(mag_q_full);
            free(ind_order);
            mexErrMsgTxt("Invalid mode. Use 0..6.");
        }

        free(mag_q_full);
    }

    for (uint64_t i = 0; i < n; i++) {
        sorted_idx_out[i] = ind_order[i] + 1U; /* MATLAB 1-based */
    }

    free(ind_order);
}

/* ----------------------------------------------------------
 * MEX interface
 *
 * MATLAB usage:
 *   sorted_idx = sorter_quantized_mex(y_soft, mode)
 * ---------------------------------------------------------- */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2) {
        mexErrMsgTxt("Usage: sorted_idx = sorter_quantized_mex(y_soft, mode)");
    }

    if (nlhs != 1) {
        mexErrMsgTxt("One output required: sorted_idx.");
    }

    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgTxt("Input y_soft must be a real double vector.");
    }

    if (mxGetN(prhs[0]) != 1 && mxGetM(prhs[0]) != 1) {
        mexErrMsgTxt("Input y_soft must be a vector.");
    }

    if (mxGetNumberOfElements(prhs[1]) != 1) {
        mexErrMsgTxt("Input mode must be a scalar.");
    }

    uint64_t n = (uint64_t)mxGetNumberOfElements(prhs[0]);
    const double *y_soft = mxGetPr(prhs[0]);

    double mode_in = mxGetScalar(prhs[1]);
    if (mode_in < 0 || mode_in > 6 || floor(mode_in) != mode_in) {
        mexErrMsgTxt("mode must be one of {0,1,2,3,4,5,6}.");
    }
    uint32_t mode = (uint32_t)mode_in;

    plhs[0] = mxCreateNumericMatrix((mwSize)n, 1, mxUINT32_CLASS, mxREAL);
    uint32_t *sorted_idx = (uint32_t *)mxGetData(plhs[0]);

    sorter_mex_core(y_soft, n, mode, sorted_idx);
}