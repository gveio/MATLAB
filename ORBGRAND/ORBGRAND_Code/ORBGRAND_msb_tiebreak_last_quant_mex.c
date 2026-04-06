/* C Functions */
#include <mex.h>
#include <math.h>
#include <stdint.h>
#include <stddef.h>
#define Inf 0x7fffffff

void copy_vector(double *v1,uint8_t *v2,uint64_t n){
    size_t i;
    for(i=0; i<n; i++){
        v1[i] = v2[i];
    }
}

int32_t findMax(int32_t a, int32_t b) {
    return !(b > a) ? a : b;
}

void hard_dec(uint8_t* y_hard,double *llr,uint64_t n){
    size_t i;

    for(i=0; i<n; i++){
        if(llr[i] >= 0){
           y_hard[i] = 0;
        }
        else{
            y_hard[i] = 1;
        }
    }
}

uint8_t memb_check(uint8_t *y, uint8_t *H,uint64_t n,uint64_t s){
    
    uint8_t syndrome;
    size_t i,j; 

    for (i = 0; i < s; i++) {
         syndrome = 0;
            for (j = 0; j < n; j++) {
                syndrome ^= (y[j] * H[i*n +j]);
            }
            if(syndrome == 1){
                return 1;
            }
    }

    return 0;
}

void mountain_build(int32_t *u, int32_t kk, int32_t w, int32_t W1, int32_t n1){
    size_t i;
    uint64_t W2,q,r;

    for(i = kk + 1; i < w; i++)
        u[i] = u[kk];
    W2 = W1;
    for( i = 0; i < w; i++)
        W2 -= u[i];
    q = floor( W2 / (n1 - u[kk]) );
    r = W2 - q*(n1 - u[kk]);
    if (q != 0){
        for(i = w-q; i < w; i++)
            u[i] = n1;
    }
    if (w > q)
        u[w-q-1] = u[w-q-1] + r;
}

static inline uint32_t msb_key(uint32_t x, const int MSB_NUM, const int B_MAG) {
    return x >> (B_MAG - MSB_NUM);
}

static inline uint32_t lsb_key(uint32_t x, const int LSB_NUM) {
    return x & ((1u << LSB_NUM) - 1u);
}

// returns 1 if (a > b) in lexicographic order (primary, secondary)
static inline int key_gt(uint32_t a, uint32_t b, int use_tiebreak, const int MSB_NUM, const int LSB_NUM, const int B_MAG)
{
    uint32_t ap = msb_key(a, MSB_NUM, B_MAG);
    uint32_t bp = msb_key(b, MSB_NUM, B_MAG);
    if (ap != bp) return (ap > bp);

    if (!use_tiebreak) return 0; // equal in MSB-space → do not force swap
    uint32_t as = lsb_key(a, LSB_NUM);
    uint32_t bs = lsb_key(b, LSB_NUM);
    return (as > bs);
}

// returns 1 if (a < b) in lexicographic order (primary, secondary)
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

/**
 * Mixed precision in stages:
 * - Always: MSB_NUM (e.g., 3) is used as primary sort key.
 * - Last 2 stages (w >= n/2): add LSB_NUM (e.g., 2) tie-break.
 *
 * Optional pruning:
 * - If prune_final_stage != 0, it SKIPS the very last stage (w == n).
 *   (This is risky for quality, but you can test it.)
 */
void bitonic_sort_msb3_tiebreak_last2(uint32_t *mag_q,
                                      uint32_t *ind_order,
                                      uint64_t n,
                                      const int B_MAG,
                                      const int MSB_NUM,
                                      const int LSB_NUM)
{

    for (uint64_t w = 2; w <= n; w <<= 1) {

        int use_tiebreak = (w == n);

        for (uint64_t j = (w >> 1); j > 0; j >>= 1) {

            for (uint64_t i = 0; i < n; i++) {

                uint64_t l = i ^ j;
                if (l <= i) continue;

                    int asc = ((i & w) == 0);

                    int do_swap = 0;
                    if (asc) {
                        // ascending: swap if key(i) > key(l)
                        do_swap = key_gt(mag_q[i], mag_q[l], use_tiebreak,
                                         MSB_NUM, LSB_NUM, B_MAG);
                    } else {
                        // descending: swap if key(i) < key(l)
                        do_swap = key_lt(mag_q[i], mag_q[l], use_tiebreak,
                                         MSB_NUM, LSB_NUM, B_MAG);
                    }

                    if (do_swap) {
                        uint32_t tmp_mag = mag_q[i];
                        mag_q[i] = mag_q[l];
                        mag_q[l] = tmp_mag;

                        uint32_t tmp_idx = ind_order[i];
                        ind_order[i] = ind_order[l];
                        ind_order[l] = tmp_idx;
                    }
            }
        }
    }
}

void ORBGRAND(double *y_decoded,double *n_guesses,double *y_soft,uint8_t *H,uint64_t n,uint64_t s,uint64_t n_guesses_max){

    size_t i,j,decoded=0,sy;
    /*vectors*/
    uint8_t* y_hard = (uint8_t* )calloc(n,sizeof(uint8_t));
    uint8_t* y_guessed = (uint8_t* )calloc(n,sizeof(uint8_t));
    uint8_t* tep = (uint8_t* )calloc(n,sizeof(uint8_t));

    /*variables for quantized values for sorter*/
    const double LLR_max = 31;
    const int B = 6;
    const int B_mag = B - 1;
    const int MSB_NUM = 3; // number of MSBs to extract
    const int LSB_NUM = 2; // number of LSBs to use as tie-break

    // round n up to next power of two
    uint64_t n_effective = 1;
    while (n_effective < n) n_effective *= 2;

    // allocate arrays of size n_effective
    uint32_t *LLR_mag_q = (uint32_t *)calloc(n_effective, sizeof(uint32_t));
    uint32_t *sorted_list_q = (uint32_t *)calloc(n_effective, sizeof(uint32_t));

    /*variables for pattern generator*/
    double temp=1+2*(double)n;
    int32_t W=0,w=0;
    int32_t k;
    int32_t* u = (int32_t*)calloc(n,sizeof(int32_t));
    int32_t* d = (int32_t*)calloc(n,sizeof(int32_t));
    int32_t* D = (int32_t*)calloc(n,sizeof(int32_t));

    /*hard y*/
    hard_dec(y_hard,y_soft,n);
    
    /*First query is hard y*/
    n_guesses[0]=1;  
    if(n_guesses_max == 0){
        n_guesses_max = Inf;
    }
    /*Codebook membership check on hard y*/
    sy = memb_check(y_hard,H,n,s);

    if (sy == 0){
        copy_vector(y_decoded, y_hard, n);
        decoded=1;
    }
    else{
        /*If codebook membership check on hard y fail sorter begins*/
        // Quantized
    for (uint64_t i = 0; i < n; i++) {
        // Clip and quantize
        double L = fmax(fmin(y_soft[i], LLR_max), -LLR_max);
        LLR_mag_q[i] = round(fabs(L) / LLR_max * (pow(2, B_mag) - 1));
        sorted_list_q[i] = i;
    }
    for (uint64_t i = n; i < n_effective; i++) {
        LLR_mag_q[i] = pow(2, B_mag) - 1; // pad = 31
        sorted_list_q[i] = i;
    }

   bitonic_sort_msb3_tiebreak_last2(LLR_mag_q, sorted_list_q, n_effective,B_mag, MSB_NUM, LSB_NUM);

    while(n_guesses[0]<n_guesses_max){
        W++; /*Increment logistic weight*/
        w = findMax(1, (int32_t)ceil((temp-sqrt(pow((double)temp, 2.0)-8*W))/2) ); /*Increment Hamming weight*/
        while ( w <= n ){  //w<=floor((sqrt(1+8*W)-1)/2))
            if(W < w*(w+1)/2)
                break;
                if(decoded) goto cleanup;
                /*landslide algorithm for generating integer partitions*/
                int32_t W1=W-w*(w+1)/2; 
                int32_t n1=n-w;
                /*Initial pattern*/
                    for (i = 0; i < w; i++)
                        u[i] = 0;
                    mountain_build(u, 0, w, W1, n1);
                    for(i=0; i<n; i++)
                        tep[i] = 0;  /*initialize tep*/
                    for(i=0; i<w; i++)
                         tep[sorted_list_q[u[i]+i]] = 1;                
                    for (i = 0; i < n; i++)
                        y_guessed[i] = y_hard[i] ^ tep[i];
                    n_guesses[0]++;
                    sy = memb_check(y_guessed,H,n,s);
                        if (sy == 0){
                            copy_vector(y_decoded,y_guessed,n);
                            decoded=1;
                            goto cleanup;
                        }
                    for (i = 0; i < w - 1; i++)
                        d[i] = u[i+1] - u[i];
                    d[w-1] = 0;
                    D[w-1] = d[w-1];
                    for (i = 1; i < w; i++)
                        D[w-i-1] = D[w-i] + d[w-i-1];
                    while( D[0] >= 2 ){
                        k = 0;
                        for (i = w-1; i > 0; i--){
                            if (D[i] >= 2){
                                k = i;
                                break;
                            }
                        }
                        u[k]++;
                        mountain_build(u,k,w,W1,n1);
                        for (i = 0; i < n; i++)
                             tep[i] = 0;
                        for(i=0; i<w; i++)
                            tep[sorted_list_q[u[i]+i]] = 1;
                        for (i = 0; i < n; i++)
                            y_guessed[i] = y_hard[i] ^ tep[i];
                        n_guesses[0]++;
                        sy = memb_check(y_guessed,H,n,s);
                        if (sy == 0){
                            copy_vector(y_decoded,y_guessed,n);
                            decoded=1;
                            break;
                        }        
                        for (size_t i = 0; i < w - 1; i++)
                            d[i] = u[i+1] - u[i];
                        d[w-1] = 0;
                        D[w-1] = d[w-1];
                        for (size_t i = 1; i < w; i++)
                            D[w-i-1] = D[w-i] + d[w-i-1];      
                    }  
                w++;
        }
        if(decoded) goto cleanup;
    }
    } 
    cleanup:
    free(y_hard);
    free(LLR_mag_q);
    free(sorted_list_q);
    free(u);  
    free(d);
    free(D);
    free(tep);
    free(y_guessed);
}

/* Mexfunction Interface */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    /* Check for proper number of arguments */
    if( nrhs != 3 )
       mexErrMsgTxt("Must have 3 input arguments.");
    if( nlhs != 2 )
        mexErrMsgTxt("Must have 3 output arguments.");

    if( mxGetN(prhs[0]) != 1 || !mxIsClass(prhs[0], "double"))
        mexErrMsgTxt("First Input (LLR) must be a column-vector of type double.");

    if( mxGetN(prhs[1]) != 1 || !mxIsClass(prhs[1], "uint8"))
        mexErrMsgTxt("Second Input (Parity check matrix) must be a column-vector of length s*n and type uint8.");

    if( mxGetNumberOfElements(prhs[2]) != 1 || !mxIsClass(prhs[2], "uint64"))
        mexErrMsgTxt("Third Input (Maximum number of guesses) must be an uint64 scalar.");

    uint64_t n = mxGetNumberOfElements(prhs[0]); // code length n
    uint64_t N = mxGetNumberOfElements(prhs[1]); // s*n
    uint64_t s = (uint64_t) N/n; 

    /* input */
    double *y_soft = mxGetPr(prhs[0]);
    uint8_t *H = (uint8_t *)mxGetData(prhs[1]);
    uint64_t n_guesses_max = mxGetScalar(prhs[2]);
    /* output */
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *y_decoded = (double *)mxGetData(plhs[0]);
    double *n_guesses = (double *)mxGetData(plhs[1]);
    /* use C functions in mexfunction */
    ORBGRAND(y_decoded,n_guesses,y_soft,H,n,s,n_guesses_max);
}
