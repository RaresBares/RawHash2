/* PELT CUDA kernel — parallel SSR cost computation across candidates per j step.
 * Built with: nvcc -O3 --shared -Xcompiler -fPIC pelt_cuda.cu -o libpelt_cuda.so
 * The .so is dlopen'd at runtime by detect_events_pelt_cuda; if absent the C
 * implementation falls back to the AVX2 path.
 *
 * API (extern "C"):
 *   int pelt_cuda_run(const float *signal, int n, float pen_mult, int min_size,
 *                     int win_size, int *out_cps, int max_cp);
 * Returns number of change-points written into out_cps, or -1 on error.
 */
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define CUDA_CHECK(x) do { cudaError_t e=(x); if(e!=cudaSuccess){ \
    fprintf(stderr,"[pelt_cuda] %s\n",cudaGetErrorString(e)); return -1; } } while(0)

/* For each (t in cands) compute SSR cost of segment [t, j) and contribute
 * to F[j] via parallel reduction. One thread per candidate. */
__global__ void pelt_step_kernel(
    const double * __restrict__ ps, const double * __restrict__ pss,
    const double * __restrict__ F, const int * __restrict__ cands, int n_cands,
    int j, int min_size, double penalty,
    double *out_best_total, int *out_best_t)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    __shared__ double s_best[256];
    __shared__ int    s_best_t[256];
    s_best[threadIdx.x]   = 1e30;
    s_best_t[threadIdx.x] = -1;

    if (tid < n_cands) {
        int t = cands[tid];
        if (j - t >= min_size) {
            double seg_sum   = ps[j]  - ps[t];
            double seg_sumsq = pss[j] - pss[t];
            double seg_len   = (double)(j - t);
            double seg_mean  = seg_sum / seg_len;
            double cost      = seg_sumsq - seg_sum * seg_mean;
            double total     = F[t] + cost + penalty;
            s_best[threadIdx.x]   = total;
            s_best_t[threadIdx.x] = t;
        }
    }
    __syncthreads();

    /* Block-level reduction (min). */
    for (int off = blockDim.x / 2; off > 0; off >>= 1) {
        if (threadIdx.x < off) {
            if (s_best[threadIdx.x + off] < s_best[threadIdx.x]) {
                s_best[threadIdx.x]   = s_best[threadIdx.x + off];
                s_best_t[threadIdx.x] = s_best_t[threadIdx.x + off];
            }
        }
        __syncthreads();
    }

    if (threadIdx.x == 0) {
        /* Atomic-min is hard for double; serialize on per-block result. */
        unsigned long long *addr = (unsigned long long*)out_best_total;
        unsigned long long old   = *addr, assumed;
        double cand_val = s_best[0];
        do {
            assumed = old;
            if (__longlong_as_double(assumed) <= cand_val) break;
            old = atomicCAS(addr, assumed,
                            __double_as_longlong(cand_val));
        } while (assumed != old);
        if (__longlong_as_double(*addr) == cand_val)
            *out_best_t = s_best_t[0];
    }
}

extern "C" int pelt_cuda_run(const float *signal, int n,
                             float pen_mult, int min_size, int win_size,
                             int *out_cps, int max_cp)
{
    if (n < 10) return 0;
    double penalty = (double)pen_mult * log((double)n);
    int n1 = n + 1;

    /* Host prefix sums first (cheap, n is small). */
    double *h_ps  = (double*)malloc(n1 * sizeof(double));
    double *h_pss = (double*)malloc(n1 * sizeof(double));
    h_ps[0] = h_pss[0] = 0;
    for (int i = 0; i < n; i++) {
        h_ps[i+1]  = h_ps[i]  + signal[i];
        h_pss[i+1] = h_pss[i] + (double)signal[i] * signal[i];
    }

    double *d_ps, *d_pss, *d_F, *d_best;
    int *d_cands, *d_prev, *d_best_t;
    CUDA_CHECK(cudaMalloc(&d_ps,  n1 * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_pss, n1 * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_F,   n1 * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_cands, n1 * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_prev,  n1 * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_best,    sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_best_t,  sizeof(int)));

    CUDA_CHECK(cudaMemcpy(d_ps,  h_ps,  n1 * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_pss, h_pss, n1 * sizeof(double), cudaMemcpyHostToDevice));
    free(h_ps); free(h_pss);

    /* Initialize F and prev on device. */
    double *h_F = (double*)malloc(n1 * sizeof(double));
    h_F[0] = -penalty;
    for (int i = 1; i <= n; i++) h_F[i] = 1e30;
    CUDA_CHECK(cudaMemcpy(d_F, h_F, n1 * sizeof(double), cudaMemcpyHostToDevice));
    free(h_F);
    CUDA_CHECK(cudaMemset(d_prev, 0, n1 * sizeof(int)));

    /* Candidate list begins with [0]. */
    int *h_cands = (int*)malloc(n1 * sizeof(int));
    h_cands[0] = 0;
    int n_cands = 1;
    CUDA_CHECK(cudaMemcpy(d_cands, h_cands, sizeof(int), cudaMemcpyHostToDevice));

    /* Main DP loop — j sequential, inner reduction parallel. */
    for (int j = min_size; j <= n; j++) {
        double init_best = 1e30;
        int    init_t    = 0;
        CUDA_CHECK(cudaMemcpy(d_best,   &init_best, sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_best_t, &init_t,    sizeof(int),    cudaMemcpyHostToDevice));

        int threads = 256;
        int blocks  = (n_cands + threads - 1) / threads;
        if (blocks < 1) blocks = 1;
        pelt_step_kernel<<<blocks, threads>>>(
            d_ps, d_pss, d_F, d_cands, n_cands,
            j, min_size, penalty, d_best, d_best_t);

        double best_total; int best_t;
        CUDA_CHECK(cudaMemcpy(&best_total, d_best,   sizeof(double), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(&best_t,     d_best_t, sizeof(int),    cudaMemcpyDeviceToHost));
        if (best_total < 1e30) {
            CUDA_CHECK(cudaMemcpy((char*)d_F + j*sizeof(double),
                                  &best_total, sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy((char*)d_prev + j*sizeof(int),
                                  &best_t, sizeof(int), cudaMemcpyHostToDevice));
        }

        /* Window-PELT pruning: keep candidates within last `win_size` samples + new j */
        int new_count = 0;
        int win_lo = (j > win_size) ? j - win_size : 0;
        for (int i = 0; i < n_cands; i++) if (h_cands[i] >= win_lo) h_cands[new_count++] = h_cands[i];
        h_cands[new_count++] = j;
        n_cands = new_count;
        CUDA_CHECK(cudaMemcpy(d_cands, h_cands, n_cands * sizeof(int), cudaMemcpyHostToDevice));
    }

    /* Backtrack on host. */
    int *h_prev = (int*)malloc(n1 * sizeof(int));
    CUDA_CHECK(cudaMemcpy(h_prev, d_prev, n1 * sizeof(int), cudaMemcpyDeviceToHost));
    int *cps = (int*)malloc(n1 * sizeof(int));
    int n_cp = 0;
    int pos = n;
    while (pos > 0 && n_cp < max_cp) {
        if (h_prev[pos] > 0) cps[n_cp++] = h_prev[pos];
        pos = h_prev[pos];
    }
    for (int i = 0; i < n_cp / 2; i++) {
        int tmp = cps[i]; cps[i] = cps[n_cp-1-i]; cps[n_cp-1-i] = tmp;
    }
    for (int i = 0; i < n_cp && i < max_cp; i++) out_cps[i] = cps[i];
    free(cps); free(h_prev); free(h_cands);

    cudaFree(d_ps); cudaFree(d_pss); cudaFree(d_F);
    cudaFree(d_cands); cudaFree(d_prev);
    cudaFree(d_best); cudaFree(d_best_t);
    return n_cp;
}
