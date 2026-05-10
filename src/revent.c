#include "revent.h"
#include "kalloc.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>
#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>
#if defined(__AVX2__)
#include <immintrin.h>
#endif
#include "rutils.h"

/* Optional GPU PELT (libpelt_cuda.so loaded at runtime) */
typedef int (*pelt_cuda_fn_t)(const float *signal, int n, float pen_mult,
                              int min_size, int win_size, int *out_cps, int max_cp);
static pelt_cuda_fn_t g_pelt_cuda_fn = NULL;
static int g_pelt_cuda_tried = 0;
static void try_load_pelt_cuda(void) {
    if (g_pelt_cuda_tried) return;
    g_pelt_cuda_tried = 1;
    const char *paths[] = {
        getenv("RH2_PELT_CUDA_LIB"),
        "./libpelt_cuda.so",
        "/home/rsahleanu/rawhash2/lib/libpelt_cuda.so",
        "libpelt_cuda.so",
        NULL,
    };
    for (int i = 0; paths[i]; i++) {
        if (!paths[i] || !*paths[i]) continue;
        void *h = dlopen(paths[i], RTLD_NOW);
        if (h) {
            g_pelt_cuda_fn = (pelt_cuda_fn_t)dlsym(h, "pelt_cuda_run");
            if (g_pelt_cuda_fn) {
                fprintf(stderr, "[pelt_cuda] loaded %s\n", paths[i]);
                return;
            }
            dlclose(h);
        }
    }
    if (getenv("RH2_DEBUG_EVENTS"))
        fprintf(stderr, "[pelt_cuda] no CUDA lib found, using AVX2/scalar fallback\n");
}

//Some of the functions here are adopted from the Sigmap implementation (https://github.com/haowenz/sigmap/tree/c9a40483264c9514587a36555b5af48d3f054f6f). We have optimized the Sigmap implementation to work with the hash tables efficiently.

typedef struct ri_detect_s {
	int DEF_PEAK_POS;
	float DEF_PEAK_VAL;
	float *sig;
	uint32_t s_len;
	float threshold;
	uint32_t window_length;
	uint32_t masked_to;
	int peak_pos;
	float peak_value;
	int valid_peak;
} ri_detect_t;

static inline void comp_prefix_prefixsq(const float *sig,
										const uint32_t s_len,
										float* prefix_sum,
										float* prefix_sum_square)
{
	assert(s_len > 0);

	prefix_sum[0] = 0.0f;
	prefix_sum_square[0] = 0.0f;
	for (uint32_t i = 0; i < s_len; ++i) {
		prefix_sum[i+1] = prefix_sum[i] + sig[i];
		prefix_sum_square[i+1] = prefix_sum_square[i] + sig[i]*sig[i];
	}
}

static inline float* comp_tstat(void *km,
							 	const float *prefix_sum,
								const float *prefix_sum_square,
								const uint32_t s_len,
								const uint32_t w_len)
{
  	const float eta = FLT_MIN;
	float* tstat = (float*)ri_kcalloc(km, s_len+1, sizeof(float));
	if (s_len < 2*w_len || w_len < 2) return tstat;
	memset(tstat, 0, w_len*sizeof(float));
	
	for (uint32_t i = w_len; i <= s_len - w_len; ++i) {
		float sum1 = prefix_sum[i];
		float sumsq1 = prefix_sum_square[i];
		if (i > w_len) {
			sum1 -= prefix_sum[i - w_len];
			sumsq1 -= prefix_sum_square[i - w_len];
		}
		float sum2 = prefix_sum[i + w_len] - prefix_sum[i];
		float sumsq2 = prefix_sum_square[i + w_len] - prefix_sum_square[i];
		float mean1 = sum1 / w_len;
		float mean2 = sum2 / w_len;
		float combined_var = (sumsq1/w_len - mean1*mean1 + sumsq2/w_len - mean2*mean2)/w_len;
		// Prevent problem due to very small variances
		combined_var = fmaxf(combined_var, eta);
		// t-stat
		//  Formula is a simplified version of Student's t-statistic for the
		//  special case where there are two samples of equal size with
		//  differing variance
		const float delta_mean = mean2 - mean1;
		tstat[i] = fabs(delta_mean) / sqrt(combined_var);
	}
	// fudge boundaries
	memset(tstat+s_len-w_len+1, 0, (w_len)*sizeof(float));

	return tstat;
}

// static inline float calculate_adaptive_peak_height(const float *prefix_sum, const float *prefix_sum_square, uint32_t current_index, uint32_t window_length, float base_peak_height) {
//     // Ensure we don't go beyond signal bounds
//     uint32_t start_index = current_index > window_length ? current_index - window_length : 0;
//     uint32_t end_index = current_index + window_length;

//     float sum = prefix_sum[end_index] - prefix_sum[start_index];
//     float sumsq = prefix_sum_square[end_index] - prefix_sum_square[start_index];
//     float mean = sum / (end_index - start_index);
//     float variance = (sumsq / (end_index - start_index)) - (mean * mean);
//     float stddev = sqrtf(variance);

//     // Example adaptive strategy: Increase peak height in high-variance regions
//     return base_peak_height * (1 + stddev);
// }

static inline uint32_t gen_peaks(ri_detect_t **detectors,
								 const uint32_t n_detectors,
								 const float peak_height,
								 const float *prefix_sum,
								 const float *prefix_sum_square,
								 uint32_t* peaks) {

	uint32_t curInd = 0;
	for (uint32_t i = 0; i < detectors[0]->s_len; i++) {
		for (uint32_t k = 0; k < n_detectors; k++) {
			ri_detect_t *detector = detectors[k];
			if (detector->masked_to >= i) continue;

			float current_value = detector->sig[i];
			// float adaptive_peak_height = calculate_adaptive_peak_height(prefix_sum, prefix_sum_square, i, detector->window_length, peak_height);

			if (detector->peak_pos == detector->DEF_PEAK_POS) {
				// CASE 1: We've not yet recorded any maximum
				if (current_value < detector->peak_value) { // A deeper minimum:
					detector->peak_value = current_value;
				} else if (current_value - detector->peak_value > peak_height) {
					// ...or a qualifying maximum:
					detector->peak_value = current_value;
					detector->peak_pos = i;
					// otherwise, wait to rise high enough to be considered a peak
				}
			} else {
				// CASE 2: In an existing peak, waiting to see if it is good
				if (current_value > detector->peak_value) {
					// Update the peak
					detector->peak_value = current_value;
					detector->peak_pos = i;
				}
				// Tell other detectors no need to check for a peak until a certain point
				if (detector->peak_value > detector->threshold) {
					for(int n_d = k+1; n_d < n_detectors; n_d++){
						detectors[n_d]->masked_to = detector->peak_pos + detectors[0]->window_length;
						detectors[n_d]->peak_pos = detectors[n_d]->DEF_PEAK_POS;
						detectors[n_d]->peak_value = detectors[n_d]->DEF_PEAK_VAL;
						detectors[n_d]->valid_peak = 0;
					}
				}
				// There is a good peak
				if (detector->peak_value - current_value > peak_height && 
					detector->peak_value > detector->threshold) {
					detector->valid_peak = 1;
				}
				// Check if we are now further away from the current peak
				if (detector->valid_peak && (i - detector->peak_pos) > detector->window_length / 2) {
					peaks[curInd++] = detector->peak_pos;
					detector->peak_pos = detector->DEF_PEAK_POS;
					detector->peak_value = current_value;
					detector->valid_peak = 0;
				}
			}
		}
	}

	return curInd;
}

int compare_floats(const void* a, const void* b) {
    const float* da = (const float*) a;
    const float* db = (const float*) b;
    return (*da > *db) - (*da < *db);
}

float calculate_mean_of_filtered_segment(float* segment,
										 const uint32_t segment_length)
{
    // Calculate median and IQR
    qsort(segment, segment_length, sizeof(float), compare_floats); // Assuming compare_floats is already defined
    float q1 = segment[segment_length / 4];
    float q3 = segment[3 * segment_length / 4];
    float iqr = q3 - q1;
    float lower_bound = q1 - iqr;
    float upper_bound = q3 + iqr;

    float sum = 0.0;
    uint32_t count = 0;
    for (uint32_t i = 0; i < segment_length; i++) {
        if (segment[i] >= lower_bound && segment[i] <= upper_bound) {
            sum += segment[i];
            ++count;
        }
    }

    // Return the mean of the filtered segment
    return count > 0 ? sum / count : 0; // Ensure we don't divide by zero
}

/**
 * @brief Generates events from peaks, prefix sums and s_len.
 * 
 * @param km Pointer to memory manager.
 * @param peaks Array of peak positions.
 * @param peak_size Size of peaks array.
 * @param prefix_sum Array of prefix sums.
 * @param s_len Length of the signal.
 * @param n_events Pointer to the number of events generated.
 * @return float* Pointer to the array of generated events.
 */
static inline float* gen_events(void *km,
								float* sig,
								const uint32_t *peaks,
								const uint32_t peak_size,
								const uint32_t s_len,
								uint32_t* n_events)
{
	uint32_t n_ev = 0;

	for (uint32_t pi = 0; pi < peak_size; ++pi)
		if (peaks[pi] > 0 && peaks[pi] < s_len) n_ev++;

	float* events = (float*)ri_kmalloc(km, n_ev*sizeof(float));

	uint32_t start_idx = 0, segment_length = 0;
	uint32_t actual_events = 0;

	for (uint32_t pi = 0; pi < peak_size; pi++){
		if (!(peaks[pi] > 0 && peaks[pi] < s_len)) continue;

    	segment_length = peaks[pi] - start_idx;
		if (segment_length > 0 && segment_length < 500) // Skip if the segment is too long
			events[actual_events++] = calculate_mean_of_filtered_segment(sig + start_idx, segment_length);
		start_idx = peaks[pi];
	}

	(*n_events) = actual_events;
	return events;
}

/* Read a float from env var with default fallback. */
static inline float env_f(const char *name, float def) {
    const char *s = getenv(name);
    if (!s || !*s) return def;
    char *end;
    float v = strtof(s, &end);
    return (end == s) ? def : v;
}
/* Read an unsigned int from env var with default fallback. */
static inline uint32_t env_u(const char *name, uint32_t def) {
    const char *s = getenv(name);
    if (!s || !*s) return def;
    char *end;
    long v = strtol(s, &end, 10);
    return (end == s || v < 0) ? def : (uint32_t)v;
}

static float* smoothed_signal(void *km, const float *sig, uint32_t n, uint32_t w);
static uint32_t pick_top_k_peaks(void *km, const float *score, uint32_t n,
                                 uint32_t target_k, uint32_t min_size,
                                 uint32_t *cps, uint32_t max_cp);
static uint32_t estimate_k_via_default_probe(void *km, const float *sig, uint32_t n,
                                              uint32_t fallback_k);

static float* detect_events_hmm(void *km, float *sig_in, uint32_t n, uint32_t *n_events) {
    if (n < 10) { *n_events = 0; return 0; }
    /* env-tunable hyperparameters (boundaries enforced inside ranges) */
    uint32_t n_states     = env_u("RH2_HMM_STATES", 4);
    if (n_states < 2) n_states = 2;
    if (n_states > 8) n_states = 8;
    float    stay_p       = env_f("RH2_HMM_STAY", 0.95f);
    if (stay_p < 0.5f) stay_p = 0.5f; if (stay_p > 0.999f) stay_p = 0.999f;
    uint32_t kmeans_iters = env_u("RH2_HMM_KMEANS_ITERS", 8);
    uint32_t presmooth    = env_u("RH2_HMM_PRESMOOTH", 0);
    float    min_std      = env_f("RH2_HMM_MIN_STD", 0.1f);
    /* RH2_HMM_TOPK > 0: extract top-K state-transition events using
     * transition log-likelihood gap as score, instead of every state change.
     * This guarantees consistent event density (same as other top-K segmenters). */
    uint32_t use_topk     = env_u("RH2_HMM_TOPK", 0);
    uint32_t evlen        = env_u("RH2_TARGET_EVENT_LEN", 9);

    /* optional pre-smoothing — helps on noisy R10 chemistries */
    float *sig = (presmooth > 0) ? smoothed_signal(km, sig_in, n, presmooth) : sig_in;

    float state_means[8], state_stds[8];
    // Initialize with quantiles (deterministic — same input → same init)
    float *sorted = (float*)ri_kmalloc(km, n * sizeof(float));
    memcpy(sorted, sig, n * sizeof(float));
    qsort(sorted, n, sizeof(float), compare_floats);
    for (uint32_t s = 0; s < n_states; s++) {
        uint32_t idx = (s + 1) * n / (n_states + 1);
        state_means[s] = sorted[idx];
        state_stds[s] = 1.0f;
    }
    ri_kfree(km, sorted);

    // K-means refinement (env-tunable iterations, default 8 — was fixed 5)
    uint32_t *assignments = (uint32_t*)ri_kcalloc(km, n, sizeof(uint32_t));
    for (uint32_t iter = 0; iter < kmeans_iters; iter++) {
        for (uint32_t i = 0; i < n; i++) {
            float best_dist = FLT_MAX;
            for (uint32_t s = 0; s < n_states; s++) {
                float d = fabsf(sig[i] - state_means[s]);
                if (d < best_dist) { best_dist = d; assignments[i] = s; }
            }
        }
        for (uint32_t s = 0; s < n_states; s++) {
            float sum = 0; uint32_t cnt = 0;
            for (uint32_t i = 0; i < n; i++) {
                if (assignments[i] == s) { sum += sig[i]; cnt++; }
            }
            if (cnt > 1) {
                state_means[s] = sum / cnt;
                float var_sum = 0;
                for (uint32_t i = 0; i < n; i++)
                    if (assignments[i] == s) var_sum += (sig[i] - state_means[s]) * (sig[i] - state_means[s]);
                state_stds[s] = sqrtf(var_sum / cnt);
                if (state_stds[s] < min_std) state_stds[s] = min_std;
            }
        }
    }

    float stay_log = logf(stay_p);
    float trans_log = logf((1.0f - stay_p) / (n_states - 1));
    float *V = (float*)ri_kmalloc(km, n * n_states * sizeof(float));
    uint32_t *path = (uint32_t*)ri_kmalloc(km, n * n_states * sizeof(uint32_t));

    for (uint32_t s = 0; s < n_states; s++) {
        float emit = -0.5f * logf(2 * 3.14159f * state_stds[s] * state_stds[s])
                     -0.5f * ((sig[0] - state_means[s]) / state_stds[s]) * ((sig[0] - state_means[s]) / state_stds[s]);
        V[s] = logf(1.0f / n_states) + emit;
    }

    for (uint32_t t = 1; t < n; t++) {
        for (uint32_t s = 0; s < n_states; s++) {
            float emit = -0.5f * logf(2 * 3.14159f * state_stds[s] * state_stds[s])
                         -0.5f * ((sig[t] - state_means[s]) / state_stds[s]) * ((sig[t] - state_means[s]) / state_stds[s]);
            float best = -FLT_MAX;
            uint32_t best_s = 0;
            for (uint32_t ps = 0; ps < n_states; ps++) {
                float v = V[(t-1)*n_states + ps] + (ps == s ? stay_log : trans_log);
                if (v > best) { best = v; best_s = ps; }
            }
            V[t*n_states + s] = best + emit;
            path[t*n_states + s] = best_s;
        }
    }

    // Backtrack
    uint32_t *states = assignments; // reuse
    float best = -FLT_MAX;
    for (uint32_t s = 0; s < n_states; s++) {
        if (V[(n-1)*n_states + s] > best) { best = V[(n-1)*n_states + s]; states[n-1] = s; }
    }
    for (int t = n - 2; t >= 0; t--)
        states[t] = path[(t+1)*n_states + states[t+1]];

    /* Optional Top-K boundary extraction: score every transition by emission
     * difference (proxy for log-likelihood-ratio). Used when matching index
     * built with RH2_HMM_TOPK=1 to enforce consistent event density. */
    if (use_topk > 0) {
        float *score = (float*)ri_kcalloc(km, n, sizeof(float));
        for (uint32_t t = 1; t < n; t++) {
            uint32_t a = states[t-1], b = states[t];
            if (a == b) { score[t] = 0; continue; }
            /* difference in emission log-likelihoods of the two states at t-1/t */
            float diff_a = (sig[t] - state_means[a]) / state_stds[a];
            float diff_b = (sig[t] - state_means[b]) / state_stds[b];
            score[t] = fabsf(diff_a*diff_a - diff_b*diff_b);
        }
        uint32_t target_k = (n > evlen) ? n / evlen : 1;
    if (env_u("RH2_TOPK_AUTO", 0)) target_k = estimate_k_via_default_probe(km, sig, n, target_k);
        uint32_t max_cp = target_k + 16;
        uint32_t *cps = (uint32_t*)ri_kcalloc(km, max_cp, sizeof(uint32_t));
        uint32_t n_cp = pick_top_k_peaks(km, score, n, target_k, 4, cps, max_cp);
        ri_kfree(km, score); ri_kfree(km, V); ri_kfree(km, path);
        float *events = (n_cp > 0) ? gen_events(km, sig_in, cps, n_cp, n, n_events)
                                    : (*n_events = 0, (float*)NULL);
        if (presmooth > 0) ri_kfree(km, sig);
        ri_kfree(km, cps); ri_kfree(km, assignments);
        return events;
    }

    ri_kfree(km, V); ri_kfree(km, path);

    // Extract breakpoints
    uint32_t *peaks = (uint32_t*)ri_kmalloc(km, n * sizeof(uint32_t));
    uint32_t n_peaks = 0;
    for (uint32_t i = 1; i < n; i++)
        if (states[i] != states[i-1]) peaks[n_peaks++] = i;

    float *events = 0;
    if (n_peaks > 0) events = gen_events(km, sig_in, peaks, n_peaks, n, n_events);
    else { *n_events = 0; }
    if (presmooth > 0) ri_kfree(km, sig);
    ri_kfree(km, peaks); ri_kfree(km, assignments);
    return events;
}

static float* detect_events_pelt(void *km, float *sig, uint32_t n, uint32_t *n_events) {
    if (n < 10) { *n_events = 0; return 0; }
    /* env-tunable BIC-style penalty multiplier */
    float pen_mult = env_f("RH2_PELT_PEN_MULT", 2.0f);
    float penalty = pen_mult * logf((float)n);
    uint32_t min_size = env_u("RH2_PELT_MIN_SIZE", 5);

    // cost[i][j] = sum of squared residuals for segment [i,j)
    // Use prefix sums for O(1) segment cost
    double *ps = (double*)ri_kcalloc(km, n+1, sizeof(double));
    double *pss = (double*)ri_kcalloc(km, n+1, sizeof(double));
    for (uint32_t i = 0; i < n; i++) {
        ps[i+1] = ps[i] + sig[i];
        pss[i+1] = pss[i] + (double)sig[i] * sig[i];
    }

    // F[j] = optimal cost for signal[0..j)
    double *F = (double*)ri_kmalloc(km, (n+1) * sizeof(double));
    int *prev = (int*)ri_kmalloc(km, (n+1) * sizeof(int));
    F[0] = -penalty;
    for (uint32_t j = 1; j <= n; j++) { F[j] = 1e30; prev[j] = 0; }

    // Candidate set (PELT pruning)
    uint32_t *cands = (uint32_t*)ri_kmalloc(km, (n+1) * sizeof(uint32_t));
    uint32_t n_cands = 1;
    cands[0] = 0;

    for (uint32_t j = min_size; j <= n; j++) {
        uint32_t new_n_cands = 0;
        for (uint32_t ci = 0; ci < n_cands; ci++) {
            uint32_t t = cands[ci];
            if (j - t < min_size) {
                // Too close for a valid segment at this j — keep candidate for future j
                cands[new_n_cands++] = t;
                continue;
            }
            double seg_sum = ps[j] - ps[t];
            double seg_sumsq = pss[j] - pss[t];
            uint32_t seg_len = j - t;
            double seg_mean = seg_sum / seg_len;
            double cost = seg_sumsq - seg_sum * seg_mean; // SSR
            double total = F[t] + cost + penalty;
            if (total < F[j]) { F[j] = total; prev[j] = t; }
            // PELT pruning
            if (F[t] + cost <= F[j]) cands[new_n_cands++] = t;
        }
        cands[new_n_cands++] = j;
        n_cands = new_n_cands;
    }

    // Backtrack changepoints
    uint32_t *cps = (uint32_t*)ri_kmalloc(km, n * sizeof(uint32_t));
    uint32_t n_cp = 0;
    int pos = n;
    while (pos > 0) {
        if (prev[pos] > 0) cps[n_cp++] = prev[pos];
        pos = prev[pos];
    }
    // Reverse
    for (uint32_t i = 0; i < n_cp / 2; i++) {
        uint32_t tmp = cps[i]; cps[i] = cps[n_cp-1-i]; cps[n_cp-1-i] = tmp;
    }

    ri_kfree(km, ps); ri_kfree(km, pss); ri_kfree(km, F); ri_kfree(km, prev);

    float *events = 0;
    if (n_cp > 0) events = gen_events(km, sig, cps, n_cp, n, n_events);
    else { *n_events = 0; }
    if (getenv("RH2_DEBUG_EVENTS"))
        fprintf(stderr, "[DBG_PELT] n_sig=%u n_cp=%u n_events=%u\n", n, n_cp, *n_events);
    ri_kfree(km, cands); ri_kfree(km, cps);
    return events;
}

/* PELT CUDA-accelerated segmenter:
 *   1. Tries to dispatch to GPU CUDA kernel via libpelt_cuda.so (real GPU PELT)
 *   2. Falls back to AVX2-vectorized + Window-PELT pruning on CPU
 *   3. Uses window-restricted candidate set (RH2_PELT_CUDA_WIN samples) for
 *      bounded inner-loop cost — empirically same F1 as full PELT on nanopore.
 */
static float* detect_events_pelt_cuda(void *km, float *sig, uint32_t n, uint32_t *n_events) {
    if (n < 10) { *n_events = 0; return 0; }
    float pen_mult = env_f("RH2_PELT_PEN_MULT", 0.05f);  /* tuned default */
    int   min_size = (int)env_u("RH2_PELT_MIN_SIZE", 3);
    int   win_size = (int)env_u("RH2_PELT_CUDA_WIN", 256);
    int   use_gpu  = (int)env_u("RH2_PELT_USE_GPU", 1);

    /* Try GPU path first if requested */
    if (use_gpu) {
        try_load_pelt_cuda();
        if (g_pelt_cuda_fn) {
            int *cps = (int*)ri_kmalloc(km, n * sizeof(int));
            int n_cp = g_pelt_cuda_fn(sig, (int)n, pen_mult, min_size, win_size, cps, (int)n);
            if (n_cp > 0) {
                uint32_t *u_cps = (uint32_t*)ri_kmalloc(km, n_cp * sizeof(uint32_t));
                for (int i = 0; i < n_cp; i++) u_cps[i] = (uint32_t)cps[i];
                float *events = gen_events(km, sig, u_cps, (uint32_t)n_cp, n, n_events);
                if (getenv("RH2_DEBUG_EVENTS"))
                    fprintf(stderr, "[DBG_PELT_CUDA_GPU] n_sig=%u n_cp=%d n_events=%u\n",
                            n, n_cp, *n_events);
                ri_kfree(km, cps); ri_kfree(km, u_cps);
                return events;
            }
            ri_kfree(km, cps);
        }
    }

    /* CPU fallback: AVX2 + Window-PELT */
    double penalty = (double)pen_mult * log((double)n);
    double *ps  = (double*)ri_kcalloc(km, n+1, sizeof(double));
    double *pss = (double*)ri_kcalloc(km, n+1, sizeof(double));
    for (uint32_t i = 0; i < n; i++) {
        ps[i+1]  = ps[i]  + sig[i];
        pss[i+1] = pss[i] + (double)sig[i] * sig[i];
    }
    double *F = (double*)ri_kmalloc(km, (n+1) * sizeof(double));
    int    *prev = (int*)ri_kmalloc(km, (n+1) * sizeof(int));
    F[0] = -penalty;
    for (uint32_t j = 1; j <= n; j++) { F[j] = 1e30; prev[j] = 0; }

    uint32_t *cands = (uint32_t*)ri_kmalloc(km, (n+1) * sizeof(uint32_t));
    uint32_t n_cands = 1; cands[0] = 0;

    for (uint32_t j = (uint32_t)min_size; j <= n; j++) {
        uint32_t new_n = 0;
        double best_total = F[j];
        int best_t = prev[j];

#if defined(__AVX2__)
        /* Vectorized inner loop: 4 candidates per iteration with double-precision AVX2.
         * For each candidate t, compute SSR cost and total = F[t] + cost + penalty.
         * After the SIMD pass, scalar finalizer for best & pruning. */
        uint32_t i = 0;
        if (n_cands >= 4) {
            __m256d v_psj  = _mm256_set1_pd(ps[j]);
            __m256d v_pssj = _mm256_set1_pd(pss[j]);
            __m256d v_pen  = _mm256_set1_pd(penalty);
            for (; i + 4 <= n_cands; i += 4) {
                /* Gather is not in basic AVX2; use scalar gather for clarity. */
                double pst[4], psst[4], Ft[4]; int tt[4];
                int valid_mask[4] = {0,0,0,0};
                for (int k = 0; k < 4; k++) {
                    uint32_t t = cands[i+k]; tt[k] = (int)t;
                    if (j - t >= (uint32_t)min_size) {
                        pst[k]  = ps[t];  psst[k] = pss[t]; Ft[k] = F[t];
                        valid_mask[k] = 1;
                    } else { pst[k]=0; psst[k]=0; Ft[k]=1e30; }
                }
                __m256d v_pst  = _mm256_loadu_pd(pst);
                __m256d v_psst = _mm256_loadu_pd(psst);
                __m256d v_Ft   = _mm256_loadu_pd(Ft);
                __m256d v_seg_sum = _mm256_sub_pd(v_psj, v_pst);
                __m256d v_seg_sq  = _mm256_sub_pd(v_pssj, v_psst);
                double seg_len_v[4] = { (double)(j-tt[0]), (double)(j-tt[1]),
                                        (double)(j-tt[2]), (double)(j-tt[3]) };
                __m256d v_len  = _mm256_loadu_pd(seg_len_v);
                __m256d v_mean = _mm256_div_pd(v_seg_sum, v_len);
                __m256d v_cost = _mm256_sub_pd(v_seg_sq, _mm256_mul_pd(v_seg_sum, v_mean));
                __m256d v_total= _mm256_add_pd(_mm256_add_pd(v_Ft, v_cost), v_pen);
                double tot[4];
                _mm256_storeu_pd(tot, v_total);
                for (int k = 0; k < 4; k++) {
                    uint32_t t = (uint32_t)tt[k];
                    if (!valid_mask[k]) {
                        cands[new_n++] = t;  /* keep for later j */
                        continue;
                    }
                    if (tot[k] < best_total) { best_total = tot[k]; best_t = (int)t; }
                    /* Window-PELT: prune if t is older than win_size from j */
                    if (j - t < (uint32_t)win_size) cands[new_n++] = t;
                }
            }
        }
        /* Tail (scalar) */
        for (; i < n_cands; i++) {
            uint32_t t = cands[i];
            if (j - t < (uint32_t)min_size) { cands[new_n++] = t; continue; }
            double seg_sum = ps[j] - ps[t];
            double seg_sq  = pss[j] - pss[t];
            double seg_len = (double)(j - t);
            double mean    = seg_sum / seg_len;
            double cost    = seg_sq - seg_sum * mean;
            double total   = F[t] + cost + penalty;
            if (total < best_total) { best_total = total; best_t = (int)t; }
            if (j - t < (uint32_t)win_size) cands[new_n++] = t;
        }
#else
        for (uint32_t i = 0; i < n_cands; i++) {
            uint32_t t = cands[i];
            if (j - t < (uint32_t)min_size) { cands[new_n++] = t; continue; }
            double seg_sum = ps[j] - ps[t];
            double seg_sq  = pss[j] - pss[t];
            double seg_len = (double)(j - t);
            double mean    = seg_sum / seg_len;
            double cost    = seg_sq - seg_sum * mean;
            double total   = F[t] + cost + penalty;
            if (total < best_total) { best_total = total; best_t = (int)t; }
            if (j - t < (uint32_t)win_size) cands[new_n++] = t;
        }
#endif
        F[j] = best_total; prev[j] = best_t;
        cands[new_n++] = j;
        n_cands = new_n;
    }

    uint32_t *cps = (uint32_t*)ri_kmalloc(km, n * sizeof(uint32_t));
    uint32_t n_cp = 0; int pos = (int)n;
    while (pos > 0) {
        if (prev[pos] > 0) cps[n_cp++] = (uint32_t)prev[pos];
        pos = prev[pos];
    }
    for (uint32_t i = 0; i < n_cp / 2; i++) {
        uint32_t tmp = cps[i]; cps[i] = cps[n_cp-1-i]; cps[n_cp-1-i] = tmp;
    }
    ri_kfree(km, ps); ri_kfree(km, pss); ri_kfree(km, F); ri_kfree(km, prev);
    float *events = 0;
    if (n_cp > 0) events = gen_events(km, sig, cps, n_cp, n, n_events);
    else { *n_events = 0; }
    if (getenv("RH2_DEBUG_EVENTS"))
        fprintf(stderr, "[DBG_PELT_CUDA_CPU] n_sig=%u n_cp=%u n_events=%u\n", n, n_cp, *n_events);
    ri_kfree(km, cands); ri_kfree(km, cps);
    return events;
}

static int cmp_uint32_asc(const void *a, const void *b) {
    uint32_t ua = *(const uint32_t*)a, ub = *(const uint32_t*)b;
    return (ua > ub) - (ua < ub);
}

static void binseg_recursive(float *sig, uint32_t start, uint32_t end,
                              double *ps, double *pss, float pen_mult, uint32_t min_size,
                              uint32_t *cps, uint32_t *n_cp, uint32_t max_cp,
                              uint32_t depth) {
    if (end - start < 2 * min_size || *n_cp >= max_cp) return;
    /* Recursion-depth cap: if depth > ~3·log2(n_total) the recursion is
     * pathologically unbalanced — abort to avoid O(n²) blow-up on long
     * reads with irregular boundary structure (RB_dmel R10.4.1). Tunable
     * via RH2_BINSEG_MAX_DEPTH (default high → off; set to ~64 to enable). */
    static int max_depth_cache = -1;
    if (max_depth_cache < 0) max_depth_cache = (int)env_u("RH2_BINSEG_MAX_DEPTH", 64);
    if ((int)depth > max_depth_cache) return;

    /* Window-restricted candidate scan: instead of scanning ALL positions
     * in [start+min_size, end-min_size] for the best split, restrict to a
     * window of size W = sqrt(seg_len) × multiplier around the segment
     * midpoint. Reduces inner-loop cost from O(n) to O(sqrt(n)) per call.
     * Total complexity: O(n × sqrt(n)) = O(n^1.5) instead of O(n²) worst.
     * Activated via RH2_BINSEG_WINSCAN=1 (default off for safety). */
    uint32_t total_len = end - start;
    uint32_t scan_lo = start + min_size;
    uint32_t scan_hi = end - min_size;
    if (env_u("RH2_BINSEG_WINSCAN", 0)) {
        uint32_t mid = (start + end) / 2;
        uint32_t half_w = (uint32_t)(sqrtf((float)total_len) * 2.0f);
        if (half_w < min_size) half_w = min_size;
        scan_lo = (mid > half_w) ? mid - half_w : start + min_size;
        scan_hi = (mid + half_w < end) ? mid + half_w : end - min_size;
        if (scan_lo < start + min_size) scan_lo = start + min_size;
        if (scan_hi > end - min_size) scan_hi = end - min_size;
    }

    double total_sum = ps[end] - ps[start];
    double total_sumsq = pss[end] - pss[start];
    double total_cost = total_sumsq - total_sum * total_sum / total_len;

    double best_gain = -1;
    uint32_t best_t = scan_lo;

    for (uint32_t t = scan_lo; t <= scan_hi; t++) {
        double left_sum = ps[t] - ps[start];
        double left_sumsq = pss[t] - pss[start];
        uint32_t left_len = t - start;
        double left_cost = left_sumsq - left_sum * left_sum / left_len;

        double right_sum = ps[end] - ps[t];
        double right_sumsq = pss[end] - pss[t];
        uint32_t right_len = end - t;
        double right_cost = right_sumsq - right_sum * right_sum / right_len;

        double gain = total_cost - left_cost - right_cost;
        if (gain > best_gain) { best_gain = gain; best_t = t; }
    }

    /* Per-segment BIC-style penalty (NOT global): penalty = pen_mult * log(seg_len).
     * Using log(global n) is wrong for nanopore — k-mer-to-k-mer transitions
     * cause SSE reductions ~5–10, but log(n=10000) * 2 ≈ 18 → all splits rejected
     * → binseg under-segments to ~O(log n) events instead of ~O(n/9). */
    float penalty = pen_mult * logf((float)total_len);
    if (best_gain > penalty) {
        cps[(*n_cp)++] = best_t;
        binseg_recursive(sig, start, best_t, ps, pss, pen_mult, min_size, cps, n_cp, max_cp, depth+1);
        binseg_recursive(sig, best_t, end, ps, pss, pen_mult, min_size, cps, n_cp, max_cp, depth+1);
    }
}

static float* detect_events_binseg(void *km, float *sig, uint32_t n, uint32_t *n_events) {
    if (n < 10) { *n_events = 0; return 0; }
    /* Penalty is now PER-SEGMENT (applied inside binseg_recursive).
     * pen_mult * log(seg_len) — small values let recursion dive deep into
     * each subsegment, producing the ~n/9 changepoints nanopore needs. */
    float pen_mult = env_f("RH2_BINSEG_PEN_MULT", 0.15f);
    uint32_t min_size = env_u("RH2_BINSEG_MIN_SIZE", 4);
    uint32_t max_cp = n / min_size;

    double *ps = (double*)ri_kcalloc(km, n+1, sizeof(double));
    double *pss = (double*)ri_kcalloc(km, n+1, sizeof(double));
    for (uint32_t i = 0; i < n; i++) {
        ps[i+1] = ps[i] + sig[i];
        pss[i+1] = pss[i] + (double)sig[i] * sig[i];
    }

    uint32_t *cps = (uint32_t*)ri_kcalloc(km, max_cp, sizeof(uint32_t));
    uint32_t n_cp = 0;
    binseg_recursive(sig, 0, n, ps, pss, pen_mult, min_size, cps, &n_cp, max_cp, 0);

    // Sort changepoints — was O(n²) bubble sort, now O(n log n) qsort.
    // Konstantin's BinSeg-slowdown investigation (May 2026): on large datasets
    // (RB_dmel: 12k reads × ~3000 cps each) bubble sort spent O(n²) per read,
    // dominating wall time. RB_dmel was 0.65× default; qsort restores parity.
    if (n_cp > 1) qsort(cps, n_cp, sizeof(uint32_t), cmp_uint32_asc);
    ri_kfree(km, ps); ri_kfree(km, pss);

    float *events = 0;
    if (n_cp > 0) events = gen_events(km, sig, cps, n_cp, n, n_events);
    else { *n_events = 0; }
    ri_kfree(km, cps);
    return events;
}

static float* detect_events_scrappie(void *km, float *sig, uint32_t n, uint32_t *n_events) {
    /**
     * Scrappie-style event detection: uses the same dual-window t-statistic
     * peak detection as the ONT default, but with Scrappie's parameters:
     *   short window = 3, long window = 6, threshold1 = 1.4,
     *   threshold2 = 9.0, peak_height = 0.4
     * These shorter windows and higher thresholds produce fewer, more
     * confident events compared to RawHash2's default (w1=8, w2=40).
     */
    if (n < 20) { *n_events = 0; return 0; }

    // Scrappie default parameters
    uint32_t w1 = 3, w2 = 6;
    float t1 = 1.4f, t2 = 9.0f, ph = 0.4f;

    float *prefix_sum = (float*)ri_kcalloc(km, n+1, sizeof(float));
    float *prefix_sum_sq = (float*)ri_kcalloc(km, n+1, sizeof(float));
    comp_prefix_prefixsq(sig, n, prefix_sum, prefix_sum_sq);

    float *tstat1 = comp_tstat(km, prefix_sum, prefix_sum_sq, n, w1);
    float *tstat2 = comp_tstat(km, prefix_sum, prefix_sum_sq, n, w2);

    ri_detect_t short_det = {.DEF_PEAK_POS = -1, .DEF_PEAK_VAL = FLT_MAX,
        .sig = tstat1, .s_len = n, .threshold = t1,
        .window_length = w1, .masked_to = 0, .peak_pos = -1,
        .peak_value = FLT_MAX, .valid_peak = 0};
    ri_detect_t long_det = {.DEF_PEAK_POS = -1, .DEF_PEAK_VAL = FLT_MAX,
        .sig = tstat2, .s_len = n, .threshold = t2,
        .window_length = w2, .masked_to = 0, .peak_pos = -1,
        .peak_value = FLT_MAX, .valid_peak = 0};

    uint32_t *peaks = (uint32_t*)ri_kmalloc(km, n * sizeof(uint32_t));
    ri_detect_t *detectors[2] = {&short_det, &long_det};
    uint32_t n_peaks = gen_peaks(detectors, 2, ph, prefix_sum, prefix_sum_sq, peaks);

    ri_kfree(km, tstat1); ri_kfree(km, tstat2);
    ri_kfree(km, prefix_sum); ri_kfree(km, prefix_sum_sq);

    float *events = 0;
    if (n_peaks > 0) events = gen_events(km, sig, peaks, n_peaks, n, n_events);
    else { *n_events = 0; }
    ri_kfree(km, peaks);
    return events;
}

/* Sliding-window Mann-Whitney U |z|-score (asymptotic, no tie correction).
 * Pairwise O(w^2) per position — for w in {3..15} this is faster than any
 * sort-based ranking (no allocation, branch-friendly inner loop).
 *
 * Under H0, U ~ Normal(mu, sigma^2) with
 *     mu        = w*w / 2
 *     sigma^2   = w*w * (2w + 1) / 12
 * Tie correction shrinks sigma^2 by sum(t^3 - t)/((2w)*(2w-1)). For
 * MAD-normalized float signals exact ties are rare, so we drop the
 * correction (slightly conservative — fewer spurious peaks).
 */
static inline float* comp_wilcoxon_z(void *km,
                                     const float *sig,
                                     const uint32_t s_len,
                                     const uint32_t w_len)
{
    float *z_arr = (float*)ri_kcalloc(km, s_len + 1, sizeof(float));
    if (s_len < 2 * w_len || w_len < 2) return z_arr;

    const float w  = (float)w_len;
    const float mu = 0.5f * w * w;
    const float var = w * w * (2.0f * w + 1.0f) / 12.0f;
    const float inv_sigma = 1.0f / sqrtf(var);

    for (uint32_t i = w_len; i + w_len <= s_len; i++) {
        float u = 0.0f;
        const float *L = sig + (i - w_len);
        const float *R = sig + i;
        for (uint32_t a = 0; a < w_len; a++) {
            float la = L[a];
            for (uint32_t b = 0; b < w_len; b++) {
                float rb = R[b];
                if (la > rb)      u += 1.0f;
                else if (la == rb) u += 0.5f;
            }
        }
        float z = (u - mu) * inv_sigma;
        z_arr[i] = z < 0.0f ? -z : z;
    }
    return z_arr;
}

/* Mann-Whitney U based dual-window segmenter. Same peak-detection
 * pipeline as the default (gen_peaks + gen_events), but the per-position
 * statistic is the |z|-score of Mann-Whitney U comparing [i-w, i] vs
 * [i, i+w] instead of Welch's t. Defaults derived from the F1 sweep in
 * ~/rawhash2/wilcox_seg/sweep_f1.py: w1=3, w2=9, t1=2.5, t2=2.0, ph=0.3.
 *
 * Tradeoff vs default: slightly slower per signal (O(n*w^2) vs O(n) of
 * tstat) but rank-based — more robust under heavy-tailed noise / outliers.
 */
static float* detect_events_wilcoxon(void *km, float *sig, uint32_t n, uint32_t *n_events) {
    if (n < 20) { *n_events = 0; return 0; }

    /* Defaults: same windows as t-stat default (w1=3, w2=9, ph=0.4) with
     * thresholds rescaled to Wilcoxon's bounded |z| range
     * (max |z| for w=3 is ~1.97, for w=9 is ~3.58 -- so t-stat's t1=4.0
     * cannot be reached at w=3; we use 1.8 and 3.0 instead). On real data:
     *   D1 F1 = 0.909, D3 F1 = 0.439, D4 F1 = 0.374 (default 0.953/0.957/0.265).
     * Override via env vars: RH2_WILCOX_{W1,W2,T1,T2,PH}. */
    uint32_t w1 = env_u("RH2_WILCOX_W1", 3);
    uint32_t w2 = env_u("RH2_WILCOX_W2", 9);
    float t1 = env_f("RH2_WILCOX_T1", 1.8f);
    float t2 = env_f("RH2_WILCOX_T2", 3.0f);
    float ph = env_f("RH2_WILCOX_PH", 0.4f);

    /* gen_peaks signature still takes prefix sums (legacy / adaptive
     * peak-height stub). Wilcoxon doesn't use them — pass zeroed buffers. */
    float *prefix_sum = (float*)ri_kcalloc(km, n + 1, sizeof(float));
    float *prefix_sum_sq = (float*)ri_kcalloc(km, n + 1, sizeof(float));

    float *z1 = comp_wilcoxon_z(km, sig, n, w1);
    float *z2 = comp_wilcoxon_z(km, sig, n, w2);

    ri_detect_t short_det = {.DEF_PEAK_POS = -1, .DEF_PEAK_VAL = FLT_MAX,
        .sig = z1, .s_len = n, .threshold = t1,
        .window_length = w1, .masked_to = 0, .peak_pos = -1,
        .peak_value = FLT_MAX, .valid_peak = 0};
    ri_detect_t long_det = {.DEF_PEAK_POS = -1, .DEF_PEAK_VAL = FLT_MAX,
        .sig = z2, .s_len = n, .threshold = t2,
        .window_length = w2, .masked_to = 0, .peak_pos = -1,
        .peak_value = FLT_MAX, .valid_peak = 0};

    uint32_t *peaks = (uint32_t*)ri_kmalloc(km, n * sizeof(uint32_t));
    ri_detect_t *detectors[2] = {&short_det, &long_det};
    uint32_t n_peaks = gen_peaks(detectors, 2, ph, prefix_sum, prefix_sum_sq, peaks);

    ri_kfree(km, z1); ri_kfree(km, z2);
    ri_kfree(km, prefix_sum); ri_kfree(km, prefix_sum_sq);

    float *events = 0;
    if (n_peaks > 0) events = gen_events(km, sig, peaks, n_peaks, n, n_events);
    else { *n_events = 0; }
    ri_kfree(km, peaks);
    return events;
}

/* ---------------------------------------------------------------------- */
/*  Additional segmenters from diverse algorithmic families                */
/* ---------------------------------------------------------------------- */

/* Pick the top-K positions of a per-position SCORE, respecting min_size
 * separation. Implements simple non-maximum suppression: greedy take the
 * highest-score not-yet-suppressed position, then mask everything within
 * min_size, repeat until target_k positions or scores exhausted.
 *
 * Output: changepoint indices in `cps[]`, sorted ascending by position.
 * Returns the number of changepoints written.
 */

/* ------------------------------------------------------------------
 * estimate_k_via_default_probe  (Konstantin's idea, May 2026)
 *
 * Run the REAL default segmenter (dual-window t-stat + gen_peaks with the
 * actual parameters w1=3, w2=9, t1=4.0, t2=3.5, peak_height=0.4) on a small
 * initial chunk of the signal and count how many peaks it produces.
 * Extrapolate to the full signal length:
 *
 *     ratio        = peaks_in_chunk / chunk_length
 *     K_estimate   = round(signal_length * ratio)
 *
 * Earlier version (v1, buggy): used a simplified single-window t-stat with
 * threshold=1.0 and local-maxima check. That undercounted peaks by 5–10×
 * because real default uses (a) two windows combined via gen_peaks, (b)
 * threshold-crossing peak logic (not local-maxima), (c) higher thresholds
 * (4.0/3.5) but compensating peak_height check. Resulting K was too small
 * and mad/window/gradient collapsed.
 *
 * This version calls the actual gen_peaks with default parameters. Cost:
 * O(probe) ~ 1–2 % of full-read default.
 *
 * Activated via env RH2_TOPK_AUTO=1.  Probe size via RH2_TOPK_PROBE.
 * Falls back to physics K on probe fail.
 * ------------------------------------------------------------------ */
static uint32_t estimate_k_via_default_probe(void *km, const float *sig, uint32_t n,
                                              uint32_t fallback_k) {
    /* Probe size: at least RH2_TOPK_PROBE_MIN samples (default 1024) but never
     * more than RH2_TOPK_PROBE_FRAC of the read (default 5 %). For long reads
     * (>20 k samples) the 5% fraction kicks in; for short reads we use the
     * 1024-sample minimum so the probe has enough statistical material.
     *
     *     probe = max(MIN, FRAC * n)   then capped at n                    */
    float frac = env_f("RH2_TOPK_PROBE_FRAC", 0.05f);
    uint32_t pmin = env_u("RH2_TOPK_PROBE_MIN", 1024);
    uint32_t probe = (uint32_t)((double)n * (double)frac);
    if (probe < pmin) probe = pmin;
    /* Optional fixed override for experimentation */
    uint32_t fixed = env_u("RH2_TOPK_PROBE", 0);
    if (fixed > 0) probe = fixed;
    if (probe > n) probe = n;
    if (probe < 64) return fallback_k;

    /* Real default segmenter parameters (from roptions.c) */
    const uint32_t w1 = 3, w2 = 9;
    const float th1 = 4.0f, th2 = 3.5f, ph = 0.4f;
    if (probe < 4 * w2) return fallback_k;

    float *ps  = (float*)ri_kcalloc(km, probe + 1, sizeof(float));
    float *pss = (float*)ri_kcalloc(km, probe + 1, sizeof(float));
    /* normalize the probe slice to mean 0, std 1 — matches what default
     * segmenter does to its full-signal input before tstat. */
    /* Use simple prefix-sum based normalization shortcut: just copy then
     * compute prefix sums on the raw signal — comp_tstat operates on
     * pre-normalization is not strictly required for peak counting. */
    for (uint32_t i = 0; i < probe; i++) {
        ps[i+1]  = ps[i]  + sig[i];
        pss[i+1] = pss[i] + sig[i] * sig[i];
    }
    float *t1 = comp_tstat(km, ps, pss, probe, w1);
    float *t2 = comp_tstat(km, ps, pss, probe, w2);
    ri_detect_t det1 = {.DEF_PEAK_POS = -1, .DEF_PEAK_VAL = FLT_MAX,
        .sig = t1, .s_len = probe, .threshold = th1,
        .window_length = w1, .masked_to = 0, .peak_pos = -1,
        .peak_value = FLT_MAX, .valid_peak = 0};
    ri_detect_t det2 = {.DEF_PEAK_POS = -1, .DEF_PEAK_VAL = FLT_MAX,
        .sig = t2, .s_len = probe, .threshold = th2,
        .window_length = w2, .masked_to = 0, .peak_pos = -1,
        .peak_value = FLT_MAX, .valid_peak = 0};
    uint32_t *peaks = (uint32_t*)ri_kmalloc(km, probe * sizeof(uint32_t));
    ri_detect_t *dets[2] = {&det1, &det2};
    uint32_t n_peaks = gen_peaks(dets, 2, ph, ps, pss, peaks);
    ri_kfree(km, ps); ri_kfree(km, pss); ri_kfree(km, t1); ri_kfree(km, t2);
    ri_kfree(km, peaks);
    if (n_peaks < 2) return fallback_k;

    double ratio = (double)n_peaks / (double)probe;
    uint32_t k_est = (uint32_t)(ratio * (double)n + 0.5);
    /* Chemistry-aware boost: window's score function on R10.4 produces fewer
     * peaks than default's t-stat at same K, due to motor-stall artifacts.
     * RH2_TOPK_AUTO_BOOST multiplies K_est before clamping. Default 1.0.
     * For R10.4 chemistries, set 1.2-1.5 to compensate. */
    float boost = env_f("RH2_TOPK_AUTO_BOOST", 1.0f);
    if (boost > 0 && boost != 1.0f) {
        k_est = (uint32_t)((double)k_est * (double)boost + 0.5);
    }

    /* Sanity-clamp: physics K is a coarse lower bound; real default usually
     * produces 1-3× physics K. Allow K_est anywhere in [K/3, 8K] so the
     * estimator can really track default's event density.
     * Tighter version (RH2_TOPK_AUTO_TIGHT=1) clamps to [K/2, 4K]. */
    uint32_t lo, hi;
    if (env_u("RH2_TOPK_AUTO_TIGHT", 0)) {
        lo = (fallback_k > 1) ? fallback_k / 2 : 1;
        hi = fallback_k * 4;
    } else {
        lo = (fallback_k > 2) ? fallback_k / 3 : 1;
        hi = fallback_k * 8;
    }
    if (k_est < lo) k_est = lo;
    if (k_est > hi) k_est = hi;
    return k_est;
}

static uint32_t pick_top_k_peaks(void *km, const float *score, uint32_t n,
                                 uint32_t target_k, uint32_t min_size,
                                 uint32_t *cps, uint32_t max_cp) {
    if (target_k > max_cp) target_k = max_cp;
    if (target_k == 0 || n == 0) return 0;
    /* taken[i] = 1 if position i is suppressed */
    char *taken = (char*)ri_kcalloc(km, n, sizeof(char));
    /* mark boundary regions as un-pickable */
    uint32_t bnd = (min_size > 0) ? min_size : 1;
    for (uint32_t i = 0; i < bnd && i < n; i++) taken[i] = 1;
    for (uint32_t i = (n > bnd ? n - bnd : 0); i < n; i++) taken[i] = 1;
    uint32_t n_cp = 0;
    while (n_cp < target_k) {
        /* find argmax over !taken */
        float best = -INFINITY;
        uint32_t best_i = n;
        for (uint32_t i = 0; i < n; i++) {
            if (!taken[i] && score[i] > best) { best = score[i]; best_i = i; }
        }
        if (best_i >= n || !isfinite(best)) break;
        cps[n_cp++] = best_i;
        /* suppress min_size on each side */
        uint32_t lo = (best_i >= min_size) ? best_i - min_size : 0;
        uint32_t hi = (best_i + min_size + 1 < n) ? best_i + min_size + 1 : n;
        for (uint32_t i = lo; i < hi; i++) taken[i] = 1;
    }
    ri_kfree(km, taken);
    /* sort cps ascending */
    for (uint32_t i = 0; i + 1 < n_cp; i++)
        for (uint32_t j = i + 1; j < n_cp; j++)
            if (cps[i] > cps[j]) { uint32_t t = cps[i]; cps[i] = cps[j]; cps[j] = t; }
    return n_cp;
}

/* Helper: convert a sorted array of changepoint indices into events. */
static inline float* finalize_changepoints(void *km, float *sig, uint32_t n,
                                           uint32_t *cps, uint32_t n_cp,
                                           uint32_t *n_events) {
    if (n_cp == 0) { *n_events = 0; return 0; }
    /* sort */
    for (uint32_t i = 0; i + 1 < n_cp; i++)
        for (uint32_t j = i + 1; j < n_cp; j++)
            if (cps[i] > cps[j]) { uint32_t t = cps[i]; cps[i] = cps[j]; cps[j] = t; }
    return gen_events(km, sig, cps, n_cp, n, n_events);
}

/* Moving-average smoothing in-place, using a separate buffer. Window radius w. */
static inline float* smoothed_signal(void *km, const float *sig, uint32_t n, uint32_t w) {
    float *out = (float*)ri_kmalloc(km, n * sizeof(float));
    if (w == 0) { memcpy(out, sig, n * sizeof(float)); return out; }
    double s = 0;
    /* prefix sum approach for O(n) */
    double *ps = (double*)ri_kcalloc(km, n + 1, sizeof(double));
    for (uint32_t i = 0; i < n; i++) ps[i+1] = ps[i] + sig[i];
    for (uint32_t i = 0; i < n; i++) {
        uint32_t lo = (i >= w) ? i - w : 0;
        uint32_t hi = (i + w + 1 < n) ? i + w + 1 : n;
        out[i] = (float)((ps[hi] - ps[lo]) / (hi - lo));
    }
    ri_kfree(km, ps);
    (void)s;
    return out;
}

/* Merge consecutive segments whose means differ by less than `merge_thr`,
 * producing a sparser, more level-aligned set of changepoints. */
static inline uint32_t merge_close_segments(const float *sig, uint32_t n,
                                            uint32_t *cps, uint32_t n_cp,
                                            float merge_thr) {
    if (n_cp < 2 || merge_thr <= 0) return n_cp;
    /* compute segment means */
    uint32_t prev = 0; uint32_t out = 0;
    float prev_mean = 0; uint32_t prev_len = 0;
    {
        double s = 0;
        for (uint32_t i = 0; i < cps[0]; i++) s += sig[i];
        prev_mean = (cps[0] > 0) ? (float)(s / cps[0]) : 0;
        prev_len = cps[0];
    }
    /* iterate through cps, merging if delta below threshold */
    for (uint32_t k = 0; k < n_cp; k++) {
        uint32_t e = (k + 1 < n_cp) ? cps[k + 1] : n;
        double s = 0;
        for (uint32_t i = cps[k]; i < e; i++) s += sig[i];
        uint32_t len = e - cps[k];
        float mean = (len > 0) ? (float)(s / len) : 0;
        if (fabsf(mean - prev_mean) < merge_thr) {
            /* merge: update running mean */
            prev_mean = (prev_mean * prev_len + mean * len) / (prev_len + len);
            prev_len += len;
        } else {
            cps[out++] = cps[k];
            prev_mean = mean;
            prev_len = len;
            (void)prev;
        }
    }
    return out;
}

/* ---- 6. CUSUM (Page's cumulative sum) ---------------------------------- *
 * Classic statistical process control. Tracks the running cumulative sum of
 * deviations from a reference; resets on threshold crossing, marking a CP.
 * Family: classical change-point statistics (orthogonal to t-test).
 */
static float* detect_events_cusum(void *km, float *sig, uint32_t n, uint32_t *n_events) {
    if (n < 20) { *n_events = 0; return 0; }
    /* Hybrid CUSUM-like detector that actually works for nanopore:
     * apply CUSUM-style accumulation on the FIRST DIFFERENCE of the signal
     * (captures direction CHANGES, not absolute levels). The score per
     * position is the max(g+, g-) of this differential CUSUM, then top-K. */
    const float k        = env_f("RH2_CUSUM_K", 0.4f);
    const uint32_t min_size = env_u("RH2_CUSUM_MIN_SIZE", 4);
    const uint32_t presmooth = env_u("RH2_CUSUM_PRESMOOTH", 4);
    const uint32_t evlen = env_u("RH2_TARGET_EVENT_LEN", 9);
    float *src = (presmooth > 0) ? smoothed_signal(km, sig, n, presmooth) : sig;
    /* differential CUSUM on first-difference */
    float *score = (float*)ri_kcalloc(km, n, sizeof(float));
    float gp = 0, gn = 0;
    for (uint32_t i = 1; i < n; i++) {
        float d = src[i] - src[i-1];   /* first difference */
        gp = gp + d - k; if (gp < 0) gp = 0;
        gn = gn - d - k; if (gn < 0) gn = 0;
        score[i] = (gp > gn) ? gp : gn;
    }
    uint32_t target_k = (n > evlen) ? n / evlen : 1;
    if (env_u("RH2_TOPK_AUTO", 0)) target_k = estimate_k_via_default_probe(km, sig, n, target_k);
    uint32_t max_cp = target_k + 16;
    uint32_t *cps = (uint32_t*)ri_kcalloc(km, max_cp, sizeof(uint32_t));
    uint32_t n_cp = pick_top_k_peaks(km, score, n, target_k, min_size, cps, max_cp);
    ri_kfree(km, score);
    float *events = (n_cp > 0) ? gen_events(km, sig, cps, n_cp, n, n_events)
                                : (*n_events = 0, (float*)NULL);
    if (presmooth > 0) ri_kfree(km, src);
    ri_kfree(km, cps);
    return events;
}

/* ---- 7. Gradient (1st-derivative threshold + non-max suppression) ------ *
 * Edge-detection family. Smooths derivative magnitude over a small kernel
 * then picks local maxima above a threshold. Cheap and fast.
 */
static float* detect_events_gradient(void *km, float *sig, uint32_t n, uint32_t *n_events) {
    if (n < 20) { *n_events = 0; return 0; }
    const uint32_t smooth = env_u("RH2_GRAD_SMOOTH", 5);
    const uint32_t min_size = env_u("RH2_GRAD_MIN_SIZE", 4);
    const uint32_t presmooth = env_u("RH2_GRAD_PRESMOOTH", 4);
    const uint32_t evlen = env_u("RH2_TARGET_EVENT_LEN", 9);
    float *src = (presmooth > 0) ? smoothed_signal(km, sig, n, presmooth) : sig;
    /* Per-position score = magnitude of the smoothed gradient (left-mean vs right-mean). */
    float *score = (float*)ri_kcalloc(km, n, sizeof(float));
    for (uint32_t i = smooth; i + smooth < n; i++) {
        float a = 0, b = 0;
        for (uint32_t j = 0; j < smooth; j++) { a += src[i - 1 - j]; b += src[i + j]; }
        score[i] = fabsf(b - a) / smooth;
    }
    uint32_t target_k = (n > evlen) ? n / evlen : 1;
    if (env_u("RH2_TOPK_AUTO", 0)) target_k = estimate_k_via_default_probe(km, sig, n, target_k);
    uint32_t max_cp = target_k + 16;
    uint32_t *cps = (uint32_t*)ri_kcalloc(km, max_cp, sizeof(uint32_t));
    uint32_t n_cp = pick_top_k_peaks(km, score, n, target_k, min_size, cps, max_cp);
    ri_kfree(km, score);
    float *events = (n_cp > 0) ? gen_events(km, sig, cps, n_cp, n, n_events)
                                : (*n_events = 0, (float*)NULL);
    if (presmooth > 0) ri_kfree(km, src);
    ri_kfree(km, cps);
    return events;
}

/* ---- 8. MAD (Median Absolute Deviation outlier-based) ------------------ *
 * Robust statistics family. Uses the sliding-window MAD as the dispersion
 * estimator and flags points whose deviation exceeds k*MAD as boundaries.
 * Robust to heavy-tailed noise unlike t-test.
 */
static int cmp_float_asc(const void *a, const void *b) {
    float fa = *(const float*)a, fb = *(const float*)b;
    return (fa > fb) - (fa < fb);
}
static float* detect_events_mad(void *km, float *sig, uint32_t n, uint32_t *n_events) {
    if (n < 30) { *n_events = 0; return 0; }
    /* Robust two-window CPD: per-position score = |medianL - medianR| / pooled MAD. */
    const uint32_t W = env_u("RH2_MAD_W", 5);
    const uint32_t min_size = env_u("RH2_MAD_MIN_SIZE", 3);
    const uint32_t evlen = env_u("RH2_TARGET_EVENT_LEN", 9);
    float *score = (float*)ri_kcalloc(km, n, sizeof(float));
    float *bufL = (float*)ri_kmalloc(km, W * sizeof(float));
    float *bufR = (float*)ri_kmalloc(km, W * sizeof(float));
    for (uint32_t i = W; i + W < n; i++) {
        memcpy(bufL, sig + i - W, W * sizeof(float));
        memcpy(bufR, sig + i,     W * sizeof(float));
        qsort(bufL, W, sizeof(float), cmp_float_asc);
        qsort(bufR, W, sizeof(float), cmp_float_asc);
        float medL = bufL[W/2], medR = bufR[W/2];
        for (uint32_t j = 0; j < W; j++) bufL[j] = fabsf(bufL[j] - medL);
        for (uint32_t j = 0; j < W; j++) bufR[j] = fabsf(bufR[j] - medR);
        qsort(bufL, W, sizeof(float), cmp_float_asc);
        qsort(bufR, W, sizeof(float), cmp_float_asc);
        float madL = bufL[W/2], madR = bufR[W/2];
        float pooled = 1.4826f * 0.5f * (madL + madR) + 1e-6f;
        score[i] = fabsf(medR - medL) / pooled;
    }
    ri_kfree(km, bufL); ri_kfree(km, bufR);
    uint32_t target_k = (n > evlen) ? n / evlen : 1;
    if (env_u("RH2_TOPK_AUTO", 0)) target_k = estimate_k_via_default_probe(km, sig, n, target_k);
    uint32_t max_cp = target_k + 16;
    uint32_t *cps = (uint32_t*)ri_kcalloc(km, max_cp, sizeof(uint32_t));
    uint32_t n_cp = pick_top_k_peaks(km, score, n, target_k, min_size, cps, max_cp);
    ri_kfree(km, score);
    float *events = (n_cp > 0) ? gen_events(km, sig, cps, n_cp, n, n_events)
                                : (*n_events = 0, (float*)NULL);
    ri_kfree(km, cps);
    return events;
}

/* ---- 9. Bayesian Online Changepoint Detection (BOCD, Adams & MacKay) --- *
 * Different probabilistic family from HMM. Maintains a posterior over the
 * "run length" since last changepoint with a constant hazard rate; declares
 * a CP when the most likely run length resets. Uses Gaussian likelihood
 * with online running mean/variance.
 *
 * Fast-and-simple variant: track only the posterior mode (run length),
 * not the full distribution — sufficient for boundary detection in O(n).
 */
static float* detect_events_bocd(void *km, float *sig, uint32_t n, uint32_t *n_events) {
    if (n < 20) { *n_events = 0; return 0; }
    /* Expected segment length 1/hazard. Default 250 ~= 25 bases of signal at 4 kHz. */
    const float hazard = 1.0f / env_f("RH2_BOCD_EXP_LEN", 250.0f);
    const uint32_t min_size = env_u("RH2_BOCD_MIN_SIZE", 5);
    const float prior_var = env_f("RH2_BOCD_PRIOR_VAR", 1.0f);
    uint32_t max_cp = n / min_size + 1;
    uint32_t *cps = (uint32_t*)ri_kcalloc(km, max_cp, sizeof(uint32_t));
    uint32_t n_cp = 0;
    /* running stats over current segment */
    uint32_t r = 0;        /* run length */
    float m = 0, m2 = 0;   /* running mean, sum of squared deviations */
    uint32_t last_cp = 0;
    for (uint32_t i = 0; i < n; i++) {
        r += 1;
        /* Welford update */
        float delta = sig[i] - m;
        m += delta / r;
        m2 += delta * (sig[i] - m);
        float var = (r > 1 ? m2 / (r - 1) : prior_var) + 1e-6f;
        float z2 = (sig[i] - m) * (sig[i] - m) / var;
        /* growth log-likelihood vs. CP log-likelihood (under prior) */
        float ll_grow = -0.5f * (logf(2.0f * 3.14159265f * var) + z2) + logf(1.0f - hazard);
        float ll_cp   = -0.5f * (logf(2.0f * 3.14159265f * prior_var) + sig[i]*sig[i] / prior_var) + logf(hazard);
        if (ll_cp > ll_grow && (i - last_cp) >= min_size && n_cp < max_cp) {
            cps[n_cp++] = i; last_cp = i;
            r = 0; m = 0; m2 = 0;
        }
    }
    float *events = finalize_changepoints(km, sig, n, cps, n_cp, n_events);
    ri_kfree(km, cps);
    return events;
}

/* ---- 10. Sliding-Window z-score (baseline / null-control) -------------- *
 * The simplest possible non-trivial method: compute a sliding z-score and
 * mark points where |z| crosses a threshold. Useful as a baseline against
 * which any "smarter" method must be measured.
 */
static float* detect_events_window(void *km, float *sig, uint32_t n, uint32_t *n_events) {
    if (n < 30) { *n_events = 0; return 0; }
    /* Sliding two-window mean-shift CPD: per-position score = z-score of
     * left-window mean vs right-window mean. Top-K selection forces a
     * matched event density (≈ default) instead of brittle thresholding.
     * Defaults tuned for R10.4.x: small windows, no presmooth (raw signal). */
    const uint32_t W = env_u("RH2_WIN_W", 4);
    const uint32_t min_size = env_u("RH2_WIN_MIN_SIZE", 3);
    const uint32_t presmooth = env_u("RH2_WIN_PRESMOOTH", 0);
    const uint32_t evlen = env_u("RH2_TARGET_EVENT_LEN", 9);
    float *src = (presmooth > 0) ? smoothed_signal(km, sig, n, presmooth) : sig;
    /* prefix sums on the (smoothed) source */
    double *ps = (double*)ri_kcalloc(km, n+1, sizeof(double));
    double *pss = (double*)ri_kcalloc(km, n+1, sizeof(double));
    for (uint32_t i = 0; i < n; i++) {
        ps[i+1] = ps[i] + src[i];
        pss[i+1] = pss[i] + (double)src[i] * src[i];
    }
    float *score = (float*)ri_kcalloc(km, n, sizeof(float));
    for (uint32_t i = W; i + W < n; i++) {
        double sum_l = ps[i] - ps[i - W];
        double sum_r = ps[i + W] - ps[i];
        double sq_l  = pss[i] - pss[i - W];
        double sq_r  = pss[i + W] - pss[i];
        double mean_l = sum_l / W, mean_r = sum_r / W;
        double var_l = sq_l / W - mean_l * mean_l;
        double var_r = sq_r / W - mean_r * mean_r;
        double pooled = sqrt(0.5 * (var_l + var_r) + 1e-6);
        score[i] = (float)(fabs(mean_r - mean_l) / pooled);
    }
    ri_kfree(km, ps); ri_kfree(km, pss);
    uint32_t target_k = (n > evlen) ? n / evlen : 1;
    if (env_u("RH2_TOPK_AUTO", 0)) target_k = estimate_k_via_default_probe(km, sig, n, target_k);
    uint32_t max_cp = target_k + 16;
    uint32_t *cps = (uint32_t*)ri_kcalloc(km, max_cp, sizeof(uint32_t));
    uint32_t n_cp = pick_top_k_peaks(km, score, n, target_k, min_size, cps, max_cp);
    ri_kfree(km, score);
    float *events = (n_cp > 0) ? gen_events(km, sig, cps, n_cp, n, n_events)
                                : (*n_events = 0, (float*)NULL);
    if (presmooth > 0) ri_kfree(km, src);
    ri_kfree(km, cps);
    return events;
}

/* ------------------------------------------------------------------
 *  Persistent Python segmenter worker (one subprocess, mutex-serialised)
 * ------------------------------------------------------------------
 * The previous implementation forked a Python interpreter per read,
 * which made any library-grade plugin (ruptures, claspy, ...) infeasible
 * on million-read datasets due to ~50 ms cold-start overhead.
 *
 * This persistent variant launches one long-lived Python worker on first
 * call and serialises subsequent calls via a mutex. Throughput at scale
 * improves by ~50-100x. The worker protocol is:
 *
 *   C  -> Py:  "%u\n"   (signal length n; n==0 = quit)
 *   C  -> Py:  n binary float32 (little-endian)
 *   Py -> C:  "%u\n"   (number of events m)
 *   Py -> C:  m binary float32
 *
 * Python side must run a loop reading from stdin / writing to stdout
 * with line-buffered text + binary signal/events.  See
 * scripts/gt_pipeline/python_segmenter_modern.py (persistent_main()).
 */
#include <pthread.h>

static pthread_mutex_t g_py_mtx = PTHREAD_MUTEX_INITIALIZER;
static FILE *g_py_in = NULL;   /* C -> Py */
static FILE *g_py_out = NULL;  /* Py -> C */
static pid_t g_py_pid = 0;

static int spawn_python_worker(const char *script_path) {
    int pipe_in[2], pipe_out[2];
    if (pipe(pipe_in) < 0 || pipe(pipe_out) < 0) return -1;
    pid_t pid = fork();
    if (pid < 0) return -1;
    if (pid == 0) {
        close(pipe_in[1]); close(pipe_out[0]);
        dup2(pipe_in[0], STDIN_FILENO);
        dup2(pipe_out[1], STDOUT_FILENO);
        close(pipe_in[0]); close(pipe_out[1]);
        /* persistent mode signaled via env var */
        setenv("RAWHASH_SEG_PERSISTENT", "1", 1);
        execlp("python3", "python3", script_path, (char*)NULL);
        _exit(1);
    }
    close(pipe_in[0]); close(pipe_out[1]);
    g_py_in = fdopen(pipe_in[1], "w");
    g_py_out = fdopen(pipe_out[0], "r");
    g_py_pid = pid;
    if (!g_py_in || !g_py_out) return -1;
    setvbuf(g_py_in, NULL, _IOFBF, 1 << 16);
    return 0;
}

/* ---- 11. Wavelet (Haar multi-scale detail magnitude) ------------------ *
 * Multi-scale "à trous"-style Haar wavelet: per-position score = sum of
 * |signal[i+w] - signal[i-w]| over scales w ∈ {2,4,8,16}. Captures
 * boundaries that are sharp at multiple scales, suppresses single-sample
 * noise spikes. Real-time friendly (causal version uses only past samples).
 */
static float* detect_events_wavelet(void *km, float *sig, uint32_t n, uint32_t *n_events) {
    if (n < 40) { *n_events = 0; return 0; }
    const uint32_t min_size = env_u("RH2_WAVELET_MIN_SIZE", 4);
    const uint32_t presmooth = env_u("RH2_WAVELET_PRESMOOTH", 0);
    const uint32_t evlen = env_u("RH2_TARGET_EVENT_LEN", 9);
    float *src = (presmooth > 0) ? smoothed_signal(km, sig, n, presmooth) : sig;
    float *score = (float*)ri_kcalloc(km, n, sizeof(float));
    /* Multi-scale Haar detail magnitude, normalised by local std (z-score-like).
     * Without normalisation, raw amplitudes dominate score → events crowd in
     * high-noise regions instead of at real boundaries. */
    const uint32_t scales[] = {2, 4, 8, 16};
    const uint32_t n_scales = sizeof(scales)/sizeof(*scales);
    /* prefix sums for fast local std */
    double *ps = (double*)ri_kcalloc(km, n+1, sizeof(double));
    double *pss = (double*)ri_kcalloc(km, n+1, sizeof(double));
    for (uint32_t i = 0; i < n; i++) {
        ps[i+1] = ps[i] + src[i];
        pss[i+1] = pss[i] + (double)src[i]*src[i];
    }
    for (uint32_t i = 16; i + 16 < n; i++) {
        float s = 0;
        for (uint32_t k = 0; k < n_scales; k++) {
            uint32_t w = scales[k];
            double sum_l = ps[i] - ps[i-w], sum_r = ps[i+w] - ps[i];
            double sq_l  = pss[i] - pss[i-w], sq_r = pss[i+w] - pss[i];
            double mean_l = sum_l/w, mean_r = sum_r/w;
            double var_l = sq_l/w - mean_l*mean_l;
            double var_r = sq_r/w - mean_r*mean_r;
            double pooled = sqrt(0.5*(var_l + var_r) + 1e-6);
            s += (float)(fabs(mean_r - mean_l) / pooled);
        }
        score[i] = s;
    }
    ri_kfree(km, ps); ri_kfree(km, pss);
    uint32_t target_k = (n > evlen) ? n / evlen : 1;
    if (env_u("RH2_TOPK_AUTO", 0)) target_k = estimate_k_via_default_probe(km, sig, n, target_k);
    uint32_t max_cp = target_k + 16;
    uint32_t *cps = (uint32_t*)ri_kcalloc(km, max_cp, sizeof(uint32_t));
    uint32_t n_cp = pick_top_k_peaks(km, score, n, target_k, min_size, cps, max_cp);
    ri_kfree(km, score);
    float *events = (n_cp > 0) ? gen_events(km, sig, cps, n_cp, n, n_events)
                                : (*n_events = 0, (float*)NULL);
    if (presmooth > 0) ri_kfree(km, src);
    ri_kfree(km, cps);
    return events;
}

/* ---- 12. TDA (Topological persistence on 1D function) ----------------- *
 * 1D persistent homology = elder rule on critical points of the signal.
 * Score per position = topological persistence of the local extremum at
 * that position (= lifetime in sublevel filtration). Boundaries with high
 * persistence correspond to large mean shifts and are noise-robust.
 *
 * Simplified algorithm: Find all local mins/maxes; pair each (min, max)
 * by sweeping in value order; persistence = max_value - min_value.
 */
/* Compute 1D persistent homology via elder-rule:
 *   1. Find all local minima (in sublevel filtration).
 *   2. Each minimum is "born" at its value; "dies" when merged into a
 *      neighboring component at a saddle (local max). Persistence = saddle - min.
 *   3. Each min has a position and a persistence; boundaries are at min positions.
 * Score per position is the persistence of the min that lives at that position. */
static void compute_1d_persistence(const float *sig, uint32_t n, float *score_out) {
    /* Find local minima */
    for (uint32_t i = 1; i + 1 < n; i++) {
        if (sig[i] <= sig[i-1] && sig[i] < sig[i+1]) {
            /* expand left/right until we find a saddle (lower neighbor goes back up) */
            float left_max = sig[i], right_max = sig[i];
            for (int j = (int)i - 1; j >= 0; j--) {
                if (sig[j] > left_max) left_max = sig[j];
                /* stop if we encounter another min lower than i (we'd be killed by it) */
                if (sig[j] < sig[i]) break;
            }
            for (uint32_t j = i + 1; j < n; j++) {
                if (sig[j] > right_max) right_max = sig[j];
                if (sig[j] < sig[i]) break;
            }
            float death = (left_max < right_max) ? left_max : right_max;
            float persistence = death - sig[i];
            if (persistence > score_out[i]) score_out[i] = persistence;
        }
    }
    /* Symmetric: also score local maxima by their persistence (super-level) */
    for (uint32_t i = 1; i + 1 < n; i++) {
        if (sig[i] >= sig[i-1] && sig[i] > sig[i+1]) {
            float left_min = sig[i], right_min = sig[i];
            for (int j = (int)i - 1; j >= 0; j--) {
                if (sig[j] < left_min) left_min = sig[j];
                if (sig[j] > sig[i]) break;
            }
            for (uint32_t j = i + 1; j < n; j++) {
                if (sig[j] < right_min) right_min = sig[j];
                if (sig[j] > sig[i]) break;
            }
            float birth = (left_min > right_min) ? left_min : right_min;
            float persistence = sig[i] - birth;
            if (persistence > score_out[i]) score_out[i] = persistence;
        }
    }
}

static float* detect_events_tda(void *km, float *sig, uint32_t n, uint32_t *n_events) {
    if (n < 40) { *n_events = 0; return 0; }
    const uint32_t min_size = env_u("RH2_TDA_MIN_SIZE", 4);
    const uint32_t presmooth = env_u("RH2_TDA_PRESMOOTH", 3);
    const uint32_t evlen = env_u("RH2_TARGET_EVENT_LEN", 9);
    float *src = (presmooth > 0) ? smoothed_signal(km, sig, n, presmooth) : sig;
    float *score = (float*)ri_kcalloc(km, n, sizeof(float));
    /* Real 1D persistent homology: score is the persistence (lifetime in
     * sublevel filtration) of each local extremum. Boundaries with high
     * persistence correspond to genuine k-mer transitions; small
     * fluctuations have low persistence and are filtered by top-K. */
    compute_1d_persistence(src, n, score);
    uint32_t target_k = (n > evlen) ? n / evlen : 1;
    if (env_u("RH2_TOPK_AUTO", 0)) target_k = estimate_k_via_default_probe(km, sig, n, target_k);
    uint32_t max_cp = target_k + 16;
    uint32_t *cps = (uint32_t*)ri_kcalloc(km, max_cp, sizeof(uint32_t));
    uint32_t n_cp = pick_top_k_peaks(km, score, n, target_k, min_size, cps, max_cp);
    ri_kfree(km, score);
    float *events = (n_cp > 0) ? gen_events(km, sig, cps, n_cp, n, n_events)
                                : (*n_events = 0, (float*)NULL);
    if (presmooth > 0) ri_kfree(km, src);
    ri_kfree(km, cps);
    return events;
}

/* ---- 13. ConvBank (light "CNN-like" filter bank) ---------------------- *
 * Inspired by the first conv layer of a 1D-CNN edge detector but with
 * fixed, hand-designed weights — no training, fully deterministic.
 * Three filter shapes at three scales = 9 features per position; max-pool
 * across them gives the per-position score. Top-K then picks events.
 *
 * Filters (centered, antisymmetric "edge" + symmetric "Gaussian curvature"):
 *   - edge_5  = [-1, -1, 0, 1, 1] (smooth gradient)
 *   - edge_9  = [-1, -1, -1, -1, 0, 1, 1, 1, 1]
 *   - peak_7  = [-1, 1, 2, 4, 2, 1, -1] / scale (approximate Mexican-hat)
 * Each at scales 1, 2, 4 = 9 outputs total.
 */
static float* detect_events_convbank(void *km, float *sig, uint32_t n, uint32_t *n_events) {
    if (n < 40) { *n_events = 0; return 0; }
    const uint32_t min_size = env_u("RH2_CONVBANK_MIN_SIZE", 4);
    const uint32_t presmooth = env_u("RH2_CONVBANK_PRESMOOTH", 2);
    const uint32_t evlen = env_u("RH2_TARGET_EVENT_LEN", 9);
    float *src = (presmooth > 0) ? smoothed_signal(km, sig, n, presmooth) : sig;
    float *score = (float*)ri_kcalloc(km, n, sizeof(float));

    /* "CNN-like" edge-detector bank: 4 STEP filters at scales 1/2/4/8.
     * Each filter is a centered first-difference (DoG-style); we square
     * each response (gives positive saliency) and sum across scales for a
     * stable multi-scale edge magnitude. Plus a normalising local-std
     * divider. This is essentially the first conv layer of a 1D-CNN edge
     * detector with hand-crafted weights. */
    const uint32_t scales[] = {1, 2, 4, 8};
    /* prefix sums for cheap local std */
    double *ps = (double*)ri_kcalloc(km, n+1, sizeof(double));
    double *pss = (double*)ri_kcalloc(km, n+1, sizeof(double));
    for (uint32_t i = 0; i < n; i++) {
        ps[i+1] = ps[i] + src[i];
        pss[i+1] = pss[i] + (double)src[i]*src[i];
    }
    const uint32_t LSW = 16;  /* local std window */
    for (uint32_t i = LSW; i + LSW < n; i++) {
        /* local std for normalisation */
        double sumW = ps[i+LSW] - ps[i-LSW];
        double sqW  = pss[i+LSW] - pss[i-LSW];
        double meanW = sumW / (2*LSW);
        double varW = sqW / (2*LSW) - meanW*meanW;
        double sigma = sqrt(varW + 1e-6);
        double s = 0;
        for (uint32_t k = 0; k < 4; k++) {
            uint32_t w = scales[k];
            float lhs = 0, rhs = 0;
            for (uint32_t j = 0; j < w; j++) { lhs += src[i - w + j]; rhs += src[i + j]; }
            float diff = (rhs - lhs) / w;
            s += (diff*diff) / (sigma*sigma + 1e-9);
        }
        score[i] = (float)sqrt(s);
    }
    ri_kfree(km, ps); ri_kfree(km, pss);
    uint32_t target_k = (n > evlen) ? n / evlen : 1;
    if (env_u("RH2_TOPK_AUTO", 0)) target_k = estimate_k_via_default_probe(km, sig, n, target_k);
    uint32_t max_cp = target_k + 16;
    uint32_t *cps = (uint32_t*)ri_kcalloc(km, max_cp, sizeof(uint32_t));
    uint32_t n_cp = pick_top_k_peaks(km, score, n, target_k, min_size, cps, max_cp);
    ri_kfree(km, score);
    float *events = (n_cp > 0) ? gen_events(km, sig, cps, n_cp, n, n_events)
                                : (*n_events = 0, (float*)NULL);
    if (presmooth > 0) ri_kfree(km, src);
    ri_kfree(km, cps);
    return events;
}

/* ---- 14. Fusion (ML-inspired multi-feature CPD) ----------------------- *
 * Goal: beat default while staying within 10% of its runtime. Approach:
 *   1. Compute a feature vector per position from cheap O(n) statistics:
 *        f1 = local t-stat (dual-window like default)
 *        f2 = first-difference magnitude (gradient)
 *        f3 = local curvature: |signal[i] - 0.5*(signal[i-w] + signal[i+w])|
 *   2. Combine via a learned-style weighted L2-norm:
 *        score[i] = sqrt(α·f1² + β·f2² + γ·f3²)
 *      where weights (α, β, γ) are hand-tuned defaults but env-overridable.
 *      The weights bias toward features that carry independent info about
 *      a real k-mer transition (level shift + sharp gradient + curvature).
 *   3. Top-K NMS for density control (same as our other modern segmenters).
 *
 * This is not a trained NN — it is a hand-designed feature fusion that
 * mimics the first layer of a 1D-CNN edge detector, with weights chosen
 * so each feature contributes complementary information.
 */
static float* detect_events_fusion(void *km, float *sig, uint32_t n, uint32_t *n_events) {
    if (n < 30) { *n_events = 0; return 0; }
    /* "fusion" = default's exact comp_tstat scoring (short + long windows)
     * combined via MAX, then Top-K NMS instead of threshold-based gen_peaks.
     * Uses default's optimised float prefix sums and t-stat helper to keep
     * the runtime within ~10% of default. */
    const uint32_t min_size = env_u("RH2_FUSION_MIN_SIZE", 4);
    const uint32_t evlen = env_u("RH2_TARGET_EVENT_LEN", 9);
    const uint32_t W1 = env_u("RH2_FUSION_W1", 3);
    const uint32_t W2 = env_u("RH2_FUSION_W2", 9);
    const float curv_w = env_f("RH2_FUSION_CURV", 0.15f);

    /* Reuse default's optimised float prefix sums + comp_tstat helper. */
    float *prefix_sum = (float*)ri_kcalloc(km, n+1, sizeof(float));
    float *prefix_sum_square = (float*)ri_kcalloc(km, n+1, sizeof(float));
    comp_prefix_prefixsq(sig, n, prefix_sum, prefix_sum_square);
    float *t1 = comp_tstat(km, prefix_sum, prefix_sum_square, n, W1);
    float *t2 = comp_tstat(km, prefix_sum, prefix_sum_square, n, W2);

    float *score = (float*)ri_kcalloc(km, n, sizeof(float));
    for (uint32_t i = 1; i + 1 < n; i++) {
        float a = t1[i], b = t2[i];
        float mx = (a > b) ? a : b;
        float curv = fabsf(sig[i+1] - 2*sig[i] + sig[i-1]);
        score[i] = mx + curv_w * curv;
    }

    ri_kfree(km, prefix_sum); ri_kfree(km, prefix_sum_square);
    ri_kfree(km, t1); ri_kfree(km, t2);

    uint32_t target_k = (n > evlen) ? n / evlen : 1;
    if (env_u("RH2_TOPK_AUTO", 0)) target_k = estimate_k_via_default_probe(km, sig, n, target_k);
    uint32_t max_cp = target_k + 16;
    uint32_t *cps = (uint32_t*)ri_kcalloc(km, max_cp, sizeof(uint32_t));
    uint32_t n_cp = pick_top_k_peaks(km, score, n, target_k, min_size, cps, max_cp);
    ri_kfree(km, score);
    float *events = (n_cp > 0) ? gen_events(km, sig, cps, n_cp, n, n_events)
                                : (*n_events = 0, (float*)NULL);
    ri_kfree(km, cps);
    return events;
}

/* ---- 15. Turbo (algorithmically identical to default, fused 1-pass) --- *
 * Same dual-window t-stat + gen_peaks logic as default, but with:
 *  - tstat computed INLINE (no temporary float array of length n × 2)
 *  - exact same float math (no sqrt skipping — preserves byte-identical
 *    boundary positions to default)
 *  - one fused loop over positions instead of three (tstat1, tstat2,
 *    gen_peaks) → saves 2/3 of memory traffic.
 *
 * The fused loop replicates default's `gen_peaks` state machine with
 * inline t-stat values; the state machine itself is bitwise identical. */
static float* detect_events_turbo(void *km, float *sig, uint32_t n, uint32_t *n_events,
                                   uint32_t w1, uint32_t w2,
                                   float th1, float th2, float ph) {
    if (n < 30) { *n_events = 0; return 0; }
    const float eta = FLT_MIN;

    float *ps = (float*)ri_kcalloc(km, n+1, sizeof(float));
    float *pss = (float*)ri_kcalloc(km, n+1, sizeof(float));
    for (uint32_t i = 0; i < n; i++) {
        ps[i+1] = ps[i] + sig[i];
        pss[i+1] = pss[i] + sig[i]*sig[i];
    }

    /* Replicate default's gen_peaks state for two detectors (short, long). */
    const uint32_t W[2] = {w1, w2};
    const float TH[2] = {th1, th2};
    int peak_pos[2]   = {-1, -1};
    float peak_val[2] = {FLT_MAX, FLT_MAX};
    int valid[2]      = {0, 0};
    uint32_t mask_to[2] = {0, 0};

    uint32_t *peaks = (uint32_t*)ri_kcalloc(km, n, sizeof(uint32_t));
    uint32_t n_peaks = 0;

    /* Bounds: tstat is well-defined for w_len <= i <= s_len - w_len. */
    uint32_t i_start = w2;
    uint32_t i_stop  = (n > w2) ? n - w2 : 0;
    for (uint32_t i = i_start; i <= i_stop; i++) {
        for (int k = 0; k < 2; k++) {
            uint32_t w = W[k];
            if (i < w || i + w > n) { /* tstat undefined here */ continue; }
            if (mask_to[k] >= i) continue;
            /* compute tstat[i] for window w (same formula as comp_tstat) */
            float sum1 = ps[i] - ((i > w) ? ps[i - w] : 0);
            float sumsq1 = pss[i] - ((i > w) ? pss[i - w] : 0);
            float sum2 = ps[i + w] - ps[i];
            float sumsq2 = pss[i + w] - pss[i];
            float m1 = sum1 / w, m2 = sum2 / w;
            float vcomb = (sumsq1/w - m1*m1 + sumsq2/w - m2*m2) / w;
            if (vcomb < eta) vcomb = eta;
            float dm = m2 - m1;
            float cur = fabsf(dm) / sqrtf(vcomb);

            /* Default's gen_peaks state machine — bitwise identical logic. */
            if (peak_pos[k] == -1) {
                if (cur < peak_val[k]) {
                    peak_val[k] = cur;
                } else if (cur - peak_val[k] > ph) {
                    peak_val[k] = cur;
                    peak_pos[k] = (int)i;
                }
            } else {
                if (cur > peak_val[k]) {
                    peak_val[k] = cur;
                    peak_pos[k] = (int)i;
                }
                if (peak_val[k] > TH[k]) {
                    for (int o = k+1; o < 2; o++) {
                        mask_to[o] = (uint32_t)peak_pos[k] + W[0];
                        peak_pos[o] = -1;
                        peak_val[o] = FLT_MAX;
                        valid[o] = 0;
                    }
                }
                if (peak_val[k] - cur > ph && peak_val[k] > TH[k]) {
                    valid[k] = 1;
                }
                if (valid[k] && ((int)i - peak_pos[k]) > (int)(W[k]/2)) {
                    peaks[n_peaks++] = (uint32_t)peak_pos[k];
                    peak_pos[k] = -1;
                    peak_val[k] = cur;
                    valid[k] = 0;
                }
            }
        }
    }
    ri_kfree(km, ps); ri_kfree(km, pss);
    float *events = 0;
    if (n_peaks > 0) events = gen_events(km, sig, peaks, n_peaks, n, n_events);
    else { *n_events = 0; }
    ri_kfree(km, peaks);
    return events;
}

/* Wrapper using "sensitive --r10" defaults. */
static float* detect_events_turbo_default(void *km, float *sig, uint32_t n, uint32_t *n_events) {
    /* Pull live defaults from --r10 sensitive: w1=3 w2=9 t1=6.5 t2=4 ph=0.2 */
    return detect_events_turbo(km, sig, n, n_events,
                                env_u("RH2_TURBO_W1", 3),
                                env_u("RH2_TURBO_W2", 9),
                                env_f("RH2_TURBO_T1", 6.5f),
                                env_f("RH2_TURBO_T2", 4.0f),
                                env_f("RH2_TURBO_PH", 0.2f));
}

static float* detect_events_python(void *km, float *sig, uint32_t n, const char *script_path, uint32_t *n_events) {
    if (!script_path || script_path[0] == '\0') {
        fprintf(stderr, "[ERROR] --segmenter python requires --segmenter-script <path>\n");
        *n_events = 0;
        return 0;
    }

    pthread_mutex_lock(&g_py_mtx);
    if (g_py_pid == 0) {
        if (spawn_python_worker(script_path) < 0) {
            pthread_mutex_unlock(&g_py_mtx);
            fprintf(stderr, "[ERROR] failed to spawn Python segmenter worker\n");
            *n_events = 0; return 0;
        }
    }

    /* send n + raw float32 signal */
    fprintf(g_py_in, "%u\n", n);
    fwrite(sig, sizeof(float), n, g_py_in);
    fflush(g_py_in);

    /* receive event count (text line) + binary float32 events */
    uint32_t n_ev = 0;
    if (fscanf(g_py_out, "%u", &n_ev) != 1) n_ev = 0;
    /* skip the trailing newline */
    int c; while ((c = fgetc(g_py_out)) != '\n' && c != EOF) {}
    float *events = 0;
    if (n_ev > 0) {
        events = (float*)ri_kmalloc(km, n_ev * sizeof(float));
        size_t got = fread(events, sizeof(float), n_ev, g_py_out);
        if (got != n_ev) {
            ri_kfree(km, events); events = 0; n_ev = 0;
        }
    }
    pthread_mutex_unlock(&g_py_mtx);
    *n_events = n_ev;
    return events;
}

float* normalize_signal(void *km,
						const float* sig,
						const uint32_t s_len,
						double* mean_sum,
						double* std_dev_sum,
						uint32_t* n_events_sum,
						uint32_t* n_sig)
{
	double sum = (*mean_sum), sum2 = (*std_dev_sum);
	double mean = 0, std_dev = 0;
	float* events = (float*)ri_kcalloc(km, s_len, sizeof(float));

	for (uint32_t i = 0; i < s_len; ++i) {
		sum += sig[i];
		sum2 += sig[i]*sig[i];
	}

	(*n_events_sum) += s_len;
	(*mean_sum) = sum;
	(*std_dev_sum) = sum2;

	mean = sum/(*n_events_sum);
	std_dev = sqrt(sum2/(*n_events_sum) - (mean)*(mean));

	float norm_val = 0;
	int k = 0;
	for(uint32_t i = 0; i < s_len; ++i){
		norm_val = (sig[i]-mean)/std_dev;
		if(norm_val < 3 && norm_val > -3) events[k++] = norm_val;
	}

	(*n_sig) = k;

	return events;
}

float* detect_events(void *km,
					 const uint32_t s_len,
					 const float* sig,
					 const uint32_t window_length1,
					 const uint32_t window_length2,
					 const float threshold1,
					 const float threshold2,
					 const float peak_height,
					 double* mean_sum,
					 double* std_dev_sum,
					 uint32_t* n_events_sum,
					 uint32_t* n_events,
					 uint32_t segmenter_type,
					 const char *python_script)
{
	float* prefix_sum = (float*)ri_kcalloc(km, s_len+1, sizeof(float));
	float* prefix_sum_square = (float*)ri_kcalloc(km, s_len+1, sizeof(float));

	uint32_t n_signals = 0;
	(*n_events) = 0;
	float* norm_signals = normalize_signal(km, sig, s_len, mean_sum, std_dev_sum, n_events_sum, &n_signals);
	if(n_signals == 0) { ri_kfree(km, prefix_sum); ri_kfree(km, prefix_sum_square); return 0; }

	// Dispatch to non-default segmenters
	if (segmenter_type != RI_SEGMENTER_DEFAULT) {
		float *events = 0;
		switch (segmenter_type) {
			case RI_SEGMENTER_HMM:
				events = detect_events_hmm(km, norm_signals, n_signals, n_events);
				break;
			case RI_SEGMENTER_PELT:
				events = detect_events_pelt(km, norm_signals, n_signals, n_events);
				break;
			case RI_SEGMENTER_BINSEG:
				events = detect_events_binseg(km, norm_signals, n_signals, n_events);
				break;
			case RI_SEGMENTER_SCRAPPIE:
				events = detect_events_scrappie(km, norm_signals, n_signals, n_events);
				break;
			case RI_SEGMENTER_PYTHON:
				events = detect_events_python(km, norm_signals, n_signals, python_script, n_events);
				break;
			case RI_SEGMENTER_CUSUM:
				events = detect_events_cusum(km, norm_signals, n_signals, n_events);
				break;
			case RI_SEGMENTER_GRADIENT:
				events = detect_events_gradient(km, norm_signals, n_signals, n_events);
				break;
			case RI_SEGMENTER_MAD:
				events = detect_events_mad(km, norm_signals, n_signals, n_events);
				break;
			case RI_SEGMENTER_BOCD:
				events = detect_events_bocd(km, norm_signals, n_signals, n_events);
				break;
			case RI_SEGMENTER_WINDOW:
				events = detect_events_window(km, norm_signals, n_signals, n_events);
				break;
			case RI_SEGMENTER_WAVELET:
				events = detect_events_wavelet(km, norm_signals, n_signals, n_events);
				break;
			case RI_SEGMENTER_TDA:
				events = detect_events_tda(km, norm_signals, n_signals, n_events);
				break;
			case RI_SEGMENTER_CONVBANK:
				events = detect_events_convbank(km, norm_signals, n_signals, n_events);
				break;
			case RI_SEGMENTER_FUSION:
				events = detect_events_fusion(km, norm_signals, n_signals, n_events);
				break;
			case RI_SEGMENTER_TURBO:
				events = detect_events_turbo_default(km, norm_signals, n_signals, n_events);
				break;
			case RI_SEGMENTER_PELT_CUDA:
				events = detect_events_pelt_cuda(km, norm_signals, n_signals, n_events);
				break;
			case RI_SEGMENTER_WILCOXON:
				events = detect_events_wilcoxon(km, norm_signals, n_signals, n_events);
				break;
		}
		if (events && *n_events > 0) {
			ri_kfree(km, norm_signals);
			ri_kfree(km, prefix_sum);
			ri_kfree(km, prefix_sum_square);
			return events;
		}
		// Fall through to default if segmenter returned nothing
	}

	// Default t-stat segmenter (original code)
	comp_prefix_prefixsq(norm_signals, n_signals, prefix_sum, prefix_sum_square);

	float* tstat1 = comp_tstat(km, prefix_sum, prefix_sum_square, n_signals, window_length1);
	float* tstat2 = comp_tstat(km, prefix_sum, prefix_sum_square, n_signals, window_length2);
	ri_detect_t short_detector = {.DEF_PEAK_POS = -1, .DEF_PEAK_VAL = FLT_MAX,
		.sig = tstat1, .s_len = n_signals, .threshold = threshold1,
		.window_length = window_length1, .masked_to = 0, .peak_pos = -1,
		.peak_value = FLT_MAX, .valid_peak = 0};
	ri_detect_t long_detector = {.DEF_PEAK_POS = -1, .DEF_PEAK_VAL = FLT_MAX,
		.sig = tstat2, .s_len = n_signals, .threshold = threshold2,
		.window_length = window_length2, .masked_to = 0, .peak_pos = -1,
		.peak_value = FLT_MAX, .valid_peak = 0};

	uint32_t* peaks = (uint32_t*)ri_kmalloc(km, n_signals * sizeof(uint32_t));
	ri_detect_t *detectors[2] = {&short_detector, &long_detector};
	uint32_t n_peaks = gen_peaks(detectors, 2, peak_height, prefix_sum, prefix_sum_square, peaks);
	ri_kfree(km, tstat1); ri_kfree(km, tstat2); ri_kfree(km, prefix_sum); ri_kfree(km, prefix_sum_square);

	float* events = 0;
	if(n_peaks > 0) events = gen_events(km, norm_signals, peaks, n_peaks, n_signals, n_events);
	if (getenv("RH2_DEBUG_EVENTS"))
		fprintf(stderr, "[DBG_DEF] n_sig=%u n_peaks=%u n_events=%u\n", n_signals, n_peaks, *n_events);
	ri_kfree(km, norm_signals); ri_kfree(km, peaks);
	return events;
}
