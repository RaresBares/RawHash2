#include "revent.h"
#include "kalloc.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>
#include <stdio.h>
#include "rutils.h"

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

	for (uint32_t pi = 0, i = 0; pi < peak_size && i < n_ev; pi++){
		if (!(peaks[pi] > 0 && peaks[pi] < s_len)) continue;

    	segment_length = peaks[pi] - start_idx;
		if (segment_length < 500) // Skip if the segment is too long
			events[i++] = calculate_mean_of_filtered_segment(sig + start_idx, segment_length);
		start_idx = peaks[pi];
	}

	(*n_events) = n_ev;
	return events;
}

static float* detect_events_hmm(void *km, float *sig, uint32_t n, uint32_t *n_events) {
    if (n < 10) { *n_events = 0; return 0; }
    uint32_t n_states = 4;
    float state_means[4], state_stds[4];
    // Initialize with quantiles
    float *sorted = (float*)ri_kmalloc(km, n * sizeof(float));
    memcpy(sorted, sig, n * sizeof(float));
    qsort(sorted, n, sizeof(float), compare_floats);
    for (uint32_t s = 0; s < n_states; s++) {
        uint32_t idx = (s + 1) * n / (n_states + 1);
        state_means[s] = sorted[idx];
        state_stds[s] = 1.0f;
    }
    ri_kfree(km, sorted);

    // K-means refinement (5 iterations)
    uint32_t *assignments = (uint32_t*)ri_kcalloc(km, n, sizeof(uint32_t));
    for (int iter = 0; iter < 5; iter++) {
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
                if (state_stds[s] < 0.1f) state_stds[s] = 0.1f;
            }
        }
    }

    // Viterbi
    float stay_log = logf(0.95f);
    float trans_log = logf(0.05f / (n_states - 1));
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

    ri_kfree(km, V); ri_kfree(km, path);

    // Extract breakpoints
    uint32_t *peaks = (uint32_t*)ri_kmalloc(km, n * sizeof(uint32_t));
    uint32_t n_peaks = 0;
    for (uint32_t i = 1; i < n; i++)
        if (states[i] != states[i-1]) peaks[n_peaks++] = i;

    float *events = 0;
    if (n_peaks > 0) events = gen_events(km, sig, peaks, n_peaks, n, n_events);
    else { *n_events = 0; }
    ri_kfree(km, peaks); ri_kfree(km, assignments);
    return events;
}

static float* detect_events_pelt(void *km, float *sig, uint32_t n, uint32_t *n_events) {
    if (n < 10) { *n_events = 0; return 0; }
    float penalty = 2.0f * logf((float)n);  // BIC penalty
    uint32_t min_size = 5;

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
            if (j - t < min_size) continue;
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
    ri_kfree(km, cands); ri_kfree(km, cps);
    return events;
}

static void binseg_recursive(float *sig, uint32_t start, uint32_t end,
                              double *ps, double *pss, float penalty, uint32_t min_size,
                              uint32_t *cps, uint32_t *n_cp, uint32_t max_cp) {
    if (end - start < 2 * min_size || *n_cp >= max_cp) return;

    double total_sum = ps[end] - ps[start];
    double total_sumsq = pss[end] - pss[start];
    uint32_t total_len = end - start;
    double total_cost = total_sumsq - total_sum * total_sum / total_len;

    double best_gain = -1;
    uint32_t best_t = start + min_size;

    for (uint32_t t = start + min_size; t <= end - min_size; t++) {
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

    if (best_gain > penalty) {
        cps[(*n_cp)++] = best_t;
        binseg_recursive(sig, start, best_t, ps, pss, penalty, min_size, cps, n_cp, max_cp);
        binseg_recursive(sig, best_t, end, ps, pss, penalty, min_size, cps, n_cp, max_cp);
    }
}

static float* detect_events_binseg(void *km, float *sig, uint32_t n, uint32_t *n_events) {
    if (n < 10) { *n_events = 0; return 0; }
    float penalty = logf((float)n); // BIC with k=1 (mean only, normalized signals)
    uint32_t min_size = 5;
    uint32_t max_cp = n / min_size;

    double *ps = (double*)ri_kcalloc(km, n+1, sizeof(double));
    double *pss = (double*)ri_kcalloc(km, n+1, sizeof(double));
    for (uint32_t i = 0; i < n; i++) {
        ps[i+1] = ps[i] + sig[i];
        pss[i+1] = pss[i] + (double)sig[i] * sig[i];
    }

    uint32_t *cps = (uint32_t*)ri_kcalloc(km, max_cp, sizeof(uint32_t));
    uint32_t n_cp = 0;
    binseg_recursive(sig, 0, n, ps, pss, penalty, min_size, cps, &n_cp, max_cp);

    // Sort changepoints
    for (uint32_t i = 0; i < n_cp; i++)
        for (uint32_t j = i+1; j < n_cp; j++)
            if (cps[i] > cps[j]) { uint32_t tmp = cps[i]; cps[i] = cps[j]; cps[j] = tmp; }

    ri_kfree(km, ps); ri_kfree(km, pss);

    float *events = 0;
    if (n_cp > 0) events = gen_events(km, sig, cps, n_cp, n, n_events);
    else { *n_events = 0; }
    ri_kfree(km, cps);
    return events;
}

static float* detect_events_window(void *km, float *sig, uint32_t n, uint32_t *n_events) {
    if (n < 10) { *n_events = 0; return 0; }
    uint32_t w = 20; // window size
    float penalty = logf((float)(2 * w)); // local window penalty, not global
    uint32_t min_size = 5;

    double *ps = (double*)ri_kcalloc(km, n+1, sizeof(double));
    double *pss = (double*)ri_kcalloc(km, n+1, sizeof(double));
    for (uint32_t i = 0; i < n; i++) {
        ps[i+1] = ps[i] + sig[i];
        pss[i+1] = pss[i] + (double)sig[i] * sig[i];
    }

    // Compute cost reduction at each point using window
    uint32_t *peaks = (uint32_t*)ri_kmalloc(km, n * sizeof(uint32_t));
    uint32_t n_peaks = 0;

    for (uint32_t i = w; i <= n - w; i += min_size) {
        uint32_t left_s = (i > w) ? i - w : 0;
        uint32_t right_e = (i + w < n) ? i + w : n;

        double total_sum = ps[right_e] - ps[left_s];
        double total_sumsq = pss[right_e] - pss[left_s];
        uint32_t total_len = right_e - left_s;
        double total_cost = total_sumsq - total_sum * total_sum / total_len;

        double left_sum = ps[i] - ps[left_s];
        double left_sumsq = pss[i] - pss[left_s];
        uint32_t left_len = i - left_s;
        double left_cost = left_sumsq - left_sum * left_sum / left_len;

        double right_sum = ps[right_e] - ps[i];
        double right_sumsq = pss[right_e] - pss[i];
        uint32_t right_len = right_e - i;
        double right_cost = right_sumsq - right_sum * right_sum / right_len;

        double gain = total_cost - left_cost - right_cost;
        if (gain > penalty) peaks[n_peaks++] = i;
    }

    ri_kfree(km, ps); ri_kfree(km, pss);

    float *events = 0;
    if (n_peaks > 0) events = gen_events(km, sig, peaks, n_peaks, n, n_events);
    else { *n_events = 0; }
    ri_kfree(km, peaks);
    return events;
}

static float* detect_events_python(void *km, float *sig, uint32_t n, const char *script_path, uint32_t *n_events) {
    if (!script_path || script_path[0] == '\0') {
        fprintf(stderr, "[ERROR] --segmenter python requires --segmenter-script <path>\n");
        *n_events = 0;
        return 0;
    }

    int pipe_in[2], pipe_out[2];
    if (pipe(pipe_in) < 0 || pipe(pipe_out) < 0) {
        fprintf(stderr, "[ERROR] pipe() failed for Python segmenter\n");
        *n_events = 0;
        return 0;
    }

    pid_t pid = fork();
    if (pid == 0) {
        close(pipe_in[1]); close(pipe_out[0]);
        dup2(pipe_in[0], STDIN_FILENO);
        dup2(pipe_out[1], STDOUT_FILENO);
        close(pipe_in[0]); close(pipe_out[1]);
        execlp("python3", "python3", script_path, (char*)NULL);
        _exit(1);
    }

    close(pipe_in[0]); close(pipe_out[1]);

    FILE *to_child = fdopen(pipe_in[1], "w");
    fprintf(to_child, "%u\n", n);
    for (uint32_t i = 0; i < n; i++) fprintf(to_child, "%.8f\n", sig[i]);
    fclose(to_child);

    FILE *from_child = fdopen(pipe_out[0], "r");
    uint32_t n_ev = 0;
    if (fscanf(from_child, "%u", &n_ev) != 1) n_ev = 0;

    float *events = 0;
    if (n_ev > 0) {
        events = (float*)ri_kmalloc(km, n_ev * sizeof(float));
        for (uint32_t i = 0; i < n_ev; i++) {
            if (fscanf(from_child, "%f", &events[i]) != 1) { events[i] = 0; }
        }
    }
    fclose(from_child);

    int status;
    waitpid(pid, &status, 0);
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
			case RI_SEGMENTER_WINDOW:
				events = detect_events_window(km, norm_signals, n_signals, n_events);
				break;
			case RI_SEGMENTER_PYTHON:
				events = detect_events_python(km, norm_signals, n_signals, python_script, n_events);
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
	ri_kfree(km, norm_signals); ri_kfree(km, peaks);
	return events;
}
