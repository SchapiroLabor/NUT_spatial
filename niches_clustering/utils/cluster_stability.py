###DETERMINE THE OPTIMAL NUMBER OF CLUSTERS FOR THE NICHE IDENTIFICATION METHOD###
import numpy as np
from sklearn.cluster import KMeans 
from sklearn.metrics import adjusted_rand_score
import pandas as pd


def cluster_stability_calculate(
    adata,
    uns_key="spatial_prop_delaunay",
    k_min=5,
    k_max=15,
    n_repeats=5,
    sample_size=30000,
    random_state=0
):
    X = adata.uns[uns_key].to_numpy()
    rng = np.random.default_rng(random_state)

    results = []

    for k in range(k_min, k_max + 1):

        labelings = []

        for r in range(n_repeats):
            idx = rng.choice(X.shape[0], size=min(sample_size, X.shape[0]), replace=False)
            X_sub = X[idx]

            km = KMeans(n_clusters=k, n_init=20, random_state=r)
            labels = km.fit_predict(X_sub)

            labelings.append((idx, labels))

        # Compute pairwise ARI between repeats
        ari_scores = []

        for i in range(len(labelings)):
            for j in range(i + 1, len(labelings)):

                idx_i, lab_i = labelings[i]
                idx_j, lab_j = labelings[j]

                # Find overlapping cells between two subsamples
                common, i_pos, j_pos = np.intersect1d(idx_i, idx_j, return_indices=True)

                if len(common) > 0:
                    ari = adjusted_rand_score(lab_i[i_pos], lab_j[j_pos])
                    ari_scores.append(ari)

        mean_ari = np.mean(ari_scores)
        results.append((k, mean_ari))

        print(f"k={k} | mean ARI={mean_ari:.4f}")

    results_df = pd.DataFrame(results, columns=["k", "mean_ARI"]).set_index("k")

    best_k = results_df["mean_ARI"].idxmax()
    best_score = results_df["mean_ARI"].max()

    return best_k, best_score, results_df