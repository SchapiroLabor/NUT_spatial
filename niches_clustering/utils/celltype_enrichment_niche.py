### FUNCTION TO CALCULATE WHICH CELL TYPES ARE ENRICHED IN A NICHE ###

import numpy as np
import pandas as pd
from tqdm import tqdm

def _proportion(obs: pd.DataFrame, id_key: str, val_key: str, normalize: bool = True) -> pd.DataFrame:
    """
    Returns a contingency table:
      rows = id_key
      cols = val_key
    If normalize=True, each row sums to 1 (row-wise proportions).
    """
    df = (
        obs[[id_key, val_key]]
        .value_counts()
        .reset_index(name="n")
        .pivot(index=id_key, columns=val_key, values="n")
        .fillna(0.0)
    )
    if normalize:
        row_sums = df.sum(axis=1).replace(0.0, np.nan)
        return df.div(row_sums, axis=0).fillna(0.0)
    return df

def _enrichment(observed: pd.DataFrame, expected: pd.DataFrame, log: bool = True, eps: float = 1e-12) -> pd.DataFrame:
    """
    Enrichment = observed / expected, optionally log2-transformed.
    eps prevents division by 0 and log2(0).
    """
    ratio = (observed + eps).div(expected + eps)
    if log:
        ratio = np.log2(ratio)
    return ratio

def _observed_permuted(obs_df: pd.DataFrame, group_key: str, label_key: str, rng: np.random.Generator) -> pd.DataFrame:
    """
    Permute group labels across cells (break association with neighborhoods),
    then compute proportions (rows=groups, cols=labels).
    """
    tmp = obs_df[[group_key, label_key]].copy()
    tmp[group_key] = rng.permutation(tmp[group_key].to_numpy())
    # proportions by label, then transpose to (groups x labels)
    return _proportion(tmp, id_key=label_key, val_key=group_key, normalize=True).T

def _empirical_pvalues_from_perm_enr(obs_enr: pd.DataFrame, perm_enr_stack: np.ndarray) -> pd.DataFrame:
    """
    One-sided empirical p-values with +1 correction (never exactly 0):
      - if obs_enr > 0: p = (#{perm >= obs} + 1) / (n_perms + 1)
      - if obs_enr < 0: p = (#{perm <= obs} + 1) / (n_perms + 1)
      - if obs_enr == 0: p = 1
    """
    obs = obs_enr.to_numpy()
    n_perms = perm_enr_stack.shape[0]
    p = np.ones_like(obs, dtype=float)

    pos = obs > 0
    neg = obs < 0

    if np.any(pos):
        k = (perm_enr_stack[:, pos] >= obs[pos]).sum(axis=0)
        p[pos] = (k + 1) / (n_perms + 1)

    if np.any(neg):
        k = (perm_enr_stack[:, neg] <= obs[neg]).sum(axis=0)
        p[neg] = (k + 1) / (n_perms + 1)

    return pd.DataFrame(p, index=obs_enr.index, columns=obs_enr.columns)

def nhood_enrichment(
    adata,
    group_key: str = "cell_category",
    label_key: str = "neigh_kmeans_count",
    log: bool = True,
    pvalues: bool = False,
    n_perms: int = 200,
    random_state: int = 0,
    eps: float = 1e-12,
    copy: bool = False,
):
    """
    Computes neighborhood enrichment of `group_key` across `label_key`.

    Outputs into:
      adata.uns[f"{label_key}_{group_key}_enrichment"] = {
        "enrichment": (groups x labels) DataFrame,
        "observed":   (groups x labels) observed proportions,
        "expected":   (groups x labels) expected proportions,
        "pvalue":     (optional) empirical p-values
      }
    """
    # Observed proportions: rows=groups, cols=labels
    observed = _proportion(adata.obs, id_key=label_key, val_key=group_key, normalize=True).T
    observed = observed.fillna(0.0)

    labels = observed.columns
    groups = observed.index

    if not pvalues:
        # Expected = global group frequency, repeated for each label
        expected_vec = adata.obs[group_key].value_counts(normalize=True)
        expected = pd.DataFrame({lab: expected_vec.reindex(groups).fillna(0.0) for lab in labels})
        expected = expected.loc[groups, labels]

        enr = _enrichment(observed, expected, log=log, eps=eps)
        result = {"enrichment": enr, "observed": observed, "expected": expected}

    else:
        rng = np.random.default_rng(random_state)
        obs_df = adata.obs[[group_key, label_key]].copy()

        # Build null distribution of ENRICHMENT (same statistic as observed)
        perm_enr_list = []
        perm_obs_list = []

        for _ in tqdm(range(n_perms)):
            perm_obs = _observed_permuted(obs_df, group_key, label_key, rng=rng)
            perm_obs = perm_obs.reindex(index=groups, columns=labels).fillna(0.0)
            perm_obs_list.append(perm_obs.to_numpy())

            # For the permuted table, expected under shuffle = mean across labels (global in that perm)
            perm_expected_vec = perm_obs.mean(axis=1)  # (groups,)
            perm_expected = pd.DataFrame({lab: perm_expected_vec for lab in labels}).loc[groups, labels]

            perm_enr = _enrichment(perm_obs, perm_expected, log=log, eps=eps)
            perm_enr_list.append(perm_enr.to_numpy())

        perm_obs_stack = np.stack(perm_obs_list, axis=0)   # (n_perms, n_groups, n_labels)
        perm_enr_stack = np.stack(perm_enr_list, axis=0)   # (n_perms, n_groups, n_labels)

        # Expected (for reporting) = mean permuted proportions
        expected = pd.DataFrame(perm_obs_stack.mean(axis=0), index=groups, columns=labels)

        # Observed enrichment vs expected
        enr = _enrichment(observed, expected, log=log, eps=eps)

        # Empirical p-values by comparing obs enrichment to perm enrichment distribution
        pval = _empirical_pvalues_from_perm_enr(enr, perm_enr_stack)

        result = {"enrichment": enr, "pvalue": pval, "observed": observed, "expected": expected}

    result["params"] = dict(
        group_key=group_key,
        label_key=label_key,
        log=log,
        pvalues=pvalues,
        n_perms=n_perms,
        random_state=random_state,
        eps=eps,
    )

    if copy:
        return result

    adata.uns[f"{label_key}_{group_key}_enrichment"] = result
    return None