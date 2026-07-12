library(ggplot2)
library(patchwork)
library(scales)
library(tidyverse)
library(ggpubr)
library(ggbeeswarm)
library(ComplexHeatmap)
library(broom)
library(circlize)
library(survminer)
library(survival)
library(showtext)
library(purrr)


# ---- builders ----
pct_patient <- function(df, celltype, niche = NULL) {
  d <- df; if (!is.null(niche)) d <- d %>% filter(neighborhood == niche)
  d %>%
    group_by(patient_ID, patient_exp) %>%
    summarise(n_total = n(),
              value = 100 * sum(cell_category == celltype) / n(), .groups = "drop") %>%
    filter(n_total >= MIN_CELLS) %>%
    group_by(patient_ID) %>% summarise(value = mean(value), .groups = "drop")
}
niche_pct_patient <- function(df, niche) {
  df %>%
    group_by(patient_ID, patient_exp) %>%
    summarise(n_total = n(),
              value = 100 * sum(neighborhood == niche) / n(), .groups = "drop") %>%
    filter(n_total >= MIN_CELLS) %>%
    group_by(patient_ID) %>% summarise(value = mean(value), .groups = "drop")
}

# ---- univariable Cox engine ----
make_km <- function(pat_df, label, split = "median", thresh = 1,
                    xlim = NULL, break.by = NULL) {
  d <- pat_df %>% inner_join(surv_patient, by = "patient_ID") %>%
    filter(!is.na(value), !is.na(survival_time_months), !is.na(event))

  if (split == "median") {
    cut <- median(d$value, na.rm = TRUE); lv <- c("Low", "High")
    d <- d %>% mutate(grp = ifelse(value >= cut, "High", "Low"))
  } else if (split == "presence") {
    d <- d %>% mutate(grp = ifelse(value >= thresh, "Present", "Absent"))
    lv <- c("Absent", "Present")
  }
  d$grp <- factor(d$grp, levels = lv)
  tab <- table(d$grp)

  if (length(tab) < 2 || any(tab < 2)) {
    return(list(plot = NULL, summary = tibble(
      variable = label, split = split, HR = NA, conf.low = NA, conf.high = NA,
      p.value = NA, ph_p = NA, n = nrow(d), events = sum(d$event),
      n_g1 = tab[[1]], n_g2 = ifelse(length(tab) > 1, tab[[2]], 0),
      note = "degenerate split -> try presence")))
  }

  fit_cox <- coxph(Surv(survival_time_months, event) ~ grp, data = d)
  s  <- summary(fit_cox); cn <- rownames(s$coefficients)[1]
  cox_p  <- s$coefficients[cn, "Pr(>|z|)"]
  cox_hr <- s$coefficients[cn, "exp(coef)"]
  ci_lo  <- s$conf.int[cn, "lower .95"]
  ci_hi  <- s$conf.int[cn, "upper .95"]
  ph_p   <- tryCatch(cox.zph(fit_cox)$table[1, "p"], error = function(e) NA_real_)

  cat(sprintf("%-42s HR=%.2f [%.2f-%.2f]  p=%.3f  PH_p=%.3f  (%s)\n",
              label, cox_hr, ci_lo, ci_hi, cox_p, ph_p,
              paste(names(tab), tab, sep="=", collapse=", ")))

  fit_km <- survfit(Surv(survival_time_months, event) ~ grp, data = d)
  xmax <- if (!is.null(xlim)) xlim[2] else max(d$survival_time_months)
  km <- ggsurvplot(fit_km, data = d,
    pval = paste0("Cox p = ", signif(cox_p, 2), "\nHR = ", signif(cox_hr, 2)),
    pval.coord = c(xmax * 0.45, 0.85), conf.int = TRUE,
    risk.table = TRUE, risk.table.height = 0.28,
    palette = c("#2E86AB", "#E74C3C"),
    xlab = "Time (months)", ylab = "Survival probability", title = label,
    legend.title = "", legend.labs = paste0(lv, " (n=", c(tab[[1]], tab[[2]]), ")"),
    xlim = xlim, break.time.by = break.by,
    ggtheme = theme_classic(base_size = 12))

  list(plot = km, summary = tibble(
    variable = label, split = split, HR = cox_hr,
    conf.low = ci_lo, conf.high = ci_hi, p.value = cox_p, ph_p = ph_p,
    n = nrow(d), events = sum(d$event), n_g1 = tab[[1]], n_g2 = tab[[2]], note = ""))
}

`%||%` <- function(a, b) if (is.null(a)) b else a
