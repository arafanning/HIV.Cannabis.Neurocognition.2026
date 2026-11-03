# ==============================================================
# HIV + CANNABIS + NEUROCOGNITION: PRIMARY ANALYSIS RECONCILIATION
# Draft recovery of the post-October-2026 within-PWH specification
# ==============================================================
#
# PROVENANCE / STATUS
# - This is NEW, assembled from available repo scripts and later
#   documented analysis decisions. It is NOT an exact recovered
#   final original script and has not been executed on study data.
# - Retains original six repository .R files unchanged.
# - Reinstates within-PWH test-level z scoring as primary (N=113).
# - Treats None / Low / High past-year cannabis use categorically;
#   planned orthogonal contrasts characterize group patterns.
# - Uses HC3 robust inference and four hierarchical model steps.
# - Does not automatically write any subject-level data or output.
# - No participant-level data are included in this repository.
#
# Review RECOVERY_NOTES.md before using for a manuscript.
# ==============================================================

# 0. Dependencies and restricted local input --------------------
needed <- c("haven", "lavaan", "sandwich", "lmtest", "car", "emmeans")
missing_packages <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required packages manually: ",
       paste(missing_packages, collapse = ", "))
}

input_file <- "data/HIV+113.sav"
if (!file.exists(input_file)) {
  stop("Restricted input is missing: ", input_file,
       ". Place authorized local data in this path before running.")
}
data <- haven::read_sav(input_file)
if (nrow(data) != 113L) {
  warning("Expected 113 PWH cases, found ", nrow(data),
          ". Review inclusion criteria and analysis sample.")
}

needed_raw <- c(
  "bv_tr_rs", "hv_tt_rs", "hv_dr_rs", "bv_dr_rs", "mst_s_ct1",
  "gpt_nd_ttime", "wc_pe_rs", "igt_nt_rs", "wai_cd_rs", "wai_ss_rs",
  "lab_thel", "mhq_4_cd4_lowest", "lab_hiv", "du_mar4_12m_a",
  "phq_2_age", "phq_7_degree", "race_eth_binary"
)
absent <- setdiff(needed_raw, names(data))
if (length(absent)) {
  stop("Missing variables: ", paste(absent, collapse = ", "))
}

# 1. Neurocognition: within-sample test standardization --------
# Primary scores reflect the mean and SD of all 113 PWH cases,
# NOT the HIV-/cannabis- comparison group. Higher is better.
raw_tests <- c(
  "bv_tr_rs", "hv_tt_rs", "hv_dr_rs", "bv_dr_rs", "mst_s_ct1",
  "gpt_nd_ttime", "wc_pe_rs", "igt_nt_rs", "wai_cd_rs", "wai_ss_rs"
)
test_ids <- c(
  "bvtr", "hvtt", "hvdr", "bvdr", "mstct",
  "gpttim", "wcpe", "igtnt", "cd", "ss"
)
reverse_ids <- c("mstct", "gpttim", "wcpe")
test_z <- paste0(test_ids, "_z")
for (i in seq_along(raw_tests)) {
  x <- as.numeric(data[[raw_tests[i]]])
  if (test_ids[i] %in% reverse_ids) x <- -x
  if (!is.finite(stats::sd(x, na.rm = TRUE)) ||
      stats::sd(x, na.rm = TRUE) == 0) {
    stop("Cannot standardize test: ", raw_tests[i])
  }
  data[[test_z[i]]] <- as.numeric(scale(x))
}

# Faithful to the existing test-to-domain mappings.
# rowMeans(..., na.rm=TRUE) permits a one-test domain;
# a global score requires at least 9 of 10 test scores.
domain_map <- list(
  learning_z = c("bvtr_z", "hvtt_z"),
  memory_z = c("hvdr_z", "bvdr_z"),
  motor_z = c("mstct_z", "gpttim_z"),
  executive_z = c("wcpe_z", "igtnt_z"),
  processing_speed_z = c("cd_z", "ss_z")
)
for (nm in names(domain_map)) {
  data[[nm]] <- rowMeans(data[, domain_map[[nm]], drop = FALSE],
                         na.rm = TRUE)
  data[[nm]][is.nan(data[[nm]])] <- NA_real_
}
data$n_neuro_tests <- rowSums(!is.na(data[, test_z, drop = FALSE]))
data$global_z <- rowMeans(data[, test_z, drop = FALSE], na.rm = TRUE)
data$global_z[is.nan(data$global_z) | data$n_neuro_tests < 9] <- NA_real_

outcomes <- c(
  "global_z", "learning_z", "memory_z", "motor_z",
  "executive_z", "processing_speed_z"
)
print(summary(data[, outcomes, drop = FALSE]))

# 2. HIV disease severity --------------------------------------
# Preserve the original three-indicator lavaan construction.
# Negative square roots: lower CD4 values imply greater severity.
# The fixed 0.001 residual for current CD4 addresses the Heywood
# case encountered in the original model; it requires reporting.
data$cd4sqrt <- -sqrt(as.numeric(data$lab_thel))
data$nadirsqrt <- -sqrt(as.numeric(data$mhq_4_cd4_lowest))
data$log_vl2 <- log1p(as.numeric(data$lab_hiv))
severity_spec <- '
  DiseaseSev =~ nadirsqrt + cd4sqrt + log_vl2
  cd4sqrt ~~ 0.001*cd4sqrt
'
severity_fit <- lavaan::cfa(
  severity_spec, data = data, estimator = "MLR", missing = "fiml"
)
print(summary(severity_fit, fit.measures = TRUE, standardized = TRUE))
severity_scores <- as.numeric(lavaan::lavPredict(severity_fit))
if (length(severity_scores) != nrow(data)) {
  stop("CFA factor score length differs from data. Check missingness.")
}
data$DiseaseSev_c <- as.numeric(scale(
  severity_scores, center = TRUE, scale = FALSE
))

# 3. Past-year cannabis-use group ------------------------------
# Median among nonusers excluded from the median calculation.
grams <- as.numeric(data$du_mar4_12m_a)
if (any(grams < 0, na.rm = TRUE)) {
  stop("Past-year cannabis amounts include negatives; check coding.")
}
median_users <- stats::median(grams[is.finite(grams) & grams > 0])
data$cu_group <- factor(
  ifelse(is.na(grams), NA_character_,
  ifelse(grams == 0, "None",
  ifelse(grams <= median_users, "Low", "High"))),
  levels = c("None", "Low", "High")
)
group_counts <- table(data$cu_group, useNA = "ifany")
print(group_counts)
if (!identical(unname(as.integer(table(data$cu_group))),
               c(62L, 27L, 24L))) {
  warning("Group counts differ from documented 62/27/24. Verify source.")
}
medians <- vapply(levels(data$cu_group), function(g) {
  stats::median(grams[data$cu_group == g], na.rm = TRUE)
}, numeric(1))
print(medians)  # Earlier documented: None=0, Low=36, High=720 g

# Categorical inference does NOT assume equal intervals.
# Orthogonal L/Q weights are secondary descriptions of the three
# adjusted group means, not a fitted continuous dose-response.
C_equal <- stats::contr.poly(3)
rownames(C_equal) <- levels(data$cu_group)
C_observed <- stats::contr.poly(3, scores = log1p(medians))
rownames(C_observed) <- levels(data$cu_group)

# 4. Demographic terms -----------------------------------------
data$phq_2_age_c <- as.numeric(scale(as.numeric(data$phq_2_age),
                                   center = TRUE, scale = FALSE))
# Existing analysis treated phq_7_degree as a numeric variable;
# confirm that the coding represents the intended education unit.
data$phq_7_degree_c <- as.numeric(scale(
  as.numeric(data$phq_7_degree), center = TRUE, scale = FALSE
))
data$race_eth_binary_covfac <- factor(data$race_eth_binary)
print(table(data$race_eth_binary_covfac, data$cu_group,
            useNA = "ifany"))

# 5. HC3 inference helper --------------------------------------
hc3 <- function(model) sandwich::vcovHC(model, type = "HC3")

# Joint Wald F-test of selected regression coefficients with
# an explicit contrast matrix to avoid formula-term ambiguity.
joint_hc3 <- function(model, terms) {
  coef_names <- names(stats::coef(model))
  if (!length(terms) || any(!terms %in% coef_names)) {
    stop("Cannot test missing terms: ",
         paste(setdiff(terms, coef_names), collapse = ", "))
  }
  H <- matrix(0, nrow = length(terms), ncol = length(coef_names),
              dimnames = list(terms, coef_names))
  for (i in seq_along(terms)) H[i, terms[i]] <- 1
  car::linearHypothesis(model, hypothesis.matrix = H,
                        vcov. = hc3(model), test = "F")
}

# 6. Four-step hierarchy, same complete-case sample per outcome
# Step 1: severity
# Step 2: + 3-level cannabis group
# Step 3: + age, education, race/ethnicity
# Step 4: + cannabis x HIV disease severity
#
# All stages use identical cases *within each outcome*, so
# nested changes in R^2 remain interpretable. OLS fit is used
# for R^2; statistical inference uses HC3 covariance.
fit_outcome <- function(y) {
  vars <- c(y, "DiseaseSev_c", "cu_group",
            "phq_2_age_c", "phq_7_degree_c",
            "race_eth_binary_covfac")
  d <- as.data.frame(data[stats::complete.cases(data[, vars]), vars])
  d$cu_group <- droplevels(d$cu_group)
  d$race_eth_binary_covfac <- droplevels(d$race_eth_binary_covfac)
  if (nlevels(d$cu_group) != 3L) {
    stop("Insufficient cannabis groups in complete-case data for ", y)
  }
  if (nlevels(d$race_eth_binary_covfac) < 2L) {
    stop("Insufficient race/ethnicity levels for ", y)
  }

  m1 <- stats::lm(stats::as.formula(paste(y, "~ DiseaseSev_c")), data = d)
  m2 <- stats::lm(stats::as.formula(paste(y, "~ DiseaseSev_c + cu_group")),
                  data = d)
  m3 <- stats::lm(stats::as.formula(paste(
    y, "~ DiseaseSev_c + cu_group + phq_2_age_c +",
    "phq_7_degree_c + race_eth_binary_covfac"
  )), data = d)
  m4 <- stats::lm(stats::as.formula(paste(
    y, "~ DiseaseSev_c * cu_group + phq_2_age_c +",
    "phq_7_degree_c + race_eth_binary_covfac"
  )), data = d)

  beta_cu_m2 <- grep("^cu_group", names(stats::coef(m2)), value = TRUE)
  beta_cu_m3 <- grep("^cu_group", names(stats::coef(m3)), value = TRUE)
  beta_int <- grep(":", names(stats::coef(m4)), value = TRUE)

  # HC3-robust 2 df block tests (not the usual OLS nested ANOVA)
  test_cu_m2 <- joint_hc3(m2, beta_cu_m2)
  test_cu_m3 <- joint_hc3(m3, beta_cu_m3)
  test_mod <- joint_hc3(m4, beta_int)

  # Adjusted group means, robust contrasts and Tukey comparisons.
  emm <- emmeans::emmeans(m3, ~ cu_group, vcov. = hc3(m3))
  pairwise <- pairs(emm, adjust = "tukey")
  planned_equal <- emmeans::contrast(
    emm, method = list(L = C_equal[, ".L"],
                       Q = C_equal[, ".Q"]), adjust = "none"
  )
  planned_observed <- emmeans::contrast(
    emm, method = list(L_observed = C_observed[, ".L"],
                       Q_observed = C_observed[, ".Q"]), adjust = "none"
  )

  # Race-removed sensitivity is performed on the same cases.
  m3_no_race <- stats::lm(stats::as.formula(paste(
    y, "~ DiseaseSev_c + cu_group + phq_2_age_c + phq_7_degree_c"
  )), data = d)

  r2 <- vapply(list(m1, m2, m3, m4), function(m) summary(m)$r.squared,
               numeric(1))
  list(
    outcome = y, n = nrow(d),
    models = list(severity = m1, cannabis = m2,
                  adjusted = m3, moderation = m4),
    r2 = r2,
    delta_r2 = diff(r2),
    hc3_coefficients = lapply(list(m1, m2, m3, m4), function(m) {
      lmtest::coeftest(m, vcov. = hc3(m))
    }),
    hc3_group_at_step2 = test_cu_m2,
    hc3_group_at_step3 = test_cu_m3,
    hc3_interaction_at_step4 = test_mod,
    adjusted_means = emm,
    pairwise_tukey = pairwise,
    contrasts_equal_spacing = planned_equal,
    contrasts_observed_spacing = planned_observed,
    no_race_sensitivity = list(
      model = m3_no_race,
      hc3_group = joint_hc3(
        m3_no_race,
        grep("^cu_group", names(stats::coef(m3_no_race)), value = TRUE)
      )
    )
  )
}

analysis_results <- stats::setNames(lapply(outcomes, fit_outcome), outcomes)
for (y in outcomes) {
  z <- analysis_results[[y]]
  cat("\n\n=== ", y, " (n = ", z$n, ") ===\n", sep = "")
  print(z$r2)
  print(z$delta_r2)
  print(z$hc3_group_at_step3)
  print(z$hc3_interaction_at_step4)
  print(summary(z$pairwise_tukey))
  print(summary(z$contrasts_equal_spacing))
  print(summary(z$contrasts_observed_spacing))
}

# REVIEW BEFORE REPORTING
# - This reconstruction has NOT been executed on the restricted data.
# - Confirm model construction, sample-specific coding, and CFA.
# - Confirm whether emmeans Tukey comparisons should be shown with
#   robust vcov and any across-domain multiplicity handling.
# - User-only continuous cannabis amount, alcohol covariate, and
#   other exploratory analyses should be recovered separately,
#   not silently added to these primary models.
