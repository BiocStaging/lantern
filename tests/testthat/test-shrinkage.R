# ============================================================================
# Shrinkage of ambiguous splits in split_diploid()/split_diploid_multi()
# (`shrink = TRUE`, always on in ancestry_split()): the w-weighted blend of
# the raw ratio toward 1/2
# ============================================================================

# ============================================================================
# split_diploid() (2-population) shrinkage
# ============================================================================

test_that("split_diploid with shrink = FALSE reproduces the pre-shrinkage estimator exactly", {
  gt <- matrix(c(1, 1, 1), nrow = 1, ncol = 3)
  an <- matrix(c(3, 2, 1), nrow = 1, ncol = 3)
  result <- split_diploid(gt, an)   # shrink = FALSE
  expect_equal(result$african[1, 2], 0.5)
  expect_equal(result$european[1, 2], 0.5)
})

test_that("split_diploid singleton case shrinks fully to 1/2 (w = 1)", {
  gt <- matrix(1, nrow = 1, ncol = 1)
  an <- matrix(2, nrow = 1, ncol = 1)   # single mixed het, no other evidence
  result <- split_diploid(gt, an, shrink = TRUE)
  expect_equal(result$african[1, 1], 0.5, tolerance = 1e-8)
  expect_equal(result$european[1, 1], 0.5, tolerance = 1e-8)
})

test_that("split_diploid without ambiguous carriers is unaffected by shrinkage (w = 0)", {
  # No mixed hets (N5=0) at all -> w=0, raw ratio unchanged
  gt <- matrix(c(1, 1, 1), nrow = 1, ncol = 3)
  an <- matrix(c(3, 3, 1), nrow = 1, ncol = 3)   # 2 pure-AFR het, 1 pure-EUR het
  result_shrunk   <- split_diploid(gt, an, shrink = TRUE)
  result_unshrunk <- split_diploid(gt, an)
  expect_equal(result_shrunk$african, result_unshrunk$african)
  expect_equal(result_shrunk$european, result_unshrunk$european)
})

test_that("split_diploid partial-w blends raw ratio and 1/2 proportionally", {
  # variant: 1 pure-AFR het (N2=1), 1 mixed het (N5=1); no hom-alt carriers
  # here, so allele-weighted total_alt happens to equal the unweighted
  # carrier count (2). w = N5/total_alt = 1/2
  # raw: denom = total_alt - N5 = (1+1) - 1 = 1 ; p1_raw = (0+1+0)/1 = 1, p2_raw = 0
  gt <- matrix(c(1, 1), nrow = 1, ncol = 2)
  an <- matrix(c(3, 2), nrow = 1, ncol = 2)
  result <- split_diploid(gt, an, shrink = TRUE)
  w <- 0.5
  expect_equal(result$african[1, 2], (1 - w) * 1.0 + w * 0.5, tolerance = 1e-8)
  expect_equal(result$european[1, 2], (1 - w) * 0.0 + w * 0.5, tolerance = 1e-8)
})

test_that("split_diploid partial-w gives hom-alt carriers 2x the confidence weight of het carriers", {
  # variant: 1 pure-AFR HOM-ALT (N1=1, worth 2 alt alleles of confidence),
  # 1 mixed het (N5=1, ambiguous, worth 1 allele of confidence)
  # total_alt = 2*N1 + N5 = 2 + 1 = 3 ; w = N5/total_alt = 1/3
  # (a naive unweighted w = N5/carriers would give 1/2, which is wrong: the
  # hom-alt carrier's evidence should outweigh the ambiguous het's 2-to-1)
  gt <- matrix(c(2, 1), nrow = 1, ncol = 2)
  an <- matrix(c(3, 2), nrow = 1, ncol = 2)
  result <- split_diploid(gt, an, shrink = TRUE)
  w <- 1 / 3
  expect_equal(result$african[1, 2], (1 - w) * 1.0 + w * 0.5, tolerance = 1e-8)
  expect_equal(result$european[1, 2], (1 - w) * 0.0 + w * 0.5, tolerance = 1e-8)
})

# ============================================================================
# split_diploid_multi() (K-population) shrinkage
# ============================================================================

test_that("split_diploid_multi with shrink = FALSE reproduces original behavior (K=3)", {
  gt  <- matrix(c(1, 1, 1, 1), nrow = 1, ncol = 4)
  anc <- matrix(c(1, 2, 4, 3), nrow = 1, ncol = 4)  # AFR pure, AFR/EUR mixed, NAT pure, AFR/NAT mixed
  pure  <- c(AFR = 1L, EUR = 2L, NAT = 3L)
  mixed <- data.frame(code = c(2L, 4L, 5L),
                      pop1 = c("AFR", "AFR", "EUR"),
                      pop2 = c("EUR", "NAT", "NAT"))
  out <- split_diploid_multi(gt, anc, pure, mixed)
  expect_named(out, c("AFR", "EUR", "NAT"))
})

test_that("split_diploid_multi singleton mixed pair shrinks fully to 1/2 (w = 1)", {
  # Single AFR/EUR mixed het, nothing else -> no evidence for this pair,
  # w = 1, split exactly 1/2 within the pair; NAT gets nothing
  gt  <- matrix(1, nrow = 1, ncol = 1)
  anc <- matrix(4, nrow = 1, ncol = 1)   # mixed code 4 = AFR/EUR pair below
  pure  <- c(AFR = 1L, EUR = 2L, NAT = 3L)
  mixed <- data.frame(code = c(4L, 5L, 6L),
                      pop1 = c("AFR", "AFR", "EUR"),
                      pop2 = c("EUR", "NAT", "NAT"))
  out <- split_diploid_multi(gt, anc, pure, mixed, shrink = TRUE)
  expect_equal(out$AFR[1, 1], 0.5, tolerance = 1e-8)
  expect_equal(out$EUR[1, 1], 0.5, tolerance = 1e-8)
  expect_equal(out$NAT[1, 1], 0.0)
})

# ============================================================================
# End-to-end: ancestry_split() shrinks toward 1/2 (not GLA) in dosage mode
# ============================================================================

test_that("ancestry_split dosage mode shrinks toward 1/2 proportionally to w", {
  td <- tempfile()
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  # 3 samples, chr19 p-arm tract: S1 pure AFR, S2 mixed AFR/EUR, S3 pure AFR
  msp <- tempfile(fileext = ".tsv.gz", tmpdir = td)
  con <- gzfile(msp, "wt")
  writeLines(c(
    "#Subpopulation order/codes:\tAFR=0\tEUR=1",
    "#chm\tspos\tepos\tsgpos\tegpos\tn snps\tS1.0\tS1.1\tS2.0\tS2.1\tS3.0\tS3.1"
  ), con)
  writeLines(c(
    "chr19\t100\t500\t0.1\t0.2\t3\t0\t0\t0\t1\t0\t0"
  ), con)
  close(con)

  # Variant: S1 het (pure AFR carrier), S2 het (mixed, ambiguous), S3 hom-ref
  vcf <- tempfile(fileext = ".vcf", tmpdir = td)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
    "chr19\t200\t.\tG\tC\t.\tPASS\t.\tGT\t0|1\t0|1\t0|0"
  ), vcf)

  result <- ancestry_split(vcf, msp, mode = "dosage", chrom = "chr19", verbose = FALSE)

  # Shrinkage target is 1/2 regardless of the tract's ancestry mix (this
  # tract's GLA would be 5/6 AFR -- the old target -- and must NOT be used).
  # w = N5/total_alt: N2=1 (S1 pure AFR het), N5=1 (S2 mixed het), no hom-alt
  # carriers here so total_alt equals the unweighted carrier count (2) -> w=1/2
  # raw p1 = (0+1+0)/((1+1)-1) = 1
  w <- 0.5
  expected_afr <- (1 - w) * 1.0 + w * 0.5      # = 0.75 (GLA would give 11/12)
  s2_idx <- match("S2", result$sample_ids)
  expect_equal(result$AFR[1, s2_idx], expected_afr, tolerance = 1e-8)
  expect_equal(result$AFR[1, s2_idx] + result$EUR[1, s2_idx], 1, tolerance = 1e-8)
})

# ============================================================================
# K-population: each pair's w uses only that pair's evidence
#   w_ij = N_ij^(1) / (n_i + n_j + N_ij^(1)),  n_k = unambiguous copies on k
#   pi~_ij = (1 - w_ij) * n_i / (n_i + n_j) + w_ij / 2
#          = (n_i + N_ij^(1) / 2) / (n_i + n_j + N_ij^(1))
# ============================================================================

# one variant from per-state carrier counts; codes: 1..3 pure, 4 = 1/2,
# 5 = 1/3, 6 = 2/3; `n` is a named vector "<code>_<genotype>" -> people
k3_variant <- function(n) {
  code <- as.integer(sub("_.*", "", names(n))); g <- as.integer(sub(".*_", "", names(n)))
  list(gt = matrix(rep(g, n), nrow = 1), anc = matrix(rep(code, n), nrow = 1))
}
k3_pure  <- c(AFR = 1L, EUR = 2L, NAT = 3L)
k3_mixed <- data.frame(code = 4:6, pop1 = c("AFR", "AFR", "EUR"), pop2 = c("EUR", "NAT", "NAT"))
k3_split <- function(v) split_diploid_multi(v$gt, v$anc, k3_pure, k3_mixed, shrink = TRUE)
first_of <- function(v, code) which(v$anc[1, ] == code & v$gt[1, ] == 1L)[1]

test_that("split_diploid_multi: other ancestries' carriers do not lower a pair's w", {
  # pair AFR/EUR: one unambiguous AFR copy, none on EUR, 10 ambiguous AFR/EUR
  # hets, plus 100 pure-NAT carriers that say nothing about AFR vs EUR.
  # Per-pair w = 10 / (1 + 0 + 10) -> pi~ = (1 + 5) / 11 = 6/11; a variant-wide
  # T (1 + 100 + 10) would give w = 10/111 and pi~ ~ 0.95.
  v <- k3_variant(c(`1_1` = 1, `4_1` = 10, `3_1` = 100))
  out <- k3_split(v)
  j <- first_of(v, 4L)
  expect_equal(out$AFR[1, j], 6 / 11, tolerance = 1e-12)
  expect_equal(out$EUR[1, j], 5 / 11, tolerance = 1e-12)
  expect_equal(out$NAT[1, j], 0)
  # and the NAT carriers really are irrelevant: same split without them
  v0 <- k3_variant(c(`1_1` = 1, `4_1` = 10))
  expect_equal(k3_split(v0)$AFR[1, first_of(v0, 4L)], out$AFR[1, j], tolerance = 1e-12)
})

test_that("split_diploid_multi: three-way worked example (collaborators report, Table 3)", {
  # n = (12, 8, 4); ambiguous carriers small (1/2: 1), medium (1/3: 8),
  # dominating (2/3: 24); one mixed hom-alt each in 1/2 and 2/3
  v <- k3_variant(c(`1_2` = 1, `1_1` = 9, `2_1` = 6, `3_1` = 3,
                    `4_2` = 1, `4_1` = 1, `5_1` = 8, `6_2` = 1, `6_1` = 24))
  out <- k3_split(v)
  expect_equal(out$AFR[1, first_of(v, 4L)], (12 + 0.5) / (12 + 8 + 1))  # 0.5952, w = 1/21
  expect_equal(out$AFR[1, first_of(v, 5L)], (12 + 4) / (12 + 4 + 8))    # 0.6667, w = 1/3
  expect_equal(out$EUR[1, first_of(v, 6L)], (8 + 12) / (8 + 4 + 24))    # 0.5556, w = 2/3
  expect_equal(sum(out$AFR) + sum(out$EUR) + sum(out$NAT), sum(v$gt))   # every copy credited once
})

test_that("split_diploid_multi matches the per-pair formula on random 3-way variants", {
  set.seed(11)
  for (rep in 1:200) {
    n <- setNames(rpois(12, c(rep(c(0.3, 2), 3), rep(c(0.2, 1.5), 3))),
                  paste(rep(1:6, each = 2), c(2, 1), sep = "_"))
    if (sum(n) == 0) next
    v <- k3_variant(n); out <- k3_split(v)
    num <- c(2 * n[["1_2"]] + n[["1_1"]] + n[["4_2"]] + n[["5_2"]],
             2 * n[["2_2"]] + n[["2_1"]] + n[["4_2"]] + n[["6_2"]],
             2 * n[["3_2"]] + n[["3_1"]] + n[["5_2"]] + n[["6_2"]])
    pairs <- list(c(1, 2, 4), c(1, 3, 5), c(2, 3, 6))
    for (pr in pairs) {
      j <- first_of(v, pr[3]); if (is.na(j)) next
      amb <- n[[paste0(pr[3], "_1")]]; ni <- num[pr[1]]; nj <- num[pr[2]]
      expected <- (ni + amb / 2) / (ni + nj + amb)
      expect_equal(out[[names(k3_pure)[pr[1]]]][1, j], expected, tolerance = 1e-12)
    }
  }
})

test_that("split_diploid_multi with K = 2 equals split_diploid (2-population formula unchanged)", {
  set.seed(12)
  for (rep in 1:100) {
    nv <- 5; ns <- 40
    gt <- matrix(sample(0:2, nv * ns, TRUE, c(0.85, 0.12, 0.03)), nv)
    an <- matrix(sample(1:3, nv * ns, TRUE), nv)          # 1 = EUR/EUR, 2 = mixed, 3 = AFR/AFR
    two <- split_diploid(gt, an, shrink = TRUE)
    multi <- split_diploid_multi(gt, an, c(AFR = 3L, EUR = 1L),
                                 data.frame(code = 2L, pop1 = "AFR", pop2 = "EUR"),
                                 shrink = TRUE)
    expect_equal(multi$AFR, two$african, tolerance = 1e-12)
    expect_equal(multi$EUR, two$european, tolerance = 1e-12)
  }
})
