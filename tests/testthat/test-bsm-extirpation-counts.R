context("BSM range-contraction counts")

bsm_count_ana_fixture <- function() {
  data.frame(
    abs_event_time = 1:5,
    event_type = c("e", "e", "e", "d", "a"),
    current_rangetxt = c("AB", "AC", "AB", "A", "B"),
    new_rangetxt = c("B", "C", "A", "AB", "C"),
    ana_dispersal_from = c("-", "-", "-", "A", "B"),
    dispersal_to = c("-", "-", "-", "B", "C"),
    extirpation_from = c("A", "A", "B", "-", "B"),
    stringsAsFactors = FALSE
  )
}

bsm_count_clado_fixture <- function() {
  data.frame(
    time_bp = 6:9,
    clado_event_type = c("sympatry (y)", "subset (s)", "vicariance (v)", "founder (j)"),
    clado_event_txt = c("A->A,A", "AB->A,AB", "AB->A,B", "A->A,B"),
    clado_dispersal_from = c("", "", "", "A"),
    clado_dispersal_to = c("", "", "", "B"),
    stringsAsFactors = FALSE
  )
}

test_that("range contractions are counted by lost area, separately from dispersal", {
  counts <- count_ana_dispersal_events(
    bsm_count_ana_fixture(), c("A", "B", "C"), c("Alpha", "Beta", "Gamma")
  )
  expect_equal(as.numeric(counts$e_counts_df[1, ]), c(2, 1, 0))
  expect_equal(sum(counts$e_counts_df), 3)
  expect_equal(sum(counts$d_counts_df), 1)
  expect_equal(sum(counts$a_counts_df), 1)
  expect_equal(sum(counts$counts_df), 2)
  expect_equal(counts$d_counts_df["Alpha", "Beta"], 1)
  expect_equal(counts$a_counts_df["Beta", "Gamma"], 1)
})

test_that("range contractions do not require a dispersal event", {
  ana <- bsm_count_ana_fixture()
  counts <- count_ana_dispersal_events(ana[ana$event_type == "e", ], c("A", "B", "C"))
  expect_equal(as.numeric(counts$e_counts_df[1, ]), c(2, 1, 0))
  expect_equal(sum(counts$counts_df), 0)
})

test_that("no contractions and empty tables retain zero contraction counts", {
  ana <- bsm_count_ana_fixture()
  for (events in list(ana[ana$event_type != "e", ], ana[FALSE, ])) {
    counts <- count_ana_dispersal_events(events, c("A", "B", "C"))
    expect_equal(sum(counts$e_counts_df), 0)
    expect_equal(sum(counts$counts_df), nrow(events))
  }
})

test_that("range switching is counted once, not as a separate contraction", {
  ana <- bsm_count_ana_fixture()
  invisible(capture.output(
    counts <- count_ana_clado_events(
      list(bsm_count_clado_fixture()), list(ana[ana$event_type == "a", ]),
      c("A", "B", "C"), c("A", "B", "C")
    )
  ))
  expect_equal(counts$e_totals_list, 0)
  expect_equal(counts$a_totals_list, 1)
  expect_equal(counts$d_totals_list, 0)
  expect_equal(counts$anagenetic_dispersals_totals_list, 1)
  expect_equal(counts$ana_totals_list, 1)
  expect_equal(counts$all_totals_list, 5)
})

test_that("BSM event totals include contractions but dispersal totals do not", {
  ana <- bsm_count_ana_fixture()
  clado <- bsm_count_clado_fixture()
  invisible(capture.output(
    counts <- count_ana_clado_events(
      list(clado, clado), list(ana, ana[ana$event_type == "e", ]),
      c("A", "B", "C"), c("Alpha", "Beta", "Gamma")
    )
  ))
  expect_equal(counts$e_totals_list, c(3, 3))
  expect_equal(counts$d_totals_list, c(1, 0))
  expect_equal(counts$a_totals_list, c(1, 0))
  expect_equal(counts$anagenetic_dispersals_totals_list, c(2, 0))
  expect_equal(counts$founder_totals_list, c(1, 1))
  expect_equal(counts$all_dispersals_totals_list, c(3, 1))
  expect_equal(counts$ana_totals_list, c(5, 3))
  expect_equal(counts$clado_totals_list, c(4, 4))
  expect_equal(counts$all_totals_list, c(9, 7))
  expect_equal(counts$summary_counts_BSMs["means", "e"], 3)
  expect_equal(counts$summary_counts_BSMs["means", "all_ana"], 4)
  expect_equal(counts$summary_counts_BSMs["means", "total_events"], 8)
})
