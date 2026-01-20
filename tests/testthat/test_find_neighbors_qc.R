path_to_counts <- system.file("extdata", "pbmc_raw.txt", package = "Seurat")

# Helper function to build test Seurat object with QC metrics
get_test_data_with_qc <- function() {
  raw_counts <- read.table(path_to_counts, sep = "\t", row.names = 1)
  counts <- as.sparse(as.matrix(raw_counts))
  assay <- CreateAssay5Object(counts)
  test_data <- CreateSeuratObject(assay)

  test_data <- NormalizeData(test_data, verbose = FALSE)
  test_data <- FindVariableFeatures(test_data, verbose = FALSE)
  test_data <- ScaleData(test_data, verbose = FALSE)
  test_data <- RunPCA(test_data, npcs = 10, verbose = FALSE)

  return(test_data)
}

context("FindNeighborsQC")

test_that("FindNeighborsQC returns Seurat object with correct graph names", {
  test_case <- get_test_data_with_qc()

  # Run FindNeighborsQC with nFeature_RNA as QC metric
  result <- FindNeighborsQC(
    test_case,
    qc.metric = "nFeature_RNA",
    k.min = 5,
    k.max = 20,
    dims = 1:10,
    verbose = FALSE
  )

  # Check that result is a Seurat object
  expect_s4_class(result, "Seurat")

  # Check that QC_nn and QC_snn graphs are present
  expect_true("QC_nn" %in% names(result))
  expect_true("QC_snn" %in% names(result))

  # Check that the graphs have correct dimensions
  n_cells <- ncol(test_case)
  expect_equal(nrow(result[["QC_nn"]]), n_cells)
  expect_equal(ncol(result[["QC_nn"]]), n_cells)
  expect_equal(nrow(result[["QC_snn"]]), n_cells)
  expect_equal(ncol(result[["QC_snn"]]), n_cells)
})

test_that("FindNeighborsQC falls back to FindNeighbors when qc.metric is NULL", {
  test_case <- get_test_data_with_qc()

  # Run with qc.metric = NULL
  result <- FindNeighborsQC(
    test_case,
    qc.metric = NULL,
    k.param = 10,
    dims = 1:10,
    verbose = FALSE
  )

  # Should use default assay naming convention
  default_assay <- DefaultAssay(test_case)
  expect_true(paste0(default_assay, "_nn") %in% names(result))
  expect_true(paste0(default_assay, "_snn") %in% names(result))
})

test_that("FindNeighborsQC stores k values in metadata", {
  test_case <- get_test_data_with_qc()

  result <- FindNeighborsQC(
    test_case,
    qc.metric = "nFeature_RNA",
    k.min = 5,
    k.max = 20,
    dims = 1:10,
    verbose = FALSE
  )

  # Check that k values are stored
  expect_true("nFeature_RNA_k" %in% colnames(result[[]]))

  # Check that k values are within bounds
  k_values <- result[["nFeature_RNA_k", drop = TRUE]]
  expect_true(all(k_values >= 5))
  expect_true(all(k_values <= 20))
})

test_that("FindNeighborsQC validates input parameters", {
  test_case <- get_test_data_with_qc()

  # Invalid qc.metric column
  expect_error(
    FindNeighborsQC(
      test_case,
      qc.metric = "nonexistent_column",
      dims = 1:10,
      verbose = FALSE
    ),
    "not found in object metadata"
  )

  # k.min must be at least 1
  expect_error(
    FindNeighborsQC(
      test_case,
      qc.metric = "nFeature_RNA",
      k.min = 0,
      k.max = 20,
      dims = 1:10,
      verbose = FALSE
    ),
    "k.min must be at least 1"
  )

  # k.max must be >= k.min
  expect_error(
    FindNeighborsQC(
      test_case,
      qc.metric = "nFeature_RNA",
      k.min = 30,
      k.max = 20,
      dims = 1:10,
      verbose = FALSE
    ),
    "k.max must be greater than or equal to k.min"
  )
})

test_that("FindNeighborsQC works with custom qc.mapping function", {
  test_case <- get_test_data_with_qc()

  # Logarithmic mapping function
  log_mapping <- function(x) log1p(x * 9) / log(10)

  result <- FindNeighborsQC(
    test_case,
    qc.metric = "nFeature_RNA",
    qc.mapping = log_mapping,
    k.min = 5,
    k.max = 20,
    dims = 1:10,
    verbose = FALSE
  )

  expect_s4_class(result, "Seurat")
  expect_true("QC_nn" %in% names(result))
  expect_true("QC_snn" %in% names(result))
})

test_that("FindNeighborsQC works with nCount_RNA metric", {
  test_case <- get_test_data_with_qc()

  result <- FindNeighborsQC(
    test_case,
    qc.metric = "nCount_RNA",
    k.min = 5,
    k.max = 20,
    dims = 1:10,
    verbose = FALSE
  )

  expect_s4_class(result, "Seurat")
  expect_true("QC_nn" %in% names(result))
  expect_true("nCount_RNA_k" %in% colnames(result[[]]))
})

test_that("FindNeighborsQC output is compatible with FindClusters", {
  test_case <- get_test_data_with_qc()

  # Build adaptive neighbor graph
  test_case <- FindNeighborsQC(
    test_case,
    qc.metric = "nFeature_RNA",
    k.min = 5,
    k.max = 20,
    dims = 1:10,
    verbose = FALSE
  )

  # Run clustering on the QC-adaptive SNN graph
  result <- FindClusters(
    test_case,
    graph.name = "QC_snn",
    resolution = 0.5,
    verbose = FALSE
  )

  # Check that clustering succeeded
  expect_true("seurat_clusters" %in% colnames(result[[]]))
  expect_false(any(is.na(result$seurat_clusters)))
})

test_that("FindNeighborsQC handles equal QC values gracefully", {
  test_case <- get_test_data_with_qc()

  # Set all nFeature_RNA values to the same value
  test_case[["constant_qc"]] <- rep(100, ncol(test_case))

  result <- FindNeighborsQC(
    test_case,
    qc.metric = "constant_qc",
    k.min = 5,
    k.max = 20,
    dims = 1:10,
    verbose = FALSE
  )

  expect_s4_class(result, "Seurat")

  # All k values should be the same (midpoint)
  k_values <- result[["constant_qc_k", drop = TRUE]]
  expect_true(length(unique(k_values)) == 1)
})

test_that("FindNeighborsQC respects custom graph names", {
  test_case <- get_test_data_with_qc()

  result <- FindNeighborsQC(
    test_case,
    qc.metric = "nFeature_RNA",
    k.min = 5,
    k.max = 20,
    dims = 1:10,
    graph.name = c("custom_nn", "custom_snn"),
    verbose = FALSE
  )

  expect_true("custom_nn" %in% names(result))
  expect_true("custom_snn" %in% names(result))
})

test_that("FindNeighborsQC works without SNN computation", {
  test_case <- get_test_data_with_qc()

  result <- FindNeighborsQC(
    test_case,
    qc.metric = "nFeature_RNA",
    k.min = 5,
    k.max = 20,
    dims = 1:10,
    compute.SNN = FALSE,
    graph.name = "QC_nn_only",
    verbose = FALSE
  )

  expect_true("QC_nn_only" %in% names(result))
})
