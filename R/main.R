#' Predict the subtypes of selected cancer type based published papers
#'
#' @param data data set to predict the subtypes which is a numeric matrix with row names of features and column names of samples
#' @param disease character string of the disease to predict subtypes, currently support 'LUSC', 'LUAD'
#' @return an object of class "SubtypeClass" with four slots: genes used for predictiong, predicted subtypes of samples, a matrix of predicting scores, and the method.
#' @export
#' @import SummarizedExperiment
#' @import e1071
#' @importFrom methods new
#' @importFrom stats cor median predict
#' @examples
#' \dontrun{
#' library(OncoSubtype)
#' data <- get_median_centered(example_fpkm)
#' data <- assays(data)$centered
#' rownames(data) <- rowData(example_fpkm)$external_gene_name
#' centroids_subtype(data, disease = 'HNSC')
#' }

centroids_subtype <- function(data, disease = 'LUSC') {
  centroids_data <- get(paste0(tolower(disease), '_centroids'))
  output <- new("SubtypeClass", features = intersect(rownames(data), rownames(lusc_centroids)))
  test_data <- data[output@features, ]
  predict_data <- centroids_data[output@features, ]
  cor_scores <- do.call(rbind.data.frame, lapply(colnames(test_data), function(x) {
    cor(test_data[,x], predict_data, use = 'complete') }))
  rownames(cor_scores) <- colnames(test_data)
  subtypes <- colnames(cor_scores)[apply(cor_scores, 1, which.max)]
  names(subtypes) <- colnames(test_data)
  output@subtypes <- as.character(subtypes)
  output@score <- as.matrix(cor_scores)
  output@train_set <- data.frame(predict_data)
  output@test_set <- data.frame(test_data)
  output@method <- 'Nearest centroids'
  return(output)
}

#' select highly variable genes from a expression matrix
#'
#' @param data a (normalized) matrix with rownames of features and colnames of samples
#' @param top number of top highly variable genes to output
#' @return subset with top ranked genes by the variances
#' @export
#' @import SummarizedExperiment
#' @importFrom stats var
#' @examples
#' \dontrun{
#' library(OncoSubtype)
#' data <- get_median_centered(example_fpkm)
#' data <- assays(data)$centered
#' get_hvg(data)
#' }
get_hvg <- function(data, top = 1000){
  varorder <- order(apply(data, 1, var, na.rm = TRUE), decreasing=T)[1:top]
  data[varorder,]
}

#' convert expression matrix to median-centered
#'
#' @param data a numeric matrix or `S4` object
#' @param log2 logical, if `TRUE` \eqn{log2(x+1)} will be applied before median centering. Defaut `TRUE`.
#' @return median-centered express matrix or an object with new slot "centered"
#' @export
#' @import SummarizedExperiment
#' @examples
#' \dontrun{
#' get_median_centered(example_fpkm)
#'}
#'
get_median_centered <- function(data, log2 = TRUE) {
  if(is.matrix(data)){
    if(log2){
      data <- log2(assay(data) + 1)
    }
    return (sweep(data, 1,  apply(data, 1, median, na.rm = TRUE)))
  }
  if(typeof(data) == 'S4'){
    edata <- assays(data)[[1]]
    if(log2){
      edata <- log2(edata + 1)
    }
    assays(data)$centered <- sweep(edata, 1,  apply(edata, 1, median, na.rm = TRUE))
    return(data)
  }
}

#' Predict the subtypes of selected cancer type
#'
#' @param train_set training set with rownames of samples, first column named 'mRNA_subtype' and the rest of features and expression values.
#' @param test_set test set with rownames of features and colnames of samples.
#' @param method character string of the method to use currently support 'rf'.
#' @param seed integer seed to use.
#' @return a matrix with column names of subtypes and predicted probabilities.
#' @import caret
#' @import randomForest

get_rf_pred <- function(train_set, test_set, method = 'rf', seed = NULL){

  if (!is.null(seed)) set.seed(seed)

  trControl = trainControl(method = 'repeatedcv',
                           number = 4,
                           repeats = 5,
                           savePredictions = 'final')

  rf_model <- train(mRNA_subtype  ~ . ,
                    data = train_set,
                    trControl = trControl,
                    method = method,
                    tuneLength = 5)

  out_pred <- predict(rf_model, newdata = t(test_set), type = 'prob')

  return(list(model = rf_model, results = out_pred))
}

#' Predict the subtypes of selected cancer type using machine learning
#'
#' @param data data set to predict the subtypes which is a numeric matrix with row names of features and column names of samples
#' @param disease character string of the disease to predict subtypes, currently support 'LUSC', 'LUAD', and 'BLCA'.
#' @param method character string of the method to use currently support 'rf'.
#' @param seed integer seed to use.
#' @param removeBatch whether do batch effect correction using \code{limma::removeBatchEffect}, default TRUE.
#' @return An object of class "SubtypeClass" with four slots: genes used for predictiong, predicted subtypes of samples, a matrix of predicting scores, and the method.
#' @references
#' \enumerate{
#' \item \insertRef{wilkerson2010}{OncoSubtype}
#'
#' \item \insertRef{wilkerson2012}{OncoSubtype}
#'
#' \item \insertRef{tcga2015}{OncoSubtype}
#'}
#' @export
#' @import SummarizedExperiment
#' @import caret randomForest e1071 Rdpack
#' @importFrom methods new
#' @importFrom limma removeBatchEffect
#' @importFrom stats cor median predict
#' @examples
#' \dontrun{
#' library(OncoSubtype)
#' data <- get_median_centered(example_fpkm)
#' data <- assays(data)$centered
#' rownames(data) <- rowData(example_fpkm)$external_gene_name
#' ml_subtype(data, disease = 'LUAD', method = 'rf', seed = 123)
#' }

ml_subtype <- function(data, disease = 'LUSC', method = 'rf', removeBatch = TRUE, seed = NULL) {
  # Load the dataset and capture the path of the downloaded dataset
  dataset_path <- load_dataset_from_github(disease = disease)

  # Assume load_dataset_from_github loads the dataset into a variable named `train_data`
  # Note: Adjust the name `train_data` based on the actual name used within your .rda files
  train_data_env <- new.env()
  load(dataset_path, envir = train_data_env)
  train_data <- get(paste0(tolower(disease), '_tcga'), envir = train_data_env)

  output <- new("SubtypeClass", features = intersect(rownames(data), colnames(train_data)[-1]))

  if (method == 'rf') {
    test_data <- data[output@features, ]
    predict_data <- train_data[, c('mRNA_subtype', output@features)]

    if (removeBatch) {
      batch <- c(rep('test', ncol(test_data)), rep('train', nrow(predict_data)))
      mat_batch <- cbind(test_data, t(predict_data[,-which(colnames(predict_data)=='mRNA_subtype')]))
      mat_rmbatch <- removeBatchEffect(mat_batch, batch)
      test_data_new <- mat_rmbatch[, batch == 'test']
      predict_data_new <- t(mat_rmbatch[, batch == 'train'])
      if (identical(rownames(predict_data_new), rownames(predict_data))) {
        predict_data_new <- cbind.data.frame(mRNA_subtype = predict_data[,'mRNA_subtype'], predict_data_new)
      }
      predict_data <- predict_data_new
      test_data <- test_data_new
    }

    out_pred <- get_rf_pred(train_set = predict_data, test_set = test_data, seed = seed)
    output@subtypes <- colnames(out_pred$results)[apply(out_pred$results, 1, which.max)]
    output@score <- out_pred$results
    output@train_set <- data.frame(predict_data)
    output@test_set <- data.frame(test_data)
    output@method <- out_pred$model
  }

  # Auto-delete the downloaded dataset file after processing
  unlink(dataset_path)

  return(output)
}

#' Plot heatmap of the train set or test set
#'
#' @param object a SubtypeClass object
#' @param set options could be 'test', 'train' or 'both'. Default 'test'.
#' @param ...  Parameters passed to \code{pheatmap}.
#' @return a pheatmap object
#' @export
#' @import pheatmap
#' @import tibble
#' @importFrom rlang .data
#' @importFrom dplyr mutate filter select relocate arrange bind_rows
#' @examples
#' \dontrun{
#' library(OncoSubtype)
#' data <- get_median_centered(example_fpkm)
#' data <- assays(data)$centered
#' rownames(data) <- rowData(example_fpkm)$external_gene_name
#' object <- MLSubtype(data, disease = 'LUSC')
#' PlotHeat(object, set = 'both', fontsize = 10, show_rownames = FALSE, show_colnames = FALSE)
#' }

PlotHeat <- function(object, set = 'test', ...) {
  test_set <- t(object@test_set) %>% data.frame %>%
    rownames_to_column('Samples') %>%
    mutate(mRNA_subtype = object@subtypes,
           set = 'test_set') %>%
    relocate(.data$Samples, .data$set, .data$mRNA_subtype, object@features)
  train_set <- object@train_set %>% data.frame %>%
    rownames_to_column('Samples') %>%
    mutate(set = 'train_set') %>%
    relocate(.data$Samples, .data$set, .data$mRNA_subtype, object@features)

  joint_set <- bind_rows(test_set, train_set) %>%
    arrange(.data$set, .data$mRNA_subtype)

  rownames(joint_set) <- make.unique(joint_set$Samples)

  if (set == 'both') {
    plotdat <- joint_set
  }

  if (set == 'test') {
    plotdat <- dplyr::filter(joint_set,  set == 'test_set')
  }

  if (set == 'train') {
    plotdat <- dplyr::filter(joint_set,  set == 'train_set')
  }

  anno <- plotdat[, c('set', 'mRNA_subtype')]

  pheatmap(t(plotdat[, object@features]),
           cluster_cols = FALSE,
           annotation_col = anno, ...)
}


#' Load Dataset from GitHub Repository
#'
#' Downloads a specified dataset from a GitHub repository if it is not already
#' present in the specified local directory, then loads the dataset into the global
#' environment. This function is designed to help manage package size by storing
#' data externally and loading it on-demand.
#'
#' @param disease A character string specifying the disease, which corresponds
#'   to the name of the dataset to be loaded (e.g., "LUSC"). The function constructs
#'   the filename as \code{tolower(disease)_tcga.rda} and attempts to load this dataset.
#' @param local_dir An optional character string specifying the path to the directory
#'   where datasets should be stored locally. If not provided, defaults to a
#'   subdirectory named \code{your_package_name_data} within the user's home directory.
#'   Users can specify their own directory path if they prefer to store data in a different
#'   location.
#'
#' @return Invisible NULL. The function is primarily used for its side effect of
#'   loading a dataset into the global environment. However, the function itself
#'   does not return the dataset directly.
#'
#' @examples
#' \dontrun{
#'   load_dataset_from_github("LUSC")
#' }
#'
#' @export
#' @importFrom utils download.file

load_dataset_from_github <- function(disease, local_dir = path.expand(getwd())) {
  # Ensure the directory exists
  if (!dir.exists(local_dir)) {
    dir.create(local_dir, recursive = TRUE)
  }

  # Construct the filename and URL for the dataset
  dataset_name <- paste0(tolower(disease), '_tcga.rda')
  # Update this URL to point to the raw content
  url <- paste0("https://raw.githubusercontent.com/DadongZ/OncoSubtype_data/main/data/", dataset_name)

  # Determine the local file path
  local_dataset_path <- file.path(local_dir, dataset_name)

  # Download the file if it does not exist locally
  if (!file.exists(local_dataset_path)) {
    message("Downloading dataset: ", dataset_name)
    download.file(url, destfile = local_dataset_path, mode = "wb")
  } else {
    message("Dataset already downloaded.")
  }

  # Load the dataset into the global environment
  load(local_dataset_path)

  return(local_dataset_path)
}
