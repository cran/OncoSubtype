#' Set the SubtypeClass
#'
#' @return an object of SubtypeClass with three empty solts
SubtypeClass <- setClass('SubtypeClass', slots = c(features = "character",
                                    subtypes = "character",
                                    score = "ANY",
                                    train_set = "data.frame",
                                    test_set = "data.frame",
                                    method = "ANY"),
                         contains = "character")
