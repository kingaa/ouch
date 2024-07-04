##' Coerce an \pkg{ouch} object to a data frame
##'
##' @name as_data_frame
##' @rdname as_data_frame
##' @family methods for ouch trees
##' @include ouchtree.R brown.R hansen.R
##' @inheritParams base::as.data.frame
##'
NULL

##' @rdname as_data_frame
##' @export
as.data.frame.ouchtree <- function (x, ...) as(x,"data.frame")

##' @rdname as_data_frame
##' @export
as.data.frame.browntree <- function (x, ...) as(x,"data.frame")

##' @rdname as_data_frame
##' @export
as.data.frame.hansentree <- function (x, ...) as(x,"data.frame")
