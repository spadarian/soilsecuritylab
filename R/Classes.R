#' An S4 class to represent a Cluster
#'
#' @slot data a \code{data.frame} with numerical data. Each column corresponds to a dimension of a \code{n}-dimensional clustering domain.
#' @slot U tmp
#' @slot W tmp
#' @slot centroids tmp
#' @slot phi a numerical value representing the degree of fuzziness or overlap of the generated clusters.
#' @slot classes numerical (integer) vector with the number of cluster to generate.
#' @slot distance used to calculate the clustering. Possible values are: 'Euclidean', 'Diagonal', and 'Mahalanobis'.
#' @slot alpha tmp
#' @slot Ue_mean-Ue_req tmp
#' @slot iterations total number of iterations
#' @slot pred_int tmp
setClass('FuzzyCluster',
         representation(data='matrix',U='matrix',W='matrix',centroids='matrix',phi='numeric',classes='integer',distance='character',
                                       alpha='numeric',`Ue_mean-Ue_req`='numeric',iterations='integer',pred_int='list'))

#' An S4 class to represent a group of FuzzyCluster objects
#'
#' @slot clusters A list containing FuzzyCluster objects
#' @slot confidence Confidence level for all the FuzzyCluster objects
setClass('FuzzyClusterGroup',representation(clusters='list',confidence='numeric'))
