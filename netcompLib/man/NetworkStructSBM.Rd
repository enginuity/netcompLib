% Generated by roxygen2 (4.0.1): do not edit by hand
\name{NetworkStructSBM}
\alias{NetworkStructSBM}
\title{Constructor for RND network structure}
\usage{
NetworkStructSBM(model_params = set_model_param(), NetM = NULL)
}
\arguments{
\item{model_params}{[list; DEFAULT = \code{\link{set_model_param}}()] :: Model parameters}

\item{NetM}{[\code{\link{NetworkModelSBM}}] :: If not NULL, this model is used, and its structure extracted}
}
\value{
[NetworkStructSBM] :: A representation of the generated model structure
}
\description{
This will generate a network model with random class assignments and class probabilities (unless otherwise specified)
}

