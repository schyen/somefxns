% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pca_loading.R
\name{pca_loading}
\alias{pca_loading}
\title{Plotting PCA Loading Plot}
\usage{
pca_loading(pca_data, PCs = 1:4, label = NULL, title = NULL)
}
\arguments{
\item{pca_data}{the PCA model output of prcomp()}

\item{PCs}{vector of Principal Components to consider}

\item{label}{vector of compounds that should be labelled on the score plot;
to label all compounds (or bins), set label = 'all'}

\item{title}{string; title of loading plot}
}
\value{
list of: \code{p} ggplot object with generated plot
                 \code{collected_data} all the data used to plot pca data
                                  using ggplot2. can use collected_data for
                                  further customization using ggplot
}
\description{
Plotting PCA Loading Plot
}
