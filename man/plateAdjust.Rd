% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plateAdjust.R
\name{plateAdjust}
\alias{plateAdjust}
\title{plateAdjust}
\usage{
plateAdjust(plateDF, metadata, well_include = NULL, blank_by = "row",
  blank_label = "blank")
}
\arguments{
\item{plateDF}{dataframe. first spreadsheet of victor file}

\item{metadata}{dataframe. well metadata. must have columns:
Well, platerow, platecol, strain, curveID, welltype, media, abx, wellconc}

\item{well_include}{string or vector of string. Default NULL. Which wells to
include in analysis}

\item{blank_by}{\code{c('row','col', 'T0')}. Default \code{'row'}.
Sets orientation of blanks.
\code{'row'} means blanks are in the same row as corresponding wells
\code{'col'} means blanks are in the same column as corresponding wells
\code{'timepoint'} means blank by reading at time 0 (aka repeat 1)}

\item{blank_label}{string. Default \code{'blank'}.
How are blanks identified in metadata}
}
\value{
full dataframe of with metadata, absorbance, adjusted absorbance in one dataframe
}
\description{
normalized absorbance values by adjusting to blank
}
