% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_heatmap.R
\name{generate_heatmap}
\alias{generate_heatmap}
\title{explore_heatmap.R}
\usage{
generate_heatmap(data, x_var, y_var, value_var, x_dist = "euclidean",
  x_hclust = "complete", y_dist = "euclidean", y_hclust = "complete",
  x_meta = NULL, y_meta = NULL, x_label = NULL, y_label = NULL,
  heatmap_w = 20, heatmap_h = 20, dendro_x_h = 8, dendro_y_w = 8,
  output)
}
\arguments{
\item{data}{dataframe. long format, with column of xdata, variable
(such as taxonomy, or compounds), and value (such as abundance or concentration)}

\item{x_var}{string. variable to be heatmap x-axis. example: sample
should be a column name in data}

\item{y_var}{string. variable to be heatmap y-axis. example: species or compound
should be column name in data}

\item{value_var}{string. column containing heatmap values}

\item{x_dist}{string; dist method for x dendrogram;
method for computing phylogenetic distance;
supplied to the dist(); default method: "euclidean"
dist method options: "euclidean", "maximum", "manhattan",
"canberra", "binary" or "minkowski"}

\item{x_hclust}{string; hclust method for x dendrogram;
 Linkage criteria to determine distance between sets of
 observations; to be used in hclust(); default method: "complete"
 hclust method options:"ward.D", "ward.D2", "single", "complete",
"average" (= UPGMA), "mcquitty" (= WPGMA),
"median" (= WPGMC) or "centroid" (= UPGMC)}

\item{y_dist}{string; dist method for ydata dendrogram;
method for computing phylogenetic distance;
supplied to the dist(); default method: "euclidean"
dist method options: "euclidean", "maximum", "manhattan",
"canberra", "binary" or "minkowski"}

\item{y_hclust}{string; hclust method for ydata dendrogram;
Linkage criteria to determine distance between sets of
observations; to be used in hclust(); default method: "complete"
hclust method options:"ward.D", "ward.D2", "single", "complete",
"average" (= UPGMA), "mcquitty" (= WPGMA),
"median" (= WPGMC) or "centroid" (= UPGMC)
When set to FALSE, xdata listed along y-axis}

\item{x_meta}{dataframe. default NULL. Any metadata on heatmap x axis
variable (such as sample) to be included
x_var should be the identifier to link metadata}

\item{y_meta}{dataframe. default NULL. any metadata  on heatmap y axis
variable (such as species) to be included
y_var should be identifier to link metadata}

\item{x_label}{string. default NULL. Should be column in x_meta.
How x-axis of heatmap should be labelled. When set to NULL, x_var is used}

\item{y_label}{string. default NULL. Should be a column in y_meta.
How y-axis of heatmap should be labelled. When set to NULL, y_var is used}

\item{heatmap_w}{numeric; default 20; width (cm) of heatmap panel}

\item{heatmap_h}{numeric; default 20; height (cm) of heatmap panel}

\item{dendro_x_h}{numeric; default 8; height (cm) of x-axis dendrogram;
helps with ensuring x-axis dendrogram labels are not cut off}

\item{dendro_y_w}{numeric; default 8; width(cm) of y-axis dendrogram;
helps with ensuring y-axis dendrogram labels are not cut off}

\item{output}{= string; path and file name of resulting heatmap png
Example: "/path/to/my_heatmap.png"}
}
\value{
\code{pdata} dataframe for heatmap
    \code{ddata_x} dataframe for x-axis dendrogram
    \code{ddata_x} dataframe for y-axis dendrogram
    \code{p_hl} ggplot of heatmap and legend
    \code{p_x} ggplot of x-axis dendrogram
    \code{p_y} ggplot of y-axis dendrogram
    \code{g} gtable grob of entire heatmap and dendrograms
}
\description{
Produce heatmap of count data. uses euclidean distances (recommdended) of clr transformed data; zero replacement completed by Bayesian-multiplicative replacement
}
