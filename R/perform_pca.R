#' Performing PCA analysis

#' @param data dataframe; data to be used for PCA. Sample metadata is in columns
#'     and the variables involved (i.e. compound or species)
#' @param variable string; the column name in which variables are stored
#                 (i.e. 'compound' or 'species')
#' @param value string; default 'value'; the column name in which data values
#'     are stored (i.e. 'concentration' or 'abundance')
#' @param convert.NA logical; where value is NA, convert value to 0;
#'     desirable when want to explicitly say missing values are true zeros
#'     (i.e. below limit of detection) and not missing from analysis
#' @param scale 'UV' for unit variance, 'UVN' for no unit variance, 'pareto';
#'     default UV
#' @param center logical; mean centering; default to TRUE
#' @param sampleID string. column name for sampleID. Default \code{'sample_name'}
#' @param print logic. Default \code{TRUE}. print preview of pca results
#'
#' @import reshape2
#' @return prcomp object; PCA-processed data
#' @export

perform_pca <- function(data, sampleID='sample_name', variable, value='value',
                       scale='UV', center = TRUE, convert.NA=FALSE,
                       print=TRUE) {

  #Getting column names of metadata and concentrations
  var <- as.character(unique(data[,variable]))
  d_colname = c(variable,value)
  meta_colname = colnames(data)[!colnames(data) %in% d_colname]

  #Building the formula needed to cast data to wide format
  form.y <- paste(meta_colname, collapse=' + ')

  form = as.formula(paste(form.y, variable, sep=' ~ '))

  #Casting data to wide format; includes qualifiers
  data.wide = reshape2::dcast(data=data, formula=form,
                    value.var=value)
  if(convert.NA==TRUE) {
    #converting NAs into 0
    data.wide[is.na(data.wide)] <- 0
  }

  #separating out concentration and qualifiers in wide format
  conc.wide = data.wide[, var]
  meta.wide = data.wide[, colnames(data.wide) %in% meta_colname]

  # check for empty compounds
  empty <- c()
  for(i in 1:ncol(conc.wide)){

    val <- unique(conc.wide[,i])
    if(length(val) == 1) {
      empty <- c(empty, colnames(conc.wide)[i])
      }
  }
  if(length(empty) > 0) {
    msg <- sprintf('Removing %s because it has has concentration of zero in every sample\n', empty)
    message(msg)

    conc.wide <- conc.wide[,!colnames(conc.wide) %in% empty]
  }
  #Converting concentration data to numeric class
  d <- as.data.frame(sapply(conc.wide, as.numeric))

  #Unit Variance Scaling, scaling weight of 1/sk
  if(scale == 'UV') {
    scale = TRUE
  }
  #No Unit Variance Scaling
  else if(scale == 'UVN') {
    scale = FALSE
  }
  #Pareto scaling, scaling weight of 1/(sk)e-2
  else if(scale == 'pareto') {
    scale =  for(column in 1:ncol(data)) {
      data[,column] = data[,column]/sqrt(sd(data[,column]))
    }
  }

  #Creating pca model
  pca_data = prcomp(d, scale=scale, center=center, na.action=na.omit)

  if(print==TRUE) print(summary(pca_data))

  #Adding qualifiers to pca score data
  pca_data$x = cbind(meta.wide, pca_data$x)

  #Melting pca score data back to long format,
  #where variable is PC number and value is score value on that PC
  pca.melt = melt(pca_data$x, id.vars=meta_colname)

  #Removing 'PC' in PC column names, so only left with a number
  pca.melt$variable = as.numeric(gsub('PC', '', pca.melt$variable))

  #Get rid of last PC if it is odd
  if(max(pca.melt$variable) %% 2 == 1) {
    pca.melt = pca.melt[pca.melt$variable != max(pca.melt$variable), ]
  }

  #Creating a new object where has columns of of odd number PCs and even number PCs
  odd.PCs = pca.melt$variable %% 2 == 1
  split = pca.melt[odd.PCs, meta_colname]
  split$x.value = pca.melt$value[odd.PCs]
  split$y.value = pca.melt$value[!odd.PCs]

  #Creating a column called PC, to give the PC number
  split$x.PC = pca.melt$variable[odd.PCs]
  split$y.PC = pca.melt$variable[!odd.PCs]

  #Adding column 'view' to describe which PCs are being plotted
  split$view = paste('PC',pca.melt$variable[odd.PCs], 'vs',
                     'PC',pca.melt$variable[!odd.PCs])

  #Replacing the pca model score data with the melted and sorted version created here
  pca_data$x = split

  return(pca_data)
  }
