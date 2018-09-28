#' Plotting PCA Loading Plot

#' @param pca_data the PCA model output of prcomp()
#' @param PCs vector of Principal Components to consider
#' @param label vector of compounds that should be labelled on the score plot;
#'     to label all compounds (or bins), set label = 'all'
#' @param title string; title of loading plot
#' @return list of: \code{p} ggplot object with generated plot
#'                  \code{collected_data} all the data used to plot pca data
#'                                   using ggplot2. can use collected_data for
#'                                   further customization using ggplot
#' @import ggplot2

pca_loading = function(pca_data, PCs=1:4, label=NULL, title=NULL)
{
  #------------------------------------------------------->
  #Input testing

  #Testing if the pca_data is indeed the result of
  #a prcomp function
  if (class(pca_data) != 'prcomp')
  {
    stop('Input "pca_data" to pca_loading must be the output of prcomp() function')
  }

  #Testing if PCs are numeric
  if (!is.numeric(PCs))
  {
    stop('Input "PCs" to pca_loading must be a vector of numeric values corresponding to target principal components')
  }

  #Testing to make sure that PCs are in the input data
  if (any(PCs > ncol(pca_data$rotation)))
  {
    stop('Input "PCs" to pca_loading must not be outside of the pca_data principal component range')
  }

  #Checking if given labels are in the loading data
  if (any(!label %in% rownames(pca_data$rotation)))
  {
    if (label != 'all') {
      stop(paste('The following labels do not correspond to variables in "pca_data":\n\t',
                 paste(label[!label %in% rownames(pca_data$rotation)], collapse=', ')))
    }

  }
  #------------------------------------------------------->

  #Extract loadings
  loading_data <- as.data.frame(pca_data$rotation)

  variance <- as.data.frame(pca_data$sdev)
  #Converting standard deviations to variance
  variance <- variance**2
  variance <- variance/sum(variance)*100
  variance <- round(variance, digits=1)

  #Initializing a new vector to store collected_data
  collected_data = c()

  #Iterating through the PCs, two at a time.
  #If there is an odd number of PCs, the last one is dropped
  for (i in seq(1, length(PCs)-1, by=2))
  {
    view_data <- loading_data[, c(i, i+1)]
    view_data$View <- paste('PC', i, '(', variance[i,], '%) ','vs.',
                           'PC', i+1,'(', variance[i+1,], '%)')

    #Changing score_column names
    colnames(view_data)[1:2] = c('xPC', 'yPC')

    #Initializing text column with blanks
    view_data$Text = ''

    if(label != 'all') {
      #Adding labels
      view_data[label, 'Text'] = label

      #Adding to collected_data
      collected_data = rbind(collected_data, view_data)
    }
    else if(label == 'all') {
      view_data$Text <- rownames(view_data)
    }
    collected_data <- rbind(collected_data, view_data)
  }

  # #Calculating plot limits
  # xlimit = max(abs(collected_data$xPC))*1.1
  # ylimit = max(abs(collected_data$yPC))*1.1

  collected_data <- unique(collected_data)

  #Generating plot
  p = ggplot(collected_data)
  p = p + geom_point(aes(x=xPC, y=yPC), colour='black', size=3, alpha=.7)

  #Adding text only if the "label" descriptor has been provided,
  #the points with text are also highlighted
  if (!is.null(label))
  {
    p = p + ggrepel::geom_text_repel(aes(x=xPC, y=yPC, label=Text),
                      hjust=-0.1, vjust=0, size=3, fontface=2)

  }

  #Adding title only if "title" descriptor has been provided
  if (!is.null(title)) p = p + ggtitle(label=title)

  #Drawing origin
  p = p + geom_vline(xintercept=0, alpha=0.3)
  p = p + geom_hline(yintercept=0, alpha=0.3)

  p = p + xlab('PC loading') + ylab('PC loading')

  #Setting limits
  # p = p + xlim(-xlimit,xlimit)
  # p = p + ylim(-ylimit,ylimit)

  # facetting through PC pairs
  p = p + facet_wrap(~ View, ncol=min(floor(length(PCs)/2), 3)) +
    theme_bw(16) +
    theme(axis.title.x = element_text(size=14, face='bold'),
          axis.title.y = element_text(size=14, face='bold'),
          axis.text = element_text(size=12, colour='black'),
          legend.title = element_text(size=14),
          legend.text = element_text(size=12),
          legend.key = element_rect(colour='white'),
          legend.spacing = unit(0, "cm"),
          panel.border=element_rect(size=2, colour='black'),
          strip.text.x = element_text(size=14, face='bold'))

  return(list('p'=p, 'collected_data'=collected_data))
}
