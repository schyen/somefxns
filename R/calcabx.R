#' calcabx
#'
#' calculate antibiotic concentration based on MIC of reference curve and MIC
#' of sample curve
#'
#' @param std_max numeric. highest abx concentration indilution series tested
#' @param std_MIC_dil numeric. number of dilutions to reach abx MIC
#' @param fw_MIC_dil numeric. number of dilutions to reach fecal water MIC
#' @param sample_wt numeric. sample weight in grams
#' @param solvent_vol numeric. volume of solvent used to maek fecal slurry -- mL
#'
#' @return concentration fo antibiltic in sample. in ug abx / g poop

calcabx <- function(std_max, std_MIC_dil, fw_MIC_dil, sample_wt, solvent_vol) {

  # abx concentration at MIC -- ug
  abx_MIC <- std_max*(0.5^std_MIC_dil)

  # fecal slurry based on weight -- g/mL
  slurry <- sample_wt / (solvent_vol+sample_wt)

  # dilution factor
  dilfactor <- 0.5 ^ fw_MIC_dil

  #calculate fw concentration at fecal MIC -- g/mL
  fw_MIC <- slurry * dilfactor

  # abx in sample -- ug abx / g poop
  abx_sample <- abx_MIC / fw_MIC

  # multiple <- 100/fw_MIC
  # out <- std_MIC * multiple
  return(abx_sample)
}
