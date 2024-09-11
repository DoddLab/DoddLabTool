################################################################################
# plot_ms2 ---------------------------------------------------------------------

#' @title plot_ms2
#' @description plot ms2 spec
#' @author Zhiwei Zhou
#' @param ms2_spec a table of spectrum. 1st column "mz", 2nd column "intensity"
#' @param is_normalize whether normalize the spectrum
#' @export
#' @examples
#' \dontrun{
#' ms2_spec <- get_ms2(lab_id = 'S00001', ce = '20', polarity = 'positive')
#' plot_ms2(ms2_spec, is_normalize = FALSE)
#' plot_ms2(ms2_spec, is_normalize = TRUE)
#' }


# ms2_spec <- get_ms2(lab_id = 'S00001', ce = '20', polarity = 'positive')
# plot_ms2(ms2_spec, is_normalize = FALSE)
# plot_ms2(ms2_spec, is_normalize = TRUE)

plot_ms2 <- function(
    ms2_spec,
    is_normalize = TRUE
) {
  if (class(ms2_spec)[1] == 'matrix') {
    ms2_spec <- ms2_spec %>% as.data.frame()
  }

  if (colnames(ms2_spec)[1] != "mz" | colnames(ms2_spec)[2] != "intensity") {
    message(crayon::red("The first 2 columns should be mz and intensity\n"))
    stop()
  }

  if (is_normalize) {
    ms2_spec$intensity <- ms2_spec$intensity/max(ms2_spec$intensity)
  }

  temp_plot <- ggplot2::ggplot(ms2_spec) +
    ggplot2::geom_segment(ggplot2::aes(x = mz, xend = mz,
                                       y = 0, yend = intensity)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::xlim(0.95*min(ms2_spec$mz),
                  1.05*max(ms2_spec$mz)) +
    ggplot2::xlab('m/z') +
    ggplot2::ylab(ifelse(is_normalize, 'Rela. Intensity', 'Intensity')) +
    ZZWtool::ZZWTheme() +
    # ggplot2::ggtitle(label = paste0(obj_ms2@info$name, '{',
    #                                 round(obj_ms2@info$scoreForward, 4), '}')) +
    ggplot2::theme(legend.position = c(0.8, 0.75),
                   title = ggplot2::element_text(vjust = 0.5))


  return(temp_plot)
}


# plot_id_ms2 ------------------------------------------------------------------
#' @title plot_id_ms2
#' @description generate ms2 mirror plot
#' @author Zhiwei Zhou
#' @param obj_ms2 object of SpectraData
#' @param obj_spec
#' @export
#' @examples

plot_id_ms2 <- function(
    obj_ms2 = NULL,
    obj_spec = NULL
) {
  # mode <- match(mode)

  if (length(obj_ms2) > 0) {
    obj_spec <- obj_ms2@matchedFragments[[1]]
  }

  if (length(obj_spec) == 0) {
    stop('Please input obj_spec')
  }

  temp_spec <- obj_spec %>%
    # dplyr::mutate(int = tidyr::replace_na(int, 0)) %>%
    dplyr::mutate(int_lib = intensity/max(intensity),
                  int_exp = intensityExp/max(intensityExp)) %>%
    dplyr::rename(mz_lib = mz,
                  mz_exp = mzExp) %>%
    dplyr::select(mz_lib, int_lib, mz_exp, int_exp, fragPrecursor) %>%
    dplyr::mutate(label = dplyr::case_when(
      int_lib > 0 & int_exp > 0 ~ 'matched',
      !(int_lib > 0 & int_exp > 0) ~ 'unmatched'
    ))

  # switch (mode,
  #         'SpecLibMatch' = {
  #
  #         },
  #         'NeighborMatch' = {
  #
  #         }
  # )

  temp_spec1 <- temp_spec %>%
    dplyr::select(mz_lib, int_lib, fragPrecursor, label) %>%
    dplyr::rename(mz = mz_lib, int = int_lib) %>%
    # dplyr::mutate(int = 0-int,
    #               tag = 'library')
    dplyr::mutate(int = 0-int,
                  tag = dplyr::case_when(
                    label == 'unmatched' ~ 'frag_unmatch',
                    label == 'matched' ~ 'library'
                  ))

  temp_spec2 <- temp_spec %>%
    dplyr::select(mz_exp, int_exp, fragPrecursor, label) %>%
    dplyr::rename(mz = mz_exp, int = int_exp) %>%
    # dplyr::mutate(tag = 'experiment')
    dplyr::mutate(tag = dplyr::case_when(
      label == 'unmatched' ~ 'frag_unmatch',
      label == 'matched' ~ 'experiment'
    ))


  temp_data <- temp_spec1 %>%
    dplyr::bind_rows(temp_spec2)


  temp_plot <- ggplot2::ggplot(temp_data) +
    ggplot2::geom_segment(ggplot2::aes(x = mz, xend = mz,
                                       y = 0, yend = int,
                                       colour = tag)) +
    ggplot2::geom_point(ggplot2::aes(x = mz,
                                     y = int,
                                     shape = label,
                                     colour = tag))+
    ggplot2::scale_colour_manual(values = c(
      'experiment' = 'dodgerblue',
      'library' = 'tomato',
      'frag_unmatch' = 'gray'
    )) +
    ggplot2::scale_shape_manual(values = c(
      'matched' = 16,
      'unmatched' = 4
    )) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::xlim(0.95*min(temp_data$mz),
                  1.05*max(temp_data$mz)) +
    ggplot2::xlab('m/z') +
    ggplot2::ylab('Relative intensity') +
    ZZWtool::ZZWTheme() +
    # ggplot2::ggtitle(label = paste0(obj_ms2@info$name, '{',
    #                                 round(obj_ms2@info$scoreForward, 4), '}')) +
    ggplot2::theme(legend.position = c(0.8, 0.75),
                   title = ggplot2::element_text(vjust = 0.5))


  return(temp_plot)
}



# plot_id_shift_ms2 -----------------------------------------------------------
#' @title plot_id_shift_ms2
#' @description generate ms2 mirror plot for shift match (bonanza, hybrid, gnps)
#' @author Zhiwei Zhou
#' @param obj_ms2 object of SpectraData
#' @export
#' @examples

# load('/home/zhouzw/Data_processing/20210224_metdna2_development_test/obj_ms2_cpd1_for_metdna2.RData')
# load('/home/zhouzw/Data_processing/20210224_metdna2_development_test/obj_ms2_cpd2_for_metdna2.RData')
# score_dp <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'dp')
# score_bonanza <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'bonanza')
# score_hybrid <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'hybrid')
# score_gnps <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'gnps')
#
# plotIdShiftMs2(obj_ms2 = score_bonanza)
# plotIdShiftMs2(obj_ms2 = score_hybrid)


plot_id_shift_ms2 <- function(
    obj_ms2 = NULL,
    obj_spec_frag_match,
    obj_spec_frag_nl
) {

  if (length(obj_ms2) > 0) {
    obj_spec_frag_match <- obj_ms2@matchedFragments[[1]]
    obj_spec_frag_nl <- obj_ms2@nlFragments[[1]]
  }

  if (length(obj_spec_frag_match) == 0 & length(obj_spec_frag_nl) == 0) {
    stop('Please input obj_spec_frag_match and obj_spec_frag_nl')
  }


  obj_spec_frag_match <- obj_spec_frag_match %>%
    dplyr::mutate(type = 'exact_match')
  obj_spec_frag_nl <- obj_spec_frag_nl %>%
    dplyr::mutate(type = 'nl_match')

  temp_spec_frag_match <- obj_spec_frag_match %>%
    dplyr::rename(mz_lib = mz,
                  mz_exp = mzExp,
                  int_lib = intensity,
                  int_exp = intensityExp) %>%
    dplyr::select(mz_lib, int_lib, mz_exp, int_exp, fragPrecursor, type) %>%
    dplyr::mutate(label = dplyr::case_when(
      int_lib > 0 & int_exp > 0 ~ 'matched',
      !(int_lib > 0 & int_exp > 0) ~ 'unmatched'
    ))

  temp_spec_frag_nl <- obj_spec_frag_nl %>%
    dplyr::rename(mz_lib = mz,
                  mz_exp = mzExp,
                  int_lib = intensity,
                  int_exp = intensityExp) %>%
    dplyr::select(mz_lib, int_lib, mz_exp, int_exp, fragPrecursor, type) %>%
    dplyr::filter(int_lib > 0 & int_exp > 0) %>%
    dplyr::mutate(mz_lib = dplyr::case_when(type == 'exact_match' ~ mz_lib,
                                            type == 'nl_match' ~ mz_exp)) %>%
    dplyr::mutate(label = dplyr::case_when(
      int_lib > 0 & int_exp > 0 ~ 'matched',
      !(int_lib > 0 & int_exp > 0) ~ 'unmatched'
    ))

  temp_spec <- temp_spec_frag_match %>%
    dplyr::bind_rows(temp_spec_frag_nl)

  temp_spec1 <- temp_spec %>%
    dplyr::select(mz_lib, int_lib, fragPrecursor, type, label) %>%
    dplyr::mutate(int_lib = int_lib/max(int_lib)) %>%
    dplyr::rename(mz = mz_lib, int = int_lib) %>%
    dplyr::mutate(int = 0-int,
                  tag = dplyr::case_when(
                    label == 'unmatched' ~ 'frag_unmatch',
                    label == 'matched' & type == 'exact_match' ~ 'library',
                    label == 'matched' & type == 'nl_match' ~ 'library_shift'
                  ))

  temp_spec2 <- temp_spec %>%
    dplyr::select(mz_exp, int_exp, fragPrecursor, type, label) %>%
    dplyr::mutate(int_exp = int_exp/max(int_exp)) %>%
    dplyr::rename(mz = mz_exp, int = int_exp) %>%
    dplyr::mutate(tag = dplyr::case_when(
      label == 'unmatched' ~ 'frag_unmatch',
      label == 'matched' ~ 'experiment'
    ))


  temp_data <- temp_spec1 %>%
    dplyr::bind_rows(temp_spec2)

  temp_plot <- ggplot2::ggplot(temp_data) +
    ggplot2::geom_segment(ggplot2::aes(x = mz, xend = mz,
                                       y = 0, yend = int,
                                       colour = tag)) +
    ggplot2::geom_point(ggplot2::aes(x = mz,
                                     y = int,
                                     shape = label,
                                     colour = tag))+
    ggplot2::scale_colour_manual(values = c(
      'experiment' = 'dodgerblue',
      'library' = 'tomato',
      'frag_unmatch' = 'gray',
      'library_shift' = 'orange'
    )) +
    ggplot2::scale_shape_manual(values = c(
      'matched' = 16,
      'unmatched' = 4
    )) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::xlim(0.95*min(temp_data$mz),
                  1.05*max(temp_data$mz)) +
    ggplot2::ylim(-1, 1) +
    ggplot2::xlab('m/z') +
    ggplot2::ylab('Relative intensity') +
    ZZWTheme() +
    # ggplot2::ggtitle(label = paste0(obj_ms2@info$name, '{',
    #                                 round(obj_ms2@info$scoreForward, 4), '}')) +
    ggplot2::theme(legend.position = c(0.8, 0.75),
                   title = ggplot2::element_text(vjust = 0.5))


  return(temp_plot)

}


