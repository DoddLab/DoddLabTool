################################################################################
# Search database function -----------------------------------------------------
  # get_compound -----------------------------------------------------------------
#' @title get_compound
#' @author Zhiwei Zhou
#' @description search compound information via lab_id
#' @param lab_id Dodd lab ID
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' get_compound(lab_id = "S00001")
#' }


# get_compound(lab_id = "S00001")

setGeneric(name = "get_compound",
           def = function(lab_id){
             if (length(lab_id) > 1) {
               message(crayon::red("Please search one compound per time\n"))
               stop()
             }

             data('cpd_stanford_lib', envir = environment())
             data('cpd_iroa_lib', envir = environment())

             cpd_table <- cpd_stanford_lib$compound_table %>%
               dplyr::bind_rows(cpd_iroa_lib$compound_table)
             rm(cpd_iroa_lib, cpd_stanford_lib);gc()

             idx <- which(cpd_table$lab_id == lab_id)
             if (length(idx) == 0){
               message(crayon::red("This compound is not found in DoddLib\n"))
               stop()
             }
             if (length(idx) > 1) {
               message(crayon::red("This compound have multiple items in DoddLib\n"))
               stop()
             }

             cpd_table %>%
               dplyr::slice(idx) %>%
               dplyr::mutate_all(as.character) %>%
               tidyr::pivot_longer(everything()) %>%
               print(n = Inf)
           })


  # get_ms2 ----------------------------------------------------------------------
#' @title get_ms2
#' @author Zhiwei Zhou
#' @description search compound ms2 spec via lab_id, ce, and polarity
#' @param lab_id Dodd lab ID
#' @param ce collision enery. "10", "20", "40"
#' @param polarity ionization polarity. "positive", "negative"
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' get_ms2(lab_id = 'S00001', ce = '10', polarity = 'positive')
#' get_ms2(lab_id = 'S00001', ce = '20', polarity = 'positive')
#' get_ms2(lab_id = 'S00001', ce = '40', polarity = 'positive')
#' get_ms2(lab_id = 'S00100', ce = '20', polarity = 'negative')
#' get_ms2(lab_id = 'I00101', ce = '20', polarity = 'negative')
#' get_ms2(lab_id = 'H00101', ce = '20', polarity = 'negative')
#' }

# get_ms2(lab_id = 'S00001', ce = '10', polarity = 'positive')
# get_ms2(lab_id = 'S00001', ce = '20', polarity = 'positive')
# get_ms2(lab_id = 'S00001', ce = '40', polarity = 'positive')
# get_ms2(lab_id = 'S00100', ce = '20', polarity = 'negative')
# get_ms2(lab_id = 'I00101', ce = '20', polarity = 'negative')
# get_ms2(lab_id = 'H00101', ce = '20', polarity = 'negative')
# get_ms2(lab_id = 'S00017', ce = '10', polarity = 'positive')

setGeneric(name = 'get_ms2',
           def = function(lab_id,
                          ce = c('10', '20', '40'),
                          polarity = c('positive', "negative")){

             ce <- match.arg(ce)
             polarity <- match.arg(polarity)


             data("ms2_dodd_lib", envir = environment())
             data("ms2_stanford_lib", envir = environment())
             data("ms2_iroa_lib", envir = environment())

             ids_pos <- names(ms2_dodd_lib$positive)
             ids_neg <- names(ms2_dodd_lib$positive)

             # search dodd lib first, if not, according id to search separated db
             if (polarity == 'positive' & (lab_id %in% ids_pos)) {
               temp_spec <- ms2_dodd_lib$positive[[lab_id]][[ce]]
               if (is.null(temp_spec)) {
                 message(crayon::red("No spectra of this compound & this ce is found!"))
               } else {
                 return(temp_spec)
               }
             }

             if (polarity == 'negative' & (lab_id %in% ids_neg)) {
               temp_spec <- ms2_dodd_lib$negative[[lab_id]][[ce]]
               if (is.null(temp_spec)) {
                 message(crayon::red("No spectra of this compound & this ce is found!"))
               } else {
                 return(temp_spec)
               }
             }

             if (stringr::str_detect(lab_id, pattern = "S\\d+")) {
                if (polarity == "positive") {
                  temp_spec1 <- ms2_stanford_lib$hilic_pos[[lab_id]][[ce]]
                  temp_spec2 <- ms2_stanford_lib$c18_pos[[lab_id]][[ce]]

                  # if two spec found, return the highest one
                  if (all(!is.null(temp_spec1), !is.null(temp_spec2))) {

                    if (nrow(temp_spec1) >= 10) {
                      temp_int1 <- temp_spec1[,2] %>% sort(decreasing = TRUE) %>% .[1:10] %>% sum()
                    } else {
                      temp_int1 <- temp_spec1[,2] %>% sum()
                    }

                    if (nrow(temp_spec2) >= 10) {
                      temp_int2 <- temp_spec2[,2] %>% sort(decreasing = TRUE) %>% .[1:10] %>% sum()
                    } else {
                      temp_int2 <- temp_spec2[,2] %>% sum()
                    }

                    if (temp_int1 > temp_int2) {
                      message(crayon::yellow("Two spectra of this compound are found, the highest spectrum (HILIC_Pos) is returned!"))
                      return(temp_spec1)
                    } else {
                      message(crayon::yellow("Two spectra of this compound are found, the highest spectrum (C18_Pos) is returned!"))
                      return(temp_spec2)
                    }

                  }

                  if (!is.null(temp_spec1) & is.null(temp_spec2)) {
                    return(temp_spec1)
                  }

                  if (is.null(temp_spec1) & !is.null(temp_spec2)) {
                    return(temp_spec2)
                  }

                  if (is.null(temp_spec1) & is.null(temp_spec2)) {
                    message(crayon::red("No spectra of this compound is found!"))
                  }

                } else {
                  temp_spec1 <- ms2_stanford_lib$hilic_neg[[lab_id]][[ce]]
                  temp_spec2 <- ms2_stanford_lib$c18_neg[[lab_id]][[ce]]

                  # if two spec found, return the highest one
                  if (all(!is.null(temp_spec1), !is.null(temp_spec2))) {

                    if (nrow(temp_spec1) >= 10) {
                      temp_int1 <- temp_spec1[,2] %>% sort(decreasing = TRUE) %>% .[1:10] %>% sum()
                    } else {
                      temp_int1 <- temp_spec1[,2] %>% sum()
                    }

                    if (nrow(temp_spec2) >= 10) {
                      temp_int2 <- temp_spec2[,2] %>% sort(decreasing = TRUE) %>% .[1:10] %>% sum()
                    } else {
                      temp_int2 <- temp_spec2[,2] %>% sum()
                    }

                    if (temp_int1 > temp_int2) {
                      message(crayon::yellow("Two spectra of this compound are found, the highest spectrum (HILIC_Pos) is returned!"))
                      return(temp_spec1)
                    } else {
                      message(crayon::yellow("Two spectra of this compound are found, the highest spectrum (C18_Pos) is returned!"))
                      return(temp_spec2)
                    }

                  }

                  if (!is.null(temp_spec1) & is.null(temp_spec2)) {
                    return(temp_spec1)
                  }

                  if (is.null(temp_spec1) & !is.null(temp_spec2)) {
                    return(temp_spec2)
                  }

                  if (is.null(temp_spec1) & is.null(temp_spec2)) {
                    message(crayon::red("No spectra of this compound is found!"))
                  }
                }
             } else if (stringr::str_detect(lab_id, pattern = "I\\d+")) {
               if (polarity == "positive") {
                 temp_spec1 <- ms2_iroa_lib$hilic_pos[[lab_id]][[ce]]
                 temp_spec2 <- ms2_iroa_lib$c18_pos[[lab_id]][[ce]]

                 # if two spec found, return the highest one
                 if (all(!is.null(temp_spec1), !is.null(temp_spec2))) {

                   if (nrow(temp_spec1) >= 10) {
                     temp_int1 <- temp_spec1[,2] %>% sort(decreasing = TRUE) %>% .[1:10] %>% sum()
                   } else {
                     temp_int1 <- temp_spec1[,2] %>% sum()
                   }

                   if (nrow(temp_spec2) >= 10) {
                     temp_int2 <- temp_spec2[,2] %>% sort(decreasing = TRUE) %>% .[1:10] %>% sum()
                   } else {
                     temp_int2 <- temp_spec2[,2] %>% sum()
                   }

                   if (temp_int1 > temp_int2) {
                     message(crayon::yellow("Two spectra of this compound are found, the highest spectrum (HILIC_Pos) is returned!"))
                     return(temp_spec1)
                   } else {
                     message(crayon::yellow("Two spectra of this compound are found, the highest spectrum (C18_Pos) is returned!"))
                     return(temp_spec2)
                   }

                 }

                 if (!is.null(temp_spec1) & is.null(temp_spec2)) {
                   return(temp_spec1)
                 }

                 if (is.null(temp_spec1) & !is.null(temp_spec2)) {
                   return(temp_spec2)
                 }

                 if (is.null(temp_spec1) & is.null(temp_spec2)) {
                   message(crayon::red("No spectra of this compound is found!"))
                 }

               } else {
                 temp_spec1 <- ms2_iroa_lib$hilic_neg[[lab_id]][[ce]]
                 temp_spec2 <- ms2_iroa_lib$c18_neg[[lab_id]][[ce]]

                 # if two spec found, return the highest one
                 if (all(!is.null(temp_spec1), !is.null(temp_spec2))) {

                   if (nrow(temp_spec1) >= 10) {
                     temp_int1 <- temp_spec1[,2] %>% sort(decreasing = TRUE) %>% .[1:10] %>% sum()
                   } else {
                     temp_int1 <- temp_spec1[,2] %>% sum()
                   }

                   if (nrow(temp_spec2) >= 10) {
                     temp_int2 <- temp_spec2[,2] %>% sort(decreasing = TRUE) %>% .[1:10] %>% sum()
                   } else {
                     temp_int2 <- temp_spec2[,2] %>% sum()
                   }

                   if (temp_int1 > temp_int2) {
                     message(crayon::yellow("Two spectra of this compound are found, the highest spectrum (HILIC_Pos) is returned!"))
                     return(temp_spec1)
                   } else {
                     message(crayon::yellow("Two spectra of this compound are found, the highest spectrum (C18_Pos) is returned!"))
                     return(temp_spec2)
                   }

                 }

                 if (!is.null(temp_spec1) & is.null(temp_spec2)) {
                   return(temp_spec1)
                 }

                 if (is.null(temp_spec1) & !is.null(temp_spec2)) {
                   return(temp_spec2)
                 }

                 if (is.null(temp_spec1) & is.null(temp_spec2)) {
                   message(crayon::red("No spectra of this compound is found!"))
                 }
               }
             } else {
               message(crayon::red("No spectra of this compound is found!"))
             }


           })



################################################################################
# spectral match function ------------------------------------------------------
  # convert_spectra_data ---------------------------------------------------------
#' @title convert_spectra_data
#' @param ms2_data
#' @importClassesFrom SpectraTools 'SpectraData'
#' @export

convert_spectra_data <- function(ms2_data) {
  options(readr.num_columns = 0)
  temp_info <- ms2_data$info %>%
    dplyr::rename(name = NAME,
                  mz = PRECURSORMZ) %>%
    dplyr::select(name:mz) %>%
    readr::type_convert()

  temp_ms2_data <- ms2_data$spec

  result <- new('SpectraData',
                info = temp_info,
                spectra = list(temp_ms2_data))

  return(result)
}



  # run_spec_match -------------------------------------------------------------

#' @title run_spec_match
#' @description a wrapper of runing SpectraTools
#' @author Zhiwei Zhou
#' @param obj_ms2_cpd1 experimental ms2 object. Note: info is a data.frame, spec is a list formed by matrix
#' @param obj_ms2_cpd2 library ms2 object.
#' @param mz_tol_ms2 Default: 35 ppm
#' @param scoring_approach 'dp', 'bonanza', 'hybrid', 'gnps'
#' @export
#' \dontrun{
#' cpd_info_zhulab_10v_HMDB0000448 <- data.frame(name = 'Adipic acid', mz = 147.0652, stringsAsFactors = FALSE)
#' spec_zhulab_10v_HMDB0000448 <- data.frame(mz = c(43.02108, 55.05529, 59.04936, 68.95148, 83.04965, 91.05245, 100.74727, 101.06027, 110.70124, 111.04364, 129.05302, 147.04354),
#'                                           intensity = c(0.02571532, 0.31513903, 0.08299478, 0.01847157, 0.54350417, 0.02804378, 0.01791196, 1.00000000, 0.01294796, 0.53553867, 0.12961211, 0.02379571),
#'                                           stringsAsFactors = FALSE)
#' obj_zhulab_10v_HMDB0000448 <- new('SpectraData',
#'                                   info = cpd_info_zhulab_10v_HMDB0000448,
#'                                   spectra = list(spec_zhulab_10v_HMDB0000448))
#'
#' cpd_info_dodd_10v_HMDB0000448 <- data.frame(name = 'S00017', mz = 147.0652, stringsAsFactors = FALSE)
#' spec_doddlib_10v_HMDB0000448 <- data.frame(mz = c(41.03830, 55.05456, 83.04922, 101.06048, 106.96222, 111.04383, 129.05473),
#'                                            intensity = c(209.8333, 453.0305, 172.4083, 980.7045, 151.1392, 1022.1850, 451.0853),
#'                                            stringsAsFactors = FALSE)
#' obj_doddlib_10v_HMDB0000448 <- new('SpectraData',
#'                                    info = cpd_info_doddlib_10v_HMDB0000448,
#'                                    spectra = list(spec_doddlib_10v_HMDB0000448))
#'
#' score_dp <- run_spec_match(obj_ms2_cpd1 = obj_doddlib_10v_HMDB0000448,
#'                          obj_ms2_cpd2 = obj_zhulab_10v_HMDB0000448,
#'                          mz_tol_ms2 = 35, scoring_approach = 'dp')
#' }
#

# cpd_info_zhulab_10v_HMDB0000448 <- data.frame(name = 'Adipic acid', mz = 147.0652, stringsAsFactors = FALSE)
# spec_zhulab_10v_HMDB0000448 <- data.frame(mz = c(43.02108, 55.05529, 59.04936, 68.95148, 83.04965, 91.05245, 100.74727, 101.06027, 110.70124, 111.04364, 129.05302, 147.04354),
#                                           intensity = c(0.02571532, 0.31513903, 0.08299478, 0.01847157, 0.54350417, 0.02804378, 0.01791196, 1.00000000, 0.01294796, 0.53553867, 0.12961211, 0.02379571),
#                                           stringsAsFactors = FALSE)
# obj_zhulab_10v_HMDB0000448 <- new('SpectraData',
#                                   info = cpd_info_zhulab_10v_HMDB0000448,
#                                   spectra = list(spec_zhulab_10v_HMDB0000448))
#
# cpd_info_dodd_10v_HMDB0000448 <- data.frame(name = 'S00017', mz = 147.0652, stringsAsFactors = FALSE)
# spec_doddlib_10v_HMDB0000448 <- data.frame(mz = c(41.03830, 55.05456, 83.04922, 101.06048, 106.96222, 111.04383, 129.05473),
#                                            intensity = c(209.8333, 453.0305, 172.4083, 980.7045, 151.1392, 1022.1850, 451.0853),
#                                            stringsAsFactors = FALSE)
# obj_doddlib_10v_HMDB0000448 <- new('SpectraData',
#                                    info = cpd_info_doddlib_10v_HMDB0000448,
#                                    spectra = list(spec_doddlib_10v_HMDB0000448))
#
# score_dp <- run_spec_match(obj_ms2_cpd1 = obj_doddlib_10v_HMDB0000448,
#                          obj_ms2_cpd2 = obj_zhulab_10v_HMDB0000448,
#                          mz_tol_ms2 = 35, scoring_approach = 'dp')
# plot_id_ms2(obj_ms2 = score_dp)


run_spec_match <- function(
    obj_ms2_cpd1,
    obj_ms2_cpd2,
    mz_tol_ms2 = 35,
    scoring_approach = c('dp', 'bonanza', 'hybrid', 'gnps'),
    ...
) {
  # browser()
  match.arg(scoring_approach)

  switch (scoring_approach,
          'dp' = {
            intensityNormedMethod <- 'maximum'
            methodScore <- 'dp'
          },
          'bonanza' = {
            intensityNormedMethod <- 'bonanza'
            methodScore <- 'bonanza'
          },
          'hybrid' = {
            intensityNormedMethod <- 'maximum'
            methodScore <- 'hybrid'
          },
          'gnps' = {
            intensityNormedMethod <- 'gnps'
            methodScore <- 'gnps'
          }
  )

  matchParam <- SpectraTools::MatchParam(ppm = mz_tol_ms2,
                                         cutoff = 0,
                                         weightIntensity = 1,
                                         weightMZ = 0,
                                         normIntensity = TRUE,
                                         tuneLibSpectra = TRUE,
                                         intensityExpNormed = TRUE,
                                         intensityLibNormed = TRUE,
                                         includePrecursor = TRUE,
                                         ppmPrecursorFilter = 30,
                                         thrIntensityAbs = 0,
                                         thrIntensityRel = 0,
                                         intensityNormedMethod = intensityNormedMethod,
                                         methodMatch = 'direct',
                                         methodScore = methodScore) %>%
    new(Class = 'MatchParam')


  result <- try(SpectraTools::MatchSpectra(dataExp = obj_ms2_cpd1,
                                           dataRef = obj_ms2_cpd2,
                                           matchParam),
                silent = TRUE)




  # # if the spectra has only one fragment, and it larger than precursor, it was removed
  # if (length(result) == 0) {
  #   n_frag_cpd1 <- nrow(obj_ms2_cpd1@spectra[[1]])
  #   n_frag_cpd2 <- nrow(obj_ms2_cpd2@spectra[[1]])
  #   n_frag_match <- 0
  #   n_nl_match <- 0
  #
  #   if (scoring_approach == 'dp') {
  #
  #   }
  #   result@info <- obj_ms2_cpd2@info %>%
  #     dplyr::mutate(scoreReverse = 0,
  #                   scoreForward = 0) %>%
  #     dplyr::mutate(n_frag_cpd1 = n_frag_cpd1,
  #                   n_frag_cpd2 = n_frag_cpd2,
  #                   n_frag_match = n_frag_match,
  #                   n_frag_nl = n_nl_match) %>%
  #     dplyr::mutate(n_frag_total = n_frag_match + n_frag_nl)
  # }

  # add matched_frag and matched_nl into the result table
  stat_matched_frag <- lapply(seq_along(result@matchedFragments), function(i){
    temp_matchedFragments <- result@matchedFragments[[i]]

    if (length(temp_matchedFragments) > 0) {
      n_frag_cpd1 <- temp_matchedFragments %>%
        dplyr::filter(intensity > 0) %>%
        dplyr::count() %>%
        dplyr::pull()

      n_frag_cpd2 <- temp_matchedFragments %>%
        dplyr::filter(intensityExp > 0) %>%
        dplyr::count() %>%
        dplyr::pull()

      n_frag_match <- temp_matchedFragments %>%
        dplyr::filter(intensity > 0 & intensityExp > 0) %>%
        dplyr::count() %>%
        dplyr::pull()

    } else {
      n_frag_cpd1 <- nrow(obj_ms2_cpd1@spectra[[1]])
      n_frag_cpd2 <- nrow(obj_ms2_cpd2@spectra[[i]])
      n_frag_match <- 0
    }

    if (scoring_approach == 'dp') {
      n_nl_match <- 0
    } else {
      n_nl_match <- result@nlFragments[[1]] %>%
        dplyr::filter(intensity > 0 & intensityExp > 0) %>%
        dplyr::count() %>%
        dplyr::pull()
    }

    temp_result <- tibble::tibble(n_frag_cpd1 = n_frag_cpd1,
                                  n_frag_cpd2 = n_frag_cpd2,
                                  n_frag_match = n_frag_match,
                                  n_frag_nl = n_nl_match) %>%
      dplyr::mutate(n_frag_total = n_frag_match + n_frag_nl)

  })

  stat_matched_frag <- stat_matched_frag %>% dplyr::bind_rows()

  result@info <- result@info %>%
    dplyr::bind_cols(stat_matched_frag)

  return(result)
}


