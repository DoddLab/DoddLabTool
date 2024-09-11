#' @title match_istd
#' @description match isotope labeling internal standards (ISTD) in peak table according to the data mode
#' @author Zhiwei Zhou
#' @param obj a object of mass_dataset
#' @param mode data modes. "hilic_pos", "hilic_neg", "c18_pos", "c18_neg"
#' @param mz_ppm mz tolerance (in ppm) to match ISTD. defaule: 10 ppm
#' @param rt_second rt tolerance (in second) to match ISTD. defaule: 25 second
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' result <- match_istd(object = object,
#'   mode = 'hilic_pos',
#'   mz_ppm = 10,
#'   rt_second = 15)
#' }

# result <- match_istd(object = object,
#   mode = 'c18_pos',
#   lib_version = 'v2',
#   mz_ppm = 10,
#   rt_second = 15)

# load('~/Project/00_IBD_project/Data/20230808_data_quality_check/B0001/hilic_neg/object.RData')
# result <- match_istd(object = object,
#   mode = 'hilic_neg',
#   lib_version = 'v2',
#   mz_ppm = 10,
#   rt_second = 25)


match_istd <- function(object,
                       mode = c('hilic_pos', 'hilic_neg', 'c18_pos', 'c18_neg'),
                       lib_version = c('v2', 'v1'),
                       mz_ppm = 10,
                       rt_second = 15) {
  message(crayon::blue("Start match ..."))

  mode <- match.arg(mode)
  lib_version <- match.arg(lib_version)
  data('istd_lib', envir = environment())
  temp_istd_lib <- istd_lib[[lib_version]][[mode]]

  variable_id_istd <- mapply(function(x, y){
    matchMzRT(object = object,
              mz = x,
              rt= y,
              mz_ppm = mz_ppm,
              rt_second = rt_second)
  },
  x = temp_istd_lib$mz,
  y = temp_istd_lib$rt,
  SIMPLIFY = FALSE)

  # names(variable_id_istd) <- temp_istd_lib$name

  istd_variable <- lapply(seq_along(variable_id_istd), function(i){
    istd_mz <- temp_istd_lib$mz[i]
    istd_rt <- temp_istd_lib$rt[i]
    istd_name <- temp_istd_lib$name[i]

    temp_variable_id <- variable_id_istd[[i]]
    if (length(temp_variable_id) == 0){
      istd_variable <- data.frame(
        name_istd = istd_name,
        mz_istd = istd_mz,
        rt_istd = istd_rt,
        stringsAsFactors = FALSE) %>%
        dplyr::mutate(variable_id = NA,
                      mz = NA,
                      rt = NA,
                      mz_error = NA,
                      rt_error = NA)

      return(istd_variable)
    } else {
      temp_variable_id_info <- object@variable_info %>%
        dplyr::select(variable_id:rt) %>%
        dplyr::filter(variable_id %in% temp_variable_id)

      istd_variable <- temp_variable_id_info %>%
        dplyr::mutate(mz_istd = istd_mz,
                      rt_istd = istd_rt,
                      name_istd = istd_name) %>%
        dplyr::mutate(mz_error = round(abs(mz - mz_istd)/mz_istd * 10^6, digits = 2),
                      rt_error = round(abs(rt - rt_istd))) %>%
        dplyr::select(name_istd, mz_istd, rt_istd, dplyr::everything())

      return(istd_variable)
    }
  })


  # check which ISTD is multiple match/missing
  temp_length <- sapply(variable_id_istd, length)

  if (all(temp_length == 1)) {
    cat('All ISTD are detected!\n')
    result <- istd_variable %>%
      dplyr::bind_rows()

    return(result)
  } else {
    idx_missing <- which(temp_length < 1)
    idx_multiple <- which(temp_length > 1)

    result <- istd_variable %>%
      dplyr::bind_rows()

    if (length(idx_missing) > 0) {
      message(crayon::red(paste('Below', length(idx_missing), 'ISTD are missing:\n')))
      temp <- istd_variable[idx_missing] %>%
        dplyr::bind_rows()
      message(crayon::red(paste(unique(temp$name_istd), collapse = '; '), '\n'))
    }

    if (length(idx_multiple) > 0) {
      message(crayon::red(paste('Below', length(idx_missing), 'ISTD have multiple variables detected:\n')))
      temp <- istd_variable[idx_missing] %>%
        dplyr::bind_rows()
      message(crayon::red(paste(unique(temp$name_istd), collapse = '; '), '\n'))
    }

    return(result)
  }

}


matchMzRT <- function(object,
                      mz,
                      rt,
                      mz_ppm = 10,
                      rt_second = 15) {
  variable_info <- object@variable_info
  idx1 <- which(abs(variable_info$mz - mz)/mz*10^6 <= mz_ppm)
  idx2 <- which(abs(variable_info$rt - rt) <= rt_second)
  idx <- intersect(idx1, idx2)

  if (length(idx) == 0) {
    return(NULL)
  } else {
    temp_variable_id <- object@variable_info$variable_id[idx]
    return(temp_variable_id)
  }
}







#' @title match_biomet
#' @description match some real metabolite in peak table according to the data mode
#' @author Zhiwei Zhou
#' @param obj a object of mass_dataset
#' @param mode data modes. "hilic_pos", "hilic_neg", "c18_pos", "c18_neg"
#' @param mz_ppm mz tolerance (in ppm) to match biomet. defaule: 10 ppm
#' @param rt_second rt tolerance (in second) to match biomet. defaule: 25 second
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' result <- match_istd(object = object,
#'   mode = 'hilic_pos',
#'   mz_ppm = 10,
#'   rt_second = 15)
#' }

# result <- match_istd(object = object,
#   mode = 'c18_pos',
#   lib_version = 'v2',
#   mz_ppm = 10,
#   rt_second = 15)

# load('~/Project/00_IBD_project/Data/20230808_data_quality_check/B0001/hilic_neg/object.RData')
# result <- match_istd(object = object,
#   mode = 'hilic_neg',
#   lib_version = 'v2',
#   mz_ppm = 10,
#   rt_second = 25)


match_biomet <- function(object,
                         mode = c('hilic_pos', 'hilic_neg', 'c18_pos', 'c18_neg'),
                         lib_version = c('v1'),
                         mz_ppm = 10,
                         rt_second = 15) {
  message(crayon::blue("Start match ..."))

  mode <- match.arg(mode)
  lib_version <- match.arg(lib_version)
  data('biomet_lib', envir = environment())
  temp_biomet_lib <- biomet_lib[[lib_version]]$lib_table

  switch (mode,
    'hilic_pos' = {
      temp_biomet_lib <- temp_biomet_lib %>% dplyr::filter(polarity == 'positive' & column == 'hilic')
    },
    'hilic_neg' = {
      temp_biomet_lib <- temp_biomet_lib %>% dplyr::filter(polarity == 'negative' & column == 'hilic')
    },
    'c18_pos' = {
      temp_biomet_lib <- temp_biomet_lib %>% dplyr::filter(polarity == 'positive' & column == 'c18')
    },
    'c18_neg' = {
      temp_biomet_lib <- temp_biomet_lib %>% dplyr::filter(polarity == 'negative' & column == 'c18')
    }
  )


  variable_id_biomet <- mapply(function(x, y){
    matchMzRT(object = object,
              mz = x,
              rt= y,
              mz_ppm = mz_ppm,
              rt_second = rt_second)
  },
  x = temp_biomet_lib$mz,
  y = temp_biomet_lib$rt,
  SIMPLIFY = FALSE)

  # names(variable_id_biomet) <- temp_biomet_lib$name

  biomet_variable <- lapply(seq_along(variable_id_biomet), function(i){
    biomet_mz <- temp_biomet_lib$mz[i]
    biomet_rt <- temp_biomet_lib$rt[i]
    biomet_name <- temp_biomet_lib$name[i]

    temp_variable_id <- variable_id_biomet[[i]]
    if (length(temp_variable_id) == 0){
      biomet_variable <- data.frame(
        name_biomet = biomet_name,
        mz_biomet = biomet_mz,
        rt_biomet = biomet_rt,
        stringsAsFactors = FALSE) %>%
        dplyr::mutate(variable_id = NA,
                      mz = NA,
                      rt = NA,
                      mz_error = NA,
                      rt_error = NA)

      return(biomet_variable)
    } else {
      temp_variable_id_info <- object@variable_info %>%
        dplyr::select(variable_id:rt) %>%
        dplyr::filter(variable_id %in% temp_variable_id)

      biomet_variable <- temp_variable_id_info %>%
        dplyr::mutate(mz_biomet = biomet_mz,
                      rt_biomet = biomet_rt,
                      name_biomet = biomet_name) %>%
        dplyr::mutate(mz_error = round(abs(mz - mz_biomet)/mz_biomet * 10^6, digits = 2),
                      rt_error = round(abs(rt - rt_biomet))) %>%
        dplyr::select(name_biomet, mz_biomet, rt_biomet, dplyr::everything())

      return(biomet_variable)
    }
  })


  # check which biomet is multiple match/missing
  temp_length <- sapply(variable_id_biomet, length)

  if (all(temp_length == 1)) {
    cat('All biomet are detected!\n')
    result <- biomet_variable %>%
      dplyr::bind_rows()

    return(result)
  } else {
    idx_missing <- which(temp_length < 1)
    idx_multiple <- which(temp_length > 1)

    result <- biomet_variable %>%
      dplyr::bind_rows()

    if (length(idx_missing) > 0) {
      message(crayon::red(paste('Below', length(idx_missing), 'biomet are missing:\n')))
      temp <- biomet_variable[idx_missing] %>%
        dplyr::bind_rows()
      message(crayon::red(paste(unique(temp$name_biomet), collapse = '; '), '\n'))
    }

    if (length(idx_multiple) > 0) {
      message(crayon::red(paste('Below', length(idx_missing), 'biomet have multiple variables detected:\n')))
      temp <- biomet_variable[idx_missing] %>%
        dplyr::bind_rows()
      message(crayon::red(paste(unique(temp$name_biomet), collapse = '; '), '\n'))
    }

    return(result)
  }

}
