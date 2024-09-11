################################################################################
# Metabolite identification Main function --------------------------------------
  # annotate_metabolite --------------------------------------------------------

#' @title annotate_metabolite
#' @description Annotate peak using MS/MS spectra in in-house database.
#' @author Zhiwei Zhou
#' @param ms1_file The name of ms1 peak table. Column 1 is "name", Column 2 is "mz" and column is "rt".
#' @param ms2_file Default: NULL
# #' @param sample_info_file Default: sample.info.csv
#' @param ms2_type "mgf", "mzXML", "msp", "cef"
#' @param metdna_version 'version1', 'version2'. Default: "version1"
#' @param path Default: '.'
#' @param instrument The instrument you used to acquire data. "AgilentQTOF", "SciexTripleTOF", "BrukerQTOF", "ThermoOrbitrap", "ThermoExploris". Default: "SciexTripleTOF"
#' @param lib Default: 'zhuMetLib'
#' @param column "hilic", "rp"
#' @param ce "10", "20", "30", "35,15", "40", "50". Default: "30"
#' @param method_lc 'Amide12min', 'Amide23min'. Default: 'Amide12min'
#' @param excluded_adduct adduct list for exclusion. Default: NULL
#' @param is_rt_calibration Default: FALSE
#' @param mz_tol ms1 match. Default: 25 ppm
#' @param pf_rt_tol penalty-free rt range. Default: 0 s
#' @param tolerance_rt_range maxmium rt tolerance. Default: 30 s
#' @param pf_ccs_tol penalty-free rt range. Default: 0 %
#' @param tolerance_ccs_range maxmium rt tolerance. Default: 2 %
#' @param is_filter Whether filter candidates whose rt/ccs exceed the maximum tolerance. Default: FALSE
#' @param is_rt_score logical vector. Default: TRUE
#' @param is_ccs_score logical vector. Default: FALSE
#' @param is_ms2_score logical vector. Default: TRUE
#' @param is_include_precursor Default: TRUE
#' @param int_ms2_min_abs Default: 30
#' @param int_ms2_min_relative Default: 0.01
#' @param mz_tol_combine_ms1_ms2 Default: 25 ppm
#' @param rt_tol_combine_ms1_ms2 Default: 10 s
#' @param ccs_tol_combine_ms1_ms2 Only appliable for IM-MS data. Default: NULL
#' @param mz_tol_ms2 Default: 35 ppm
#' @param dp_cutoff The cutoff of dot product. Default is 0.8.
#' @param matched_frag_cutoff The cutoff of matched fragment number. Default: 1
#' @param direction Reserve which direct dot product. 'reverse', 'forward'. Default: 'reverse'
#' @param scoring_approach spectral match approach, including 'dp', 'bonanza', 'hybrid', 'gnps'. Default: 'dp'.
#' @param is_plot_ms2 Output MS2 match plot or not. Default: TRUE
#' @return Return the annotation result.
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @importClassesFrom SpectraTools 'SpectraData'
#' @export




# load('~/Project/00_IBD_project/Data/20230207_data_cleaning_IBD/HILIC_pos/02_data_cleaning/05_object_hilic_pos_outlier_removal.RData')
# object <- object_hilic_pos;rm(object_hilic_pos);gc()
# path <- '~/Project/00_IBD_project/Data/20230227_develop_metabolite_ID_workflow/'
#
# object = NULL
# ms1_file = "data.csv"
# ms2_file = NULL
# ms2_type = c('mzML')
#
# # parameters of loadDB
# polarity = c("positive")
# lib = c('dodd')
# column = c("hilic")
# ce = "20"
# adduct_list = c('[M+H]+')
#
# mz_tol = 15
# mz_ppm_thr = 150
# pf_rt_range = 0
# tolerance_rt_range = 30
# is_rt_score = TRUE
# is_ms2_score = TRUE
#
# is_include_precursor = TRUE
# int_ms2_min_abs = 50
# int_ms2_min_relative = 0.01
# mz_tol_combine_ms1_ms2 = 25 # ppm
# rt_tol_combine_ms1_ms2 = 10 # s
# mz_tol_ms2 = 35
# dp_cutoff = 0.8
# matched_frag_cutoff = 1
# direction = c('forward')
# scoring_approach = c('dp')


setGeneric(name = "annotate_metabolite",
           def = function(object = NULL,
                          ms1_file = "data.csv",
                          ms2_file = NULL,
                          ms2_type = c('mzML', "mzXML", "mgf", "msp"),
                          path = ".",

                          # parameters of loadDB
                          lib_db,
                          polarity = c("positive", "negative"),
                          lib = c('dodd'),
                          column = c("hilic", "c18"),
                          ce = c("20", "10", "40"),
                          adduct_list = c('[M+H]+'),

                          # parameters of matchMs1WithSpecLib
                          mz_tol = 15,
                          mz_ppm_thr = 150,
                          pf_rt_range = 10,
                          tolerance_rt_range = 20,
                          is_rt_score = TRUE,
                          is_ms2_score = TRUE,

                          # parameters of matchMs2WithSpecLib
                          is_include_precursor = TRUE,
                          int_ms2_min_abs = 50,
                          int_ms2_min_relative = 0.01,
                          mz_tol_combine_ms1_ms2 = 20, # ppm
                          rt_tol_combine_ms1_ms2 = 20, # s
                          mz_tol_ms2 = 35,
                          dp_cutoff = 0.8,
                          matched_frag_cutoff = 1,
                          direction = c('reverse', 'forward'),
                          scoring_approach = c('dp', 'bonanza', 'hybrid', 'gnps'),

                          # plots
                          is_plot_ms2 = TRUE
           ){
             polarity <- match.arg(polarity)
             column <- match.arg(column)
             ms2_type <- match.arg(ms2_type)
             ce <- match.arg(ce)
             direction <- match.arg(direction)
             scoring_approach <- match.arg(scoring_approach)

             path_output <- file.path(path, "01_metabolite_annotation")
             dir.create(file.path(path_output, "00_intermediate_data"), showWarnings = FALSE, recursive = TRUE)

             # check ms1_file and ms2_file
             if (is.null(object)) {
               temp <- list.files(path = path, pattern = ms1_file, recursive = TRUE)
               if (length(temp) == 0) {stop("There is no", ms1_file, '\n')}
               temp <- list.files(path = path, pattern = paste0("\\.", ms2_type), recursive = TRUE)
               if (length(temp) == 0) {stop("There is no", ms2_type, '\n')}
             }

             # get ms2 file names
             if (length(ms2_file) == 0) {
               ms2_file <- list.files(path = path, pattern = paste0("\\.", ms2_type), recursive = TRUE)
             }

             message(crayon::blue('Load database ...\n'))

             if (missing(lib_db)) {
               lib_db <- load_spec_db(lib = lib,
                                      column = column,
                                      ce = ce,
                                      polarity = polarity,
                                      adduct_list = adduct_list)
             }

             lib_meta <- lib_db$lib_meta
             lib_spec <- lib_db$lib_spec

             save(lib_meta,
                  file = file.path(path_output, "00_intermediate_data", 'lib_meta'),
                  compress = 'gzip',
                  version = 2)

             save(lib_spec,
                  file = file.path(path_output, "00_intermediate_data", 'lib_spec'),
                  compress = 'gzip',
                  version = 2)

             rm(lib_db);gc()

             # generate ms1 data
             if (is.null(object)) {
               ms1_data <- readMs1(filename = ms1_file,
                                   instrument = instrument,
                                   path = path)
             } else {
               variable_data <- object@variable_info %>%
                 dplyr::select(variable_id:rt) %>%
                 dplyr::rename(name = variable_id)
               expression_profile_data <- object@expression_data
               ms1_data <- list(info = variable_data, subject = expression_profile_data)
               rm(variable_data, expression_profile_data);gc()
             }

             save(ms1_data,
                  file = file.path(path_output, "00_intermediate_data", 'ms1_data'),
                  compress = 'gzip',
                  version = 2)

             # mz, rt match for ms1 data
             message(crayon::blue('MS1 & RT match...\n'))
             ms1_result <- match_ms1_rt(ms1_data = ms1_data,
                                        lib_meta = lib_meta,
                                        mz_tol = mz_tol,
                                        mz_ppm_thr = mz_ppm_thr,
                                        pf_rt_range = pf_rt_range,
                                        tolerance_rt_range = tolerance_rt_range,
                                        is_rt_score = is_rt_score)

             save(ms1_result,
                  file = file.path(path_output, "00_intermediate_data", 'ms1_result'),
                  compress = 'gzip',
                  version = 2)

             # if is_ms2_score
             #    check input object whether contain ms2
             #        true: extract ms2 from obj data, and save as ms2
             #        false:
             #          read & clean ms2 data
             #          select ms2 data
             #          if object existed, add into object

             if (is_ms2_score) {
               cat("\n")

               # process or load ms2
               if ('ms2' %in% list.files(file.path(path_output, "00_intermediate_data"))) {
                 cat('Read, purify ms2 spec, and combine with ms1 & ms2...\n')
                 cat('Note: load the existed ms2\n\n')
                 load(file.path(path_output, "00_intermediate_data", 'ms2'))

                 # adjust has_ms2 slot in the class
                 ms1_result <- add_has_ms2_to_SpecAnnotationClass(obj_class = ms1_result, ms2 = ms2)
               } else {
                 if (!is.null(object) & length(object@ms2_data) > 0) {
                   info <- data.frame(mz = object@ms2_data[[1]]@ms2_mz,
                                      rt = object@ms2_data[[1]]@ms2_rt,
                                      filename = object@ms2_data[[1]]@ms2_file,
                                      stringsAsFactors = FALSE)
                   spec <- object@ms2_data[[1]]@ms2_spectra
                   ms2_data_combined <- list(info = info, spec = spec)
                   rm(info, spec);gc()

                   save(ms2_data_combined,
                        file = file.path(path_output, "00_intermediate_data", 'ms2_data_combined'),
                        compress = 'gzip',
                        version = 2)

                   # generate ms2 for ms2 match
                   ms2 <- split_ms2(ms2_data_combined)
                   save(ms2, file = file.path(path_output, "00_intermediate_data", 'ms2'))

                   # adjust has_ms2 slot in the class
                   ms1_result <- add_has_ms2_to_SpecAnnotationClass(obj_class = ms1_result, ms2 = ms2)

                   # export ms2 as MSP file
                   purrr::walk(ms2, function(x){
                     generate_msp(file_name = file.path(path, 'ms2.msp'),
                                  cmp_name = x$info[1],
                                  precusormz = x$info[2],
                                  rt = x$info[3],
                                  polarity = polarity,
                                  spec = x$spec)
                   })

                 } else {
                   ms2_data <- read_ms2(ms2_file = file.path(path, ms2_file),
                                        ms2_type = ms2_type)
                   ms2_data <- integrate_ms2(ms2_data = ms2_data,
                                             ms2_file = file.path(path, ms2_file),
                                             ms2_type = ms2_type,
                                             is_include_precursor = is_include_precursor,
                                             is_deisotope = FALSE,
                                             int_ms2_min_abs = int_ms2_min_abs,
                                             int_ms2_min_relative = int_ms2_min_relative)
                   save(ms2_data,
                        file = file.path(path_output, "00_intermediate_data", 'ms2_data'),
                        compress = 'gzip',
                        version = 2)

                   ms2_data_combined <- combine_ms1_ms2(ms1_data = ms1_data,
                                                        ms2_data = ms2_data,
                                                        ms2_type = ms2_type,
                                                        mz_tol_combine_ms1_ms2 = mz_tol_combine_ms1_ms2,
                                                        rt_tol_combine_ms1_ms2 = rt_tol_combine_ms1_ms2)

                   save(ms2_data_combined,
                        file = file.path(path_output, "00_intermediate_data", 'ms2_data_combined'),
                        compress = 'gzip',
                        version = 2)

                   # generate ms2 for ms2 match
                   ms2 <- split_ms2(ms2_data_combined)
                   save(ms2, file = file.path(path_output, "00_intermediate_data", 'ms2'))

                   # adjust has_ms2 slot in the class
                   ms1_result <- add_has_ms2_to_SpecAnnotationClass(obj_class = ms1_result, ms2 = ms2)

                   # export ms2 as MSP file
                   purrr::walk(ms2, function(x){
                     generate_msp(file_name = file.path(path, 'ms2.msp'),
                                  cmp_name = x$info[1],
                                  precusormz = x$info[2],
                                  rt = x$info[3],
                                  polarity = polarity,
                                  spec = x$spec)
                   })

                   # if object existed, add ms2_data to tidymass object
                   if (!is.null(object) & length(object@ms2_data) < 1) {
                     # convert to ms2_data class for tidymass
                     temp_ms2_data <-
                       new(
                         Class = "ms2_data",
                         column = column,
                         polarity = polarity,
                         variable_id = names(ms2_data_combined$spec),
                         ms2_spectrum_id = names(ms2_data_combined$spec),
                         ms2_mz = ms2_data_combined$info$mz,
                         ms2_rt = ms2_data_combined$info$rt,
                         ms2_file = ms2_data_combined$info$filename,
                         ms2_spectra = ms2_data_combined$spec,
                         mz_tol = mz_tol_combine_ms1_ms2,
                         rt_tol = rt_tol_combine_ms1_ms2
                       )

                     temp_ms2_data = list(name = temp_ms2_data)
                     temp_name <- paste(sort(basename(ms2_file)), collapse = ";")
                     names(temp_ms2_data) <- temp_name
                     object@ms2_data <- temp_ms2_data
                     save(object, file = file.path(path_output, "00_intermediate_data", 'object_with_ms2'))
                     rm(temp_ms2_data, temp_name);gc()
                   }
                 }
               }

               # calculate ms2 score
               if (!('result_annotation' %in% list.files(file.path(path_output, "00_intermediate_data")))) {
                 cat('\n');cat('start ms2 matching...\n')
                 # ms2 order and name matched with ms1
                 ms2_name <- unname(unlist(lapply(ms2, function(x) x[[1]][1,])))
                 ms1_name <- ms1_data$info$name
                 temp_idx <- match(ms2_name, ms1_name)
                 exp_spec <- vector(mode = "list", length = nrow(ms1_data$info))
                 names(exp_spec) <- ms1_name
                 exp_spec[temp_idx] <- ms2

                 exp_spec <- lapply(exp_spec, function(x){
                   if(is.null(x)) return(x)
                   x[[2]]
                 })

                 result_annotation <- match_ms2(ms1_result = ms1_result,
                                                exp_spec = exp_spec,
                                                lib_meta = lib_meta,
                                                lib_spec = lib_spec,
                                                mz_tol_ms2 = mz_tol_ms2,
                                                dp_cutoff = dp_cutoff,
                                                matched_frag_cutoff = matched_frag_cutoff,
                                                direction = direction,
                                                scoring_approach = scoring_approach,
                                                path = path_output,
                                                is_include_precursor = is_include_precursor)

                 save(result_annotation,
                      file = file.path(path_output, '00_intermediate_data', 'result_annotation'),
                      version = 2)
               } else {
                 cat('start ms2 matching...\n')
                 cat('Note: load the existed result_annotation\n')

                 load(file.path(path_output, '00_intermediate_data', 'result_annotation'))
               }

             } else {
               result_annotation <- ms1_result

               save(result_annotation,
                    file = file.path(path_output, '00_intermediate_data', 'result_annotation'),
                    version = 2)
             }

             # generate annotation table for export
             table_annotation <- convertSpecAnnotationClass2Table(ms1_data = ms1_data,
                                                                  result_annotation = result_annotation)

             readr::write_csv(table_annotation,
                              file.path(path_output, "ms2_match_annotation_result.csv"))

             save(table_annotation,
                  file = file.path(path_output, '00_intermediate_data', 'table_annotation'),
                  version = 2)

             # generate annotation summary table: contain more details for annotation
             annot_table <- generate_summary_table(result_annotation = result_annotation)
             writexl::write_xlsx(annot_table,
                                 path = file.path(path_output, "annotation_summary.xlsx"),
                                 format_headers = FALSE)

             save(annot_table,
                  file = file.path(path_output, '00_intermediate_data', 'annot_table'),
                  version = 2)


             # ms2 plots -------------------------------------------------------
             is_ms2_candidates <- any(nrow(annot_table) > 0)

             if (is_plot_ms2 & is_ms2_score & is_ms2_candidates) {
               load(file.path(path_output, '00_intermediate_data', 'ms2_result'))

               cat('\n')

               if (scoring_approach == 'dp') {
                 purrr::walk(c('forward', 'reverse'), function(x){
                   cat('Plot', x, 'spec match...\n')

                   if (x == 'forward') {
                     temp_result <- annot_table %>%
                       dplyr::rename(ms2_score = msms_score_forward)
                   } else {
                     temp_result <- annot_table %>%
                       dplyr::rename(ms2_score = msms_score_reverse)
                   }

                   plot_list <- pbapply::pblapply(seq_along(temp_result$name), function(i){
                     # cat(i, ' ')
                     temp_feature_name <- temp_result$feature_name[i]

                     cpd_id <- temp_result$id[i]
                     cpd_name <- temp_result$name[i]
                     cpd_score <- temp_result$ms2_score[i]

                     idx <- which(names(ms2_result) == temp_feature_name)
                     temp_ms2_obj <- ms2_result[[idx]]

                     # need to be fix
                     idx_id <- match(cpd_id, temp_ms2_obj@info$name)

                     plot_list <- purrr::map(seq_along(idx_id), function(j){
                       temp_idx <- idx_id[j]
                       suppressMessages(
                         temp_plot <- plot_id_ms2(obj_spec = temp_ms2_obj@matchedFragments[[temp_idx]]) +
                           ggplot2::scale_colour_manual(
                             name = 'Attribute',
                             labels= c(paste0('Experiment ', '(', temp_feature_name, ')'),
                                       'Unmatched fragments',
                                       paste0('Library ', '(', cpd_id, ')')),
                             values = c(
                               'experiment' = 'black',
                               'library' = 'red',
                               'frag_unmatch' = 'gray'
                             )
                           ) +
                           ggplot2::scale_shape_manual(
                             name = 'Label',
                             labels= c('matched' = "Matched",
                                       'unmatched' = "Unmatched"),
                             values = c(
                               'matched' = 16,
                               'unmatched' = 4
                             )
                           ) +
                           ggplot2::ggtitle(label = paste0(temp_feature_name,
                                                           ': ', cpd_name[j],
                                                           ' (ID: ', cpd_id[j],
                                                           ', DP: ', cpd_score[j],
                                                           ')')) +
                           ggplot2::theme(legend.position = c(0.85, 0.85),
                                          title = ggplot2::element_text(vjust = 0.5))
                       )
                       return(temp_plot)
                     })

                     return(plot_list)

                   })

                   # export plot list as one pdf
                   plot_list <- do.call(c, plot_list)
                   dir.create(file.path(path_output, '02_experimental_ms2_spec_plot'),
                              showWarnings = FALSE, recursive = TRUE)
                   ggplot2::ggsave(
                     filename = file.path(path_output,
                                          '02_experimental_ms2_spec_plot',
                                          paste0('ms2_match_plot_', x, '.pdf')),
                     plot = gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1),
                     width = 15, height = 9
                   )

                   cat('\n\n')
                 })


                 # purrr::walk(c('forward', 'reverse'), function(x){
                 #   cat('Plot', x, 'spec match...\n')
                 #
                 #   if (x == 'forward') {
                 #     temp_result <- table_annotation %>%
                 #       dplyr::filter(!is.na(id_forward_summary)) %>%
                 #       dplyr::rename(id = id_forward_summary)
                 #   } else {
                 #     temp_result <- table_annotation %>%
                 #       dplyr::filter(!is.na(id_reverse_summary)) %>%
                 #       dplyr::rename(id = id_reverse_summary)
                 #   }
                 #
                 #   plot_list <- purrr::map(seq_along(temp_result$name), function(i){
                 #     cat(i, ' ')
                 #     temp_feature_name <- temp_result$name[i]
                 #
                 #     temp_id <- temp_result %>%
                 #       dplyr::filter(name == temp_feature_name) %>%
                 #       dplyr::select(name:rt, id) %>%
                 #       tidyr::separate_rows(id, sep = '\\};') %>%
                 #       dplyr::pull(id)
                 #     cpd_id <- temp_id %>%
                 #       stringr::str_extract(pattern = 'labid\\{[A-Za-z0-9]+\\}') %>%
                 #       gsub(pattern = 'labid\\{', replacement = '') %>%
                 #       gsub(pattern = '\\}', replacement = '')
                 #     cpd_name <- temp_id %>%
                 #       stringr::str_extract(pattern = "name\\{[^\\{]+\\}") %>%
                 #       gsub(pattern = 'name\\{', replacement = '') %>%
                 #       gsub(pattern = '\\}', replacement = '')
                 #     cpd_score <- temp_id %>%
                 #       stringr::str_extract(pattern = 'score\\{0\\.[0-9]+\\}|score\\{1\\}') %>%
                 #       gsub(pattern = 'score\\{', replacement = '') %>%
                 #       gsub(pattern = '\\}', replacement = '')
                 #
                 #     idx <- which(names(ms2_result) == temp_feature_name)
                 #     temp_ms2_obj <- ms2_result[[idx]]
                 #
                 #     # need to be fix
                 #     idx_id <- match(cpd_id, temp_ms2_obj@info$name)
                 #
                 #     plot_list <- purrr::map(seq_along(idx_id), function(j){
                 #       temp_idx <- idx_id[j]
                 #       suppressMessages(
                 #         temp_plot <- plot_id_ms2(obj_spec = temp_ms2_obj@matchedFragments[[temp_idx]]) +
                 #           ggplot2::scale_colour_manual(
                 #             name = 'Attribute',
                 #             labels= c(paste0('Experiment ', '(', temp_feature_name, ')'),
                 #                       'Unmatched fragments',
                 #                       paste0('Library ', '(', cpd_id, ')')),
                 #             values = c(
                 #               'experiment' = 'black',
                 #               'library' = 'red',
                 #               'frag_unmatch' = 'gray'
                 #             )
                 #           ) +
                 #           ggplot2::scale_shape_manual(
                 #             name = 'Label',
                 #             labels= c('matched' = "Matched",
                 #                       'unmatched' = "Unmatched"),
                 #             values = c(
                 #               'matched' = 16,
                 #               'unmatched' = 4
                 #             )
                 #           ) +
                 #           ggplot2::ggtitle(label = paste0(temp_feature_name,
                 #                                           ': ', cpd_name[j],
                 #                                           ' (ID: ', cpd_id[j],
                 #                                           ', DP: ', cpd_score[j],
                 #                                           ')')) +
                 #           ggplot2::theme(legend.position = c(0.85, 0.85),
                 #                          title = ggplot2::element_text(vjust = 0.5))
                 #       )
                 #       return(temp_plot)
                 #     })
                 #
                 #     return(plot_list)
                 #
                 #   })
                 #
                 #   # export plot list as one pdf
                 #   plot_list <- do.call(c, plot_list)
                 #   dir.create(file.path(path_output, '02_experimental_ms2_spec_plot'),
                 #              showWarnings = FALSE, recursive = TRUE)
                 #   ggplot2::ggsave(
                 #     filename = file.path(path_output,
                 #                          '02_experimental_ms2_spec_plot',
                 #                          paste0('ms2_match_plot_', x, '.pdf')),
                 #     plot = gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1),
                 #     width = 15, height = 9
                 #   )
                 #
                 #   cat('\n\n')
                 # })
               } else {
                 cat('Plot shifted spec match...\n')

                 temp_result <- annot_table %>%
                   dplyr::rename(ms2_score = id_reverse_summary)

                 plot_list <- purrr::map(seq_along(temp_result$name), function(i){

                   temp_feature_name <- temp_result$feature_name[i]

                   cpd_id <- temp_result$id[i]
                   cpd_name <- temp_result$name[i]
                   cpd_score <- temp_result$ms2_score[i]

                   idx <- which(names(ms2_result) == temp_feature_name)
                   temp_ms2_obj <- ms2_result[[idx]]

                   idx_id <- match(cpd_id, temp_ms2_obj@info$name)

                   plot_list <- purrr::map(seq_along(idx_id), function(j){
                     temp_idx <- idx_id[j]
                     suppressMessages(
                       temp_plot <- plot_id_shift_ms2(obj_spec_frag_match = temp_ms2_obj@matchedFragments[[temp_idx]],
                                                      obj_spec_frag_nl = temp_ms2_obj@nlFragments[[temp_idx]]) +
                         ggplot2::scale_colour_manual(
                           name = 'Attribute',
                           labels= c(paste0('Experiment ', '(', temp_feature_name, ')'),
                                     'Unmatched fragments',
                                     paste0('Library ', '(', temp_ms2_obj@info$name[temp_idx], ')'),
                                     paste0('Library shift', '(', temp_ms2_obj@info$name[temp_idx], ')')),
                           values = c(
                             'experiment' = 'black',
                             'library' = 'red',
                             'frag_unmatch' = 'gray',
                             'library_shift' = 'blue'
                           )
                         ) +
                         ggplot2::scale_shape_manual(
                           name = 'Label',
                           labels= c('matched' = "Matched",
                                     'unmatched' = "Unmatched"),
                           values = c(
                             'matched' = 16,
                             'unmatched' = 4
                           )
                         ) +
                         ggplot2::ggtitle(label = paste0(temp_feature_name,
                                                         ': ', cpd_name[j],
                                                         ' (ID: ', cpd_id[j],
                                                         ', DP: ', cpd_score[j],
                                                         ')')) +
                         ggplot2::theme(legend.position = c(0.85, 0.85),
                                        title = ggplot2::element_text(vjust = 0.5))

                     )

                     return(temp_plot)
                   })

                   return(plot_list)
                 })

                 # export plot list as one pdf
                 plot_list <- do.call(c, plot_list)
                 dir.create(file.path(path_output, '02_experimental_ms2_spec_plot'),
                            showWarnings = FALSE, recursive = TRUE)
                 ggplot2::ggsave(
                   filename = file.path(path_output,
                                        '02_experimental_ms2_spec_plot',
                                        paste0('ms2_match_plot_', scoring_approach, '.pdf')),
                   plot = gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1),
                   width = 15, height = 9
                 )

                 # temp_result <- table_annotation %>%
                 #   dplyr::filter(!is.na(id_reverse_summary)) %>%
                 #   dplyr::rename(id = id_reverse_summary)
                 #
                 # plot_list <- purrr::map(seq_along(temp_result$name), function(i){
                 #   temp_feature_name <- temp_result$name[i]
                 #   temp_id <- temp_result %>%
                 #     dplyr::filter(name == temp_feature_name) %>%
                 #     dplyr::select(name:rt, id) %>%
                 #     tidyr::separate_rows(id, sep = ';score') %>%
                 #     dplyr::pull(id)
                 #
                 #   cpd_id <- temp_id %>%
                 #     stringr::str_extract(pattern = 'labid\\{[A-Za-z0-9]+\\}') %>%
                 #     gsub(pattern = 'labid\\{', replacement = '') %>%
                 #     gsub(pattern = '\\}', replacement = '')
                 #   cpd_name <- temp_id %>%
                 #     stringr::str_extract(pattern = "name\\{[^\\{]+\\}") %>%
                 #     gsub(pattern = 'name\\{', replacement = '') %>%
                 #     gsub(pattern = '\\}', replacement = '')
                 #   cpd_score <- temp_id %>%
                 #     stringr::str_extract(pattern = 'score\\{0\\.[0-9]+\\}|score\\{1\\}') %>%
                 #     gsub(pattern = 'score\\{', replacement = '') %>%
                 #     gsub(pattern = '\\}', replacement = '')
                 #
                 #   idx <- which(names(ms2_result) == temp_feature_name)
                 #   temp_ms2_obj <- ms2_result[[idx]]
                 #
                 #   idx_id <- match(cpd_id, temp_ms2_obj@info$name)
                 #
                 #   plot_list <- purrr::map(seq_along(idx_id), function(j){
                 #     temp_idx <- idx_id[j]
                 #     suppressMessages(
                 #       temp_plot <- plot_id_shift_ms2(obj_spec_frag_match = temp_ms2_obj@matchedFragments[[temp_idx]],
                 #                                      obj_spec_frag_nl = temp_ms2_obj@nlFragments[[temp_idx]]) +
                 #         ggplot2::scale_colour_manual(
                 #           name = 'Attribute',
                 #           labels= c(paste0('Experiment ', '(', temp_feature_name, ')'),
                 #                     'Unmatched fragments',
                 #                     paste0('Library ', '(', temp_ms2_obj@info$name[temp_idx], ')'),
                 #                     paste0('Library shift', '(', temp_ms2_obj@info$name[temp_idx], ')')),
                 #           values = c(
                 #             'experiment' = 'black',
                 #             'library' = 'red',
                 #             'frag_unmatch' = 'gray',
                 #             'library_shift' = 'blue'
                 #           )
                 #         ) +
                 #         ggplot2::scale_shape_manual(
                 #           name = 'Label',
                 #           labels= c('matched' = "Matched",
                 #                     'unmatched' = "Unmatched"),
                 #           values = c(
                 #             'matched' = 16,
                 #             'unmatched' = 4
                 #           )
                 #         ) +
                 #         ggplot2::ggtitle(label = paste0(temp_feature_name,
                 #                                         ': ', cpd_name[j],
                 #                                         ' (ID: ', cpd_id[j],
                 #                                         ', DP: ', cpd_score[j],
                 #                                         ')')) +
                 #         ggplot2::theme(legend.position = c(0.85, 0.85),
                 #                        title = ggplot2::element_text(vjust = 0.5))
                 #
                 #     )
                 #
                 #     return(temp_plot)
                 #   })
                 #
                 #   return(plot_list)
                 # })
                 #
                 # # export plot list as one pdf
                 # plot_list <- do.call(c, plot_list)
                 # dir.create(file.path(path_output, '02_experimental_ms2_spec_plot'),
                 #            showWarnings = FALSE, recursive = TRUE)
                 # ggplot2::ggsave(
                 #   filename = file.path(path_output,
                 #                        '02_experimental_ms2_spec_plot',
                 #                        paste0('ms2_match_plot_', scoring_approach, '.pdf')),
                 #   plot = gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1),
                 #   width = 15, height = 9
                 # )

                 cat('\n\n')

               }

             }
             if (is_plot_ms2 & is_ms2_score & (!is_ms2_candidates)) {message(crayon::red('No ms2 annotated ...\n'))}

             message(crayon::blue('Job done!\n'))
           })



################################################################################
# load_spec_db -----------------------------------------------------------------

#' @title load_spec_db
#' @author Zhiwei Zhou
#' @param lib database name, 'zhuMetLib', 'zhuMetLib_orbitrap', 'fiehnHilicLib'. Default: 'zhuMetLib'
#' @param instrument "SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap", "ThermoExploris", "AgilentDTIMMS", "BrukerTIMS". Default: "SciexTripleTOF"
#' @param column 'hilic', 'rp'. Default: 'hilic'
#' @param method_lc 'Amide12min' or 'Amide23min'
#' @param ce "10", "20", "30", "35,15", "40", "50"; Default: '30'
#' @param polarity 'positive' or 'negative'. Default: 'positive'
#' @param adduct_list NULL
# #' @param is_rt_score whether only reserve compounds with RT
# #' @param is_ccs_score whether only reserve compounds with CCS
#' @param is_rt_calibration TRUE
#' @param path '.'
#' @export

# test <- load_spec_db(lib = 'dodd', column = 'hilic', ce = '20', polarity = 'positive', adduct_list = '[M+H]+')

setGeneric(name = 'load_spec_db',
           def = function(
    lib = c('dodd'),
    column = c('hilic', 'c18'),
    ce = c('10', '20', '40'),
    polarity = c('positive', 'negative'),
    adduct_list = NULL
    # file_rt_ref = 'RT_recalibration_table.csv',
    # is_rt_calibration = FALSE,
    # path = NULL
           ){
             lib <- match.arg(lib)
             column <- match.arg(column)
             ce <- match.arg(ce)
             polarity <- match.arg(polarity)

             # check adduct type
             if (length(adduct_list) == 0) {
               stop('Please input adduct_list\n')
             }

             if (polarity == 'positive') {
               temp <- all(adduct_list %in% lib_adduct_nl$positive$adduct)
               if (!temp) stop('Please check the selected adducts!\n')
             } else {
               temp <- all(adduct_list %in% lib_adduct_nl$negative$adduct)
               if(!temp) stop('Please check the selected adducts!\n')
             }

             switch (lib,
                     'dodd' = {
                       data('cpd_dodd_lib', envir = environment())
                       data('ms2_dodd_lib', envir = environment())
                       cpd_lib <- cpd_dodd_lib
                       ms2_lib <- ms2_dodd_lib
                       rm(cpd_dodd_lib, ms2_dodd_lib);gc()
                     }
             )

             # modify lib_meta
             lib_meta <- cpd_lib$compound_table
             switch(polarity,
                    'positive' = {
                      temp_adduct_table <- lib_adduct_nl$positive %>% dplyr::filter(adduct %in% adduct_list)
                    },
                    'negative' = {
                      temp_adduct_table <- lib_adduct_nl$negative %>% dplyr::filter(adduct %in% adduct_list)
                    }
             )

             # calculate mz
             lib_mz <- lapply(seq_along(lib_meta$lab_id), function(i){
               temp_mz <- lib_meta$monoisotopic_mass[i]

               result <- sapply(seq_along(temp_adduct_table$adduct), function(j){
                 calculateMz(exact_mass = temp_mz,
                             adduct = temp_adduct_table$adduct[j],
                             delta_mz = temp_adduct_table$delta_mz[j])
               })
               result
             })
             lib_mz <- lib_mz %>% do.call(rbind, .) %>% tibble::as_tibble()
             colnames(lib_mz) <- temp_adduct_table$adduct

             # extract compound info
             lib_meta <- lib_meta %>%
               dplyr::bind_cols(lib_mz) %>%
               tidyr::pivot_longer(cols = temp_adduct_table$adduct,
                                   names_to = 'adduct',
                                   values_to = 'mz') %>%
               dplyr::rename(id = lab_id, name = compound_name) %>%
               dplyr::select(id:formula, adduct:mz, dplyr::everything())

             # select mode responding RT
             if (column == 'hilic') {
               if (polarity == 'positive') {
                 lib_meta <- lib_meta %>%
                   dplyr::mutate(rt = rt_hilic_pos) %>%
                   dplyr::mutate(rt = round(rt*60)) %>%
                   dplyr::select(id:mz, rt, dplyr::everything())
               } else {
                 lib_meta <- lib_meta %>%
                   dplyr::mutate(rt = rt_hilic_neg) %>%
                   dplyr::mutate(rt = round(rt*60)) %>%
                   dplyr::select(id:mz, rt, dplyr::everything())
               }
             }

             if (column == 'c18') {
               if (polarity == 'positive') {
                 lib_meta <- lib_meta %>%
                   dplyr::mutate(rt = rt_c18_pos) %>%
                   dplyr::mutate(rt = round(rt*60)) %>%
                   dplyr::select(id:mz, rt, dplyr::everything())
               } else {
                 lib_meta <- lib_meta %>%
                   dplyr::mutate(rt = rt_c18_neg) %>%
                   dplyr::mutate(rt = round(rt*60)) %>%
                   dplyr::select(id:mz, rt, dplyr::everything())
               }
             }

             # modify lib_spec
             if (polarity == 'positive') {
               lib_spec <- ms2_lib$positive
             } else {
               lib_spec <- ms2_lib$negative
             }

             idx_spec <- match(unique(lib_meta$id), names(lib_spec))
             idx_spec <- which(!is.na(idx_spec))
             ms2_id <- unique(lib_meta$id)[idx_spec]
             lib_spec <- match(ms2_id, names(lib_spec)) %>%
               lib_spec[.] %>%
               lapply(., function(x){
                 x[[ce]]
               })

             # rm null spec
             idx_null <- sapply(lib_spec, is.null) %>% which()
             lib_spec <- lib_spec[-idx_null]
             result <- list(lib_meta = lib_meta,
                            lib_spec = lib_spec)

             return(result)
           }
)




#   calibrateRT ----------------------------------------------------------------

#' @title calibrateRT
#' @description RT calibration according to the RTQC
#' @author Zhiwei Zhou, Mingdu Luo
#' \email{zhouzw@@sioc.ac.cn}
#' @param file_rt_ref Default: 'RT_recalibration_table.csv'
#' @param lib_rt
#' @param is_rt_calibration
#' @param is_plot Default: TRUE
#' @param method_lc 'Amide12min', 'Amide23min'. Default:  'Amide12min'
#' @param column Default: 'hilic'
#' @param path '.'

setGeneric(name = 'calibrateRT',
           def = function(file_rt_ref = 'RT_recalibration_table.csv',
                          lib_rt,
                          is_rt_calibration = TRUE,
                          is_plot = TRUE,
                          method_lc = c('Amide12min', 'Amide23min'),
                          column = 'hilic',
                          path = '.'){

             # Load in house RT library according to LC method
             switch(method_lc,
                    "Amide12min" = {
                      data('rt_ref', envir = environment())
                      ref_rtqc_table <- rt_ref[[1]]
                      lc_start <- 0
                      lc_end <- 720},

                    "Amide23min" = {
                      data('rt_ref', envir = environment())
                      ref_rtqc_table <- rt_ref[[2]]
                      lc_start <- 0
                      lc_end <- 1380}
             )

             # Check whether to do RT recalibration, only hilic method were allowed ###
             if (any(!is_rt_calibration, column != 'hilic')){
               cat('RT recalibration was turned off.\n')
             } else {

               ### Check existence of RT calibration table ###
               if (!("RT_recalibration_table.csv" %in% dir(path))){
                 stop("There was no 'RT_recalibration_table.csv' table in the file folder.\n")
               } else {

                 exp_rtqc_rt <- read.csv(file.path(path, file_rt_ref), stringsAsFactors = F)

                 ### Check the format of exp_rtqc_rt ###
                 if (!identical(colnames(exp_rtqc_rt)[1:6], c("compound.name", "id.zhulab", "id.pubchem", "ref.mz", "rt", "polarity"))){
                   stop("Please check the format of 'RT_recalibration_table.csv' table and correct it according to our tutorial.\n")
                 } else {

                   ### Begin RT calibration ###
                   if (max(exp_rtqc_rt$rt) < 60) {
                     exp_rtqc_rt$rt <- round(exp_rtqc_rt$rt*60, digits = 4)
                   } else {
                     exp_rtqc_rt$rt <- round(exp_rtqc_rt$rt, digits = 4)
                   }

                   # the index of match rt compound and the warning
                   idx <- match(toupper(ref_rtqc_table$name),
                                toupper(exp_rtqc_rt$id.zhulab))

                   if (sum(!is.na(idx)) < 7){
                     warning(paste0("The number of used compouds for RT recalibration with LOESS was ",
                                    sum(!is.na(idx)),
                                    " and might be insufficient for a good performance.\n\n"))
                   } else {
                     cat("The number of used compouds for RT recalibration with LOESS was ",
                         sum(!is.na(idx)),
                         ".\n\n",
                         sep ='')
                   }

                   training.data <- data.frame(ref.rt = c(lc_start, ref_rtqc_table$rt[!is.na(idx)], lc_end),
                                               exp.rt = c(lc_start, exp_rtqc_rt$rt[idx[!is.na(idx)]], lc_end))

                   rownames(training.data) <- c('Start', exp_rtqc_rt$id.zhulab[idx[!is.na(idx)]], 'End')

                   ### rt.recalibration.model <- lm(exp.rt~ref.rt, data = training.data)
                   rt.recalibration.model <- loess(exp.rt~ref.rt,
                                                   data = training.data,
                                                   span = 0.75, degree = 2)

                   new.data <- data.frame(ref.rt = lib_rt$rt, stringsAsFactors = FALSE)
                   lib.calibrated.rt <- round(predict(object =  rt.recalibration.model,
                                                      newdata = new.data),
                                              digits = 2)

                   lib_rt$rt <- lib.calibrated.rt

                   result <- list(lib_rt = lib_rt,
                                  training.data = training.data,
                                  rt.recalibration.model = rt.recalibration.model)

                   dir.create(file.path(path, "01_result_initial_seed_annotation", '00_intermediate_data'), showWarnings = FALSE, recursive = TRUE)
                   save(result,
                        file = file.path(path, "01_result_initial_seed_annotation", '00_intermediate_data', 'rt_calibration_result'),
                        version = 2)

                   cat("RT recalibration was done.\n")

                   if (is_plot) {
                     # browser()
                     dir.create(file.path(path, "01_result_initial_seed_annotation", '01_rt_calibration_plot'),
                                showWarnings = FALSE, recursive = TRUE)
                     plotRtCalibration(result_rt_calibration = result[[2]],
                                       rt_recalibration_model = result[[3]],
                                       path = file.path(path, "01_result_initial_seed_annotation", '01_rt_calibration_plot'))
                   }

                   return(result)
                   rm(c(idx, lc_start, lc_end))

                 }
               }
             }
           })





#   calculateMz ----------------------------------------------------------------
setGeneric(name = 'calculateMz',
           def = function(
    exact_mass,
    adduct,
    delta_mz,
    nmol = NULL,
    ncharge = NULL
           ){

             if (length(nmol) == 0) {
               if (stringr::str_detect(adduct, pattern = '2M')) {
                 mz <- exact_mass*2 + delta_mz
               } else if (stringr::str_detect(adduct, pattern = '3M')) {
                 mz <- exact_mass*3 + delta_mz
               } else {
                 mz <- exact_mass + delta_mz
               }
             } else {
               mz <- exact_mass*nmol + delta_mz
             }


             if (length(ncharge) == 0) {
               if (stringr::str_detect(adduct, pattern = '\\]2\\-|\\]2\\+')) {
                 mz <- mz/2
               } else if (stringr::str_detect(adduct, pattern = '\\]3\\-|\\]3\\+')) {
                 mz <- mz/3
               } else {
                 mz
               }
             } else {
               mz <- mz/ncharge
             }

             mz
           })

################################################################################
# read_ms1 ---------------------------------------------------------------------

#' @title read_ms1
#' @author Zhiwei Zhou
#' @param filename The name of csv file.
#' @param instrument "SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap", 'ThermoExploris', "AgilentDTIMMS", "BrukerTIMS", 'WatersQTOF','WatersTWIMMS'
#' @param path '.'
#' @return A list. info: a table of feature information, column 1_feature name, column 2_mz, column3_rt in seconds, column4_ccs. subject: a table of peak area.
#' @importFrom magrittr '%>%'
#' @export

setGeneric('read_ms1',
           def = function(filename,
                          path = '.') {

             options(readr.num_columns = 0)
             temp <- readr::read_csv(file.path(path, filename),
                                     col_types = readr::cols())

             if (sum(tolower(colnames(temp)[1:3])==c('name', "mz", "rt"))!=3){
               stop("Please check the format of 'MS1_table.csv'. [The name of column 1-3 must be named as 'name', 'mz', 'rt']")
             }

             colnames(temp)[1:3] <- c("name", "mz", "rt")
             info <- data.frame(name = temp$name,
                                mz = as.numeric(temp$mz),
                                rt = as.numeric(temp$rt),
                                ccs = NA,
                                stringsAsFactors = F)

             if (max(info$rt) < 60) {
               info$rt <- info$rt*60
             }

             subject <- temp[, -c(1:3), drop = F]
             rownames(subject) <- info$feature

             result <- list(info=info, subject=subject)

             return(result)

           })


#   featureReName --------------------------------------------------------------
#' @title featureReName
#' @author Zhiwei Zhou
#' @description define feature name of replicates.
#' @param name a vector. character.

setGeneric('featureReName',
           def = function(name){
             sapply(seq(length(name)), function(i){
               temp.name <- name[i]
               temp.idx <- which(name==temp.name)
               if (length(temp.idx) > 1) {
                 temp.name <- paste(temp.name, which(temp.idx==i), sep = "_")
               }
               temp.name
             })
           })



#   getMzRange -----------------------------------------------------------------
setGeneric(name = 'getMzRange',
           def = function(mz,
                          ppm = 25,
                          mz_ppm_thr = 400){
             result <- sapply(mz, function(x) {
               if (x >= mz_ppm_thr) {
                 x * (1 + c(-1, 1) * ppm * 1e-6)
               } else {
                 temp1 <- x + mz_ppm_thr * c(-1, 1) * ppm * 1e-6
               }
             })

             t(result)
           }
)

#   SXTMTmatch -----------------------------------------------------------------
#' @title SXTMTmatch
#' @description Match two data according to mz and RT.
#' @author Xiaotao Shen, Zhiwei Zhou
#' @param data1 First data for matching, first column must be mz and seconod column must be rt.
#' @param data2 Second data for matching, first column must be mz and seconod column must be rt.
#' @param mz.tol mz tol for ms1 and ms2 data matching. ppm
#' @param rt.tol RT tol for ms1 and ms2 data matching. s or %
#' @return Return a result which give the matching result of data1 and database.
# export

setGeneric(name = "SXTMTmatch",
           def = function(data1,
                          data2,
                          mz.tol = 25,
                          #rt.tol is relative
                          rt.tol = 10,
                          ccs.tol = NULL,
                          rt.error.type = c("abs", "relative")){
             rt.error.type <- match.arg(rt.error.type)

             if (nrow(data1) == 0 | nrow(data2) == 0) {
               result <- NULL
               return(result)
             }

             if (length(ccs.tol) == 0) {
               info1 <- data1[,c(1,2)]
               info1 <- apply(info1, 1, list)

               mz2 <- as.numeric(data2[, 1])
               rt2 <- as.numeric(data2[, 2])
             } else {
               info1 <- data1[,c(1,2,3)]
               info1 <- apply(info1, 1, list)

               mz2 <- as.numeric(data2[, 1])
               rt2 <- as.numeric(data2[, 2])
               ccs2 <- as.numeric(data2[, 3])
             }


             result <- pbapply::pblapply(info1, function(x) {
               temp.mz1 <- x[[1]][[1]]
               temp.rt1 <- x[[1]][[2]]
               mz.error <- abs(temp.mz1 - mz2) * 10 ^ 6 / temp.mz1
               if (rt.error.type == "relative") {
                 rt.error <- abs(temp.rt1 - rt2) * 100 / temp.rt1
               } else {
                 rt.error <- abs(temp.rt1 - rt2)
               }

               if (length(ccs.tol) == 0) {
                 j <- which(mz.error <= mz.tol & rt.error <= rt.tol)

                 if (length(j) == 0) {
                   matrix(NA, ncol = 7)
                 } else {
                   cbind(j, temp.mz1, mz2[j], mz.error[j], temp.rt1, rt2[j], rt.error[j])
                 }

               } else {
                 temp.ccs1 <- x[[1]][[3]]
                 ccs.error <- abs(temp.ccs1 - ccs2)*100/temp.ccs1

                 j <- which(mz.error <= mz.tol & rt.error <= rt.tol & ccs.error <= ccs.tol)

                 if (length(j) == 0) {
                   matrix(NA, ncol = 10)
                 } else {
                   cbind(j, temp.mz1, mz2[j], mz.error[j], temp.rt1, rt2[j], rt.error[j], temp.ccs1, ccs2[j], ccs.error[j])
                 }
               }

             })

             # add a column of number
             if (length(result) == 1) {
               result <- cbind(1,result[[1]])
             } else {
               result <- mapply(function(x,y){list(cbind(x,y))},
                                x <- 1:length(info1),
                                y = result)
               result <- do.call(rbind, result)
             }


             if (length(ccs.tol) == 0) {
               result <- matrix(result[which(!apply(result,1,function(x) any(is.na(x)))),], ncol = 8)
               if(nrow(result) == 0) return(NULL)
               colnames(result) <-
                 c("Index1",
                   "Index2",
                   "mz1",
                   "mz2",
                   "mz error",
                   "rt1",
                   "rt2",
                   "rt error")

               return(result)
             } else {
               result <- matrix(result[which(!apply(result,1,function(x) any(is.na(x)))),], ncol = 11)
               if(nrow(result) == 0) return(NULL)
               colnames(result) <-
                 c("Index1",
                   "Index2",
                   "mz1",
                   "mz2",
                   "mz error",
                   "rt1",
                   "rt2",
                   "rt error",
                   'ccs1',
                   'ccs2',
                   'ccs error')

               return(result)

             }

           })



################################################################################
# match_ms1_rt -----------------------------------------------------------------

#' @title match_ms1_rt
#' @author Zhiwei Zhou
#' @param ms1_data ms1 data, including info and subject.
#' @param lib_meta the meta information of library
#' @param mz_tol m/z tolerance. Default: 10ppm
#' @param mz_ppm_thr
#' @param pf_rt_range Default: 0
#' @param tolerance_rt_range Defaule: 20
#' @param is_combined_rt Defaule: FALSE
#' @param is_rt_score Defaule: FALSE

# temp <- match_ms1_rt(ms1_data = ms1_data,
#                      lib_meta = lib_meta,
#                      mz_tol = 15,
#                      mz_ppm_thr = 150,
#                      pf_rt_range = 0,
#                      tolerance_rt_range = 15,
#                      is_combined_rt = FALSE,
#                      is_rt_score = FALSE)


setGeneric(name = 'match_ms1_rt',
           def = function(ms1_data,
                          lib_meta,
                          mz_tol = 15,
                          mz_ppm_thr = 150,
                          pf_rt_range = 0,
                          tolerance_rt_range = 20,
                          is_combined_rt = FALSE,
                          # is_filter = TRUE,
                          is_rt_score = FALSE,
                          is_check_cation = FALSE
           ){
             cat("m/z & RT & CCS match...\n\n")
             ms1_info <- ms1_data$info
             result_annotation <- convertMs1Data2SpecAnnotationClass(ms1_info = ms1_info)

             # browser()
             exp_mz_range <- getMzRange(mz = ms1_info$mz,
                                        ppm = mz_tol,
                                        mz_ppm_thr = mz_ppm_thr)

             idx_match <- lapply(seq_along(result_annotation), function(i){
               temp_idx <- which(lib_meta$mz >= exp_mz_range[i,1] & lib_meta$mz <= exp_mz_range[i,2])
               temp_idx
             })

             # generate indexs
             if (is_rt_score) {
               # cat("rt match\n")
               exp_rt_range <- getRtRange(data = ms1_info$rt,
                                          abs_rt_dev = tolerance_rt_range,
                                          is_combined = is_combined_rt)

               temp_idx <- lapply(seq_along(result_annotation), function(i){
                 temp_idx <- which(lib_meta$rt >= exp_rt_range[i,1] & lib_meta$rt <= exp_rt_range[i,2])
                 # # include no experimental RT candidates
                 # temp_idx2 <- which(is.na(lib_meta$rt))
                 # result <- unique(c(temp_idx, temp_idx2))
               })
               idx_match <- mapply(function(x, y){
                 intersect(x,y)
               },
               x=idx_match,
               y=temp_idx)
             }


             # add scores
             result <- pbapply::pblapply(seq_along(result_annotation), function(i){
               # cat(i);cat(' ')
               if (length(idx_match[[i]]) > 0) {
                 temp_idx <- idx_match[[i]]
                 lib_mz <- lib_meta$mz[temp_idx]
                 mz_error <- round(abs(ms1_info$mz[i]-lib_mz)/lib_mz*10^6,
                                   digits = 1)


                 if (is_rt_score) {
                   temp_rt <- ms1_data$info$rt[i]
                   temp_lib_rt <- lib_meta$rt[temp_idx]
                   delta_rt <- abs(temp_rt-temp_lib_rt)
                   rt_score <- getTrapezoidalScore(delta = delta_rt,
                                                   pf_range = pf_rt_range,
                                                   range = tolerance_rt_range)
                   rt_score <- round(rt_score, digits = 2)
                   rt_error <- round(delta_rt, digits = 1)
                 } else {
                   rt_score <- -1
                   rt_error <- -1
                 }

                 if (is_check_cation) {
                   idx_check_cation <- mapply(function(x, y){
                     if (x %in% c('[M]+', '[M-2H]-')) {
                       !is.na(y)
                     } else {
                       return(TRUE)
                     }
                   },
                   x = lib_meta$adduct[temp_idx],
                   y = lib_meta$note[temp_idx]) %>%
                     which()

                   if (length(idx_check_cation) == 0) {
                     return(NULL)
                   } else {
                     temp_idx <- temp_idx[idx_check_cation]
                   }

                   result <- tibble::tibble(idx = temp_idx,
                                            id = lib_meta$id[temp_idx],
                                            name = lib_meta$name[temp_idx],
                                            formula = lib_meta$formula[temp_idx],
                                            smiles = lib_meta$smiles[temp_idx],
                                            inchikey = lib_meta$inchikey[temp_idx],
                                            inchikey1 = lib_meta$inchikey1[temp_idx],
                                            adduct = lib_meta$adduct[temp_idx],
                                            mz_lib = lib_meta$mz[temp_idx],
                                            rt_lib = lib_meta$rt[temp_idx],
                                            ccs_lib = lib_meta$ccs[temp_idx],
                                            mz_error = mz_error[idx_check_cation],
                                            # mz_score=mz_score,
                                            rt_error = rt_error[idx_check_cation],
                                            rt_score = rt_score[idx_check_cation],
                                            ccs_error = ccs_error[idx_check_cation],
                                            ccs_score = ccs_score[idx_check_cation],
                                            msms_score_forward = -1,
                                            msms_score_reverse = -1,
                                            msms_matched_frag = -1)

                   return(result)
                 }

                 result <- tibble::tibble(idx = temp_idx,
                                          id = lib_meta$id[temp_idx],
                                          name = lib_meta$name[temp_idx],
                                          formula = lib_meta$formula[temp_idx],
                                          smiles = lib_meta$smiles[temp_idx],
                                          inchikey = lib_meta$inchikey[temp_idx],
                                          inchikey1 = lib_meta$inchikey1[temp_idx],
                                          adduct = lib_meta$adduct[temp_idx],
                                          mz_lib = lib_meta$mz[temp_idx],
                                          rt_lib = lib_meta$rt[temp_idx],
                                          ccs_lib = lib_meta$ccs[temp_idx],
                                          mz_error = mz_error,
                                          # mz_score=mz_score,
                                          rt_error = rt_error,
                                          rt_score = rt_score,
                                          msms_score_forward = -1,
                                          msms_score_reverse = -1,
                                          msms_matched_frag = -1)

                 return(result)
               }

             })


             result_annotation <- mapply(function(x, y){
               if (length(y)>0) {
                 x@annotation_result <- y
                 return(x)
               } else {
                 return(x)
               }
             },
             x = result_annotation,
             y = result,
             SIMPLIFY = FALSE)

             return(result_annotation)

           }
)


  # getRtRange -----------------------------------------------------------------
#' @title getRtRange
#' @author Zhiwei Zhou
#' @param data Numeric.
#' @param abs_rt_dev Numeric. Absolute rt deviation tolerance in seconds.
#' @param rel_rt_dev Numeric. Relative rt deviation tolerance in percent %.
#' @return a table of range. column 1: minimum value; column 2: maxmium value.

# getRtRange(data = 300, abs_rt_dev = 30)
# getRtRange(data = 600, rel_rt_dev = 6)
# getRtRange(data = 600, abs_rt_dev = 30, rel_rt_dev = 6, is_combined = TRUE)
# getRtRange(data = 300, abs_rt_dev = 30, rel_rt_dev = 6, is_combined = TRUE)

setGeneric(name = 'getRtRange',
           def = function(data,
                          abs_rt_dev,
                          rel_rt_dev,
                          rt_threshold = 500,
                          is_combined = FALSE){

             if (is_combined) {
               if (!any(missing(abs_rt_dev), missing(rel_rt_dev))) {
                 rt_range <- t(sapply(data, function(x){
                   if (x <= rt_threshold) {
                     result <- x+c(-1, 1)*abs_rt_dev
                   } else {
                     result <- x*(1+c(-1, 1)*rel_rt_dev/100)
                   }

                   result

                 }))
               }

             } else {
               if (!missing(abs_rt_dev)) {
                 rt_range <- t(sapply(data, function(x){
                   x+c(-1, 1)*abs_rt_dev
                 }))
               }

               if (!missing(rel_rt_dev)) {
                 rt_range <- t(sapply(data, function(x){
                   x*(1+c(-1, 1)*rel_rt_dev/100)
                 }))
               }
             }



             return(rt_range)

           })


  # getTrapezoidalScore --------------------------------------------------------
#' @title getTrapezoidalScore
#' @author Zhiwei Zhou
#' @description Trapezodial scoring function
#' @param delta Numeric. The difference between experiment value and library value.
#' @param pf_range Numeric. Penalty free range.
#' @param range Numeric. The maxmium tolerance.

setGeneric('getTrapezoidalScore',
           def = function(delta,
                          pf_range,
                          range){
             delta <- abs(delta)
             delta[delta <= pf_range] <- pf_range
             score <- 1-((delta-pf_range)/(range-pf_range))
             score[score<0] <- 0
             return(score)
           })


  # getLinerScore --------------------------------------------------------------
#' @title getLinerScore
#' @author Zhiwei Zhou
#' @description calculate liner score
#' @param delta Numeric. The difference between experiment value and library value.
#' @param tolerance The maxmium tolerance

setGeneric(name = 'getLinerScore',
           def = function(delta,
                          tolerance){
             score <- 1 - abs(delta) / tolerance
             score
           })


  # SpecAnnotationClass --------------------------------------------------------
setClass(Class = "SpecAnnotationClass",
         representation(peak_info = "list",
                        annotation_result = "data.frame")
)

setMethod(f = "show",
          signature = "SpecAnnotationClass",
          definition = function(object) {
            cat("-----------Meta information------------\n")
            cat("Feature:", object@peak_info$name, "\n")
            cat("m/z:", object@peak_info$mz, "\n")
            cat("RT (s):", object@peak_info$rt, "\n")
            cat("With ms2:", ifelse(object@peak_info$has_ms2 > 0, 'Yes', 'No'), "\n")
            # cat("CCS (A2):",object@peak_info$ccs, "\n")
            # cat("Median peak area:",object@peak_info$intmed, '\n\n')

            cat("-----------Annotation result------------\n")

            if(nrow(object@annotation_result) == 0) {
              cat("No annotation.\n\n")
            }else{
              for(i in 1:nrow(object@annotation_result)){
                cat("\n")
                cat("Annotation:", i, '\n')
                cat("ID:", object@annotation_result$id[i], "\n")
                cat("Name:", object@annotation_result$name[i], "\n")
                cat("Formula:", object@annotation_result$formula[i], "\n")
                cat("SMILES:", object@annotation_result$smiles[i], "\n")
                cat("InChIKey:", object@annotation_result$inchikey[i], "\n")
                cat("InChIKey1:", object@annotation_result$inchikey1[i], "\n")
                cat("Adduct:", object@annotation_result$adduct[i], "\n")
                cat("Lib m/z:", object@annotation_result$mz_lib[i], "\n")
                cat("Lib RT:", object@annotation_result$rt_lib[i], "\n")
                # cat("Lib CCS:", object@annotation_result$ccs_lib[i], "\n")
                cat("m/z error (ppm):", object@annotation_result$mz_error[i], "\n")
                cat("RT error (s):", object@annotation_result$rt_error[i], "\n")
                cat("RT score:", object@annotation_result$rt_score[i], "\n")
                # cat("CCS error (%):", object@annotation_result$ccs_error[i], "\n")
                # cat("CCS score:", object@annotation_result$ccs_score[i], "\n")
                cat("MS2 score (forward):", object@annotation_result$msms_score_forward[i], "\n")
                cat("MS2 score (reverse):", object@annotation_result$msms_score_reverse[i], "\n")
                cat("Number of matched fragments:", object@annotation_result$msms_matched_fra[i], "\n\n")
              }

            }
          }
)

  # convertMs1Data2SpecAnnotationClass -----------------------------------------
#' @title convertMs1Data2SpecAnnotationClass
#' @description Convert Ms1DataInfor to AnnotationResult class
#' @author Zhiwei Zhou
#' @param ms1_info

# test <- Ms1DataInfo2AnnotationResult(ms1_info = ms1_info)

setGeneric(name = 'convertMs1Data2SpecAnnotationClass',
           def = function(
    ms1_info
           ){
             cat('Convert to SpecAnnotationClass...\n')
             temp_annotation_result <- tibble::tibble(idx = numeric(),
                                                      id = character(),
                                                      name = character(),
                                                      formula = character(),
                                                      smiles = character(),
                                                      inchikey = character(),
                                                      inchikey1 = character(),
                                                      adduct = character(),
                                                      mz_lib = numeric(),
                                                      rt_lib = numeric(),
                                                      ccs_lib = numeric(),
                                                      mz_error = numeric(),
                                                      # mz_score=mz_score,
                                                      rt_error = numeric(),
                                                      rt_score = numeric(),
                                                      msms_score_forward = numeric(),
                                                      msms_score_reverse = numeric(),
                                                      msms_matched_frag = numeric())


             result <- pbapply::pblapply(seq_along(ms1_info$name), function(i){
               peak_info <- ms1_info[i,] %>% as.list()
               peak_info$has_ms2 <- -1
               new(Class = "SpecAnnotationClass",
                   peak_info = peak_info,
                   annotation_result = temp_annotation_result)
             })

             names(result) <- ms1_info$name

             return(result)
           }
)



################################################################################
# read_ms2 ----------------------------------------------------------------------

#' @title read_ms2
#' @author Zhiwei Zhou
#' @param ms2_file MS/MS file names
#' @param ms2_type 'mgf', 'cef', 'msp', 'mzXML', 'mzML'
#' @export

# path <- 'H:/00_projects/03_MetDNA2/00_data/20200922_metdna2_development/pos_201007_metdna2/'
# temp_files <- dir(path) %>% stringr::str_detect('\\.mgf') %>% dir(path)[.]
# ms2_data <- readMs2(ms2_file = file.path(path, temp_files),
#                     ms2_type = 'mgf')
#
# path <- '/home/zhouzw/Data_processing/20210315_tims_data_processing_metdna/'
# temp_files <- dir(path) %>% stringr::str_detect('\\.mgf') %>% dir(path)[.]
# test <- readMs2(ms2_file = file.path(path, temp_files),
#                 ms2_type = 'mgf',
#                 instrument = 'IMMS')

setGeneric(name = 'read_ms2',
           def = function(
    ms2_file,
    ms2_type = c('mgf', 'cef', 'msp', 'mzXML', 'mzML'),
    ...
           ){
             ms2_type <- match.arg(ms2_type)

             cat('Read MS2...\n\n')

             switch(ms2_type,
                    'mgf' = {
                      pbapply::pboptions(type='timer', char='+')
                      ms2_data <- pbapply::pblapply(ms2_file, function(x){
                        readMGF(x, ...)
                      })
                    },

                    'cef' = {
                      pbapply::pboptions(type='timer', char='+')
                      ms2_data <- pbapply::pblapply(ms2_file, function(x){
                        readCEF(x, ...)
                      })
                    },

                    'msp' = {
                      pbapply::pboptions(type='timer', char='+')
                      ms2_data <- pbapply::pblapply(ms2_file, function(x){
                        readMSP(x, ...)
                      })
                    },

                    'mzXML' = {
                      pbapply::pboptions(type='timer', char='+')
                      ms2_data <- pbapply::pblapply(ms2_file, function(x){
                        readMzXML(x, ...)
                      })

                    },

                    'mzML' = {
                      pbapply::pboptions(type='timer', char='+')
                      ms2_data <- pbapply::pblapply(ms2_file, function(x){
                        readMzML(x, ...)
                      })

                    }
             )

             file_name <- basename(ms2_file)
             ms2_data <- lapply(seq_along(ms2_data), function(i){
               temp_ms2 <- ms2_data[[i]]
               result <- lapply(temp_ms2, function(x){
                 x$info <- c(x$info, 'file_name_idx' = i)
                 return(x)
               })
               return(result)
             })

             ms2_data <- do.call(c, ms2_data)
           }
)



#   readMGF ----------------------------------------------------------------------

#' @title readMGF
#' @author Zhiwei Zhou
#' @param file the name of mgf file
#' @export

# path <- 'H:/00_projects/03_MetDNA2/00_data/20200922_metdna2_development/pos_201007_metdna2'
# temp_files <- dir(path) %>% stringr::str_detect('\\.mgf') %>% dir(path)[.]
# test <- readMGF(file.path(path, temp_files[1]))

# path <- '/home/zhouzw/Data_processing/20210315_tims_data_processing_metdna/'
# temp_files <- dir(path) %>% stringr::str_detect('\\.mgf') %>% dir(path)[.]
# test <- readMGF(file.path(path, temp_files[1]))

setGeneric('readMGF',
           def = function(file,
                          ...){

             mgf.data <- ListMGF(file)
             whether_has_ccs <- stringr::str_detect(mgf.data[[1]], pattern = '^CCS') %>% any()

             if (!whether_has_ccs) {
               # mgf.data <- ListMGF(file)
               # nl.spec <- grep('^\\d', mgf.data)
               nl.spec <- lapply(mgf.data, function(x) grep('^\\d', x)) # select spec index
               info.mz <- lapply(mgf.data, function(x) grep('^PEPMASS', x, value = T))
               info.rt <- lapply(mgf.data, function(x) grep('^RTINSECONDS', x, value = T))

               info.mz <- sapply(info.mz, function(x){
                 temp <- gsub(pattern = "\\w+=", "", x)
                 if (stringr::str_detect(temp, pattern = "\ ")) {
                   temp <- as.numeric(substr(temp, start = 1, stop = regexpr("\ ", temp)))
                 }

                 temp <- as.numeric(temp)
                 # temp <- round(as.numeric(temp), digits = 4)
               })

               info.rt <- unlist(info.rt)
               # info.rt <- round(as.numeric(gsub(pattern = "\\w+=", "", info.rt)), digits = 0)
               info.rt <- as.numeric(gsub(pattern = "\\w+=", "", info.rt))

               spec <- mapply(function(x, y){
                 list(do.call(rbind, strsplit(x[y], split = " ")))
               },
               x = mgf.data,
               y = nl.spec)

               spec <- lapply(seq(length(info.mz)), function(i){
                 # data.frame(mz=as.numeric(spec[1:10, i]), intensity=as.numeric(spec[11:20, i]), stringsAsFactors = F)
                 data.frame(mz=as.numeric(spec[[i]][, 1]), intensity=as.numeric(spec[[i]][, 2]), stringsAsFactors = F)
               })

               ms2 <- lapply(seq(length(info.mz)), function(i){
                 temp.info <- c(info.mz[i], info.rt[i])
                 names(temp.info) <- c("mz", "rt")
                 temp.spec <- spec[[i]]
                 temp.result <- list(info=temp.info, spec=temp.spec)
               })

               return(ms2)
             } else {
               # mgf.data <- ListMGF(file)
               # nl.spec <- grep('^\\d', mgf.data)
               nl.spec <- lapply(mgf.data, function(x) grep('^\\d', x)) # select spec index
               info.mz <- lapply(mgf.data, function(x) grep('^PEPMASS', x, value = T))
               info.rt <- lapply(mgf.data, function(x) grep('^RTINSECONDS', x, value = T))
               info.ccs <- lapply(mgf.data, function(x) grep('^CCS', x, value = T))

               info.mz <- sapply(info.mz, function(x){
                 temp <- gsub(pattern = "\\w+=", "", x)
                 if (stringr::str_detect(temp, pattern = "\ ")) {
                   temp <- as.numeric(substr(temp, start = 1, stop = regexpr("\ ", temp)))
                 }

                 temp <- as.numeric(temp)
                 # temp <- round(as.numeric(temp), digits = 4)
               })

               info.rt <- unlist(info.rt)
               # info.rt <- round(as.numeric(gsub(pattern = "\\w+=", "", info.rt)), digits = 0)
               info.rt <- as.numeric(gsub(pattern = "\\w+=", "", info.rt))

               info.ccs <- unlist(info.ccs)
               info.ccs <- as.numeric(gsub(pattern = "\\w+=", "", info.ccs))


               spec <- mapply(function(x, y){
                 list(do.call(rbind, strsplit(x[y], split = " ")))
               },
               x = mgf.data,
               y = nl.spec)

               spec <- lapply(seq(length(info.mz)), function(i){
                 # data.frame(mz=as.numeric(spec[1:10, i]), intensity=as.numeric(spec[11:20, i]), stringsAsFactors = F)
                 data.frame(mz=as.numeric(spec[[i]][, 1]), intensity=as.numeric(spec[[i]][, 2]), stringsAsFactors = F)
               })

               ms2 <- lapply(seq(length(info.mz)), function(i){
                 temp.info <- c(info.mz[i], info.rt[i], info.ccs[i])
                 names(temp.info) <- c("mz", "rt", 'ccs')
                 temp.spec <- spec[[i]]
                 temp.result <- list(info=temp.info, spec=temp.spec)
               })

               return(ms2)
             }

           })


#' @title ListMGF
#' @author Zhiwei Zhou
#' @param file the name of mgf file

setGeneric('ListMGF',
           def = function(file){
             mgf.data <- readLines(file)
             nl.rec.new <- 1
             idx.rec <- 1
             rec.list <- list()
             for(nl in 1:length(mgf.data))
             {
               if(mgf.data[nl]=="END IONS")
               {
                 rec.list[idx.rec] <- list(Compound = mgf.data[nl.rec.new : nl])
                 nl.rec.new <- nl + 1
                 idx.rec <- idx.rec + 1
               }
             }
             rec.list
           })



#   readMSP ----------------------------------------------------------------------
#' #' @title readMSP
#' #' @description  Read a MSP file and return a list of spectra for all feature with feature information
#' #' @param file path of the msp file
#' #' @export
#'
#' setGeneric('readMSP', function(file) {
#'   msp.data.list <- ListDB(file)
#'   nr.num.pk <- grep('Num Peaks', msp.data.list[[1]])
#'   info.spec <- lapply(msp.data.list, function(msp.data) {
#'     info.list <- strsplit(msp.data[1:nr.num.pk], split = ': ')
#'     info <- do.call(cbind, info.list)
#'     colnames(info) <- info[1, ]
#'     mz <- round(as.numeric(info[2,3]), digits = 4)
#'     rt <- round(as.numeric(strsplit(info[2,"Comment"], split = "_")[[1]][1])*60, digits = 0)
#'     name <- info[2,"Comment"]
#'
#'     # info <- data.frame(mz=mz, rt=rt)
#'     rownames(info) <- name
#'
#'     spec.list <- strsplit(msp.data[(nr.num.pk+1):length(msp.data)], ' |\\t')
#'     spec.list <- lapply(spec.list, as.numeric)
#'     spec <- do.call(rbind, spec.list)[, c(1, 2), drop = FALSE]
#'     colnames(spec) <- c('mz', 'intensity')
#'
#'     list('info' = info,
#'          'spec' = spec)
#'
#'   })
#'
#'   info.spec
#' })
#'
#' setGeneric('ListDB', function(file) {
#'   msp.data <- readLines(file)
#'   nl.db.new <- 1
#'   idx.db <- 1
#'   db.list <- list()
#'   len.data <- length(msp.data)
#'   for(nl in 1:len.data)
#'   {
#'     if(msp.data[nl]=="") {
#'       db.list[idx.db] <- list(Compound=msp.data[nl.db.new:(nl-1)])
#'       nl.db.new <- nl + 1
#'       idx.db <- idx.db + 1
#'     } else if (nl == len.data){
#'       db.list[idx.db] <- list(Compound=msp.data[nl.db.new:nl])
#'     }
#'   }
#'   db.list
#' })



#' @title readMSP
#' @description read MSP spectra files
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param file the file name
#' @param mode standard: extract name, mz and RT from the file, which fit for MSP data exported from QI software; all: extract all information in the file
#' @return
#' A list object. info: MS1 information, including ms1, rt etc.; spec: the MS/MS spectrum
#' @examples
#' test <- readMSP(file = 'F:/01 MetIMMS/00 data processing/190515 external validation msms data extraction/zhumetlib_validation_pos_20v_190520.msp', mode = 'all')

setGeneric('readMSP', function(file,
                               mode = c('all', 'standard'),
                               source = c('MetAnalyzer', 'MSDIAL', 'Other')) {
  # devtools::use_package('dplyr')

  mode <- match.arg(mode)
  source <- match.arg(source)
  msp.data.list <- ListDB(file)
  nr.num.pk <- grep('Num Peaks', stringr::str_to_title(msp.data.list[[1]]))
  info.spec <- lapply(msp.data.list, function(msp.data) {
    info.list <- strsplit(msp.data[1:nr.num.pk], split = ': ')
    info <- do.call(cbind, info.list)
    colnames(info) <- info[1, ]

    if (mode=='standard') {
      mz <- round(as.numeric(info[2,3]), digits = 4)
      rt <- round(as.numeric(strsplit(info[2, "Comment"], split = "_")[[1]][1])*60, digits = 0)
      name <- info[2, "Comment"]
      # info <- matrix(c(mz, rt), ncol = 2)
      info <- data.frame(mz=mz, rt=rt)
      rownames(info) <- name
      # colnames(info) <- c("mz", "rt")
    } else {
      info <- as.data.frame(tibble::as.tibble(info))
      info <- info[-1,,drop=F]

      if ('comment' %in% tolower(colnames(info))) {
        info$Comment <- stringr::str_replace(info$Comment, pattern = '\t', replacement = '')
      }

      if (info$NAME == 'Unknown') {
        source <- 'MSDIAL'
      }

      # change colnames for MSP file from MetAnalyzer
      switch (source,
              'MetAnalyzer' = {
                info <- info %>%
                  dplyr::rename(name = NAME,
                                mz = PRECURSORMZ,
                                npeaks = `Num Peaks`) %>%
                  dplyr::mutate(mz = as.numeric(mz))
              },
              'MSDIAL' = {
                info <- info %>%
                  dplyr::rename(name = NAME,
                                mz = PRECURSORMZ,
                                npeaks = `Num Peaks`) %>%
                  dplyr::mutate(mz = as.numeric(mz),
                                name = Comment)
              }
      )

      rownames(info) <- NULL
    }

    # if NULL spectra exported, return NULL (changed for MSDIAL)
    if (length(msp.data) <= nr.num.pk) {
      return(NULL)
    }

    spec.list <- strsplit(msp.data[(nr.num.pk+1):length(msp.data)], ' |\\t')
    spec.list <- lapply(spec.list, as.numeric)
    spec <- do.call(rbind, spec.list)[, c(1, 2), drop = FALSE]
    colnames(spec) <- c('mz', 'intensity')
    # spec <- list('spec' = spec)

    # list('info' = info[-1, , drop = FALSE],
    #      'spec' = spec)

    list('info' = info,
         'spec' = spec)

  })

  # remove NULL spectra
  idx_rm <- sapply(info.spec, function(x){length(x) == 0}) %>% which()
  if (length(idx_rm) > 0) {
    info.spec <- info.spec[-idx_rm]
  }

  return(info.spec)
})



setGeneric('ListDB', function(file) {
  msp.data <- readLines(file)
  nl.db.new <- 1
  idx.db <- 1
  db.list <- list()
  len.data <- length(msp.data)
  for(nl in 1:len.data)
  {
    if(msp.data[nl]=="") {
      db.list[idx.db] <- list(Compound=msp.data[nl.db.new:(nl-1)])
      nl.db.new <- nl + 1
      idx.db <- idx.db + 1
    } else if (nl == len.data){
      db.list[idx.db] <- list(Compound=msp.data[nl.db.new:nl])
    }
  }
  db.list
})


#   readCEF ----------------------------------------------------------------------
#' @title readCEF
#' @description  Read a CEF file and return a list of spectra for all feature
#'  with feature information
#' @param file path of the CEF file
#' @export

# file <- 'H:/00_projects/03_MetDNA2/00_data/20201012_metdna2_developmemt/CEF/mice_liver_pos/QC_01_MA-d3-c3_SR.cef'
# test <- readCEF(file)

setGeneric('readCEF', function(file) {
  cef.data.list <- ListCEF(file)
  # nr.num.pk: the number of
  # nr.num.pk <- grep('Num Peaks', cef.data.list[[1]])
  info.spec <- lapply(cef.data.list, function(cef.data) {

    # Extract basic information
    #  line 13, extract mz
    #  line 2, extract string rt, ccs "rt=\"\\d+.\\d+", example "rt=\"5.934"
    #    rt in seconds

    mz <- stringr::str_extract(string = cef.data[13], "\\d+.\\d+")
    mz <- round(as.numeric(mz), 4)

    rt <- stringr::str_extract(string = cef.data[2], "rt=\"\\d+.\\d+")
    ccs <- stringr::str_extract(string = cef.data[2], "ccs=\"\\d+.\\d+")

    rt <- round(as.numeric(stringr::str_extract(string = rt, pattern = "\\d+.\\d+"))*60, 0)
    ccs <- round(as.numeric(stringr::str_extract(string = ccs, pattern = "\\d+.\\d+")), 1)

    # info <- data.frame(mz=mz, rt=rt, ccs=ccs)
    info <- c(mz, rt, ccs)
    names(info) <- c('mz', 'rt', 'ccs')

    # extract MS/MS spectrum -------------------
    idx.spec.start <- stringr::str_locate(string = cef.data, pattern = "\\s+<MSPeaks>")
    idx.spec.end <- stringr::str_locate(string = cef.data, pattern = "\\s+</MSPeaks>")

    # search the label "<MSPeaks>"
    # 1st is the product ion, 2nd is the isopote

    idx.spec.start <- which(idx.spec.start[,1]>0)[1]+1
    idx.spec.end <- which(idx.spec.end[,1]>0)[1]-1

    spec.list.mz<- stringr::str_extract(cef.data[idx.spec.start:idx.spec.end],
                                        pattern = "x=\"\\d+.\\d+")
    spec.list.int <- stringr::str_extract(cef.data[idx.spec.start:idx.spec.end],
                                          pattern = "y=\"\\d+.\\d+")

    spec.list.mz <- round(as.numeric(stringr::str_extract(spec.list.mz, "\\d+.\\d+")),4)
    spec.list.int <- round(as.numeric(stringr::str_extract(spec.list.int, "\\d+.\\d+")),0)


    spec <- data.frame(mz=spec.list.mz, intensity=spec.list.int)
    spec <- as.matrix(spec)

    list('info' = info,
         'spec' = spec)

  })

  return(info.spec)
})



#' @title ListCEF
#' @author Zhiwei Zhou
#' @param file the name of CEF file

setGeneric('ListCEF',
           def = function(file){
             # cef.data <- readLines(file)
             cef.data <- readr::read_lines(file)

             idx.start <- stringr::str_locate(string = cef.data, "<Compound mppid=\"")
             idx.end <- stringr::str_locate(string = cef.data, "</Compound>")

             idx.start <- which(as.numeric(idx.start[,1])>0)
             idx.end <- which(as.numeric(idx.end[,1])>0)

             rec.list <- mapply(function(x, y){
               result <- cef.data[x:y]
               return(result)
             },
             x = idx.start,
             y = idx.end)

             return(rec.list)
           })


#   readMzXML ----------------------------------------------------------------------

#' @title readMzXML
#' @author Zhiwei Zhou, Xiaotao Shen
#' @param file
#' @export
# file <- 'H:/00_projects/03_MetDNA2/00_data/20201012_metdna2_developmemt/mzXML/POS/Sample2-W30 POS06-POS-W30 POS06.mzXML'
# test <- readMzXML(file)

setGeneric('readMzXML',
           def = function(file){
             mzxml.data <- mzR::openMSfile(file)
             mzxml.info <- mzR::header(mzxml.data)
             mzxml.peak <- mzR::peaks(mzxml.data)

             ms2.idx <- which(mzxml.info$msLevel == 2)
             ms2.info <- mzxml.info[ms2.idx, c("precursorMZ", "retentionTime")]
             ms2.info <- apply(ms2.info, 1, list)
             ms2.spec <- mzxml.peak[ms2.idx]

             result.ms2 <- mapply(function(x, y){
               # temp_info <- data.frame(mz = x[[1]][1],
               #                         rt = x[[1]][2],
               #                         stringsAsFactors = FALSE)

               temp_info <- as.numeric(c(x[[1]][1], x[[1]][2]))
               names(temp_info) <- c("mz", "rt")

               temp_spec <- y
               colnames(temp_spec) <- c("mz", "intensity")

               result <- list('info' = temp_info,
                              'spec' = temp_spec)

               return(result)
             },
             x = ms2.info,
             y = ms2.spec,
             SIMPLIFY = FALSE)

             return(result.ms2)

           })


#   readMzML -------------------------------------------------------------------

#' @title readMzML
#' @author Zhiwei Zhou, Xiaotao Shen
#' @param file
#' @export
#' file <- 'G:/00_projects/03_MetDNA2/10_project/MetDNA2_project/Data/20220704_mzML_test/mzml/NIST_urine01_pos-NIST_urine01.mzML'
# test <- readMzXML(file)

setGeneric('readMzML',
           def = function(file){
             mzxml.data <- mzR::openMSfile(file)
             mzxml.info <- mzR::header(mzxml.data)
             mzxml.peak <- mzR::peaks(mzxml.data)

             ms2.idx <- which(mzxml.info$msLevel == 2)
             ms2.info <- mzxml.info[ms2.idx, c("precursorMZ", "retentionTime")]
             ms2.info <- apply(ms2.info, 1, list)
             ms2.spec <- mzxml.peak[ms2.idx]

             result.ms2 <- mapply(function(x, y){
               # temp_info <- data.frame(mz = x[[1]][1],
               #                         rt = x[[1]][2],
               #                         stringsAsFactors = FALSE)

               temp_info <- as.numeric(c(x[[1]][1], x[[1]][2]))
               names(temp_info) <- c("mz", "rt")

               temp_spec <- y
               colnames(temp_spec) <- c("mz", "intensity")

               result <- list('info' = temp_info,
                              'spec' = temp_spec)

               return(result)
             },
             x = ms2.info,
             y = ms2.spec,
             SIMPLIFY = FALSE)

             return(result.ms2)

           })



################################################################################
# integrate_ms2 ----------------------------------------------------------------

#' @title integrate_ms2
#' @description Purify and integrate multiple MS/MS files
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param ms2_data ms2 data
#' @param is_include_precursor remove precursor peaks in raw spectra or not. Default: FALSE
#' @param is_deisotope remove isotope peaks in raw spectra or not. Default: TRUE
#' @param int_ms2_min_abs the minmium intensity of ms2 fragment. Default: 10
#' @param int_ms2_min_relative the minmium relative of ms2 fragment. Default: 0.01 ((1%))
#' @param ppm_precursor_filter the m/z tolerance to rexmove precursor ion in spectra. Default: 20
#' @param extract_ms2_file extract ms2 file. Default: NULL
#' @param mz_range_ms2 the range of ms2 data.
#' @export

# path <- 'H:/00_projects/03_MetDNA2/00_data/20200922_metdna2_development/pos_201007_metdna2/'
# temp_files <- dir(path) %>% stringr::str_detect('\\.mgf') %>% dir(path)[.]
# ms2_data <- read(ms2_file = file.path(path, temp_files),
#                     ms2_type = 'mgf')
# test <- inteMs2(ms2_data = ms2_data,
#                 ms2_type = 'mgf',
#                 is_include_precursor = TRUE,
#                 is_deisotope = FALSE,
#                 int_ms2_min_abs = 0,
#                 int_ms2_min_relative = 0.01,
#                 ppm_precursor_filter = 25,
#                 mz_range_ms2 = NULL)

setGeneric(name = 'integrate_ms2',
           def = function(
    ms2_data,
    ms2_file,
    ms2_type = c("mgf", "msp", "cef", 'mzXML', 'mzML'),
    is_include_precursor = TRUE,
    is_deisotope = FALSE,
    int_ms2_min_abs = 0,
    int_ms2_min_relative = 0.01,
    ppm_precursor_filter = 20,
    mz_range_ms2=NULL,
    ...
           ){
             match.arg(ms2_type)

             cat('Purify and integrate MS/MS spectra\n')
             tune_MS2_spec <- pbapply::pblapply(seq_along(ms2_data), function(i){
               # cat(i, ' ')
               result <- purifyMs2(spec = ms2_data[[i]]$spec,
                                   mz_precursor = ms2_data[[i]]$info['mz'],
                                   is_include_precursor = is_include_precursor,
                                   is_deisotope = is_deisotope,
                                   int_ms2_min_abs = int_ms2_min_abs,
                                   int_ms2_min_relative = int_ms2_min_relative,
                                   ppm_precursor_filter = ppm_precursor_filter,
                                   mz_range_ms2 = mz_range_ms2)


               return(result)
             })

             # add file names
             info <- lapply(seq(length(ms2_data)), function(i){
               ms2_data[[i]]$info
             })

             info <- do.call(rbind, info)

             if (ms2_type == "mgf") {
               info <- data.frame(mz=as.numeric(info[,1]),
                                  rt=as.numeric(info[,2]),
                                  filename = as.numeric(info[,3]),
                                  stringsAsFactors = F)
             }

             if (ms2_type == "msp") {
               temp <- rownames(info)
               info <- data.frame(name=as.character(info[,1]),
                                  mz=as.numeric(info[,2]),
                                  filename = as.numeric(info[,3]),
                                  stringsAsFactors = F)
               rownames(info) <- temp
             }

             if (ms2_type == 'mzXML') {
               temp <- rownames(info)
               info <- data.frame(mz=as.numeric(info[,1]),
                                  rt=as.numeric(info[,2]),
                                  filename = as.numeric(info[,3]),
                                  stringsAsFactors = F)
             }

             if (ms2_type == 'mzML') {
               temp <- rownames(info)
               info <- data.frame(mz=as.numeric(info[,1]),
                                  rt=as.numeric(info[,2]),
                                  filename = as.numeric(info[,3]),
                                  stringsAsFactors = F)
             }

             idx_null <- which(sapply(tune_MS2_spec, is.null))

             if (length(idx_null) > 0) {
               info <- info[-idx_null,,drop=F]
               tune_MS2_spec <- tune_MS2_spec[-idx_null]
             }

             file_name_list <- basename(ms2_file)
             info <- info %>%
               dplyr::mutate(filename = file_name_list[info$filename])
             result <- list(info=info, spec=tune_MS2_spec)
             return(result)

           }
)


  # purifyMs2 -------------------------------------------------------------------


#' @title purifyMs2
#' @author Zhiwei Zhou
#' @param spec the spectra matrix. column1: mz; column2: intensity
#' @param mz_precursor the precursor m/z. Numeric.
#' @param is_include_precursor remove precursor peak in raw spectra or not_ Default: False
#' @param mz_range_ms2 the range of m/z. Default: NULL
#' @param int_ms2_min_abs the minmium intensity of ms2 fragment_ Default: 10
#' @param int_ms2_min_relative the minmium relative of ms2 fragment_ Default: 0.03 ((3%))
#' @param ppm_precursor_filter the m/z tolerance to remove precursor ion in spectra_ Default: 20
#' @export

setGeneric(name = 'purifyMs2',
           def = function(
    spec,
    mz_precursor,
    is_include_precursor = FALSE,
    is_deisotope = FALSE,
    is_remove_ring_effect = TRUE,
    mz_range_ms2 = NULL,
    int_ms2_min_abs = 10,
    int_ms2_min_relative = 0.03,
    ppm_precursor_filter = 20
           ){
             # browser()
             spec <- spec[order(spec[, 'mz']), , drop = FALSE]

             # considering precursor ion
             if (missing(mz_precursor)) {
               mz_precursor <- max(spec[, 'mz'])
               mz_precursor <- as.numeric(mz_precursor)
             }

             mz_precursor_range <- getMzRange(mz_precursor, ppm_precursor_filter, mz_ppm_thr = 400)
             idx_mz_precursor_range <- ifelse(is_include_precursor, 2, 1)

             #change mz range depend precusor include or not
             mz_cutoff <- mz_precursor_range[idx_mz_precursor_range]
             spec <- spec[spec[,'mz'] < mz_cutoff, , drop = FALSE]

             if (!is.null(mz_range_ms2)) {
               nr_keep <- which(spec[, 'mz'] >= mz_range_ms2[1] &
                                  spec[, 'mz'] <= mz_range_ms2[2])
               if (length(nr_keep) > 0) {
                 spec <- spec[nr_keep, , drop = FALSE]
               }
               else {
                 return()
               }
             }

             # discarding low intensity spec (1% highest int and int.ms2.min.abs)
             int_cutoff <- max(max(spec[, 'intensity']) *
                                 int_ms2_min_relative,
                               int_ms2_min_abs)
             spec <- spec[spec[, 'intensity'] >= int_cutoff, , drop = FALSE]
             if (nrow(spec) == 0) {
               return()
             }

             if (is_deisotope) {
               spec <- removeIsotopes(spec)
             }

             # discarding ring effects
             if (is_remove_ring_effect) {
               spec <- removeRingEffect(spec)
             }

             return(spec)

           }
)



#' @title removeIsotopes
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param spec the matrix spectra. 1st column: "mz", 2nd column: "intensity"

setGeneric('removeIsotopes',
           def = function(
    spec
           ){
             # transform spec to data.frame
             temp_spec <- data.frame(spec, annotation='M', stringsAsFactors = F)

             is_filter_isotope <- TRUE

             while (is_filter_isotope) {
               idx_max <- which.max(temp_spec[,'intensity'])
               # cat(idx.max); cat(' ')

               mono_mz <- temp_spec[idx_max, 'mz']
               mono_int <- temp_spec[idx_max, 'intensity']

               temp_spec$intensity[idx_max] <- 0

               # generate potential isotope list
               isotopes_list <- mono_mz + c(1.0003, 2.0044)

               isotopes_list <- lapply(isotopes_list, function(x){
                 c(x, x + c(-1, 1) * 0.003)
               })

               isotopes_list <- do.call(rbind, isotopes_list)
               colnames(isotopes_list) <- c('mz', 'min', 'max')

               # m/z match
               # if multiple spectra were

               temp_idx_isotope <- lapply(seq(nrow(isotopes_list)), function(i){
                 idx_isotope <- which(temp_spec[,'mz'] >= isotopes_list[i,'min'] & temp_spec[,'mz'] <= isotopes_list[i,'max'])

                 if (length(idx_isotope) > 1) {
                   temp_mz_error <- abs(temp_spec[idx_isotope, 'mz'] - isotopes_list[i,'mz'])
                   idx_isotope <- idx_isotope[which.min(temp_mz_error)]
                 }

                 idx_isotope
               })


               # judge the intensity relationship
               # if matched M+1 and M+2, intensity relationship: M > M+1 > M+2
               # if matched M+1, intensity relationship: M > M+1

               if (length(temp_idx_isotope[[1]]) > 0 & length(temp_idx_isotope[[2]]) > 0) {
                 temp_m_1_int <- temp_spec$intensity[temp_idx_isotope[[1]]]
                 temp_m_2_int <- temp_spec$intensity[temp_idx_isotope[[2]]]

                 # intensity decrease
                 if (mono_int > temp_m_1_int & temp_m_1_int > temp_m_2_int) {
                   idx_isotope <- c(temp_idx_isotope[[1]], temp_idx_isotope[[2]])
                   temp_spec$annotation[idx_isotope] <- c('M+1', 'M+2')
                   temp_spec$intensity[idx_isotope] <- 0

                   # names(idx_isotope) <- c('M+1', 'M+2')
                 }

                 if (mono_int > temp_m_1_int & temp_m_1_int <= temp_m_2_int) {
                   idx_isotope <- c(temp_idx_isotope[[1]])
                   temp_spec$annotation[idx_isotope] <- c('M+1')
                   temp_spec$intensity[idx_isotope] <- 0
                   # names(idx_isotope) <- c('M+1', 'M+2')
                 }

               }

               if (length(temp_idx_isotope[[1]]) > 0 & length(temp_idx_isotope[[2]]) == 0) {
                 temp_m_1_int <- temp_spec$intensity[temp_idx_isotope[[1]]]

                 if (mono_int > temp_m_1_int) {
                   idx_isotope <- c(temp_idx_isotope[[1]])
                   temp_spec$annotation[idx_isotope] <- c('M+1')
                   temp_spec$intensity[idx_isotope] <- 0
                 }

               }


               if (all(temp_spec$intensity==0)) {
                 is_filter_isotope <- FALSE
               }

             }

             idx_remove <- which(temp_spec$annotation != 'M')
             spec <- spec[-idx_remove, , drop=F]
             spec
           }
)




#' @title removeRingEffect
#' @author Yandong Yin; Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param spec 2-column table. column 1: mz; column 2: intensity
#' @param mz_diff_thr the threold of mz difference. Default: 0.3
#' @param int_rel_thr the threold of relative intensity. Default: 0.2
#' @export

setGeneric(name = 'removeRingEffect',
           def = function(spec,
                          mz_diff_thr = 0.3,
                          int_rel_thr = 0.2){
             nr_ring <- nrow(spec) + 1
             mz <- spec[, 'mz']

             mz_diff <- diff(mz)
             idx_mzdiff <- which(mz_diff <= mz_diff_thr)
             if (length(idx_mzdiff) == 0) {
               return(spec)
             }

             nr_ring_possible <- unique(c(idx_mzdiff, idx_mzdiff + 1))

             # remove ringeffect loop
             while (TRUE) {

               idx_int_max <- which.max(spec[nr_ring_possible, 2])
               nr_int_max <- nr_ring_possible[idx_int_max] # the index of possible Ringeffect ions with maxium intensity
               int_thr <- spec[nr_int_max, 2] * int_rel_thr # the threshold = 0_2*max_int (possible ring)

               mz_diff <- abs(mz[nr_ring_possible[-idx_int_max]] - mz[nr_int_max])
               int <- spec[nr_ring_possible[-idx_int_max], 2]
               nr_ring <- append(nr_ring, nr_ring_possible[-idx_int_max][which(mz_diff <= mz_diff_thr & int <= int_thr)])
               nr_ring_possible <- nr_ring_possible[!nr_ring_possible %in% c(nr_ring, nr_int_max)]
               if (length(nr_ring_possible) == 0) {
                 break # break loop untill satisfy the nr_ring_possible==0
               }
             }

             return(spec[-nr_ring, , drop = FALSE])
           }
)


################################################################################
# combine_ms1_ms2 --------------------------------------------------------------
#' @title combine_ms1_ms2
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param ms1_data The name of ms1 peak table. Column 1 is "name", Column 2 is "mz" and column 3 is "rt".
#' @param ms2_data The vector of names of ms2 files. MS2 file must be mzXML.
#' @param mz_tol_combine_ms1_ms2 mz tol for ms1 and ms2 data matching.
#' @param rt_tol_combine_ms1_ms2 RT tol for ms1 and ms2 data matching.
#' @param ms2_type The type of MS2 file, default is mzXML.
#' @return Return ms1 and ms2 data.

setGeneric(name = "combine_ms1_ms2",
           function(ms1_data,
                    ms2_data,
                    mz_tol_combine_ms1_ms2 = 15,
                    rt_tol_combine_ms1_ms2 = 15,
                    ms2_type = c("mgf", "mzXML", 'mzML')){

             ms2_info <- ms2_data$info
             ms2_spec <- ms2_data$spec

             if (ms2_type %in% c("mzXML", 'mzML', 'mgf')) {
               cat("\n"); cat("Match MS1 and MS2 data.(mz tolerance ",
                              mz_tol_combine_ms1_ms2," ppm, RT tolerance ",
                              rt_tol_combine_ms1_ms2, " second).\n",
                              sep = "")

               temp_ms1_info <- ms1_data$info %>% select(mz, rt)
               match_result <- SXTMTmatch(data1 = ms2_info,
                                          data2 = temp_ms1_info,
                                          mz.tol = mz_tol_combine_ms1_ms2,
                                          rt.tol = rt_tol_combine_ms1_ms2,
                                          rt.error.type = "abs")

               match_result <- match_result %>%
                 tibble::as_tibble() %>%
                 dplyr::mutate(peak_name = ms1_data$info$name[Index2])

               # remove duplicated MS2 spectrum, if one peak has more than 1 ms2 spectrum,
               #    select the biggest of the sum intensity of top 10 fragments.

               unique_name <- match_result$peak_name %>% unique()
               cat("\n")
               cat("Select the most intense MS2 spectrum for one peak.\n")
               remain_idx <- pbapply::pblapply(seq_along(unique_name), function(i){
                 name <- unique_name[i]
                 temp_idx <- match_result %>% dplyr::filter(peak_name == name) %>% dplyr::pull(Index1)
                 # temp_idx <- which(ms2_name == name)
                 if(length(temp_idx) == 1) return(temp_idx)
                 temp_ms2 <- ms2_spec[temp_idx]
                 # temp_ms2 <- lapply(temp_ms2, function(x) x[[2]])
                 temp_int <- lapply(temp_ms2, function(x) {
                   x <- as.data.frame(x)
                   x <- x[order(x[,2], decreasing = TRUE),]
                   if(nrow(x) >= 10) {sum(x[1:10,2])}else{sum(x[,2])}
                 })
                 temp_int <- unlist(temp_int)
                 temp_idx <- temp_idx[which.max(temp_int)]
                 temp_idx
               })

               names(remain_idx) <- unique_name
               remain_idx <- unlist(remain_idx)

               ms2_name <- names(remain_idx)
               ms2_info <- ms2_info[remain_idx,]
               ms2_spec <- ms2_spec[remain_idx]
               names(ms2_spec) <- ms2_name

               ms2_list <- list(info = ms2_info, spec = ms2_spec)
               return(ms2_list)
             }

             if(ms2_type == "msp"){
               cat("\n")
               cat("Match MS1 and MS2 data according to the name.\n")

               # remove spectra not in peak table
               ms1_info <- ms1_data$info
               idx <- match(ms2_info$name, ms1_info$name) %>% is.na() %>% which()

               if (length(idx) > 0) {
                 ms2_info <- ms2_info[-idx,,drop = FALSE]
                 ms2_spec <- ms2_spec[-idx]
               }

               idx <- match(ms2_info$name, ms1_info$name)
               ms2_name <- ms2_info$name
               ms2_mz <- ms1_info$mz[idx]
               ms2_rt <- ms1_info$rt[idx]
               ms2_info <- data.frame(mz = ms2_mz,
                                      rt = ms2_rt,
                                      filename = ms2_name,
                                      stringsAsFactors = FALSE)
               names(ms2_spec) <- ms2_name

               ms2_list <- list(info = ms2_info, spec = ms2_spec)
               return(ms2_list)

               # # # change ms2 info to MetDNA1 default format
               # ms2 <- lapply(seq_along(ms2_name), function(i){
               #   temp_info <- data.frame('NAME' = ms2_name[i],
               #                           'PRECURSORMZ' = ms2_mz[i],
               #                           'PRECURSORRT' = ms2_rt[i],
               #                           stringsAsFactors = FALSE)
               #   temp_info <- t(temp_info)
               #   rownames(temp_info) <- c('NAME', 'PRECURSORMZ', 'PRECURSORRT')
               #   result <- list(info = temp_info,
               #                  spec = ms2_spec[[i]])
               #
               #   return(result)
               # })
               #
               # names(ms2) <- ms2_name
               # # save(ms2, file = file.path(path, "ms2"), compress = "xz", version = 2)
               #
               # return(ms2)
             }

           })






  # split_ms2 ------------------------------------------------------------------

setGeneric(name = 'split_ms2',
           def = function(ms2_data_combined){
             ms2_name <- ms2_data_combined$spec %>% names()
             ms2 <- pbapply::pblapply(seq_along(ms2_name), function(i){
               NAME <- ms2_name[i]
               PRECURSORMZ <- ms2_data_combined$info$mz[i]
               PRECURSORRT <- ms2_data_combined$info$rt[i]
               result_info <- data.frame(NAME, PRECURSORMZ,
                                         PRECURSORRT, stringsAsFactors = FALSE)
               result_info <- t(result_info)
               rownames(result_info) <- c('NAME', 'PRECURSORMZ', 'PRECURSORRT')

               result <- list(info = result_info,
                              spec = ms2_data_combined$spec[[i]])

               return(result)
             })

             names(ms2) <- ms2_name
             return(ms2)
           })



  # add_has_ms2_to_SpecAnnotationClass
#' @title add_has_ms2_to_SpecAnnotationClass
#' @author Zhiwei Zhou
#' @param obj_class a pecAnnotationClass object
#' @param ms2
#' @export

setGeneric(name = "add_has_ms2_to_SpecAnnotationClass",
           def = function(obj_class, ms2){
             peaks_with_ms2 <- names(ms2)
             result <- lapply(obj_class, function(x){
               if (x@peak_info$name %in% peaks_with_ms2) {
                 x@peak_info$has_ms2 <- 1
               } else {
                 x@peak_info$has_ms2 <- -1
               }

               return(x)
             })

             names(result) <- names(obj_class)
             return(result)
           })

################################################################################
# match_ms2 --------------------------------------------------------------------
#' @title match_ms2
#' @author Zhiwei Zhou
#' @param ms1_result ms1 result
#' @param exp_spec correspond msms data
#' @param lib_meta library of meta information
#' @param lib_spec library of spectrum
#' @param mz_tol_ms2 Default: 35 ppm
#' @param dp_cutoff Default: 0.8
#' @param matched_frag_cutoff Default: 1
#' @param direction 'forward', 'reverse'
#' @param scoring_approach 'dp', 'bonanza', 'hybrid', 'gnps'. Default: 'dp'
#' @export

# test <- matchMs2WithSpecLib(ms1_result = ms1_result,
#                             exp_spec = exp_spec,
#                             lib_meta = lib_meta,
#                             lib_spec = lib_spec,
#                             mz_tol_ms2 = 35,
#                             dp_cutoff = 0.8,
#                             direction = 'reverse',
#                             path = path)

setGeneric(name = 'match_ms2',
           def = function(
    ms1_result,
    exp_spec,
    lib_meta,
    lib_spec,
    mz_tol_ms2 = 35,
    dp_cutoff = 0.8,
    matched_frag_cutoff = 1,
    path = '.',
    is_include_precursor = TRUE,
    direction = c('reverse', 'forward'),
    scoring_approach = c('dp', 'bonanza', 'hybrid', 'gnps'),
    ...
           ){
             cat('Match ms2 with Spec library...\n')

             direction <- match.arg(direction)
             scoring_approach <- match.arg(scoring_approach)


             if (!('ms2_result' %in% list.files(file.path(path, "00_intermediate_data")))) {
               pbapply::pboptions(type='timer', char='+')
               ms2_result <- pbapply::pblapply(seq_along(ms1_result), function(i){
                 # modify the experimental ms2
                 #    if there have the exp ms2 & the candidate have the lib ms2, then run match
                 if (length(exp_spec[[i]]) > 0) {
                   if (nrow(ms1_result[[i]]@annotation_result) > 0) {
                     # modify the msms format for SpectraTools
                     temp_peak_name <- ms1_result[[i]]@peak_info$name
                     temp_peak_mz <- ms1_result[[i]]@peak_info$mz
                     temp_peak_spec <- exp_spec[[i]]

                     temp_peak_data <- list(info = tibble::tibble(NAME = temp_peak_name,
                                                                  PRECURSORMZ = temp_peak_mz),
                                            spec = temp_peak_spec)
                     exp_ms2 <- convertSpectraData(ms2_data = temp_peak_data)

                     # modify the library ms2
                     temp_lib_info <- ms1_result[[i]]@annotation_result %>%
                       dplyr::select(id, mz_lib) %>%
                       dplyr::rename(name = id,
                                     mz = mz_lib)

                     temp_lib_spec <- match(temp_lib_info$name, names(lib_spec)) %>% lib_spec[.] # if the candidate don't have the ms2, the spec is NULL
                     names(temp_lib_spec) <- c(temp_lib_info$name)

                     lib_ms2 <- new('SpectraData',
                                    info = temp_lib_info,
                                    spectra = temp_lib_spec)

                     result <- try(runSpecMatch(obj_ms2_cpd1 = exp_ms2,
                                                obj_ms2_cpd2 = lib_ms2,
                                                mz_tol_ms2 = mz_tol_ms2,
                                                scoring_approach = scoring_approach),
                                   silent = TRUE)


                     if ((class(result) == 'try-error') | (length(result) == 0)) {
                       return(NULL)
                     }

                     # result@info[is.na(result@info)] <- 0
                     idx_zero <- which(result@info$n_frag_cpd2 < 1)
                     result@info$scoreReverse[idx_zero] <- -1
                     result@info$scoreForward[idx_zero] <- -1

                     return(result)
                   } else {
                     return(NULL)
                   }
                 } else {
                   return(NULL)
                 }
               })

               names(ms2_result) <- names(ms1_result)

               dir.create(file.path(path, "00_intermediate_data"),
                          recursive = TRUE, showWarnings = FALSE)

               save(ms2_result,
                    file = file.path(path, "00_intermediate_data", 'ms2_result'),
                    compress = 'gzip',
                    version = 2)
             } else {
               cat('Note: load the ms2_result\n')
               load(file.path(path, '00_intermediate_data', 'ms2_result'))
             }



             # for DP score, record ms2 score as its direction
             # for other score, record ms2 score as reverse
             cat('Add the ms2 score into SpecAnnotationClass...\n')
             ms2_result_annotation <- mapply(function(x, y){
               if (nrow(x@annotation_result) > 0) {
                 if (scoring_approach == 'dp') {
                   if(length(y) > 0) {
                     temp <- y@info

                     # if (lib == 'gnpsLib') {
                     #   # if one metabolites have multiple spec, select the highest scoreReverse as its result
                     #   unique_name <- unique(temp$name)
                     #   idx <- lapply(unique_name, function(z){
                     #     idx <- which(temp$name == z)
                     #     if (length(idx) > 1) {
                     #       idx <- idx[which.max(temp$scoreReverse[idx])]
                     #     }
                     #     idx
                     #   }) %>% unlist()
                     #   temp <- temp[idx,]
                     # }

                     temp_forward <- temp$scoreForward
                     temp_reverse <- temp$scoreReverse
                     temp_n_frag <- temp$n_frag_total
                   } else {
                     temp_forward <- temp_reverse <- temp_n_frag <- 0
                   }
                   x@annotation_result$msms_score_forward <- temp_forward
                   x@annotation_result$msms_score_reverse <- temp_reverse
                   x@annotation_result$msms_matched_frag <- temp_n_frag
                   return(x)

                 } else {
                   if(length(y) > 0) {

                     # if (lib == 'gnpsLib') {
                     #   # if one metabolites have multiple spec, select the highest scoreReverse as its result
                     #   unique_name <- unique(temp$name)
                     #   idx <- lapply(unique_name, function(z){
                     #     idx <- which(temp$name == z)
                     #     if (length(idx) > 1) {
                     #       idx <- idx[which.max(temp$score[idx])]
                     #     }
                     #     idx
                     #   }) %>% unlist()
                     #   temp <- temp[idx,]
                     # }

                     temp <- y@info
                     temp_forward <- 0
                     temp_reverse <- temp$score
                     temp_n_frag <- temp$n_frag_total
                   } else {
                     temp_forward <- temp_reverse <- temp_n_frag <- 0
                   }
                   x@annotation_result$msms_score_forward <- temp_forward
                   x@annotation_result$msms_score_reverse <- temp_reverse
                   x@annotation_result$msms_matched_frag <- temp_n_frag
                   return(x)
                 }
               }

               return(x)

             },
             x = ms1_result,
             y = ms2_result,
             SIMPLIFY = FALSE)

             ms2_result_annotation <- pbapply::pblapply(ms2_result_annotation, function(x){
               if (nrow(x@annotation_result) > 0) {
                 if (scoring_approach == 'dp') {
                   switch (direction,
                           'reverse' = {
                             x@annotation_result <- x@annotation_result %>%
                               dplyr::filter(msms_score_reverse >= dp_cutoff & msms_matched_frag >= matched_frag_cutoff) %>%
                               dplyr::mutate(msms_score_forward = round(msms_score_forward, 4),
                                             msms_score_reverse = round(msms_score_reverse, 4))
                           },
                           'forward' = {
                             x@annotation_result <- x@annotation_result %>%
                               dplyr::filter(msms_score_forward >= dp_cutoff & msms_matched_frag >= matched_frag_cutoff) %>%
                               dplyr::mutate(msms_score_forward = round(msms_score_forward, 4),
                                             msms_score_reverse = round(msms_score_reverse, 4))
                           }
                   )
                 } else {
                   # for other score, the score was assigned as reverse score
                   x@annotation_result <- x@annotation_result %>%
                     dplyr::filter(msms_score_reverse >= dp_cutoff & msms_matched_frag >= matched_frag_cutoff) %>%
                     dplyr::mutate(msms_score_forward = round(msms_score_forward, 4),
                                   msms_score_reverse = round(msms_score_reverse, 4))
                 }

               }
               return(x)
             })


             return(ms2_result_annotation)
           }
)


#   runSpecMatch ---------------------------------------------------------------
#' @title runSpecMatch
#' @description a interphace of runing SpectraTools
#' @author Zhiwei Zhou
#' @param obj_ms2_cpd1 experimental ms2 object
#' @param obj_ms2_cpd2 library ms2 object
#' @param mz_tol_ms2 Default: 35 ppm
#' @param scoring_approach 'dp', 'bonanza', 'hybrid', 'gnps'
#' @export

# load('/home/zhouzw/Data_processing/20210224_metdna2_development_test/obj_ms2_cpd1_for_metdna2.RData')
# load('/home/zhouzw/Data_processing/20210224_metdna2_development_test/obj_ms2_cpd2_for_metdna2.RData')
# score_dp <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'dp')
# score_bonanza <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'bonanza')
# score_hybrid <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'hybrid')
# score_gnps <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'gnps')

setGeneric(name = 'runSpecMatch',
           def = function(
    obj_ms2_cpd1,
    obj_ms2_cpd2,
    mz_tol_ms2 = 35,
    scoring_approach = c('dp', 'bonanza', 'hybrid', 'gnps'),
    ...
           ){

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

                 if (is.null(obj_ms2_cpd2@spectra[[i]])) {
                   n_frag_cpd2 <- 0
                 } else {
                   n_frag_cpd2 <- nrow(obj_ms2_cpd2@spectra[[i]])
                 }

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
           })



#   convertSpectraData ---------------------------------------------------------

#' @title convertSpectraData
#' @param ms2_data
#' @importClassesFrom SpectraTools 'SpectraData'
#' @export

setGeneric(name = 'convertSpectraData',
           def = function(
    ms2_data
           ){
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
           })





################################################################################
# convertSpecAnnotationClass2Table ---------------------------------------------

#' @title convertSpecAnnotationClass2Table
#' @param ms1_data
#' @param result_annotation
#' @param lib_meta
#' @param instrument
#' @param is_rt_score
#' @param is_ccs_score
#' @param is_msms_score
#' @param direction
#' @export

setGeneric(name = 'convertSpecAnnotationClass2Table',
           def = function(
    ms1_data,
    result_annotation,
    ...
           ){
             cat('\n');cat('Generate metabolite annotation report\n')

             pbapply::pboptions(type='timer', char='+')
             report_table <- pbapply::pblapply(seq_along(result_annotation), function(i){
               # cat(i, ' ')
               x <- result_annotation[[i]]

               if (nrow(x@annotation_result) > 0) {
                 result <- x@annotation_result %>%
                   dplyr::arrange(desc(msms_score_forward),
                                  desc(msms_score_reverse),
                                  desc(rt_score)) %>%
                   dplyr::mutate(feature = x@peak_info$name) %>%
                   dplyr::select(feature, dplyr::everything())
                 return(result)
               } else {
                 return(NULL)
               }

             })

             report_table <- report_table %>% dplyr::bind_rows()

             # browser()

             if (length(report_table) == 0) {
               report_result <- ms1_data$info %>%
                 dplyr::bind_cols(ms1_data$subject) %>%
                 dplyr::mutate(id_reverse_summary = NA,
                               id_forward_summary = NA)

               return(report_result)
             }

             id_summary_forward <- report_table %>%
               dplyr::mutate(id = paste0('score{', msms_score_forward, '}',
                                         'frag{', msms_matched_frag, '}',
                                         'adduct{', adduct, '}',
                                         'name{', name, '}',
                                         'labid{', id, '}')) %>%
               dplyr::group_by(feature) %>%
               dplyr::summarise(id_forward_summary = paste(id, collapse = ';')) %>%
               dplyr::ungroup()

             id_summary_reverse <- report_table %>%
               dplyr::mutate(id = paste0('score{', msms_score_reverse, '}',
                                         'frag{', msms_matched_frag, '}',
                                         'adduct{', adduct, '}',
                                         'name{', name, '}',
                                         'labid{', id, '}')) %>%
               dplyr::group_by(feature) %>%
               dplyr::summarise(id_reverse_summary = paste(id, collapse = ';')) %>%
               dplyr::ungroup()

             report_result <- ms1_data$info %>%
               # dplyr::bind_cols(ms1_data$subject) %>%
               dplyr::left_join(id_summary_reverse, by = c('name' = 'feature')) %>%
               dplyr::left_join(id_summary_forward, by = c('name' = 'feature'))

             return(report_result)

           })



################################################################################
# generate_summary_table -------------------------------------------------------

#' @title generate_summary_table
#' @description generate a summary table of metabolite annotation
#' @author Zhiwei Zhou
#' @param result_annotation
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract

# load('~/Project/00_IBD_project/Data/20230227_develop_metabolite_ID_workflow/01_metabolite_annotation_ms2/00_intermediate_data/result_annotation')
# load('~/Project/00_IBD_project/Data/20230227_develop_metabolite_ID_workflow/01_metabolite_annotation_ms2/00_intermediate_data/lib_meta')
# load('~/Project/00_IBD_project/Data/20230227_develop_metabolite_ID_workflow/01_metabolite_annotation_ms2/00_intermediate_data/ms1_data')


setGeneric(name = 'generate_summary_table',
           def = function(result_annotation){

             idx_id <- sapply(result_annotation, function(x){
               x@annotation_result %>% nrow()
             })

             if (any(idx_id> 0)) {
               annot_table <- lapply(result_annotation, function(x){
                 if (nrow(x@annotation_result) > 0) {
                   temp_feature_name <- x@peak_info$name
                   temp_feature_mz <- x@peak_info$mz
                   temp_feature_rt <- x@peak_info$rt
                   temp_feature_has_ms2 <- x@peak_info$has_ms2

                   result <- x@annotation_result %>%
                     dplyr::mutate(feature_name = temp_feature_name,
                                   mz = temp_feature_mz,
                                   rt = temp_feature_rt,
                                   with_ms2 = temp_feature_has_ms2) %>%
                     dplyr::select(feature_name:with_ms2, dplyr::everything()) %>%
                     dplyr::select(-idx)

                   return(result)
                 } else {
                   return(NULL)
                 }
               }) %>% bind_rows()

             } else {
               annot_table <- result_annotation[[1]]@annotation_result %>%
                 dplyr::mutate(feature_name = character(),
                               mz = numeric(),
                               rt = numeric(),
                               with_ms2 = numeric()) %>%
                 dplyr::select(feature_name:with_ms2, dplyr::everything()) %>%
                 dplyr::select(-idx)
             }

             return(annot_table)
           })




