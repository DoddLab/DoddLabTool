################################################################################
# modify_sample_info -----------------------------------------------------------

#' @title modify_sample_info
#' @description modify the sample information for data processing (Aedan extraction table)
#' @author Zhiwei Zhou
#' @param file_name
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' setwd('~/Project/00_IBD_project/Data/00_IBD_worklists & sample_id/')
#' file_name <- '~/Project/00_IBD_project/Data/00_IBD_worklists & sample_id/sampleinfo/2023-08-25 IBD enrollment, NewBatch_004 samples.csv'
#' modify_sample_info(file_name)
#' }


# setwd('~/Project/00_IBD_project/Data/00_IBD_worklists & sample_id/')
# file_name <- './sampleinfo/2023-08-25 IBD enrollment, NewBatch_004 samples.csv'
# modify_sample_info(file_name)

modify_sample_info <- function(file_name){
  raw_data <- readr::read_csv(file_name)
  batch_number <- file_name %>% stringr::str_extract('Batch_\\d+|batch_\\d+') %>% stringr::str_extract('\\d+')
  batch_number <- paste0('B', batch_number)

  result <- raw_data %>%
    dplyr::rename('sample_name' = '#',
                  'originating_id' = 'Originating Id',
                  'deidentified_master_patient_id' = 'DEIDENTIFIED_MASTER_PATIENT_ID',
                  'phenotype_group' = 'Phenotype Group',
                  'visit_number' = 'Visit Number',
                  'patient_ids' = 'Samples',
                  'plasma_location' = 'Plasma Location',
                  'ms_ids' = 'MS #',
                  'extract_location' = 'Extract location',
                  'note' = 'Note',
                  'extract_confirmation' = 'Extract Confrimation')


  temp_sample_id <- result$ms_ids %>% stringr::str_extract(pattern = 'S\\d+')

  result <- result %>%
    dplyr::mutate(sample_id = paste0(batch_number, '_', temp_sample_id)) %>%
    dplyr::select(sample_id, dplyr::everything()) %>%
    dplyr::mutate(phenotype_group2 = dplyr::case_when(phenotype_group == 'Non-IBD' ~ 'non_IBD',
                                                      phenotype_group != 'Non-IBD' ~ 'CD')) %>%
    dplyr::select(sample_id, originating_id:phenotype_group, phenotype_group2, visit_number:note)

  return(result)
}



################################################################################
# modify_worklist --------------------------------------------------------------

#' @title modify_worklist
#' @description modify the worklist for data processing (LC-MS worklist)
#' @author Zhiwei Zhou
#' @param file_name
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' setwd('~/Project/00_IBD_project/Data/00_IBD_worklists & sample_id/')
#' file_name <- '~/Project/00_IBD_project/Data/00_IBD_worklists & sample_id/worklists/B003_Worklist_c18_neg.xlsx'
#' modify_worklist(file_name)
#' }

# file_name <- './worklists/B003_Worklist_c18_neg.xlsx'
modify_worklist <- function(file_name){
  raw_data <- readxl::read_xlsx(file_name, col_names = FALSE)
  # batch_number <- file_name %>% stringr::str_extract('Batch_\\d+|batch_\\d+|B\\d+') %>% stringr::str_extract('\\d+')
  # batch_number <- paste0('B', batch_number)

  colnames(raw_data) <- c('sample_name', 'position', 'method', 'file_name', 'type')

  result <- raw_data %>%
    dplyr::select(1:5) %>%
    dplyr::filter(!is.na(file_name)) %>%
    dplyr::mutate(sample_type = dplyr::case_when(sample_name == 'Blank' ~ 'blank',
                                                 sample_name == 'PoolQC-equilibrate' ~ 'poolQC_equilibrate',
                                                 sample_name == "ProcedureBlank" ~ "procedureBlank",
                                                 sample_name == 'NIST-QC' ~ 'nistQC',
                                                 sample_name == 'PoolQC' ~ 'poolQC',
                                                 sample_name == 'PoolQC-DDA' ~ 'poolQC_dda',
                                                 stringr::str_detect(sample_name, 'S\\d+') ~ 'sample'))

  batch_id <- result$file_name %>% stringr::str_extract(pattern = 'B\\d+') %>% unique()
  mode <- result$method %>% stringr::str_extract(pattern = '\\\\Untargeted_.+') %>% stringr::str_extract(pattern = 'C18_Pos|C18_neg|HILIC_pos|HILIC_neg')
  column <- mode %>% stringr::str_extract(pattern = 'C18|HILIC')
  polarity <- mode %>% stringr::str_extract(pattern = 'Pos|Neg')
  new_sample_name <- result$file_name %>% stringr::str_replace(pattern = '\\.d', replacement = '')

  result <- result %>%
    dplyr::mutate(batch_id = batch_id, column = column, polarity = polarity) %>%
    dplyr::mutate(sample_name = new_sample_name) %>%
    dplyr::select(sample_name, sample_type:polarity) %>%
    dplyr::mutate(injection_order = seq(n()))

  return(result)
}



################################################################################
# generate one table for peak extraction ---------------------------------------

#' @title merge_worklist_sample_info
#' @description modify the worklist for data processing (LC-MS worklist)
#' @author Zhiwei Zhou
#' @param file_name
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' setwd('~/Project/00_IBD_project/Data/00_IBD_worklists & sample_id/')
#' file_name <- '~/Project/00_IBD_project/Data/00_IBD_worklists & sample_id/worklists/B003_Worklist_c18_neg.xlsx'
#' modify_worklist(file_name)
#' }

# for each mode, read and modify all sampleinfo and worklists
# library(tidyverse)
#
# setwd('~/Project/00_IBD_project/Data/00_IBD_worklists & sample_id/')
# sample_info_files <- list.files('./sampleinfo')
# worklist_files <- list.files('./worklists')

merge_worklist_sample_info <- function(sample_info_files,
                                       worklist_files,
                                       path = '.',
                                       mode = c('c18_pos', 'c18_neg', 'hilic_pos', 'hilic_neg')) {

  # modify sample_info table
  list_sample_info <- lapply(sample_info_files, function(x){
    modify_sample_info(file_name = file.path(path, 'sampleinfo', x))
  }) %>%
    dplyr::bind_rows()


  # modify worklist table
  temp_worklist_files <- worklist_files %>% stringr::str_detect(mode) %>% worklist_files[.]
  list_worklist <- lapply(temp_worklist_files, function(x){
    cat(x, ' ')
    modify_worklist(file_name = file.path(path, 'worklists', x))
  }) %>% dplyr::bind_rows()

  result_table <- list_worklist %>%
    dplyr::left_join(list_sample_info, by = c('sample_name' = 'sample_id')) %>%
    dplyr::mutate(sample_type2 = dplyr::case_when(
      sample_type == 'sample' ~ phenotype_group2,
      sample_type != 'sample' ~ sample_type
    )) %>%
    dplyr::select(sample_name:sample_type, sample_type2, dplyr::everything())

  # result_table %>%
  #   dplyr::filter(!(sample_type %in% c('blank', 'poolQC_equilibrate')))

  dir.create(file.path(path, 'merged_worklist'), showWarnings = FALSE, recursive = TRUE)
  writexl::write_xlsx(result_table, path = file.path(path, 'merged_worklist', paste0('merged_worklist_', mode, '.xlsx')), format_headers = FALSE)
}



################################################################################
# create_ibd_object ------------------------------------------------------------

#' @title create_ibd_object
#' @author Zhiwei Zhou
#' @param peak_table
#' @param merged_worklist
#' @param dir_path
#' @param mode c('c18_pos', 'c18_neg', 'hilic_pos', 'hilic_neg')
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @importClassesFrom massdataset mass_dataset
#' @export
#' @examples
#' \dontrun{
#' peak_table <- 'Peak-table.csv'
#' dir_path <- '~/Project/00_IBD_project/Data/20230906_B001_B004_batch_merge/C18_pos/00_raw_data_processing'
#' merged_worklist <- '~/Project/00_IBD_project/Data/00_IBD_worklists & sample_id/merged_worklist/merged_worklist_c18_pos.xlsx'
#' mode <- 'c18_pos'
#' create_ibd_object(peak_table = 'Peak-table.csv',
#'                   merged_worklist = '~/Project/00_IBD_project/Data/00_IBD_worklists & sample_id/merged_worklist/merged_worklist_c18_pos.xlsx',
#'                   dir_path = '~/Project/00_IBD_project/Data/20230906_B001_B004_batch_merge/C18_pos/00_raw_data_processing',
#'                   mode = 'c18_pos')
#' }



# peak_table <- 'Peak-table.csv'
# dir_path <- '~/Project/00_IBD_project/Data/20230906_B001_B004_batch_merge/C18_pos/00_raw_data_processing'
# merged_worklist <- '~/Project/00_IBD_project/Data/00_IBD_worklists & sample_id/merged_worklist/merged_worklist_c18_pos.xlsx'
# mode <- 'c18_pos'

# create_ibd_object(peak_table = 'Peak-table.csv',
#                   merged_worklist = '~/Project/00_IBD_project/Data/00_IBD_worklists & sample_id/merged_worklist/merged_worklist_c18_pos.xlsx',
#                   dir_path = '~/Project/00_IBD_project/Data/20230906_B001_B004_batch_merge/C18_pos/00_raw_data_processing',
#                   mode = 'c18_pos')


create_ibd_object <- function(peak_table = 'Peak-table.csv',
                              merged_worklist = 'merged_worklist_c18_pos.xlsx',
                              dir_path = '.',
                              mode = c('c18_pos', 'c18_neg', 'hilic_pos', 'hilic_neg')) {

  mode <- match.arg(mode)

  cat('Modify the worklist for tidymass worklist ...\n')
  # modify the worklist format for object
  worklist_data <- readxl::read_xlsx(merged_worklist)
  sample_info <- worklist_data %>%
    dplyr::filter(!(sample_type %in% c("blank", "poolQC_equilibrate", "procedureBlank", "poolQC_dda"))) %>%
    dplyr::rename(sample_id = sample_name,
                  injection.order = injection_order,
                  class = sample_type,
                  batch = batch_id,
                  subject.id = deidentified_master_patient_id,
                  group = sample_type2) %>%
    dplyr::select(sample_id, injection.order, class, batch, subject.id, group, everything()) %>%
    dplyr::mutate(batch = stringr::str_extract(batch, pattern = '\\d+') %>% as.numeric()) %>%
    dplyr::mutate(sample_id = stringr::str_replace(sample_id, pattern = '\\-', replacement = '_'))  %>%
    dplyr::group_by(batch) %>%
    dplyr::mutate(injection.order = seq(dplyr::n())) %>%
    dplyr::ungroup()

  # modify the peak table
  cat('Modify the expression_data & variable_info ...\n')
  sample_id <- sample_info$sample_id
  peak_table <- readr::read_csv(file.path(dir_path, peak_table), show_col_types = FALSE)
  colnames(peak_table) <- colnames(peak_table) %>%
    stringr::str_replace(pattern = '-', replacement = '_')

  peak_table <- peak_table %>%
    dplyr::select(name, mzmed, rtmed, sample_id) %>%
    dplyr::rename(mz = mzmed,
                  rt = rtmed,
                  variable_id = name) %>%
    dplyr::mutate(variable_id = paste0(variable_id, '_', mode))

  # expression_data_pos
  expression_data <- peak_table %>%
    dplyr::select(-c(variable_id:rt)) %>%
    as.data.frame()

  # variable_info_pos
  variable_info <- peak_table %>%
    dplyr::select(variable_id:rt) %>%
    as.data.frame()

  rownames(expression_data) <- variable_info$variable_id

  # check sample_id in sample_info whether same with expression_data, if not, reorganize
  colnames(expression_data) == sample_info$sample_id
  expression_data <- expression_data[, sample_info$sample_id]

  # creat mass_data object -------------------------------------------------------
  cat('Create the mass dataset object ...\n')
  object <- massdataset::create_mass_dataset(expression_data = expression_data,
                                             sample_info = sample_info,
                                             variable_info = variable_info)

  dir.create(path = file.path(dir_path, '01_input_data_cleaning'), showWarnings = FALSE, recursive = TRUE)
  save(object,
       file = file.path(dir_path, '/01_input_data_cleaning/object.RData'))

  massdataset::export_mass_dataset(object = object,
                                   file_type = "xlsx",
                                   path = file.path(dir_path, "/01_input_data_cleaning/mass_data_tables"))

  cat('Done!\n')
}



