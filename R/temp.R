# ms2_spec <- get_ms2(lab_id = 'I00092', ce = '20', polarity = 'positive')
# plot_ms2(ms2_spec, is_normalize = FALSE)
#
# ms2_spec <- get_ms2(lab_id = 'I00634', ce = '40', polarity = 'negative')
# plot_ms2(ms2_spec, is_normalize = FALSE)
#
# get_compound(lab_id = 'I00092')
# get_compound(lab_id = 'I00089')
# get_compound(lab_id = 'I00634')



# export_db_msdial(lib = 'stanford',
#                  column = 'hilic',
#                  polarity = 'positive',
#                  msp_name = 'stanford_lib_hilic_pos_230227.msp',
#                  dir_path = '~/Project/00_IBD_project/Data/20230215_MS2_lib_processing/04_final_tables')
#
#
# export_db_msdial(lib = 'stanford',
#                  column = 'hilic',
#                  polarity = 'negative',
#                  msp_name = 'stanford_lib_hilic_neg_230227.msp',
#                  dir_path = '~/Project/00_IBD_project/Data/20230215_MS2_lib_processing/04_final_tables')
#
#
# export_db_msdial(lib = 'stanford',
#                  column = 'c18',
#                  polarity = 'positive',
#                  msp_name = 'stanford_lib_c18_pos_230227.msp',
#                  dir_path = '~/Project/00_IBD_project/Data/20230215_MS2_lib_processing/04_final_tables')
#
#
# export_db_msdial(lib = 'stanford',
#                  column = 'c18',
#                  polarity = 'negative',
#                  msp_name = 'stanford_lib_c18_neg_230227.msp',
#                  dir_path = '~/Project/00_IBD_project/Data/20230215_MS2_lib_processing/04_final_tables')
#
#
# export_db_msdial(lib = 'iroa',
#                  column = 'hilic',
#                  polarity = 'positive',
#                  msp_name = 'iroa_lib_hilic_pos_230227.msp',
#                  dir_path = '~/Project/00_IBD_project/Data/20230215_MS2_lib_processing/04_final_tables')
#
#
# export_db_msdial(lib = 'iroa',
#                  column = 'hilic',
#                  polarity = 'negative',
#                  msp_name = 'iroa_lib_hilic_neg_230227.msp',
#                  dir_path = '~/Project/00_IBD_project/Data/20230215_MS2_lib_processing/04_final_tables')
#
#
# export_db_msdial(lib = 'iroa',
#                  column = 'c18',
#                  polarity = 'positive',
#                  msp_name = 'iroa_lib_c18_pos_230227.msp',
#                  dir_path = '~/Project/00_IBD_project/Data/20230215_MS2_lib_processing/04_final_tables')
#
#
# export_db_msdial(lib = 'iroa',
#                  column = 'c18',
#                  polarity = 'negative',
#                  msp_name = 'iroa_lib_c18_neg_230227.msp',
#                  dir_path = '~/Project/00_IBD_project/Data/20230215_MS2_lib_processing/04_final_tables')


# export_db_tidymass(lib = )

# load('~/Project/00_IBD_project/Data/20230227_develop_metabolite_ID_workflow/05_object_hilic_pos_outlier_removal.RData')
# load('~/Downloads/01_result_initial_seed_annotation/00_intermediate_data/lib_meta')
# load('~/Downloads/01_result_initial_seed_annotation/00_intermediate_data/lib_spec')
# load('~/Downloads/01_result_initial_seed_annotation/00_intermediate_data/ms1_data')
# load('~/Downloads/01_result_initial_seed_annotation/00_intermediate_data/ms2')
#
# test <- sapply(ms1_result, function(x){
#   nrow(x@annotation_result)
# })
#
# test_idx <- which(test > 0)
#
# result_annotation[test_idx]
#
# sum(unlist())
#
# lapply(object@ms2_data, function(x) {
#   length(unique(x@ms2_spectra))
# })


# whole workflows
# load('~/Project/00_IBD_project/Data/20230227_develop_metabolite_ID_workflow/05_object_hilic_pos_outlier_removal.RData')
# annotate_metabolite(object = object_hilic_pos,
#                     ms2_type = 'mzML',
#                     path = '~/Project/00_IBD_project/Data/20230227_develop_metabolite_ID_workflow',
#                     polarity = 'positive',
#                     lib = 'dodd',
#                     column = 'hilic',
#                     ce = '20',
#                     adduct_list = '[M+H]+',
#                     mz_tol = 15,
#                     mz_ppm_thr = 150,
#                     pf_rt_range = 10,
#                     tolerance_rt_range = 20,
#                     is_rt_score = TRUE,
#                     is_ms2_score = TRUE)
#
#
# load('~/Project/00_IBD_project/Data/20230227_develop_metabolite_ID_workflow/01_metabolite_annotation_ms2/00_intermediate_data/ms2')

# load('~/Project/00_IBD_project/Data/20230227_develop_metabolite_ID_workflow/01_metabolite_annotation/00_intermediate_data/ms1_result')
# load('~/Project/00_IBD_project/Data/20230227_develop_metabolite_ID_workflow/01_metabolite_annotation/00_intermediate_data/ms1_data')
# load('~/Project/00_IBD_project/Data/20230227_develop_metabolite_ID_workflow/01_metabolite_annotation/00_intermediate_data/result_annotation')
#
# load('~/Project/00_IBD_project/Data/20230227_develop_metabolite_ID_workflow/05_object_hilic_pos_outlier_removal.RData')
# annotate_metabolite(object = object_hilic_pos,
#                     ms2_type = 'mzML',
#                     path = '~/Project/00_IBD_project/Data/20230227_develop_metabolite_ID_workflow',
#                     polarity = 'positive',
#                     lib = 'dodd',
#                     column = 'hilic',
#                     ce = '20',
#                     adduct_list = '[M+H]+',
#                     mz_tol = 15,
#                     mz_ppm_thr = 150,
#                     pf_rt_range = 10,
#                     tolerance_rt_range = 20,
#                     is_rt_score = TRUE,
#                     is_ms2_score = FALSE)
#
# load('~/Project/00_IBD_project/Data/20230301_annotate_metabolites/HILIC_pos/01_metabolite_annotation_mz_rt/00_intermediate_data/ms1_result')
# load('~/Project/00_IBD_project/Data/20230301_annotate_metabolites/HILIC_pos/01_metabolite_annotation_mz_rt_ms2/00_intermediate_data/ms2')
# load('~/Project/00_IBD_project/Data/20230301_annotate_metabolites/DoddLib/HILIC_pos/01_metabolite_annotation_mz_rt/00_intermediate_data/lib_meta')
# load('~/Project/00_IBD_project/Data/20230301_annotate_metabolites/DoddLib/HILIC_pos/01_metabolite_annotation_mz_rt/00_intermediate_data/lib_spec')

# mz_tol = 15
# mz_ppm_thr = 150
# pf_rt_range = 10
# tolerance_rt_range = 20
# is_rt_score = TRUE
# is_ms2_score = FALSE


# load('~/Project/00_IBD_project/Data/20230227_develop_metabolite_ID_workflow/05_object_hilic_pos_outlier_removal.RData')
# load('~/Project/04_package/00_Database/MSDIAL_library/msdial_db_pos_230302.RData')
# annotate_metabolite(object = object_hilic_pos,
#                     ms2_type = 'mzML',
#                     path = '~/Project/00_IBD_project/Data/20230227_develop_metabolite_ID_workflow',
#                     polarity = 'positive',
#                     lib_db = msdial_db_pos,
#                     mz_tol = 15,
#                     mz_ppm_thr = 150,
#                     pf_rt_range = 10,
#                     tolerance_rt_range = 20,
#                     dp_cutoff = 0.8,
#                     is_rt_score = FALSE,
#                     is_ms2_score = TRUE)

# path_output <- '~/Project/00_IBD_project/Data/20230227_develop_metabolite_ID_workflow/01_metabolite_annotation/'
# load('~/Project/00_IBD_project/Data/20230227_develop_metabolite_ID_workflow/01_metabolite_annotation/00_intermediate_data/result_annotation')
# load('~/Project/00_IBD_project/Data/20230227_develop_metabolite_ID_workflow/01_metabolite_annotation/00_intermediate_data/table_annotation')


################################################################################
# # 20230809 ---------------------------------------------------------------------
# load('~/Project/04_package/00_Database/DoddLib/03_isotope_standard_ISTD/v230215/istd_lib_230215.RData')
# istd_lib_v1 <- istd_lib
#
# load('~/Project/04_package/00_Database/DoddLib/03_isotope_standard_ISTD/v230808/istd_lib_v2_230809.RData')
# istd_lib_v2 <- istd_lib
#
# istd_lib <- list('v1' = istd_lib_v1, 'v2' = istd_lib_v2)
#
# usethis::use_data(istd_lib, overwrite = TRUE)


################################################################################
# 20230816 biomet --------------------------------------------------------------
#
#
# biomet_lib <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20230808_data_quality_check/biomet_lib_230816.xlsx')
#
#
# lib_info <- list('db_name' = 'biomet_lib',
#                  'version' = 'v1.0.0',
#                  'date' = Sys.time(),
#                  'submitter' = 'Zhiwei Zhou')
#
# biomet_lib <- list(lib_info = lib_info,
#                    'lib_table' = biomet_lib)
#
# save(biomet_lib, file = '~/Project/00_IBD_project/Data/20230808_data_quality_check/biomet_lib_230816.RData')
#
#
# biomet_lib <- list('v1' = biomet_lib)
#
# dir.create('~/Project/04_package/00_Database/DoddLib/04_biomet/230816', showWarnings = FALSE, recursive = TRUE)
# save(biomet_lib,
#      file = '~/Project/04_package/00_Database/DoddLib/04_biomet/230816/biomet_lib_230816.RData')
#
# usethis::use_data(biomet_lib, overwrite = TRUE)
#
# # update istd
# load('~/Project/04_package/00_Database/DoddLib/03_isotope_standard_ISTD/v230215/istd_lib_230215.RData')
# istd_lib_v1 <- istd_lib
#
# load('~/Project/04_package/00_Database/DoddLib/03_isotope_standard_ISTD/v230808/istd_lib_v2_230816.RData')
# istd_lib_v2 <- istd_lib
#
# istd_lib <- list('v1' = istd_lib_v1, 'v2' = istd_lib_v2)
#
# usethis::use_data(istd_lib, overwrite = TRUE)
