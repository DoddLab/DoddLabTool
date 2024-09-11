################################################################################
# startup massage --------------------------------------------------------------
.onAttach <- function(libname, pkgname){
  packageStartupMessage("
Version 0.1.5
-------------
Authors: Zhiwei Zhou
Maintainer: Zhiwei Zhou

Updates
-------------
o Add functions (modify_sample_info, modify_worklist, merge_worklist_sample_info, create_ibd_object) for IBD project
")
}
