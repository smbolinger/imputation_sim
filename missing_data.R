#############################################################################
#' List missing data values
#############################################################################
dat_miss <- function(ndGLM_){
  camFateNA    <- ndGLM_$nest[which(is.na(ndGLM_$cam_fate))]
  fdateNA  <- ndGLM_$nest[which(is.na(ndGLM_$fdate))]
  # fdateNA      <- ndGLM_$nest[which(is.na(ndGLM_$k))]
  # is fate date when fate was assigned, or when fate occurred?
  nest_ageNA   <- ndGLM_$nest[which(is.na(ndGLM_$nest_age))]

  # what is the difference between base R union() and dplyr union()?

  cfate_fdateNA <- union(camFateNA, fdateNA) # all NAs except nest age
  allNA <- union(cfate_fdateNA, nest_ageNA)  # all NAs
  justFateNA <- setdiff(cfate_fdateNA, camFateNA)
  justCamNA <- setdiff(cfate_fdateNA, fdateNA)
  numNA <- length(allNA)                     # how many nests are missing data?
  ndGLM_$numNA <- numNA                      # export the counts as part of the data
  # ndGLM1 <- ndGLM_ %>% dplyr::select(!c(.data$status))    # remove status (list of lists) & rename
  cat("Nests missing age (", length(nest_ageNA),"total):\n\n", nest_ageNA)
  # cat("\n\nNests missing camera fate (", length(camFateNA), "total):\n\n", camFateNA)
  cat("\n\nNests missing camera fate (", length(camFateNA), "total):\n\n", camFateNA)
  cat("\n\nNests missing k (", length(fdateNA),"total):\n\n", fdateNA)
  cat("\n\nTotal rows that would be removed:\n\n", numNA)
  cat("\n\nAll the missing values:\n\n")
  ndGLM1 %>%
    # dplyr::select(!c(fate, i, j, k, c.fate, numNA))%>%
    sapply(function(x) sum(is.na(x)))

  return(allNA)

}

#############################################################################
#' Create a table of missingness
#############################################################################
miss_tabs <- function(dName, expl, dep){
  # datDiag <- get(dName)

  ndDiag <- get(dName)
  # creeate an index of rows with at least one NA:
  exp_index <- lapply(expl, function(x) sum(is.na(ndDiag[x]))) > 0
  exp_miss <- expl[exp_index]

  tabNames <- list()
  for(v in seq_along(exp_miss)){
    missTab <- ndDiag %>%
      finalfit::missing_compare(dependent = expl[v],
                                # explanatory = expl[expl!=v])
                                explanatory = expl[-v])
    print(missTab)
    titl <- names(missTab)[1]
    names(missTab) <- c("var", "val", "not_missing", "missing", "p")
    missTabFancy <- missTab %>%
      gt::gt(rowname_col = "var") %>%
      gt::tab_header(title=titl) %>%
      gt::cols_label(val ~ "",
                 not_missing ~ "Not missing (percent)",
                 missing ~ "Missing (percent)",
                 p ~ "p-value") %>%
      # nestGLM::save_tab(suffix=paste(dName,dep,expl[v],sep="_"), rtf=TRUE)
      save_tab(suffix=paste(dName,dep,expl[v],sep="_"), rtf=TRUE)
    tabNames[v] <- missTab
  }
    # return(missTabFancy)
  return(tabNames)
}
#############################################################################

#' Remove all rows containing NAs
#`############################################################################
#'
#' @param ndGLM1 nest data; must include all predictor vars (nest_age, fate_date, cam_fate, species, & obs_int)
#' @param allNA list of nest numbers for all nests with missing values (NA)
#'
#' @return nest data with NAs removed
#' @export
#'
remove_na <- function(ndGLM1, allNA){
  ndGLM1$remove <- ndGLM1$nest %in% allNA # use %in% instead of ==
  true1 <- sum(ndGLM1$remove)

  # write.csv(ndGLM1, sprintf("nestdata_NAremoved_%s.csv", now))

  ndGLM1 <- ndGLM1[ndGLM1$remove == FALSE,]
  # in this subset we are not excluding nests with unknown age
  # AKA keep only rows that aren't flagged with remove == TRUE
  # fates <- ndGLM2 %>% select(nest,final_fate,cam_fate,misclass,how_mis,hatchfail,c_hatchfail,HF_mis, k)
  # table(ndGLM2$c.fate)
  cat("fates, all NAs removed:",table(ndGLM1$cam_fate))

  return(ndGLM1)
}

homeDir <- "C://Users/sarah/Dropbox/Models/fate_glm/"

missing_tab <- function(datName, prVars, dep = c("HF_mis", "is_u") ){
  ndDiag <- get(datName)
  expl = prVars
  
  # creeate an index of variables with at least one NA:
  exp_index <- lapply(expl, function(x) sum(is.na(ndDiag[x]))) > 0
  exp_miss <- expl[exp_index] 
  exp_miss
  
  tabNames <- list()
  # for(v in seq_along(exp_miss)){
  for(v in which(exp_index)){ # use the index for the vars withint expl
    # v=1
    # expl[v]
    # expl[-v]
    missTab <- ndDiag %>%
      dplyr::mutate(across(c(species, cam_fate), as.factor)) 
      
    # for this table, "dependent" is the independent var that is being compared
    missTab <- missTab %>%
      # finalfit::missing_compare(dependent = exp_miss[v], # this should be the table
      finalfit::missing_compare(dependent = expl[v], # this should be the table
      # finalfit::missing_pairs(dependent = expl[v],  # this is the plot
                                # explanatory = expl[expl!=v])
                                explanatory = expl[-v])
    # print(missTab)
    titl <- names(missTab)[1]
    names(missTab) <- c("var", "val", "not_missing", "missing", "p")
    
    missTabFancy <- missTab %>%
      gt::gt(rowname_col = "var") %>%
      gt::tab_header(title=titl) %>%
      gt::cols_label(val ~ "",
                 not_missing ~ "Not missing (percent)",
                 missing ~ "Missing (percent)",
                 p ~ "p-value") %>%
      # nestGLM::save_tab(suffix=paste(dName,dep,expl[v],uVar,sep="_"), rtf=TRUE)
      # nestGLM::save_tab(suffix=paste(datName,expl[v],uVar,sep="_"), dir = homeDir,rtf=TRUE)
      save_tab(suffix=paste(datName,expl[v],uVar,sep="_"), dir = homeDir,rtf=TRUE)
    
    tabNames[v] <- missTab
    # missPP <- list()
    for(i in seq_along(dep)){
      # missPP[i] <- ndDiag %>%
      missPP <- ndDiag %>%
        finalfit::missing_pairs(dependent = dep[i],
                                explanatory = expl,
                                position = "fill")
      now = format(Sys.time(), "_%m%d_%H%M_")
      # fname <- paste0(homeDir,"glm_script/figures/pp_",datName,dep[i],expl[v],uVar,now,".svg")
      fname <- paste0(homeDir,"figures/pp_",datName,dep[i],expl[v],uVar,now,".svg")
      ggsave(plot=missPP, filename = fname, device = "svg")
    }
    
    return(missTabFancy) # this returns the filename?
    # now = format(Sys.time(), "_%m%d_%H%M_")
    # fname <- paste0(homeDir,"glm_script/figures/pp_",datName,dep[i],expl[v],uVar,now,".svg")
    # ggsave(plot=missPP, filename = fname)
  }
}


#  > visdat:::vis_create_
vis_create_ <- function (x) 
{
  ggplot2::ggplot(data = x, ggplot2::aes(x = variable, y = rows,  text = value)) + 
    ggplot2::geom_raster(ggplot2::aes(fill = valueType)) + 
    ggplot2::theme_minimal() + 
    # ggplot2::theme(axis.text.x = ggplot2::element_text(size = 14, angle = 45,  vjust = 1, hjust = 1)) +
    # ggplot2::theme(axis.text.x = ggplot2::element_text(size = 22, angle = 55, vjust=0, hjust=0),
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 20, angle = 55, vjust=0, hjust=0),
                   axis.text.y  = ggplot2::element_text(size=22),
                   axis.title.y = ggplot2::element_text(size=24),
                   title        = ggplot2::element_text(size=22, hjust=0.5))+
    ggplot2::labs(x = "", y = "Observation (Row) Number") + 
    ggplot2::scale_y_reverse() +
    # ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 0.5)) + 
    # ggplot2::theme(plot.margin = margin(0,2,0,0, "cm" ))
    # ggplot2::theme(plot.margin = margin(0,20,0,0)) +
    # ggplot2::coord_cartesian(clip="off") +
    ggplot2::guides(colour = "none")
}
# <bytecode: 0x0000020079de1210>
#   <environment: namespace:visdat>

#   > visdat:::miss_guide_label
miss_guide_label <- function (x) 
{
  p_miss <- (mean(is.na(x)) * 100)
  if (p_miss == 0) {
    p_miss_lab <- "No Missing Values"
    p_pres_lab <- "Present (100%)"
  }
  else if (p_miss < 0.1) {
    p_miss_lab <- "Missing (< 0.1%)"
    p_pres_lab <- "Present (> 99.9%)"
  }
  else {
    p_miss <- round(p_miss, 1)
    p_pres <- 100 - p_miss
    p_miss_lab <- glue::glue("Missing \n({p_miss}%)")
    p_pres_lab <- glue::glue("Present \n({p_pres}%)")
  }
  label_frame <- tibble::tibble(p_miss_lab, p_pres_lab)
  return(label_frame)
}
# <bytecode: 0x000002007c5602d8>
#   <environment: namespace:visdat>
#   

# > visdat:::fingerprint
fingerprint <-  function (x) 
{
  if (!is.list(x)) {
    ifelse(is.na(x), yes = NA, no = glue::glue_collapse(class(x), 
                                                        sep = "\n"))
  }
  else {
    ifelse(purrr::map_lgl(x, ~length(.x) == 0), yes = NA, 
           no = glue::glue_collapse(class(x), sep = "\n"))
  }
}
# <bytecode: 0x00000200748f7100>
#   <environment: namespace:visdat>
# ff <- fingerprint(ndGLM1$nest_age)
# ff
# 
# any_miss <- function (x)
# {
#   if(!is.list(x)){  
#     ifelse(any(is.na(x)), )
#   }
# }

#   > visdat:::fingerprint_df
fingerprint_df <- function (x) 
{
  purrr::map_df(x, fingerprint)
}
# <bytecode: 0x000002007c80e260>
#   <environment: namespace:visdat>
#   

# > visdat:::label_col_missing_pct
label_col_missing_pct <-  function (x, col_order_index, labs) 
{
  labelled_pcts <- colMeans(is.na(x))[col_order_index] %>% 
    purrr::map_chr(function(x) {
      dplyr::case_when(x == 0 ~ "0%",
                       x >= 0.001 ~ scales::percent(x,  accuracy = 1),
                       x < 0.001 ~ "<0.1%")
    })
  # glue::glue("{col_order_index} ({labelled_pcts})")
  glue::glue("{labs[match(col_order_index, names(labs))]} ({labelled_pcts})")
}
# <bytecode: 0x0000020077133050>
#   <environment: namespace:visdat>

# > visdat::vis_miss
vis_miss <- function (x, cluster = FALSE, sort_miss = FALSE, show_perc = TRUE, 
          # show_perc_col = TRUE, large_data_size = 9e+05, warn_large_data = TRUE, 
          show_perc_col = TRUE, large_data_size = 9e+05, warn_large_data = FALSE, 
          facet, pTitle, labs) 
{
  # test_if_dataframe(x)
  # test_if_large_data(x, large_data_size, warn_large_data)
  # ------------------------------------------------------------------------------------
  # sort columns?
  # x <- ndGLM1
  # if (sort_miss) {
  #   col_order_index <- names(n_miss_col(x, sort = TRUE))
  # }
  # else if (!sort_miss) {
  #   col_order_index <- names(x)
  # }
  # ------------------------------------------------------------------------------------
  
  # facet?
  # if (!missing(facet)) {
  #   vis_miss_data <- x %>% dplyr::group_by({
  #     {
  #       facet
  #     }
  #   }) %>% data_vis_miss(cluster)
  #   col_order_index <- update_col_order_index(col_order_index, 
  #                                             facet, environment())
  # }
  # else {
  #   vis_miss_data <- visdat::data_vis_miss(x, cluster)
  # }
  # ------------------------------------------------------------------------------------
  
  col_order_index <- names(x)
  vis_miss_data <- visdat::data_vis_miss(x, cluster) # creates tidy dataframe for plotting
  x_fingerprinted <- fingerprint_df(x) # creates a "fingerprint" of the data - type where value is present, NA where missing
  
  if (show_perc) {
    temp <- miss_guide_label(x_fingerprinted)
    p_miss_lab <- temp$p_miss_lab
    p_pres_lab <- temp$p_pres_lab
  }
  else {
    p_miss_lab <- "Missing"
    p_pres_lab <- "Present"
  }
  
  vis_miss_plot <- vis_create_(vis_miss_data) + 
    # ggplot2::scale_x_discrete(labels=labs) +
    ggplot2::scale_fill_manual(name = "",                           # fill color specifications & names (for legend)
                               values = c("grey80", "grey20"), 
                               labels = c(p_pres_lab,  p_miss_lab)) +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) + 
    ggplot2::theme(legend.position = "bottom",
                   legend.text     = element_text(size=24),
                   legend.title    = element_text(size=20)) + 
    # ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 0),
    # ggplot2::theme(axis.text.y  = ggplot2::element_text(size=22),
    #                axis.title.y = ggplot2::element_text(size=24),
    #                title        = ggplot2::element_text(size=22, hjust=0.5))+
    ggplot2::ggtitle(pTitle) 
    # ggplot2::lege
    
  
  # ------------------------------------------------------------------------------------
  if (ncol(x) == 1) {
    if (show_perc_col) {
      return(vis_miss_plot <- vis_miss_plot +
               ggplot2::theme(plot.margin = margin(10,60,10,10)) +
               # ggplot2::theme(plot.margin = margin(0,20,0,0)) +
               ggplot2::coord_cartesian(clip="off") +
               ggplot2::scale_x_discrete(position = "top", 
                                         labels = label_col_missing_pct(x_fingerprinted, 
                                                                        col_order_index,
                                                                        labs=labs )))
    }
    else if (!show_perc_col) {
      return(vis_miss_plot <- vis_miss_plot +
               ggplot2::scale_x_discrete(position = "top", 
                                         labels = col_order_index))
    }
  }
  
  # ------------------------------------------------------------------------------------
  if (!missing(facet)) {
    vis_miss_plot <- vis_miss_plot + ggplot2::facet_wrap(facets = dplyr::vars({
      {
        facet
      }
    }))
  }
  
  # ------------------------------------------------------------------------------------
  if (show_perc_col && missing(facet)) {
    vis_miss_plot <- vis_miss_plot + 
      ggplot2::theme(plot.margin = margin(10,60,10,10)) +
      ggplot2::coord_cartesian(clip="off") +
      ggplot2::scale_x_discrete(position = "top",
                                limits = col_order_index,
                                labels = label_col_missing_pct(x_fingerprinted, 
                                                               col_order_index, 
                                                               labs=labs))
  }
  else {
    vis_miss_plot <- vis_miss_plot + 
      ggplot2::scale_x_discrete(position = "top", 
                                limits = col_order_index)
  }
  
  # ------------------------------------------------------------------------------------
  return(vis_miss_plot)
}
# <bytecode: 0x0000020072144190>
#   <environment: namespace:visdat>
vNames <- c("nest_age" = "NEST_AGE", 
            "obs_int" = "FINAL_INTERVAL", 
            "fdate" = "END_DATE",
            "species" = "SPECIES",
            "is_u" = "MARKED_UNKNOWN",
            "HF_mis" = "MISCLASSIFIED",
            "cam_fate" = "TRUE_FATE",
            "any_missing" = "ANY_MISSING"
            )
ndGLM1 <- readRDS("dataframes/ndGLM1.rds")

# vis_miss(dat4imp)
# ndGLM1$anyM <- ifelse(is.na(ndGLM1), NA, "no") # this is colwise, not rowwise
allPl <- ndGLM1 %>%
  # select(species, cam_fate, fdate, nest_age, obs_int) %>%
  select(is_u, HF_mis, species, cam_fate,  obs_int, fdate, nest_age) %>%
#   rowwise %>%
#   mutate(anyMiss = ifelse(any(is.na(.)), NA, "numeric")) %>%
#   # mutate(anyMiss = ifelse(if_any(is.na(.)), NA, "numeric")) %>%
  vis_miss(pTitle="Missingness -\nBy Variable",labs=vNames)
ggsave("figures/allPl.svg", allPl, device="svg", width=11, height=9, unit="in")

# allPljj
# nd_anymiss <- function(x) ifelse(any(is.na(x)), yes=NA, no="no") 
ndGLM1$any_missing <-ndGLM1 %>%
  # select(species, cam_fate, fdate, nest_age, obs_int) %>%
  select(is_u, HF_mis, species, cam_fate,  obs_int, fdate, nest_age) %>%
  # select(is_u, HF_mis, species, cam_fate, fdate, nest_age, obs_int) %>%
  apply(1, function(x) ifelse(any(is.na(x)), yes=NA, no="no") )
# nd_anymiss
# anyPl <- vis_miss(as.data.frame(any_missing), pTitle="Missingness -\nBy Observation", labs=vNames)
anyPl <- ndGLM1 %>% select(any_missing) %>% as.data.frame() %>% vis_miss(pTitle="Missingness-\nBy Observation", labs=vNames)
# anyPl <- vis_miss(as.data.frame(ndGLM1$any_missing), pTitle="Missingness -\nBy Observation", labs=vNames)

ggsave("figures/anyPl.svg", anyPl,device = "svg",width = 4, height=9, unit="in")


# ggsave("figures/allPl.svg", allPl, device="svg", width=11, height=9, unit="in")
# nd_anymiss <-function(nd){
#   lapply(nd, function(x) ifelse(any(is.na(x)), yes=NA, no="no") )

nd_anymiss(ndGLM1)
ndGLM1 %>% 
  select(species, cam_fate, fdate, nest_age, obs_int) %>%
  rowwise %>%
  nd_anymiss()
ndMiss <- purrr::map_df(ndGLM1, nd_anymiss) # goes across cols to see if there is any NA
