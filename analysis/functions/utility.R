# rounds up, then centers on (integer) midpoint of the rounding points
roundmid_any <- function(x, to=1){
  ceiling(x/to)*to - (floor(to/2)*(x!=0))
}


# function for printing summaries of datasets to txt files
my_skim <- function(
    .data, # dataset to be summarised
    path,
    id_suffix = "_id" # (set to NULL if no id columns)
) {
  
  # don't run when in an interactive session as for some reason is causes R to crash
  if(Sys.getenv("OPENSAFELY_BACKEND") != "") {
    
    # specify summary function for each class
    my_skimmers <- list(
      logical = skimr::sfl(
      ),
      # numeric applied to numeric and integer
      numeric = skimr::sfl(
        mean = ~ mean(.x, na.rm=TRUE),
        sd = ~ sd(.x, na.rm=TRUE),
        min = ~ min(.x, na.rm=TRUE),
        p10 = ~ quantile(.x, p=0.1, na.rm=TRUE, type=1),
        p25 = ~ quantile(.x, p=0.25, na.rm=TRUE, type=1),
        p50 = ~ quantile(.x, p=0.5, na.rm=TRUE, type=1),
        p75 = ~ quantile(.x, p=0.75, na.rm=TRUE, type=1),
        p90 = ~ quantile(.x, p=0.9, na.rm=TRUE, type=1),
        max = ~ max(.x, na.rm=TRUE)
      ),
      character = skimr::sfl(),
      factor = skimr::sfl(),
      Date = skimr::sfl(
        # wrap in as.Date to avoid errors when all missing
        min = ~ as.Date(min(.x, na.rm=TRUE)),
        p50 = ~ as.Date(quantile(.x, p=0.5, na.rm=TRUE, type=1)),
        max = ~ as.Date(max(.x, na.rm=TRUE))
      ),
      POSIXct = skimr::sfl(
        # wrap in as.POSIXct to avoid errors when all missing
        min = ~ as.POSIXct(min(.x, na.rm=TRUE)),
        p50 = ~ as.POSIXct(quantile(.x, p=0.5, na.rm=TRUE, type=1)),
        max = ~ as.POSIXct(max(.x, na.rm=TRUE))
      )
    )
    
    my_skim_fun <- skimr::skim_with(
      !!!my_skimmers,
      append = FALSE
    )
    
    # summarise factors as the printing is not very nice or flexible in skim
    summarise_factor <- function(var) {
      
      out <- .data %>%
        group_by(across(all_of(var))) %>%
        count() %>%
        ungroup() %>%
        mutate(across(n, ~roundmid_any(.x, to = 7))) %>%
        mutate(percent = round(100*n/sum(n),2)) %>%
        arrange(!! sym(var)) 
      
      total <- nrow(out)
      
      out %>%
        slice(1:min(total, 10)) %>% 
        knitr::kable(
          format = "pipe",
          caption = glue::glue("{min(total, 10)} of {total} factor levels printed")
        ) %>% 
        print()
      
    }
    
    vars <- .data %>% 
      select(-ends_with(id_suffix)) %>% 
      select(where(~ is.factor(.x) | is.character(.x))) %>%
      names()
    
    options(width = 120)
    capture.output(
      {
        cat("The following id variables are removed from this summary:\n")
        print(.data %>% select(ends_with(id_suffix)) %>% names())
        cat("\n")
        print(my_skim_fun(.data, -ends_with(id_suffix)))
        cat("\n")
        cat("--- counts for factor and character variables ---")
        for (v in vars) {
          summarise_factor(v)
        }
      },
      file = path,
      append = FALSE
    )
    
  } else {
    
    capture.output(
      cat("Don't run skim in interactice sessions."),
      file = path,
      append = FALSE
    )
    
  }
}

flow_stats_rounded <- function(.data, to) {
  .data %>%
    mutate(
      n = roundmid_any(n, to = to),
      n_exclude = lag(n) - n,
      pct_exclude = n_exclude/lag(n),
      pct_all = n / first(n),
      pct_step = n / lag(n),
    ) %>%
    mutate(across(starts_with("pct_"), ~round(.x, 3)))
}
