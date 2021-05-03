# Function to fit the various models of interest to a single dataset
# df = the name of the dataset
# B = the number of bootstrap samples to use in esimating the variance for the interval censoring methods. Set to 3 by default here since interest in these simulations is only in the parameter estimate and not the variance so we use a low number for speed. You would want a much higher number if the variance/confidence intervals/p-values were of itnerest
fitsimmod <-
  function(df, B = 3) {
    
    tibble(
      # Cox model for RFS - censored at lastobs1
      tryCatch({
        tab <- coxph(
          Surv(rfs_months, rfs_event) ~ z,
          data = df
        ) %>% 
          broom::tidy(conf.int = TRUE) %>% 
          select(estimate, conf.low, conf.high) %>% 
          mutate(model = "Cox RFS - lastobs1")
        tab$issue <- NA
        tab
      }, warning = function(w){
        tab <- coxph(
          Surv(rfs_months, rfs_event) ~ z,
          data = df
        ) %>% 
          broom::tidy(conf.int = TRUE) %>% 
          select(estimate, conf.low, conf.high) %>% 
          mutate(model = "Cox RFS - lastobs1")
        tab$issue <- "w"
        tab
      }, error = function(e){
        tab <- coxph(
          Surv(rfs_months, rfs_event) ~ z,
          data = df
        ) %>% 
          broom::tidy(conf.int = TRUE) %>% 
          select(estimate, conf.low, conf.high) %>% 
          mutate(model = "Cox RFS - lastobs1")
        tab$issue <- "e"
        tab
      })
    ) %>% 
      add_row(
        # Cox model for RFS - censored at dtime
        tryCatch({
          tab <- coxph(
            Surv(rfs_months_2, rfs_event) ~ z,
            data = df
          ) %>% 
            broom::tidy(conf.int = TRUE) %>% 
            select(estimate, conf.low, conf.high) %>% 
            mutate(model = "Cox RFS - dtime")
          tab$issue <- NA
          tab
        }, warning = function(w){
          tab <- coxph(
            Surv(rfs_months_2, rfs_event) ~ z,
            data = df
          ) %>% 
            broom::tidy(conf.int = TRUE) %>% 
            select(estimate, conf.low, conf.high) %>% 
            mutate(model = "Cox RFS - dtime")
          tab$issue <- "w"
          tab
        }, error = function(e){
          tab <- coxph(
            Surv(rfs_months_2, rfs_event) ~ z,
            data = df
          ) %>% 
            broom::tidy(conf.int = TRUE) %>% 
            select(estimate, conf.low, conf.high) %>% 
            mutate(model = "Cox RFS - dtime")
          tab$issue <- "e"
          tab
        })
      ) %>% 
      add_row(
        # Cox interval censoring model for RFS - censored at lastobs1
        tryCatch({
          tab <- ic_sp(
            Surv(left_time, right_time, type = "interval2") ~ z,
            model = "ph",
            bs_samples = B,
            data = df
          ) %>% 
            tidy_ic_sp() %>% 
            select(estimate, conf.low, conf.high) %>% 
            mutate(model = "Interval RFS - lastobs1")
          tab$issue <- NA
          tab
        }, warning = function(w){
          tab <- ic_sp(
            Surv(left_time, right_time, type = "interval2") ~ z,
            model = "ph",
            bs_samples = B,
            data = df
          ) %>% 
            tidy_ic_sp() %>% 
            select(estimate, conf.low, conf.high) %>% 
            mutate(model = "Interval RFS - lastobs1")
          tab$issue <- "w"
          tab
        }, error = function(e){
          tab <- ic_sp(
            Surv(left_time, right_time, type = "interval2") ~ z,
            model = "ph",
            bs_samples = B,
            data = df
          ) %>% 
            tidy_ic_sp() %>% 
            select(estimate, conf.low, conf.high) %>% 
            mutate(model = "Interval RFS - lastobs1")
          tab$issue <- "e"
          tab
        })
      ) %>% 
      add_row(
        # Cox interval censoring model for RFS - censored at dtime
        tryCatch({
          tab <- ic_sp(
            Surv(left_time_2, right_time, type = "interval2") ~ z,
            model = "ph",
            bs_samples = B,
            data = df
          ) %>% 
            tidy_ic_sp() %>% 
            select(estimate, conf.low, conf.high) %>% 
            mutate(model = "Interval RFS - dtime")
          tab$issue <- NA
          tab
        }, warning = function(w){
          tab <- ic_sp(
            Surv(left_time_2, right_time, type = "interval2") ~ z,
            model = "ph",
            bs_samples = B,
            data = df
          ) %>% 
            tidy_ic_sp() %>% 
            select(estimate, conf.low, conf.high) %>% 
            mutate(model = "Interval RFS - dtime")
          tab$issue <- "w"
          tab
        }, error = function(e){
          tab <- ic_sp(
            Surv(left_time_2, right_time, type = "interval2") ~ z,
            model = "ph",
            bs_samples = B,
            data = df
          ) %>% 
            tidy_ic_sp() %>% 
            select(estimate, conf.low, conf.high) %>% 
            mutate(model = "Interval RFS - dtime")
          tab$issue <- "e"
          tab
        })
      ) %>%
      add_row(
        # Traditional cause-specific for recurrence
        tryCatch({
          tab <- coxph(
            Surv(rcr_months, rcr_event) ~ z,
            data = df
          )%>% 
            broom::tidy(conf.int = TRUE) %>% 
            select(estimate, conf.low, conf.high) %>% 
            mutate(model = "CS Cox recurrence")
          tab$issue <- NA
          tab
        }, warning = function(w){
          tab <- coxph(
            Surv(rcr_months, rcr_event) ~ z,
            data = df
          )%>% 
            broom::tidy(conf.int = TRUE) %>% 
            select(estimate, conf.low, conf.high) %>% 
            mutate(model = "CS Cox recurrence")
          tab$issue <- "w"
          tab
        }, error = function(e){
          tab <- coxph(
            Surv(rcr_months, rcr_event) ~ z,
            data = df
          )%>% 
            broom::tidy(conf.int = TRUE) %>% 
            select(estimate, conf.low, conf.high) %>% 
            mutate(model = "CS Cox recurrence")
          tab$issue <- "e"
          tab
        })
      ) %>% 
      add_row(
        # Traditional cause-specific for death
        tryCatch({
          tab <- coxph(
            Surv(death_time_1, death_event_1) ~ z,
            data = df
          )%>% 
            broom::tidy(conf.int = TRUE) %>% 
            select(estimate, conf.low, conf.high) %>% 
            mutate(model = "CS Cox death")
          tab$issue <- NA
          tab
        }, warning = function(w){
          tab <- coxph(
            Surv(death_time_1, death_event_1) ~ z,
            data = df
          )%>% 
            broom::tidy(conf.int = TRUE) %>% 
            select(estimate, conf.low, conf.high) %>% 
            mutate(model = "CS Cox death")
          tab$issue <- "w"
          tab
        }, error = function(e){
          tab <- coxph(
            Surv(death_time_1, death_event_1) ~ z,
            data = df
          )%>% 
            broom::tidy(conf.int = TRUE) %>% 
            select(estimate, conf.low, conf.high) %>% 
            mutate(model = "CS Cox death")
          tab$issue <- "e"
          tab
        })
      ) %>% 
    add_row(
      # Cox interval censoring model for recurrence
      tryCatch({
        tab <- ic_sp(
          Surv(left_time2, right_time2, type = "interval2") ~ z,
          model = "ph",
          bs_samples = B,
          data = df
        ) %>% 
          tidy_ic_sp() %>% 
          select(estimate, conf.low, conf.high) %>% 
          mutate(model = "CS Interval recurrence")
        tab$issue <- NA
        tab
      }, warning = function(w){
        tab <- ic_sp(
          Surv(left_time2, right_time2, type = "interval2") ~ z,
          model = "ph",
          bs_samples = B,
          data = df
        ) %>% 
          tidy_ic_sp() %>% 
          select(estimate, conf.low, conf.high) %>% 
          mutate(model = "CS Interval recurrence")
        tab$issue <- "w"
        tab
      }, error = function(e){
        tab <- ic_sp(
          Surv(left_time2, right_time2, type = "interval2") ~ z,
          model = "ph",
          bs_samples = B,
          data = df
        ) %>% 
          tidy_ic_sp() %>% 
          select(estimate, conf.low, conf.high) %>% 
          mutate(model = "CS Interval recurrence")
        tab$issue <- "e"
        tab
      })
    ) 
  }