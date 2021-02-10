# NPDR --------------
prevalence_df = tab2_df %>% 
  mutate_at(vars(d_retinopathy,
                 d_nephropathy,
                 d_retinopathy_yn,
                 d_nephropathy_yn,
                 dkd
  ),function(x) case_when(x == "Missing" ~ NA_character_,
                          TRUE ~ as.character(x))) %>% 
  mutate_at(vars(d_retinopathy_yn,
                 d_nephropathy_yn,
                 dkd),
            function(x) case_when( x  == "0, No" ~ 0,
                                   x == "1, Yes" ~ 1,
                                   TRUE ~ NA_real_))

# NPDR,  PDR -------------------

multinom_retinopathy = prevalence_df %>% 
  nnet::multinom(d_retinopathy ~ cluster + ageonset + d_duration + gender + hba,
                 data=.)
glm_retinopathy = prevalence_df %>% 
  glm(d_retinopathy_yn ~ cluster + ageonset + d_duration + gender+ hba,
      family=binomial(),data=.)

# Micro, Macro ---------------------
multinom_nephropathy = prevalence_df %>% 
  nnet::multinom(d_nephropathy ~ cluster + ageonset + d_duration + gender+ hba,
                 data=.)
glm_nephropathy = prevalence_df %>% 
  glm(d_nephropathy_yn ~ cluster + ageonset + d_duration + gender+ hba,
      family=binomial(),data=.)
glm_dkd = prevalence_df %>% 
  glm(dkd ~ cluster + ageonset + d_duration + gender+ hba,
      family=binomial(),data=.)

# SUMMARY ---------------

summary_df = bind_rows(broom::tidy(multinom_nephropathy) %>% 
                         mutate(type = "multinom nephropathy"),
                       broom::tidy(multinom_retinopathy)%>% 
                         mutate(type = "multinom retinopathy"),
                       broom::tidy(glm_retinopathy) %>% 
                         mutate(type = "glm retinopathy"),
                       broom::tidy(glm_nephropathy) %>% 
                         mutate(type = "glm nephropathy"),
                       broom::tidy(glm_dkd) %>% 
                         mutate(type = "glm dkd")
) %>% 
  mutate(estimate = case_when(type %in% c("multinom nephropathy","multinom retinopathy") ~ log(estimate),
                              TRUE ~ estimate)) %>% 
  mutate(output = paste0("OR = ",exp(estimate) %>% round(.,2)," (",
                         exp(estimate - 1.96*std.error) %>% round(.,2),", ",
                         exp(estimate + 1.96*std.error) %>% round(.,2),")"
  )) %>% 
  dplyr::filter(str_detect(term,"cluster")) 