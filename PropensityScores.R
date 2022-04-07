## ----setup, include=FALSE---------------------------
library(tidyverse); 
library(broom);
library(DiagrammeR);
library(sandwich)
library(haven)
library(mgcv)
library(MatchIt)
knitr::opts_chunk$set(echo = T, warning = F, message = F);
knitr::knit_hooks$set(mysize = function(before, options, envir) {
  if (before) 
    return(options$size);
})
options(digits = 3);
figure_scaler = 5/8;#1/2 for ioslides; ~1/3 for word, pdf
text_scaler = 1;#1 for ioslides; 2/3 for word, pdf
fig.x = 16 * figure_scaler;
fig.y = 9 * figure_scaler;
cex_scale = 1.75 * figure_scaler;
theme_set(theme_bw())



## ---- echo = F--------------------------------------
create_graph() %>%
  add_node(label = "Z",node_aes = node_aes(color = "black", fontname = "serif")) %>%
  add_node(label = "X",node_aes = node_aes(color = "black", fontname = "serif")) %>%
  add_node(label = "Y",node_aes = node_aes(color = "black", fontname = "serif")) %>%
  add_edge(from = "X", to = "Y",edge_aes = edge_aes(color = "black")) %>%
  add_edge(from = "X", to = "Z",edge_aes = edge_aes(color = "black")) %>%
  add_edge(from = "Z", to = "Y",edge_aes = edge_aes(color = "black")) %>%
  render_graph()


## ---------------------------------------------------

rhc <- 
  read_csv("rhc.csv", show_col_types = FALSE) %>%
  mutate(income = as_factor(income), 
         cat1 = as_factor(cat1), 
         cat2 = as_factor(cat2)) %>%
  rename(subjid = id)

predictors <- setdiff(colnames(rhc),
                      c("subjid",
                        "swang1", #treatment
                        "dth30" #outcome
                      ))

ps_formula <-  
  as.formula(paste0("swang1 ~",paste0(predictors, collapse = "+")))

ps_model <- 
  glm(ps_formula, family = "binomial", data = rhc) 

rhc_ps <- 
  bind_cols(rhc %>% select(subjid, swang1, dth30), 
            model.matrix(ps_model)[,-1],
            ps_swang1 = predict(ps_model, type = "response", newdata = rhc))

predictors <- colnames(model.matrix(ps_model)[,-1])



## ---------------------------------------------------
tidy(ps_model) %>% 
  arrange(desc(abs(statistic))) %>% 
  print(n = Inf)



## ---- echo = FALSE, fig.width=fig.x, fig.height = fig.y----

ggplot() + 
  geom_histogram(data = filter(rhc_ps, 
                               swang1 == 1), 
                 mapping = aes(x = ps_swang1, 
                               y = after_stat(density),
                               fill = "RHC"),
                 bins = 40) + 
  geom_histogram(data = filter(rhc_ps, 
                               swang1 == 0), 
                 mapping = aes(x = ps_swang1,
                               y = -1 * after_stat(density),
                               fill = "No RHC"),
                 bins = 40) +
  scale_x_continuous(name = "Estimated propensity for RHC", 
                     expand = expansion(add = 0.003)) +
  scale_y_continuous(name = "Unweighted frequency", 
                     breaks = NULL,
                     expand = expansion(mult = 0.035)) + 
  scale_fill_manual(name = "Actual treatment delivered", 
                    breaks = c("RHC", "No RHC"),
                    values = c("#4DAF4A", "#377EB8")) + 
  theme(text = element_text(size = 22, family = "serif"),
        legend.position = c(0.8, 0.2));



## ---------------------------------------------------

# Use 0.025-width propensity bins
rhc_ps <- 
  rhc_ps %>% 
  mutate(ps_quantile = cut_width(ps_swang1, width = 0.025, center = 0.0125)) 

rhc_ps %>% 
  count(swang1, ps_quantile) %>% 
  pivot_wider(id_cols = ps_quantile, 
              names_from = swang1, 
              values_from = n,
              values_fill = 0) %>% 
  print(n = Inf)

# Calculate weights based upon strata size, treatment size, and strata-treatment size
rhc_ps <- 
  left_join(rhc_ps, 
            rhc_ps %>% 
              count(ps_quantile, name = "stratum_size"), 
            by = "ps_quantile") %>%
  left_join(rhc_ps %>%
              count(swang1, name = "treatment_size"), 
            by = "swang1") %>%
  left_join(rhc_ps %>%
              count(swang1, ps_quantile, name = "stratum_treatment_size"),
            by = c("swang1", "ps_quantile")) %>%
  mutate(stratified_weight = stratum_size * treatment_size / stratum_treatment_size / n(), 
         .keep = "unused")

# weight given to each observation in each group
rhc_ps %>%
  group_by(swang1, ps_quantile) %>%
  summarize(num = n(), 
            wgt_per_obs = first(stratified_weight), 
            wgt_total  = num * wgt_per_obs) %>% 
  pivot_wider(id_cols = ps_quantile, 
              names_from = swang1, 
              values_from = c(num, wgt_per_obs, wgt_total))



## ---------------------------------------------------
library(MatchIt)
matched_rhc <- matchit(ps_formula, data = rhc, caliper = 0.025, std.caliper = FALSE)

matched_rhc 

matched_rhc_pairs <- 
  get_matches(matched_rhc, weights = "matched_weight")

head(matched_rhc_pairs)

rhc_ps <- 
  left_join(rhc_ps, 
            matched_rhc_pairs %>%
              select(subjid, matched_weight), 
            by = "subjid") %>%
  mutate(matched_weight = ifelse(is.na(matched_weight), 0, matched_weight))



## ---------------------------------------------------
rhc_ps <- 
  rhc_ps %>%
  mutate(ip_weight = swang1 / ps_swang1 + (1 - swang1) / (1 - ps_swang1))


## ---- echo = FALSE, fig.width=fig.x, fig.height = 1.2*fig.y----

weighted.var <- function(x, w, ...) {
  sum(w) * sum(w * (x - weighted.mean(x, w, ...))^2) / (sum(w)^2  - sum(w^2))
}

equal_weighted <- 
  full_join(
    rhc_ps %>%
      group_by(swang1) %>%
      summarize(across(all_of(predictors), mean)) %>%
      arrange(desc(swang1)) %>%
      select(-swang1) %>%
      summarize(across(where(is.numeric), ~ first(.x) - last(.x))) %>%
      pivot_longer(everything(), values_to = "mean_diff"),
    rhc_ps %>%
      group_by(swang1) %>%
      summarize(across(all_of(predictors), var)) %>%
      arrange(desc(swang1)) %>%
      select(-swang1) %>%
      summarize(across(where(is.numeric), ~ sqrt(mean(.x)))) %>%
      pivot_longer(everything(), values_to = "pooled_sd")) %>%
  mutate(std_diff = 100 * mean_diff / pooled_sd)

ip_weighted <-
  full_join(
    rhc_ps %>%
      group_by(swang1) %>%
      summarize(across(all_of(predictors), ~weighted.mean(.x, ip_weight))) %>%
      arrange(desc(swang1)) %>%
      select(-swang1) %>%
      summarize(across(where(is.numeric), ~ first(.x) - last(.x))) %>%
      pivot_longer(everything(), values_to = "mean_diff_ipw"),
    rhc_ps %>%
      group_by(swang1) %>%
      summarize(across(all_of(predictors), ~weighted.var(.x, ip_weight))) %>%
      arrange(desc(swang1)) %>%
      select(-swang1) %>%
      summarize(across(where(is.numeric), ~ sqrt(mean(.x)))) %>%
      pivot_longer(everything(), values_to = "pooled_sd_ipw")) %>%
  mutate(std_diff_ipw = 100 * mean_diff_ipw / pooled_sd_ipw)


stratified_weighted <-
  full_join(
    rhc_ps %>%
      group_by(swang1) %>%
      summarize(across(all_of(predictors), ~weighted.mean(.x, stratified_weight))) %>%
      arrange(desc(swang1)) %>%
      select(-swang1) %>%
      summarize(across(where(is.numeric), ~ first(.x) - last(.x))) %>%
      pivot_longer(everything(), values_to = "mean_diff_stratified"),
    rhc_ps %>%
      group_by(swang1) %>%
      summarize(across(all_of(predictors), ~weighted.var(.x, stratified_weight))) %>%
      arrange(desc(swang1)) %>%
      select(-swang1) %>%
      summarize(across(where(is.numeric), ~ sqrt(mean(.x)))) %>%
      pivot_longer(everything(), values_to = "pooled_sd_stratified")) %>%
  mutate(std_diff_stratified = 100 * mean_diff_stratified / pooled_sd_stratified)

matched_weighted <-
  full_join(
    rhc_ps %>%
      group_by(swang1) %>%
      summarize(across(all_of(predictors), ~weighted.mean(.x, matched_weight))) %>%
      arrange(desc(swang1)) %>%
      select(-swang1) %>%
      summarize(across(where(is.numeric), ~ first(.x) - last(.x))) %>%
      pivot_longer(everything(), values_to = "mean_diff_matched"),
    rhc_ps %>%
      group_by(swang1) %>%
      summarize(across(all_of(predictors), ~weighted.var(.x, matched_weight))) %>%
      arrange(desc(swang1)) %>%
      select(-swang1) %>%
      summarize(across(where(is.numeric), ~ sqrt(mean(.x)))) %>%
      pivot_longer(everything(), values_to = "pooled_sd_matched")) %>%
  mutate(std_diff_matched = 100 * mean_diff_matched / pooled_sd_matched)

ggplot(data = full_join(equal_weighted, 
                        ip_weighted, 
                        by = "name") %>%
         full_join(matched_weighted, 
                   by = "name") %>%
         full_join(stratified_weighted, 
                   by = "name") %>%
         arrange(abs(std_diff)) %>%
         mutate(name = factor(name) %>% fct_inorder())) + 
  geom_vline(xintercept = 0) + 
  geom_point(aes(x = abs(std_diff), y = name, color = "Equal", shape = "Equal"), size = 3) + 
  geom_point(aes(x = abs(std_diff_ipw), y = name, color = "IP", shape = "IP"), size = 3) + 
  geom_point(aes(x = abs(std_diff_stratified), y = name, color = "Stratified", shape = "Stratified"), size = 3) + 
  geom_point(aes(x = abs(std_diff_matched), y = name, color = "Matched", shape = "Matched"), size = 3) + 
  scale_y_discrete(name = NULL) + 
  scale_x_continuous(name = "Absolute Weighted Standardized Difference", 
                     breaks = seq(-30, 50, by = 10), 
                     expand = expansion(mult = 0.02)) + 
  scale_color_brewer(name = NULL, palette ="Dark2") + 
  scale_shape_discrete(name = NULL) + 
  theme(text = element_text(size = 22, family = "serif"),
        legend.position = c(0.8, 0.2));



## ---------------------------------------------------
model_obs <- gam(dth30 ~ swang1, family = "binomial", data = rhc_ps)

obs_diff <- 
  predict(model_obs, type = "response", newdata = rhc_ps %>% slice(1) %>% mutate(swang1 = 1)) - 
  predict(model_obs, type = "response", newdata = rhc_ps %>% slice(1) %>% mutate(swang1 = 0)) 

obs_diff


## ---------------------------------------------------

# Parametric component for treatment; non-parametric component for ps
model_adjust_ps <- gam(dth30 ~ swang1 + s(ps_swang1), family = "binomial", data = rhc_ps)

# Estimated treatment effect is average predicted response if all were treated
resp_rhc <- predict(model_adjust_ps, type = "response", newdata = rhc_ps %>% mutate(swang1 = 1))
resp_not_rhc <- predict(model_adjust_ps, type = "response", newdata = rhc_ps %>% mutate(swang1 = 0))
adjusted_diff <-
  mean(resp_rhc) - 
  mean(resp_not_rhc)

adjusted_diff



## ---- echo = FALSE, fig.width=fig.x, fig.height = 0.9*fig.y----
ggplot() + 
  geom_point(aes(x = rhc_ps %>% pull(ps_swang1), y = resp_rhc - resp_not_rhc),
             size = 1, 
             alpha = 0.5) + 
  scale_x_continuous(name = "Estimated Propensity") + 
  scale_y_continuous(name = "Estimated Treatment Effect") + 
  theme(text = element_text(size = 22, family = "serif"),
        legend.position = c(0.8, 0.2));


## ---------------------------------------------------


# Fit GLM weighting by stratum weights
model_stratified <- glm(dth30 ~ swang1, family = "binomial", weights = stratified_weight, data = rhc_ps)

# Estimated treatment effect is average predicted response if all were treated
# versus if all were untreated
stratified_diff <- 
  predict(model_stratified, type = "response", newdata = rhc_ps %>% slice(1) %>% mutate(swang1 = 1)) - 
  predict(model_stratified, type = "response", newdata = rhc_ps %>% slice(1) %>% mutate(swang1 = 0)) 
stratified_diff

# Same as taking difference in weighted means
rhc_ps %>% 
  group_by(swang1) %>%
  summarize(mean_dth30 = weighted.mean(dth30, stratified_weight)) %>%
  arrange(swang1) %>%
  pull(mean_dth30) %>%
  diff()


## ---------------------------------------------------
matched_diff <- 
  rhc_ps %>% 
  group_by(swang1) %>%
  summarize(mean_dth30 = weighted.mean(dth30, matched_weight)) %>%
  arrange(swang1) %>%
  pull(mean_dth30) %>%
  diff()

matched_diff


## ---------------------------------------------------
model_ipw <- glm(dth30 ~ swang1, family = "binomial", weights = ip_weight, data = rhc_ps)

# Estimated treatment effect is average predicted response if all were treated
# versus if all were untreated
ipw_diff <- 
  predict(model_ipw, type = "response", newdata = rhc_ps %>% slice(1) %>% mutate(swang1 = 1)) - 
  predict(model_ipw, type = "response", newdata = rhc_ps %>% slice(1) %>% mutate(swang1 = 0)) 

ipw_diff


## ---- echo = FALSE----------------------------------
tribble(~"Method", ~"Result",
        "Observed", obs_diff,
        "Adjusted", adjusted_diff,
        "Stratified", stratified_diff,
        "Matched", matched_diff, 
        "IPW", ipw_diff) %>% knitr::kable()


## ---------------------------------------------------
# regular variance estimator
sqrt(diag(vcov(model_ipw)))
# vcovHC from sandwich package for robust variance estimator
sqrt(diag(vcovHC(model_ipw, type = "HC")))


