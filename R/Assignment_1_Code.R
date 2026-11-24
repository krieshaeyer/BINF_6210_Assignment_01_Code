rm(list = ls())

### YS - The code "rm(list = ls())" recommended by Dr. Cottenie is wonderful. Clears your variables before you rerun, and I find it really helps you find errors.

install.packages("tidyverse")

### YS - completely personal preference but I personally would have liked to see the install packages at the top and commented out. Just so that when you have codes that use more niche packages, people are able to fully execute your code from your code alone without needing to search! 

# 01 Load libraries + datasets ----

library(tidyverse)
library(viridis)
library(ggtext)
library(ggplot2)

t_spiralis <- read_csv("../data/GCA_001447595.1.csv")
t_britovi <- read_csv("../data/GCA_001447585.1.csv")

#excluding non-TE matching classes, and adding row numbers to create unique identifier for each row, enabling exclusion of outliers

ls_non_te <- c("Simple_repeat", "Low_complexity", "Satellite", "snRNA", "rRNA", "ARTEFACT", "Unknown", "simple_repeat")

'%notin%' <- Negate('%in%')

### YS - %notin% is throwing an error. Please either define it at the top of the script or use standard negation (!) instead. I've done negation as defining can be a one time at the beginning. 

t_spiralis <- t_spiralis |> mutate(row = 1:n()) |> filter(matching_class %notin% ls_non_te)
dim(t_spiralis) #1037

t_britovi <- t_britovi |> mutate(row = 1:n()) |> filter(matching_class %notin% ls_non_te)
dim(t_britovi) #1177


# 02 Identify, Plot, and Remove Outliers (using interquartile range method) ----
Q1_spiralis <- quantile(t_spiralis$perc_div, .25)
Q3_spiralis <- quantile(t_spiralis$perc_div, .75)
IQR_spiralis <- IQR(t_spiralis$perc_div)
outliers_spiralis <- subset(t_spiralis, t_spiralis$perc_div<(Q1_spiralis - 1.5*IQR_spiralis) | t_spiralis$perc_div>(Q3_spiralis + 1.5*IQR_spiralis))
dim(outliers_spiralis) #7

t_spiralis |> ggplot(aes(matching_class,y = perc_div)) + geom_boxplot() + geom_hline(aes(yintercept = Q1_spiralis - 1.5*IQR_spiralis, color = "Q1")) + geom_hline(aes(yintercept = Q3_spiralis + 1.5*IQR_spiralis,  colour = "Q3")) + scale_color_manual(name = "Outlier Thresholds", values = c("red", "blue")) + guides(x = guide_axis(angle = 45)) + labs(title = "Outlier Detection by TE Matching Class - *T. spiralis*", x = "TE Matching Class", y = "Percent Divergence") + theme(plot.title = element_markdown())

Q1_britovi <- quantile(t_britovi$perc_div, .25)
Q3_britovi <- quantile(t_britovi$perc_div, .75)
IQR_britovi <- IQR(t_britovi$perc_div)
outliers_britovi <- subset(t_britovi, t_britovi$perc_div<(Q1_britovi - 1.5*IQR_britovi) | t_britovi$perc_div>(Q3_britovi + 1.5*IQR_britovi)) 
dim(outliers_britovi) #1

t_britovi |> ggplot(aes(matching_class,y = perc_div)) + geom_boxplot() + geom_hline(aes(yintercept = Q1_britovi - 1.5*IQR_britovi, color = "Q1")) + geom_hline(aes(yintercept = Q3_britovi + 1.5*IQR_britovi,  colour = "Q3" )) + scale_color_manual(name = "Outlier Threshold", values = c("red", "blue"))+ guides(x = guide_axis(angle = 45)) + labs(title = "Outlier Detection by TE Matching Class - *T. britovi*", x = "TE Matching Class", y = "Percent Divergence") + theme(plot.title = element_markdown())

#Remove outliers from dataset
t_spiralis_out_removed <- t_spiralis |> filter(row %notin% outliers_spiralis$row)

dim(t_spiralis_out_removed) #1030

t_britovi_out_removed <- t_britovi |> filter(row %notin% outliers_britovi$row)
dim(t_britovi_out_removed) #1176


# 03 Perform Wilcoxon tests: percent divergence between species, within each TA mathcing_class ----

### YS - something I think worth considering are satisfying the assumptions/caveats of Wilcoxin tests = dependency on sample size. We don't want that to skew/misrepresent your data so perhaps a good idea is to include a short summary where you can see sample sizes/stats of the two compared groups?

te_combined <- bind_rows(
  t_spiralis |> mutate(species = "spiralis"),
  t_britovi  |> mutate(species = "britovi")
) |>
  filter(
    !matching_class %in% ls_non_te,
    !is.na(matching_class),
    is.finite(perc_div)
  )

te_check <- te_combined |>
  group_by(species, matching_class) |>
  summarise(
    count = n(),
    med_div = median(perc_div, na.rm = TRUE),
    mis_div = sum(is.na(perc_div)),
    .groups = "drop"
  ) |>
  arrange(matching_class, species)
print(te_check)

te_combined_withoutoutliers <- bind_rows(
  t_spiralis_out_removed |> mutate(species = "spiralis"),
  t_britovi_out_removed  |> mutate(species = "britovi")
) |>
  filter(
    !matching_class %in% ls_non_te,
    !is.na(matching_class),
    is.finite(perc_div)
  )

te_check_no_outliers <- te_combined_withoutoutliers |>
  group_by(species, matching_class) |>
  summarise(
    count = n(),
    med_div = median(perc_div, na.rm = TRUE),
    mis_div = sum(is.na(perc_div)),
    .groups = "drop"
  ) |>
  arrange(matching_class, species)
print(te_check_no_outliers)

### YS - consider the large differences in sample size between the species in some classes (such as Harbinger and Mariner - could potentially need to look into consequences on the statistics and whether wilcoxin is best in this case! I think this is very informative for further statistical testing!)

### YS - hmmm I'm wondering why the outlier-corrected sets weren't used for further analysis as outliers were detected? Personally, I would have found a lot of value in generating a comparison plot/image to compare the dataset with and without the outlier! I've included it below just in case it is something you'd be curious of too! (te_combined_withoutoutliers)

wilcox_results <- te_combined |>
  group_by(matching_class) |>
  group_modify(~ {
    x <- .x$perc_div[.x$species == "spiralis"]
    y <- .x$perc_div[.x$species == "britovi"]
    
    # Skip classes with fewer than two observations in either species
    # or where every value is identical
    can_test <- length(x) >= 2 &&
      length(y) >= 2 &&
      n_distinct(c(x, y)) > 1
    
    tibble(
      n_spiralis = length(x),
      n_britovi = length(y),
      median_spiralis = if (length(x) > 0) median(x) else NA_real_,
      median_britovi = if (length(y) > 0) median(y) else NA_real_,
      p_value = if (can_test) {
        wilcox.test(
          x, y,
          alternative = "two.sided",
          paired = FALSE,
          exact = FALSE
        )$p.value
      } else NA_real_
    )
  }) |>
  ungroup() |>
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significant = p_adjusted < 0.05
  ) |>
  arrange(p_adjusted)

print(wilcox_results, n = Inf)

wilcox_results_wo <- te_combined_withoutoutliers |>
  group_by(matching_class) |>
  group_modify(~ {
    x <- .x$perc_div[.x$species == "spiralis"]
    y <- .x$perc_div[.x$species == "britovi"]
    
    # Skip classes with fewer than two observations in either species
    # or where every value is identical
    can_test <- length(x) >= 2 &&
      length(y) >= 2 &&
      n_distinct(c(x, y)) > 1
    
    tibble(
      n_spiralis = length(x),
      n_britovi = length(y),
      median_spiralis = if (length(x) > 0) median(x) else NA_real_,
      median_britovi = if (length(y) > 0) median(y) else NA_real_,
      p_value = if (can_test) {
        wilcox.test(
          x, y,
          alternative = "two.sided",
          paired = FALSE,
          exact = FALSE
        )$p.value
      } else NA_real_
    )
  }) |>
  ungroup() |>
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significant = p_adjusted < 0.05
  ) |>
  arrange(p_adjusted)

print(wilcox_results_wo, n = Inf)

# pretty cool! It looks like there were some changes, especially for those lined past ranking 10! We may be seeing a large impact on smaller groups by outliers so I personally would continue to visualize both together! 

# 04 Plot average percent divergence within matching_class between t_spiralis & t_britovi using boxplots ----


# find distribution of %div of classes in britovi
shared_classes <- intersect(t_spiralis$matching_class, t_britovi$matching_class)

te_shared <- bind_rows(
  t_britovi  |> 
    select(matching_class, perc_div) |> 
    mutate(species = "T. britovi"),
  t_spiralis |> 
    select(matching_class, perc_div) |> 
    mutate(species = "T. spiralis")) |>
  filter(matching_class %in% shared_classes)

te_shared |> filter(species == "T. britovi") |>
  ggplot(aes(x = perc_div, y = matching_class)) + geom_boxplot()

te_shared |> filter(species == "T. spiralis") |>
  ggplot(aes(x = perc_div, y = matching_class)) + geom_boxplot()

# combine datasets to make combined graph
te_div <- ggplot(data = te_shared, mapping = aes(x = perc_div, y = matching_class, fill = species)) +
  geom_boxplot() +
  labs(
    title = "% Divergence in Shared TE Classes of T. britovi and T. spiralis",
    x = "% Divergence",
    y = "TE Class")

### YS - I personally like to code into but comment out how you export the figures! It really helps with organizing file names, figure numbers and sizing! Depending on how you size your image, you can make some things easier and harder to see -- for example a figure that is super busy horizontally benefits from a longer width than height - it'll help with resolving busy datapoints! I like to make sure the reader will see exactly what I see before I make conclusions, but it would be up to your personal preference!

#making an outlier-free version to compare/replace (I think it is beneficial as you already had set up the dataset without outlier already, and it is useful for supporting your conclusions or further steps): 

shared_classes_no_outliers <- intersect(t_spiralis_out_removed$matching_class, t_britovi_out_removed$matching_class)

te_shared_no_outliers <- bind_rows(
  t_britovi_out_removed  |> 
    select(matching_class, perc_div) |> 
    mutate(species = "T. britovi"),
  t_spiralis_out_removed |> 
    select(matching_class, perc_div) |> 
    mutate(species = "T. spiralis")) |>
  filter(matching_class %in% shared_classes_no_outliers)

te_shared_no_outliers |> filter(species == "T. britovi") |>
  ggplot(aes(x = perc_div, y = matching_class)) + geom_boxplot()

te_shared_no_outliers |> filter(species == "T. spiralis") |>
  ggplot(aes(x = perc_div, y = matching_class)) + geom_boxplot()

# combine datasets to make combined graph
te_div_no_out <- ggplot(data = te_shared_no_outliers, mapping = aes(x = perc_div, y = matching_class, fill = species)) +
  geom_boxplot() +
  labs(
    title = "% Divergence in Shared TE Classes of T. britovi and T. spiralis (Outlier-free)",
    x = "% Divergence",
    y = "TE Class")

### YS - combining a side-by-side panel to show the significance of the outliers

library(gridExtra)
comparisondiv_outnoout <- grid.arrange(te_div, te_div_no_out, ncol=2)

ggsave(
  filename = "DivPlot_OutNoOut.png",
  plot = comparisondiv_outnoout,
  width = 16,
  height = 8,
  units = "in",
  dpi = 300
)

### YS - The results seem super interesting! The specific impact of outliers on Mariner could give us an interesting avenue of "next" steps as to the differences in divergence (of your main question).

