## Packages --------------------------------------------------------------------
# Data analysis
library(data.table)
library(tidyverse)

# MSstats packages
library(MSstats)
library(MSstatsConvert)
library(MSstatsBig)

annotation = fread("MetaData_2_5_2025.csv")

# Process V Instrument ---------------------------------------------------------
v_instrument_input = "spectronaut_input_v.csv"

v_converted_data <- bigSpectronauttoMSstatsFormat(
  v_instrument_input,
  "v_output_file.csv",
  intensity = "F.NormalizedPeakArea",
  filter_by_excluded = TRUE,
  filter_by_identified = FALSE,
  filter_by_qvalue = TRUE,
  qvalue_cutoff = 0.01,
  filter_unique_peptides = TRUE,
  aggregate_psms = TRUE,
  filter_few_obs = FALSE,
  remove_annotation = FALSE,
  calculateAnomalyScores=TRUE,
  anomalyModelFeatures=c("FG.ShapeQualityScore..MS1.",
                         "FG.ShapeQualityScore..MS2.",
                         "EG.ApexRT"),
  backend="arrow")

v_converted_data = collect(v_converted_data)
v_converted_data = as.data.table(v_converted_data)

v_converted_data = MSstatsConvert::MSstatsBalancedDesign(
  v_converted_data, 
  c("PeptideSequence", "PrecursorCharge", "FragmentIon", "ProductCharge"), 
  fill_incomplete = TRUE, handle_fractions = TRUE, fix_missing = FALSE,
  remove_few = FALSE, 
  anomaly_metrics = c("FGShapeQualityScoreMS1", "FGShapeQualityScoreMS2", "EGApexRT"))

v_converted_data$Intensity = ifelse(
  v_converted_data$Intensity == 0,
  NA,
  v_converted_data$Intensity
)

anomalyModelFeatures=c("FGShapeQualityScoreMS1", 
                       "FGShapeQualityScoreMS2", 
                       "EGApexRT")
anomalyModelFeatureTemporal=c("mean_decrease",
                              "mean_decrease",
                              "dispersion_increase")
runOrder=temporal
n_trees=100
max_depth="auto"
numberOfCores=8

anomalyModelFeatures=MSstatsConvert:::.standardizeColnames(anomalyModelFeatures)

msstats_data_v = MSstatsConvert::MSstatsAnomalyScores(
  v_converted_data, anomalyModelFeatures,
  anomalyModelFeatureTemporal, .5, 100, runOrder, n_trees,
  max_depth, numberOfCores)

summarized_v = dataProcess(msstats_data_v,
                           normalization="equalizeMedians",
                           featureSubset = "topN",
                           n_top_feature = 100,
                           MBimpute=TRUE,
                           summaryMethod="linear",
                           maxQuantileforCensored = 0.99,
                           numberOfCores = 12)

summarized_v$AnomalyScores=NULL
summarized_v_msstats = dataProcess(summarized_v,
                                   normalization="equalizeMedians",
                                   featureSubset = "topN",
                                   n_top_feature = 100,
                                   MBimpute=TRUE,
                                   summaryMethod="TMP",
                                   numberOfCores = 12)

# Process T Instrument ---------------------------------------------------------
t_instrument_input = "spectronaut_input_t.csv"

t_converted_data <- bigSpectronauttoMSstatsFormat(
  t_instrument_input,
  "t_output_file.csv",
  intensity = "F.NormalizedPeakArea",
  filter_by_excluded = TRUE,
  filter_by_identified = FALSE,
  filter_by_qvalue = TRUE,
  qvalue_cutoff = 0.01,
  filter_unique_peptides = TRUE,
  aggregate_psms = TRUE,
  filter_few_obs = FALSE,
  remove_annotation = FALSE,
  calculateAnomalyScores=TRUE,
  anomalyModelFeatures=c("FG.ShapeQualityScore..MS1.",
                         "FG.ShapeQualityScore..MS2.",
                         "EG.ApexRT"),
  backend="arrow")

t_converted_data = collect(t_converted_data)
t_converted_data = as.data.table(t_converted_data)

t_converted_data = MSstatsConvert::MSstatsBalancedDesign(
  t_converted_data, 
  c("PeptideSequence", "PrecursorCharge", "FragmentIon", "ProductCharge"), 
  fill_incomplete = TRUE, handle_fractions = TRUE, fix_missing = FALSE,
  remove_few = FALSE, 
  anomaly_metrics = c("FGShapeQualityScoreMS1", "FGShapeQualityScoreMS2", "EGApexRT"))

t_converted_data$Intensity = ifelse(
  t_converted_data$Intensity == 0, NA, t_converted_data$Intensity
)

anomalyModelFeatures=c("FGShapeQualityScoreMS1", 
                       "FGShapeQualityScoreMS2", 
                       "EGApexRT")
anomalyModelFeatureTemporal=c("mean_decrease",
                              "mean_decrease",
                              "dispersion_increase")
runOrder=temporal
n_trees=100
max_depth="auto"
numberOfCores=8

anomalyModelFeatures=MSstatsConvert:::.standardizeColnames(anomalyModelFeatures)

msstats_data_t = MSstatsConvert::MSstatsAnomalyScores(
  t_converted_data, anomalyModelFeatures,
  anomalyModelFeatureTemporal, .5, 100, runOrder, n_trees,
  max_depth, numberOfCores)

summarized_t = dataProcess(msstats_data_t,
                           normalization="equalizeMedians",
                           featureSubset = "topN",
                           n_top_feature = 100,
                           MBimpute=TRUE,
                           summaryMethod="linear",
                           maxQuantileforCensored = 0.99,
                           numberOfCores = 12)

msstats_data_t$AnomalyScores=NULL
summarized_t_msstats = dataProcess(msstats_data_t,
                                   normalization="equalizeMedians",
                                   featureSubset = "topN",
                                   n_top_feature = 100,
                                   MBimpute=TRUE,
                                   summaryMethod="TMP",
                                   numberOfCores = 12)

## Correct for instrument effects ----------------------------------------------
correct_instrument_effect <- function(dt) {
  dt[, LogIntensities_adjusted := {
    if (all(is.na(LogIntensities)) || length(unique(
      INSTRUMENT[!is.na(LogIntensities)])) == 1) {
      LogIntensities  # can't adjust if only one level
    } else {
      fit = lm(LogIntensities ~ INSTRUMENT)
      predicted = predict(fit, newdata = .SD)
      LogIntensities - (predicted - mean(LogIntensities, na.rm = TRUE))
    }
  }, by = Protein]
  return(dt)
}

# Proposed
proposed_summarized = list()
proposed_summarized$ProteinLevelData = rbindlist(
  list(summarized_t$ProteinLevelData, summarized_v$ProteinLevelData),
  use.names=TRUE)
proposed_summarized$FeatureLevelData = rbindlist(
  list(summarized_t$FeatureLevelData, summarized_v$FeatureLevelData),
)
proposed_summarized$ProteinLevelData = correct_instrument_effect(
  as.data.table(proposed_summarized$ProteinLevelData))
proposed_summarized$ProteinLevelData$LogIntensities = proposed_summarized$ProteinLevelData$LogIntensities_adjusted

# MSstats
msstats_summarized = list()
msstats_summarized$ProteinLevelData = rbindlist(
  list(msstats_summarized_t$ProteinLevelData, msstats_summarized_v$ProteinLevelData),
  use.names=TRUE)
msstats_summarized$FeatureLevelData = rbindlist(
  list(msstats_summarized_t$FeatureLevelData, msstats_summarized_v$FeatureLevelData),
)
msstats_summarized$ProteinLevelData = correct_instrument_effect(
  as.data.table(msstats_summarized$ProteinLevelData))
msstats_summarized$ProteinLevelData$LogIntensities = msstats_summarized$ProteinLevelData$LogIntensities_adjusted

## Differential analysis -------------------------------------------------------
model_instrument = function(data){
  
  data$ProteinLevelData = data$ProteinLevelData %>%
    filter(SEX != "" & TIME == "BASELINE" &
             GROUP != "STUDY 3")
  
  data$ProteinLevelData = data$ProteinLevelData %>%
    mutate(
      condition_bucket = fifelse(
        CDRSB < 4.5,
        "Low", "High"
      )
    )
  
  data$ProteinLevelData$GROUP = data$ProteinLevelData$condition_bucket

  msstats_model = groupComparison("pairwise", data,
                                  numberOfCores=4)
  return(msstats_model)
}

# MSstats+
proposed_model = model_instrument(proposed_summarized)
proposed_model = proposed_model$ComparisonResult

#MSstats
msstats_model = model_instrument(msstats_summarized)
msstats_model = msstats_model$ComparisonResult

# Analysis plots ---------------------------------------------------------------
# Create a dummy data frame for vertical lines
vline_data <- data.frame(
  cutoff = c(4.5, 9.5),
  severity = c("Mild", "Moderate")
)

annotation %>%
  filter(SEX != "", TIME == "BASELINE", STUDY != "STUDY 3") %>%
  ggplot(aes(x = CDRSB)) +
  geom_histogram(
    bins = 25,
    fill = "steelblue",
    color = "black",
    alpha = 0.7
  ) +
  geom_vline(
    data = vline_data,
    aes(xintercept = cutoff, color = severity),
    linetype = "dashed",
    linewidth = 2,
    show.legend = TRUE
  ) +
  scale_color_manual(
    name = "Severity Cutoff",
    values = c("Mild" = "darkorange", "Moderate" = "red")
  ) +
  annotate("text", x = 4.5, y = Inf, label = "4.5", vjust = -0.5, hjust = -0.1,
           size = 6, fontface = "bold", color = "darkorange") +
  annotate("text", x = 9.5, y = Inf, label = "9.5", vjust = -0.5, hjust = -0.1,
           size = 6, fontface = "bold", color = "red") +
  theme_minimal(base_size = 16) +
  labs(
    x = "CDR-SB Score",
    y = "Participant Count",
    title = "Baseline CDR-SB Distribution"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
    axis.title = element_text(face = "bold", size = 26),
    axis.text = element_text(size = 24),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20),
    panel.grid.minor = element_blank()
  )


proposed_model %>% ggplot(aes(x = log2FC, y = -log10(adj.pvalue))) +
  geom_point(aes(
    color = ifelse(adj.pvalue > 0.05, "insignificant",
                   ifelse(log2FC < 0, "down", "up"))
  ), size = 3) +
  scale_color_manual(
    values = c("significant" = "gray", 
               "down" = "blue", 
               "up" = "red")
  ) +
  geom_vline(xintercept = c(-.5, .5), linetype = "dashed") +  # typical log2FC thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +  # p-value threshold
  labs(
    x = "log2 Fold Change",
    y = "-log10(adj. p-value)",
    title = ""
  ) +
  theme_minimal(base_size = 24) +
  theme(legend.position = "none") + 
  xlim(-2.75,2.75) + ylim(0, 5)

msstats_model %>% ggplot(aes(x = log2FC, y = -log10(adj.pvalue))) +
  geom_point(aes(
    color = ifelse(adj.pvalue > 0.05, "insignificant",
                   ifelse(log2FC < 0, "down", "up"))
  ), size = 3) +
  scale_color_manual(
    values = c("significant" = "gray", 
               "down" = "blue", 
               "up" = "red")
  ) +
  geom_vline(xintercept = c(-.5, .5), linetype = "dashed") +  # typical log2FC thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +  # p-value threshold
  labs(
    x = "log2 Fold Change",
    y = "-log10(adj. p-value)",
    title = ""
  ) +
  theme_minimal(base_size = 24) +
  theme(legend.position = "none") + 
  xlim(-2.75,2.75) + ylim(0, 5)

cbPalette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Compute sets
only_proposed <- setdiff(proposed_proteins, msstats_proteins)
only_msstats <- setdiff(msstats_proteins, proposed_proteins)
overlap <- intersect(proposed_proteins, msstats_proteins)

proposed_model %>% filter(log2FC > .4 & Protein %in% only_proposed)

# Actual counts (for labeling)
actual_counts <- c(
  "MSstats+" = length(only_proposed),
  "MSstats" = length(only_msstats),
  "MSstats+&MSstats" = length(overlap)
)

# Transformed counts (for area)
transformed_counts <- sqrt(actual_counts)+2.5

# Fit eulerr with transformed sizes
fit <- euler(transformed_counts)

# Plot with real counts shown
plot(
  fit,
  fills = list(fill = c(cbPalette[[1]], cbPalette[[2]])),
  quantities = list(labels = actual_counts, cex = 3),
  labels = list(cex = 2.75),
  edges = list(lwd = 1.5)
)

# Profile plots 
plot_protein = function(data, protein){
  plot = ggplot() +
    geom_point(data=data$FeatureLevelData %>% filter(PROTEIN == protein),
               aes(x=Order, y=newABUNDANCE, group=FEATURE), color="grey", alpha=.6) +
    geom_line(data=data$FeatureLevelData %>% filter(PROTEIN == protein),
              aes(x=Order, y=newABUNDANCE, group=FEATURE), color="grey", alpha=.6) +
    geom_point(data=data$ProteinLevelData %>% filter(Protein == protein),
               aes(x=Order, y=LogIntensities, group=Protein), color="darkred", size=1.5) +
    geom_line(data=data$ProteinLevelData %>% filter(Protein == protein),
              aes(x=Order, y=LogIntensities, group=Protein), color="darkred", linewidth=1) +
    geom_errorbar(data=data$ProteinLevelData %>% filter(Protein == protein),
                  aes(x=Order, 
                      ymin = LogIntensities-sqrt(Variance), 
                      ymax = LogIntensities+sqrt(Variance)),
                  width = 2, linewidth=1, color="darkred") +
    theme_bw() + theme(legend.position="none",
                       # axis.text.x = element_blank(),
                       axis.title = element_text(size = 30),
                       axis.text = element_text(size = 26),
                       plot.title = element_text(size = 32)) +
    xlim(0, 125) + #ylim(5, 22) + 
    labs(title="Summarized Intensities", x="Run", y=expression(log[2] * " Intensity"))
  print(plot)
}

plot_protein(summarized_v, "P11597")
plot_protein(summarized_t, "P11597")
proposed_model %>% filter(Protein == "P11597")
msstats_model %>% filter(Protein == "P11597")

plot = ggplot() +
  geom_point(data=summarized_t$ProteinLevelData %>% filter(Protein == "P11597") ,
             aes(x=Order, y=Variance, group=Protein), color="darkred", size=3) +
  geom_line(data=summarized_t$ProteinLevelData %>% filter(Protein == "P11597") ,
            aes(x=Order, y=Variance, group=Protein), color="darkred", linewidth=1) +
  theme_bw() + theme(legend.position="none",
                     # axis.text.x = element_blank(),
                     axis.title = element_text(size = 30),
                     axis.text = element_text(size = 26),
                     plot.title = element_text(size = 32)) +
  xlim(0, 125) +
  labs(title="Summarized Intensities", x="Run", y="Variance")
print(plot)

summarized_t$ProteinLevelData %>% filter(Protein == "P11597") %>% 
  ggplot(
    aes(x=Order, y=Variance)
  ) + geom_line() + geom_point()

plot_protein(summarized_v, "Q14766")
plot_protein(summarized_t, "Q14766")
proposed_model$ComparisonResult %>% filter(Protein == "Q14766")
msstats_model$ComparisonResult %>% filter(Protein == "Q14766")
summarized_t$ProteinLevelData %>% filter(Protein == "Q14766")


msstats_data_t = msstats_data_t %>% merge(temporal_t, by.x="Run", by.y="Run")
msstats_data_t$PSM = paste(msstats_data_t$PeptideSequence,
                           msstats_data_t$PrecursorCharge,
                           sep="_")
msstats_data_t$FEATURE = paste(msstats_data_t$PeptideSequence,
                               msstats_data_t$PrecursorCharge,
                               msstats_data_t$FragmentIon,
                               msstats_data_t$ProductCharge,
                               sep="_")

psms = msstats_data_t %>% filter(ProteinName == "P11597") %>% distinct(PSM) %>% unlist()

# Intensities
msstats_data_t %>% 
  filter(PSM == psms[[6]]) %>% 
  ggplot() + 
  geom_line(aes(x=Order, y=AnomalyScores), 
            linewidth=1.25, color="darkorchid4") + 
  geom_point(aes(x=Order, y=AnomalyScores), color="darkorchid4", size=3) + 
  theme_bw()  + theme(
    legend.position = "bottom",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16), 
    # axis.text.x = element_blank(),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    plot.title = element_text(size = 32)
  ) + 
  xlim(0, 125) + ylim(0, .7) + 
  labs(title="Intensities", x="Run", y="")

anomaly_cols = c("FGShapeQualityScore(MS2)", 
                 "FGShapeQualityScore(MS1)", "EGApexRT", 
                 "FGShapeQualityScore(MS2).mean_decrease",
                 "FGShapeQualityScore(MS1).mean_decrease",
                 "EGApexRT.dispersion_increase")

for (y in anomaly_cols){
  
  plot = msstats_data_t %>% 
    filter(PSM == psms[[6]]) %>% 
    ggplot() + 
    geom_line(aes(x=Order, y=.data[[y]], group=PSM), 
              color="grey", linewidth=1.25) + 
    geom_point(aes(x=Order, y=.data[[y]]), color="grey", size=3) + 
    theme_bw()  + theme(legend.position="none",
                        axis.text.x = element_blank(),
                        axis.title = element_text(size = 18),
                        axis.text = element_text(size = 16),
                        plot.title = element_text(size = 30)) + 
    labs(title=y, x="Run", y=y)
  print(plot)
}

# Imputed values
plot = msstats_data_v %>% 
  filter(PSM == "_HPPEASVQIHQVSR__3") %>% 
  merge(summarized_v$FeatureLevelData, all.x = TRUE, all.y = FALSE, 
        by.x = c("FEATURE", "Run"), by.y = c("FEATURE", "originalRUN")) %>% 
  ggplot() + 
  geom_line(aes(x = Order.x, y = newABUNDANCE, group = FragmentIon), 
            linewidth = 1.25, color = "grey") + 
  geom_errorbar(aes(x = Order.x, ymin = newABUNDANCE - AnomalyScores,
                    ymax = newABUNDANCE + AnomalyScores), 
                width = 1, linewidth = 1.1, color = "red", show.legend = FALSE) +
  geom_point(aes(x = Order.x, y = newABUNDANCE, group = FragmentIon), 
             size = 3, color = "grey") + 
  geom_point(aes(x = Order.x, y = predicted, group = FragmentIon), 
             color = "white", size = 2) +
  scale_color_viridis_d(option = "D") + 
  theme_bw() + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16), 
    # axis.text.x = element_blank(),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    plot.title = element_text(size = 32)
  ) + 
  xlim(0, 125) + ylim(5, 22) + 
  labs(title = "dsafa", x = "Run", y=expression(log[2] * " Intensity"), color = "Fragment")
print(plot)
