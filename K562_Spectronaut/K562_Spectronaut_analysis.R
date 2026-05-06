## Packages --------------------------------------------------------------------
# Data analysis
library(data.table)
library(tidyverse)

# Comparison packages
library(MSstats)
library(MSstatsConvert)
library(QFeatures)
library(msqrob2)
library(iq)
library(limma)
library(limpa)
library(DEqMS)

## Load + prepare inputs -------------------------------------------------------
# Load data
data_folder = 
annotation = fread(paste(data_folder, "K562_annotation.csv", sep="/"))

# Prepare annotation info
annotation = annotation[annotation$Condition != "Blank",]
annotation$Run = annotation$R.FileName
run_order = unique(annotation[, .(Run, Order)])

# Predefined true positives/negatives
protein_swap_list = fread("K562_protein_swap_list.csv")

true_positives = unlist(
  protein_swap_list[protein_swap_list$Label == "Positive", "Proteins"])
true_negatives = unlist(
  protein_swap_list[protein_swap_list$Label == "Negative", "Proteins"])

color_p = c("#E69F00", "#56B4E9", "#009E73", 
            "#F0E442", "#0072B2", "#D55E00", 
            "#CC79A7")

## Analysis --------------------------------------------------------------------
## MSstats+ analysis -----------------------------------------------------------
# Load data
msstats_input = fread(file=paste(data_folder, "MSstats+", 
                                 "MSstats+_input.csv", sep="/"))
load(file=paste(data_folder, "MSstats+", "MSstats+_summarized.rda", sep="/"))
msstatsplus_model = fread(file=paste(data_folder, "MSstats+",
                                  "MSstats+_model.csv", sep="/"))

# Prepare input
msstats_input = merge(msstats_input, run_order, by="Run", all.x=TRUE, all.y=FALSE)
msstats_input$PSM = paste0(msstats_input$PeptideSequence, "_", msstats_input$PrecursorCharge)
msstats_input$Feature = paste(msstats_input$PeptideSequence, msstats_input$PrecursorCharge,
                              msstats_input$FragmentIon, msstats_input$ProductCharge, sep="_")

summarized$FeatureLevelData$Order = as.integer(str_split_i(
  summarized$FeatureLevelData$originalRUN, "Seq", 2))
summarized$ProteinLevelData$Order = as.integer(str_split_i(
  summarized$ProteinLevelData$originalRUN, "Seq", 2))

# Experiment-wide plots
# Skewness analysis
health_info = MSstatsConvert::CheckDataHealth(msstats_input)
skew_score = health_info[[2]]

ggplot(skew_score, aes(x = skew)) + 
  geom_histogram(fill = "#E69F00", color = "black", binwidth = 0.2) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1.5) + 
  theme_minimal(base_size = 16) + 
  labs(
    x = "Pearson's moment coefficient of skewness",
    y = "Count"
  ) +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)
  ) +
  xlim(-1, 7) +
  scale_y_continuous(expand = c(0, 0))

msstats_input$Anomaly = ifelse(
  msstats_input$Order > 30, TRUE, FALSE
)

# Anomaly vs non Anomaly plot
msstats_input %>%
  ggplot(aes(x = AnomalyScores, fill = Anomaly)) +
  geom_histogram(
    position = "identity", 
    alpha = 0.6,
    bins = 50,
    color = "black"
  ) +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8")) +
  labs(
    title = "Dataset 1: K562 
    Benchmark (Spectronaut)",
    x = "Anomaly Score",
    y = "Density",
    fill = "Anomaly"
  ) +
  xlim(0, .82) + 
  scale_y_continuous(labels = scales::scientific) +
  theme_minimal(base_size = 28) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

# Anomaly distribution
msstats_input %>%
  merge(
    summarized$FeatureLevelData,
    all.x = TRUE, all.y = FALSE,
    by.x = c("Feature", "Run"),
    by.y = c("FEATURE", "originalRUN")
  ) %>%
  filter(!is.na(predicted)) %>%
  ggplot(aes(x = AnomalyScores)) +
  geom_histogram(
    bins = 50,
    fill = "#4C72B0",
    color = "black",
    alpha = 0.7
  ) +
  labs(
    title = "Dataset 1: K562 Benchmark (Spectronaut)",
    x = "Anomaly Score",
    y = "Count"
  ) +
  theme_minimal(base_size = 28) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

# Anomaly vs Intensity
msstats_input %>%
  ggplot(aes(x = log2(Intensity), y = AnomalyScores)) +
  geom_hex(bins = 50) +
  scale_fill_gradientn(
    colours = c("grey90", "#fcae91", "#fb6a4a", "#cb181d"), 
    values = scales::rescale(c(0, 0.3, 0.6, 1)),
    trans = "sqrt",
    guide = "none"
  ) +
  theme_minimal(base_size = 16) +
  labs(
    x = expression(log[2] * " Intensity"),
    y = "Anomaly Score",
    title = ""
  ) +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 26),
    panel.grid.major = element_line(size = 0.3),
    panel.grid.minor = element_blank()
  )

# Individual profile plot (Update protein/PSM names to plot different cases)
# A4D126
anomaly_cols = c("FGShapeQualityScore(MS2)",
                 "FGShapeQualityScore(MS2).mean_decrease", 
                 "FGShapeQualityScore(MS1)",
                 "FGShapeQualityScore(MS1).mean_decrease", 
                 "EGDeltaRT", "EGDeltaRT.dispersion_increase")
new_cols = c("MS2 Shape Quality", "MS2 Shape Quality Trend",
             "MS1 Shape Quality", "MS1 Shape Quality Trend",
             "Delta Retention Time", "Delta RT Trend")
y_label = c("Shape Quality", "Negative CUSUM",
            "Shape Quality", "Negative CUSUM",
            "Delta RT", "Dispersion CUSUM")

setnames(msstats_input, anomaly_cols, new_cols)

for (y in seq_along(new_cols)){
  
  plot = msstats_input %>% 
    filter(PeptideSequence == "VAGQDGSVVQFK" & PrecursorCharge == 2) %>% 
    ggplot() + 
    geom_line(aes(x=Order, y=.data[[new_cols[[y]]]], group=PSM), 
              color="darkgrey", linewidth=1.25) + 
    geom_point(aes(x=Order, y=.data[[new_cols[[y]]]]), color="darkgrey", size=3) + 
    theme_bw()  + theme(legend.position="none",
                        # axis.text.x = element_blank(),
                        axis.title = element_text(size = 30),
                        axis.text = element_text(size = 28),
                        plot.title = element_text(size = 38)) + 
    labs(title=new_cols[[y]], x="Run (Ordered by acquisition time)", 
         y=y_label[[y]])
  print(plot)
}

msstats_input %>% 
  filter(PeptideSequence == "VAGQDGSVVQFK" & PrecursorCharge == 2) %>% 
  ggplot() + 
  geom_line(aes(x=Order, y=`AnomalyScores`, group=PrecursorCharge), 
            color="darkorchid4", linewidth=1.25) + 
  geom_point(aes(x=Order, y=`AnomalyScores`), color="darkorchid4", size=3) + 
  theme_bw()  + theme(legend.position="none",
                      # axis.text.x = element_blank(),
                      axis.title = element_text(size = 30),
                      axis.text = element_text(size = 28),
                      plot.title = element_text(size = 38)) + 
  labs(title="Precusor Anomaly Scores", x="Run (Ordered by acquisition time)", 
       y="Anomaly Score")

msstats_input %>% 
  filter(PeptideSequence == "VAGQDGSVVQFK" & PrecursorCharge == 2) %>% 
  ggplot() + 
  geom_line(aes(x=Order, y=log2(Intensity), group=FragmentIon), 
            color="darkred", linewidth=1.25) + 
  geom_point(aes(x=Order, y=log2(Intensity)), color="darkred", size=3) + 
  theme_bw()  + theme(legend.position="none",
                      # axis.text.x = element_blank(),
                      axis.title = element_text(size = 30),
                      axis.text = element_text(size = 28),
                      plot.title = element_text(size = 38)) + 
  labs(title="Fragment Intensities", x="Run (Ordered by acquisition time)", 
       y=expression(log[2] * " Intensity"))



msstats_input %>% 
  filter(ProteinName == "P55854") %>% 
  ggplot() + 
  geom_line(aes(x=Order, y=`AnomalyScores`, group=PSM, color=PSM), linewidth=1.25) + 
  geom_point(aes(x=Order, y=`AnomalyScores`, group=PSM, color=PSM), size=3) + 
  # scale_color_manual(values = color_p) +
  scale_color_viridis_d() + 
  theme_bw()  + theme(legend.position="none",
                      # axis.text.x = element_blank(),
                      axis.title = element_text(size = 30),
                      axis.text = element_text(size = 28),
                      plot.title = element_text(size = 38)) + 
  labs(title="Anomaly Scores Across Precursors", x="Run (Ordered by acquisition time)", 
       y="Anomaly Score")

protein = "P55854"
ggplot() +
  geom_point(data=summarized$FeatureLevelData %>%
               filter(PROTEIN == protein),
             aes(x=Order, y=newABUNDANCE, group=FEATURE), color="grey", alpha=.6) +
  geom_line(data=summarized$FeatureLevelData %>%
              filter(PROTEIN == protein),
            aes(x=Order, y=newABUNDANCE, group=FEATURE), color="grey", alpha=.6) +
  geom_point(data=summarized$ProteinLevelData %>%
               filter(Protein == protein),
             aes(x=Order, y=LogIntensities, group=Protein), color="darkred", size=3) +
  geom_line(data=summarized$ProteinLevelData %>%
              filter(Protein == protein),
            aes(x=Order, y=LogIntensities, group=Protein), color="darkred", linewidth=1.25) +
  geom_errorbar(data=summarized$ProteinLevelData %>%
                  filter(Protein == protein),
                aes(x=Order, ymin = LogIntensities-sqrt(Variance),
                    ymax = LogIntensities+sqrt(Variance)),
                width = 1, linewidth=1.1, color="darkred") + 
  theme_bw()  + theme(legend.position="none",
                      # axis.text.x = element_blank(),
                      axis.title = element_text(size = 30),
                      axis.text = element_text(size = 28),
                      plot.title = element_text(size = 38)) + 
  labs(title="Summarized Protein-level Intensities", 
       x="Run (Ordered by acquisition time)", 
       y=expression(log[2] * " Intensity"))

# Missing value plot
msstats_psm = "EAGGAFGK_1"
plot = msstats_input %>% 
  filter(PSM == msstats_psm) %>% 
  merge(summarized$FeatureLevelData, all.x = TRUE, all.y = FALSE, 
        by.x = c("Feature", "Run"), by.y = c("FEATURE", "originalRUN")) %>% 
  ggplot() + 
  geom_line(aes(x = Order.x, y = newABUNDANCE, group = FragmentIon, color = FragmentIon), 
            linewidth = 1.25) + 
  geom_point(aes(x = Order.x, y = newABUNDANCE, group = FragmentIon, color = FragmentIon), 
             size = 3) + 
  geom_point(aes(x = Order.x, y = predicted, group = FragmentIon), 
             color = "white", size = 2) +
  scale_color_manual(values=c("b6" = "#009E73",
                              "y3" = "#D55E00",
                              "y7" = "#0072B2")) + 
  theme_bw() + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16), 
    axis.text.x = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 32)
  ) + 
  labs(title = msstats_psm, x = "Run", y = "Log2 Intensity", color = "Fragment")
print(plot)

plot = msstats_input %>% 
  filter(PSM == msstats_psm) %>% 
  ggplot() + 
  geom_line(aes(x = Order, y = AnomalyScores, group = PSM), color = "darkorchid4", linewidth = 1.25) + 
  geom_point(aes(x = Order, y = AnomalyScores, group = PSM), color = "darkorchid4", 
             size = 3) +
  theme_bw() + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16), 
    axis.text.x = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 32)
  ) + 
  labs(title = msstats_psm, x = "Run", y = "Log2 Intensity", color = "Fragment")
print(plot)

## Compare to other models -----------------------------------------------------
msstats_model = fread(file=paste("MSstats", 
                                 "MSstats_model.csv", sep="/"))

msqrob2_model = fread(file=paste("msqrob2", 
                                 "msqrob2_model.csv", sep="/"))

limma_model = fread(file=paste("limma", 
                               "limma_model.csv", sep="/"))

limpa_model = fread(file=paste("limpa",
                               "limpa_model.csv", sep="/"))

deqms_model = fread(file=paste("DEqMS", 
                               "deqms_model.csv", sep="/"))

mapdia_model = fread(file=paste("mapDIA", 
                                "mapdia_model.csv", sep="/"))

# number of proteins tested
msstatsplus_model %>% filter(Protein != "") %>% 
  group_by(Label) %>% summarize(tested = n_distinct(Protein))
msstats_model %>% filter(Protein != "") %>% 
  group_by(Label) %>% summarize(tested = n_distinct(Protein))
msqrob2_model %>% filter(Protein != "") %>% 
  group_by(Label) %>% summarize(tested = n_distinct(Protein))
limma_model %>% filter(Protein != "") %>% 
  group_by(Label) %>% summarize(tested = n_distinct(Protein))
limpa_model %>% filter(Protein != "") %>% 
  group_by(Label) %>% summarize(tested = n_distinct(Protein))
deqms_model %>% filter(Protein != "") %>% 
  group_by(Label) %>% summarize(tested = n_distinct(Protein))
mapdia_model %>% filter(Protein != "") %>% 
  group_by(Label) %>% summarize(tested = n_distinct(Protein))

# Calculate TPR/FDR
# TPR
msstatsplus_model %>% filter(pvalue < .05 & is.finite(log2FC) &
                                    Label == "Positive") %>% nrow() /
  msstatsplus_model %>% filter(is.finite(log2FC) & Label == "Positive") %>% nrow()

msstats_model %>% filter(pvalue < .05 & is.finite(log2FC) & Label == "Positive") %>% nrow() /
  msstats_model %>% filter(is.finite(log2FC) & Label == "Positive") %>% nrow()

limpa_model %>% filter(pvalue < .05 & is.finite(logFC) & Label == "Positive") %>% nrow() /
  limpa_model %>% filter(is.finite(logFC) & Label == "Positive") %>% nrow()

limma_model %>% filter(pvalue < .05 & is.finite(logFC) & Label == "Positive") %>% nrow() /
  limma_model %>% filter(is.finite(logFC) & Label == "Positive") %>% nrow()

msqrob2_model %>% filter(pval < .05 & is.finite(logFC) & Label == "Positive") %>% nrow() /
  msqrob2_model %>% filter(is.finite(logFC) & Label == "Positive") %>% nrow()

deqms_model %>% filter(pvalue < .05 & is.finite(logFC) & Label == "Positive") %>% nrow() /
  deqms_model %>% filter(is.finite(logFC) & Label == "Positive") %>% nrow()

# FDR is adjusted (non comparable with the other models)
mapdia_model %>% filter(FDR < .05 & is.finite(log2FC) & Label == "Positive") %>% nrow() /
  mapdia_model %>% filter(is.finite(log2FC) & Label == "Positive") %>% nrow()

# PPV
1 - msstatsplus_model %>% filter(pvalue < .05 & is.finite(log2FC) &
                                    Label == "Negative") %>% nrow() /
  msstatsplus_model %>% filter(pvalue < .05 & is.finite(log2FC)) %>% nrow()

1 - msstats_model %>% filter(pvalue < .05 & is.finite(log2FC) & Label == "Negative") %>% nrow() /
  msstats_model %>% filter(pvalue < .05 & is.finite(log2FC)) %>% nrow()

1 - limpa_model %>% filter(pvalue < .05 & is.finite(logFC) & Label == "Negative") %>% nrow() /
  limpa_model %>% filter(pvalue < .05 & is.finite(logFC)) %>% nrow()

1 - limma_model %>% filter(pvalue < .05 & is.finite(logFC) & Label == "Negative") %>% nrow() /
  limma_model %>% filter(pvalue < .05 & is.finite(logFC)) %>% nrow()

1 - msqrob2_model %>% filter(pval < .05 & is.finite(logFC) & Label == "Negative") %>% nrow() /
  msqrob2_model %>% filter(pval < .05 & is.finite(logFC)) %>% nrow()

1 - deqms_model %>% filter(pvalue < .05 & is.finite(logFC) & Label == "Negative") %>% nrow() /
  deqms_model %>% filter(pvalue < .05 & is.finite(logFC)) %>% nrow()

1 - mapdia_model %>% filter(FDR < .05 & is.finite(log2FC) & Label == "Negative") %>% nrow() /
  mapdia_model %>% filter(FDR < .05 & is.finite(log2FC)) %>% nrow()

# Comparison plots
plot_dt = rbindlist(list(
  data.table(model = "Proposed",
             adj_p = msstatsplus_model$adj.pvalue,
             pval = msstatsplus_model$pvalue,
             logFC = msstatsplus_model$log2FC),
  
  data.table(model = "MSstats",
             adj_p = msstats_model$adj.pvalue,
             pval = msstats_model$pvalue,
             logFC = msstats_model$log2FC),
  
  data.table(model = "msqrob2",
             adj_p = msqrob2_model$adjPval,
             pval = msqrob2_model$pval,
             logFC = msqrob2_model$logFC),
  
  data.table(model = "MaxLFQ + limma",
             adj_p = limma_model$adj.pvalue,
             pval = limma_model$pvalue,
             logFC = limma_model$logFC),
  
  data.table(model = "DEqMS",
             adj_p = deqms_model$adj.pvalue,
             pval = deqms_model$pvalue,
             logFC = deqms_model$logFC),
  
  data.table(model = "limpa",
             adj_p = limpa_model$adj.pvalue,
             pval = limpa_model$pvalue,
             logFC = limpa_model$logFC)
  ), fill = TRUE)

model_colors = c(
  "Proposed" = color_p[[1]],
  "MSstats" = color_p[[2]],
  "msqrob2" = color_p[[3]],
  "MaxLFQ + limma" = color_p[[4]],
  "DEqMS" = color_p[[5]],
  "limpa" = color_p[[6]]
)

plot_dt[, model := factor(
  model,
  levels = c(
    "DEqMS",
    "limpa",
    "MaxLFQ + limma",
    "msqrob2",
    "MSstats",
    "Proposed"
  )
)]

ggplot(plot_dt[!is.na(adj_p)]) +
  geom_histogram(
    aes(x = adj_p, fill = model),
    bins = 40,
    color = "black"
  ) +
  geom_vline(
    xintercept = 0.05,
    linetype = "dashed",
    color = "red",
    linewidth = 2
  ) +
  facet_grid(model ~ ., scales = "free_y") +
  coord_cartesian(xlim = c(-0.05, 1.05)) +
  scale_fill_manual(values = model_colors, drop = FALSE) +
  scale_y_log10() +
  theme_bw(base_size = 18) +
  labs(
    title = "Adjusted p-value distribution",
    x = "Adjusted p-value",
    y = expression(log[10]~count)
  ) + 
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 24),
    strip.text.y = element_text(angle = 0, size = 22, face = "bold"),
    axis.title.x = element_text(size = 24),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_blank(),
    text = element_text(size = 22),
    strip.background = element_blank(),
    title = element_text(size = 26)
  )

ggplot(plot_dt[!is.na(adj_p)]) +
  geom_histogram(
    aes(x = pval, fill = model),
    bins = 40,
    color = "black"
  ) +
  geom_vline(
    xintercept = 0.05,
    linetype = "dashed",
    color = "red",
    linewidth = 2
  ) +
  facet_grid(model ~ ., scales = "free_y") +
  coord_cartesian(xlim = c(-0.05, 1.05)) +
  scale_fill_manual(values = model_colors, drop = FALSE) +
  scale_y_log10() +
  theme_bw(base_size = 18) +
  labs(
    title = "Nominal p-value distribution",
    x = "Nominal p-value",
    y = expression(log[10]~count)
  ) + 
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 24),
    strip.text.y = element_text(angle = 0, size = 22, face = "bold"),
    axis.title.x = element_text(size = 24),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_blank(),
    text = element_text(size = 22),
    strip.background = element_blank(),
    title = element_text(size = 26)
  )


ggplot() +
  geom_boxplot(data=msstatsplus_model %>% filter(Label == "Positive"),
               aes(x = "Proposed", y = log2FC),
               fill = color_p[[1]], color = "black") +
  geom_boxplot(data=msstats_model %>% filter(Label == "Positive"),
               aes(x = "MSstats", y = log2FC),
               fill = color_p[[2]], color = "black") +
  geom_boxplot(data=msqrob2_model %>% filter(Label == "Positive"),
               aes(x = "msqrob2", y = logFC),
               fill = color_p[[3]], color = "black") +
  geom_boxplot(data=limma_model %>% filter(Label == "Positive"),
               aes(x = "MaxLFQ + limma", y = logFC),
               fill = color_p[[4]], color = "black") +
  geom_boxplot(data=deqms_model %>% filter(Label == "Positive"),
               aes(x = "DEqMS", y = logFC),
               fill = color_p[[5]], color = "black") +
  geom_boxplot(data=limpa_model %>% filter(Label == "Positive"),
               aes(x = "limpa", y = logFC),
               fill = color_p[[6]], color = "black") +
  # geom_boxplot(data=mapdia_model %>% filter(Label == "Positive"),
  #              aes(x = "mapDIA", y = -log2FC),
  #              fill = color_p[[7]], color = "black") +
  geom_hline(aes(yintercept = -1), color="red", linetype="dashed", size=1.5) + 
  scale_fill_viridis_d(option = "plasma") + 
  ylim(-2.1, 1) + 
  theme_minimal(base_size = 28) +  # Make text readable
  labs(
    title = expression("K562 benchmark - TP" ~ log[2] ~ FC),
    x = "",
    y = expression(log[2] * " Fold Change")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size=44),
    axis.title.y = element_text( size=32),
    axis.text.x = element_text(angle = 45, hjust = 1, size=30),
    axis.text.y = element_text(size=30),
    legend.position = "none",
    legend.title = element_text(face = "bold")
  )

# Fold change negative
ggplot() +
  geom_boxplot(data=msstatsplus_model %>% filter(Label == "Negative"), 
               aes(x = "Proposed", y = log2FC), 
               fill = color_p[[1]], color = "black") +
  geom_boxplot(data=msstats_model %>% filter(Label == "Negative"),
               aes(x = "MSstats", y = log2FC),
               fill = color_p[[2]], color = "black") +
  geom_boxplot(data=msqrob2_model %>% filter(Label == "Negative"),
               aes(x = "msqrob2", y = logFC),
               fill = color_p[[3]], color = "black") +
  geom_boxplot(data=limma_model %>% filter(Label == "Negative"),
               aes(x = "MaxLFQ + limma", y = logFC),
               fill = color_p[[4]], color = "black") +
  geom_boxplot(data=deqms_model %>% filter(Label == "Negative"), 
               aes(x = "DEqMS", y = logFC), 
               fill = color_p[[5]], color = "black") +
  geom_boxplot(data=limpa_model %>% filter(Label == "Negative"),
               aes(x = "limpa", y = logFC),
               fill = color_p[[6]], color = "black") +
  # geom_boxplot(data=mapdia_model %>% filter(Label == "Negative"),
  #              aes(x = "mapDIA", y = -log2FC),
  #              fill = color_p[[7]], color = "black") +
  geom_hline(aes(yintercept = 0), color="red", linetype="dashed", size=1.5) + 
  scale_fill_viridis_d(option = "plasma") + 
  # ylim(-1, 1) + 
  theme_minimal(base_size = 28) +  # Make text readable
  labs(
    title = expression("K562 benchmark - FP" ~ log[2] ~ FC),
    x = "",
    y = expression(log[2] * " Fold Change")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size=44),
    axis.title.y = element_text( size=32),
    axis.text.x = element_text(angle = 45, hjust = 1, size=30),
    axis.text.y = element_text(size=30),
    legend.position = "none",
    legend.title = element_text(face = "bold")
  )

# Standard error
ggplot() +
  geom_boxplot(data=msstatsplus_model %>% filter(Label == "Positive"),
               aes(x = "Proposed", y = SE),
               fill = color_p[[1]], color = "black") +
  geom_boxplot(data=msstats_model %>% filter(Label == "Positive"),
               aes(x = "MSstats", y = SE),
               fill = color_p[[2]], color = "black") +
  geom_boxplot(data=msqrob2_model %>% filter(Label == "Positive"),
               aes(x = "msqrob2", y = se),
               fill = color_p[[3]], color = "black") +
  geom_boxplot(data=limma_model %>% filter(Label == "Positive"),
               aes(x = "MaxLFQ + limma", y = SE),
               fill = color_p[[4]], color = "black") +
  geom_boxplot(data=deqms_model %>% filter(Label == "Positive"),
               aes(x = "DEqMS", y = SE),
               fill = color_p[[5]], color = "black") +
  geom_boxplot(data=limpa_model %>% filter(Label == "Positive"),
               aes(x = "limpa", y = SE),
               fill = color_p[[6]], color = "black") +
  ylim(0, 1) + 
  theme_minimal(base_size = 28) +  # Make text readable
  labs(
    title = "K562 benchmark - SE",
    x = "",
    y = expression(log[2] * " Fold Change")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size=44),
    axis.title.y = element_text( size=32),
    axis.text.x = element_text(angle = 45, hjust = 1, size=30),
    axis.text.y = element_text(size=30),
    legend.position = "none",
    legend.title = element_text(face = "bold")
  )
