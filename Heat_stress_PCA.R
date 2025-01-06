library(ggfortify)  
library(tidyverse)
library(lme4)
library(brms)
library(missMDA)
library(FactoMineR)
library(ggpubr)
library(ggrepel)
library(factoextra)
library(forcats)

# Set working directory
setwd("~/Documents/Rice_Research/")

# Load trait data
solidity <- read.csv("~/Downloads/USDA_seed_scan_phenotypes.csv") %>%
  select(., c(GSOR, Solidity))
names(solidity) <- c("Samples", "Solidity")

biomass <- read.csv("~/Documents/Rice_Research/biomass_data.csv") %>%
  select(., c(X, weight_after..g.))
names(biomass) <- c("Samples", "Biomass")

seed <- read.csv("~/Documents/Rice_Research/Seeds-count.csv", header = T) 

fl <- read.csv("~/Documents/Rice_Research/flowering_data_USDA.csv", header = F)
names(fl) <- c("Samples", "flowering_time", "sub_pop", "chr3", "chr7", "chromosome")

HS_meta_data <- fl %>%
  separate(Samples, into = c("Samples", "Treatment"), sep = "-") %>%
  separate(Treatment, into = c("Block", "Treatment"), sep = "(?<=\\d)(?=\\D)") %>%
  select(Samples, Block, Treatment, sub_pop, chr3, chr7, chromosome)

write.csv(HS_meta_data, file = "~/Documents/Rice_Research/HS_meta_data.csv")

# Perform the joins
merged_df <- seed %>%
  full_join(fl, by = "Samples") %>%
  #inner_join(solidity, by = "Samples") %>%
  full_join(biomass, by = "Samples")

merged_df <- merged_df %>%
  separate(Samples, into = c("Samples", "Treatment"), sep = "-") %>%
  separate(Treatment, into = c("Block", "Treatment"), sep = "(?<=\\d)(?=\\D)")

#solidity$Samples <- as.character(solidity$Samples)

#merged_df <- merged_df %>%
 # inner_join(solidity, by = "Samples")

# Drop the overall averages to avoid conflicts
seed_clean <- merged_df %>%
  select(-Avg_seed_per_pan, -Avg_spikelet_fertility, -Avg_spikelet_no, -Avg_seed_weight)

seed_clean$second_Avg_spikelet_fertility <- as.numeric(seed_clean$second_Avg_spikelet_fertility)
seed_clean$third_Avg_spikelet_fertility <- as.numeric(seed_clean$third_Avg_spikelet_fertility)


# Reshape the dataframe into long format
seed_long <- seed_clean %>%
  pivot_longer(
    cols = c(first_Avg_spikelet_fertility, second_Avg_spikelet_fertility, third_Avg_spikelet_fertility,
             first_Avg_spikelet_no, second_Avg_spikelet_no, third_Avg_spikelet_no,
             first_Avg_seed_per_pan, second_Avg_seed_per_pan, third_Avg_seed_per_pan,
             first_avg_seed_weight, second_avg_seed_weight, third_avg_seed_weight),
    names_to = c("group", ".value"),
    names_pattern = "(.+?)_(.+)"
  )

# seed_meta <- seed_clean %>%
  # select(Samples,Block, Treatment, sub_pop, chr3, chr7, chromosome)
# write.csv(seed_meta, file = "./Documents/Rice_Research/HS_meta_data.csv")

# View the transformed dataframe
head(seed_long)

seed_long <- as.data.frame(seed_long)

seed_long <- seed_long %>%
  dplyr::rename(
    Spikelet_fertility = Avg_spikelet_fertility,
    Seed_per_panicle = Avg_seed_per_pan,
    Spikelet_no = Avg_spikelet_no,
    Seed_weight = avg_seed_weight,
    panicle = group
  )


# Convert column A to numeric
seed_long$Spikelet_fertility <- as.numeric(seed_long$Spikelet_fertility)

# Rename Treatment column 
seed_long <- seed_long %>%
  mutate(Treatment = recode(Treatment, "C" = "Control", "H" = "Heat"))

# Convert chr3 to a factor and set levels to ensure "GG" appears before "AA"
seed_long <- seed_long %>%
  mutate(chr3 = factor(chr3, levels = c("GG", "AA")))

#------------------------------------------------------------------------
# mixed model regression with PCs
# Select relevant columns for PCA for chr3 for GG allele
chr3_data_GG <- seed_long %>%
  select(Biomass, flowering_time, Spikelet_fertility, Seed_per_panicle, Spikelet_no, Seed_weight, Treatment, chromosome, chr3, chr7, sub_pop, Samples, Block, panicle) %>%
  filter(Treatment %in% c("Control", "Heat")) %>%
  filter(chromosome %in% "chr3") %>%
  filter(chr3 %in% "GG")


pca_data <- chr3_data_GG %>%
  select(., Biomass, flowering_time, Spikelet_fertility, Seed_per_panicle, Spikelet_no, Seed_weight)
  
traits_standardized <- scale(pca_data)

# Perform PCA with missing values
pca_result <- PCA(traits_standardized, ncp = 5, scale.unit = TRUE, graph = FALSE)

# View results
summary(pca_result)

# Extract the first 3 PCs from the individual coordinates
pca_scores <- pca_result$ind$coord[, 1:3]

# Preview the extracted scores
head(pca_scores)

# Create a dataframe with PC scores
pca_data <- as.data.frame(pca_scores)
pca_data$Treatment <- chr3_data_GG$Treatment
pca_data$sub_pop <- chr3_data_GG$sub_pop
pca_data$Samples <- chr3_data_GG$Samples
pca_data$allele <- chr3_data_GG$chr3
pca_data$panicle <- chr3_data_GG$panicle
pca_data$Block <- chr3_data_GG$Block

pca_data$allele <- as.character(pca_data$allele)

pca_data <- cbind(pca_data, chr3_data_GG[, 1:6]) 


# Add the treatment data to pca_result
pca_result$ind$Treatment <- pca_data$Treatment
pca_result$ind$allele <- pca_data$allele


# Now create the plot
pc <- fviz_pca_biplot(pca_result, 
                axes = c(1, 2),  # Plot PC1 and PC2
                col.var = "red",  # Color arrows
                col.ind = pca_data$Treatment,
                geom = c("point"),  # Plot only points, no text labels
                pointsize = 3,  # Increase point size
                repel = TRUE,  # Avoid text overlap
                addEllipses = TRUE,  # Add confidence ellipses
                legend.title = "Treatment") +  # Map color to Treatment
  theme_minimal() +
  labs(
    title = "GG",
    x = paste0("PC1 (", round(pca_result$eig[1, 2], 1), "%)"),
    y = paste0("PC2 (", round(pca_result$eig[2, 2], 1), "%)")
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Title customization
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),   # Axis text size
    legend.position = "right"              # Position the legend on the right
  ) +
  scale_color_manual(values = c("Control" = "darkgreen", "Heat" = "tomato3")) #+
  #facet_wrap(~ pca_data$allele)  # Facet by allele (GG and AA)# Custom color scale

#-------------------------------------------------------------------------------
# Select relevant columns for PCA for chr3 for GG allele
chr3_data_AA <- seed_long %>%
  select(Biomass, flowering_time, Spikelet_fertility, Seed_per_panicle, Spikelet_no, Seed_weight, Treatment, chromosome, chr3, chr7, sub_pop, Samples, Block, panicle) %>%
  filter(Treatment %in% c("Control", "Heat")) %>%
  filter(chromosome %in% "chr3") %>%
  filter(chr3 %in% "AA")


pca_data <- chr3_data_AA %>%
  select(., Biomass, flowering_time, Spikelet_fertility, Seed_per_panicle, Spikelet_no, Seed_weight)

traits_standardized <- scale(pca_data)

# Perform PCA with missing values
pca_result <- PCA(traits_standardized, ncp = 5, scale.unit = TRUE, graph = FALSE)

# View results
summary(pca_result)

# Extract the first 3 PCs from the individual coordinates
pca_scores <- pca_result$ind$coord[, 1:3]

# Preview the extracted scores
head(pca_scores)

# Create a dataframe with PC scores
pca_data <- as.data.frame(pca_scores)
pca_data$Treatment <- chr3_data_AA$Treatment
pca_data$sub_pop <- chr3_data_AA$sub_pop
pca_data$Samples <- chr3_data_AA$Samples
pca_data$allele <- chr3_data_AA$chr3
pca_data$panicle <- chr3_data_AA$panicle
pca_data$Block <- chr3_data_AA$Block

pca_data$allele <- as.character(pca_data$allele)

pca_data <- cbind(pca_data, chr3_data_AA[, 1:6]) 


# Add the treatment data to pca_result
pca_result$ind$Treatment <- pca_data$Treatment
pca_result$ind$allele <- pca_data$allele


# Now create the plot
pc <- fviz_pca_biplot(pca_result, 
                      axes = c(1, 2),  # Plot PC1 and PC2
                      col.var = "red",  # Color arrows
                      col.ind = pca_data$Treatment,
                      geom = c("point"),  # Plot only points, no text labels
                      pointsize = 3,  # Increase point size
                      repel = TRUE,  # Avoid text overlap
                      addEllipses = TRUE,  # Add confidence ellipses
                      legend.title = "Treatment") +  # Map color to Treatment
  theme_minimal() +
  labs(
    title = "AA",
    x = paste0("PC1 (", round(pca_result$eig[1, 2], 1), "%)"),
    y = paste0("PC2 (", round(pca_result$eig[2, 2], 1), "%)")
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Title customization
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),   # Axis text size
    legend.position = "right"              # Position the legend on the right
  ) +
  scale_color_manual(values = c("Control" = "darkgreen", "Heat" = "tomato3")) #+
#facet_wrap(~ pca_data$allele)  # Facet by allele (GG and AA)# Custom color scale


#------------------------------------------------------------------------
chr3_data <- seed_long %>%
  select(Biomass, flowering_time, Spikelet_fertility, Seed_per_panicle, Spikelet_no, Seed_weight, Treatment, chromosome, chr3, chr7, sub_pop, Samples, Block, panicle) %>%
  filter(Treatment %in% c("Control", "Heat")) %>%
  filter(chromosome %in% "chr3") 

pca_data <- chr3_data %>%
  select(., Biomass, flowering_time, Spikelet_fertility, Seed_per_panicle, Spikelet_no, Seed_weight)

traits_standardized <- scale(pca_data)

# Perform PCA with missing values
pca_result <- PCA(traits_standardized, ncp = 5, scale.unit = TRUE, graph = FALSE)

# View results
summary(pca_result)

# Extract the first 3 PCs from the individual coordinates
pca_scores <- pca_result$ind$coord[, 1:3]

# Preview the extracted scores
head(pca_scores)

# Create a dataframe with PC scores
pca_data <- as.data.frame(pca_scores)
pca_data$Treatment <- chr3_data$Treatment
pca_data$sub_pop <- chr3_data$sub_pop
pca_data$Samples <- chr3_data$Samples
pca_data$allele <- chr3_data$chr3
pca_data$panicle <- chr3_data$panicle
pca_data$Block <- chr3_data$Block

pca_data$allele <- as.character(pca_data$allele)

pca_data <- cbind(pca_data, chr3_data[, 1:6]) 


# Fit the multivariate model PCs in Chr3
pca_data$allele <- factor(pca_data$allele, levels = c("GG", "AA"))
pca_data$sub_pop <- factor(pca_data$sub_pop, levels = c("TRJ", "AUS"))

model_brms <- brm(
  formula = mvbind(Dim.1, Dim.2, Dim.3) ~ Treatment * allele * sub_pop + (1 | Block) + (1 | Samples) + (1 | panicle),
  data = pca_data,
  family = gaussian()
)

# Summarize the model
summary(model_brms)


# Extract fixed effects (estimates, standard errors, and confidence intervals)
fixed_effects <- fixef(model_brms, summary = TRUE)

# Convert to a dataframe for ggplot2
fixed_effects_df <- as.data.frame(fixed_effects)
fixed_effects_df$Term <- rownames(fixed_effects_df)

# Rename columns for clarity
colnames(fixed_effects_df) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Term")

# Filter and separate the Term column using the correct separator
fixed_effects_df <- fixed_effects_df %>%
  filter((`l-95% CI` > 0 & `u-95% CI` > 0) | (`l-95% CI` < 0 & `u-95% CI` < 0)) #%>%
  #separate(Term, into = c("Response", "Interaction"), sep = "_", fill = "right")

# Remove "PC_" from the Term column
fixed_effects_df$Term <- gsub("PC\\d_", "", fixed_effects_df$Term)


# Plot fixed effects with confidence intervals
p1 <- ggplot(fixed_effects_df, aes(x = Estimate, y = Term)) +
  geom_point(colour = "blue", size = 3) +  # Change point color to blue and increase point size
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`), height = 0.2, colour = "blue", size = 1.5) +  # Change error bar color to blue and increase line thickness
  theme_minimal() +
  #facet_grid(.~ Response, scales = "free") + 
  xlab("Effect Estimate") + 
  ylab("Interaction") +
  ggtitle("Fixed Effects with 95% Confidence Intervals for Chr3") +
  theme(
    axis.text.x = element_text(size = 15),  # Increase x-axis text size
    axis.text.y = element_text(size = 15),  # Increase y-axis text size
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16)   # Increase y-axis title size
  )



#------------------------------------------------------------------------
# Fit the multivariate model for each reproductive trait for chromosome 3
model_flowering <- brm(
  formula = flowering_time ~ Treatment * allele * sub_pop + (1 | Block) + (1 | Samples) + (1 | panicle),
  data = pca_data,
  family = gaussian()
)

# Summarize the model
summary(model_flowering)

# Extract fixed effects (estimates, standard errors, and confidence intervals)
fixed_effects <- fixef(model_flowering, summary = TRUE)

# Convert to a dataframe for ggplot2
fixed_effects_df <- as.data.frame(fixed_effects)
fixed_effects_df$Term <- rownames(fixed_effects_df)

# Rename columns for clarity
colnames(fixed_effects_df) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Term")

# Filter and separate the Term column using the correct separator
fixed_effects_df_1 <- fixed_effects_df %>%
  filter((`l-95% CI` > 0 & `u-95% CI` > 0) | (`l-95% CI` < 0 & `u-95% CI` < 0)) #%>%
#separate(Term, into = c("Response", "Interaction"), sep = "_", fill = "right")

# Nothing significant for flowering time

#-----------------------------------------------------------------------------
# Fit the multivariate model for trait and not PC
model_biomass <- brm(
  formula = Biomass ~ Treatment * allele * sub_pop + (1 | Block) + (1 | Samples) + (1 | panicle),
  data = pca_data,
  family = gaussian()
)

# Summarize the model
summary(model_biomass)

# Extract fixed effects (estimates, standard errors, and confidence intervals)
fixed_effects <- fixef(model_biomass, summary = TRUE)

# Convert to a dataframe for ggplot2
fixed_effects_df <- as.data.frame(fixed_effects)
fixed_effects_df$Term <- rownames(fixed_effects_df)

# Rename columns for clarity
colnames(fixed_effects_df) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Term")

# Filter and separate the Term column using the correct separator
fixed_effects_df_1 <- fixed_effects_df %>%
  filter((`l-95% CI` > 0 & `u-95% CI` > 0) | (`l-95% CI` < 0 & `u-95% CI` < 0)) #%>%
#separate(Term, into = c("Response", "Interaction"), sep = "_", fill = "right")

write.table(fixed_effects_df_1, file = "~/Documents/Rice_Research/biomass_response_to_heat_chr3.txt", sep = "\t", quote = F)

#----------------------------------------------------------------------------
# Fit the multivariate model for trait and not PC
model_spikelet_no <- brm(
  formula = Spikelet_no ~ Treatment * allele * sub_pop + (1 | Block) + (1 | Samples) + (1 | panicle),
  data = pca_data,
  family = gaussian()
)

# Summarize the model
summary(model_spikelet_no)

# Extract fixed effects (estimates, standard errors, and confidence intervals)
fixed_effects <- fixef(model_spikelet_no, summary = TRUE)

# Convert to a dataframe for ggplot2
fixed_effects_df <- as.data.frame(fixed_effects)
fixed_effects_df$Term <- rownames(fixed_effects_df)

# Rename columns for clarity
colnames(fixed_effects_df) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Term")

# Filter and separate the Term column using the correct separator
fixed_effects_df_1 <- fixed_effects_df %>%
  filter((`l-95% CI` > 0 & `u-95% CI` > 0) | (`l-95% CI` < 0 & `u-95% CI` < 0)) #%>%
#separate(Term, into = c("Response", "Interaction"), sep = "_", fill = "right")

# Nothing significant for spikelet number

#------------------------------------------------------------------------
# Fit the multivariate model for trait and not PC
model_seed_per_panicle <- brm(
  formula = Seed_per_panicle ~ Treatment * allele * sub_pop + (1 | Block) + (1 | Samples) + (1 | panicle),
  data = pca_data,
  family = gaussian()
)

# Summarize the model
summary(model_seed_per_panicle)

# Extract fixed effects (estimates, standard errors, and confidence intervals)
fixed_effects <- fixef(model_seed_per_panicle, summary = TRUE)

# Convert to a dataframe for ggplot2
fixed_effects_df <- as.data.frame(fixed_effects)
fixed_effects_df$Term <- rownames(fixed_effects_df)

# Rename columns for clarity
colnames(fixed_effects_df) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Term")

# Filter and separate the Term column using the correct separator
fixed_effects_df_1 <- fixed_effects_df %>%
  filter((`l-95% CI` > 0 & `u-95% CI` > 0) | (`l-95% CI` < 0 & `u-95% CI` < 0)) #%>%
#separate(Term, into = c("Response", "Interaction"), sep = "_", fill = "right")

write.table(fixed_effects_df_1, file = "~/Documents/Rice_Research/seed_per_panicle_response_to_heat_chr3.txt", sep = "\t", quote = F)

#------------------------------------------------------------------------
# Fit the multivariate model for trait and not PC
model_spikelet_fertility <- brm(
  formula = Spikelet_fertility ~ Treatment * allele * sub_pop + (1 | Block) + (1 | Samples) + (1 | panicle),
  data = pca_data,
  family = gaussian()
)

# Summarize the model
summary(model_spikelet_fertility)

# Extract fixed effects (estimates, standard errors, and confidence intervals)
fixed_effects <- fixef(model_spikelet_fertility, summary = TRUE)

# Convert to a dataframe for ggplot2
fixed_effects_df <- as.data.frame(fixed_effects)
fixed_effects_df$Term <- rownames(fixed_effects_df)

# Rename columns for clarity
colnames(fixed_effects_df) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Term")

# Filter and separate the Term column using the correct separator
fixed_effects_df_1 <- fixed_effects_df %>%
  filter((`l-95% CI` > 0 & `u-95% CI` > 0) | (`l-95% CI` < 0 & `u-95% CI` < 0)) #%>%
#separate(Term, into = c("Response", "Interaction"), sep = "_", fill = "right")

# nothing significant for spikelet fertility

#-------------------------------------------------------------------------
# Fit the multivariate model for trait and not PC
model_seed_weight <- brm(
  formula = Seed_weight ~ Treatment * allele * sub_pop + (1 | Block) + (1 | Samples) + (1 | panicle),
  data = pca_data,
  family = gaussian()
)

# Summarize the model
summary(model_seed_weight)

# Extract fixed effects (estimates, standard errors, and confidence intervals)
fixed_effects <- fixef(model_seed_weight, summary = TRUE)

# Convert to a dataframe for ggplot2
fixed_effects_df <- as.data.frame(fixed_effects)
fixed_effects_df$Term <- rownames(fixed_effects_df)

# Rename columns for clarity
colnames(fixed_effects_df) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Term")

# Filter and separate the Term column using the correct separator
fixed_effects_df_1 <- fixed_effects_df %>%
  filter((`l-95% CI` > 0 & `u-95% CI` > 0) | (`l-95% CI` < 0 & `u-95% CI` < 0)) #%>%
#separate(Term, into = c("Response", "Interaction"), sep = "_", fill = "right")

write.table(fixed_effects_df_1, file = "~/Documents/Rice_Research/seed_weight_response_to_heat_chr3.txt", sep = "\t", quote = F)




#----------------------------------------------------------------------------
# mixed model regression with PCs
# Select relevant columns for PCA for chr7 for AA allele
chr7_data_AA <- seed_long %>%
  select(Biomass, flowering_time, Spikelet_fertility, Seed_per_panicle, Spikelet_no, Seed_weight, Treatment, chromosome, chr3, chr7, sub_pop, Samples, Block, panicle) %>%
  filter(Treatment %in% c("Control", "Heat")) %>%
  filter(chromosome %in% "chr7") %>%
  filter(chr7 %in% "AA")

pca_data <- chr7_data_AA %>%
  select(., Biomass, flowering_time, Spikelet_fertility, Seed_per_panicle, Spikelet_no, Seed_weight)

traits_standardized <- scale(pca_data)

# Perform PCA with missing values
pca_result <- PCA(traits_standardized, ncp = 5, scale.unit = TRUE, graph = FALSE)

# View results
summary(pca_result)

# Extract the first 3 PCs from the individual coordinates
pca_scores <- pca_result$ind$coord[, 1:3]

# Preview the extracted scores
head(pca_scores)

# Create a dataframe with PC scores
pca_data <- as.data.frame(pca_scores)
pca_data$Treatment <- chr7_data_AA$Treatment
pca_data$sub_pop <- chr7_data_AA$sub_pop
pca_data$Samples <- chr7_data_AA$Samples
pca_data$allele <- chr7_data_AA$chr7
pca_data$panicle <- chr7_data_AA$panicle
pca_data$Block <- chr7_data_AA$Block

pca_data$allele <- as.character(pca_data$allele)

pca_data <- cbind(pca_data, chr7_data_AA[, 1:6]) 


# Add the treatment data to pca_result
pca_result$ind$Treatment <- pca_data$Treatment
pca_result$ind$allele <- pca_data$allele


# Now create the plot
pc <- fviz_pca_biplot(pca_result, 
                      axes = c(1, 2),  # Plot PC1 and PC2
                      col.var = "red",  # Color arrows
                      col.ind = pca_data$Treatment,
                      geom = c("point"),  # Plot only points, no text labels
                      pointsize = 3,  # Increase point size
                      repel = TRUE,  # Avoid text overlap
                      addEllipses = TRUE,  # Add confidence ellipses
                      legend.title = "Treatment") +  # Map color to Treatment
  theme_minimal() +
  labs(
    title = "AA",
    x = paste0("PC1 (", round(pca_result$eig[1, 2], 1), "%)"),
    y = paste0("PC2 (", round(pca_result$eig[2, 2], 1), "%)")
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Title customization
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),   # Axis text size
    legend.position = "right"              # Position the legend on the right
  ) +
  scale_color_manual(values = c("Control" = "darkgreen", "Heat" = "tomato3")) #+
#facet_wrap(~ pca_data$allele)  # Facet by allele (GG and AA)# Custom color scale

#-------------------------------------------------------------------------------
# Select relevant columns for PCA for chr3 for GG allele
chr7_data_GG <- seed_long %>%
  select(Biomass, flowering_time, Spikelet_fertility, Seed_per_panicle, Spikelet_no, Seed_weight, Treatment, chromosome, chr3, chr7, sub_pop, Samples, Block, panicle) %>%
  filter(Treatment %in% c("Control", "Heat")) %>%
  filter(chromosome %in% "chr7") %>%
  filter(chr7 %in% "GG")


pca_data <- chr7_data_GG %>%
  select(., Biomass, flowering_time, Spikelet_fertility, Seed_per_panicle, Spikelet_no, Seed_weight)

traits_standardized <- scale(pca_data)

# Perform PCA with missing values
pca_result <- PCA(traits_standardized, ncp = 5, scale.unit = TRUE, graph = FALSE)

# View results
summary(pca_result)

# Extract the first 3 PCs from the individual coordinates
pca_scores <- pca_result$ind$coord[, 1:3]

# Preview the extracted scores
head(pca_scores)

# Create a dataframe with PC scores
pca_data <- as.data.frame(pca_scores)
pca_data$Treatment <- chr7_data_GG$Treatment
pca_data$sub_pop <- chr7_data_GG$sub_pop
pca_data$Samples <- chr7_data_GG$Samples
pca_data$allele <- chr7_data_GG$chr3
pca_data$panicle <- chr7_data_GG$panicle
pca_data$Block <- chr7_data_GG$Block

pca_data$allele <- as.character(pca_data$allele)

pca_data <- cbind(pca_data, chr7_data_GG[, 1:6]) 


# Add the treatment data to pca_result
pca_result$ind$Treatment <- pca_data$Treatment
pca_result$ind$allele <- pca_data$allele


# Now create the plot
pc <- fviz_pca_biplot(pca_result, 
                      axes = c(1, 2),  # Plot PC1 and PC2
                      col.var = "red",  # Color arrows
                      col.ind = pca_data$Treatment,
                      geom = c("point"),  # Plot only points, no text labels
                      pointsize = 3,  # Increase point size
                      repel = TRUE,  # Avoid text overlap
                      addEllipses = TRUE,  # Add confidence ellipses
                      legend.title = "Treatment") +  # Map color to Treatment
  theme_minimal() +
  labs(
    title = "GG",
    x = paste0("PC1 (", round(pca_result$eig[1, 2], 1), "%)"),
    y = paste0("PC2 (", round(pca_result$eig[2, 2], 1), "%)")
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Title customization
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),   # Axis text size
    legend.position = "right"              # Position the legend on the right
  ) +
  scale_color_manual(values = c("Control" = "darkgreen", "Heat" = "tomato3")) #+
#facet_wrap(~ pca_data$allele)  # Facet by allele (GG and AA)# Custom color scale


#------------------------------------------------------------------------
chr7_data <- seed_long %>%
  select(Biomass, flowering_time, Spikelet_fertility, Seed_per_panicle, Spikelet_no, Seed_weight, Treatment, chromosome, chr3, chr7, sub_pop, Samples, Block, panicle) %>%
  filter(Treatment %in% c("Control", "Heat")) %>%
  filter(chromosome %in% "chr7") 

pca_data <- chr7_data %>%
  select(., Biomass, flowering_time, Spikelet_fertility, Seed_per_panicle, Spikelet_no, Seed_weight)

traits_standardized <- scale(pca_data)

# Perform PCA with missing values
pca_result <- PCA(traits_standardized, ncp = 5, scale.unit = TRUE, graph = FALSE)

# View results
summary(pca_result)

# Extract the first 3 PCs from the individual coordinates
pca_scores <- pca_result$ind$coord[, 1:3]

# Preview the extracted scores
head(pca_scores)

# Create a dataframe with PC scores
pca_data <- as.data.frame(pca_scores)
pca_data$Treatment <- chr7_data$Treatment
pca_data$sub_pop <- chr7_data$sub_pop
pca_data$Samples <- chr7_data$Samples
pca_data$allele <- chr7_data$chr7
pca_data$panicle <- chr7_data$panicle
pca_data$Block <- chr7_data$Block

pca_data <- cbind(pca_data, chr7_data[, 1:6]) 

model_brms <- brm(
  formula = mvbind(Dim.1, Dim.2, Dim.3) ~ Treatment * allele + (1 | Block) + (1 | Samples) + (1 | panicle),
  data = pca_data,
  family = gaussian()
)

# Summarize the model
summary(model_brms)


# Extract fixed effects (estimates, standard errors, and confidence intervals)
fixed_effects <- fixef(model_brms, summary = TRUE)

# Convert to a dataframe for ggplot2
fixed_effects_df <- as.data.frame(fixed_effects)
fixed_effects_df$Term <- rownames(fixed_effects_df)

# Rename columns for clarity
colnames(fixed_effects_df) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Term")

# Filter and separate the Term column using the correct separator
fixed_effects_df <- fixed_effects_df %>%
  filter((`l-95% CI` > 0 & `u-95% CI` > 0) | (`l-95% CI` < 0 & `u-95% CI` < 0)) #%>%
#separate(Term, into = c("Response", "Interaction"), sep = "_", fill = "right")

# Remove "PC_" from the Term column
fixed_effects_df$Term <- gsub("PC\\d_", "", fixed_effects_df$Term)


# Plot fixed effects with confidence intervals
p1 <- ggplot(fixed_effects_df, aes(x = Estimate, y = Term)) +
  geom_point(colour = "blue", size = 3) +  # Change point color to blue and increase point size
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`), height = 0.2, colour = "blue", size = 1.5) +  # Change error bar color to blue and increase line thickness
  theme_minimal() +
  #facet_grid(.~ Response, scales = "free") + 
  xlab("Effect Estimate") + 
  ylab("Interaction") +
  ggtitle("Fixed Effects with 95% Confidence Intervals for Chr3") +
  theme(
    axis.text.x = element_text(size = 15),  # Increase x-axis text size
    axis.text.y = element_text(size = 15),  # Increase y-axis text size
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16)   # Increase y-axis title size
  )

#-----------------------------------------------------------------------
# Fit the multivariate model for traits only
model_biomass <- brm(
  formula = Biomass ~ Treatment * allele + (1 | Block) + (1 | Samples) + (1 | panicle),
  data = pca_data,
  family = gaussian()
)

# Summarize the model
summary(model_biomass)

# Extract fixed effects (estimates, standard errors, and confidence intervals)
fixed_effects <- fixef(model_biomass, summary = TRUE)

# Convert to a dataframe for ggplot2
fixed_effects_df <- as.data.frame(fixed_effects)
fixed_effects_df$Term <- rownames(fixed_effects_df)

# Rename columns for clarity
colnames(fixed_effects_df) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Term")

# Filter and separate the Term column using the correct separator
fixed_effects_df_1 <- fixed_effects_df %>%
  filter((`l-95% CI` > 0 & `u-95% CI` > 0) | (`l-95% CI` < 0 & `u-95% CI` < 0)) #%>%
#separate(Term, into = c("Response", "Interaction"), sep = "_", fill = "right")

write.table(fixed_effects_df_1, file = "~/Documents/Rice_Research/biomass_response_to_heat_chr7.txt", sep = "\t", quote = F)

#----------------------------------------------------------------------------
# Fit the multivariate model for trait and not PC
model_spikelet_no <- brm(
  formula = Spikelet_no ~ Treatment * allele + (1 | Block) + (1 | Samples) + (1 | panicle),
  data = pca_data,
  family = gaussian()
)

# Summarize the model
summary(model_spikelet_no)

# Extract fixed effects (estimates, standard errors, and confidence intervals)
fixed_effects <- fixef(model_spikelet_no, summary = TRUE)

# Convert to a dataframe for ggplot2
fixed_effects_df <- as.data.frame(fixed_effects)
fixed_effects_df$Term <- rownames(fixed_effects_df)

# Rename columns for clarity
colnames(fixed_effects_df) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Term")

# Filter and separate the Term column using the correct separator
fixed_effects_df_1 <- fixed_effects_df %>%
  filter((`l-95% CI` > 0 & `u-95% CI` > 0) | (`l-95% CI` < 0 & `u-95% CI` < 0)) #%>%
#separate(Term, into = c("Response", "Interaction"), sep = "_", fill = "right")

# Nothing significant for spikelet number

#------------------------------------------------------------------------
# Fit the multivariate model for trait and not PC
model_seed_per_panicle <- brm(
  formula = Seed_per_panicle ~ Treatment * allele + (1 | Block) + (1 | Samples) + (1 | panicle),
  data = pca_data,
  family = gaussian()
)

# Summarize the model
summary(model_seed_per_panicle)

# Extract fixed effects (estimates, standard errors, and confidence intervals)
fixed_effects <- fixef(model_seed_per_panicle, summary = TRUE)

# Convert to a dataframe for ggplot2
fixed_effects_df <- as.data.frame(fixed_effects)
fixed_effects_df$Term <- rownames(fixed_effects_df)

# Rename columns for clarity
colnames(fixed_effects_df) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Term")

# Filter and separate the Term column using the correct separator
fixed_effects_df_1 <- fixed_effects_df %>%
  filter((`l-95% CI` > 0 & `u-95% CI` > 0) | (`l-95% CI` < 0 & `u-95% CI` < 0)) #%>%
#separate(Term, into = c("Response", "Interaction"), sep = "_", fill = "right")

# no response of seed_per_panicle to heat stress

#------------------------------------------------------------------------
# Fit the multivariate model for trait and not PC
model_spikelet_fertility <- brm(
  formula = Spikelet_fertility ~ Treatment * allele + (1 | Block) + (1 | Samples) + (1 | panicle),
  data = pca_data,
  family = gaussian()
)

# Summarize the model
summary(model_spikelet_fertility)

# Extract fixed effects (estimates, standard errors, and confidence intervals)
fixed_effects <- fixef(model_spikelet_fertility, summary = TRUE)

# Convert to a dataframe for ggplot2
fixed_effects_df <- as.data.frame(fixed_effects)
fixed_effects_df$Term <- rownames(fixed_effects_df)

# Rename columns for clarity
colnames(fixed_effects_df) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Term")

# Filter and separate the Term column using the correct separator
fixed_effects_df_1 <- fixed_effects_df %>%
  filter((`l-95% CI` > 0 & `u-95% CI` > 0) | (`l-95% CI` < 0 & `u-95% CI` < 0)) #%>%
#separate(Term, into = c("Response", "Interaction"), sep = "_", fill = "right")

# nothing significant for spikelet fertility

#-------------------------------------------------------------------------
# Fit the multivariate model for trait and not PC
model_seed_weight <- brm(
  formula = Seed_weight ~ Treatment * allele + (1 | Block) + (1 | Samples) + (1 | panicle),
  data = pca_data,
  family = gaussian()
)

# Summarize the model
summary(model_seed_weight)

# Extract fixed effects (estimates, standard errors, and confidence intervals)
fixed_effects <- fixef(model_seed_weight, summary = TRUE)

# Convert to a dataframe for ggplot2
fixed_effects_df <- as.data.frame(fixed_effects)
fixed_effects_df$Term <- rownames(fixed_effects_df)

# Rename columns for clarity
colnames(fixed_effects_df) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Term")

# Filter and separate the Term column using the correct separator
fixed_effects_df_1 <- fixed_effects_df %>%
  filter((`l-95% CI` > 0 & `u-95% CI` > 0) | (`l-95% CI` < 0 & `u-95% CI` < 0)) #%>%
#separate(Term, into = c("Response", "Interaction"), sep = "_", fill = "right")

# Nothing significant for weight

#----------------------------------------------------------------------
# Fit the multivariate model for trait and not PC
model_flowering_time <- brm(
  formula = flowering_time ~ Treatment * allele + (1 | Block) + (1 | Samples) + (1 | panicle),
  data = pca_data,
  family = gaussian()
)

# Summarize the model
summary(model_flowering_time)

# Extract fixed effects (estimates, standard errors, and confidence intervals)
fixed_effects <- fixef(model_flowering_time, summary = TRUE)

# Convert to a dataframe for ggplot2
fixed_effects_df <- as.data.frame(fixed_effects)
fixed_effects_df$Term <- rownames(fixed_effects_df)

# Rename columns for clarity
colnames(fixed_effects_df) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Term")

# Filter and separate the Term column using the correct separator
fixed_effects_df_1 <- fixed_effects_df %>%
  filter((`l-95% CI` > 0 & `u-95% CI` > 0) | (`l-95% CI` < 0 & `u-95% CI` < 0)) #%>%
#separate(Term, into = c("Response", "Interaction"), sep = "_", fill = "right")

write.table(fixed_effects_df_1, file = "~/Documents/Rice_Research/flowering_time_response_to_heat_chr7.txt", sep = "\t", quote = F)

#-----------------------------------------------------------------------
# Plot interaction plots for all significant responses for Chr3
chr3_response <- read.table("~/Documents/Rice_Research/Chr3_reproductive_traits_response_heat.txt", skip = 1)
colnames(chr3_response) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Term", "Trait")


# Plot fixed effects with confidence intervals
p1 <- ggplot(chr3_response, aes(x = Estimate, y = Term, color = Trait)) +
  geom_point(size = 3) +  # Color points based on Trait and set size
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`, color = Trait), 
                 height = 0.2, size = 1.5) +  # Color error bars based on Trait
  theme_minimal() +
  facet_grid(. ~ Trait, scales = "free") + 
  xlab("Effect Estimate") + 
  ylab("Interaction") +
  ggtitle("Effect Estimates by Trait on Chr3") +
  theme(
    axis.text.x = element_text(size = 15),  # Increase x-axis text size
    axis.text.y = element_text(size = 15),  # Increase y-axis text size
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16)   # Increase y-axis title size
  )


# Plot interaction plots for all significant responses for Chr7
chr7_response <- read.table("~/Documents/Rice_Research/chr7_reproductive_trait_response_heat.txt", 
                            sep = "\t", header = TRUE, fill = TRUE)

colnames(chr7_response) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Term", "Trait")

# Plot fixed effects with confidence intervals
p1 <- ggplot(chr7_response, aes(x = Estimate, y = Term, color = Trait)) +
  geom_point(size = 3) +  # Color points based on Trait and set size
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`, color = Trait), 
                 height = 0.2, size = 1.5) +  # Color error bars based on Trait
  theme_minimal() +
  facet_grid(. ~ Trait, scales = "free") + 
  xlab("Effect Estimate") + 
  ylab("Interaction") +
  ggtitle("Effect Estimates by Trait on Chr7") +
  theme(
    axis.text.x = element_text(size = 15),  # Increase x-axis text size
    axis.text.y = element_text(size = 15),  # Increase y-axis text size
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16)   # Increase y-axis title size
  )



# Create a summary of the means for multiple traits
interaction_data <- chr3_data %>%
  pivot_longer(
    cols = c(flowering_time, Spikelet_no, Spikelet_fertility, Seed_per_panicle, Biomass, Seed_weight),
    names_to = "trait", 
    values_to = "value"
  ) %>%
  group_by(Treatment, chr3, trait, sub_pop) %>%  # Keep sub_pop as a grouping variable
  summarize(mean_value = mean(value, na.rm = TRUE), .groups = 'drop')


# Create the interaction plot with faceting for each trait
p7 <- ggplot(interaction_data, aes(x = Treatment, y = mean_value, color = chr3, group = chr3)) +
  geom_line(size = 1) +            # Line for each group
  geom_point(size = 3) +           # Points for mean values
  labs(title = "Interaction Plot of Treatment and Chr3 Allele on Traits",
       x = "Treatment",
       y = "Mean Value") +
  theme_minimal() +                # Clean theme
  scale_color_manual(values = c("red", "blue")) + # Custom colors
  theme(legend.title = element_blank()) +  # Remove legend title
  facet_wrap(trait ~ sub_pop, scales = "free_y")  # Facet by trait and sub_pop

# Print the plot
print(p7)

#------------------------------------------------------------------------
# Select relevant columns for PCA for chr7
chr7_data <- seed_long %>%
  select(Biomass, flowering_time, Spikelet_fertility, Seed_per_panicle, Spikelet_no, Treatment, chromosome, chr3, chr7, sub_pop, Samples, Block, panicle) %>%
  filter(Treatment %in% c("Control", "Heat")) %>%
  filter(chromosome %in% "chr7") 

pca_data <- chr7_data %>%
  select(., Biomass, flowering_time, Spikelet_fertility, Seed_per_panicle, Spikelet_no)


traits_standardized <- scale(pca_data)

# Perform PCA with missing values
pca_result <- PCA(traits_standardized, ncp = 5, scale.unit = TRUE, graph = FALSE)

# View results
summary(pca_result)

# Extract the first 3 PCs from the individual coordinates
pca_scores <- pca_result$ind$coord[, 1:3]

# Preview the extracted scores
head(pca_scores)

# Create a dataframe with PC scores
pca_data <- as.data.frame(pca_scores)
pca_data$Treatment <- chr7_data$Treatment
pca_data$sub_pop <- chr7_data$sub_pop
pca_data$Samples <- chr7_data$Samples
pca_data$allele <- chr7_data$chr7
pca_data$panicle <- chr7_data$panicle
pca_data$Block <- chr7_data$Block

# Extract the explained variance (percentage) for PC1 and PC2
explained_variance <- pca_result$eig[1:2, 2]  # Extract the percentages for PC1 and PC2

# Print the explained variance percentages
explained_variance


# Plot PCA with custom points and larger labels/titles
ggplot(pca_data, aes(x = Dim.1, y = Dim.2, color = Treatment)) +
  geom_point(size = 4) +  # Increase point size
  scale_colour_manual(values = c("Control" = "darkgreen", "Heat" = "tomato3")) +  # Custom color scale
  facet_wrap(~ allele) +  # Facet by allele
  theme_minimal() +
  labs(
    title = "PCA of Reproductive Traits (Chromosome 7)",
    x = paste0("PC1 (", round(explained_variance[1], 1), "%)"),
    y = paste0("PC2 (", round(explained_variance[2], 1), "%)")
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Increase title size and make it bold
    axis.title = element_text(size = 14),  # Increase axis title size
    axis.text = element_text(size = 12)   # Increase axis label size
  )

# Fit the multivariate model
model_brms <- brm(
  formula = mvbind(Dim.1, Dim.2, Dim.3) ~ Treatment * allele + (1 | Block) + (1 | Samples) + (1 | panicle),
  data = pca_data,
  family = gaussian()
)

# Summarize the model
summary(model_brms)

# Extract fixed effects (estimates, standard errors, and confidence intervals)
fixed_effects <- fixef(model_brms, summary = TRUE)

# Convert to a dataframe for ggplot2
fixed_effects_df <- as.data.frame(fixed_effects)
fixed_effects_df$Term <- rownames(fixed_effects_df)

# Rename columns for clarity
colnames(fixed_effects_df) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Term")

# Filter and separate the Term column using the correct separator
fixed_effects_df <- fixed_effects_df %>%
  filter((`l-95% CI` > 0 & `u-95% CI` > 0) | (`l-95% CI` < 0 & `u-95% CI` < 0)) %>%
  separate(Term, into = c("Response", "Interaction"), sep = "_", fill = "right")

# Plot fixed effects with confidence intervals
p1 <- ggplot(fixed_effects_df, aes(x = Estimate, y = Interaction)) +
  geom_point(colour = "blue", size = 3) +  # Change point color to blue and increase point size
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`), height = 0.2, colour = "blue", size = 1.5) +  # Change error bar color to blue and increase line thickness
  theme_minimal() +
  facet_wrap(~ Response, scales = "free") + 
  xlab("Effect Estimate") + 
  ylab("Interaction") +
  ggtitle("Fixed Effects with 95% Confidence Intervals for Chr7") +
  theme(
    axis.text.x = element_text(size = 15),  # Increase x-axis text size
    axis.text.y = element_text(size = 15),  # Increase y-axis text size
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16)   # Increase y-axis title size
  )


# Create a summary of the means for multiple traits
interaction_data <- chr7_data %>%
  pivot_longer(cols = c(flowering_time, Spikelet_no, Spikelet_fertility, Seed_per_panicle, Biomass),  # Add other traits here
               names_to = "trait", 
               values_to = "value") %>%
  group_by(Treatment, chr7, trait) %>%
  summarize(mean_value = mean(value, na.rm = TRUE), .groups = 'drop')

# Create the interaction plot with faceting for each trait
p7 <- ggplot(interaction_data, aes(x = Treatment, y = mean_value, color = chr7, group = chr7)) +
  geom_line(size = 1) +            # Line for each group
  geom_point(size = 3) +           # Points for mean values
  labs(title = "Interaction Plot of Treatment and Chr7 Allele on Traits",
       x = "Treatment",
       y = "Mean Value") +
  theme_minimal() +                 # Clean theme
  scale_color_manual(values = c("red", "blue")) +  # Custom colors
  theme(legend.title = element_blank()) +  # Remove legend title
  facet_wrap(~ trait, scales = "free_y")  # Facet by trait

# Print the plot
print(p7)

#------------------------------------------------------------------------
# Select relevant columns for PCA for Kitaake
kx_data <- seed_long %>%
  select(Biomass, flowering_time, Spikelet_fertility, Seed_per_panicle, Spikelet_no, Treatment, sub_pop, Samples, Block, panicle) %>%
  filter(Treatment %in% c("Control", "Heat")) %>%
  filter(sub_pop %in% c("Kitaake"))

kx_data$Samples <- factor(kx_data$Samples)
kx_data$Samples <- relevel(kx_data$Samples, ref = "WT")

pca_data <- kx_data %>%
  select(., Biomass, flowering_time, Spikelet_fertility, Seed_per_panicle, Spikelet_no)


traits_standardized <- scale(pca_data)

# Perform PCA with missing values
pca_result <- PCA(traits_standardized, ncp = 5, scale.unit = TRUE, graph = FALSE)

# View results
summary(pca_result)

# Extract the first 3 PCs from the individual coordinates
pca_scores <- pca_result$ind$coord[, 1:3]

# Preview the extracted scores
head(pca_scores)

# Create a dataframe with PC scores
pca_data <- as.data.frame(pca_scores)
pca_data$Treatment <- kx_data$Treatment
pca_data$sub_pop <- kx_data$sub_pop
pca_data$Samples <- kx_data$Samples
pca_data$panicle <- kx_data$panicle
pca_data$Block <- kx_data$Block

# Extract the explained variance (percentage) for PC1 and PC2
explained_variance <- pca_result$eig[1:2, 2]  # Extract the percentages for PC1 and PC2

# Print the explained variance percentages
explained_variance


# Plot PCA with custom points and larger labels/titles
ggplot(pca_data, aes(x = Dim.1, y = Dim.2, color = Treatment)) +
  geom_point(size = 4) +  # Increase point size
  scale_colour_manual(values = c("Control" = "darkgreen", "Heat" = "tomato3")) +  # Custom color scale
  facet_wrap(~ Samples) +  # Facet by allele
  theme_minimal() +
  labs(
    title = "PCA of Reproductive Traits (Chromosome 3)",
    x = paste0("PC1 (", round(explained_variance[1], 1), "%)"),
    y = paste0("PC2 (", round(explained_variance[2], 1), "%)")
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Increase title size and make it bold
    axis.title = element_text(size = 14),  # Increase axis title size
    axis.text = element_text(size = 12)   # Increase axis label size
  )


# Fit the multivariate model
model_brms <- brm(
  formula = mvbind(Dim.1, Dim.2, Dim.3) ~ Treatment * Samples + (1 | Block) + (1 | panicle),
  data = pca_data,
  family = gaussian()
)

# Summarize the model
summary(model_brms)

# Extract fixed effects (estimates, standard errors, and confidence intervals)
fixed_effects <- fixef(model_brms, summary = TRUE)

# Convert to a dataframe for ggplot2
fixed_effects_df <- as.data.frame(fixed_effects)
fixed_effects_df$Term <- rownames(fixed_effects_df)

# Rename columns for clarity
colnames(fixed_effects_df) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Term")

# Filter and separate the Term column using the correct separator
fixed_effects_df <- fixed_effects_df %>%
  filter((`l-95% CI` > 0 & `u-95% CI` > 0) | (`l-95% CI` < 0 & `u-95% CI` < 0)) %>%
  separate(Term, into = c("Response", "Interaction"), sep = "_", fill = "right")

# Plot fixed effects with confidence intervals
p1 <- ggplot(fixed_effects_df, aes(x = Estimate, y = Interaction)) +
  geom_point(colour = "cyan4", size = 3) +  # Change point color to blue and increase point size
  geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u-95% CI`), height = 0.2, colour = "cyan4", size = 1.5) +  # Change error bar color to blue and increase line thickness
  theme_minimal() +
  facet_wrap(~ Response, scales = "free") + 
  xlab("Effect Estimate") + 
  ylab("Interaction") +
  ggtitle("Fixed Effects with 95% Confidence Intervals for Kitaake Mutants") +
  theme(
    axis.text.x = element_text(size = 15),  # Increase x-axis text size
    axis.text.y = element_text(size = 15),  # Increase y-axis text size
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16)   # Increase y-axis title size
  )


# Create a summary of the means for multiple traits
interaction_data <- kx_data %>%
  pivot_longer(cols = c(flowering_time, Spikelet_no, Spikelet_fertility, Seed_per_panicle, Biomass),  # Add other traits here
               names_to = "trait", 
               values_to = "value") %>%
  group_by(Treatment, Samples, trait) %>%
  summarize(mean_value = mean(value, na.rm = TRUE), .groups = 'drop')

# Create the interaction plot with faceting for each trait
p7 <- ggplot(interaction_data, aes(x = Treatment, y = mean_value, color = Samples, group = Samples)) +
  geom_line(size = 1) +            # Line for each group
  geom_point(size = 3) +           # Points for mean values
  labs(title = "Interaction Plot of Treatment and Kitaake Samples on Traits",
       x = "Treatment",
       y = "Mean Value") +
  theme_minimal() +                 # Clean theme
  scale_color_manual(values = c("violet", "cyan4", "yellow4")) +  # Custom colors
  theme(legend.title = element_blank()) +  # Remove legend title
  facet_wrap(~ trait, scales = "free_y")  # Facet by trait

# Print the plot
print(p7)

#------------------------------------------------------------------------

seed_long %>%
  filter(Treatment %in% c("Control", "Heat")) %>%
  filter(chromosome == "chr3") %>%
  #filter(sub_pop == "TRJ") %>%
  ggplot(aes(x = Treatment, y = Spikelet_no, fill = Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = c("Control" = "palegreen1", "Heat" = "tomato3")) +  # Custom fill colors
  stat_compare_means(aes(group = Treatment), label = "p.signif") +
  labs(title = "Spikelet Number by Treatment (Chromosome 3)",
       x = "Treatment",
       y = "Spikelet Number") +
  theme_minimal() +
  facet_wrap(~sub_pop)  # Facet by two variables
