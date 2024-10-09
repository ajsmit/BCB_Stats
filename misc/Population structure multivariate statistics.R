# Population structure multivariate analysis
# 01 October 2024



# Load the data and packages ----------------------------------------------
library(tidyverse)
library(vegan)
library(viridis)

pop <- read_csv("misc/population_structure_final_2024.csv") # loaded the data for the pcoa

pop$year <- as.factor(pop$year) # made year as a factor
pop$species <- as.factor(pop$species)

pop1 <- pop %>%
  select(-estuary, -species, -year) # removed the non-numeric columns to standardize and do the pcoa
pop1 <- sqrt(pop1) # square root transform the data
pop1_std <- decostand(pop1, method = "standardize") # did the standardization

pop_pcoa <- capscale(pop1_std ~ 1)
summary(pop_pcoa)

col <- c("black", "firebrick2")
pch <- c(0, 1, 2)

opar <- par(no.readonly = TRUE)
ordiplot(pop_pcoa, type = "n", scaling = 2,
         main = "PCoA of mangrove stuff")
points(pop_pcoa, "sites", pch = pch[pop$species], cex = 1.75, col = col[pop$year], bg = "white")
text(pop_pcoa, "sites", col = col[pop$year], cex = 0.5)

# attempt at adding a custom legend corresponding to the colours and shapes
leg.txt <- c("AM", "BG", "RM", "2018", "2019")
legend("bottomleft",leg.txt, col = col[pop1$year], pch = pch[pop1$species])


# Alternative: extract sores for plotting via ggplot


data.scores = as.data.frame(scores(pop_pcoa)$sites) # extracting scores and making it a data frame

data.scores$year = pop1$year # adding Year and site back to the data frame
data.scores$species = pop1$species

ggplot(data.scores, aes(x = MDS1, y = MDS2)) +
  geom_point(size = 4, aes( shape = species, colour = year))+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 12),
        legend.text = element_text(size = 12, face ="bold", colour ="black"),
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.background = element_blank(),panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank()) +
  labs(x = "MDS1", colour = "year", y = "MDS2", shape = "species") +
  scale_colour_manual(values = c("firebrick1", "gray27"))




# PCA -------------------------------------------------

pop_pca <- rda(pop1_std)
summary(pop_pca)

scores = as.data.frame(scores(pop_pca)$species)

scores$year = pop$year # adding Year and site back to the dataframe
scores$species = pop$species

ggplot(scores, aes(x = PC1, y = PC2)) +
  geom_point(size = 4, aes(shape = species, colour = year))+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 12),
        legend.text = element_text(size = 12, face ="bold", colour ="black"),
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.background = element_blank(),panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank()) +
  labs(x = "PC1", colour = "species", y = "PC2", shape = "yeaer") +
  # scale_colour_manual(values = c("firebrick1", "gray27", "royalblue"))
  scale_colour_manual(values = c("firebrick1", "gray27"))




# GLM's -------------------------------------------------------------------


pop <- read_csv("Data/Pop_struct_final.csv")

# Questions to ask

# Within an estuary, do species differ in terms of the response variables (8 columns)
# For only BW:
# glm(adult_dens ~ species) is adult density different across species in beachwood estuary, repeat for all 3 estuaries


# Does AM response differ as a function of estuary
 #filter out only AM
# glm (adult_dens ~ estuary), repeat for next species

# Can't do interactions between things, does adult density depend upon estuary and is it modif
# glm(adult_dens ~ species + estuary)


# Hypothesis 1 ------------------------------------------------------------
# Within an estuary, do species differ in terms of the response variables?



# Beachwood estuary -------------------------------------------------------

# Adults

pop_bw <- pop %>%
  filter(estuary == "BW") # filtered out Beachwood estuary

adult_bw <- glm(adult_dens ~ estuary * species, family = poisson, data = pop) # run the glm and make it an object to summarize

summary(adult_bw) # summary output of the model

# visualizing the results
ggplot(data = pop_bw, aes( y = adult_dens, fill = species)) +
  geom_boxplot() +
  labs(y = "Adult density", x = "", title = " Beachwood ") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# Saplings

sap_bw <- glm(sap_dens ~ species, family = poisson, data = pop_bw) # run the glm and make it an object to summarize

summary(sap_bw) # summary output of the model

# visualizing the results
ggplot(data = pop_bw, aes( y = sap_dens, fill = species)) +
  geom_boxplot() +
  labs(y = "Sapling density", x = "", title = " Beachwood ") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# Seedlings

seed_bw <- glm(seed_dens ~ species, family = poisson, data = pop_bw) # run the glm and make it an object to summarize

summary(seed_bw)

ggplot(data = pop_bw, aes( y = seed_dens, fill = species)) +
  geom_boxplot() +
  labs(y = "Seedling density", x = "", title = " Beachwood ") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())




# CCA ----------------------------------------------------------------

# PCA -------------------------------------------------



pop_cca <- cca(pop1)
summary(pop_cca)



