# Data analysis: Profile data ISRaD #
# Relationship between 14C and depth/SOC #
# Sophie von Fromm #
# 08/03/2023 #

## Derivatives of carbon with depth (based on Carlos' idea)

library(tidyverse)
library(tikzDevice)

# Load splined and filtered data
lyr_data <- read_csv("../Data/ISRaD_flat_splined_filled_2023-02-08.csv")

# Calculate derivatives for 14C and SOC
lyr_data_deri <- lyr_data %>% 
  dplyr::select(id, UD, LD, CORG_msp, lyr_14c_msp, ClimateZoneAnd, MineralType,
                site_name) %>% 
  group_by(id) %>% 
  mutate(CORG_d1 = lead(CORG_msp, 1) - CORG_msp,
         CORG_d2 = lead(CORG_d1, 1) - CORG_d1) %>% 
  ungroup()

# Calculate median SOC profile for each climate zone
climate_c <- lyr_data_deri %>%
  dplyr::group_by(ClimateZoneAnd, UD) %>% 
  dplyr::mutate(Derivative_1 = median(CORG_d1, na.rm = TRUE),
                Derivative_2 = median(CORG_d2, na.rm = TRUE),
                median_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$estimate,
                lci_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[1],
                uci_c = wilcox.test(CORG_msp, conf.level = 0.95, conf.int = TRUE)$conf.int[2],
                n = n(),
                n_site = n_distinct(site_name)) %>% 
  distinct(ClimateZoneAnd, UD, .keep_all = TRUE) %>%
  ungroup(UD) %>%
  dplyr::mutate(n_rel = n * 100 / max(n))

climate_c$ClimateZoneAnd <- factor(climate_c$ClimateZoneAnd,
                                   levels = c("volcanic soils", "tundra/polar", "cold temperate", 
                                              "warm temperate", "arid", "tropical"))

# Plot splined data: Figure 6
tikz(file="../Figures/ISRaD_Profiles.tex", standAlone = TRUE, width=14, height=7)
climate_c %>% 
  #remove data with less than 5 profiles,less than 3 sites and less than 33% of the total data
  dplyr::filter(n > 4 & n_rel > 33 & n_site > 2) %>% 
  ggplot() + 
  geom_path(aes(y = UD, x = median_c, color = ClimateZoneAnd), linewidth = 2) +
  geom_errorbarh(aes(xmin = lci_c, xmax = uci_c, y = UD, 
                     color = ClimateZoneAnd), alpha = 0.3) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position = "top",
        legend.background = element_blank()) +
  scale_y_reverse("Depth (cm)", expand = c(0,0), limits = c(100,0)) +
  scale_x_continuous("Soil organic carbon (wt-\\%)", expand = c(0,0),
                     limits = c(0,30)) +
  scale_color_discrete("") +
  theme(axis.line = element_line(color = 'black'))
dev.off()

ggsave(file = paste0("./Figure/ISRaD_msp_SOC_climate_depth_profile_deriv_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)


#Plot first derivative
c_depth_1 <- climate_c %>% 
  drop_na() %>% 
  filter(n > 4 & n_rel > 33 & n_site > 3) %>% 
  ggplot() + 
  geom_path(aes(x = Derivative_1, y = UD, color = ClimateZoneAnd), linewidth = 1.5) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position = "top",
        legend.background = element_blank(), 
        plot.margin = margin(l = 25, t = 25, b = 25)) +
  scale_x_continuous(expression(paste("Change in soil organic carbon (wt-\\% c", m^-1,")")), 
                     limits = c(-0.6, 0.05), expand = c(0,0),
                     breaks = seq(-0.6,0,0.2)) +
  scale_y_reverse("Depth (cm)", limits = c(100,0), expand = c(0,0)) +
  scale_color_discrete("") +
theme(axis.line = element_line(color = 'black'))

# Plot second derivative
c_depth_2 <- climate_c %>% 
  drop_na() %>% 
  filter(n > 4 & n_rel > 33 & n_site > 2) %>% 
  ggplot() + 
  geom_path(aes(x = Derivative_2, y = UD, color = ClimateZoneAnd), linewidth = 1.5) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(color = "black"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position = "top",
        legend.background = element_blank(), 
        plot.margin = margin(r = 25, t = 25, b = 25)) +
  scale_x_continuous(expression(paste("Change in 1st derivative (wt-\\% c", m^-2,")")), 
                     expand = c(0,0), limits = c(-0.09, 0.03)) +
  scale_y_reverse("", limits = c(100,0), expand = c(0,0)) +
  scale_color_discrete("") +
theme(axis.line = element_line(color = 'black'))

# Plot both together: Figure 7
tikz("../Figures/ISRaD_derivatives.tex", standAlone = TRUE, width=14)
ggpubr::ggarrange(c_depth_1, c_depth_2, common.legend = TRUE, 
                  labels = c("1st derivative", "2nd derivative"), 
                  font.label = list(face = "plain"), hjust = c(-0.8, -0.5))
dev.off()

ggsave(file = paste0("./Figure/ISRaD_msp_14C_SOC_climate_depth_deriv_", Sys.Date(),
                     ".jpeg"), width = 12, height = 6)


