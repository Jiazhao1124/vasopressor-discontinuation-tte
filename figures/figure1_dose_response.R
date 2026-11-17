## Figure 1: standardised dose-response, MIMIC-IV and eICU-CRD
suppressPackageStartupMessages({library(ggplot2); library(ragg); library(svglite)})
fam <- "Liberation Sans"   # metric-identical to Arial; SVG/PDF name the face Arial
m <- read.csv("../results/dose_response_mimic.csv"); m$db <- "MIMIC-IV (n = 11,518)"
e <- read.csv("../results/dose_response_eicu.csv");  e$db <- "eICU-CRD (n = 1,852)"
d <- rbind(m, e); d$db <- factor(d$db, levels = unique(d$db))
d[c("risk","lo","hi")] <- d[c("risk","lo","hi")] * 100
cols <- c("#0072B5", "#BC3C29")

p <- ggplot(d, aes(end_dose, risk, colour = db, fill = db)) +
  geom_vline(xintercept = 0.05, linewidth = 0.25, linetype = "22", colour = "grey55") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.45) +
  geom_point(aes(shape = db), size = 1.3, stroke = 0.35, fill = "white") +
  scale_shape_manual(values = c(21, 22)) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  scale_x_continuous(limits = c(0, 0.205), breaks = seq(0, 0.20, 0.05),
                     labels = c("0", "0.05", "0.10", "0.15", "0.20"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(10, 70), breaks = seq(10, 70, 10), expand = c(0, 0)) +
  labs(x = "Norepinephrine-equivalent dose at discontinuation (µg/kg/min)",
       y = "Vasopressor re-initiation within 24 h (%)") +
  theme_classic(base_size = 8, base_family = fam) +
  theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        axis.ticks.length = unit(1.2, "mm"), axis.text = element_text(colour = "black", size = 7),
        axis.title = element_text(size = 8),
        legend.position = c(0.98, 0.04), legend.justification = c(1, 0),
        legend.title = element_blank(), legend.text = element_text(size = 7),
        legend.key.width = unit(5, "mm"), legend.key.height = unit(3.2, "mm"),
        legend.background = element_blank(), plot.margin = margin(2, 4, 1, 1, "mm"))

w <- 120; h <- 85
ggsave("output/Figure1.png",  p, device = agg_png,  width = w, height = h, units = "mm", dpi = 600)
ggsave("output/Figure1.tiff", p, device = agg_tiff, width = w, height = h, units = "mm", dpi = 600, compression = "lzw")
ggsave("output/Figure1.pdf",  p, device = cairo_pdf, width = w, height = h, units = "mm")
svglite("output/Figure1.svg", width = w/25.4, height = h/25.4,
        system_fonts = list(sans = "Liberation Sans"), web_fonts = list())
print(p); invisible(dev.off())
