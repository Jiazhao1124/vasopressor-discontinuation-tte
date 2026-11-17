## Figure S2: subgroup forest plot, MIMIC-IV (overlap-weighted risk differences)
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(ragg); library(svglite)})
fam <- "Liberation Sans"
fmt <- function(x) sub("^-", "−", sprintf("%.1f", x))
ci  <- function(lo, hi) if (lo < 0) paste0(fmt(lo), " to ", fmt(hi)) else paste0(fmt(lo), "–", fmt(hi))

res <- read.csv("../results/main_results.csv", stringsAsFactors = FALSE)
pick <- function(key) res[res$analysis == key, ]
spec <- list(c("All patients", "M primary OW"), c("Sepsis-3", NA), c("  Yes", "M sub: Sepsis-3"),
             c("  No", "M sub: No sepsis-3"), c("ICU type", NA), c("  Cardiovascular", "M sub: Cardiovascular ICU"),
             c("  Other", "M sub: Other ICUs"), c("Vasopressin at discontinuation", NA),
             c("  Yes", "M sub: Vasopressin at stop"), c("  No", "M sub: No vasopressin at stop"),
             c("Single agent in last 2 h", "M sub: Single agent, last 2 h"), c("Episode duration", NA),
             c("  \u226524 h", "M sub: Episode >=24 h"), c("  <24 h", "M sub: Episode <24 h"),
             c("Invasive ventilation at discontinuation", NA), c("  Yes", "M sub: IMV at stop"),
             c("  No", "M sub: No IMV at stop"))
r <- do.call(rbind, lapply(spec, function(s) {
  if (is.na(s[2])) return(data.frame(lab = s[1], n = NA, rd = NA, lo = NA, hi = NA))
  x <- pick(s[2])
  data.frame(lab = s[1], n = x$n, rd = round(100 * x$rd, 1), lo = round(100 * x$rd_lo, 1), hi = round(100 * x$rd_hi, 1))
}))
r$y <- rev(seq_len(nrow(r)))
r$head <- is.na(r$n) | r$lab == "All patients"
r$ntxt <- ifelse(is.na(r$n), "", formatC(r$n, big.mark = ",", format = "d"))
r$etxt <- mapply(function(a, b, c) if (is.na(a)) "" else paste0(fmt(a), " (", ci(b, c), ")"), r$rd, r$lo, r$hi)
r$face <- ifelse(r$head, "bold", "plain")
top <- max(r$y) + 1.1

base <- theme_void(base_family = fam) + theme(plot.margin = margin(0, 0, 0, 0))
txt <- function(col, x, hj, header) {
  ggplot(r, aes(y = y)) +
    geom_text(aes(x = x, label = .data[[col]], fontface = face), hjust = hj, size = 7/.pt, family = fam) +
    annotate("text", x = x, y = top, label = header, hjust = hj, fontface = "bold", size = 7/.pt, family = fam) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.4, top + 0.5), expand = c(0, 0)) + base
}
left  <- txt("lab",  0,    0, "Subgroup")
mid_n <- txt("ntxt", 0.95, 1, "ICU stays")
right <- txt("etxt", 1,    1, "Risk difference (95% CI)")

f <- subset(r, !is.na(rd))
forest <- ggplot(f, aes(y = y)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = r$rd[1], linewidth = 0.25, linetype = "22", colour = "grey55") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.35) +
  geom_point(aes(x = rd, shape = lab == "All patients", size = lab == "All patients"), fill = "black", stroke = 0) +
  scale_shape_manual(values = c(15, 18), guide = "none") +
  scale_size_manual(values = c(1.4, 2.6), guide = "none") +
  scale_x_continuous(limits = c(-5, 25), breaks = seq(-5, 25, 5),
                     labels = function(x) sub("^-", "−", x), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, top + 0.5), expand = c(0, 0)) +
  labs(x = "Abrupt minus tapered (percentage points)") +
  theme_classic(base_size = 7, base_family = fam) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        axis.title.y = element_blank(), axis.line.x = element_line(linewidth = 0.3),
        axis.ticks.x = element_line(linewidth = 0.3), axis.ticks.length = unit(1, "mm"),
        axis.text.x = element_text(colour = "black", size = 7), axis.title.x = element_text(size = 7),
        plot.margin = margin(0, 2, 0, 2, "mm"))

p <- left + mid_n + forest + right + plot_layout(widths = c(4.4, 1.3, 4.0, 2.6))
w <- 170; h <- 100
ggsave("output/FigureS2.png",  p, device = agg_png,  width = w, height = h, units = "mm", dpi = 600)
ggsave("output/FigureS2.tiff", p, device = agg_tiff, width = w, height = h, units = "mm", dpi = 600, compression = "lzw")
ggsave("output/FigureS2.pdf",  p, device = cairo_pdf, width = w, height = h, units = "mm")
svglite("output/FigureS2.svg", width = w/25.4, height = h/25.4,
        system_fonts = list(sans = "Liberation Sans"), web_fonts = list())
print(p); invisible(dev.off())
