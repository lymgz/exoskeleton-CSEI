# =============================================================================
# arm_network_visualization.R
# -----------------------------------------------------------------------------
# Association-rule mining (ARM) and dependency-network visualization for:
#
#   "Validation Logic and Transfer Limits in Construction Exoskeleton
#    Research: A Bibliometric and Association-Rule Analysis"
#
# This script reproduces the ARM results and the overall and setting-specific
# dependency networks reported in the manuscript.
#
# Input
#   ./data/step14_case_dataset.csv
#     Coded matrix of 97 retained studies x 10 ARM dimensions (D1-D10).
#     Deposited together with this script.
#

#
# Usage
#   1. Place the deposited CSV at ./data/step14_case_dataset.csv
#   2. Set the working directory to the folder containing this script.
#   3. source("step16_arm_network_visualization.R")
#
# Tested with R 4.3+ on Windows, macOS, and Linux.
# =============================================================================

rm(list = ls())

# -----------------------------------------------------------------------------
# 1. Packages
# -----------------------------------------------------------------------------
required_packages <- c(
  "arules", "arulesViz", "igraph", "tidygraph", "ggraph",
  "ggplot2", "ggrepel", "dplyr", "stringr", "tibble",
  "purrr", "patchwork", "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Please install the following packages first: ",
       paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(arules);   library(arulesViz); library(igraph)
  library(tidygraph); library(ggraph);   library(ggplot2);  library(ggrepel)
  library(dplyr);    library(stringr);   library(tibble);   library(purrr)
  library(patchwork); library(scales)
})

# -----------------------------------------------------------------------------
# 2. Paths and configuration
# -----------------------------------------------------------------------------
INPUT_CSV  <- file.path("data", "step14_case_dataset.csv")
OUTPUT_DIR <- "outputs"
PNG_DIR    <- file.path(OUTPUT_DIR, "png")
PDF_DIR    <- file.path(OUTPUT_DIR, "pdf")

for (d in c(OUTPUT_DIR, PNG_DIR, PDF_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ARM dimensions
ARM_COLS <- c(
  "ARM_D1_study_setting",
  "ARM_D2_task_dynamics",
  "ARM_D3_target_region",
  "ARM_D4_structural_form",
  "ARM_D5_actuation_type",
  "ARM_D6_control_strategy",
  "ARM_D7_primary_objective_evidence",
  "ARM_D8_objective_coverage",
  "ARM_D9_subjective_focus",
  "ARM_D10_embodied_sensing"
)
ARM_ALIAS_COLS <- paste0("D", seq_along(ARM_COLS))
names(ARM_ALIAS_COLS) <- ARM_COLS

ARM_DIMENSION_FAMILY <- c(
  D1 = "Setting / Context",       D2 = "Setting / Context",
  D3 = "Design / Configuration",  D4 = "Design / Configuration",
  D5 = "Design / Configuration",  D6 = "Design / Configuration",
  D7 = "Objective / Metric",      D8 = "Objective / Metric",
  D9 = "Objective / Metric",      D10 = "Sensing / Reporting"
)

EDGE_PALETTE <- c(
  "Setting / Context"      = "#5B7BD5",
  "Design / Configuration" = "#E6B84A",
  "Objective / Metric"     = "#E97B72",
  "Sensing / Reporting"    = "#6FB7D6"
)
NODE_FILL <- c(
  "Setting / Context"      = "#C9D2F0",
  "Design / Configuration" = "#E9E1C6",
  "Objective / Metric"     = "#E7D6D6",
  "Sensing / Reporting"    = "#D4E8EC"
)
NODE_BORDER <- c(
  "Setting / Context"      = "#7D8CC4",
  "Design / Configuration" = "#B99F55",
  "Objective / Metric"     = "#B89A9A",
  "Sensing / Reporting"    = "#6C9CAD"
)

# Apriori thresholds (match the manuscript)
APR_SUPPORT       <- 0.10
APR_CONFIDENCE    <- 0.90
APR_LIFT_MIN      <- 1.10
SUBGROUP_FALLBACK <- 0.08   # relaxed support used only if a subgroup is too sparse

# Network rendering
NETWORK_TOP_EDGES_OVERALL  <- 36
NETWORK_TOP_EDGES_SUBGROUP <- 20
MIN_RETAINED_RULES_NETWORK <- 4
MIN_NETWORK_EDGES          <- 3
NODE_LABEL_TOP_N           <- 22
NETWORK_LAYOUT_SEED        <- 20260322

PANEL_ORDER  <- c("Overall", "Laboratory", "Simulated", "On-Site")
PANEL_TITLES <- c(
  Overall = "(a) Overall",  Laboratory = "(b) Laboratory",
  Simulated = "(c) Simulated", "On-Site" = "(d) On-site"
)

LAYOUT_CENTERS <- tibble::tribble(
  ~family,                   ~x,    ~y,
  "Setting / Context",      -2.3,   1.8,
  "Design / Configuration", -0.4,   0.6,
  "Objective / Metric",      2.0,  -1.7,
  "Sensing / Reporting",     2.0,   1.6
)

# -----------------------------------------------------------------------------
# 3. Helper functions
# -----------------------------------------------------------------------------
parse_rule_side <- function(side) {
  side <- sub("^\\{", "", sub("\\}$", "", side))
  if (!nzchar(side)) return(character(0))
  trimws(unlist(strsplit(side, ",", fixed = TRUE)))
}

item_alias  <- function(x) sub("^(D[0-9]+)=.*$", "\\1", x)
item_family <- function(x) {
  fam <- ARM_DIMENSION_FAMILY[[item_alias(x)]]
  if (is.null(fam) || is.na(fam)) "Design / Configuration" else fam
}

apriori_filtered <- function(trans_obj, supp, conf, lift_min) {
  r <- arules::apriori(
    trans_obj,
    parameter = list(supp = supp, conf = conf, target = "rules")
  )
  if (length(r) == 0) return(r)
  subset(r, subset = lift >= lift_min)
}

retain_nonredundant <- function(rules_obj) {
  if (length(rules_obj) == 0) return(rules_obj)
  rules_obj[!arules::is.redundant(rules_obj)]
}

build_edges <- function(rules_obj, top_n, label) {
  if (length(rules_obj) == 0) return(tibble())
  rdf <- as(rules_obj, "data.frame")

  edges <- purrr::map_dfr(seq_len(nrow(rdf)), function(i) {
    parts <- strsplit(rdf$rules[i], " => ", fixed = TRUE)[[1]]
    if (length(parts) != 2) return(tibble())
    lhs <- parse_rule_side(parts[1]); rhs <- parse_rule_side(parts[2])
    if (!length(lhs) || !length(rhs)) return(tibble())
    expand.grid(from = lhs, to = rhs, stringsAsFactors = FALSE) |>
      tibble::as_tibble() |>
      dplyr::mutate(
        support     = rdf$support[i],
        confidence  = rdf$confidence[i],
        lift        = rdf$lift[i],
        count       = rdf$count[i],
        network_label = label,
        from_family = vapply(from, item_family, character(1)),
        to_family   = vapply(to,   item_family, character(1)),
        edge_group  = to_family
      )
  })
  if (nrow(edges) == 0) return(tibble())

  edges |>
    dplyr::group_by(from, to, edge_group, from_family, to_family, network_label) |>
    dplyr::summarise(
      rule_count  = dplyr::n(),
      weight      = sum(count, na.rm = TRUE),
      mean_lift   = mean(lift, na.rm = TRUE),
      max_support = max(support, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(edge_score = weight * mean_lift) |>
    dplyr::arrange(dplyr::desc(edge_score)) |>
    dplyr::slice_head(n = top_n)
}

build_nodes <- function(edge_df) {
  if (nrow(edge_df) == 0) return(tibble())
  dplyr::bind_rows(
    edge_df |> dplyr::transmute(name = from, family = from_family),
    edge_df |> dplyr::transmute(name = to,   family = to_family)
  ) |>
    dplyr::distinct() |>
    dplyr::left_join(
      edge_df |> dplyr::group_by(name = from) |>
        dplyr::summarise(out_w = sum(weight), .groups = "drop"),
      by = "name") |>
    dplyr::left_join(
      edge_df |> dplyr::group_by(name = to) |>
        dplyr::summarise(in_w = sum(weight), .groups = "drop"),
      by = "name") |>
    dplyr::mutate(
      out_w       = dplyr::coalesce(out_w, 0),
      in_w        = dplyr::coalesce(in_w, 0),
      node_weight = out_w + in_w,
      label = dplyr::if_else(
        rank(-node_weight, ties.method = "first") <= NODE_LABEL_TOP_N,
        name, NA_character_
      )
    )
}

make_layout <- function(node_df, edge_df, mode = "standard") {
  g <- igraph::graph_from_data_frame(
    edge_df |> dplyr::select(from, to, edge_group, weight, edge_score),
    directed = TRUE, vertices = node_df
  )

  if (identical(mode, "circular")) {
    ordered <- node_df |>
      dplyr::mutate(
        family = factor(family, levels = names(EDGE_PALETTE)),
        ow = -node_weight
      ) |>
      dplyr::arrange(family, ow, name) |>
      dplyr::mutate(
        idx = dplyr::row_number(),
        ang = seq(0, 2 * pi, length.out = dplyr::n() + 1)[idx],
        x   = cos(ang),
        y   = sin(ang)
      ) |>
      dplyr::select(name, x, y)
    coords <- tibble::tibble(name = igraph::V(g)$name) |>
      dplyr::left_join(ordered, by = "name")
  } else {
    seed_coords <- tibble::tibble(name = igraph::V(g)$name) |>
      dplyr::left_join(node_df, by = "name") |>
      dplyr::left_join(LAYOUT_CENTERS, by = "family") |>
      dplyr::mutate(
        x = x + dplyr::row_number() * 0.03,
        y = y - dplyr::row_number() * 0.02
      )
    set.seed(NETWORK_LAYOUT_SEED)
    m <- igraph::layout_with_fr(
      g, coords = as.matrix(seed_coords[, c("x", "y")]), niter = 1500
    )
    colnames(m) <- c("x", "y")
    coords <- tibble::as_tibble(m)
  }

  ggraph::create_layout(
    tidygraph::as_tbl_graph(g), layout = "manual",
    x = coords$x, y = coords$y
  )
}

network_theme <- function() {
  ggplot2::theme_void(base_family = "sans") +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = "#FFFFFF", colour = NA),
      panel.background = ggplot2::element_rect(fill = "#FFFFFF", colour = NA),
      plot.title       = ggplot2::element_text(face = "bold", size = 14, colour = "#222222"),
      plot.subtitle    = ggplot2::element_text(size = 10, colour = "#3A3A3A"),
      legend.position  = "right",
      legend.title     = ggplot2::element_text(size = 10, face = "bold"),
      legend.text      = ggplot2::element_text(size = 9),
      plot.margin      = ggplot2::margin(10, 10, 10, 10)
    )
}

placeholder_plot <- function(title_text, subtitle_text = NULL,
                             msg = "Insufficient retained rules for a meaningful network") {
  ggplot() +
    annotate("text", x = 0, y =  0.10, label = title_text, fontface = "bold", size = 5) +
    annotate("text", x = 0, y = -0.05, label = msg, size = 4.2, colour = "#3A3A3A") +
    coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1), expand = FALSE) +
    labs(title = title_text, subtitle = subtitle_text) +
    network_theme()
}

render_network <- function(layout_obj, edge_df, title_text, subtitle_text = NULL) {
  if (is.null(layout_obj) || nrow(edge_df) < MIN_NETWORK_EDGES) {
    return(placeholder_plot(title_text, subtitle_text))
  }
  ggraph::ggraph(layout_obj) +
    ggraph::geom_edge_link(
      aes(edge_colour = edge_group, edge_width = edge_score),
      alpha = 0.34,
      arrow = grid::arrow(length = grid::unit(2, "mm"), type = "closed"),
      end_cap   = ggraph::circle(2.5, "mm"),
      start_cap = ggraph::circle(2.5, "mm")
    ) +
    ggraph::geom_node_point(
      aes(fill = family, colour = family, size = node_weight),
      shape = 21, stroke = 0.8, alpha = 0.95
    ) +
    ggraph::geom_node_text(
      aes(label = label),
      repel = TRUE, family = "sans", size = 4.0, colour = "#3A3A3A",
      bg.colour = alpha("#FFFFFF", 0.7), max.overlaps = Inf, na.rm = TRUE
    ) +
    ggraph::scale_edge_colour_manual(values = EDGE_PALETTE, name = "Edge family") +
    ggplot2::scale_fill_manual(values = NODE_FILL, name = "Node family") +
    ggplot2::scale_colour_manual(values = NODE_BORDER, guide = "none") +
    ggplot2::scale_size_continuous(range = c(3.2, 8.5), name = "Node weight") +
    ggraph::scale_edge_width_continuous(range = c(0.4, 1.8), guide = "none") +
    labs(title = title_text, subtitle = subtitle_text) +
    network_theme()
}

save_plot <- function(p, name, w = 13, h = 9, dpi = 600) {
  ggplot2::ggsave(
    file.path(PNG_DIR, paste0(name, ".png")),
    p, width = w, height = h, units = "in", dpi = dpi, bg = "#FFFFFF"
  )
  ggplot2::ggsave(
    file.path(PDF_DIR, paste0(name, ".pdf")),
    p, width = w, height = h, units = "in",
    device = grDevices::cairo_pdf, bg = "#FFFFFF"
  )
}

build_network <- function(rules_obj, top_n, label, mode = "standard") {
  edge_df <- build_edges(rules_obj, top_n, label)
  if (nrow(edge_df) == 0) {
    return(list(edge_df = edge_df, node_df = tibble(), layout = NULL))
  }
  node_df <- build_nodes(edge_df)
  layout  <- make_layout(node_df, edge_df, mode = mode)
  list(edge_df = edge_df, node_df = node_df, layout = layout)
}

sanitize_name <- function(x) {
  x |> stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]+", "_") |>
    stringr::str_replace_all("^_|_$", "")
}

# -----------------------------------------------------------------------------
# 4. Load and prepare data
# -----------------------------------------------------------------------------
if (!file.exists(INPUT_CSV)) {
  stop("Input CSV not found: ", INPUT_CSV,
       "\nPlace the deposited coded matrix at this path before running.")
}

case_records <- read.csv(
  INPUT_CSV, stringsAsFactors = FALSE,
  fileEncoding = "UTF-8", check.names = FALSE
)

missing_cols <- setdiff(ARM_COLS, names(case_records))
if (length(missing_cols) > 0) {
  stop("Missing ARM columns in input: ", paste(missing_cols, collapse = ", "))
}
case_records <- case_records[complete.cases(case_records[, ARM_COLS]), , drop = FALSE]
cat(sprintf("Loaded %d complete cases from %s\n", nrow(case_records), INPUT_CSV))

trans_input <- case_records[, ARM_COLS, drop = FALSE]
names(trans_input) <- ARM_ALIAS_COLS[ARM_COLS]
trans_input[] <- lapply(trans_input, factor)
trans <- arules::transactions(trans_input)

# -----------------------------------------------------------------------------
# 5. Overall Apriori and retained-rule set
# -----------------------------------------------------------------------------
rules <- apriori_filtered(trans, APR_SUPPORT, APR_CONFIDENCE, APR_LIFT_MIN)
rules_sorted <- sort(rules, by = "confidence", decreasing = TRUE)
write.csv(
  as(rules_sorted, "data.frame"),
  file.path(OUTPUT_DIR, "step14_all_association_rules.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

retained_rules <- retain_nonredundant(rules_sorted)
cat(sprintf("All rules: %d  |  Retained non-redundant: %d\n",
            length(rules_sorted), length(retained_rules)))

# -----------------------------------------------------------------------------
# 6. Overall dependency network (standard + circular layouts)
# -----------------------------------------------------------------------------
build_and_save_overall <- function(mode, file_stub) {
  net <- build_network(retained_rules, NETWORK_TOP_EDGES_OVERALL, "Overall", mode = mode)
  subtitle <- sprintf(
    "%s retained-rule dependency network | top %d edges | retained rules = %d",
    if (mode == "circular") "Circular" else "Aggregated",
    NETWORK_TOP_EDGES_OVERALL, length(retained_rules)
  )
  p <- render_network(net$layout, net$edge_df, PANEL_TITLES[["Overall"]], subtitle)
  save_plot(p, file_stub)
  list(plot = p, net = net)
}
overall_std <- build_and_save_overall("standard", "network_overall_full")
overall_cir <- build_and_save_overall("circular", "network_overall_full_circular")

# -----------------------------------------------------------------------------
# 7. Subgroup rules grouped by study setting (CSV export only)
# -----------------------------------------------------------------------------
sub_rule_list <- list()
for (ss in unique(case_records$ARM_D1_study_setting)) {
  sdf <- subset(case_records, ARM_D1_study_setting == ss)
  sin <- sdf[, ARM_COLS, drop = FALSE]
  names(sin) <- ARM_ALIAS_COLS[ARM_COLS]
  sin[] <- lapply(sin, factor)
  st <- arules::transactions(sin)
  sr <- apriori_filtered(st, APR_SUPPORT, APR_CONFIDENCE, APR_LIFT_MIN)
  if (length(sr) == 0) next
  d <- as(sort(sr, by = "confidence", decreasing = TRUE), "data.frame")
  d$ARM_D1_study_setting <- ss
  sub_rule_list[[ss]] <- d
}
if (length(sub_rule_list) > 0) {
  out <- do.call(rbind, sub_rule_list)
  out <- out[order(out$ARM_D1_study_setting, -out$confidence, -out$lift), ]
  write.csv(
    out, file.path(OUTPUT_DIR, "step14_rules_sorted_by_study_setting.csv"),
    row.names = FALSE, fileEncoding = "UTF-8"
  )
}

# -----------------------------------------------------------------------------
# 8. Setting-specific subgroup networks (D1 excluded inside each subgroup)
# -----------------------------------------------------------------------------
build_subgroup_net <- function(setting_label, mode) {
  sdf <- subset(case_records, ARM_D1_study_setting == setting_label)
  if (nrow(sdf) == 0) {
    return(list(plot = placeholder_plot(
      PANEL_TITLES[[setting_label]], NULL, "No rows in this setting"
    )))
  }
  cols  <- setdiff(ARM_COLS, "ARM_D1_study_setting")
  alias <- ARM_ALIAS_COLS[cols]
  sin <- sdf[, cols, drop = FALSE]
  names(sin) <- alias
  sin[] <- lapply(sin, factor)
  st <- arules::transactions(sin)

  used_supp <- APR_SUPPORT
  fallback  <- FALSE
  sr <- retain_nonredundant(apriori_filtered(st, used_supp, APR_CONFIDENCE, APR_LIFT_MIN))
  if (length(sr) < MIN_RETAINED_RULES_NETWORK) {
    used_supp <- SUBGROUP_FALLBACK
    fallback  <- TRUE
    warning(sprintf("Subgroup '%s' sparse; retrying with support %.2f.",
                    setting_label, used_supp))
    sr <- retain_nonredundant(apriori_filtered(st, used_supp, APR_CONFIDENCE, APR_LIFT_MIN))
  }

  net <- build_network(sr, NETWORK_TOP_EDGES_SUBGROUP, setting_label, mode = mode)
  parts <- c(
    "Within-setting dependency network",
    sprintf("support >= %.2f", used_supp),
    sprintf("retained = %d", length(sr))
  )
  if (fallback) parts <- c(parts, "fallback support used")
  p <- render_network(net$layout, net$edge_df, PANEL_TITLES[[setting_label]],
                      paste(parts, collapse = " | "))
  list(plot = p, net = net, retained_rules = sr)
}

subgroup_std <- list(); subgroup_cir <- list()
for (s in c("Laboratory", "Simulated", "On-Site")) {
  cat(sprintf("\nBuilding subgroup network: %s\n", s))
  subgroup_std[[s]] <- build_subgroup_net(s, "standard")
  subgroup_cir[[s]] <- build_subgroup_net(s, "circular")
  stub <- sanitize_name(s)
  save_plot(subgroup_std[[s]]$plot, paste0("network_setting_", stub))
  save_plot(subgroup_cir[[s]]$plot, paste0("network_setting_", stub, "_circular"))
}

# -----------------------------------------------------------------------------
# 9. 2x2 panel plots (overall + three settings)
# -----------------------------------------------------------------------------
assemble_panel <- function(plots, title_text) {
  ordered <- lapply(PANEL_ORDER, function(n) plots[[n]])
  patchwork::wrap_plots(ordered, ncol = 2) +
    patchwork::plot_annotation(
      title = title_text,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 15,
                                           colour = "#222222", hjust = 0.5),
        plot.background = ggplot2::element_rect(fill = "#FFFFFF", colour = NA)
      )
    )
}

panel_std <- assemble_panel(
  c(list(Overall = overall_std$plot),
    lapply(subgroup_std, function(x) x$plot)),
  "Overall and setting-specific ARM dependency networks"
)
panel_cir <- assemble_panel(
  c(list(Overall = overall_cir$plot),
    lapply(subgroup_cir, function(x) x$plot)),
  "Overall and setting-specific ARM dependency networks (circular layout)"
)
save_plot(panel_std, "network_2x2_overall_plus_settings",          w = 14, h = 10)
save_plot(panel_cir, "network_2x2_overall_plus_settings_circular", w = 14, h = 10)

# -----------------------------------------------------------------------------
# 10. Summary
# -----------------------------------------------------------------------------
cat("\n=== Summary ===\n")
cat(sprintf("Input file       : %s\n", INPUT_CSV))
cat(sprintf("Cases            : %d\n", nrow(case_records)))
cat(sprintf("Transactions     : %d\n", length(trans)))
cat(sprintf("All rules        : %d\n", length(rules_sorted)))
cat(sprintf("Retained rules   : %d\n", length(retained_rules)))
cat(sprintf("Outputs written  : %s\n", normalizePath(OUTPUT_DIR)))

# End of script
