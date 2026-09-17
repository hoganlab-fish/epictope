library(shiny)
library(bslib)      # tooltip() UI helper
library(plotly)     # linked combined view -- install.packages("plotly") if missing
library(epictope)
library(NGLVieweR)  # 3D structure viewer -- install.packages("NGLVieweR") if missing

conda_lib <- file.path(Sys.getenv("CONDA_PREFIX"), "lib")
Sys.setenv(LD_LIBRARY_PATH = paste0(conda_lib, ":", Sys.getenv("LD_LIBRARY_PATH")))
print(Sys.getenv("LD_LIBRARY_PATH")); print(Sys.which("mkdssp"))

# --- Patch: epictope:::calculate_scores has a tie-handling bug. ---
# The original does:
#   c(...)[which(.x == min(.x))]
# which returns EVERY feature tied for the minimum, producing a
# variable-length (list-column) result instead of one name per row.
# That list column then crashes renderTable()'s xtable formatting with
# "Error in is.finite: default method not implemented for type 'list'".
#
# Fix: swap to which.min() (always returns exactly one index) without
# touching the installed package -- we copy the real function and patch
# only the offending line, keeping its original namespace environment
# so internal objects (max_sasa, ss_convert, weights, etc.) still resolve.
calculate_scores <- epictope:::calculate_scores
body_text <- paste(deparse(body(calculate_scores)), collapse = "\n")  # preserves original statement structure
# \\s tolerates an internal line-wrap (deparse() can split long expressions
# across lines) without touching anything else in the function body.
target_pattern <- "which\\(\\s*\\.x\\s*==\\s*min\\(\\s*\\.x\\s*\\)\\s*\\)"
if (!grepl(target_pattern, body_text)) {
  stop("calculate_scores patch failed: target pattern not found in installed epictope source.")
}
body_text <- sub(target_pattern, "which.min(.x)", body_text)
body(calculate_scores) <- parse(text = body_text)[[1]]

# --- Protein reference list for the searchable dropdown ---------------------
# zebrafish_proteins.csv must sit next to app.R. It's a starter reference
# list (UniProt reviewed Danio rerio entries: uniprot_id, gene_name,
# protein_name, length, alphafold_id) -- not necessarily the complete
# zebrafish proteome. The dropdown still accepts free-text UniProt IDs for
# anything not in the list (options = list(create = TRUE) below), so this
# never blocks a query -- it's just a convenience list.
protein_ref <- tryCatch(
  read.csv("zebrafish_proteins.csv", stringsAsFactors = FALSE),
  error = function(e) NULL
)

protein_choices <- if (!is.null(protein_ref) && nrow(protein_ref) > 0) {
  labels <- paste0(
    ifelse(nzchar(protein_ref$gene_name), protein_ref$gene_name, protein_ref$uniprot_id),
    " — ", protein_ref$protein_name,
    " (", protein_ref$uniprot_id, ")"
  )
  setNames(protein_ref$uniprot_id, labels)
} else {
  NULL
}

# --- Small plotting/analysis helpers -----------------------------------------

# Sliding-window average matching the paper's Fig 6C smoothing: a 7-residue
# window, shrinking to 4-6 residues near the sequence termini rather than
# padding/wrapping.
sliding_avg <- function(x, window = 7) {
  n <- length(x)
  half <- window %/% 2
  vapply(seq_len(n), function(i) {
    lo <- max(1, i - half)
    hi <- min(n, i + half)
    mean(x[lo:hi], na.rm = TRUE)
  }, numeric(1))
}

# Greedy local-maxima picker: finds up to k points that are local maxima of
# `val` -- strictly the highest of their immediate neighbors -- at least
# `min_sep` positions apart, chosen best-score-first so several picks never
# end up crowded on the same hump.
find_extrema <- function(pos, val, k, min_sep) {
  n <- length(val)
  if (n < 3 || k <= 0) return(data.frame(position = numeric(0), score = numeric(0)))
  is_extreme <- vapply(seq_len(n), function(i) {
    lo <- max(1, i - 1); hi <- min(n, i + 1)
    val[i] == max(val[lo:hi])
  }, logical(1))
  cand <- which(is_extreme)
  cand <- cand[order(val[cand], decreasing = TRUE)]
  picked <- integer(0)
  for (idx in cand) {
    if (length(picked) >= k) break
    if (length(picked) == 0 || all(abs(pos[picked] - pos[idx]) >= min_sep)) {
      picked <- c(picked, idx)
    }
  }
  data.frame(position = pos[picked], score = val[picked], stringsAsFactors = FALSE)
}

# Pick candidate tag sites: N-term, C-term, plus the top-scoring local
# maxima ("peaks" -- good tagging candidates) of the smoothed min-score
# curve. These are the app's default "regions of interest" -- highlighted
# everywhere (score plot, alignment, 3D structure) before the user clicks
# anything, and listed at the top of the Combined view tab for quick jumps.
pick_tag_sites <- function(final_df, edge_margin = 10, n_extrema = 5, min_sep = 15) {
  final_df <- final_df[order(final_df$position), ]
  smoothed <- sliding_avg(final_df$min, window = 7)
  pos <- final_df$position
  n <- nrow(final_df)
  interior <- which(pos > edge_margin & pos < (max(pos) - edge_margin))
  if (length(interior) == 0) interior <- seq_len(n)

  peaks <- find_extrema(pos[interior], smoothed[interior], n_extrema, min_sep)

  list(
    n_term = pos[1],
    c_term = pos[n],
    peaks = peaks,         # data.frame(position, score), best first
    position_vec = pos,    # aligned with `smoothed`, for later lookups
    smoothed = smoothed
  )
}

# Standard Clustal-style amino-acid coloring so the alignment panel reads
# like a real MSA viewer (Jalview/MView/muscle-style): hydrophobic = blue,
# positive = red, negative = magenta, polar = green, glycine = orange,
# proline = yellow, aromatic (His/Tyr) = cyan, gaps = light gray. Anything
# outside the 20 standard residues (e.g. "X") falls back to neutral gray.
aa_colors <- c(
  W = "#3B7FCC", L = "#3B7FCC", V = "#3B7FCC", I = "#3B7FCC", M = "#3B7FCC",
  A = "#3B7FCC", F = "#3B7FCC", C = "#3B7FCC",
  K = "#E6194B", R = "#E6194B",
  E = "#B10DC9", D = "#B10DC9",
  N = "#2ECC71", Q = "#2ECC71", S = "#2ECC71", T = "#2ECC71",
  G = "#FF851B",
  P = "#FFDC00",
  H = "#39CCCC", Y = "#39CCCC",
  "-" = "#F2F2F2"
)

# Shared "meaning of color" across every panel (score plot, alignment, 3D):
# amber = sequence termini, green = top scoring peaks (good tag candidates),
# red = curated UniProt binding-site residues (avoid tagging over these),
# blue = whatever the user last clicked or range-selected on any panel.
# Every one of these is used as the ONE canonical hex for that meaning
# everywhere it appears (score plot bands/lines/text, MSA outline boxes,
# 3D structure, downloadable report, legends) -- no separate "darker" or
# "lighter" variant per panel, so the same meaning never looks like two
# different colors depending on which plot you're looking at. Chosen to
# also stay clear of every `aa_colors` hue (peak green used to be the exact
# same hex as the "polar" amino acid color; binding red used to be very
# close to the "positive" amino acid color) and of the feature-plot trace
# colors (black/brown/forestgreen/teal).
COLOR_DEFAULT_SITE <- "#FFC107"
COLOR_PEAK          <- "#00C853"
COLOR_BINDING       <- "#C62828"
COLOR_USER_TAG      <- "#2979FF"

# A small colored square + label, for building color-key legends.
color_swatch <- function(color, label) {
  tags$span(
    style = "display:inline-flex; align-items:center; margin-right:14px; margin-bottom:4px; font-size:12px; color:#333; white-space:nowrap;",
    tags$span(style = paste0("display:inline-block; width:12px; height:12px; background:", color,
                              "; border:1px solid #333; border-radius:2px; margin-right:5px; flex-shrink:0;")),
    label
  )
}

# Site-highlight legend (termini/peaks/binding-site/user selection) -- the
# same four colors and meanings on the score plot, the alignment, and the
# 3D structure, so one definition is reused everywhere it applies.
site_legend_ui <- function() {
  tags$div(style = "margin: 4px 0 8px 0;",
    color_swatch(COLOR_DEFAULT_SITE, "N-term / C-term"),
    color_swatch(COLOR_PEAK, "Top peak (good tag candidate)"),
    color_swatch(COLOR_BINDING, "UniProt binding site (avoid)"),
    color_swatch(COLOR_USER_TAG, "Your current selection")
  )
}

# Amino-acid biochemical-property legend for the alignment panel's coloring.
aa_legend_ui <- function() {
  tags$div(style = "margin: 4px 0 8px 0;",
    color_swatch("#3B7FCC", "Hydrophobic (A,V,L,I,M,F,W,C)"),
    color_swatch("#E6194B", "Positive (K,R)"),
    color_swatch("#B10DC9", "Negative (D,E)"),
    color_swatch("#2ECC71", "Polar (N,Q,S,T)"),
    color_swatch("#FF851B", "Glycine (G)"),
    color_swatch("#FFDC00", "Proline (P)"),
    color_swatch("#39CCCC", "Aromatic (H,Y)"),
    color_swatch("#F2F2F2", "Gap (-)")
  )
}

# Collapse a sorted vector of local (in-row) indices into contiguous
# [start, end] runs -- so a stretch of adjacent highlighted residues (e.g.
# a multi-residue binding motif, or a dragged range selection) draws as one
# outline box spanning the whole run, instead of one small square per
# residue.
contiguous_runs <- function(idx) {
  if (length(idx) == 0) return(list())
  idx <- sort(unique(idx))
  breaks <- c(0, which(diff(idx) > 1), length(idx))
  starts <- idx[breaks[-length(breaks)] + 1]
  ends   <- idx[breaks[-1]]
  Map(function(s, e) c(start = s, end = e), starts, ends)
}

# Classic alignment "wrap" width (matches typical Clustal/EMBOSS text output)
# -- used only in "show full alignment" mode, where long alignments are
# shown as stacked blocks of this many residues instead of one endlessly
# wide row, so everything fits on screen at once.
MSA_WRAP_WIDTH <- 60

# Default (non-full) alignment view: a fixed 50-residue-wide strip centered
# on the current point of interest -- deliberately NOT zoomable (see the
# scrollable single-row renderer below); the user scrolls sideways to see
# more instead of zooming out.
MSA_DEFAULT_HALF_WIDTH <- 25

# Build the (label, id, position, kind) list for every default site --
# termini and top peaks -- shared by the top-level nav bar, the score plot,
# the alignment, the 3D view, and the downloadable report, so they never
# drift apart. `kind` drives which color a site gets.
site_list_for <- function(ts) {
  term_pos <- c(ts$n_term, ts$c_term)
  term_df <- data.frame(
    id = c("term_n", "term_c"), label = c("N-term", "C-term"),
    position = term_pos, score = ts$smoothed[match(term_pos, ts$position_vec)],
    kind = "term", stringsAsFactors = FALSE
  )
  peak_df <- if (nrow(ts$peaks) > 0) {
    data.frame(id = paste0("peak_", seq_len(nrow(ts$peaks))),
               label = paste0("Peak ", seq_len(nrow(ts$peaks))),
               position = ts$peaks$position, score = ts$peaks$score,
               kind = "peak", stringsAsFactors = FALSE)
  } else NULL
  do.call(rbind, c(list(term_df), list(peak_df)))
}

# Parse a UniProt TSV "Binding site" cell (field ft_binding) into a sorted,
# de-duplicated vector of residue positions. UniProt packs multiple
# features into one string like:
#   "BINDING 45; /ligand=\"Zn(2+)\"; ...; BINDING 48..52; /ligand=\"ATP\"; ..."
# -- we only need the position/range after each "BINDING" keyword.
parse_uniprot_binding <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) return(integer(0))
  hits <- regmatches(x, gregexpr("BINDING\\s+\\d+(\\.\\.\\d+)?", x, perl = TRUE))[[1]]
  if (length(hits) == 0) return(integer(0))
  ranges <- lapply(hits, function(h) {
    nums <- as.integer(regmatches(h, gregexpr("\\d+", h))[[1]])
    if (length(nums) == 1) nums else seq(nums[1], nums[2])
  })
  sort(unique(unlist(ranges)))
}

# Collapse a sorted vector of positions into contiguous [start, end] runs,
# so a stretch of adjacent binding-site residues draws as one shaded band
# instead of one sliver per residue.
positions_to_runs <- function(pos) {
  if (length(pos) == 0) return(data.frame(start = integer(0), end = integer(0)))
  pos <- sort(unique(pos))
  breaks <- c(0, which(diff(pos) > 1), length(pos))
  data.frame(
    start = pos[breaks[-length(breaks)] + 1],
    end   = pos[breaks[-1]]
  )
}

# Human-readable "45, 48-52, 90" summary of a set of binding-site residues,
# shared by the nav bar and the downloadable HTML report.
format_binding_summary <- function(binding_sites) {
  if (length(binding_sites) == 0) return("None found")
  runs <- positions_to_runs(binding_sites)
  paste(apply(runs, 1, function(r) {
    if (r["start"] == r["end"]) as.character(r["start"]) else sprintf("%d-%d", r["start"], r["end"])
  }), collapse = ", ")
}

# A jump-to-position button styled like an actionButton but with a
# dynamically-set target -- click sends `pos` to a single shared
# `input$jump_pos`, so an arbitrary, run-time-determined number of buttons
# (N-term/C-term/peaks) can all be wired up without pre-declaring a
# static input id (and matching observer) for each one.
make_jump_button <- function(label, pos, color, tip) {
  btn <- tags$button(
    label, type = "button", class = "btn btn-default action-button",
    style = paste0("background-color:", color,
                   "; border: 1px solid #333; margin-right: 8px; margin-bottom: 6px;"),
    onclick = sprintf("Shiny.setInputValue('jump_pos', %d, {priority: 'event'});", pos)
  )
  tooltip(btn, tip)
}

# One-line, consistently-worded banner describing the active user tag --
# repeated at the top of every panel so the current selection is never
# ambiguous, on top of the color highlighting itself.
tag_banner_text <- function(rng) {
  if (is.null(rng)) return("No residue tagged yet -- click or box-select any panel to tag one.")
  if (rng[1] == rng[2]) sprintf("Tagged: residue %d", rng[1])
  else sprintf("Tagged range: residues %d–%d", rng[1], rng[2])
}

ui <- fluidPage(
  theme = bs_theme(version = 5),
  tags$head(tags$script(HTML(
    "
    // NGLVieweR (and other WebGL/canvas widgets) initialize at zero size when
    // their tab starts hidden, and don't repaint on their own when you switch
    // to that tab -- force a window resize event whenever a Bootstrap tab
    // becomes visible, since NGL's canvas listens for that.
    $(document).on('shown.bs.tab', 'a[data-toggle=\"tab\"], a[data-bs-toggle=\"tab\"], button[data-bs-toggle=\"tab\"]', function (e) {
       setTimeout(function() { window.dispatchEvent(new Event('resize')); }, 100);
    });

    // 'Jump to site' navigation: scroll so the alignment block containing
    // the clicked site's block is centered in view.
    Shiny.addCustomMessageHandler('scrollToOffset', function(message) {
      var el = document.getElementById(message.elementId);
      if (!el) return;
      var rect = el.getBoundingClientRect();
      var targetY = window.pageYOffset + rect.top + rect.height * message.fraction - 120;
      window.scrollTo({ top: Math.max(targetY, 0), behavior: 'smooth' });
    });

    // Smoother scroll-vs-zoom: plain mouse wheel always scrolls the page
    // normally (scrollZoom is off in config() below), and Ctrl/Cmd+wheel
    // zooms that plot in place instead -- the same convention as Google
    // Maps and most embedded chart tools. That way hovering a plot while
    // scrolling the page never unexpectedly hijacks the wheel; zooming is
    // an explicit, deliberate action (Ctrl+wheel, drag-box, or the
    // toolbar's +/- buttons) instead of an accidental one.
    function setupWheelZoom(id) {
      var gd = document.getElementById(id);
      if (!gd) return;
      $(gd).off('wheel.epictope').on('wheel.epictope', function(e) {
        if (!(e.ctrlKey || e.metaKey)) return; // let the page scroll normally
        if (!gd._fullLayout || !gd._fullLayout.xaxis || !gd._fullLayout.xaxis.range) return;
        e.preventDefault();
        var range = gd._fullLayout.xaxis.range;
        var factor = e.originalEvent.deltaY < 0 ? 0.85 : (1 / 0.85); // wheel up = zoom in
        var mid = (range[0] + range[1]) / 2;
        var halfWidth = (range[1] - range[0]) / 2 * factor;
        Plotly.relayout(gd, { 'xaxis.range': [mid - halfWidth, mid + halfWidth] });
      });
    }
    $(document).on('shiny:value', function(event) {
      // msa_view intentionally excluded: it's a fixed-width, drag-to-scroll
      // strip, not a zoomable plot (see the 'not zoomable' note on that panel).
      if (event.name === 'scores_view') {
        setTimeout(function() { setupWheelZoom(event.target.id); }, 50);
      }
    });
    "
  ))),
  titlePanel("EpicTope"),
  sidebarLayout(
    sidebarPanel(
      tooltip(
        selectizeInput(
          "query", "Protein (zebrafish)",
          choices = NULL,
          options = list(
            placeholder = "Search by gene name, protein name, or UniProt ID...",
            create = TRUE,        # still allow typing a raw UniProt ID not in the list
            maxOptions = 2000,
            render = I(
              "{ option: function(item, escape) {
                   return '<div title=\"' + escape(item.label) + '\">' + escape(item.label) + '</div>';
                 } }"
            )
          )
        ),
        "Start typing a gene name (e.g. smad5, hdac1) or paste a UniProt ID directly. Hover any option to see its full protein name."
      ),
      tooltip(
        fileInput("custom_structure", "Custom AlphaFold .cif (optional)", accept = ".cif"),
        "Upload your own AlphaFold structure file instead of fetching one automatically from UniProt/AlphaFoldDB."
      ),
      tooltip(
        numericInput("start_resnum", "Start residue index (if custom structure)", value = 1, min = 1),
        "If your custom structure doesn't start numbering at residue 1, set the true starting residue number here."
      ),
      tooltip(
        actionButton("run", "Run EpicTope"),
        "Runs the full pipeline: BLAST homolog search, MUSCLE alignment (Shannon entropy), DSSP secondary structure/RSA, IUPred2A/ANCHOR2 disorder, and combined scoring."
      ),
      tooltip(
        numericInput("n_extrema", "Top peaks to highlight", value = 5, min = 1, max = 15, step = 1),
        "How many top-scoring local maxima (good tag candidates, green) to mark on every panel, in addition to the N-term/C-term. Updates instantly without re-running the pipeline."
      ),
      tooltip(
        downloadButton("download_csv", "Download results CSV"),
        "Download just the full per-residue feature and score table."
      ),
      tooltip(
        downloadButton("download_report", "Download full report (.zip)"),
        "Download everything: the score table, the alignment used (FASTA), the AlphaFold structure file, a static min-score plot, and a summary HTML report listing the top candidate tag sites."
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          tooltip(span("Score table ⓘ"), "Full per-residue table: position, normalized entropy/secondary-structure/RSA/disorder scores, combined minimum score, DSSP and IUPred2A raw values."),
          tableOutput("score_table")
        ),
        tabPanel(
          tooltip(
            span("Combined view ⓘ"),
            paste(
              "3D structure, min-score/feature plots, and the raw alignment, all on one page and linked together.",
              "Amber = sequence termini, green = top scoring peaks (good tag candidates), red = curated UniProt binding-site residues (avoid tagging these) -- use the buttons at the top to jump straight to any of them.",
              "Scrolling the page always works normally over any plot; hold Ctrl (Cmd on Mac) and scroll to zoom that plot instead, or drag a box / use the toolbar to zoom or select a range.",
              "A click (or a range selection) on any panel highlights that same residue or range -- in blue, flashing -- clearly on every other panel, including the 3D structure (with a residue-number label)."
            )
          ),
          uiOutput("site_nav_ui"),
          site_legend_ui(),
          tags$hr(),
          tooltip(
            tags$strong("3D structure"),
            "AlphaFold model, cartoon colored by chain position. Spacefill blobs mark the same sites as the plots below (see legend above); a residue-number label is drawn on each. Drag to rotate, scroll to zoom, right-drag to pan."
          ),
          div(style = "position: relative; width: 100%; height: 500px;",
              NGLVieweROutput("structure_view", height = "500px")),
          tags$hr(),
          tooltip(
            tags$strong("Tagging-score & feature plot"),
            paste(
              "Top: 'Min score' -- the minimum of the four normalized feature scores below, averaged over a 7-residue window; ranges 0-1, higher = better candidate for inserting an epitope tag without disrupting the protein.",
              "Bottom: the four underlying features, each normalized 0-1. Entropy: sequence variability across homologs (higher = less conserved = safer). Secondary structure: 1 = loop/coil, 0 = helix/sheet (higher = more tolerant of insertion). RSA: relative solvent accessibility (higher = more surface-exposed). Disorder (DBR): inverted ANCHOR2 disordered-binding-region score (higher = less likely to be a protein-binding interface).",
              "A bright magenta vertical line always marks wherever the alignment strip below is centered -- it moves live as you scroll that strip sideways, even before you've clicked anything."
            )
          ),
          tags$div(style = "font-size: 12px; color: #666; margin-bottom: 4px;",
                   "Scroll normally to move the page. Hold Ctrl (⌘ on Mac) + scroll to zoom this plot, drag to box-zoom, or use the toolbar to switch to box/lasso select for tagging a range. Trace names in the plot's own legend (top-right) toggle each line on/off."),
          plotlyOutput("scores_view", height = "570px"),
          tags$hr(),
          tooltip(
            tags$strong("Multiple sequence alignment"),
            "Homologous protein sequences (from the configured species list) aligned to your query. Letter colors group amino acids by biochemical property (see legend below) so conserved biochemical character, not just identical letters, is visible at a glance. '-' = gap."
          ),
          tooltip(
            checkboxInput("msa_full_view", "Show full alignment (disable the scrollable strip)", value = FALSE),
            "By default the alignment below is a fixed 50-residue-wide strip centered on the top candidate peak (or wherever you last clicked/selected) -- drag it sideways to scroll; it's deliberately not zoomable. Check this box to render the entire alignment as stacked wrapped blocks instead."
          ),
          textOutput("msa_window_label"),
          site_legend_ui(),
          aa_legend_ui(),
          tags$div(style = "font-size: 12px; color: #666; margin: 4px 0;",
                   "Default view: click-and-drag the strip to scroll sideways (it won't zoom -- there's no need, it's a fixed 50-residue window). Vertical page scroll still works normally. 'Show full alignment' switches to a tall, wrapped view of everything at once, with a position ruler (residue number every 10 columns) in each block's header row."),
          uiOutput("msa_view_ui")
        )
      )
    )
  )
)

server <- function(input, output, session) {

  updateSelectizeInput(session, "query", choices = protein_choices, server = TRUE, selected = character(0))

  results <- eventReactive(input$run, {
    withProgress(message = "Running EpicTope pipeline...", {
      setup_files()
      check_config()

      uniprot_fields <- c("accession", "id", "gene_names", "xref_alphafolddb",
                           "sequence", "organism_name", "organism_id", "ft_binding")
      uniprot_data <- query_uniProt(query = input$query, fields = uniprot_fields)
      binding_sites <- parse_uniprot_binding(uniprot_data$Binding.site)

      alphafold_file <- if (!is.null(input$custom_structure)) {
        input$custom_structure$datapath
      } else {
        af_id <- uniprot_data$AlphaFoldDB
        if (is.na(af_id)) af_id <- input$query
        fetch_alphafold(gsub(";", "", af_id))
      }

      dssp_res <- dssp_command(alphafold_file)
      dssp_df  <- parse_dssp(dssp_res)
      iupred_df <- iupredAnchor(input$query)

      seq <- Biostrings::AAStringSet(uniprot_data$Sequence, use.names = TRUE)
      aa_files <- list.files(cds_folder, pattern = "\\.all\\.fa$",
                              ignore.case = TRUE, full.names = TRUE, recursive = TRUE)
      names(aa_files) <- aa_files
      blast_results <- lapply(aa_files, function(.x) protein_blast(seq, .x))
      blast_best_match <- lapply(blast_results, function(.x) head(.x[order(.x$E), ], 1))
      blast_seqs <- lapply(blast_best_match, fetch_sequences)
      blast_seqs[[input$query]] <- seq
      blast_stringset <- Biostrings::AAStringSet(unlist(lapply(blast_seqs, function(.x) .x[[1]])))

      msa_res <- tryCatch(
        Biostrings::AAStringSet(muscle(blast_stringset)),
        error = function(e) stop("Alignment failed: ", conditionMessage(e))
      )

      shannon_df <- shannon_reshape(msa_res, input$query)

      features_df <- Reduce(function(x, y) merge(x, y, all = TRUE),
                             list(shannon_df, dssp_df, iupred_df))
      norm_feats_df <- calculate_scores(features_df)
      final_df <- merge(norm_feats_df, features_df)

      # Everything the various panels need, computed once per run. Peak
      # count depends on a live UI input (n_extrema), so that's computed
      # separately in the `tag_sites` reactive below instead of baked in
      # here -- that lets the user tweak the count without re-running the
      # whole (slow) pipeline.
      list(
        final_df = final_df,
        msa_res = msa_res,
        alphafold_file = alphafold_file,
        binding_sites = binding_sites,
        query_id = input$query
      )
    })
  })

  # Default candidate sites (termini + top peaks), recomputed live whenever
  # the pipeline reruns or the user changes "Top peaks to highlight" --
  # cheap to recompute, so no need to gate it behind `run`.
  tag_sites <- reactive({
    res <- results(); req(res)
    n_extrema <- input$n_extrema
    if (is.null(n_extrema) || is.na(n_extrema) || n_extrema < 1) n_extrema <- 5
    pick_tag_sites(res$final_df, n_extrema = n_extrema, min_sep = 15)
  })

  # --- Cross-panel tagging state -----------------------------------------
  # A single shared value: NULL (nothing user-tagged yet) or c(start, end)
  # (start == end for a plain click). Every panel (score/feature plot,
  # alignment, 3D structure) reads this and redraws its own highlight.
  tag_range <- reactiveVal(NULL)
  observeEvent(results(), { tag_range(NULL) })

  # Where the user is currently scrolled to in the (non-full) alignment
  # strip -- updated as they drag it sideways (see the plotly_relayout
  # observer below). Deliberately NOT read by the alignment panel itself
  # (that would fight the user's own drag with a server round-trip on every
  # scroll tick); it only drives the "current position of interest" line on
  # the score/feature plot, so that line tracks wherever they've scrolled.
  pan_center <- reactiveVal(NULL)
  observeEvent(results(), { pan_center(NULL) })
  observeEvent(event_data("plotly_relayout", source = "msa_view"), {
    ev <- event_data("plotly_relayout", source = "msa_view")
    r0 <- ev[["xaxis.range[0]"]]; r1 <- ev[["xaxis.range[1]"]]
    if (!is.null(r0) && !is.null(r1)) pan_center(mean(c(as.numeric(r0), as.numeric(r1))))
  })

  # All ungapped query-sequence columns in the full alignment (i.e. every
  # residue position 1..N).
  full_query_cols <- reactive({
    res <- results(); req(res)
    mat_full <- as.matrix(res$msa_res)
    which(mat_full[res$query_id, ] != "-")
  })

  # The single "current position of interest" -- drives the dynamic line on
  # the score/feature plot AND the alignment strip's default scroll
  # position. Priority: an explicit click/range-selection anywhere, else
  # wherever the user has scrolled the alignment to, else the single best
  # interior peak (the global maximum that isn't on a terminus).
  focus_position <- reactive({
    rng <- tag_range()
    if (!is.null(rng)) return(mean(rng))
    pc <- pan_center()
    if (!is.null(pc)) return(pc)
    ts <- tryCatch(tag_sites(), error = function(e) NULL)
    if (is.null(ts)) return(NULL)
    if (nrow(ts$peaks) > 0) ts$peaks$position[1] else ts$n_term
  })

  # The fixed-width [lo, hi] view the (non-full) alignment strip opens to.
  # Always exactly `2 * MSA_DEFAULT_HALF_WIDTH` residues wide -- the strip
  # is deliberately not zoomable, so this only ever changes by re-centering,
  # never by resizing (see the scrollable single-row renderer below).
  msa_window <- reactive({
    full_cols <- full_query_cols()
    n_full <- length(full_cols)
    if (isTRUE(input$msa_full_view) || n_full <= 2 * MSA_DEFAULT_HALF_WIDTH) {
      return(c(1, n_full))
    }
    center <- focus_position()
    if (is.null(center)) center <- 1 + MSA_DEFAULT_HALF_WIDTH
    lo <- round(center - MSA_DEFAULT_HALF_WIDTH); hi <- round(center + MSA_DEFAULT_HALF_WIDTH)
    if (lo < 1) { hi <- hi + (1 - lo); lo <- 1 }
    if (hi > n_full) { lo <- lo - (hi - n_full); hi <- n_full }
    c(max(1, lo), min(n_full, hi))
  })

  output$msa_window_label <- renderText({
    res <- tryCatch(results(), error = function(e) NULL)
    if (is.null(res)) return("")
    w <- msa_window(); n <- length(full_query_cols())
    if (isTRUE(input$msa_full_view) || (w[1] == 1 && w[2] == n)) {
      sprintf("Showing the full alignment: residues 1-%d.", n)
    } else {
      sprintf("Showing residues %d-%d of %d -- drag the strip sideways to scroll (check the box above for the full alignment instead).", w[1], w[2], n)
    }
  })

  register_tag_events <- function(source_name) {
    click <- event_data("plotly_click", source = source_name)
    if (!is.null(click) && !is.null(click$customdata)) {
      pos <- round(click$customdata[1])
      tag_range(c(pos, pos))
    }
  }
  register_range_events <- function(source_name) {
    sel <- event_data("plotly_selected", source = source_name)
    if (!is.null(sel) && !is.null(sel$customdata) && length(sel$customdata) > 0) {
      vals <- round(as.numeric(sel$customdata))
      vals <- vals[!is.na(vals)]
      if (length(vals) > 0) tag_range(c(min(vals), max(vals)))
    }
  }
  observeEvent(event_data("plotly_click", source = "scores_view"), register_tag_events("scores_view"))
  observeEvent(event_data("plotly_selected", source = "scores_view"), register_range_events("scores_view"))
  observeEvent(event_data("plotly_click", source = "msa_view"), register_tag_events("msa_view"))
  observeEvent(event_data("plotly_selected", source = "msa_view"), register_range_events("msa_view"))

  output$score_table <- renderTable(head(results()$final_df, 20))

  # --- Alignment layout (shared by the MSA panel and the nav "jump to") ---
  # Always the full, ungapped query-aligned matrix -- "show full alignment"
  # and the default scrollable strip both read from this and just choose a
  # different way to lay the same data out (see `output$msa_view`).
  msa_layout <- reactive({
    res <- results()
    mat_full <- as.matrix(res$msa_res)
    query_cols <- full_query_cols()
    mat <- mat_full[, query_cols, drop = FALSE]
    seq_names <- sub("\\..*$", "", sub(".*/", "", rownames(mat)))
    list(mat = mat, seq_names = seq_names, n_positions = length(query_cols), n_seqs = nrow(mat))
  })

  # --- Top-of-tab nav bar: jump straight to any candidate site ------------
  output$site_nav_ui <- renderUI({
    res <- tryCatch(results(), error = function(e) NULL)
    if (is.null(res)) return(tags$em("Run EpicTope to see candidate tag sites here."))
    ts <- tag_sites()
    sites <- site_list_for(ts)
    kind_color <- c(term = COLOR_DEFAULT_SITE, peak = COLOR_PEAK)
    has_binding <- length(res$binding_sites) > 0
    binding_box <- tags$div(
      style = paste0("margin-top: 6px; padding: 6px 10px; border-radius: 4px; font-size: 12px;",
                      if (has_binding) paste0(" background:", COLOR_BINDING, "22; border: 1px solid ", COLOR_BINDING, ";")
                      else " background:#f0f0f0; border: 1px solid #ccc; color:#555;"),
      if (has_binding) {
        tagList(tags$strong("UniProt binding-site residues found (shown in red everywhere): "),
                format_binding_summary(res$binding_sites), ".")
      } else {
        "No curated UniProt binding-site (ft_binding) annotations exist for this protein -- this is a real UniProt lookup, not every protein has one curated."
      }
    )
    tagList(
      tags$strong("Candidate / notable sites: "),
      tagList(lapply(seq_len(nrow(sites)), function(i) {
        s <- sites[i, ]
        make_jump_button(
          sprintf("%s — pos %d (%.3f)", s$label, s$position, s$score),
          s$position, kind_color[[s$kind]],
          sprintf("Jump to %s (residue %d): highlights it in blue everywhere and centers the alignment window on it.",
                  s$label, s$position)
        )
      })),
      binding_box,
      tags$div(style = "margin-top: 4px; font-style: italic; color: #444;",
               textOutput("tag_banner", inline = TRUE))
    )
  })

  output$tag_banner <- renderText(tag_banner_text(tag_range()))

  jump_to_position <- function(pos) {
    res <- results(); req(res)
    tag_range(c(pos, pos))

    n_positions <- nrow(res$final_df)
    half_window <- max(20, round(n_positions * 0.05))
    plotlyProxy("scores_view", session) %>%
      plotlyProxyInvoke("relayout", list(xaxis.range = list(pos - half_window, pos + half_window)))

    # The alignment re-windows itself around the new tag_range() on its own
    # (see `msa_window`); just scroll the page to bring it into view.
    session$sendCustomMessage("scrollToOffset", list(elementId = "msa_view", fraction = 0.05))
  }
  observeEvent(input$jump_pos, jump_to_position(input$jump_pos))

  # --- Panel: min-score + feature tracks (Fig 6C / 6B), zoomable ----------
  output$scores_view <- renderPlotly({
    res <- results()
    df <- res$final_df[order(res$final_df$position), ]
    ts <- tag_sites()
    smoothed <- ts$smoothed

    sites <- site_list_for(ts)
    term_sites <- sites[sites$kind == "term", ]
    peak_sites <- sites[sites$kind == "peak", ]

    add_kind_markers <- function(p, sdf, color, symbol) {
      if (nrow(sdf) == 0) return(p)
      p %>% add_markers(data = sdf, x = ~position, y = ~score, showlegend = FALSE,
                         marker = list(color = color, size = 14, symbol = symbol,
                                       line = list(color = "#333", width = 2)),
                         customdata = ~position,
                         text = ~paste0(label, ": ", sprintf("%.3f", score)),
                         hovertemplate = "%{text}<extra></extra>")
    }

    # mode = "lines+markers" with invisible markers so the toolbar's box/lasso
    # select tool has actual points to hit-test against (a pure "lines" trace
    # isn't selectable in plotly.js) -- the line still *looks* like a plain line.
    p1 <- plot_ly(source = "scores_view") %>%
      add_trace(data = df, x = ~position, y = smoothed, type = "scatter", mode = "lines+markers",
                name = "Min score", line = list(color = "black"),
                marker = list(size = 8, opacity = 0), customdata = ~position,
                hovertemplate = paste0(
                  "Position %{x}<br>Min score: %{y:.3f}",
                  "<br><i>Minimum of entropy/structure/RSA/disorder (0-1);<br>higher = better tag-insertion candidate</i>",
                  "<extra></extra>")) %>%
      add_kind_markers(term_sites, COLOR_DEFAULT_SITE, "diamond") %>%
      add_kind_markers(peak_sites, COLOR_PEAK, "triangle-up") %>%
      layout(yaxis = list(title = "Min score (7-res avg)"))

    p2 <- plot_ly(source = "scores_view")
    feat_specs <- list(
      list(col = "normalized_entropy", name = "Entropy", color = "black",
           desc = "Sequence variability across homologs (0=conserved, 1=variable); higher = safer to tag"),
      list(col = "ss_score", name = "Secondary structure", color = "#8B4513",
           desc = "0=helix/sheet, 1=loop/coil; higher = more tolerant of an insertion"),
      list(col = "rsa", name = "RSA", color = "forestgreen",
           desc = "Relative solvent accessibility (0=buried, 1=exposed); higher = more surface-exposed"),
      list(col = "inv_anchor2", name = "Disorder (DBR)", color = "#008080",
           desc = "Inverted ANCHOR2 binding-region score; higher = less likely a protein-binding interface")
    )
    for (spec in feat_specs) {
      p2 <- p2 %>% add_trace(
        data = df, x = ~position, y = df[[spec$col]], type = "scatter", mode = "lines+markers",
        name = spec$name, line = list(color = spec$color),
        marker = list(size = 8, opacity = 0), customdata = ~position,
        hovertemplate = paste0("Position %{x}<br>", spec$name, ": %{y:.3f}<br><i>", spec$desc, "</i><extra></extra>")
      )
    }
    p2 <- p2 %>% layout(yaxis = list(title = "Normalized score", range = c(0, 1)))

    # Default "regions of interest" -- always shown as a shaded band (not
    # just a hairline) plus a bold dashed center line, so they read as
    # regions at a glance, not something you have to spot a thin line for.
    # One canonical `color` per kind for both the band and the line/text --
    # a separate darker "line color" per kind used to exist here, which just
    # meant the same site looked like two different colors within one plot.
    build_site_shapes <- function(sdf, color) {
      if (nrow(sdf) == 0) return(list(bands = list(), lines = list(), annotations = list()))
      bands <- lapply(sdf$position, function(pos) {
        list(type = "rect", x0 = pos - 1.5, x1 = pos + 1.5, y0 = 0, y1 = 1, xref = "x",
             fillcolor = color, opacity = 0.30, line = list(width = 0))
      })
      lines <- lapply(sdf$position, function(pos) {
        list(type = "line", x0 = pos, x1 = pos, y0 = 0, y1 = 1, xref = "x",
             line = list(color = color, dash = "dash", width = 2))
      })
      annotations <- lapply(seq_len(nrow(sdf)), function(i) {
        list(x = sdf$position[i], y = sdf$score[i], xref = "x", yref = "y",
             text = sprintf("%.3f", sdf$score[i]), showarrow = TRUE, arrowhead = 0,
             ax = 0, ay = -25, font = list(size = 11, color = color, family = "Arial Black"))
      })
      list(bands = bands, lines = lines, annotations = annotations)
    }
    term_vis <- build_site_shapes(term_sites, COLOR_DEFAULT_SITE)
    peak_vis <- build_site_shapes(peak_sites, COLOR_PEAK)

    # Curated UniProt binding-site residues -- drawn as wide, low-opacity
    # red bands (behind everything else) so tagging near a known functional
    # site is obvious without needing per-residue diamonds/labels.
    binding_runs <- positions_to_runs(res$binding_sites)
    binding_shapes <- if (nrow(binding_runs) > 0) {
      lapply(seq_len(nrow(binding_runs)), function(i) {
        list(type = "rect", x0 = binding_runs$start[i] - 0.5, x1 = binding_runs$end[i] + 0.5,
             y0 = 0, y1 = 1, xref = "x", fillcolor = COLOR_BINDING, opacity = 0.18, line = list(width = 0))
      })
    } else list()

    # Whatever the user last clicked/selected -- on this panel or any other --
    # drawn last (on top) with a solid outline so it's unmistakable even when
    # it overlaps a default site.
    user_shapes <- if (!is.null(tag_range())) {
      rng <- tag_range()
      list(list(type = "rect", x0 = rng[1] - 1, x1 = rng[2] + 1, y0 = 0, y1 = 1, xref = "x",
                fillcolor = COLOR_USER_TAG, opacity = 0.35,
                line = list(color = COLOR_USER_TAG, width = 3)))
    } else list()

    # "Current position of interest" -- a thin line that always tracks
    # wherever the alignment strip below is centered on: a click/selection,
    # wherever the user has scrolled that strip to, or (with nothing picked
    # yet) the default top peak. Distinct from the site markers/bands above
    # (which mark fixed candidate sites) and from the user-tag box (which
    # only appears after a click) -- this one is always present and moves
    # live as the alignment strip is scrolled.
    # Bright magenta, on top of everything else (including the user-tag
    # box) -- deliberately a color used nowhere else on this plot (the
    # "Min score" and "Entropy" traces are both black, so black was
    # invisible here -- it blended straight into those curves).
    focus_pos <- focus_position()
    focus_shapes <- if (!is.null(focus_pos)) {
      list(list(type = "line", x0 = focus_pos, x1 = focus_pos, y0 = 0, y1 = 1, xref = "x",
                line = list(color = "#FF00FF", width = 3)))
    } else list()

    all_shapes <- c(binding_shapes, term_vis$bands, peak_vis$bands,
                     term_vis$lines, peak_vis$lines, user_shapes, focus_shapes)
    all_annotations <- c(term_vis$annotations, peak_vis$annotations)
    banner <- list(list(x = 0, y = 1.16, xref = "paper", yref = "paper", xanchor = "left",
                         showarrow = FALSE, font = list(size = 12, color = COLOR_USER_TAG),
                         text = tag_banner_text(tag_range())))

    combined <- subplot(p1, p2, nrows = 2, shareX = TRUE, titleY = TRUE, heights = c(0.4, 0.6))
    combined$x$source <- "scores_view"
    combined %>%
      layout(
        title = paste("EpicTope analysis:", res$query_id),
        dragmode = "zoom",   # drag to zoom; use the toolbar to switch to box/lasso select for tagging a range
        legend = list(orientation = "h", y = 1.1),
        margin = list(t = 90),
        # shapes/annotations need a yref per row (y for row1, y2 for row2) --
        # duplicate each shape onto both rows so it visibly spans the full view.
        shapes = c(
          lapply(all_shapes, function(s) { s$yref <- "y"; s }),
          lapply(all_shapes, function(s) { s$yref <- "y2"; s })
        ),
        annotations = c(lapply(all_annotations, function(a) { a$yref <- "y"; a }), banner)
      ) %>%
      event_register("plotly_click") %>%
      event_register("plotly_selected") %>%
      config(scrollZoom = FALSE, displaylogo = FALSE)  # wheel scrolls the page; Ctrl/Cmd+wheel zooms (see setupWheelZoom JS)
  })

  # --- Panel: raw, colored alignment (Fig 2 / Fig 6A style) ---------------
  # Two render modes share the same colored-cell + grouped-run-outline
  # scheme (built once below): "show full alignment" stacks the whole thing
  # into fixed-width wrapped blocks (unchanged from before); the default is
  # a fixed 50-residue-wide strip on a real position axis that the user
  # drags sideways to scroll -- deliberately not zoomable.
  output$msa_view_ui <- renderUI({
    req(results())
    ml <- msa_layout()
    if (isTRUE(input$msa_full_view)) {
      n_blocks <- ceiling(ml$n_positions / MSA_WRAP_WIDTH)
      total_rows <- n_blocks * (ml$n_seqs + 1)  # +1 header/ruler row per block
      msa_height <- max(280, total_rows * 20 + 70)
    } else {
      msa_height <- max(200, ml$n_seqs * 24 + 90)
    }
    plotlyOutput("msa_view", height = paste0(msa_height, "px"), width = "100%")
  })

  output$msa_view <- renderPlotly({
    res <- results()
    ml <- msa_layout()
    mat <- ml$mat; seq_names <- ml$seq_names
    n_positions <- ml$n_positions; n_seqs <- ml$n_seqs

    # Cell fill is ALWAYS the true amino-acid color -- never overridden by a
    # site highlight. Two color scales (termini/peaks/binding-site/gap and
    # amino-acid property) sharing the same cells would inevitably collide
    # somewhere (e.g. the peak green used to be the exact same hex as the
    # "polar" amino-acid color, so a highlighted peak was indistinguishable
    # from an ordinary polar residue elsewhere). Sites are marked purely by
    # colored outline boxes (below), which never touch the letter/fill color.
    letters_present <- sort(unique(as.vector(mat)))
    n_base <- length(letters_present)
    letter_idx <- setNames(seq_along(letters_present) - 1, letters_present)
    n_total_bins <- n_base

    bin_colors <- vapply(letters_present, function(l) {
      c <- unname(aa_colors[l]); if (is.na(c)) "#BBBBBB" else c
    }, character(1))
    msa_colorscale <- list()
    for (i in seq_along(bin_colors)) {
      msa_colorscale[[length(msa_colorscale) + 1]] <- list((i - 1) / n_total_bins, bin_colors[i])
      msa_colorscale[[length(msa_colorscale) + 1]] <- list(i / n_total_bins, bin_colors[i])
    }

    ts <- tag_sites()
    sites <- site_list_for(ts)
    term_positions    <- sites$position[sites$kind == "term"]
    peak_positions    <- sites$position[sites$kind == "peak"]
    binding_positions <- res$binding_sites
    user_positions <- if (!is.null(tag_range())) {
      rng <- tag_range(); seq(rng[1], rng[2])
    } else integer(0)

    run_color <- c(term = COLOR_DEFAULT_SITE, peak = COLOR_PEAK, binding = COLOR_BINDING, user = COLOR_USER_TAG)
    run_width <- c(term = 3, peak = 3, binding = 3, user = 4)
    # Contiguous highlighted runs -- one entry per (row, kind, run), so an
    # entire highlighted stretch (a binding motif, a dragged range) draws as
    # ONE outline box spanning the whole run instead of one square per residue.
    runs <- list()
    add_runs <- function(row_key, hit_local_idx, kind) {
      for (r in contiguous_runs(hit_local_idx)) {
        runs[[length(runs) + 1]] <<- list(row_key = row_key, start = r["start"], end = r["end"], kind = kind)
      }
    }
    build_run_shapes <- function(row_index) {
      shapes <- lapply(runs, function(r) {
        yi <- row_index[[r$row_key]]
        list(type = "rect", xref = "x", yref = "y",
             x0 = r$start - 0.5, x1 = r$end + 0.5, y0 = yi - 0.45, y1 = yi + 0.45,
             line = list(color = run_color[[r$kind]], width = run_width[[r$kind]]),
             fillcolor = "rgba(0,0,0,0)")
      })
      # User-selection boxes drawn last (on top) so they're never hidden
      # behind a default site's outline when they overlap.
      shapes[order(vapply(runs, function(r) r$kind == "user", logical(1)))]
    }
    hover_suffix <- paste0(
      "<br><i>Fill color = amino-acid property; an outlined box = a highlighted site (see legend below)</i>",
      "<extra></extra>")

    if (isTRUE(input$msa_full_view)) {
      # --- Full alignment: stacked MSA_WRAP_WIDTH-wide blocks (unchanged) ---
      n_blocks <- ceiling(n_positions / MSA_WRAP_WIDTH)
      rows_key <- character(0); rows_label <- character(0)
      z_list <- list(); text_list <- list(); customdata_list <- list()

      for (b in seq_len(n_blocks)) {
        start_col <- (b - 1) * MSA_WRAP_WIDTH + 1
        end_col   <- min(b * MSA_WRAP_WIDTH, n_positions)
        block_len <- end_col - start_col + 1
        block_abs <- start_col:end_col

        # Ruler header: residue number every 10 positions, plus the block's
        # very first column -- so a position is never more than ~10 columns
        # from a visible number (classic Clustal/EMBOSS ruler convention).
        header_text <- rep("", MSA_WRAP_WIDTH)
        tick_local <- which(block_abs %% 10 == 0)
        if (length(tick_local) == 0 || tick_local[1] != 1) tick_local <- c(1, tick_local)
        header_text[tick_local] <- as.character(block_abs[tick_local])

        rows_key   <- c(rows_key, paste0("__hdr", b))
        rows_label <- c(rows_label, "")
        z_list[[length(z_list) + 1]] <- rep(NA_real_, MSA_WRAP_WIDTH)
        text_list[[length(text_list) + 1]] <- header_text
        customdata_list[[length(customdata_list) + 1]] <- rep(NA_real_, MSA_WRAP_WIDTH)

        for (i in seq_len(n_seqs)) {
          row_key <- paste0(seq_names[i], "___b", b)
          rows_key   <- c(rows_key, row_key)
          rows_label <- c(rows_label, seq_names[i])

          row_letters <- rep("", MSA_WRAP_WIDTH)
          row_z <- rep(NA_real_, MSA_WRAP_WIDTH)
          row_customdata <- rep(NA_real_, MSA_WRAP_WIDTH)
          idx <- seq_len(block_len)
          cols <- start_col:end_col
          abs_cols <- block_abs[idx]
          row_letters[idx] <- mat[i, cols]
          row_z[idx] <- letter_idx[mat[i, cols]]
          row_customdata[idx] <- abs_cols

          hit_term <- which(abs_cols %in% term_positions)
          if (length(hit_term) > 0) add_runs(row_key, idx[hit_term], "term")
          hit_peak <- which(abs_cols %in% peak_positions)
          if (length(hit_peak) > 0) add_runs(row_key, idx[hit_peak], "peak")
          hit_binding <- which(abs_cols %in% binding_positions)
          if (length(hit_binding) > 0) add_runs(row_key, idx[hit_binding], "binding")
          hit_user <- which(abs_cols %in% user_positions)
          if (length(hit_user) > 0) add_runs(row_key, idx[hit_user], "user")

          z_list[[length(z_list) + 1]] <- row_z
          text_list[[length(text_list) + 1]] <- row_letters
          customdata_list[[length(customdata_list) + 1]] <- row_customdata
        }
      }

      z_mat <- do.call(rbind, z_list)
      text_mat <- do.call(rbind, text_list)
      customdata_mat <- do.call(rbind, customdata_list)
      row_index <- setNames(seq_along(rows_key) - 1, rows_key)

      plot_ly(source = "msa_view") %>%
        add_trace(
          x = seq_len(MSA_WRAP_WIDTH), y = rows_key, z = z_mat, text = text_mat,
          customdata = customdata_mat,
          texttemplate = "%{text}", textfont = list(size = 11, family = "monospace", color = "black"),
          type = "heatmap", showscale = FALSE, zmin = 0, zmax = n_total_bins,
          colorscale = msa_colorscale, xgap = 1, ygap = 1,
          hovertemplate = paste0("Position %{customdata}<br>Residue %{text}", hover_suffix)
        ) %>%
        layout(
          xaxis = list(title = "", showticklabels = FALSE),
          yaxis = list(title = "", autorange = "reversed", tickmode = "array",
                       tickvals = rows_key, ticktext = rows_label),
          dragmode = "zoom",
          margin = list(t = 40),
          shapes = build_run_shapes(row_index),
          annotations = list(list(x = 0, y = 1.05, xref = "paper", yref = "paper", xanchor = "left",
                                   showarrow = FALSE, font = list(size = 12, color = COLOR_USER_TAG),
                                   text = tag_banner_text(tag_range())))
        ) %>%
        event_register("plotly_click") %>%
        event_register("plotly_selected") %>%
        config(scrollZoom = FALSE, displaylogo = FALSE)  # wheel scrolls the page; Ctrl/Cmd+wheel zooms (see setupWheelZoom JS)

    } else {
      # --- Default: one continuous row per sequence, real position axis,
      # fixed-width and drag-to-scroll -- intentionally not zoomable, so
      # scrolling sideways (not zooming out) is the only way to see more. ---
      rows_key <- paste0(seq_names, "___r", seq_len(n_seqs))
      rows_label <- seq_names
      z_mat <- matrix(letter_idx[mat], nrow = n_seqs, ncol = n_positions)
      text_mat <- mat
      customdata_mat <- matrix(rep(seq_len(n_positions), each = n_seqs), nrow = n_seqs, ncol = n_positions)

      cols_all <- seq_len(n_positions)
      hit_term <- which(cols_all %in% term_positions)
      hit_peak <- which(cols_all %in% peak_positions)
      hit_binding <- which(cols_all %in% binding_positions)
      hit_user <- which(cols_all %in% user_positions)
      for (i in seq_len(n_seqs)) {
        row_key <- rows_key[i]
        if (length(hit_term) > 0)    add_runs(row_key, hit_term, "term")
        if (length(hit_peak) > 0)    add_runs(row_key, hit_peak, "peak")
        if (length(hit_binding) > 0) add_runs(row_key, hit_binding, "binding")
        if (length(hit_user) > 0)    add_runs(row_key, hit_user, "user")
      }
      row_index <- setNames(seq_along(rows_key) - 1, rows_key)
      view_range <- unname(msa_window())

      plot_ly(source = "msa_view") %>%
        add_trace(
          x = cols_all, y = rows_key, z = z_mat, text = text_mat,
          customdata = customdata_mat,
          texttemplate = "%{text}", textfont = list(size = 11, family = "monospace", color = "black"),
          type = "heatmap", showscale = FALSE, zmin = 0, zmax = n_total_bins,
          colorscale = msa_colorscale, xgap = 1, ygap = 1,
          hovertemplate = paste0("Position %{x}<br>Residue %{text}", hover_suffix)
        ) %>%
        layout(
          xaxis = list(title = "Residue position", range = view_range, dtick = 10, fixedrange = FALSE),
          yaxis = list(title = "", autorange = "reversed", tickmode = "array",
                       tickvals = rows_key, ticktext = rows_label, fixedrange = TRUE),
          dragmode = "pan",
          margin = list(t = 40),
          shapes = build_run_shapes(row_index),
          annotations = list(list(x = 0, y = 1.08, xref = "paper", yref = "paper", xanchor = "left",
                                   showarrow = FALSE, font = list(size = 12, color = COLOR_USER_TAG),
                                   text = tag_banner_text(tag_range())))
        ) %>%
        event_register("plotly_click") %>%
        event_register("plotly_selected") %>%
        event_register("plotly_relayout") %>%
        config(scrollZoom = FALSE, displaylogo = FALSE, doubleClick = FALSE,
               modeBarButtonsToRemove = c("zoom2d", "zoomIn2d", "zoomOut2d", "autoScale2d", "resetScale2d"))
    }
  })

  # --- Panel: 3D structure (Fig 6D), same default + user highlighting -----
  output$structure_view <- renderNGLVieweR({
    res <- results()
    req(res$alphafold_file)
    ts <- tag_sites()
    sites <- site_list_for(ts)

    # NGLVieweR reads local files itself (readLines) and hands NGL.js the raw
    # text as a Blob -- so a relative-vs-absolute path was never the real
    # issue (the console's "loading file ''" is just NGL's log line for an
    # unnamed Blob). The crash was inside NGL's parser: an experimental PDB
    # deposition's CIF has biological-assembly records
    # (_pdbx_struct_assembly_gen etc); an AlphaFold *computed model* CIF has
    # none (it's a single predicted chain), and NGL indexed into that
    # missing array. AlphaFoldDB also publishes the same model as a plain
    # .pdb, and NGL's PDB parser doesn't expect assembly records, so we
    # prefer that format when we can derive its URL from the .cif we fetched.
    struct_file   <- res$alphafold_file
    struct_format <- tolower(tools::file_ext(struct_file))

    if (struct_format %in% c("cif", "mmcif", "mcif") &&
        grepl("^AF-", basename(struct_file))) {
      pdb_url   <- sub("\\.(cif|mmcif|mcif)$", ".pdb",
                        file.path("https://alphafold.ebi.ac.uk/files", basename(struct_file)),
                        ignore.case = TRUE)
      pdb_local <- sub("\\.(cif|mmcif|mcif)$", ".pdb", struct_file, ignore.case = TRUE)
      got_pdb <- tryCatch({
        if (!file.exists(pdb_local)) {
          download.file(pdb_url, pdb_local, quiet = TRUE, mode = "wb")
        }
        file.exists(pdb_local) && file.info(pdb_local)$size > 0
      }, error = function(e) FALSE)
      if (isTRUE(got_pdb)) {
        struct_file   <- pdb_local
        struct_format <- "pdb"
      }
    }

    abs_path <- normalizePath(struct_file, mustWork = FALSE)
    cat("structure_view: loading", abs_path, "(format:", struct_format,
        ") exists =", file.exists(abs_path), "\n")

    # Default candidate sites, all highlighted as bold spacefill blobs (much
    # harder to miss than thin ball-and-stick) plus a residue-number label,
    # named so the proxy below can add/replace a separate "userTag"
    # selection independently. Curated UniProt binding-site residues get
    # their own translucent red spacefill so functional sites are obvious
    # right on the structure, not just in the plots.
    kind_color <- c(term = COLOR_DEFAULT_SITE, peak = COLOR_PEAK)
    viewer <- NGLVieweR(abs_path, format = struct_format) %>%
      addRepresentation("cartoon", param = list(colorScheme = "residueindex"))

    for (i in seq_len(nrow(sites))) {
      s <- sites[i, ]
      viewer <- viewer %>%
        addRepresentation("spacefill", param = list(
          name = paste0("site_", s$id), sele = as.character(s$position),
          colorValue = kind_color[[s$kind]], opacity = 0.9)) %>%
        addRepresentation("label", param = list(
          sele = as.character(s$position), labelType = "format", labelFormat = "%(resno)s",
          labelGrouping = "residue", color = kind_color[[s$kind]],
          showBackground = TRUE, backgroundColor = "black", backgroundOpacity = 0.5))
    }

    if (length(res$binding_sites) > 0) {
      viewer <- viewer %>%
        addRepresentation("spacefill", param = list(
          name = "uniprot_binding", sele = paste(res$binding_sites, collapse = " or "),
          colorValue = COLOR_BINDING, opacity = 0.75))
    }

    viewer %>% stageParameters(backgroundColor = "white")
  })

  # Whatever gets clicked/selected on the score, feature, or alignment panel
  # also highlights on the 3D structure, in the same blue used everywhere
  # else -- as a spacefill blob that flashes (opacity pulses) so a selected
  # residue buried inside the cartoon still catches the eye, not just a
  # static, easy-to-miss recoloring.
  flash_on <- reactiveVal(TRUE)
  observe({
    if (is.null(tag_range())) return()
    invalidateLater(450, session)
    isolate(flash_on(!flash_on()))
  })

  observe({
    req(results())
    rng <- tag_range()
    proxy <- NGLVieweR_proxy("structure_view")
    proxy %>% removeSelection("userTag") %>% removeSelection("userTagLabel")
    if (!is.null(rng)) {
      opacity <- if (flash_on()) 0.95 else 0.25  # flash_on() drives the pulse timing
      sele_str <- if (rng[1] == rng[2]) as.character(rng[1]) else paste0(rng[1], "-", rng[2])
      proxy %>%
        addSelection("spacefill", param = list(name = "userTag", sele = sele_str,
                                                colorValue = COLOR_USER_TAG, opacity = opacity)) %>%
        addSelection("label", param = list(name = "userTagLabel", sele = sele_str,
                                            labelType = "format", labelFormat = "%(resno)s",
                                            labelGrouping = "residue", color = COLOR_USER_TAG,
                                            showBackground = TRUE, backgroundColor = "black",
                                            backgroundOpacity = 0.5))
    }
  })

  output$download_csv <- downloadHandler(
    filename = function() paste0(results()$query_id, "_score.csv"),
    content = function(file) write.csv(results()$final_df, file, row.names = FALSE)
  )

  # --- Full report: CSV + FASTA alignment + structure + plot + HTML summary ---
  output$download_report <- downloadHandler(
    filename = function() paste0(results()$query_id, "_epictope_report.zip"),
    content = function(file) {
      res <- results()
      ts <- tag_sites()
      df <- res$final_df[order(res$final_df$position), ]
      sites <- site_list_for(ts)
      kind_color <- c(term = COLOR_DEFAULT_SITE, peak = COLOR_PEAK)
      site_colors <- kind_color[sites$kind]

      tmp_dir <- tempfile("epictope_report_")
      dir.create(tmp_dir)

      # 1. Full per-residue score table
      write.csv(res$final_df, file.path(tmp_dir, paste0(res$query_id, "_scores.csv")), row.names = FALSE)

      # 2. Alignment used, as FASTA
      Biostrings::writeXStringSet(res$msa_res, file.path(tmp_dir, paste0(res$query_id, "_alignment.fasta")))

      # 3. Structure file used for DSSP/RSA and the 3D view
      struct_ext <- tools::file_ext(res$alphafold_file)
      struct_name <- paste0(res$query_id, "_structure.", struct_ext)
      file.copy(res$alphafold_file, file.path(tmp_dir, struct_name), overwrite = TRUE)

      # 4. Static min-score plot (base R graphics only -- no extra dependency)
      png(file.path(tmp_dir, "min_score_plot.png"), width = 1000, height = 400)
      plot(df$position, ts$smoothed, type = "l", xlab = "Residue position",
           ylab = "Min score (7-res avg)", main = paste("EpicTope min-score profile:", res$query_id))
      if (length(res$binding_sites) > 0) {
        binding_runs <- positions_to_runs(res$binding_sites)
        rect(binding_runs$start - 0.5, par("usr")[3], binding_runs$end + 0.5, par("usr")[4],
             col = grDevices::adjustcolor(COLOR_BINDING, alpha.f = 0.15), border = NA)
      }
      abline(v = sites$position, col = site_colors, lty = 2)
      points(sites$position, sites$score, pch = 18, col = site_colors, cex = 1.6)
      text(sites$position, sites$score, labels = sites$label, pos = 3, col = site_colors, cex = 0.8)
      dev.off()

      # 5. Self-contained HTML summary report (relative image path -> works
      #    once the zip is extracted, no bundling/knitting dependency needed)
      kind_label <- c(term = "Terminus", peak = "Peak (candidate)")
      tag_rows <- paste0(
        "<tr><td>", sites$label, "</td><td>", kind_label[sites$kind], "</td><td>", sites$position, "</td><td>",
        sprintf("%.3f", sites$score), "</td></tr>", collapse = "\n"
      )
      binding_summary <- format_binding_summary(res$binding_sites)
      n_homologs <- nrow(as.matrix(res$msa_res)) - 1
      html <- sprintf(
        '<!DOCTYPE html><html><head><meta charset="utf-8">
<title>EpicTope report: %s</title>
<style>
body { font-family: -apple-system, Arial, sans-serif; margin: 2em; color: #222; }
h1 { margin-bottom: 0; } .meta { color: #666; margin-top: 0.2em; }
table { border-collapse: collapse; margin: 1em 0; }
th, td { border: 1px solid #ccc; padding: 6px 12px; text-align: left; }
th { background: #f4f4f4; }
img { max-width: 100%%; border: 1px solid #ddd; margin: 1em 0; }
</style></head><body>
<h1>EpicTope analysis report</h1>
<p class="meta">Query: <b>%s</b> &middot; Generated %s</p>

<h2>Candidate / notable sites</h2>
<table><tr><th>Site</th><th>Type</th><th>Residue position</th><th>Min score (7-res avg)</th></tr>
%s
</table>

<h2>UniProt binding-site residues</h2>
<p>%s</p>

<h2>Min-score profile</h2>
<img src="min_score_plot.png" alt="Min-score profile">

<h2>Run summary</h2>
<ul>
<li>Protein length analyzed: %d residues</li>
<li>Homologous sequences used in the alignment: %d</li>
<li>Structure file: %s (included in this archive)</li>
</ul>

<h2>Included files</h2>
<ul>
<li><b>%s_scores.csv</b> — full per-residue feature/score table</li>
<li><b>%s_alignment.fasta</b> — the multiple sequence alignment used for the entropy calculation</li>
<li><b>%s</b> — the AlphaFold structure used for DSSP/RSA and the 3D view</li>
<li><b>min_score_plot.png</b> — static copy of the min-score profile above</li>
</ul>
</body></html>',
        res$query_id, res$query_id, format(Sys.time(), "%Y-%m-%d %H:%M"),
        tag_rows, binding_summary, nrow(df), n_homologs, struct_name,
        res$query_id, res$query_id, struct_name
      )
      writeLines(html, file.path(tmp_dir, paste0(res$query_id, "_report.html")))

      # Bundle everything, using the "zip" R package if available (pure R,
      # no external zip binary needed -- safer given this box's history of
      # missing/mismatched native tools); falls back to base R + system zip.
      if (requireNamespace("zip", quietly = TRUE)) {
        zip::zip(zipfile = file, files = list.files(tmp_dir), root = tmp_dir)
      } else {
        old_wd <- setwd(tmp_dir)
        on.exit(setwd(old_wd), add = TRUE)
        utils::zip(zipfile = file, files = list.files("."))
      }
    },
    contentType = "application/zip"
  )
}

shinyApp(ui, server)
