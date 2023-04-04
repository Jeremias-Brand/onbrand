#' Plot Genome Alignment
#'
#' This function creates a visual representation of genome alignments using ggplot2.
#' It allows for customization of the plot, including filtering by mapping quality, displaying gaps, and setting axis limits.
#'
#' @param q_name A character string representing the query name.
#' @param r_name A character string representing the reference name.
#' @param fq A vector of query sequences.
#' @param fr A vector of reference sequences.
#' @param ref_Ns A data frame containing reference sequence information.
#' @param query_Ns A data frame containing query sequence information.
#' @param plot_df A data frame containing the plot data.
#' @param sz A numeric value representing the size of the points or lines in the plot (default: 0.009).
#' @param line A logical value indicating whether to display alignments as lines (default: FALSE).
#' @param ns A logical value indicating whether to display reference and query sequence lines (default: TRUE).
#' @param gaps A data frame containing gap information, or NULL if not provided (default: NULL).
#' @param min_mapQ A numeric value representing the minimum mapping quality for filtering (default: 0).
#' @param gapwidth A numeric value representing the minimum gap width (default: 1e5).
#' @param xlim A numeric vector representing the x-axis limits, or NULL if not provided (default: NULL).
#' @param ylim A numeric vector representing the y-axis limits, or NULL if not provided (default: NULL).
#'
#' @return A ggplot object representing the genome alignment plot.
#' @export
#'
#' @examples
#' # Load sample data and call the plot_genome_alignment() function
#' # (Note: Replace with your own data when using this function)
#' plot_genome_alignment(q_name = "Query", r_name = "Target", fq = fq_data, fr = fr_data,
#' ref_Ns = ref_Ns_data, query_Ns = query_Ns_data, plot_df = plot_df_data)

plot_genome_alignment <- function(q_name, r_name, fq, fr, ref_Ns, query_Ns, plot_df,
                                  sz = 0.009, line = FALSE, ns = TRUE, gaps = NULL,
                                  min_mapQ = 0, gapwidth = 1e5,
                                  xlim = NULL, ylim = NULL) {
  plot_df <- plot_df |>
    filter(mapQ >= min_mapQ)

  p_df <- plot_df %>%
    filter(refID %in%fr,
           queryID %in% fq)

  order_size <-  p_df |>
    group_by(queryID, refID) |>
    summarise(max_aln = max(lenAln),
              query_max = max(queryEnd )) |>
    group_by(refID) |>
    top_n(1, max_aln)

  plot <- p_df %>%
    ggplot(aes(color = mapQ))

  if(!is.null(gaps)) {
    df_lgap <- gaps  |>
      data.frame() |>
      filter(width >= gapwidth) |>
      mutate(refId = seqnames,
             refStart = start,
             refEnd = end ) |>
      filter(refId %in% fr)

    df_lgap <- left_join(df_lgap, order_size, by = c('refId' = 'refID'))
    # we adjust the maximum values depending on if there are some ylim
    if (!is.null(ylim)) {
      plot <- plot +
        geom_rect(data = df_lgap, aes(xmin = refStart, ymin = ylim[1], xmax = refEnd, ymax = ylim[2]), fill = 'grey70', color = NA)
    } else {
      plot <- plot +
        geom_rect(data = df_lgap, aes(xmin = refStart, ymin = 0, xmax = refEnd, ymax = query_max), fill = 'grey70', color = NA)
    }
  }

  if (ns) {
    pr_Ns <- ref_Ns %>%
      mutate(refStart = start,
             refEnd = end,
             refID = chr) %>%
      filter(refID %in% fr)

    pq_Ns <- query_Ns %>%
      mutate(queryStart = start,
             refEnd = end,
             queryID = chr) %>%
      filter(queryID %in% fq)

    plot <- plot +
      geom_hline(data = pq_Ns, aes(yintercept = queryStart), alpha = 0.2, color = '#0074D9') +
      geom_vline(data = pr_Ns, aes(xintercept = refStart), alpha = 0.2, color = '#FF4136')
  }

  if (line) {
    plot <- plot +
      geom_segment(mapping = aes(x = refStart, y = queryStart, xend = refEnd, yend = queryEnd, color = mapQ),
                   size = sz)
  } else {
    plot <- plot +
      geom_point(mapping = aes(x = refStart, y = queryStart),
                 size = sz) +
      geom_point(mapping = aes(x = refEnd, y = queryEnd),
                 size = sz)
  }

  plot <- plot +
    facet_grid(queryID ~ refID, scales = 'free', space = 'free') +
    theme_bw() +
    theme(text = element_text(size = 8)) +
    theme(
      axis.text.y = element_text(size = 10, angle = 15),
      axis.text.x = element_text(size = 10, angle = 15),
      legend.position = 'bottom',
      strip.background = element_rect(fill = NA),
      text = element_text(size = 12)
    ) +
    scale_color_viridis_b(direction = -1) +
    labs(color = "Mapping Quality",
         title = paste0(   paste0("Post-filtering number of aln: ", nrow(p_df), " "),
                           paste0("minimum alignment length (-m): ", opt$min_align)
         )) +
    xlab(r_name) +
    ylab(q_name)

  if (ns) {
    plot <- plot +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()
      )
  } else {
    plot <- plot +
      theme(
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()
      )
  }

  if (!is.null(xlim)) {
    plot <- plot + scale_x_continuous(breaks = br, labels = br/1e6, expand = c(0, 0), limits = xlim)
  } else {
    plot <- plot + scale_x_continuous(breaks = br, labels = br/1e6, expand = c(0, 0))
  }

  if (!is.null(ylim)) {
    plot <- plot +  scale_y_continuous(breaks = br, labels = br/1e6, expand = c(0, 0), limits = ylim)
  } else {
    plot <- plot +  scale_y_continuous(breaks = br, labels = br/1e6, expand = c(0, 0))
  }

  return(plot)
}
