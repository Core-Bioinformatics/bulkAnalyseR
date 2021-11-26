volcano_plot <- function(
  genes.de.results,
  pval.threshold = 0.05, 
  lfc.threshold = 0.5,
  alpha = 0.1,
  xlims = NULL,
  log10pval.cap = TRUE,
  add.colours = TRUE,
  add.expression.colour.gradient = TRUE,
  add.guide.lines = TRUE,
  add.labels.auto = TRUE,
  add.labels.custom = FALSE,
  ...
){
  df = genes.de.results %>%
    dplyr::mutate(gene = gene_name, log10pval = log10(pvalAdj)) %>%
    dplyr::filter(!is.na(log10pval))
  
  if(all(df$log10pval >= -10)) log10pval.cap <- FALSE
  if(log10pval.cap) df$log10pval[df$log10pval < -10] <- -10
  
  vp <- ggplot(data = df, 
               mapping = aes(x = log2FC, y = -log10pval)) +
    ggplot2::theme_minimal() +
    xlab("log2(FC)") +
    ylab("-log10(pval)") 
  
  if(is.null(xlims)){
    max.abs.lfc = max(abs(df[df$log10pval > -Inf,]$log2FC))
    vp <- vp + xlim(-max.abs.lfc, max.abs.lfc)
  }else{
    vp <- vp + xlim(-abs(xlims), abs(xlims))
  }
  
  if(log10pval.cap){
    vp <- vp + scale_y_continuous(labels=c("0.0", "2.5", "5.0", "7.5", ">10"))
  }
  
  if(any(add.colours, 
         add.expression.colour.gradient,
         add.guide.lines,
         add.labels.auto,
         add.labels.custom)){
    vp <- volcano_enhance(
      vp = vp,
      df = df,
      pval.threshold = pval.threshold,
      lfc.threshold = lfc.threshold,
      alpha = alpha,
      add.colours = add.colours,
      add.expression.colour.gradient = add.expression.colour.gradient,
      add.guide.lines = add.guide.lines,
      add.labels.auto = add.labels.auto,
      add.labels.custom = add.labels.custom,
      ...
    )
  }
  
  return(vp)
  
}

volcano_enhance <- function(
  vp,
  df,
  pval.threshold,
  lfc.threshold,
  alpha,
  add.colours,
  point.colours = c("#bfbfbf", "orange", "red", "blue"),
  raster = FALSE,
  add.expression.colour.gradient,
  colour.gradient.scale = list(left  = c("#99e6ff", "#000066"),
                               right = c("#99e6ff", "#000066")),
  colour.gradient.breaks = waiver(),
  colour.gradient.limits = NULL,
  add.guide.lines,
  guide.line.colours = c("green", "blue"),
  add.labels.auto,
  add.labels.custom,
  annotation = NULL,
  n.labels.auto = c(5, 5),
  genes.to.label = NULL,
  seed = 0,
  label.force = 1
){
  
  logp.threshold = log10(pval.threshold)
  
  if(add.colours){
    colours = vector(length=nrow(df))
    colours[] = point.colours[1]
    colours[abs(df$log2FC) > lfc.threshold] = point.colours[2]
    colours[df$log10pval < logp.threshold] = point.colours[3]
    colours[abs(df$log2FC) > lfc.threshold & df$log10pval < logp.threshold] = point.colours[4]
    df$colours <- colours
    
    if(raster){
      vp <- vp + ggrastr::rasterise(geom_point(alpha = alpha, colour = colours))
    }else{
      vp <- vp + geom_point(alpha = alpha, colour = colours)
    }
  }
  
  if(add.expression.colour.gradient){
    df.colour.gradient <- df %>%
      dplyr::filter(abs(log2FC) > lfc.threshold & log10pval < logp.threshold) %>%
      dplyr::arrange(log2exp)
    if(identical(colour.gradient.scale$left, colour.gradient.scale$right)){
      vp <- vp +
        geom_point(data = df.colour.gradient,
                   mapping = aes(x = log2FC, y = -log10pval, colour = log2exp)) +
        scale_color_gradient(low = colour.gradient.scale$left[1], 
                             high = colour.gradient.scale$left[2],
                             breaks = colour.gradient.breaks,
                             limits = colour.gradient.limits) +
        labs(colour = "log2(exp)")
    }else{
      vp <- vp +
        geom_point(data = dplyr::filter(df.colour.gradient, log2FC < 0),
                   mapping = aes(x = log2FC, y = -log10pval, colour = log2exp)) +
        scale_color_gradient(low = colour.gradient.scale$left[1], 
                             high = colour.gradient.scale$left[2],
                             breaks = colour.gradient.breaks,
                             limits = colour.gradient.limits) +
        labs(colour = "log2(exp)") +
        ggnewscale::new_scale_colour() +
        geom_point(data = dplyr::filter(df.colour.gradient, log2FC > 0),
                   mapping = aes(x = log2FC, y = -log10pval, colour = log2exp)) +
        scale_colour_gradient(low = colour.gradient.scale$right[1], 
                              high = colour.gradient.scale$right[2],
                              breaks = colour.gradient.breaks,
                              limits = colour.gradient.limits) +
        labs(colour = "log2(exp)")
    }
  }
  
  if(add.guide.lines){
    vp <- vp +
      geom_vline(xintercept =      lfc.threshold,  colour = guide.line.colours[1]) +
      geom_vline(xintercept =     -lfc.threshold,  colour = guide.line.colours[1]) +
      geom_vline(xintercept =  2 * lfc.threshold,  colour = guide.line.colours[2]) +
      geom_vline(xintercept = -2 * lfc.threshold,  colour = guide.line.colours[2]) +
      geom_hline(yintercept =     -logp.threshold, colour = guide.line.colours[1]) +
      geom_hline(yintercept = -2 * logp.threshold, colour = guide.line.colours[2])
  }
  
  if(add.labels.auto | add.labels.custom){
    if(!is.null(annotation)){
      df <- df %>% 
        dplyr::mutate(
          symbol = annotation$SYMBOL[match(gene, annotation$ENSEMBL)],
          name = ifelse(is.na(symbol), gene, symbol)
        ) %>%
        dplyr::select(-symbol)
    }else{
      df <- df %>% dplyr::mutate(name = gene)
    }
    
    df.label <- tibble::tibble()
    if(add.labels.custom){
      genes.to.rename <- genes.to.label[names(genes.to.label) != ""]
      genes.to.label <- df$name[(match(genes.to.label, c(df$name, df$gene)) - 1) %% nrow(df) + 1]
      genes.to.label <- unique(genes.to.label[!is.na(genes.to.label)])
      genes.to.rename <- genes.to.rename[genes.to.rename %in% genes.to.label]
      df.label <- dplyr::filter(df, name %in% genes.to.label)
      df.label$name[match(genes.to.rename, df.label$name)] <- names(genes.to.rename)
      if(nrow(df.label) == 0){
        message(paste0("add.labels.custom was TRUE but no genes specified; ",
                       "did you forget to supply genes.to.label or annotation?"))
      }
    }
    
    if(add.labels.auto){
      if(length(n.labels.auto) == 1) n.labels.auto <- rep(n.labels.auto, 2)
      df.significant <- dplyr::filter(df, 
                                      abs(log2FC) > lfc.threshold,
                                      log10pval < logp.threshold,
                                      !(name %in% genes.to.label))
      df.lowest.p.vals <- head(df.significant, n.labels.auto[1])
      df.rest <- tail(df.significant, nrow(df.significant) - n.labels.auto[1])
      df.highest.abn <- head(df.rest[order(df.rest$log2exp, decreasing=TRUE),], n.labels.auto[2])
      df.label <- rbind(df.lowest.p.vals, df.highest.abn, df.label) %>%
        dplyr::distinct(name, .keep_all = TRUE)
    }
    
    set.seed(seed = seed)
    vp <- vp +
      ggrepel::geom_label_repel(data = df.label, 
                                mapping = aes(x = log2FC, y = -log10pval, label = name),
                                max.overlaps = nrow(df.label),
                                force = label.force,
                                point.size = NA)
  }
  
  return(vp)
  
}

ma_plot <- function(
  genes.de.results,
  pval.threshold = 0.05, 
  lfc.threshold = 0.5,
  alpha = 0.1,
  ylims = NULL,
  add.colours = TRUE,
  add.expression.colour.gradient = TRUE,
  add.guide.lines = TRUE,
  add.labels.auto = TRUE,
  add.labels.custom = FALSE,
  ...
){
  df = genes.de.results %>%
    dplyr::mutate(gene = gene_name, log10pval = log10(pvalAdj)) %>%
    dplyr::filter(!is.na(log10pval))
  
  p <- ggplot(data = df, 
                mapping = aes(x = log2exp, y = log2FC)) +
    ggplot2::theme_minimal() +
    xlab("Average log2(exp)") +
    ylab("log2(FC)")
  
  if(is.null(ylims)){
    max.abs.lfc = max(abs(df$log2FC))
    p <- p + ylim(-max.abs.lfc, max.abs.lfc)
  }else{
    p <- p + ylim(-abs(ylims), abs(ylims))
  }
  
  if(any(add.colours, 
         add.expression.colour.gradient,
         add.guide.lines,
         add.labels.auto,
         add.labels.custom)){
    p <- ma_enhance(
      p = p,
      df = df,
      pval.threshold = pval.threshold,
      lfc.threshold = lfc.threshold,
      alpha = alpha,
      add.colours = add.colours,
      add.expression.colour.gradient = add.expression.colour.gradient,
      add.guide.lines = add.guide.lines,
      add.labels.auto = add.labels.auto,
      add.labels.custom = add.labels.custom,
      ...
    )
  }
  
  return(p)
  
}

ma_enhance <- function(
  p,
  df,
  pval.threshold,
  lfc.threshold,
  alpha,
  add.colours,
  point.colours = c("#bfbfbf", "orange", "red", "blue"),
  raster = FALSE,
  add.expression.colour.gradient,
  colour.gradient.scale = list(left  = c("#99e6ff", "#000066"),
                               right = c("#99e6ff", "#000066")),
  colour.gradient.breaks = waiver(),
  colour.gradient.limits = NULL,
  add.guide.lines,
  guide.line.colours = c("green", "blue"),
  add.labels.auto,
  add.labels.custom,
  annotation = NULL,
  n.labels.auto = c(5, 5),
  genes.to.label = NULL,
  seed = 0,
  label.force = 1
){
  
  logp.threshold = log10(pval.threshold)
  
  if(add.colours){
    colours = vector(length=nrow(df))
    colours[] = point.colours[1]
    colours[abs(df$log2FC) > lfc.threshold] = point.colours[2]
    colours[df$log10pval < logp.threshold] = point.colours[3]
    colours[abs(df$log2FC) > lfc.threshold & df$log10pval < logp.threshold] = point.colours[4]
    df$colours <- colours
    
    if(raster){
      p <- p + ggrastr::rasterise(geom_point(alpha = alpha, colour = colours))
    }else{
      p <- p + geom_point(alpha = alpha, colour = colours)
    }
  }
  
  if(add.expression.colour.gradient){
    df.colour.gradient <- df %>%
      dplyr::filter(abs(log2FC) > lfc.threshold & log10pval < logp.threshold) %>%
      dplyr::arrange(log2exp)
    if(identical(colour.gradient.scale$left, colour.gradient.scale$right)){
      p <- p +
        geom_point(data = df.colour.gradient,
                   mapping = aes(x = log2exp, y = log2FC, colour = log2exp)) +
        scale_color_gradient(low = colour.gradient.scale$left[1], 
                             high = colour.gradient.scale$left[2],
                             breaks = colour.gradient.breaks,
                             limits = colour.gradient.limits) +
        labs(colour = "log2(exp)")
    }else{
      p <- p +
        geom_point(data = dplyr::filter(df.colour.gradient, log2FC < 0),
                   mapping = aes(x = log2exp, y = log2FC, colour = log2exp)) +
        scale_color_gradient(low = colour.gradient.scale$left[1], 
                             high = colour.gradient.scale$left[2],
                             breaks = colour.gradient.breaks,
                             limits = colour.gradient.limits) +
        labs(colour = "log2(exp)") +
        ggnewscale::new_scale_colour() +
        geom_point(data = dplyr::filter(df.colour.gradient, log2FC > 0),
                   mapping = aes(x = log2exp, y = log2FC, colour = log2exp)) +
        scale_colour_gradient(low = colour.gradient.scale$right[1], 
                              high = colour.gradient.scale$right[2],
                              breaks = colour.gradient.breaks,
                              limits = colour.gradient.limits) +
        labs(colour = "log2(exp)")
    }
  }
  
  if(add.guide.lines){
    p <- p +
      geom_hline(yintercept =      lfc.threshold,  colour = guide.line.colours[1]) +
      geom_hline(yintercept =     -lfc.threshold,  colour = guide.line.colours[1]) +
      geom_hline(yintercept =  2 * lfc.threshold,  colour = guide.line.colours[2]) +
      geom_hline(yintercept = -2 * lfc.threshold,  colour = guide.line.colours[2])
  }
  
  if(add.labels.auto | add.labels.custom){
    if(!is.null(annotation)){
      df <- df %>% 
        dplyr::mutate(
          symbol = annotation$SYMBOL[match(gene, annotation$ENSEMBL)],
          name = ifelse(is.na(symbol), gene, symbol)
        ) %>%
        dplyr::select(-symbol)
    }else{
      df <- df %>% dplyr::mutate(name = gene)
    }
    
    df.label <- tibble::tibble()
    if(add.labels.custom){
      genes.to.rename <- genes.to.label[names(genes.to.label) != ""]
      genes.to.label <- df$name[(match(genes.to.label, c(df$name, df$gene)) - 1) %% nrow(df) + 1]
      genes.to.label <- unique(genes.to.label[!is.na(genes.to.label)])
      genes.to.rename <- genes.to.rename[genes.to.rename %in% genes.to.label]
      df.label <- dplyr::filter(df, name %in% genes.to.label)
      df.label$name[match(genes.to.rename, df.label$name)] <- names(genes.to.rename)
      if(nrow(df.label) == 0){
        message(paste0("add.labels.custom was TRUE but no genes specified; ",
                       "did you forget to supply genes.to.label or annotation?"))
      }
    }
    
    if(add.labels.auto){
      if(length(n.labels.auto) == 1) n.labels.auto <- rep(n.labels.auto, 2)
      df.significant <- dplyr::filter(df, 
                                      abs(log2FC) > lfc.threshold,
                                      log10pval < logp.threshold,
                                      !(name %in% genes.to.label))
      df.lowest.p.vals <- head(df.significant, n.labels.auto[1])
      df.rest <- tail(df.significant, nrow(df.significant) - n.labels.auto[1])
      df.highest.abn <- head(df.rest[order(df.rest$log2exp, decreasing=TRUE),], n.labels.auto[2])
      df.label <- rbind(df.lowest.p.vals, df.highest.abn, df.label) %>%
        dplyr::distinct(name, .keep_all = TRUE)
    }
    
    set.seed(seed = seed)
    p <- p +
      ggrepel::geom_label_repel(data = df.label, 
                                mapping = aes(x = log2exp, y = log2FC, label = name),
                                max.overlaps = nrow(df.label),
                                force = label.force,
                                point.size = NA)
  }
  
  return(p)
  
}