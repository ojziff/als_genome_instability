# iPSN ALS meta-analysis Gene Expression Browser
library(here)
library(tidyverse)
library(circlize)
library(GetoptLong)
library(shiny)
library(shinydashboard)
library(shinythemes)
library(bs4Dash)
library(DT)
library(DESeq2)
library(bslib)
library(InteractiveComplexHeatmap)
library(ComplexHeatmap)
library(DESeq2)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)
library(vsn)
library(ggrepel)
library(ggpubr)
library(ggsci)
library(ggplotify)
library(decoupleR)

# Load data ----------------
camp_path = here("/camp/lab/luscomben/home/users/ziffo")
proj_path = here::here("/camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta")
load(here::here(proj_path,"expression/deseq2/shiny_dds.RData"))
gene2ens <- readRDS(here(camp_path,"/genomes/ensembl/Homo_sapiens.GRCh38.99.gene2ens.RDS")) %>% unique
net_progeny <- get_progeny(organism = 'human', top = 100)
net_dorothea <- get_dorothea(organism='human', levels=c('A', 'B', 'C'))

data_list = list("ipsc_mn_als_datasets.dds"=ipsc_mn_als_datasets.dds,"ipsc_mn_sporadic_datasets.dds"=ipsc_mn_sporadic_datasets.dds,"ipsc_mn_c9orf72_datasets.dds"=ipsc_mn_c9orf72_datasets.dds,"ipsc_mn_fus_datasets.dds"=ipsc_mn_fus_datasets.dds,
                 "ipsc_mn_tardbp_datasets.dds"=ipsc_mn_tardbp_datasets.dds,"ipsc_mn_sod1_datasets.dds"=ipsc_mn_sod1_datasets.dds)

# Functions ----------------------
theme_oz <- function () { 
  theme_bw(base_size=12) %+replace% # theme_bw(base_size=8, base_family="Helvetica") %+replace% 
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      #panel.border = element_blank(),
      axis.line = element_line(),
      #text = element_text(color = "black"), 
      strip.text = element_text(color = "black"),
      axis.text = element_text(colour = "black"),
      axis.ticks = element_line(colour = "black"),
      panel.background  = element_blank(),
      plot.background = element_rect(fill="white", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA), 
      axis.line.y = element_line(), strip.text.x = element_text(face = "bold", margin = margin(t = 2,r = 0,b = 2,l=0))
    )
}

tibble_to_matrix <- function(tbl, ..., row_names=NULL){
  cols <- rlang::enquos(...)
  mat <- as.matrix(dplyr::select(tbl, !!! cols))
  if (!is.null(row_names)){
    if (length(row_names) == 1 & row_names[1] %in% colnames(tbl)){
      row_names <- tbl[[row_names]]
    }
    rownames(mat) <- row_names
  }
  return(mat)
}

ma_plot <- function(deseq_res, ymax = 5, xmax = 100000, significance_threshold = "padj"){
  deseq_res <- deseq_res %>% mutate(sig = case_when(!!sym(significance_threshold) < 0.05 & abs(log2FoldChange) < 1 ~ "sig", !!sym(significance_threshold) < 0.05 & abs(log2FoldChange) >= 1 ~ "sig_strong", !!sym(significance_threshold) >= 0.05 ~ "non_sig"), 
                                    direction = ifelse(log2FoldChange > 0, "up", "down"), class = paste(sig, direction), log2FoldChange = case_when(log2FoldChange > xmax ~ Inf, log2FoldChange < -xmax ~ -Inf, TRUE ~ log2FoldChange), pvalue = case_when(-log10(pvalue) > ymax ~ 10^-ymax, TRUE ~ pvalue))
  plot <- ggplot(data = deseq_res, aes(x = baseMean, y = log2FoldChange)) +  # NB if pvalue == NA then will not show
    ggrastr::rasterise(geom_point(aes(colour = class), size = 0.5), dpi = 72) +
    scale_colour_manual(values = c("non_sig up" = "gray", "non_sig down" = "gray", 
                                   "sig up" = "#F36F6F",
                                   "sig_strong up" = "#EB4445",
                                   "sig down" = "#A6CEE3",#"#4F8FC4",
                                   "sig_strong down" = "#79B1D3")) +
    labs(x = expression(log[2]~base~mean~expression), y = paste0("log2 fold change ALS vs CTRL")) + guides(colour = "none") +
    scale_y_continuous(expand = c(0,0), limits = c(-ymax,ymax)) +scale_x_continuous(limits = c(0,xmax)) +
    theme_oz() +  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), panel.border = element_blank(), axis.ticks = element_line(colour = "black") ) +
    geom_hline(yintercept = 0, linetype = 'dotted') + 
    geom_text_repel(fontface = "italic", data = top_n(deseq_res, n = 30, wt = abs(stat)), aes(x = baseMean, y = log2FoldChange, label = gene_name), max.overlaps = Inf, min.segment.length = unit(0, "lines"),size = 2.8)

  return(plot)
}


volcano_plot <- function(deseq_res, ymax = 16.5, xmax = 3, significance_threshold = "padj"){
  deseq_res <- deseq_res %>% mutate(sig = case_when(!!sym(significance_threshold) < 0.05 & abs(log2FoldChange) < 1 ~ "sig", !!sym(significance_threshold) < 0.05 & abs(log2FoldChange) >= 1 ~ "sig_strong", !!sym(significance_threshold) >= 0.05 ~ "non_sig"), direction = ifelse(log2FoldChange > 0, "up", "down"), class = paste(sig, direction), 
                                    log2FoldChange = case_when(log2FoldChange > xmax ~ Inf, log2FoldChange < -xmax ~ -Inf, TRUE ~ log2FoldChange), pvalue = case_when(-log10(pvalue) > ymax ~ 10^-ymax, TRUE ~ pvalue))
  plot <- ggplot(deseq_res, aes(x = log2FoldChange, y = -log10(pvalue))) +  
    ggrastr::rasterise(geom_point(aes(colour = class), size = 0.5), dpi = 72) +
    scale_colour_manual(values = c("non_sig up" = "gray", "non_sig down" = "gray", 
                                   "sig up" = "#F36F6F",
                                   "sig_strong up" = "#EB4445",
                                   "sig down" = "#A6CEE3",#"#4F8FC4",
                                   "sig_strong down" = "#79B1D3")) +
    labs(y = expression(-log[10]~P~value), x = paste0("log2 fold change ALS vs CTRL")) + guides(colour = "none") +
    scale_y_continuous(expand = c(0,0), limits = c(0,ymax)) +
    theme_oz() +  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), panel.border = element_blank(), axis.ticks = element_line(colour = "black") ) +
    scale_x_continuous(limits = c(-xmax,xmax)) +
    geom_text_repel(fontface = "italic", data = top_n(deseq_res, n = 30, wt = abs(stat)), aes(x = log2FoldChange, y = -log10(pvalue), label = gene_name), max.overlaps = Inf, min.segment.length = unit(0, "lines"),size = 2.8)
  return(plot)
}

make_heatmap = function(deseq_dds, deseq_res, fdr = 0.05, base_mean = 0, log2fc = 0) {
  rownames(deseq_dds) = gene2ens$gene_name[gene2ens$gene_id %in% rownames(deseq_dds)]
  heatmap_sample_column_order <- colData(deseq_dds) %>% as_tibble %>% arrange(condition)
  l = deseq_res$padj <= fdr & deseq_res$baseMean >= base_mean &  abs(deseq_res$log2FoldChange) >= log2fc; l[is.na(l)] = FALSE
    if(sum(l) == 0) return(NULL)
    m = counts(deseq_dds, normalized = TRUE)
    m = m[l, ]

    ht = as.ggplot(
      Heatmap(t(scale(t(m))), name = "normalised counts", border = TRUE,
              column_order = heatmap_sample_column_order$dataset_sample,
              top_annotation = HeatmapAnnotation(condition = colData(deseq_dds)$condition, col = list(condition = c("ctrl"="dodgerblue2", "als"="firebrick2")), dataset = colData(deseq_dds)$dataset),
              row_names_gp = gpar(fontsize = 6), show_column_names = FALSE, #row_km = 2,
              column_title = paste0(sum(l), " significant genes with FDR < ", fdr), cluster_columns = FALSE, cluster_rows = FALSE, show_row_dend = FALSE, col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))))
    return(ht)
}


# ui user interface ---------------------------
ui <- dashboardPage(
  
  dashboardHeader(title = "iPSN ALS Meta-analysis", titleWidth = 250),
  
  dashboardSidebar(
    # fluidPage(theme = bs_theme(version = 4, bootswatch = "minty")),   # tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")),
    width = 300, skin = "light",
    bs4SidebarMenu(      
      bs4SidebarMenuItem("Select comparison",  icon = icon("cog"), startExpanded = TRUE,
                         radioButtons("Select comparison", NULL, list("Pan-ALS vs CTRL" = "ipsc_mn_als_datasets.dds", "Sporadic vs CTRL" = "ipsc_mn_sporadic_datasets.dds", "C9orf72 vs CTRL" = "ipsc_mn_c9orf72_datasets.dds", "FUS vs CTRL" = "ipsc_mn_fus_datasets.dds", "TARDBP vs CTRL" = "ipsc_mn_tardbp_datasets.dds","SOD1 vs CTRL" = "ipsc_mn_sod1_datasets.dds"), 
                                      inline = FALSE)),
      bs4SidebarMenuItem("Results Table", tabName = "table", icon = icon("table")),
      bs4SidebarMenuItem("Visualisation", tabName = "vis", icon = icon("chart-area")),
      bs4SidebarMenuItem("Pathways & TFs", tabName = "pathway_tf", icon = icon("chart-bar")),
      bs4SidebarMenuItem("Heatmap", tabName = "heatmap", icon = icon("th")))),
  
  dashboardBody(
    bs4TabItems(
      bs4TabItem(tabName="table", DTOutput("resTable") ),
      bs4TabItem(tabName = "vis",
                 fluidRow( column(10, box(title = "MA plot", plotOutput("ma"), width = 12, style='height:35vw' ) ) ),
                 fluidRow( column(10, box(title = "Volcano plot", plotOutput("volcano"), width = 12, style='height:35vw' ) ) )
                 ),
      bs4TabItem(tabName = "pathway_tf",
                 fluidRow( column(10, box(title = "Signalling Pathways", plotOutput("progeny"), width = 12, style='height:35vw' ) ) ),
                 fluidRow( column(10, box(title = "Transcription Factors", plotOutput("dorothea"), width = 12, style='height:35vw' ) ) )
      ),
      bs4TabItem(tabName = "heatmap", plotOutput("ht", height = 1000px) ) 
      )
    ) 
  )



# server --------------------------
server <- function(input, output, session) { #
  
  datasetInput <- reactive({
    df <- data_list[[input$`Select comparison`]]
  })
  
  output$resTable <- renderDT({
    dds <- datasetInput() ## Note: added assignment
    res = DESeq2::results(dds) %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue) #%>% column_to_rownames(var = "gene_id")
    formatRound(datatable(res, rownames = FALSE), columns = 2:7, digits = 3)
  })
  
  output$ma <- renderPlot({
    dds <- datasetInput() ## Note: added assignment
    res = DESeq2::results(dds) %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue) #%>% column_to_rownames(var = "gene_id")
    ma = ma_plot(res)
    ma
    })
  
  output$volcano <- renderPlot({
    dds <- datasetInput() ## Note: added assignment
    res = DESeq2::results(dds) %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue) #%>% column_to_rownames(var = "gene_id")
    volcano = volcano_plot(res)
    volcano
  })
  
  output$ht <- renderPlot({
    dds <- datasetInput() ## Note: added assignment
    res = DESeq2::results(dds) %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue) #%>% column_to_rownames(var = "gene_id")
    ht = make_heatmap(dds, res)
    ht
  })
  
  output$progeny <- renderPlot({
    dds <- datasetInput() ## Note: added assignment
    res.matrix = DESeq2::results(dds) %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue) %>% select(gene_name, stat) %>% drop_na(stat) %>% distinct(gene_name, .keep_all = TRUE) %>% 
      tibble_to_matrix(stat, row_names = "gene_name")
    progeny.contrast_acts <- run_wmean(mat=res.matrix, net=net_progeny, .source='source', .target='target', .mor='weight', times = 100, minsize = 5) %>% filter(statistic == 'norm_wmean') %>% 
      mutate(p.signif = case_when(p_value < 0.003 ~ "***", p_value < 0.009 ~ "***",p_value < 0.01 ~ "**",p_value < 0.05 ~ "*", TRUE ~ ""), group1 = source, group2 = source)
    progeny <- ggplot(progeny.contrast_acts, aes(x=reorder(source, score), y = score)) +
      geom_bar(aes(fill = score), stat = "identity") +
      scale_fill_gradient2(low = "darkblue", high = "indianred", mid = "whitesmoke", midpoint = 0) + 
      theme_oz() + theme(axis.text.x = element_text(angle = 90), legend.position = "none", plot.title = element_text(hjust = 0.5)) +
      labs(title = "Signaling Pathways", x = "", y = "ALS vs CTRL Normalised enrichment") +
      stat_pvalue_manual(progeny.contrast_acts, label = "p.signif", y.position = 14, xmin = "source", xmax = NULL, size = 3, hide.ns = TRUE) 
    progeny
  })
  
  output$dorothea <- renderPlot({
    dds <- datasetInput() ## Note: added assignment
    res.matrix = DESeq2::results(dds) %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue) %>% select(gene_name, stat) %>% drop_na(stat) %>% distinct(gene_name, .keep_all = TRUE) %>% 
      tibble_to_matrix(stat, row_names = "gene_name")
    dorothea.contrast_acts <- run_wmean(mat=res.matrix, net=net_dorothea, .source='source', .target='target', .mor='mor', times = 100, minsize = 5) %>% filter(statistic == 'norm_wmean') %>% 
      mutate(p.signif = case_when(p_value < 0.05 ~ "*", TRUE ~ ""), sig = case_when(p_value < 0.135 & abs(score) < 4 ~ "sig", p_value < 0.135 & abs(score) >= 4 ~ "sig_strong", p_value >= 0.135 ~ "non_sig"), 
             direction = ifelse(score > 0, "up", "down"), class = paste(sig, direction))
    dorothea.contrast_acts.labels.up = dorothea.contrast_acts %>% top_n(score, n = 10) %>% pull(source)
    dorothea.contrast_acts.labels.down = dorothea.contrast_acts %>% top_n(-score, n = 10) %>% pull(source)
    dorothea = ggplot(dorothea.contrast_acts, aes(x = score, y = -log10(p_value))) + geom_point(aes(colour = class)) +
      scale_colour_manual(values = c("non_sig up" = "gray", "non_sig down" = "gray",  "sig up" = "#F36F6F", "sig_strong up" = "#EB4445", "sig down" = "#A6CEE3", "sig_strong down" = "#79B1D3")) +
      labs(y = expression(-log[10]~P~value), x = "Normalised enrichment ALS vs CTRL", title = "Transcription Factors") +
      guides(colour = "none") +
      theme_oz() +  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), panel.border = element_blank(), axis.ticks = element_line(colour = "black") ) +
      geom_text_repel(fontface = "italic", data = filter(dorothea.contrast_acts, source %in% c(dorothea.contrast_acts.labels.up, dorothea.contrast_acts.labels.down)), aes(label = source), max.overlaps = Inf, min.segment.length = unit(0, "lines"), size = 2.8) +  
      geom_vline(xintercept = 0, linetype = 'dotted') + ylim(c(0,2.8))
    dorothea
  })
  
}


# generate app -----------------------
shinyApp(ui, server)

